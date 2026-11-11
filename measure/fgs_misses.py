"""What FGS's one pass misses against hmmsearch on all ORFs, and what a second pass could recover.

fgs_misses.py prep READS.fasta CALLS.faa CALLS.ffn ORFS.domtbl DECOYS.domtbl CALLS.domtbl MODELS.tsv OUT_DIR
    Classifies each miss. Writes, per removal variant, the read pieces left for FGS, the best case
    (the miss's own frame, translated from what survives) and the models to search it against.
fgs_misses.py score ...same arguments...
    Second pass: why each miss is missed again, given OUT_DIR/{variant}.calls.faa, .calls.domtbl and .oracle.domtbl.
Variants: call = the whole FGS call removed; ga = its alignments at >= GA removed; any = all its alignments removed.
"""

import csv
import sys
from collections import Counter, defaultdict
from statistics import median

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from orf_filters import NAMES, bin_of, fasta, hexamers, log_odds, placed  # noqa: E402
from real_orfs import RC  # noqa: E402

VARIANTS = ("call", "ga", "any")
FGS_MIN = 60  # FGS writes a call only when longer than this, viterbi.rs:577


def domains(path):
    with open(path) as f:
        for line in f:
            if line[0] != "#":
                c = line.split()
                yield c[0], c[4].split(".")[0], float(c[13]), float(c[12]), int(c[17]), int(c[18])


def call_name(name):
    r, a, b, st = name.rsplit("_", 3)
    return r, int(a), int(b), st == "-"


def call_nt(a, b, minus, f, t):
    """A call's residues f..t on the read, half open."""
    return (b - 3 * t, b - 3 * (f - 1)) if minus else (a - 1 + 3 * (f - 1), a - 1 + 3 * t)


def frames_of(a, b, minus, off=0):
    return {(minus, (off + a - 1) % 3), (minus, (off + b) % 3)}


def remainder(n, cut):
    """What is left of 0..n after removing the intervals."""
    out, at = [], 0
    for s, e in sorted(cut):
        if s > at:
            out.append((at, s))
        at = max(at, e)
    if at < n:
        out.append((at, n))
    return out


def overlap(a, b, c, d):
    return max(0, min(b, d) - max(a, c))


def load(reads_fa, calls_faa, calls_ffn, table, decoys, calls_tbl, models):
    d = {"reads": dict(fasta(reads_fa))}
    with open(models) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    d["ga"] = {r["acc"]: float(r["ga"]) for r in rows}
    d["name"] = {r["acc"]: r["name"] for r in rows}
    d["orfs"] = {}
    for r, s in d["reads"].items():
        for j, (p, a, b, st) in enumerate(placed(s)):
            d["orfs"][f"{r}_{j}"] = (r, p, a, b, bool(st))
    d["calls"] = {}
    for name, seq in fasta(calls_faa):
        r, a, b, minus = call_name(name)
        d["calls"][r] = (a, b, minus, seq)
    truth = {}
    for q, m, s, e, f, t in domains(table):
        k = (d["orfs"][q][0], m)
        if k not in truth or s > truth[k][0]:
            truth[k] = (s, e, q, f, t)
    d["truth"] = truth
    d["rev"] = {}
    for q, m, s, *_ in domains(decoys):
        k = (q, m)
        d["rev"][k] = max(d["rev"].get(k, -1e9), s)
    d["call_doms"] = defaultdict(list)
    for q, m, s, e, f, t in domains(calls_tbl):
        d["call_doms"][q.rsplit("_", 3)[0]].append((m, s, e, f, t))
    d["found"] = {}
    for r, doms in d["call_doms"].items():
        for m, s, *_ in doms:
            d["found"][(r, m)] = max(d["found"].get((r, m), -1e9), s)
    # Hexamer log-odds, cross-fitted by read number, as orf_filters.py frames
    fold = lambda r: int(r[1:]) % 2  # noqa: E731
    coding, other = [Counter(), Counter()], [Counter(), Counter()]
    for name, seq in fasta(calls_ffn):
        coding[fold(name.rsplit("_", 3)[0])].update(hexamers(seq.upper()))
    for _, (r, _, a, b, st) in d["orfs"].items():
        if r in d["calls"] and (st, a % 3) not in frames_of(*d["calls"][r][:3]):
            s = d["reads"][r][a:b]
            other[fold(r)].update(hexamers(s[::-1].translate(RC) if st else s))
    tables = [log_odds(coding[1 - i], other[1 - i]) for i in (0, 1)]

    def coding_score(r, a, b, st):
        s = d["reads"][r][a:b]
        return sum(tables[fold(r)].get(h, 0.0) for h in hexamers(s[::-1].translate(RC) if st else s))
    d["coding"] = coding_score
    return d


def hit_nt(d, hit):
    _, _, q, f, t = hit
    _, _, a, b, st = d["orfs"][q]
    return (b - 3 * t, b - 3 * (f - 1)) if st else (a + 3 * (f - 1), a + 3 * t)


def misses(d):
    out = []
    for (r, m), hit in d["truth"].items():
        b = bin_of(hit[0] - d["ga"][m])
        got = d["found"].get((r, m))
        if got is None or (b == 0 and got < d["ga"][m]):
            out.append((r, m, b, hit))
    return out


def classify(d, r, m, hit):
    q = hit[2]
    _, prot, a, b, st = d["orfs"][q]
    if r not in d["calls"]:
        return "read not called"
    ca, cb, minus, seq = d["calls"][r]
    s, e = hit_nt(d, hit)
    if (st, a % 3) in frames_of(ca, cb, minus):
        if overlap(s, e, ca - 1, cb) < 0.8 * (e - s):
            return "right frame: call stops short of the hit"
        if seq not in prot:
            return "right frame: call differs from ORF (repeated end codon, or frameshift)"
        return "right frame: score lower" if (r, m) in d["found"] else "right frame: no hit on call"
    side = "opposite strand" if st != minus else "same strand"
    doms = d["call_doms"].get(r, [])
    if any(s >= d["ga"][x] for x, s, *_ in doms):
        rival = "call has a GA hit"
    elif doms:
        rival = "call has a weak hit"
    else:
        rival = "call has no hit"
    return f"other frame, {side}: {rival}"


def report_first(d, ms):
    rows = Counter()
    for r, m, b, hit in ms:
        rows[classify(d, r, m, hit), b] += 1
    tot = Counter(bin_of(h[0] - d["ga"][m]) for (r, m), h in d["truth"].items())
    miss = Counter(b for *_, b, _ in ms)
    print(f"{'':<70}" + "".join(f"{n:>9}" for n in NAMES))
    print(f"{'truth hits':<70}" + "".join(f"{tot[i]:>9}" for i in range(5)))
    print(f"{'missed':<70}" + "".join(f"{miss[i]:>9}" for i in range(5)))
    for c in sorted({k[0] for k in rows}, key=lambda c: -sum(rows[c, i] for i in range(5))):
        print(f"{c:<70}" + "".join(f"{rows[c, i]:>9}" for i in range(5)))

    print("\nestimated precision by where the hit sits (1 - reversed/real)")
    real, fake = Counter(), Counter()

    def where(q):
        r, _, a, _, st = d["orfs"][q]
        if r not in d["calls"]:
            return "uncalled read"
        return "called frame" if (st, a % 3) in frames_of(*d["calls"][r][:3]) else "other frame"
    for (r, m), h in d["truth"].items():
        real[where(h[2]), bin_of(h[0] - d["ga"][m])] += 1
    best = {}
    for (q, m), s in d["rev"].items():
        k = (d["orfs"][q][0], m)
        if k not in best or s > best[k][0]:
            best[k] = (s, q)
    for (r, m), (s, q) in best.items():
        fake[where(q), bin_of(s - d["ga"][m])] += 1
    for w in ("called frame", "other frame", "uncalled read"):
        print(f"{w:<16}" + " ".join(f"{NAMES[i]} {real[w, i]}/{fake[w, i]}={1 - fake[w, i] / max(real[w, i], 1):.2f}" for i in range(5)))

    print("\nmissed against kept: median per bin")
    missed = {(r, m) for r, m, *_ in ms}
    for i in range(4):
        group = {True: [], False: []}
        for k, h in d["truth"].items():
            if bin_of(h[0] - d["ga"][k[1]]) == i:
                group[k in missed].append(k)
        cols = []
        for was in (False, True):
            ks = group[was]
            ali = [d["truth"][k][4] - d["truth"][k][3] + 1 for k in ks]
            orf = [len(d["orfs"][d["truth"][k][2]][1]) for k in ks]
            gc = [sum(c in "GC" for c in d["reads"][k[0]]) / len(d["reads"][k[0]]) for k in ks]
            cols.append(f"{'missed' if was else 'kept'} n={len(ks)} alignment {median(ali):.0f} aa, ORF {median(orf):.0f} aa, GC {median(gc):.2f}")
        print(f"{NAMES[i]:<8}" + " | ".join(cols))

    print("\nother-frame misses: the missed frame against the called frame")
    for i in range(4):
        dl, longer, higher = [], 0, 0
        for r, m, b, hit in ms:
            if b != i or r not in d["calls"]:
                continue
            _, p, a, e, st = d["orfs"][hit[2]]
            ca, cb, minus, seq = d["calls"][r]
            if (st, a % 3) in frames_of(ca, cb, minus):
                continue
            x = d["coding"](r, a, e, st)
            y = d["coding"](r, ca - 1, cb, minus)
            dl.append(x - y)
            longer += len(p) > len(seq)
            higher += x > y
        print(f"{NAMES[i]:<8} n={len(dl)} median coding score difference {median(dl):+.1f}; missed frame scores higher {higher / len(dl):.2f}, is longer {longer / len(dl):.2f}")

    print("\nGA recall by read GC")
    ok, n = Counter(), Counter()
    for k, h in d["truth"].items():
        if h[0] >= d["ga"][k[1]]:
            s = d["reads"][k[0]]
            g = int(10 * sum(c in "GC" for c in s) / len(s))
            n[g] += 1
            ok[g] += d["found"].get(k, -1e9) >= d["ga"][k[1]]
    print(" ".join(f"{g * 10}-{g * 10 + 10}%: {ok[g] / n[g]:.3f} (n={n[g]})" for g in sorted(n)))
    called = Counter()
    for r, s in d["reads"].items():
        g = int(10 * sum(c in "GC" for c in s) / len(s))
        called[g, r in d["calls"]] += 1
    print("reads called by GC: " + " ".join(f"{g * 10}-{g * 10 + 10}%: {called[g, True] / (called[g, True] + called[g, False]):.3f}" for g in sorted({k[0] for k in called})))

    print("\nGA misses by family, 5+ GA hits, lowest recall")
    fam, lost = Counter(), Counter()
    kinds = defaultdict(Counter)
    for (r, m), h in d["truth"].items():
        if h[0] >= d["ga"][m]:
            fam[m] += 1
            if d["found"].get((r, m), -1e9) < d["ga"][m]:
                lost[m] += 1
                kinds[m][classify(d, r, m, h).split(":")[0]] += 1
    for m in sorted((m for m in fam if fam[m] >= 5), key=lambda m: (fam[m] - lost[m]) / fam[m])[:15]:
        print(f"{m} {d['name'][m]:<24} {fam[m] - lost[m]}/{fam[m]}  {dict(kinds[m])}")


def cuts(d, variant, r):
    """Intervals a second pass removes from read r, or nothing."""
    if r not in d["calls"]:
        return []
    a, b, minus, _ = d["calls"][r]
    if variant == "call":
        return [(a - 1, b)]
    doms = d["call_doms"].get(r, [])
    if variant == "ga":
        doms = [x for x in doms if x[1] >= d["ga"][x[0]]]
    return [call_nt(a, b, minus, x[3], x[4]) for x in doms]


def prep(d, out):
    ms = misses(d)
    report_first(d, ms)
    with open(f"{out}/models.txt", "w") as f:
        f.write("\n".join(sorted({d["name"][m] for _, m, _, _ in ms})) + "\n")
    for v in VARIANTS:
        with open(f"{out}/{v}.pieces.fasta", "w") as f:
            for r in d["calls"]:
                cut = cuts(d, v, r)
                if cut:
                    for s, e in remainder(len(d["reads"][r]), cut):
                        if e - s >= 30:
                            f.write(f">{r}p{s}\n{d['reads'][r][s:e]}\n")
        with open(f"{out}/{v}.oracle.faa", "w") as f:
            for i, (r, m, _, hit) in enumerate(ms):
                cut = cuts(d, v, r)
                if not cut:
                    continue
                _, p, a, e, st = d["orfs"][hit[2]]
                run = best = ""
                for i_aa, c in enumerate(p):
                    x, y = (e - 3 * (i_aa + 1), e - 3 * i_aa) if st else (a + 3 * i_aa, a + 3 * (i_aa + 1))
                    run = "" if any(overlap(x, y, *k) for k in cut) else run + c
                    best = max(best, run, key=len)
                if len(best) >= 5:
                    f.write(f">M{i}\n{best}\n")


def score(d, out):
    ms = misses(d)
    print(f"{len(ms)} misses")
    for v in VARIANTS:
        new_frames = defaultdict(set)
        ncalls = res = 0
        for name, seq in fasta(f"{out}/{v}.calls.faa"):
            pr, a, b, st = call_name(name)
            r, off = pr.split("p")
            new_frames[pr] |= frames_of(a, b, st, int(off))
            ncalls += 1
            res += len(seq)
        new = {}
        with open(f"{out}/{v}.calls.domtbl") as f:
            for line in f:
                if line[0] != "#":
                    c = line.split()
                    k = (c[0].rsplit("_", 3)[0].split("p")[0], c[4].split(".")[0])
                    new[k] = max(new.get(k, -1e9), float(c[13]))
        oracle = {}
        with open(f"{out}/{v}.oracle.domtbl") as f:
            for line in f:
                if line[0] != "#":
                    c = line.split()
                    k = (c[0], c[4].split(".")[0])
                    oracle[k] = max(oracle.get(k, -1e9), float(c[13]))
        pieces = sum(1 for _ in fasta(f"{out}/{v}.pieces.fasta"))
        print(f"\n== removing {v}: {pieces} pieces >= 30 nt, {ncalls} new calls, {res / 50000:.1f} new residues per read")
        rows = Counter()
        best_case = Counter()
        gained = sum(1 for k, s in new.items() if s >= d["ga"].get(k[1], 1e9) and k not in d["truth"])
        for i, (r, m, b, hit) in enumerate(ms):
            ga = d["ga"][m]
            cut = cuts(d, v, r)
            o = oracle.get((f"M{i}", m))
            best_case[b] += o is not None and (b != 0 or o >= ga)
            n = new.get((r, m))
            if n is not None and (b != 0 or n >= ga):
                why = "RECOVERED"
            elif r not in d["calls"]:
                why = "read not called: same input"
            elif not cut:
                why = "nothing removed: same call"
            else:
                s, e = hit_nt(d, hit)
                left = remainder(len(d["reads"][r]), cut)
                inside = max(left, key=lambda p: overlap(s, e, *p), default=None)
                kept = overlap(s, e, *inside) if inside else 0
                _, _, a, _, st = d["orfs"][hit[2]]
                if kept == 0:
                    why = "hit region removed entirely"
                elif inside[1] - inside[0] <= FGS_MIN:
                    why = f"piece <= {FGS_MIN} nt: FGS cannot report it"
                elif f"{r}p{inside[0]}" not in new_frames:
                    why = "piece long enough: no call"
                elif (st, a % 3) not in new_frames[f"{r}p{inside[0]}"]:
                    why = "piece long enough: other frame called"
                else:
                    why = "right frame called: hit too weak"
            rows[why, b] += 1
        print(f"{'':<44}" + "".join(f"{n:>9}" for n in NAMES))
        for c in sorted({k[0] for k in rows}, key=lambda c: -sum(rows[c, i] for i in range(5))):
            print(f"{c:<44}" + "".join(f"{rows[c, i]:>9}" for i in range(5)))
        print(f"{'best case: right frame, no length limit':<44}" + "".join(f"{best_case[i]:>9}" for i in range(5)))
        print("new GA hits not in truth:", gained)


if __name__ == "__main__":
    {"prep": prep, "score": score}[sys.argv[1]](load(*sys.argv[2:9]), sys.argv[9])
