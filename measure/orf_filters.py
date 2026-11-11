"""Filters on six-frame ORFs and gene calling, against hmmsearch on all ORFs of the same reads.

orf_filters.py reads SUBSET R1 R2 ORFS.faa OUT.fasta OUT_COORDS.tsv
    Rebuilds real_orfs.py's subset reads and each ORF's place on its read; checks the ORFs match.
orf_filters.py score COORDS.tsv ORFS.domtbl DECOYS.domtbl MODELS.tsv NREADS [NAME CALLS.faa CALLS.domtbl [CALLS_DECOY.domtbl]]...
    Recall per match-strength bin for each rule, and for each gene caller's searched output.
orf_filters.py frames READS.fasta CALLS.faa ORFS.domtbl DECOYS.domtbl CALLS.domtbl CALLS_DECOY.domtbl MODELS.tsv NREADS
    The same, choosing ORFs by coding score. Hexamer log-odds, ORFs in the called frame against the other
    frames' ORFs, trained on half the reads and applied to the other half. No new search: ORF hits are reused.
    Plain ORFs stand in for FGS's proteins, which repeat a codon at a read end on most calls.
"""

import csv
import gzip
import random
import sys
from collections import Counter, defaultdict
from math import log

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from real_orfs import CODON, RC  # noqa: E402

LOWS = (0, -5, -10, -20, float("-inf"))
NAMES = ("≥GA", "-5-GA", "-10--5", "-20--10", "<-20")


def placed(seq, min_len=15):
    """ORFs in real_orfs.py's order, with start and end on the read's forward strand."""
    out = []
    n = len(seq)
    for strand, s in enumerate((seq, seq[::-1].translate(RC))):
        for frame in range(3):
            prot = "".join(CODON.get(s[i:i + 3], "X") for i in range(frame, len(s) - 2, 3))
            pos = 0
            for p in prot.split("*"):
                if len(p) >= min_len:
                    a, b = frame + pos * 3, frame + (pos + len(p)) * 3
                    out.append((p, a, b, 0) if strand == 0 else (p, n - b, n - a, 1))
                pos += len(p) + 1
    return out


def reads(subset, r1, r2, orfs_path, out, coords):
    seqs = []
    for path in (r1, r2):
        with gzip.open(path, "rt") as f:
            seqs += [line.strip() for i, line in enumerate(f) if i % 4 == 1]
    random.Random(5).shuffle(seqs)
    with open(orfs_path) as f:
        known = dict(zip((x[1:].strip() for x in f), f))
    n = 0
    with open(out, "w") as fa, open(coords, "w") as c:
        for r, s in enumerate(seqs[:int(subset)]):
            fa.write(f">R{r}\n{s}\n")
            for j, (p, a, b, _) in enumerate(placed(s)):
                assert known[f"R{r}_{j}"].strip() == p, (r, j)
                c.write(f"R{r}\t{j}\t{a}\t{b}\n")
                n += 1
    assert n == len(known), (n, len(known))
    print("matched", n, "ORFs")


def domtbl(path):
    """(query, model) -> best domain score."""
    best = {}
    with open(path) as f:
        for line in f:
            if line[0] != "#":
                c = line.split()
                k = (c[0], c[4].split(".")[0])
                best[k] = max(best.get(k, -1e9), float(c[13]))
    return best


def bin_of(d):
    return next(i for i, low in enumerate(LOWS) if d >= low)


def per_read(hits):
    out = defaultdict(float)
    for (q, m), s in hits.items():
        k = (q.split("_")[0], m)
        out[k] = max(out[k], s) if k in out else s
    return out


def longest_per_read(by_read, keep):
    for r, rows in by_read.items():
        yield from sorted(rows, key=lambda x: x[3] - x[2], reverse=True)[:keep]


def non_overlapping(by_read, allowed):
    for rows in by_read.values():
        kept = []
        for row in sorted(rows, key=lambda x: x[3] - x[2], reverse=True):
            if all(min(row[3], k[3]) - max(row[2], k[2]) <= allowed for k in kept):
                kept.append(row)
        yield from kept


def report(name, seqs, residues, found, decoy_n, truth, ga, nreads):
    tot, ok = [0] * 5, [0] * 5
    for (r, m), s in truth.items():
        d = s - ga[m]
        b = bin_of(d)
        tot[b] += 1
        got = found.get((r, m))
        ok[b] += got is not None and (b != 0 or got >= ga[m])
    gained = sum(1 for (r, m), s in found.items() if s >= ga[m] and truth.get((r, m), -1e9) < ga[m])
    print(f"{name:<26} {seqs / nreads:5.2f} {residues / nreads:7.1f} " + " ".join(f"{o / t:.3f}" for o, t in zip(ok, tot))
          + f"  {len(found) / nreads:6.3f} "
          + ("     -" if decoy_n is None else f"{decoy_n / nreads:6.3f}") + f" {gained:5d}")


def score(coords, table, decoys, models, nreads, *calls):
    nreads = int(nreads)
    with open(models) as f:
        ga = {r["acc"]: float(r["ga"]) for r in csv.DictReader(f, delimiter="\t")}
    by_read = defaultdict(list)
    with open(coords) as f:
        for line in f:
            r, j, a, b = line.split()
            by_read[r].append((r, f"{r}_{j}", int(a), int(b)))
    hits = domtbl(table)
    fake = domtbl(decoys)
    truth = per_read(hits)
    print("truth hits per bin", dict(zip(NAMES, (sum(bin_of(s - ga[m]) == i for (_, m), s in truth.items()) for i in range(5)))))
    print(f"{'rule':<26} seq/rd res/rd  " + " ".join(f"{n:>5}" for n in NAMES) + "  hits/rd rev/rd newGA")
    every = [row for rows in by_read.values() for row in rows]
    rules = [("all ORFs", every)]
    for m in (20, 25, 30, 35):
        rules.append((f"length >= {m} aa", [x for x in every if x[3] - x[2] >= 3 * m]))
    for k in (1, 2, 3):
        rules.append((f"{k} longest per read", list(longest_per_read(by_read, k))))
    for allowed in (0, 15, 30, 60):
        rules.append((f"longest, overlap <= {allowed} nt", list(non_overlapping(by_read, allowed))))
    for name, rows in rules:
        kept = {x[1] for x in rows}
        found = per_read({k: s for k, s in hits.items() if k[0] in kept})
        n_fake = sum(1 for k in per_read({k: s for k, s in fake.items() if k[0] in kept}))
        report(name, len(rows), sum(x[3] - x[2] for x in rows) / 3, found, n_fake, truth, ga, nreads)
    i = 0
    while i < len(calls):
        name, faa, tbl, dec = calls[i:i + 4]
        i += 4
        seqs = residues = 0
        with open(faa) as f:
            for line in f:
                if line[0] == ">":
                    seqs += 1
                else:
                    residues += len(line.strip())
        norm = {(q.rsplit("_", 3)[0] + "_0", m): s for (q, m), s in domtbl(tbl).items()}
        found = per_read(norm)
        n_fake = None if dec == "-" else len(per_read({(q.rsplit("_", 3)[0] + "_0", m): s for (q, m), s in domtbl(dec).items()}))
        report(name, seqs, residues, found, n_fake, truth, ga, nreads)


def fasta(path):
    with open(path) as f:
        name = None
        for line in f:
            if line[0] == ">":
                name = line[1:].split()[0]
            else:
                yield name, line.strip()


def hexamers(s):
    return (s[i:i + 6] for i in range(0, len(s) - 5, 3))


def log_odds(coding, other):
    """Log-odds of each hexamer, with one pseudocount per hexamer."""
    nc, no = sum(coding.values()) + 4096, sum(other.values()) + 4096
    return {h: log((coding[h] + 1) / nc) - log((other[h] + 1) / no) for h in set(coding) | set(other)}


def by_family(found, truth, ga):
    tot, ok = Counter(), Counter()
    for (r, m), s in truth.items():
        if s >= ga[m]:
            tot[m] += 1
            ok[m] += found.get((r, m), -1e9) >= ga[m]
    rec = sorted(ok[m] / t for m, t in tot.items() if t >= 5)
    return f"families {len(rec)}: mean {sum(rec) / len(rec):.3f} p10 {rec[len(rec) // 10]:.2f} min {rec[0]:.2f} below .8 {sum(x < .8 for x in rec)}"


def frames(reads_fa, calls_faa, table, decoys, calls_tbl, calls_dec, models, nreads):
    nreads = int(nreads)
    with open(models) as f:
        ga = {r["acc"]: float(r["ga"]) for r in csv.DictReader(f, delimiter="\t")}
    orfs = {}  # read -> [(id, a, b, strand, nucleotides)]
    for r, s in fasta(reads_fa):
        rows = orfs[r] = []
        for j, (_, a, b, strand) in enumerate(placed(s)):
            rows.append((f"{r}_{j}", a, b, strand, s[a:b] if not strand else s[a:b][::-1].translate(RC)))
    # A call's frame at both ends, since a frameshift moves it. Names are read_start_end_strand.
    call_frames, residues = defaultdict(set), 0
    for name, seq in fasta(calls_faa):
        r, a, b, strand = name.rsplit("_", 3)
        st = strand == "-"
        call_frames[r] |= {(st, (int(a) - 1) % 3), (st, int(b) % 3)}
        residues += len(seq)
    n_calls = sum(1 for _ in fasta(calls_faa))

    coding, other = [Counter(), Counter()], [Counter(), Counter()]
    fold = lambda r: int(r[1:]) % 2  # noqa: E731
    for r, called in call_frames.items():
        for _, a, _, st, nt in orfs[r]:
            ((coding if (bool(st), a % 3) in called else other)[fold(r)]).update(hexamers(nt))
    tables = [log_odds(coding[1 - i], other[1 - i]) for i in (0, 1)]
    scored = {}  # read -> rows with score, best first
    for r, rows in orfs.items():
        t = tables[fold(r)]
        scored[r] = sorted(((sum(t.get(h, 0.0) for h in hexamers(x[4])), x) for x in rows), key=lambda v: -v[0])
    top = [v[0][1] for r, v in scored.items() if v and r in call_frames]
    print(f"best-scoring ORF is in the called frame on {sum((bool(x[3]), x[1] % 3) in call_frames[x[0].split('_')[0]] for x in top) / len(top):.3f} of called reads")

    hits, fake, truth = domtbl(table), domtbl(decoys), None
    truth = per_read(hits)
    norm = lambda d: per_read({(q.rsplit("_", 3)[0] + "_0", m): s for (q, m), s in d.items()})  # noqa: E731
    call_hits, call_fake = norm(domtbl(calls_tbl)), norm(domtbl(calls_dec))
    print(f"{'rule':<26} seq/rd res/rd  " + " ".join(f"{n:>5}" for n in NAMES) + "  hits/rd rev/rd newGA")

    def run(name, rows, with_calls):
        kept = {x[0] for x in rows}
        found = per_read({k: s for k, s in hits.items() if k[0] in kept})
        rev = set(per_read({k: s for k, s in fake.items() if k[0] in kept}))
        seqs, res = len(rows), sum(x[2] - x[1] for x in rows) / 3
        if with_calls:
            for k, s in call_hits.items():
                found[k] = max(found[k], s) if k in found else s
            rev |= set(call_fake)
            seqs, res = seqs + n_calls, res + residues
        report(name, seqs, res, found, len(rev), truth, ga, nreads)
        print(" " * 27 + by_family(found, truth, ga))

    def outside(r):
        return [v for v in scored[r] if (bool(v[1][3]), v[1][1] % 3) not in call_frames.get(r, ())]

    run("FGS alone", [], True)
    for k in (1, 2, 3):
        run(f"{k} best by coding score", [x for v in scored.values() for _, x in v[:k]], False)
    run("FGS + best other", [o[0][1] for o in map(outside, scored) if o], True)
    run("FGS + longest other", [max((x for _, x in o), key=lambda x: x[2] - x[1]) for o in map(outside, scored) if o], True)
    run("FGS + best 2 others", [x for o in map(outside, scored) for _, x in o[:2]], True)
    for t in (20, 10, 5, 0, -5, -10):
        run(f"FGS + others score >= {t}", [x for r in scored for v, x in outside(r) if v >= t], True)

    print("\nplain ORFs in FGS's frame in place of its proteins")
    def inside(r):
        return [x for v, x in scored[r] if (bool(x[3]), x[1] % 3) in call_frames.get(r, ())]

    base = [x for r in scored for x in inside(r)]
    run("frame ORFs", base, False)
    run("frame + best other", base + [o[0][1] for o in map(outside, scored) if o], False)
    run("frame + longest other", base + [max((x for _, x in o), key=lambda x: x[2] - x[1]) for o in map(outside, scored) if o], False)
    run("frame + best 2 others", base + [x for o in map(outside, scored) for _, x in o[:2]], False)
    for t in (20, 15, 10, 5, 0, -5, -10):
        run(f"frame + others score >= {t}", base + [x for r in scored for v, x in outside(r) if v >= t], False)
    run("all ORFs", [x for v in scored.values() for _, x in v], False)


if __name__ == "__main__":
    {"reads": reads, "score": score, "frames": frames}[sys.argv[1]](*sys.argv[2:])
