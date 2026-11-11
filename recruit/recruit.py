"""Recruit reads with weak HMMER hits to targets the sample confirms, keeping ambiguity as sets.

Inputs: HMMER pairs (recruit/pairs.py output, or an hmmsearch --domtblout) for the sample's ORFs and for
the reversed ORFs (`recruit.py reverse`), optionally a fraction of the reads. Steps (docs/recruit-plan.md):

1. Confirmed targets: models with >= MIN_GA reads at GA.
2. Per-pair likelihood ratio, true against null, from a logistic fit of real against reversed pairs to
   confirmed targets on FEATURES: E-value; alignment position (cut / inner / whole); score against the
   score a member is expected to reach over the same columns; ORF length.
3. Per locus (a read's pairs whose alignments overlap on the read, in any frame; separate stretches of
   one read, as two genes or two domains, are separate loci), P(true) from
   the best pair, then P(t | true) over the locus's targets as prior(t) * LR, prior(t) = GA reads of t
   plus its recruited mass (EM_ITERS rounds). Reads are kept at the per-read FDR target, estimated by
   running the reversed reads through the same model.
4. Unconfirmed targets whose summed log LR beats the reversed reads' at TARGET_FDR join the confirmed
   set (needs the whole null, fraction 1), and step 3 runs again.
5. Status: 'ga' (a pair at GA), 'one' (top P(t | true) >= TAU), 'set' (the smallest set summing to TAU),
   'none' (below the FDR cut). Every locus keeps its set; `top` is the single best target.
6. Depth check: recruited mass per target against rho * GA reads, negative binomial.

Usage:
  recruit.py run --real REAL --null NULL --models MODELS.tsv --prep PREP --orfs ORFS.faa --out DIR
                 [--null-fraction 1] [--fdr 0.02] [--tau 0.8] [--target-fdr 0.05] [--min-ga 1]
                 [--features e,pos,ratio,len] [--em-iters 1] [--orf-coords ORFS.coords.tsv]
  recruit.py reverse ORFS.faa OUT.faa [FRACTION]   reversed ORFs, of a fraction of the reads
Outputs in DIR: loci.tsv, hits.tsv, targets.tsv, candidates.faa, summary.json (docs/recruit.md).
"""

import argparse
import csv
import json
import math
import os
import sys
import zlib
from collections import defaultdict
from dataclasses import dataclass

import numpy as np

Z = 20_000_000
# Pfam/HMMER background, as measure/prep.py
Q = np.array([.0787945, .0151600, .0535222, .0668298, .0397062, .0695071, .0229198, .0590092,
              .0594422, .0963728, .0237718, .0414386, .0482904, .0395639, .0540978, .0683364,
              .0540687, .0673417, .0114135, .0304133])
Q /= Q.sum()


@dataclass
class Pair:
    orf: str
    read: str
    acc: str
    L: int
    score: float
    evalue: float
    hf: int
    ht: int
    af: int
    at: int
    pos: str = ""
    ratio: float = 0.0
    lr: float = 1.0
    rf: int = -1  # alignment on the forward read, nucleotides [rf, rt); -1 without ORF coordinates
    rt: int = -1
    extra: tuple = ()  # env_from, env_to, c_evalue, i_evalue, cigar, when pairs.py wrote them


def read_of(orf):
    return orf.rsplit("_", 1)[0]


def read_pairs(path, models):
    """Strongest domain per (ORF, model), from pairs.py TSV or an hmmsearch domtbl."""
    best = {}
    with open(path) as f:
        first = f.readline()
        if first.startswith("orf\t"):
            rows = ((c[0], c[1], int(c[3]), float(c[4]), float(c[5]), *map(int, c[6:10]), tuple(c[10:15]))
                    for c in (line.rstrip("\n").split("\t") for line in f))
        else:
            f.seek(0)
            rows = ((c[0], c[4].split(".")[0], int(c[2]), float(c[13]), float(c[6]), *map(int, c[15:19]),
                     (c[19], c[20], c[11], c[12], ""))
                    for c in (line.split() for line in f if line[0] != "#"))
        for orf, acc, L, score, e, hf, ht, af, at, extra in rows:
            if acc in models and ((orf, acc) not in best or score > best[orf, acc].score):
                best[orf, acc] = Pair(orf, read_of(orf), acc, L, score, e, hf, ht, af, at, extra=extra)
    return list(best.values())


def expected_bits(prep, accs):
    """{acc: cumulative expected bits per column}: a member's expected log-odds at each match state."""
    lo = np.fromfile(f"{prep}/profile.f32", dtype=np.float32).reshape(-1, 21)[:, :20].astype(np.float64)
    offs = np.fromfile(f"{prep}/offs.u32", dtype=np.uint32)
    names = open(f"{prep}/accs.txt").read().split()
    e = Q * np.exp2(lo)
    e /= e.sum(1, keepdims=True)
    col = (e * lo).sum(1)
    want = set(accs)
    return {a: np.concatenate(([0.0], np.cumsum(col[offs[i]:offs[i + 1]])))
            for i, a in enumerate(names) if a in want}


def annotate(pairs, models, cum):
    for p in pairs:
        M = models[p.acc]["M"]
        if p.hf <= 3 and p.ht >= M - 2:
            p.pos = "whole"
        elif (p.af <= 3 and p.hf > 3) or (p.at >= p.L - 2 and p.ht < M - 2):
            p.pos = "cut"
        else:
            p.pos = "inner"
        c = cum[p.acc]
        p.ratio = p.score / max(c[p.ht] - c[p.hf - 1], 1.0)


def design(ps, features):
    """Columns for the real-against-null logistic fit: log10 E with hinges at E = 1 and 100, position
    (cut, inner against whole) and its slope in E, score / expected score, log ORF length."""
    x = np.array([math.log10(min(max(p.evalue, 1e-30), 1e4)) for p in ps])
    cols = [np.ones(len(ps))]
    if "e" in features:
        cols += [x, np.maximum(x, 0), np.maximum(x - 2, 0)]
    if "pos" in features:
        for c in ("cut", "inner"):
            d = np.array([p.pos == c for p in ps], dtype=float)
            cols += [d] + ([d * x] if "e" in features else [])
    if "ratio" in features:
        r = np.array([min(p.ratio, 1.5) for p in ps])
        cols += [r, r * r] + ([r * x] if "e" in features else [])
    if "len" in features:
        cols.append(np.log([max(p.L, 1) for p in ps]))
    return np.column_stack(cols)


def logistic(X, y, w, ridge=1e-2, iters=50):
    """Weighted logistic regression by IRLS on standardised columns."""
    mu, sd = X.mean(0), X.std(0)
    mu[0], sd[0] = 0.0, 1.0
    sd[sd == 0] = 1.0
    Xs = (X - mu) / sd
    beta = np.zeros(X.shape[1])
    R = ridge * np.eye(X.shape[1])
    R[0, 0] = 0.0
    for _ in range(iters):
        eta = np.clip(Xs @ beta, -30, 30)
        pr = 1 / (1 + np.exp(-eta))
        W = w * pr * (1 - pr)
        g = Xs.T @ (w * (y - pr)) - R @ beta
        H = Xs.T @ (Xs * W[:, None]) + R
        step = np.linalg.solve(H, g)
        beta += step
        if np.abs(step).max() < 1e-8:
            break
    return beta, mu, sd


def fit_lr(real, null, w_null, features, apply_to):
    """Fit real against null pairs to confirmed targets; set p.lr = f1 / f0 on every pair in apply_to.
    With o the fitted odds of real, h / f0 = o * Nn / Nr and h = pi f0 + (1 - pi) f1 with pi = Nn / Nr,
    so LR = pi (o - 1) / (1 - pi)."""
    ps = real + null
    X = design(ps, features)
    y = np.array([1.0] * len(real) + [0.0] * len(null))
    w = np.array([1.0] * len(real) + [w_null] * len(null))
    beta, mu, sd = logistic(X, y, w)
    pi = min(w_null * len(null) / max(len(real), 1), 0.999)
    eta = np.clip(((design(apply_to, features) - mu) / sd) @ beta, -30, 30)
    lr = np.clip(pi * (np.exp(eta) - 1) / (1 - pi), 1e-6, 1e9)
    for p, v in zip(apply_to, lr):
        p.lr = float(v)
    return {"beta": beta.round(4).tolist(), "pairs_real": len(real), "pairs_null": len(null), "pi_pairs": pi}


def place(pairs, coords_path, reversed_orfs=False):
    """Alignment positions on the forward read, from ORF coordinates (measure/real_orfs.py coords.tsv:
    orf, strand, first base on its strand, read length). Reversed ORFs (the null) keep their ORF's
    footprint, with positions mirrored."""
    want = {p.orf for p in pairs}
    at = {}
    with open(coords_path) as f:
        next(f)
        for line in f:
            orf, sign, start, n = line.split("\t")
            if orf in want:
                at[orf] = (sign, int(start), int(n))
    for p in pairs:
        if p.orf not in at:
            continue
        sign, start, n = at[p.orf]
        af, a_t = (p.L + 1 - p.at, p.L + 1 - p.af) if reversed_orfs else (p.af, p.at)
        s0, s1 = start + 3 * (af - 1), start + 3 * a_t
        p.rf, p.rt = (s0, s1) if sign == "+" else (n - s1, n - s0)


def overlaps(p, q):
    """Alignments that share at least half of the shorter one: on the read when placed, else on one ORF."""
    if p.rf >= 0 and q.rf >= 0:
        return min(p.rt, q.rt) - max(p.rf, q.rf) >= 0.5 * min(p.rt - p.rf, q.rt - q.rf)
    return p.orf == q.orf and min(p.at, q.at) - max(p.af, q.af) + 1 >= 0.5 * min(p.at - p.af, q.at - q.af) + 1


def loci(pairs):
    """{(read, n): [pairs]}: a read's pairs grouped where their alignments overlap, strongest first.
    Placed on the read (ORF coordinates), pairs in any frame over the same bases compete and separate
    stretches (two genes, two domains) are separate loci. Without coordinates, alignments that don't
    overlap on the best ORF split, and the read's other ORFs join its best locus."""
    by_read = defaultdict(list)
    for p in pairs:
        by_read[p.read].append(p)
    out = {}
    for read, ps in by_read.items():
        ps = sorted(ps, key=lambda p: -p.score)
        placed = all(p.rf >= 0 for p in ps)
        groups = []
        for p in (ps if placed else [p for p in ps if p.orf == ps[0].orf]):
            for g in groups:
                if any(overlaps(p, q) for q in g):
                    g.append(p)
                    break
            else:
                groups.append([p])
        if not placed:
            groups[0] += [p for p in ps if p.orf != ps[0].orf]
        for n, g in enumerate(groups):
            out[read, n] = g
    return out


def score_loci(lc, pi0, prior):
    """{locus: (P(true), {acc: P(t | true)}, best pair per acc)}."""
    out = {}
    for k, ps in lc.items():
        best = {}
        for p in ps:
            if p.acc not in best or p.score > best[p.acc].score:
                best[p.acc] = p
        lrs = {a: p.lr for a, p in best.items()}
        top = max(lrs.values())
        p_true = (1 - pi0) * top / (pi0 + (1 - pi0) * top)
        w = {a: prior.get(a, 0.0) * v for a, v in lrs.items()}
        tot = sum(w.values())
        if tot <= 0:
            w, tot = lrs, sum(lrs.values())
        out[k] = (p_true, {a: v / tot for a, v in w.items()}, best)
    return out


def choose(ws, tau):
    order = sorted(ws, key=ws.get, reverse=True)
    acc, keep = 0.0, []
    for a in order:
        acc += ws[a]
        keep.append(a)
        if acc >= tau - 1e-12:
            break
    return ("one" if len(keep) == 1 else "set"), keep


def fdr_cut(real_p, null_p, w_null, target):
    """Lowest P(true) cut whose estimated FDR (scaled null reads kept / real reads kept) is <= target,
    and per-read q-values. Reads are scored by their best locus."""
    pts = sorted([(p, 0) for p in real_p] + [(p, 1) for p in null_p], reverse=True)
    q, kept_r, kept_d, cut = {}, 0, 0.0, 1.0 + 1e-9
    fdrs = []
    i = 0
    while i < len(pts):
        j = i
        while j < len(pts) and pts[j][0] == pts[i][0]:
            kept_r += pts[j][1] == 0
            kept_d += w_null * (pts[j][1] == 1)
            j += 1
        fdrs.append((pts[i][0], kept_d / max(kept_r, 1)))
        i = j
    best = math.inf
    for p, f in reversed(fdrs):  # q-value: the lowest FDR at this cut or any looser one
        best = min(best, f)
        q[p] = best
        if best <= target:
            cut = min(cut, p)
    return cut, q


def run_model(real, null, w_null, conf, ga_reads, models, features, em_iters):
    """Fit and score with confirmed set conf. Returns the fit, pi0, real loci, real and null locus scores."""
    rc = [p for p in real if p.acc in conf]
    nc = [p for p in null if p.acc in conf]
    lr = fit_lr(rc, nc, w_null, features, real + null)
    real_loci, null_loci = loci(rc), loci(nc)
    n_real = len({k[0] for k in real_loci})
    pi0 = min(w_null * len({k[0] for k in null_loci}) / max(n_real, 1), 0.999)
    prior = {a: float(ga_reads.get(a, 0)) for a in conf}
    for _ in range(em_iters + 1):
        S = score_loci(real_loci, pi0, {a: v + 1e-3 for a, v in prior.items()})
        mass = defaultdict(float)
        for k, (pt, ws, _) in S.items():
            if not any(p.score >= models[p.acc]["ga"] for p in real_loci[k]):
                for a, w in ws.items():
                    mass[a] += pt * w
        prior = {a: ga_reads.get(a, 0) + mass[a] for a in conf}
    N = score_loci(null_loci, pi0, prior)
    return lr, pi0, real_loci, S, N


def target_evidence(real, null, conf, bin_cols=15):
    """Summed log LR per unconfirmed target, counted once per region of the HMM: the best read in each
    window of BIN_COLS columns (by the alignment's midpoint). Overlapping reads of one deep region, as
    repeats and biased composition give, then add no more than one read; a present target's reads
    spread along the model."""
    out = []
    for ps in (real, null):
        best = {}
        for p in ps:
            if p.acc not in conf:
                k = (p.acc, (p.hf + p.ht) // 2 // bin_cols)
                best[k] = max(best.get(k, -math.inf), math.log(max(p.lr, 1e-9)))
        s = defaultdict(float)
        for (a, _), v in best.items():
            s[a] += max(v, 0.0)
        out.append(s)
    return out


def depth_flags(conf, ga_reads, mass):
    """rho = recruited mass per GA read; NB dispersion by moments; p-values for too many / too few."""
    from scipy.stats import nbinom

    ts = [a for a in conf if ga_reads.get(a, 0) > 0]
    rho = sum(mass.get(a, 0) for a in ts) / max(sum(ga_reads[a] for a in ts), 1)
    mu = {a: rho * ga_reads[a] for a in ts}
    num = sum((mass.get(a, 0) - mu[a]) ** 2 - mu[a] for a in ts)
    inv_k = max(num / max(sum(m * m for m in mu.values()), 1e-9), 1e-6)
    out = {}
    for a in ts:
        m, x = mu[a], round(mass.get(a, 0))
        n = 1 / inv_k
        pr = n / (n + m) if m > 0 else 1.0
        hi, lo = nbinom.sf(x - 1, n, pr), nbinom.cdf(x, n, pr)
        out[a] = (m, hi, lo, "high" if hi < 0.01 else "low" if lo < 0.01 else "")
    return rho, 1 / inv_k, out


def run(args):
    features = args.features.split(",")
    rows = list(csv.DictReader(open(args.models), delimiter="\t"))
    models = {r["acc"]: {"ga": float(r["ga"]), "M": int(r["length"]), "name": r["name"], "clan": r["clan"]}
              for r in rows}
    real, null = read_pairs(args.real, models), read_pairs(args.null, models)
    cum = expected_bits(args.prep, {p.acc for p in real + null})
    annotate(real + null, models, cum)
    if args.orf_coords:
        place(real, args.orf_coords)
        place(null, args.orf_coords, reversed_orfs=True)
    w_null = 1 / args.null_fraction
    ga_reads = defaultdict(set)
    for p in real:
        if p.score >= models[p.acc]["ga"]:
            ga_reads[p.acc].add(p.read)
    ga_n = {a: len(v) for a, v in ga_reads.items()}
    conf = {a for a, n in ga_n.items() if n >= args.min_ga}
    summary = {"features": features, "confirmed_ga": len(conf), "null_fraction": args.null_fraction}

    lr, pi0, rl, S, N = run_model(real, null, w_null, conf, ga_n, models, features, args.em_iters)
    weak_conf = {}
    if args.null_fraction >= 1 and args.target_fdr > 0:
        sr, sn = target_evidence(real, null, conf)
        vals = sorted(set(v for v in sr.values() if v > 0), reverse=True)
        thr = math.inf
        for v in vals:
            f = sum(x >= v for x in sn.values()) / max(sum(x >= v for x in sr.values()), 1)
            if f <= args.target_fdr:
                thr = v
        weak_conf = {a: v for a, v in sr.items() if v >= thr}
        summary["target_evidence_cut"] = thr if thr < math.inf else None
        if weak_conf:
            conf = conf | set(weak_conf)
            lr, pi0, rl, S, N = run_model(real, null, w_null, conf, ga_n, models, features, args.em_iters)
    summary["confirmed_weak"] = len(weak_conf)
    summary["pi0"] = pi0

    def read_best(L):
        out = defaultdict(float)
        for k, (pt, _, _) in L.items():
            out[k[0]] = max(out[k[0]], pt)
        return out

    rb, nb = read_best(S), read_best(N)
    ga_read_set = {k[0] for k, ps in rl.items() if any(p.score >= models[p.acc]["ga"] for p in ps)}
    # the FDR is over recruited reads: reads with a GA call are kept regardless
    cut, q = fdr_cut([v for k, v in rb.items() if k not in ga_read_set], list(nb.values()), w_null, args.fdr)
    summary.update(fdr_target=args.fdr, p_true_cut=cut, tau=args.tau)

    mass = defaultdict(float)
    counts = defaultdict(lambda: defaultdict(float))
    out_rows, hit_rows = [], []
    for k in sorted(S):
        pt, ws, best = S[k]
        is_ga = any(p.score >= models[p.acc]["ga"] for p in rl[k])
        kind, keep = choose(ws, args.tau)
        status = "ga" if is_ga else (kind if pt >= cut else "none")
        qv = 0.0 if is_ga else q.get(rb[k[0]], 1.0)
        top = keep[0]
        tp = best[top]
        for a in keep:
            p = best[a]
            env = p.extra if len(p.extra) == 5 else ("", "", "", "", "")
            w = ws[a] * (1.0 if status == "ga" else pt)
            hit_rows.append((p.orf, k[0], k[1], status, int(a == top), f"{w:.4f}", f"{pt:.4f}", f"{qv:.4f}",
                             int(a[2:]), p.evalue, p.score, env[2], env[3], p.hf, p.ht, env[0], env[1],
                             p.af, p.at, env[4], int(p.score >= models[a]["ga"]), p.rf, p.rt, p.L))
            counts[a][status] += ws[a]
            if status in ("one", "set"):
                mass[a] += ws[a]
        out_rows.append((k[0], k[1], status, f"{pt:.4f}", f"{qv:.4f}", top, f"{ws[top]:.4f}",
                         ",".join(f"{a}:{ws[a]:.3f}" for a in keep),
                         ",".join(sorted({models[a]["clan"] or a for a in keep})),
                         tp.orf, tp.score, tp.evalue, tp.hf, tp.ht, tp.af, tp.at, tp.rf, tp.rt, tp.pos, f"{tp.ratio:.3f}"))
    os.makedirs(args.out, exist_ok=True)
    with open(f"{args.out}/loci.tsv", "w") as o:
        o.write("read\tlocus\tstatus\tp_true\tq\ttop\ttop_p\tset\tclans\torf\tscore\tevalue\thmm_from\thmm_to\t"
                "ali_from\tali_to\tread_from\tread_to\tpos\tratio\n")
        for row in out_rows:
            o.write("\t".join(map(str, row)) + "\n")
    # for alignment against MGnify: pfam-reannotate's hits columns (query = the ORF), one row per set member
    with open(f"{args.out}/hits.tsv", "w") as o:
        o.write("query\tread\tlocus\tstatus\ttop\tweight\tp_true\tq\tpfam_accession\tseq_evalue\tscore\t"
                "c_evalue\ti_evalue\thmm_from\thmm_to\tenv_from\tenv_to\tali_from\tali_to\tcigar\tabove_ga\t"
                "read_from\tread_to\torf_len\n")
        for row in hit_rows:
            o.write("\t".join(map(str, row)) + "\n")
    rho, k_nb, flags = depth_flags(conf, ga_n, mass)
    with open(f"{args.out}/targets.tsv", "w") as o:
        o.write("acc\tname\tclan\tconfirmed\tga_reads\tone\tset_mass\trecruited_mass\texpected\tp_high\tp_low\tflag\n")
        for a in sorted(conf, key=lambda a: -(ga_n.get(a, 0) + mass[a])):
            m, hi, lo, fl = flags.get(a, (0.0, 1.0, 1.0, "weak-only"))
            o.write(f"{a}\t{models[a]['name']}\t{models[a]['clan']}\t{'ga' if a in ga_n else 'weak'}\t{ga_n.get(a, 0)}\t"
                    f"{counts[a]['one']:.0f}\t{counts[a]['set']:.2f}\t{mass[a]:.2f}\t{m:.2f}\t{hi:.3g}\t{lo:.3g}\t{fl}\n")
    # candidates for alignment: every kept locus, each set member, with its weight
    need = {}
    for row in out_rows:
        if row[2] != "none":
            for item in row[7].split(","):
                a, w = item.split(":")
                p = S[row[0], row[1]][2][a]
                need.setdefault(p.orf, []).append((a, row[2], float(w) * (1.0 if row[2] == "ga" else float(row[3])), a == row[5]))
    seqs, name = {}, None
    for line in open(args.orfs):
        if line[0] == ">":
            name = line[1:].split()[0]
        elif name in need:
            seqs[name] = seqs.get(name, "") + line.strip()
    by_target = defaultdict(list)
    for orf, items in need.items():
        for a, st, w, top in items:
            by_target[a].append((orf, st, w, top))
    with open(f"{args.out}/candidates.faa", "w") as o:
        for a in sorted(by_target):
            for orf, st, w, top in sorted(by_target[a]):
                o.write(f">{orf} target={a} status={st} weight={w:.3f} top={int(top)}\n{seqs.get(orf, '')}\n")
    st = defaultdict(int)
    for row in out_rows:
        st[row[2]] += 1
    nst = defaultdict(int)
    for k, (pt, ws, _) in N.items():
        nst["kept" if pt >= cut else "none"] += 1
    summary.update(loci=dict(st), null_loci=dict(nst),
                   reads_ga=len(ga_read_set),
                   reads_recruited=sum(1 for k, v in rb.items() if k not in ga_read_set and v >= cut),
                   null_reads_kept_scaled=w_null * sum(v >= cut for v in nb.values()),
                   rho=rho, nb_size=k_nb, depth_flags={f: sum(x[3] == f for x in flags.values()) for f in ("high", "low")},
                   lr=lr)
    json.dump(summary, open(f"{args.out}/summary.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in summary.items() if k != "lr"}))


def reverse(path, out, fraction=1.0):
    """Reversed ORFs (the null), of the reads whose name hashes under FRACTION."""
    keep, name = False, None
    with open(out, "w") as o:
        for line in open(path):
            if line[0] == ">":
                name = line[1:].split()[0]
                keep = zlib.crc32(read_of(name).encode()) / 2**32 < fraction
                if keep:
                    o.write(line)
            elif keep:
                o.write(line.strip()[::-1] + "\n")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "reverse":
        reverse(sys.argv[2], sys.argv[3], float(sys.argv[4]) if len(sys.argv) > 4 else 1.0)
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run"])
    for a in ("--real", "--null", "--models", "--prep", "--orfs", "--out"):
        ap.add_argument(a, required=True)
    ap.add_argument("--null-fraction", type=float, default=1.0)
    ap.add_argument("--fdr", type=float, default=0.02)
    ap.add_argument("--tau", type=float, default=0.8)
    ap.add_argument("--target-fdr", type=float, default=0.05)
    ap.add_argument("--min-ga", type=int, default=1)
    ap.add_argument("--features", default="e,pos,ratio,len")
    ap.add_argument("--em-iters", type=int, default=1)
    ap.add_argument("--orf-coords", help="measure/real_orfs.py coords.tsv: place alignments on the read")
    run(ap.parse_args())


if __name__ == "__main__":
    main()
