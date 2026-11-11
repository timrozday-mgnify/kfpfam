"""Recruit weak hits to confirmed targets: a mixture over reads, fitted by EM on the sample itself.

Each read with a weak pair to a confirmed target (one with a GA read) comes from the null or from one of
its targets. Step 1, P(true): pair likelihood ratio by E bin, null from reversed reads, true = real minus
reversed; pi0 = reversed weak reads / real weak reads. Step 2, P(t | true) over the read's targets by
prior * LR; prior 'depth' = GA reads of t, 'flat' = equal. A read kept at P(true) >= CUT goes to one target
when P(t | true) >= TAU, else to the smallest set with sum >= TAU (ambiguous). Checked on reversed reads (false recruitment) and on windows cut
from GA reads (known origin, recruit_trunc.py).

Usage: recruit_model.py REAL.domtbl DECOY.domtbl MODELS.tsv TRUNC.tsv
"""

import csv
import math
import sys
from collections import defaultdict

Z, EMAX = 20_000_000, 1000
PMAX = EMAX / Z


def pairs(path, ga):
    best = {}
    for line in open(path):
        if line[0] != "#":
            c = line.split()
            a = c[4].split(".")[0]
            if a in ga:
                k = (int(c[0][1:].split("_")[0]), a)
                best[k] = min(best.get(k, (9e9, 0)), (float(c[6]), -float(c[13])))
    return {k: (e, -s) for k, (e, s) in best.items()}


def by_read(P, conf, skip, ga):
    out = defaultdict(dict)
    for (q, a), (e, s) in P.items():
        if a in conf and q not in skip and s < ga[a]:
            out[q][a] = max(e / Z, 1e-300)
    return out


EDGES = (-3, -2, -1, 0, 1, 2, 3)  # log10 E bins: <= 1e-3, ..., 100-1000


def ebin(p):
    le = math.log10(p * Z)
    return next((i for i, x in enumerate(EDGES) if le <= x), len(EDGES) - 1)


def fit_lr(real_reads, rev_reads):
    """Empirical per-pair likelihood ratio by E bin: null from reversed reads, true = real - reversed."""
    r, d = [0] * len(EDGES), [0] * len(EDGES)
    for reads, h in ((real_reads, r), (rev_reads, d)):
        for hs in reads.values():
            for p in hs.values():
                h[ebin(p)] += 1
    ex = [max(a - b, 0.5) for a, b in zip(r, d)]
    f0 = [max(b, 0.5) / sum(d) for b in d]
    f1 = [x / sum(ex) for x in ex]
    pi0 = len(rev_reads) / len(real_reads)
    return [a / b for a, b in zip(f1, f0)], pi0, r, d


def post(reads, lr, pi0, prior):
    """Per read: P(from some confirmed target) from the best pair, then split over the read's targets by
    prior * LR. Returns {read: (p_true, {target: P(t | true)})}."""
    out = {}
    for q, hs in reads.items():
        best = max(lr[ebin(p)] for p in hs.values())
        p_true = (1 - pi0) * best / (pi0 + (1 - pi0) * best)
        w = {t: prior(t) * lr[ebin(p)] for t, p in hs.items()}
        tot = sum(w.values()) or 1
        out[q] = (p_true, {t: v / tot for t, v in w.items()})
    return out


def assign(ws, tau):
    """'one', 'set' or 'none', and the targets."""
    if not ws:
        return "none", []
    order = sorted(ws, key=ws.get, reverse=True)
    if ws[order[0]] >= tau:
        return "one", order[:1]
    acc, keep = 0.0, []
    for t in order:
        acc += ws[t]
        keep.append(t)
        if acc >= tau:
            return "set", keep
    return "none", []


def main():
    real, decoy, models, trunc = sys.argv[1:]
    rows = list(csv.DictReader(open(models), delimiter="\t"))
    ga = {r["acc"]: float(r["ga"]) for r in rows}
    R, D = pairs(real, ga), pairs(decoy, ga)
    g = defaultdict(int)
    called = set()
    for (q, a), (e, s) in R.items():
        if s >= ga[a]:
            g[a] += 1
            called.add(q)
    conf = set(g)
    real_reads, rev_reads = by_read(R, conf, called, ga), by_read(D, conf, set(), ga)
    # cut windows: name -> {model: p}, origin set
    win, origin = defaultdict(dict), {}
    for line in open(trunc):
        n, o, a, s, e = line.split()
        win[n][a] = max(float(e) / Z, 1e-300)
        origin[n] = set(o.split(",")) if o != "-" else set()
    lr, pi0, r, d = fit_lr(real_reads, rev_reads)
    print(f"confirmed {len(conf)}; weak reads: real {len(real_reads)}, reversed {len(rev_reads)}; pi0 {pi0:.3f}")
    print("pairs by E bin real / reversed / LR: " + "  ".join(f"<=1e{e}: {a}/{b}/{x:.2f}" for e, a, b, x in zip(EDGES, r, d, lr)))
    priors = {"depth": lambda t: g[t], "flat": lambda t: 1.0}
    for name, prior in priors.items():
        Pr, Pd, Pw = (post(x, lr, pi0, prior) for x in (real_reads, rev_reads, win))
        print(f"\nprior {name}; FDR cut on P(true), then TAU on P(t | true)")
        print(" P(true)>=  tau | real one/set  reversed  est.FDR | cut windows w: one right / set holds origin / wrong ; reversed windows kept")
        for cut in (0.5, 0.7, 0.9):
            for tau in (0.5, 0.8, 0.95):
                def res(P, n):
                    pt, ws = P[n]
                    return assign(ws, tau) if pt >= cut else ("none", [])
                cr = [res(Pr, q)[0] for q in Pr]
                cd = [res(Pd, q)[0] for q in Pd]
                kept = sum(k != "none" for k in cr)
                line = f"   {cut:4}  {tau:4} | {cr.count('one'):5} {cr.count('set'):4}  {sum(k != 'none' for k in cd):5}  {sum(k != 'none' for k in cd) / max(kept, 1):.3f} |"
                for w in (10, 15, 20, 25):
                    names = [n for n in origin if n.startswith(f"T{w}_")]
                    out = [res(Pw, n) if n in Pw else ("none", []) for n in names]
                    right = sum(k == "one" and t[0] in origin[n] for (k, t), n in zip(out, names))
                    inset = sum(k == "set" and bool(set(t) & origin[n]) for (k, t), n in zip(out, names))
                    wrong = sum(k != "none" and not set(t) & origin[n] for (k, t), n in zip(out, names))
                    dn = [n for n in origin if n.startswith(f"D{w}_")]
                    fake = sum(res(Pw, n)[0] != "none" for n in dn if n in Pw)
                    line += f" {w}: {right / len(names):.3f}/{inset / len(names):.3f}/{wrong / len(names):.3f} {fake / len(dn):.4f};"
                print(line)
    # reported (any pair) share of cut windows, the ceiling
    for w in (10, 15, 20, 25):
        names = [n for n in origin if n.startswith(f"T{w}_")]
        print(f"cut windows w={w}: {len(names)}, with a pair to the origin {sum(bool(origin[n] & set(win.get(n, {}))) for n in names) / len(names):.3f}, "
              f"with any pair {sum(n in win for n in names) / len(names):.3f}; with a pair to >1 target {sum(len(win.get(n, {})) > 1 for n in names) / len(names):.3f}")

if __name__ == "__main__":
    main()
