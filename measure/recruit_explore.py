"""Exploration: weak hits on real reads against targets confirmed at GA, with reversed reads as the null.

Usage: recruit_explore.py REAL.domtbl DECOY.domtbl MODELS.tsv
Pair = (read, model), strongest domain. Confirmed target = model with >= 1 read at GA in the real reads.
"""

import csv
import math
import sys
from collections import Counter, defaultdict


def pairs(path):
    best = {}
    for line in open(path):
        if line[0] != "#":
            c = line.split()
            k = (int(c[0][1:].split("_")[0]), c[4].split(".")[0])
            v = (float(c[13]), float(c[6]), int(c[16]) - int(c[15]) + 1)
            if k not in best or v[0] > best[k][0]:
                best[k] = v
    return best


real, decoy, models = sys.argv[1:]
rows = {r["acc"]: r for r in csv.DictReader(open(models), delimiter="\t")}
ga = {a: float(r["ga"]) for a, r in rows.items()}
clan = {a: r["clan"] or a for a, r in rows.items()}
R, D = ({k: v for k, v in pairs(p).items() if k[1] in ga} for p in (real, decoy))
gareads = defaultdict(set)
for (q, a), v in R.items():
    if v[0] >= ga[a]:
        gareads[a].add(q)
conf = set(gareads)
called = {q for s in gareads.values() for q in s}
print(f"confirmed targets {len(conf)}, reads with a GA call {len(called)}")

# 1. Weak pairs to confirmed targets on reads without any GA call, real against reversed, by score
print("\nweak pairs (read has no GA call) to confirmed targets, real / reversed, by best-domain score - GA")
for lo, hi in ((-5, 0), (-10, -5), (-15, -10), (-20, -15), (-99, -20)):
    r = sum(1 for (q, a), v in R.items() if a in conf and q not in called and lo <= v[0] - ga[a] < hi)
    d = sum(1 for (q, a), v in D.items() if a in conf and lo <= v[0] - ga[a] < hi)
    print(f"  {lo:4}..{hi:3}: real {r:6}  reversed {d:6}  excess {r - d:6}  FDR {d / max(r, 1):.2f}")
print("same, by E (Z = 20 M)")
for e in (1e-3, 1e-2, 0.1, 1, 10, 100, 1000):
    r = sum(1 for (q, a), v in R.items() if a in conf and q not in called and v[1] <= e and v[0] < ga[a])
    d = sum(1 for (q, a), v in D.items() if a in conf and v[1] <= e and v[0] < ga[a])
    print(f"  E <= {e:6g}: real {r:6}  reversed {d:6}  FDR {d / max(r, 1):.2f}")
print("and to unconfirmed targets, E <= 1 / 10 / 100:",
      [(sum(1 for (q, a), v in R.items() if a not in conf and v[1] <= e),
        sum(1 for (q, a), v in D.items() if a not in conf and v[1] <= e)) for e in (1, 10, 100)])

# 2. Depth: excess weak reads per target against GA reads
print("\nper confirmed target, by GA reads: targets, weak real / reversed per target (E <= 10, no GA call)")
wr, wd = Counter(), Counter()
for (q, a), v in R.items():
    if a in conf and q not in called and v[1] <= 10:
        wr[a] += 1
for (q, a), v in D.items():
    if a in conf and v[1] <= 10:
        wd[a] += 1
for lo, hi in ((1, 2), (2, 4), (4, 11), (11, 31), (31, 10**9)):
    ts = [a for a in conf if lo <= len(gareads[a]) < hi]
    n = sum(len(gareads[a]) for a in ts)
    r, d = sum(wr[a] for a in ts), sum(wd[a] for a in ts)
    print(f"  {lo:3}-{hi - 1 if hi < 10**9 else '':4}: {len(ts):4} targets, {n:5} GA reads, weak {r / len(ts):6.2f} / {d / len(ts):5.2f}, "
          f"excess per GA read {(r - d) / n:.2f}")
print("excess per GA read by model length (E <= 10):")
for lo, hi in ((0, 50), (50, 100), (100, 200), (200, 400), (400, 10**5)):
    ts = [a for a in conf if lo <= int(rows[a]["length"]) < hi]
    n = sum(len(gareads[a]) for a in ts)
    print(f"  length {lo}-{hi}: {len(ts)} targets, excess per GA read {(sum(wr[a] for a in ts) - sum(wd[a] for a in ts)) / max(n, 1):.2f}")

# 3. Ambiguity: per weak read, how many confirmed targets it hits (E <= 10), and how many clans
for name, P, skip in (("real", R, called), ("reversed", D, set())):
    per = defaultdict(list)
    for (q, a), v in P.items():
        if a in conf and q not in skip and v[1] <= 10 and v[0] < ga[a]:
            per[q].append((v[0], a))
    n = Counter(min(len(x), 4) for x in per.values())
    cl = Counter(min(len({clan[a] for _, a in x}), 4) for x in per.values())
    print(f"\n{name}: {len(per)} reads with a weak pair to a confirmed target; targets per read {dict(sorted(n.items()))}; clans {dict(sorted(cl.items()))}")
    gaps = [sorted(x)[-1][0] - sorted(x)[-2][0] for x in per.values() if len(x) > 1]
    if gaps:
        gaps.sort()
        print(f"  score gap best - second: median {gaps[len(gaps) // 2]:.1f}, share < 2 bits {sum(g < 2 for g in gaps) / len(gaps):.2f}")

# 4. Do the GA reads themselves have competing weak targets (how ambiguous are confirmed hits)?
comp = []
for q in called:
    pass
