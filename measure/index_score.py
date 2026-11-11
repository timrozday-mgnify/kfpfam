"""Recall and candidates of a real index search against hmmsearch, on the same reads.

Usage: index_score.py SEARCH.tsv DECOY_SEARCH.tsv HITS.domtbl DECOY_HITS.domtbl MODELS.tsv NREADS
"""

import csv
import sys
from collections import defaultdict

LOWS = (30, 10, 0, -5, -10, -20, float("-inf"))
NAMES = ("GA+30", "+10-30", "GA-+10", "-5-GA", "-10--5", "-20--10", "<-20")


def hits(path, index, ga):
    """Per (read, model), the strongest domain's bin."""
    best = {}
    with open(path) as f:
        for line in f:
            if line[0] != "#":
                c = line.split()
                k = (int(c[0][1:].split("_")[0]), index[c[4].split(".")[0]])
                best[k] = max(best.get(k, float("-inf")), float(c[13]))
    return {k: next(i for i, low in enumerate(LOWS) if s - ga[k[1]] >= low) for k, s in best.items()}


def search(path):
    """Per (read, model), the best diagonal and gapped scores."""
    best = {}
    with open(path) as f:
        for line in f:
            q, m, d, g = line.split()
            k = (int(q[1:].split("_")[0]), int(m))
            old = best.get(k, (0.0, 0.0))
            best[k] = (max(old[0], float(d)), max(old[1], float(g)))
    return best


def main():
    found, fake, table, fake_table, models, nreads = sys.argv[1:7]
    nreads = int(nreads)
    with open(models) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    index = {r["acc"]: i for i, r in enumerate(rows)}
    ga = [float(r["ga"]) for r in rows]
    truth, truth_fake = hits(table, index, ga), hits(fake_table, index, ga)
    found, fake = search(found), search(fake)
    total = defaultdict(int)
    for b in truth.values():
        total[b] += 1
    print("hits per bin", {NAMES[i]: total[i] for i in range(7)})

    def show(label, metric, t):
        ok = [0] * 7
        for k, v in found.items():
            if v[metric] >= t and k in truth:
                ok[truth[k]] += 1
        cand = sum(v[metric] >= t for v in found.values()) / nreads
        decoy = sum(v[metric] >= t for v in fake.values()) / nreads
        above = sum(ok[:3]) / max(sum(total[i] for i in range(3)), 1)
        by = " ".join(f"{ok[i] / max(total[i], 1):.3f}" for i in range(7))
        print(f"{label:8} >= {t:5.2f}  cand {cand:7.3f}  decoy {decoy:7.3f}  >=GA {above:.3f} | {by}")

    print(" " * 50, " ".join(NAMES))
    for t in (22, 24, 26, 28):
        show("diagonal", 0, t)
    import os
    gapped_bars = [float(x) for x in os.environ.get("GAPPED_BARS", "22.25 24.25 26.25 28.25").split()]
    for t in gapped_bars:
        show("gapped", 1, t)
    assert all(v[1] >= v[0] - 0.02 for v in found.values()), "gaps can only raise a score"


if __name__ == "__main__":
    main()
