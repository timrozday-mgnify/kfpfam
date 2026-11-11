"""Recall of nail and the two modes by GA bin and by how much of the HMM the hit covers.

Truth: hmmsearch's strongest domain per (read, model); coverage = (hmm_to - hmm_from + 1) / model length.
Usage: cov_bins.py DOMTBL MODELS.tsv NAIL_SEEDS NAIL_TBL MODE_A.tsv BAR_A MODE_B.tsv BAR_B
"""

import csv
import sys

GA_LOWS = (0, -5, -10, -20, float("-inf"))
GA_NAMES = (">=GA", "-5-GA", "-10--5", "-20--10", "<-20")
COV_LOWS = (0.7, 0.4, 0.2, 0.1, 0)
COV_NAMES = (">=70%", "40-70%", "20-40%", "10-20%", "<10%")
DEPTH_LOWS = (101, 31, 11, 4, 2, 1)
DEPTH_NAMES = (">100", "31-100", "11-30", "4-10", "2-3", "1")


def read_key(q):
    return int(q[1:].split("_")[0])


def main():
    domtbl, models, seeds, tbl, a_tsv, a_bar, b_tsv, b_bar = sys.argv[1:]
    rows = list(csv.DictReader(open(models), delimiter="\t"))
    ga = {r["acc"]: float(r["ga"]) for r in rows}
    by_name = {r["name"]: r["acc"] for r in rows}
    accs = [r["acc"] for r in rows]
    best = {}
    for line in open(domtbl):
        if line[0] != "#":
            c = line.split()
            acc = c[4].split(".")[0]
            k, s = (read_key(c[0]), acc), float(c[13])
            if acc in ga and (k not in best or s > best[k][0]):
                best[k] = (s, (int(c[16]) - int(c[15]) + 1) / int(c[5]))
    truth = {k: (next(i for i, lo in enumerate(GA_LOWS) if s - ga[k[1]] >= lo),
                 next(i for i, lo in enumerate(COV_LOWS) if f >= lo)) for k, (s, f) in best.items()}

    def tsv(path, bar):
        return {(read_key(q), accs[int(m)]) for q, m, _, g in map(str.split, open(path)) if float(g) >= bar}

    tools = {
        "nail seeds": {(read_key(c[1]), by_name[c[0]]) for c in map(str.split, open(seeds)) if c[0] in by_name},
        "nail report": {(read_key(c[0]), c[1].split(".")[0]) for c in map(str.split, open(tbl)) if c[0][0] != "#"},
        "Mode A": tsv(a_tsv, float(a_bar)),
        "Mode B": tsv(b_tsv, float(b_bar)),
    }

    def table(title, names, axis):
        print(f"\n{title}\n{'':12}" + "".join(f"{n:>9}" for n in names))
        print(f"{'pairs':12}" + "".join(f"{sum(v[axis] == i for v in truth.values()):9d}" for i in range(len(names))))
        for t, found in tools.items():
            cells = []
            for i in range(len(names)):
                ks = [k for k, v in truth.items() if v[axis] == i]
                cells.append(f"{sum(k in found for k in ks) / max(len(ks), 1):9.3f}")
            print(f"{t:12}" + "".join(cells))

    table("Recall by GA bin", GA_NAMES, 0)
    table("Recall by HMM coverage of hmmsearch's best domain", COV_NAMES, 1)
    # Depth: reads per target (model) in the truth. Per depth bin: recall pooled over the targets' reads,
    # the mean of per-target recall, and the share of targets with at least one read found.
    for label, keep in (("GA hits", lambda v: v[0] == 0), ("all hmmsearch hits, E <= 1000", lambda _: True)):
        reads = {}
        for k, v in truth.items():
            if keep(v):
                reads.setdefault(k[1], []).append(k)
        print(f"\nRecall by target depth (reads per model in the truth, {label}): pooled / mean per target / detected")
        print(f"{'':12}" + "".join(f"{n:>20}" for n in DEPTH_NAMES))
        groups = [[a for a, ks in reads.items() if next(i for i, lo in enumerate(DEPTH_LOWS) if len(ks) >= lo) == b]
                  for b in range(len(DEPTH_LOWS))]
        print(f"{'targets':12}" + "".join(f"{f'{len(g)} ({sum(len(reads[a]) for a in g)} reads)':>20}" for g in groups))
        for t, found in tools.items():
            cells = []
            for g in groups:
                per = [sum(k in found for k in reads[a]) / len(reads[a]) for a in g]
                pooled = sum(sum(k in found for k in reads[a]) for a in g) / max(sum(len(reads[a]) for a in g), 1)
                cells.append(f"{pooled:.3f} {sum(per) / max(len(per), 1):.3f} {sum(p > 0 for p in per) / max(len(per), 1):.3f}")
            print(f"{t:12}" + "".join(f"{c:>20}" for c in cells))
    print("\nRecall by GA bin (rows) x coverage (columns), counts in brackets")
    for t, found in tools.items():
        print(f"\n{t}\n{'':9}" + "".join(f"{n:>16}" for n in COV_NAMES))
        for gi, gn in enumerate(GA_NAMES):
            cells = []
            for ci in range(len(COV_NAMES)):
                ks = [k for k, v in truth.items() if v == (gi, ci)]
                cells.append(f"{sum(k in found for k in ks) / max(len(ks), 1):8.3f} ({len(ks):5d})")
            print(f"{gn:9}" + "".join(cells))


if __name__ == "__main__":
    main()
