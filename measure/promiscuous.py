"""What rejecting promiscuous k-mers would cost and save, for k4 T7 at the diagonal bars.

Usage: promiscuous.py DATA_DIR OUT.json

Rules, each sweeping hits removed: longest index lists (index), list length x frequency under HMMER's
background (expected), hits measured on this sample (sample), and masking tandem repeats in the reads:
a k-mer seen again within p residues in its ORF (repeat). A pair survives when any seed on a diagonal
at the bar survives. Only `sample` and `repeat` look at the reads, and they could at query time.
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np

from sketch_score import BIN_LABELS, GA_BINS, read_truth

K, SPACE = 4, 160000
A = b"ACDEFGHIKLMNPQRSTVWY"
# HMMER's amino acid background, in A's order.
BG = np.array([.0787945, .0151600, .0535222, .0668298, .0397062, .0695071, .0229198, .0590092, .0594422, .0965125,
               .0237718, .0414386, .0482904, .0395639, .0540978, .0683364, .0540687, .0673417, .0114135, .0304133])
LEVELS = (0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
PERIODS = (2, 4, 8, 16)
BARS = (22, 24, 26)
TOP_FAMILIES = 8


def scan_orfs(path):
    """Per 4-mer occurrence: k-mer frequency, and distance to the nearest same k-mer in the ORF."""
    table = {c: i for i, c in enumerate(A)}
    freq = np.zeros(SPACE, dtype=np.int64)
    near, reads = {}, set()
    with open(path, "rb") as f:
        for head, seq in zip(f, f):
            name = head[1:].split()[0].decode()
            reads.add(name.split("_")[0])
            codes, last, dist = [], {}, {}
            seq = seq.strip()
            for j in range(len(seq) - K + 1):
                digits = [table.get(c) for c in seq[j : j + K]]
                if None in digits:
                    continue
                code = digits[0] * 8000 + digits[1] * 400 + digits[2] * 20 + digits[3]
                freq[code] += 1
                if code in last:
                    gap = j - last[code]
                    dist[j] = min(dist.get(j, 255), gap)
                    dist[last[code]] = min(dist.get(last[code], 255), gap)
                last[code] = j
                codes.append((j, code))
            if dist:
                near[name] = (dist, [(j, c) for j, c in codes if j in dist])
    return freq, near, len(reads)


def repeat_hits(near, lens, period):
    """Hits from occurrences a repeat mask of this period would remove."""
    return sum(int(lens[c]) for dist, occ in near.values() for j, c in occ if dist[j] <= period)


def read_seeds(path, near):
    """Arrays over logged hits: pair id, k-mer, diagonal, nearest repeat; and the pair keys."""
    pairs, pid, code, diag, rep = {}, [], [], [], []
    with open(path) as f:
        for line in f:
            q, j, c, m, d = line.split("\t")
            key = (q.split("_")[0], int(m))
            pid.append(pairs.setdefault(key, len(pairs)))
            code.append(int(c))
            diag.append(float(d))
            entry = near.get(q)
            rep.append(entry[0].get(int(j), 255) if entry else 255)
    keys = list(pairs)
    return np.array(pid), np.array(code), np.array(diag, dtype=np.float32), np.array(rep, dtype=np.uint8), keys


def survivors(seeds, bar, rejected, period):
    pid, code, diag, rep, _ = seeds
    keep = diag >= bar
    if rejected is not None:
        keep &= ~rejected[code]
    if period is not None:
        keep &= rep > period
    return np.unique(pid[keep])


def rate(seeds, seeds_decoy, truth_bins, ga_ok, names, n_reads, bar, rule):
    """One rule at one bar. rule: (rejected k-mers or None, repeat period or None) for reads, then reversed."""
    base = survivors(seeds, bar, None, None)
    now = survivors(seeds, bar, *rule[0])
    keys = seeds[4]
    per_bin = np.zeros(len(GA_BINS))
    family = np.zeros(len(names), dtype=np.int64)
    for p in now:
        k = keys[p]
        b = truth_bins.get(k)
        if b is not None:
            per_bin[b] += 1
            if b == 0:
                family[k[1]] += 1
    baseline = np.zeros(len(names), dtype=np.int64)
    for p in base:
        k = keys[p]
        if truth_bins.get(k) == 0:
            baseline[k[1]] += 1
    totals = np.bincount(list(truth_bins.values()), minlength=len(GA_BINS))
    common = ga_ok >= 5
    recall = family[common] / ga_ok[common]
    drop = baseline - family
    order = np.argsort(-drop)[:TOP_FAMILIES]
    bray = np.abs(baseline - family).sum() / max((baseline + family).sum(), 1)
    return {
        "recall": (per_bin / totals).tolist(),
        "candidates": len(now) / n_reads,
        "decoy": len(survivors(seeds_decoy, bar, *rule[1])) / n_reads,
        "p10": float(np.quantile(recall, 0.1)),
        "bray": float(bray),
        "families_lost": int(((baseline > 0) & (family == 0)).sum()),
        "worst": [[names[m], int(baseline[m]), int(family[m])] for m in order if drop[m] > 0],
    }


def main():
    data, out = Path(sys.argv[1]), sys.argv[2]
    with open(data / "prep/models.tsv") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    index = {r["acc"]: i for i, r in enumerate(rows)}
    ga = np.array([float(r["ga"]) for r in rows])
    names = [r["name"] for r in rows]
    real = data / "real/ERR10889147"
    truth_bins = {}
    for key, s in read_truth(real / "subset.e1000.domtbl", index).items():
        for i, low in enumerate(GA_BINS):
            if s - ga[key[1]] >= low:
                truth_bins[key] = i
                break
    ga_ok = np.zeros(len(rows), dtype=np.int64)
    for key, b in truth_bins.items():
        if b == 0:
            ga_ok[key[1]] += 1

    with open(data / "index/k4T7.idx", "rb") as f:
        np.fromfile(f, np.uint64, 2)
        lens = np.diff(np.fromfile(f, np.uint64, SPACE + 1).astype(np.int64))
    freq, near, n_reads = scan_orfs(real / "orfs.subset.faa")
    freq_d, near_d, _ = scan_orfs(real / "decoy.faa")
    seeds = read_seeds(data / "sketch/seeds.real.tsv", near)
    seeds_decoy = read_seeds(data / "sketch/seeds.decoy.tsv", near_d)

    prob = np.ones(SPACE)
    digits = np.arange(SPACE)
    for _ in range(K):
        prob *= BG[digits % 20]
        digits //= 20
    hits, hits_d = freq * lens, freq_d * lens
    total = int(hits.sum())
    ranks = {"index": lens.astype(float), "expected": lens * prob, "sample": hits.astype(float)}
    result = {"hits_per_read": total / n_reads, "bins": BIN_LABELS, "bars": BARS, "periods": PERIODS, "levels": LEVELS, "rules": []}
    for bar in BARS:
        result.setdefault("base", {})[bar] = rate(seeds, seeds_decoy, truth_bins, ga_ok, names, n_reads, bar, ((None, None), (None, None)))
    for name, score in ranks.items():
        order = np.argsort(-score)
        removed = np.cumsum(hits[order]) / total
        decoy_order = np.argsort(-hits_d) if name == "sample" else order
        removed_d = np.cumsum(hits_d[decoy_order]) / hits_d.sum()
        for level in LEVELS:
            take = int(np.searchsorted(removed, level)) + 1
            rejected = np.zeros(SPACE, dtype=bool)
            rejected[order[:take]] = True
            decoy_rejected = np.zeros(SPACE, dtype=bool)
            decoy_rejected[decoy_order[: int(np.searchsorted(removed_d, level)) + 1]] = True
            result["rules"].append({
                "rule": name, "kmers": take, "removed": float(removed[take - 1]),
                "bars": {bar: rate(seeds, seeds_decoy, truth_bins, ga_ok, names, n_reads, bar, ((rejected, None), (decoy_rejected, None))) for bar in BARS},
            })
    for period in PERIODS:
        result["rules"].append({
            "rule": "repeat", "period": period, "removed": repeat_hits(near, lens, period) / total,
            "bars": {bar: rate(seeds, seeds_decoy, truth_bins, ga_ok, names, n_reads, bar, ((None, period), (None, period))) for bar in BARS},
        })
    # Rejecting nothing must return the search's own pairs.
    assert len(survivors(seeds, 22, np.zeros(SPACE, dtype=bool), None)) == len(survivors(seeds, 22, None, None))
    Path(out).write_text(json.dumps(result, separators=(",", ":")))


if __name__ == "__main__":
    main()
