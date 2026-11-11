"""Recall and precision of whole-sample k-mer sets as a model pre-filter, for the dashboard.

Usage: sketch_score.py DATA_DIR OUT.json

Reads DATA_DIR/sketch/*.fwd.tsv and *.rev.tsv (`filter sketch`) and their timings in runs.log. A model
is present when any read has a domain at or above its GA: hmmsearch on 50,000 reads, nail on all reads.
Rules: best = top k-mer score in the set; best-null = that less the same model's score on reversed
reads. The bound lumps the read-level diagonal and gapped scores per model, which no k-mer set can beat.
"""

import csv
import json
import re
import sys
from pathlib import Path

import numpy as np

GA_BINS = (0, -5, -10, -20)
BIN_LABELS = ("≥GA", "GA−5–GA", "GA−10–−5", "GA−20–−10")
MISSING = -99.0


def read_sketch(path, n):
    """Top k-mer score per model; MISSING where the set holds none at or above the sketch's T."""
    best = np.full(n, MISSING)
    with open(path) as f:
        for line in f:
            m, b = line.split("\t")[:2]
            best[int(m)] = max(float(b), MISSING)
    return best


def read_truth(path, index):
    """Best domain score per (read, model), from an hmmsearch domtblout or a nail tbl."""
    top = {}
    with open(path) as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.split()
            q, acc, s = (c[0], c[4], float(c[13])) if len(c) > 13 else (c[0], c[1], float(c[6]))
            key = (q.split("_")[0], index[acc.split(".")[0]])
            top[key] = max(top.get(key, MISSING), s)
    return top


def bin_hits(truth, ga):
    """Hits per GA bin per model (bins by columns)."""
    hits = np.zeros((len(GA_BINS), len(ga)), dtype=np.int64)
    for (_, m), s in truth.items():
        for i, low in enumerate(GA_BINS):
            if s - ga[m] >= low:
                hits[i, m] += 1
                break
    return hits


def sweep(score, null, hits, bars):
    """At each bar: models kept, on the null too, and recall of models and of hits."""
    present = hits[0] > 0
    common = hits[0] >= 5
    rows = []
    for bar in bars:
        kept = score >= bar
        n_kept = int(kept.sum())
        rows.append({
            "bar": float(bar),
            "kept": n_kept,
            "decoy": int((null >= bar).sum()),
            "recall": present[kept].sum() / max(present.sum(), 1),
            "recall5": common[kept].sum() / max(common.sum(), 1),
            "hits": [hits[i][kept].sum() / max(hits[i].sum(), 1) for i in range(len(GA_BINS))],
            "precision": present[kept].sum() / max(n_kept, 1),
        })
    assert all(a["kept"] >= b["kept"] and a["recall"] >= b["recall"] for a, b in zip(rows, rows[1:]))
    return rows


def read_timings(path):
    """From runs.log, by output path."""
    out = {}
    for line in Path(path).read_text().splitlines():
        if " out=" in line:
            fields = dict(re.findall(r"(\w+)=(\S+)", line))
            out[Path(fields["out"]).name] = fields
    return out


def read_pairs(path, n):
    """Best diagonal and gapped score per model, from `filter search`."""
    d, g = np.full(n, MISSING), np.full(n, MISSING)
    with open(path) as f:
        for line in f:
            _, m, a, b = line.split()
            m = int(m)
            d[m] = max(d[m], float(a))
            g[m] = max(g[m], float(b))
    return d, g


def main():
    data, out = Path(sys.argv[1]), sys.argv[2]
    with open(data / "prep/models.tsv") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    index = {r["acc"]: i for i, r in enumerate(rows)}
    ga = np.array([float(r["ga"]) for r in rows])
    n = len(rows)
    real = data / "real/ERR10889147"
    hits = {
        "sub": bin_hits(read_truth(real / "subset.e1000.domtbl", index), ga),
        "all": bin_hits(read_truth(real / "nail_all.tbl", index), ga),
    }
    timings = read_timings(data / "sketch/runs.log")
    runs = []
    for path in sorted((data / "sketch").glob("*.fwd.tsv")):
        m = re.match(r"(\w+)\.k(\d+)T([\d.]+)\.s(\d+)\.fwd\.tsv", path.name)
        assert m, path
        tag, k, t, scaled = m[1], int(m[2]), float(m[3]), int(m[4])
        rev = path.with_name(path.name.replace(".fwd.", ".rev."))
        if not rev.exists():
            continue
        best, null = read_sketch(path, n), read_sketch(rev, n)
        floored, floored_null = np.maximum(best, t - 0.5), np.maximum(null, t - 0.5)
        time = timings.get(path.name, {})
        runs.append({
            "depth": tag, "k": k, "T": t, "scaled": scaled,
            "best": sweep(best, null, hits[tag], np.arange(t, 61.0, 0.5)),
            "null": sweep(floored - floored_null, floored_null - floored, hits[tag], np.arange(-4.0, 16.5, 0.5)),
            **{key: float(time[key]) if key in time else None for key in ("sample_core_s", "score_core_s", "residues", "filled")},
            "absent": int(((best == MISSING) & (hits[tag][0] > 0)).sum()),
        })
    d, g = read_pairs(data / "index/real.tsv", n)
    dd, dg = read_pairs(data / "index/decoy.tsv", n)
    bars = np.arange(22.0, 61.0)
    bound = {"diagonal": sweep(d, dd, hits["sub"], bars), "gapped": sweep(g, dg, hits["sub"], bars)}
    result = {
        "models": n,
        "median_len": float(np.median([float(r["length"]) for r in rows])),
        "bins": BIN_LABELS,
        "reads": {"sub": 50000, "all": 1070914},
        "present": {t: int((h[0] > 0).sum()) for t, h in hits.items()},
        "ga_hits": {t: int(h[0].sum()) for t, h in hits.items()},
        "singletons": {t: int((h[0] == 1).sum()) for t, h in hits.items()},
        "runs": runs,
        "bound": bound,
    }
    Path(out).write_text(json.dumps(result, default=float, separators=(",", ":")))


if __name__ == "__main__":
    main()
