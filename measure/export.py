"""Scan outputs to the small files the dashboard loads.

Usage: export.py SCAN_DIR PREP_DIR ABUNDANCE.tsv OUT
"""

import csv
import gzip
import json
import shutil
import sys
from collections import Counter
from pathlib import Path

import numpy as np

SHAPES = {"k4": (0, 1, 2, 3), "k5": (0, 1, 2, 3, 4), "k6": (0, 1, 2, 3, 4, 5),
          "k7": (0, 1, 2, 3, 4, 5, 6), "s5": (0, 1, 3, 5, 6), "s6": (0, 1, 2, 4, 6, 7),
          "p33": (0, 1, 2, 6, 7, 8)}
NT, NN = 21, 9
NG = NT + 3 * NN
FRACS = [1.0, 0.5, 0.25, 0.1]


def thin_background(src, groups_src, dest, groups_dest, bar=18.0):
    """Keep anchors, and seeds on diagonals of at least `bar` bits: lower bars pass thousands of
    models and double the download. Layout is scan's write_rows."""
    b = gzip.open(src).read()
    ng, nr = np.frombuffer(b[:8], dtype=np.uint32)
    offs = np.frombuffer(b, np.uint32, ng + 1, 8)
    p = 8 + 4 * (int(ng) + 1)
    wide = {}
    for name in ("col", "hash"):
        wide[name] = np.frombuffer(b, np.uint16, nr, p)
        p += 2 * int(nr)
    cols = {}
    for name in ("j", "score", "diag", "rank_lr", "rank_p", "logp", "info", "word"):
        cols[name] = np.frombuffer(b, np.uint8, nr, p)
        p += int(nr)
    keep = (cols["diag"] >= int(bar * 2)) | (cols["word"] & 2 > 0)
    group = np.repeat(np.arange(ng), np.diff(offs))[keep]
    kept, inv = np.unique(group, return_inverse=True)
    new = np.zeros(len(kept) + 1, np.uint32)
    np.add.at(new, inv + 1, 1)
    new = np.cumsum(new, dtype=np.uint32)
    with gzip.open(dest, "wb", compresslevel=9) as f:
        f.write(np.array([len(kept), keep.sum()], np.uint32).tobytes())
        f.write(new.tobytes())
        for v in wide.values():
            f.write(v[keep].tobytes())
        for v in cols.values():
            f.write(v[keep].tobytes())
    pairs = np.frombuffer(gzip.open(groups_src).read(), np.uint16).reshape(-1, 2)
    with gzip.open(groups_dest, "wb", compresslevel=9) as f:
        f.write(pairs[kept].tobytes())


def profile_chunks(prep, out, chunks=200):
    """Log-odds per model in quarter bits, in pieces fetched on demand by the family view."""
    offs = np.fromfile(prep / "offs.u32", dtype=np.uint32)
    lo = np.memmap(prep / "profile.f32", dtype=np.float32, mode="r").reshape(-1, 21)[:, :20]
    (out / "profile").mkdir(exist_ok=True)
    n = len(offs) - 1
    for c in range(chunks):
        a, b = c * n // chunks, (c + 1) * n // chunks
        q = np.clip(np.round(lo[offs[a]:offs[b]] * 4), -127, 127).astype(np.int8)
        with gzip.open(out / "profile" / f"{c}.bin.gz", "wb", compresslevel=9) as f:
            f.write(q.tobytes())


def main():
    scan, prep, abundance, out = (Path(a) for a in sys.argv[1:5])
    out.mkdir(parents=True, exist_ok=True)
    with open(prep / "models.tsv") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    count = Counter()
    with open(abundance) as f:
        for line in f:
            acc, n = line.split()
            count[acc] = int(n)
    models = {k: [r[k] for r in rows] for k in ("acc", "name", "clan")}
    for k in ("length", "ga", "nseq", "effn", "mu", "lambda", "mean_info", "top_info", "conserved"):
        models[k] = [float(r[k]) for r in rows]
    models["domains"] = [count[r["acc"]] for r in rows]
    with gzip.open(out / "models.json.gz", "wt") as f:
        json.dump(models, f, separators=(",", ":"))

    for shape in SHAPES:
        t = np.frombuffer(gzip.open(scan / "tables" / f"{shape}.tab.f32.gz").read(), dtype=np.float32)
        t = t.reshape(len(rows), NG, 3)
        # log2 in 1/4-bit steps as signed 16 bits: precise enough, and gzip takes most of it
        log = np.clip(np.round(np.log2(t[:, :, :2] + 1e-12) * 8), -32000, 32000).astype(np.int16)
        cov = np.round(np.clip(t[:, :, 2], 0, 65000)).astype(np.uint16)
        with gzip.open(out / f"{shape}.tab.gz", "wb", compresslevel=9) as f:
            f.write(log.tobytes())
            f.write(cov.tobytes())
        strata = np.fromfile(scan / "tables" / f"{shape}.strata.f64").reshape(len(FRACS), NT, 3)
        ratio = strata[:, :, :2] / np.maximum(strata[0, :, :2], 1e-30)
        np.round(ratio, 4).astype(np.float32).tofile(out / f"{shape}.strata.bin")
        shutil.copy(scan / "bg" / f"{shape}.words.bin", out / f"{shape}.words.bin")
        shutil.copy(scan / "truth" / f"{shape}.rows.gz", out / f"{shape}.rows.gz")
        thin_background(scan / "bg" / f"{shape}.bg.gz", scan / "bg" / f"{shape}.bggroups.gz",
                        out / f"{shape}.bg.gz", out / f"{shape}.bggroups.gz")

    # pairs, with each query's best diagonal over all models
    with open(scan / "pairs.meta.tsv") as f:
        pairs = list(csv.DictReader(f, delimiter="\t"))
    best = [line.split("\t") for line in open(scan / "pairs.best.tsv").read().splitlines()]
    assert len(best) == len(pairs)
    doc = {
        "set": [p["set"] for p in pairs],
        "model": [int(p["model"]) for p in pairs],
        "weight": [float(p["weight"]) for p in pairs],
        "evalue": [float(p["evalue"]) for p in pairs],
        "score": [float(p["score"]) for p in pairs],
        "ga": [int(p["ga_pass"]) for p in pairs],
        "best": [float(b[1]) for b in best],
        "bestModel": [int(b[2]) for b in best],
    }
    with gzip.open(out / "pairs.json.gz", "wt") as f:
        json.dump(doc, f, separators=(",", ":"))
    # background reads: best over each read's ORFs
    byread = {}
    for read, score, model in (line.split("\t") for line in open(scan / "bg.best.tsv").read().splitlines()):
        if read not in byread or float(score) > byread[read][0]:
            byread[read] = (float(score), int(model))
    orfs = [line.split("\t") for line in open(scan / "bg.orfs.tsv").read().splitlines()]
    reads = max(int(o[0]) for o in orfs) + 1
    meta = {
        "reads": reads,
        "residuesPerRead": sum(len(o[1]) for o in orfs) / reads,
        "orfsPerRead": len(orfs) / reads,
        "readBest": [byread.get(str(i), (0, 0))[0] for i in range(reads)],
        "readBestModel": [byread.get(str(i), (0, 0))[1] for i in range(reads)],
        "shapes": {k: list(v) for k, v in SHAPES.items()},
        "fracs": FRACS,
        "anchorFloor": {"k4": 9, "k5": 9.5, "k6": 10.5, "k7": 11.5, "s5": 9.5, "s6": 10.5, "p33": 10.5},
        "maxWindow": 6,
    }
    (out / "meta.json").write_text(json.dumps(meta))
    profile_chunks(prep, out)
    print(sum(p.stat().st_size for p in out.iterdir()) / 1e6, "MB")


if __name__ == "__main__":
    main()
