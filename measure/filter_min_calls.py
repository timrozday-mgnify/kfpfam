"""How few called reads can train select_frames.py's tables: kept residues and GA pairs by training size.

Usage: filter_min_calls.py PIPELINE_DIR V2_DIR. Tables train on N called reads, split in two, and score all
reads, as a sample of N called reads would train and score itself, with more reads to measure on.
"""

import random
import sys
from collections import Counter

import pyarrow.parquet as pq

pipe, v = sys.argv[1:3]
sys.path.insert(0, pipe + "/modules/local/selectframes/bin")
import select_frames as sf  # noqa: E402

reads, orfs = v + "/reads.fasta", v + "/orfs.faa"
calls = sf.read_calls(v + "/fgs.out")
rows = list(sf.by_read(reads, orfs))
items = [(r, s, [(h, *sf.nucleotides(s, h)[1:], len(p)) for h, p in o]) for r, s, o in rows]
ga = {f"{r['read']}_frame={r['frame']}_begin={r['orf_begin']}_end={r['orf_end']}" for r in pq.read_table(v + "/pairs.parquet").to_pylist()}
total = sum(x[3] for _, _, o in items for x in o)
called_ids = [i for i, (r, _, _) in enumerate(items) if r in calls]
print("training called reads\tseed\tresidues kept\tGA ORFs kept\tkept outside FGS frame")
for n in (100, 250, 500, 1000, 2500, 5000, 10000, 25000, len(called_ids)):
    for seed in range(3 if n < len(called_ids) else 1):
        chosen = random.Random(seed).sample(called_ids, n)
        coding, other = [Counter(), Counter()], [Counter(), Counter()]
        fold = {}
        for j, i in enumerate(chosen):
            fold[i] = j % 2
            r, _, o = items[i]
            for _, frame, nt, _ in o:
                (coding if frame in calls[r] else other)[j % 2].update(sf.hexamers(nt))
        tables = [sf.log_odds(coding[1 - k], other[1 - k]) for k in (0, 1)]
        kept = out = 0
        ga_kept = 0
        for i, (r, _, o) in enumerate(items):
            f = fold.get(i, i % 2)
            for h, frame, nt, length in o:
                inframe = frame in calls.get(r, ())
                keep = inframe or sum(tables[1 - f].get(x, 0.0) for x in sf.hexamers(nt)) >= -5
                if keep:
                    kept += length
                    out += not inframe
                    ga_kept += h.split()[0] in ga
        print(f"{n}\t{seed}\t{kept / total:.3f}\t{ga_kept / len(ga):.4f}\t{out}", flush=True)
