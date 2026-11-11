"""Score a tool's (query, model) candidates against hmmsearch on the same read ORFs.

Usage: baseline_score.py CANDIDATES.tsv READS_DOMTBL MODELS.tsv NREADS TOOL SETTINGS SECONDS_PER_READ >> baselines.jsonl
Candidates: query and target name or accession in the first two columns. ORF names are R<read>_<orf>.
"""

import csv
import json
import sys
from collections import defaultdict

import numpy as np


def main():
    cand_path, table, models, nreads, tool, settings, ms = sys.argv[1:8]
    nreads = int(nreads)
    with open(models) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    ga = {r["acc"]: float(r["ga"]) for r in rows}
    name2acc = {r["name"]: r["acc"] for r in rows}
    truth = {}
    with open(table) as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.split()
            acc, e, s = c[4].split(".")[0], float(c[6]), float(c[13])
            key = (c[0], acc)
            old = truth.get(key, (1e9, -1e9))
            truth[key] = (min(old[0], e), max(old[1], s))
    cands = set()
    with open(cand_path) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 2:
                c = line.split()
            q, t = c[0], c[1]
            cands.add((q, name2acc.get(t, t.split(".")[0])))
    per_read = len({(q.split("_")[0], t) for q, t in cands}) / nreads
    out = {"tool": tool, "settings": settings, "candidates": per_read, "ms_per_read": float(ms)}
    for name, keep in (("e3", lambda v, a: v[0] <= 1e-3), ("ga", lambda v, a: v[1] >= ga[a])):
        wanted = [k for k, v in truth.items() if keep(v, k[1])]
        out[f"recall_{name}"] = sum(k in cands for k in wanted) / max(len(wanted), 1)
        fam = defaultdict(list)
        for k in wanted:
            fam[k[1]].append(k in cands)
        vals = np.array([np.mean(v) for v in fam.values() if len(v) >= 3])
        out[f"family_p10_{name}"] = float(np.percentile(vals, 10)) if len(vals) else None
        out[f"families_{name}"] = len(vals)
    print(json.dumps(out))


if __name__ == "__main__":
    main()
