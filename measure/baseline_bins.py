"""Recall and candidates of the baseline tools by match strength, as the dashboard bins them.

Usage: baseline_bins.py TRUTH_DIR BASELINE_DIR MODELS.tsv NREADS REAL_DIR NREAL OUT.json
Streams each tool's output once and keeps only hits hmmsearch reported, plus one integer per other
candidate: relaxed nail settings give tens of millions of seeds.
"""

import csv
import json
import sys
from array import array
from collections import defaultdict
from pathlib import Path

import numpy as np

E_TOPS = (1e-10, 1e-5, 1e-3, 1, 10, 1000, float("inf"))
SCORE_LOWS = (60, 40, 30, 25, 20, 15, -900, float("-inf"))
GA_LOWS = (30, 10, 0, -5, -10, -20, -900, float("-inf"))
SCHEMES = {"ga": len(GA_LOWS), "evalue": len(E_TOPS), "score": len(SCORE_LOWS)}
GA = {}


def bin_of(scheme, e, s, acc=None):
    if scheme == "ga":
        d = float("-inf") if s <= -900 else s - GA[acc]
        return next(i for i, low in enumerate(GA_LOWS) if d >= low)
    if scheme == "evalue":
        return next(i for i, t in enumerate(E_TOPS) if e <= t)
    return next(i for i, low in enumerate(SCORE_LOWS) if s >= low)


def hits(path):
    best = {}
    with open(path) as f:
        for line in f:
            if line[0] != "#":
                c = line.split()
                k = (c[0], c[4].split(".")[0])
                v = (float(c[6]), float(c[13]))
                if k not in best or v[1] > best[k][1]:
                    best[k] = v
    return best


class Side:
    """What one tool passed on one set of sequences."""

    def __init__(self, table):
        self.table, self.found, self.other = table, set(), array("q")

    def add(self, q, a, index):
        if (q, a) in self.table:
            self.found.add((q, a))
        elif q[0] != "F" and a in index:
            self.other.append(int(q[1:].split("_")[0]) * 40000 + index[a])

    def candidates(self, scheme, nreads, index):
        """Per read: each (read, model) once, in the bin of its strongest ORF."""
        seen = {}
        for q, a in self.found:
            k = (int(q[1:].split("_")[0]), a)
            v = self.table[(q, a)]
            if k not in seen or v[1] > seen[k][1]:
                seen[k] = v
        out = [0.0] * SCHEMES[scheme]
        for (_, a), (e, s) in seen.items():
            out[bin_of(scheme, e, s, a)] += 1 / nreads
        rated = np.array([r * 40000 + index[a] for r, a in seen], dtype=np.int64)
        out[-1] += np.setdiff1d(np.frombuffer(self.other, dtype=np.int64), rated).size / nreads
        return out

    def recall(self, scheme):
        tot, ok = [0] * SCHEMES[scheme], [0] * SCHEMES[scheme]
        for k, (e, s) in self.table.items():
            b = bin_of(scheme, e, s, k[1])
            tot[b] += 1
            ok[b] += k in self.found
        return [o / x if x else None for o, x in zip(ok, tot)]


def read_pairs(path, swap, acc):
    with open(path) as f:
        for line in f:
            if line[0] != "#":
                c = line.split(None, 2)
                q, t = (c[1], c[0]) if swap else (c[0], c[1])
                yield q, acc(t)


def main():
    truth, base, models, nreads, real, nreal, out = sys.argv[1:8]
    truth, base, real, nreads, nreal = Path(truth), Path(base), Path(real), int(nreads), int(nreal)
    with open(models) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    name2acc = {r["name"]: r["acc"] for r in rows}
    index = {r["acc"]: i for i, r in enumerate(rows)}
    GA.update({r["acc"]: float(r["ga"]) for r in rows})
    acc = lambda t: name2acc.get(t, t.split(".")[0])  # noqa: E731

    with open(truth / "frag.tsv") as f:
        frags = list(csv.DictReader(f, delimiter="\t"))
    tables = {
        "frag": hits(truth / "frag.e1000.domtbl"),
        "reads": hits(truth / "reads.e1000.domtbl"),
        "decoys": hits(truth / "decoy.e1000.domtbl"),
        "real": hits(real / "subset.e1000.domtbl"),
        "realdecoys": hits(real / "decoy.e1000.domtbl") if (real / "decoy.e1000.domtbl").exists() else {},
    }
    # A fragment's source family counts even when hmmsearch misses it
    for r in frags:
        tables["frag"].setdefault((r["id"], r["acc"]), (1e9, -999.0))
    scale = {"reads": nreads, "decoys": nreads, "real": nreal, "realdecoys": nreal}
    # Tool name, then per set the output and whether its first column is the model
    tools = [
        ("nail seeds", {"frag": (base / "nail_frag.seeds.tsv", True), "reads": (base / "nail.seeds.tsv", True), "decoys": (base / "nail_decoy.seeds.tsv", True), "real": (real / "nail_subset.seeds.tsv", True), "realdecoys": (real / "nail_decoy.seeds.tsv", True)}),
        ("nail report", {"frag": (base / "nail_frag.tbl", False), "reads": (base / "nail.tbl", False), "decoys": (base / "nail_decoy.tbl", False), "real": (real / "nail_subset.tbl", False), "realdecoys": (real / "nail_decoy.tbl", False)}),
        ("MMseqs2 s5.7 k5", {"reads": (base / "pref_s5.7k5.tsv", False)}),
        ("MMseqs2 s7.5 k5", {"reads": (base / "pref_s7.5k5.tsv", False)}),
    ]
    results = []

    def finish(name, sides):
        t = {"tool": name, "recall": {}, "candidates": {}}
        for scheme, nb in SCHEMES.items():
            if "frag" in sides:
                groups = defaultdict(lambda: ([0.0] * nb, [0.0] * nb))
                for r in frags:
                    k = (r["id"], r["acc"])
                    e, s = sides["frag"].table.get(k, (1e9, -999.0))
                    b, w = bin_of(scheme, e, s, k[1]), float(r["weight"])
                    for g in (f"f{r['length']}", "fragments"):
                        groups[g][0][b] += w
                        groups[g][1][b] += w * (k in sides["frag"].found)
                for g, (tot, ok) in groups.items():
                    t["recall"].setdefault(g, {})[scheme] = [o / x if x else None for o, x in zip(ok, tot)]
            for s in ("reads", "real"):
                if s in sides:
                    t["recall"].setdefault(s, {})[scheme] = sides[s].recall(scheme)
            t["candidates"][scheme] = {s: sides[s].candidates(scheme, scale[s], index) for s in scale if s in sides}
        results.append(t)
        print(name, flush=True)
        for s in ("reads", "real"):
            if s in t["recall"]:
                print(" ", s, [None if x is None else round(x, 3) for x in t["recall"][s]["ga"]], flush=True)
        if "fragments" in t["recall"]:
            print("  fragments", [None if x is None else round(x, 3) for x in t["recall"]["fragments"]["ga"]], flush=True)
        for s, c in t["candidates"]["ga"].items():
            print("  cand", s, [round(x, 3) for x in c], flush=True)

    for name, sources in tools:
        sides = {}
        for s, (path, swap) in sources.items():
            if path.exists():
                sides[s] = Side(tables[s])
                for q, a in read_pairs(path, swap, acc):
                    sides[s].add(q, a, index)
        finish(name, sides)
    # Relaxed nail settings: one run over reads, decoys (renamed D...) and fragments, split here.
    route = {"F": "frag", "R": "reads", "D": "decoys"}
    for table in sorted((base / "sweep").glob("*.tbl")):
        tag = table.stem
        for suffix, swap, label in ((".seeds.tsv", True, "seeds"), (".tbl", False, "report")):
            sides = {s: Side(tables[s]) for s in route.values()}
            for q, a in read_pairs(table.with_suffix(suffix), swap, acc):
                if q[0] in route:
                    sides[route[q[0]]].add("R" + q[1:] if q[0] == "D" else q, a, index)
            # The same setting on the real reads, when it was run there
            for s, name in (("real", "subset"), ("realdecoys", "decoy")):
                path = real / f"nail_{tag}_{name}{suffix}"
                if path.exists():
                    sides[s] = Side(tables[s])
                    for q, a in read_pairs(path, swap, acc):
                        sides[s].add(q, a, index)
            finish(f"nail {tag} {label}", sides)
    Path(out).write_text(json.dumps({"tools": results}))


if __name__ == "__main__":
    main()
