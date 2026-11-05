"""Recall and candidates of the baseline tools by match strength, as the dashboard bins them.

Usage: baseline_bins.py TRUTH_DIR BASELINE_DIR PAIRS_META MODELS.tsv NREADS OUT.json
"""

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

E_TOPS = (1e-10, 1e-5, 1e-3, 1, 10, 1000, float("inf"))
SCORE_LOWS = (60, 40, 30, 25, 20, 15, -900, float("-inf"))
GA_LOWS = (30, 10, 0, -5, -10, -20, -900, float("-inf"))
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


def main():
    truth, base, meta_path, models, nreads, out = sys.argv[1:7]
    truth, base, nreads = Path(truth), Path(base), int(nreads)
    with open(models) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    name2acc = {r["name"]: r["acc"] for r in rows}
    GA.update({r["acc"]: float(r["ga"]) for r in rows})
    acc = lambda t: name2acc.get(t, t.split(".")[0])  # noqa: E731

    def load(path, swap):
        got = set()
        if path.exists():
            with open(path) as f:
                for line in f:
                    if line[0] == "#":
                        continue
                    c = line.split()
                    q, t = (c[1], c[0]) if swap else (c[0], c[1])
                    got.add((q, acc(t)))
        return got

    with open(truth / "frag.tsv") as f:
        frags = list(csv.DictReader(f, delimiter="\t"))
    frag_hits = hits(truth / "frag.e1000.domtbl")
    read_hits = hits(truth / "reads.e1000.domtbl")
    decoy_hits = hits(truth / "decoy.e1000.domtbl")
    tools = [
        ("nail seeds", base / "nail_frag.seeds.tsv", True, base / "nail.seeds.tsv", True, base / "nail_decoy.seeds.tsv", True),
        ("nail report", base / "nail_frag.tbl", False, base / "nail.tbl", False, base / "nail_decoy.tbl", False),
        ("MMseqs2 s5.7 k5", None, False, base / "pref_s5.7k5.tsv", False, None, False),
        ("MMseqs2 s7.5 k5", None, False, base / "pref_s7.5k5.tsv", False, None, False),
    ]
    # Relaxed nail settings: one run over reads, decoys (renamed D...) and fragments, split here.
    sweep = base / "sweep"
    if sweep.exists():
        for table in sorted(sweep.glob("*.tbl")):
            tag = table.stem
            for suffix, swap, label in ((".seeds.tsv", True, "seeds"), (".tbl", False, "report")):
                src = sweep / f"{tag}{suffix}"
                parts = {"F": [], "R": [], "D": []}
                with open(src) as f:
                    for line in f:
                        if line[0] == "#":
                            continue
                        c = line.split()
                        q = c[1] if swap else c[0]
                        if q[0] in parts:
                            parts[q[0]].append(f"{c[0]}\t{c[1]}\n")
                paths = {}
                for k, lines in parts.items():
                    paths[k] = sweep / f".{tag}.{label}.{k}"
                    if k == "D":
                        lines = [x.replace("\tD", "\tR", 1) if swap else "R" + x[1:] for x in lines]
                    paths[k].write_text("".join(lines))
                tools.append((f"nail {tag} {label}", paths["F"], swap, paths["R"], swap, paths["D"], swap))
    results = []
    for name, fpath, fswap, rpath, rswap, dpath, dswap in tools:
        t = {"tool": name, "recall": {}, "candidates": {}}
        for scheme, nb in (("ga", len(GA_LOWS)), ("evalue", len(E_TOPS)), ("score", len(SCORE_LOWS))):
            if fpath:
                cands = load(fpath, fswap)
                groups = {"f50": frags and [], "f33": []}
                for r in frags:
                    k = (r["id"], r["acc"])
                    e, s = frag_hits.get(k, (1e9, -999.0))
                    groups[f"f{r['length']}"].append((k, e, s, float(r["weight"])))
                groups["fragments"] = groups["f50"] + groups["f33"]
                for setname, items in groups.items():
                    tot, ok = [0.0] * nb, [0.0] * nb
                    for k, e, s, w in items:
                        b = bin_of(scheme, e, s, k[1])
                        tot[b] += w
                        ok[b] += w * (k in cands)
                    t["recall"].setdefault(setname, {})[scheme] = [o / x if x else None for o, x in zip(ok, tot)]
            cands = load(rpath, rswap)
            tot, ok = [0] * nb, [0] * nb
            for k, (e, s) in read_hits.items():
                b = bin_of(scheme, e, s, k[1])
                tot[b] += 1
                ok[b] += k in cands
            t["recall"].setdefault("reads", {})[scheme] = [o / x if x else None for o, x in zip(ok, tot)]
            spread = defaultdict(lambda: [0.0] * nb)
            per_read = lambda q: q.split("_")[0]  # noqa: E731
            for label, cs, table in (("reads", cands, read_hits), ("decoys", load(dpath, dswap) if dpath else None, decoy_hits)):
                if cs is None:
                    continue
                seen = {}
                for q, a in cs:
                    e, s = table.get((q, a), (1e9, -999.0))
                    k = (per_read(q), a)
                    if k not in seen or s > seen[k][1]:
                        seen[k] = (e, s)
                for (_, a), (e, s) in seen.items():
                    spread[label][bin_of(scheme, e, s, a)] += 1 / nreads
            t["candidates"][scheme] = dict(spread)
        results.append(t)
    Path(out).write_text(json.dumps({"tools": results}))
    for t in results:
        print(t["tool"], "reads", [None if x is None else round(x, 3) for x in t["recall"]["reads"]["ga"]])
        if "fragments" in t["recall"]:
            print(" fragments", [None if x is None else round(x, 3) for x in t["recall"]["fragments"]["ga"]])
        print(" cand", {k: [round(x, 3) for x in v] for k, v in t["candidates"]["ga"].items()})


if __name__ == "__main__":
    main()
