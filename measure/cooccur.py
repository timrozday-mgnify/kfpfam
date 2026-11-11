"""Do anchor k-mers co-occur in real family members beyond what their separate frequencies predict,
and does having been seen together separate true hits from chance?

Per family, k-mers come from Pfam-A.seed only (MGnify stays out of training). An anchor is a k-mer of
the query that occurs in the family's seed alignment and scores at least T at some model window. A
pair is two non-overlapping anchors of one query. It is co-seen when one seed sequence holds both, at
the query's distance give or take W residues.

Usage: cooccur.py PREP SEQS.tsv META.tsv ORFS.faa NREADS OUT.json
Truth is keyed by score against GA: 0 is GA + 30 or more, through 6, below GA - 20; 7 is not reported.
ORFS may be decoys, so real background has no reported hits to confuse it.
"""

import csv
import gzip
import json
import sys
from collections import defaultdict

import numpy as np

SHAPES = {"k5": (0, 1, 2, 3, 4)}  # ponytail: 4-mers at low T ran for hours; add them if 5-mers show a signal
THRESHOLDS = (6.0, 8.0)
WINDOW = 6
AA = "ACDEFGHIKLMNPQRSTVWY"
INDEX = {c: i for i, c in enumerate(AA)}
NSEQ_BINS = ((0, 10), (10, 30), (30, 100), (100, 10**9))


def families(path):
    acc, seqs = None, []
    with gzip.open(path, "rb") as f:
        for line in f:
            if line.startswith(b"#=GF AC"):
                acc = line.split()[2].decode().split(".")[0]
            elif line.startswith(b"//"):
                yield acc, seqs
                seqs = []
            elif line.strip() and line[0] != 35:
                s = line.split()[-1].decode("ascii", "replace").upper()
                seqs.append("".join(c for c in s if c.isalpha()))


def words(seq, cols):
    span = cols[-1] + 1
    for j in range(len(seq) - span + 1):
        k = "".join(seq[j + o] for o in cols)
        if all(c in INDEX for c in k):
            yield j, k


def best_score(lo, kmer, cols):
    length = lo.shape[0] - cols[-1]
    if length <= 0:
        return -1e9
    s = np.zeros(length, dtype=np.float32)
    for o, c in zip(cols, kmer):
        s += lo[o:o + length, INDEX[c]]
    return float(s.max())


class Tally:
    def __init__(self):
        self.n = self.single = self.pair = self.coseen = 0
        self.observed = self.expected = 0.0

    def row(self):
        return {"queries": self.n, "one_anchor": self.single, "pair": self.pair, "co_seen": self.coseen,
                "observed_co_seen": self.observed, "expected_if_independent": self.expected}


def group(seen):
    out = defaultdict(list)
    for s, p in seen:
        out[s].append(p)
    return out


def evaluate(anchors, occ, nseq, span, tally, positions):
    tally.single += 1
    has_pair = co = False
    for a in range(len(anchors)):
        ja, ka = anchors[a]
        for b in range(a + 1, len(anchors)):
            jb, kb = anchors[b]
            if abs(jb - ja) < span:
                continue
            has_pair = True
            seen_a, seen_b = occ[ka], occ[kb]
            where_a = positions.setdefault(ka, group(seen_a))
            both = set()
            for s, p in seen_b:
                if s in where_a and any(abs((p - q) - (jb - ja)) <= WINDOW for q in where_a[s]):
                    both.add(s)
            tally.observed += len(both)
            tally.expected += len(where_a) * len(positions.setdefault(kb, group(seen_b))) / max(nseq, 1)
            if both:
                co = True
    tally.pair += has_pair
    tally.coseen += co


def main():
    prep, seqs_path, meta_path, orfs_path, nreads, out = sys.argv[1:7]
    nreads = int(nreads)
    with open(f"{prep}/models.tsv") as f:
        models = list(csv.DictReader(f, delimiter="\t"))
    index = {r["acc"]: i for i, r in enumerate(models)}
    offs = np.fromfile(f"{prep}/offs.u32", dtype=np.uint32)
    profile = np.memmap(f"{prep}/profile.f32", dtype=np.float32, mode="r").reshape(-1, 21)
    with open(meta_path) as f:
        meta = list(csv.DictReader(f, delimiter="\t"))
    with open(seqs_path) as f:
        seqs = [line.split("\t")[0] for line in f]
    assert len(meta) == len(seqs)
    by_model = defaultdict(list)
    for r, s in zip(meta, seqs):
        m = int(r["model"])
        delta = float("-inf") if float(r["score"]) <= -900 else float(r["score"]) - float(models[m]["ga"])
        strength = next(i for i, low in enumerate((30, 10, 0, -5, -10, -20, float("-inf"))) if delta >= low)
        by_model[m].append((r["set"], strength if delta > float("-inf") else 7, s))
    orfs = []
    with open(orfs_path) as f:
        for line in f:
            if line[0] == ">":
                read = int(line[2:].split("_")[0])
            elif read < nreads:
                orfs.append((read, line.strip()))
    results = {}
    for shape, cols in SHAPES.items():
        span = cols[-1] + 1
        universe = defaultdict(list)  # k-mer -> background ORFs holding it
        for i, (_, s) in enumerate(orfs):
            for j, k in words(s, cols):
                universe[k].append((i, j))
        truth_words = {m: [(st, ga, list(words(s, cols))) for st, ga, s in items] for m, items in by_model.items()}
        wanted = set(universe)
        for items in truth_words.values():
            for _, _, ws in items:
                wanted.update(k for _, k in ws)
        tallies = defaultdict(Tally)
        candidates = defaultdict(lambda: defaultdict(set))
        for acc, family in families(f"{prep}/seed.gz"):
            m = index.get(acc)
            if m is None:
                continue
            occ = defaultdict(list)
            for s, seq in enumerate(family):
                for j, k in words(seq, cols):
                    if k in wanted:
                        occ[k].append((s, j))
            if not occ:
                continue
            lo = np.asarray(profile[offs[m]:offs[m + 1]])
            scores = {}
            score = lambda k: scores.setdefault(k, best_score(lo, k, cols))  # noqa: B023
            nseq = len(family)
            positions = {}
            nb = next(i for i, (a, z) in enumerate(NSEQ_BINS) if a <= nseq < z)
            for st, ga, ws in truth_words.get(m, []):
                hits = [(j, k) for j, k in ws if k in occ]
                for t in THRESHOLDS:
                    anchors = [(j, k) for j, k in hits if score(k) >= t]
                    keys = [("truth", st, ga, t, "all")] + [("truth", st, ga, t, f"nseq{nb}")]
                    for key in keys:
                        tallies[key].n += 1
                    if anchors:
                        for key in keys:
                            evaluate(anchors, occ, nseq, span, tallies[key], positions)
            per_orf = defaultdict(list)
            for k in occ:
                for i, j in universe.get(k, ()):
                    per_orf[i].append((j, k))
            for i, hits in per_orf.items():
                for t in THRESHOLDS:
                    anchors = sorted((j, k) for j, k in hits if score(k) >= t)
                    if not anchors:
                        continue
                    t1 = Tally()
                    evaluate(anchors, occ, nseq, span, t1, positions)
                    for key in (("bg", "reads", 0, t, "all"), ("bg", "reads", 0, t, f"nseq{nb}")):
                        tallies[key].n += 1
                        tallies[key].single += t1.single
                        tallies[key].pair += t1.pair
                        tallies[key].coseen += t1.coseen
                        tallies[key].observed += t1.observed
                        tallies[key].expected += t1.expected
                    read = orfs[i][0]
                    candidates[t]["one"].add((read, m))
                    if t1.pair:
                        candidates[t]["pair"].add((read, m))
                    if t1.coseen:
                        candidates[t]["co_seen"].add((read, m))
        results[shape] = {
            "tallies": [{"key": list(k), **v.row()} for k, v in sorted(tallies.items(), key=str)],
            "candidates_per_read": {str(t): {k: len(v) / nreads for k, v in c.items()} for t, c in candidates.items()},
        }
        with open(out, "w") as o:
            json.dump(results, o)
        print(shape, "done", flush=True)


if __name__ == "__main__":
    main()
