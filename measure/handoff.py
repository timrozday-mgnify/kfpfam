"""Speed plan step 5: HMMER on the filter's candidate pairs only, against hmmsearch on every pair.

pyhmmer per model on its candidate ORFs, Z fixed at 20 million and E <= 1000 as the truth was run, so
scores and E-values are comparable. Reports, per bar on the filter's gapped score: candidates and
reversed-read candidates per read, pyhmmer time, recall by bin of hmmsearch's pairs, recall at GA of
pairs pyhmmer calls at GA, and whether shared pairs' best-domain scores agree.

Usage (in the pyhmmer container): handoff.py PFAM.hmm MODELS.tsv REAL.faa REAL.tsv DECOY.faa DECOY.tsv
                                  TRUTH.domtbl DECOY_TRUTH.domtbl NREADS CPUS BAR...
"""

import os
import sys
import time
from collections import defaultdict

import pyhmmer

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "recruit"))
from pairs import search_each  # noqa: E402

LOWS = (0, -5, -10, -20, float("-inf"))
NAMES = (">=GA", "-5-GA", "-10--5", "-20--10", "<-20")
Z = 20_000_000


def fasta(path, keep):
    out, name = {}, None
    for line in open(path):
        if line[0] == ">":
            name = line[1:].split()[0]
        elif name in keep:
            out[name] = line.strip()
    return out


def candidates(path, bar):
    """{model index: [orf]} with a gapped score >= bar."""
    out = defaultdict(set)
    for line in open(path):
        q, m, _, g = line.split()
        if float(g) >= bar:
            out[int(m)].add(q)
    return out


def truth(path):
    best = {}
    for line in open(path):
        if line[0] != "#":
            c = line.split()
            k = (c[0], c[4].split(".")[0])
            best[k] = max(best.get(k, -1e9), float(c[13]))
    return best


def search(hmms, accs, cands, seqs, cpus, ga=False):
    """{(orf, acc): best domain score} from pyhmmer on each model's candidates (recruit/pairs.py's search_each)."""
    def keep(hmm, hits):
        acc = (hmm.accession.decode() if isinstance(hmm.accession, bytes) else hmm.accession).split(".")[0]
        return [((h.name.decode() if isinstance(h.name, bytes) else h.name, acc), max(d.score for d in h.domains))
                for h in hits if h.domains]

    t = time.time()
    out = dict(search_each({m: hmms[accs[m]] for m in cands}, cands, seqs, cpus, keep,
                           **({"bit_cutoffs": "gathering"} if ga else {})))
    return out, time.time() - t


def main():
    hmm_path, models, real_faa, real_tsv, decoy_faa, decoy_tsv, tt, dt, nreads, cpus, *bars = sys.argv[1:]
    nreads, cpus = int(nreads), int(cpus)
    rows = [line.rstrip("\n").split("\t") for line in open(models)][1:]
    accs = [r[0] for r in rows]
    ga = {r[0]: float(r[4]) for r in rows}
    lowest = min(float(b) for b in bars)
    real_c, decoy_c = candidates(real_tsv, lowest), candidates(decoy_tsv, lowest)
    need = {accs[m] for m in list(real_c) + list(decoy_c)}
    t = time.time()
    hmms = {}
    with pyhmmer.plan7.HMMFile(hmm_path) as f:
        for hmm in f:
            acc = (hmm.accession.decode() if isinstance(hmm.accession, bytes) else hmm.accession).split(".")[0]
            if acc in need:
                hmms[acc] = hmm
    print(f"{len(hmms)} HMMs read in {time.time() - t:.0f} s", flush=True)
    real_s = fasta(real_faa, {q for v in real_c.values() for q in v})
    decoy_s = fasta(decoy_faa, {q for v in decoy_c.values() for q in v})
    real_t, decoy_t = truth(tt), truth(dt)
    bins = defaultdict(set)
    for k, s in real_t.items():
        if k[1] in ga:
            bins[next(i for i, low in enumerate(LOWS) if s - ga[k[1]] >= low)].add(k)
    for bar in sorted(float(b) for b in bars):
        rc = {m: {q for q in v if True} for m, v in candidates(real_tsv, bar).items()}
        dc = candidates(decoy_tsv, bar)
        found, secs = search(hmms, accs, rc, real_s, cpus)
        found_d, secs_d = search(hmms, accs, dc, decoy_s, cpus)
        n_c = sum(len(v) for v in rc.values())
        n_d = sum(len(v) for v in dc.values())
        # recall: hmmsearch's pairs that pyhmmer also reports, by bin; at GA: called at GA by both
        rec = [sum(k in found for k in bins[i]) / max(len(bins[i]), 1) for i in range(5)]
        ga_true = {k for k, s in real_t.items() if k[1] in ga and s >= ga[k[1]]}
        ga_found = {k for k, s in found.items() if s >= ga[k[1]]}
        shared = [k for k in found if k in real_t]
        diff = sorted(abs(found[k] - real_t[k]) for k in shared)
        extra = len([k for k in found if k not in real_t])
        missed = sorted(ga_true - ga_found)
        print("  GA pairs not called: " + ", ".join(f"{q} {a} truth {real_t[q, a]:.1f} found {found.get((q, a), float('nan')):.1f} ga {ga[a]}" for q, a in missed[:8]))
        print(f"\ngapped >= {bar}: candidates {n_c / nreads:.2f} per read, reversed {n_d / nreads:.2f}; "
              f"pyhmmer {secs:.1f} + {secs_d:.1f} s wall on {cpus} cpus ({(secs + secs_d) * cpus / nreads * 1e3:.3f} core-ms per read)")
        print("  recall of hmmsearch's pairs by bin " + " ".join(f"{n} {r:.3f}" for n, r in zip(NAMES, rec)))
        print(f"  GA calls: {len(ga_found & ga_true)} of {len(ga_true)} ({len(ga_found & ga_true) / len(ga_true):.4f}), "
              f"{len(ga_found - ga_true)} not called at GA by hmmsearch; reversed reads called at GA: "
              f"{sum(s >= ga[k[1]] for k, s in found_d.items())} (hmmsearch {sum(s >= ga[k[1]] for k, s in decoy_t.items() if k[1] in ga)})")
        print(f"  best-domain scores on {len(shared)} shared pairs: within 0.1 bits {sum(d <= 0.1 for d in diff) / max(len(diff), 1):.4f}, "
              f"largest {diff[-1] if diff else 0:.2f}; {extra} pairs hmmsearch didn't report")


if __name__ == "__main__":
    main()
