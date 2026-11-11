"""Speed plan step 6: pyhmmer time on the filter's candidates at scale (no truth at this depth).
Usage: handoff_time.py PFAM.hmm MODELS.tsv ORFS.faa FILTER.tsv NREADS CPUS BAR...
"""
import sys
import time

import pyhmmer

from handoff import candidates, fasta, search

hmm_path, models, faa, tsv, nreads, cpus, *bars = sys.argv[1:]
accs = [line.split("\t")[0] for line in open(models)][1:]
cands = {float(b): candidates(tsv, float(b)) for b in bars}
need = {accs[m] for c in cands.values() for m in c}
t = time.time()
hmms = {}
with pyhmmer.plan7.HMMFile(hmm_path) as f:
    for hmm in f:
        acc = (hmm.accession.decode() if isinstance(hmm.accession, bytes) else hmm.accession).split(".")[0]
        if acc in need:
            hmms[acc] = hmm
print(f"{len(hmms)} HMMs read in {time.time() - t:.0f} s", flush=True)
seqs = fasta(faa, {q for c in cands.values() for v in c.values() for q in v})
import os
ga = os.environ.get("CUT") == "ga"
for bar, c in sorted(cands.items()):
    found, secs = search(hmms, accs, c, seqs, int(cpus), ga)
    n = sum(len(v) for v in c.values())
    print(f"{'--cut_ga' if ga else 'E <= 1000'}, gapped >= {bar}: {n / int(nreads):.2f} candidates per read, {len(c)} models, {len(found)} reported; "
          f"{secs:.1f} s wall on {cpus} cpus, {secs * int(cpus) / int(nreads) * 1e3:.3f} core-ms per read, {secs * int(cpus) / n * 1e6:.1f} core-us per pair", flush=True)
