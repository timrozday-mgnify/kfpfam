"""Where pyhmmer's per-model cost goes: one thread, each piece of a model's search timed apart.
Usage: permodel.py H3M MODELS.tsv ORFS.faa FILTER.tsv BAR NMODELS
"""
import os
import random
import resource
import sys
import threading
import time
from collections import defaultdict

import pyhmmer

Z = 20_000_000


def watchdog(gb=6):
    def run():
        while True:
            if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > gb * 2**30:  # bytes on macOS
                print("watchdog: over", gb, "GB", flush=True)
                os._exit(3)
            time.sleep(0.5)
    threading.Thread(target=run, daemon=True).start()


def main():
    h3m, models, faa, tsv, bar, nmodels = sys.argv[1:]
    watchdog()
    accs = [line.split("\t")[0] for line in open(models)][1:]
    cands = defaultdict(set)
    for line in open(tsv):
        q, m, _, g = line.split()
        if float(g) >= float(bar):
            cands[accs[int(m)]].add(q)
    random.seed(1)
    pick = sorted(random.sample(sorted(cands), min(int(nmodels), len(cands))))
    need_q = {q for a in pick for q in cands[a]}
    seqs, name = {}, None
    for line in open(faa):
        if line[0] == ">":
            name = line[1:].split()[0]
        elif name in need_q:
            seqs[name] = seqs.get(name, "") + line.strip()
    alphabet = pyhmmer.easel.Alphabet.amino()
    want = set(pick)
    hmms = {}
    with pyhmmer.plan7.HMMFile(h3m) as f:
        for h in f:
            a = (h.accession.decode() if isinstance(h.accession, bytes) else h.accession).split(".")[0]
            if a in want:
                hmms[a] = h
    bg = pyhmmer.plan7.Background(alphabet)
    pipe = pyhmmer.plan7.Pipeline(alphabet, Z=Z, E=1000, domE=1000)
    nofilt = pyhmmer.plan7.Pipeline(alphabet, Z=Z, E=1000, domE=1000, F1=1.0, F2=1.0, F3=1.0, bias_filter=False)
    one = pyhmmer.easel.TextSequenceBlock([pyhmmer.easel.TextSequence(name=b"x", sequence="ACDEFGHIKLMNPQRSTVWY" * 3)]).digitize(alphabet)
    print("acc\tM\tn\tblock\tprofile\topt\tsearch_hmm\tsearch_opt\tsearch_opt_nofilt\tone_seq_hmm\tone_seq_opt\trows\thits")
    tp = time.perf_counter
    for a in pick:
        hmm = hmms[a]
        t0 = tp()
        block = pyhmmer.easel.TextSequenceBlock(
            pyhmmer.easel.TextSequence(name=q.encode(), sequence=seqs[q]) for q in sorted(cands[a])).digitize(alphabet)
        t1 = tp()
        prof = pyhmmer.plan7.Profile(hmm.M, alphabet)
        prof.configure(hmm, bg, L=400)
        t2 = tp()
        opt = prof.to_optimized()
        t3 = tp()
        hits = pipe.search_hmm(hmm, block)
        t4 = tp()
        n_hits = len(hits)
        pipe.clear()
        t5 = tp()
        pipe.search_hmm(opt, block)
        pipe.clear()
        t6 = tp()
        nofilt.search_hmm(opt, block)
        nofilt.clear()
        t7 = tp()
        pipe.search_hmm(hmm, one)
        pipe.clear()
        t8 = tp()
        pipe.search_hmm(opt, one)
        pipe.clear()
        t9 = tp()
        rows = [(h.name, d.env_from, d.env_to, d.score, d.alignment.hmm_from, d.alignment.target_sequence)
                for h in hits for d in h.domains]
        t10 = tp()
        ms = [x * 1e3 for x in (t1 - t0, t2 - t1, t3 - t2, t4 - t3, t6 - t5, t7 - t6, t8 - t7, t9 - t8, t10 - t9)]
        print(a, hmm.M, len(block), *(f"{x:.4f}" for x in ms), n_hits, sep="\t", flush=True)
        del rows


if __name__ == "__main__":
    main()
