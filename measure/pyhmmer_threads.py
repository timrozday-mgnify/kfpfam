"""handoff.search's pattern against variants, at several thread counts: wall, process CPU, per model / per pair.
Usage: threads.py H3M MODELS.tsv ORFS.faa FILTER.tsv BAR THREADS,...
"""
import sys
import threading
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import pyhmmer

sys.path.insert(0, sys.path[0])
from pyhmmer_permodel import watchdog  # noqa: E402

Z = 20_000_000
alphabet = pyhmmer.easel.Alphabet.amino()


def run(variant, hmms, opts, cands, seqs, blocks, cpus):
    local = threading.local()

    def one(acc):
        if not hasattr(local, "pipe"):
            local.pipe = pyhmmer.plan7.Pipeline(alphabet, Z=Z, E=1000, domE=1000)
        if variant == "current":  # handoff.search / pairs.py today
            block = pyhmmer.easel.TextSequenceBlock(
                pyhmmer.easel.TextSequence(name=q.encode(), sequence=seqs[q]) for q in sorted(cands[acc])).digitize(alphabet)
            hits = local.pipe.search_hmm(hmms[acc], block)
        elif variant == "preblock":
            hits = local.pipe.search_hmm(hmms[acc], blocks[acc])
        else:  # preblock + optimized profile built once
            hits = local.pipe.search_hmm(opts[acc], blocks[acc])
        n = sum(len(h.domains) for h in hits)
        local.pipe.clear()
        return n

    order = sorted(cands)
    if variant.endswith("big1st"):  # longest jobs first, so no thread is left with the tail
        order = sorted(cands, key=lambda a: -len(cands[a]) * hmms[a].M)
    w, c = time.perf_counter(), time.process_time()
    with ThreadPoolExecutor(cpus) as ex:
        n = sum(ex.map(one, order, chunksize=64 if not variant.endswith("big1st") else 1))
    return time.perf_counter() - w, time.process_time() - c, n


def main():
    h3m, models, faa, tsv, bar, threads = sys.argv[1:]
    watchdog()
    accs = [line.split("\t")[0] for line in open(models)][1:]
    cands = defaultdict(set)
    for line in open(tsv):
        q, m, _, g = line.split()
        if float(g) >= float(bar):
            cands[accs[int(m)]].add(q)
    need = {q for v in cands.values() for q in v}
    seqs, name = {}, None
    for line in open(faa):
        if line[0] == ">":
            name = line[1:].split()[0]
        elif name in need:
            seqs[name] = seqs.get(name, "") + line.strip()
    hmms = {}
    with pyhmmer.plan7.HMMFile(h3m) as f:
        for h in f:
            a = (h.accession.decode() if isinstance(h.accession, bytes) else h.accession).split(".")[0]
            if a in cands:
                hmms[a] = h
    t = time.perf_counter()
    bg = pyhmmer.plan7.Background(alphabet)
    opts = {}
    for a, h in hmms.items():
        p = pyhmmer.plan7.Profile(h.M, alphabet)
        p.configure(h, bg, L=400)
        opts[a] = p.to_optimized()
    t_opt = time.perf_counter() - t
    t = time.perf_counter()
    blocks = {a: pyhmmer.easel.TextSequenceBlock(
        pyhmmer.easel.TextSequence(name=q.encode(), sequence=seqs[q]) for q in sorted(v)).digitize(alphabet)
        for a, v in cands.items()}
    t_blk = time.perf_counter() - t
    pairs = sum(len(v) for v in cands.values())
    print(f"{len(cands)} models, {pairs} pairs; optimized profiles {t_opt:.1f} s, blocks {t_blk:.1f} s (one thread)", flush=True)
    for cpus in map(int, threads.split(",")):
        for v in ("current", "preblock", "preopt", "preopt_big1st"):
            wall, cpu, n = run(v, hmms, opts, cands, seqs, blocks, cpus)
            print(f"{v:14s} {cpus:2d} threads: wall {wall:6.2f} s, cpu {cpu:6.2f} s, wall*threads {wall * cpus:6.2f} core-s; "
                  f"{wall * cpus / len(cands) * 1e3:.3f} core-ms per model or {wall * cpus / pairs * 1e6:.1f} core-us per pair; {n} domains",
                  flush=True)


if __name__ == "__main__":
    main()
