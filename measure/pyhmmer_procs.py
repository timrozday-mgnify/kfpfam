"""End to end: fork workers build each model's block and optimized profile themselves; pool start included.
Usage: e2e.py H3M MODELS.tsv ORFS.faa FILTER.tsv BAR CPUS"""
import multiprocessing as mp
import sys
import time
from collections import defaultdict

import pyhmmer

sys.path.insert(0, sys.path[0])
from pyhmmer_permodel import watchdog  # noqa: E402

A = pyhmmer.easel.Alphabet.amino()
G = {}


def work(part):
    c = time.process_time()
    pipe = pyhmmer.plan7.Pipeline(A, Z=20_000_000, E=1000, domE=1000)
    bg = pyhmmer.plan7.Background(A)
    out = []
    for a in part:
        h = G["hmms"][a]
        p = pyhmmer.plan7.Profile(h.M, A)
        p.configure(h, bg, L=400)
        block = pyhmmer.easel.TextSequenceBlock(
            pyhmmer.easel.TextSequence(name=q.encode(), sequence=G["seqs"][q]) for q in sorted(G["cands"][a])).digitize(A)
        for hit in pipe.search_hmm(p.to_optimized(), block):
            for d in hit.domains:
                out.append((hit.name, a, d.env_from, d.env_to, d.score))
        pipe.clear()
    return out, time.process_time() - c


def main():
    h3m, models, faa, tsv, bar, cpus = sys.argv[1:]
    cpus = int(cpus)
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
    G.update(hmms=hmms, seqs=seqs, cands=cands)
    order = sorted(cands, key=lambda a: -len(cands[a]) * hmms[a].M)
    parts = [order[i::cpus * 8] for i in range(cpus * 8)]
    w = time.perf_counter()
    with mp.get_context("fork").Pool(cpus) as pool:
        res = list(pool.imap_unordered(work, parts))
    wall = time.perf_counter() - w
    n, cpu = sum(len(r[0]) for r in res), sum(r[1] for r in res)
    pairs = sum(len(v) for v in cands.values())
    print(f"{cpus} procs: wall {wall:.2f} s, worker CPU {cpu:.2f} s, wall*procs {wall * cpus:.2f} core-s, "
          f"{wall * cpus / pairs * 1e6:.1f} core-us per pair, {wall * cpus / 50000 * 1e3:.3f} core-ms per read; {n} domains", flush=True)


if __name__ == "__main__":
    main()
