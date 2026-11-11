"""Promiscuous-word lists on top of per-profile seed shapes.

Policies: k4 T7 everywhere (the earlier measurement's base), and per_k.py's per-model choice at mu 1e-3
(about 1x index) and 1e-4 (about 5x). A word is the read's residues at a shape's offsets. Rules, each
removing words until a share of the policy's hits is gone:
    sample    the words with most hits on odd reads, applied to even reads (a list from one block, used
              on the next)
    index     the words with the longest posting lists under the policy (the index alone)
    expected  posting length x the word's background frequency (the index and a fixed composition)
Hits removed are counted on even reads; recall on hmmsearch pairs, up to MAX_PER_BIN per bin: a pair
survives when a seed on a diagonal >= 22 bits has a word not on the list.

per_k_promisc.py policies PREP OUT
per_k_promisc.py scan PREP ORFS.faa N_READS WORKERS OUT
per_k_promisc.py pairs PREP ORFS.faa DOMTBL MAX_PER_BIN OUT
per_k_promisc.py postings PREP OUT WORKERS       posting lengths of the pairs' seed words
per_k_promisc.py report OUT
"""

import multiprocessing as mp
import pickle
import random
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

import gor
import per_k
from per_k import NB, SHAPES, offsets, span, ts
from spaced import blocks, diag_blocks, pair_scores, pattern_scores
from ss_kmers import BAR, BIN_LOW, BIN_NAMES, Q, diagonals, domtbl, fasta, load_prep

S = list(SHAPES)
POLICIES = {"k4 T7": None, "per model mu 1e-3": 1e-3, "per model mu 1e-4": 1e-4}


def words(x, s):
    """Each read position's word for shape s; -1 where it runs off or holds an unknown residue."""
    off = offsets(s)
    sp = span(s)
    out = np.full(len(x), -1, dtype=np.int64)
    if len(x) < sp:
        return out
    w = np.stack([x[o:len(x) - sp + 1 + o] for o in off], 1).astype(np.int64)
    c = w @ (20 ** np.arange(len(off)))
    c[(w == 20).any(1)] = -1
    out[:len(c)] = c
    return out


def policies(prep, out):
    out = Path(out)
    no = pickle.load(open(out / "per_k_noise.pkl", "rb"))
    pr = pickle.load(open(out / "per_k_profile.pkl", "rb"))
    opts = [(s, t) for s in S for t in range(NB)]
    h1 = np.stack([no["halves"][s][1][:, ::-1].cumsum(1)[:, ::-1][:, t] for s, t in opts], 1)
    ent = np.stack([pr["entries"][s][:, t] for s, t in opts], 1)
    base = opts.index(("k4", 0))
    e0, M = ent[:, base].sum(), len(ent)
    det = np.stack([pr["det"][s, 1.0][:, t] for s, t in opts], 1)
    ok = det >= det[:, [base]] - 1e-12
    ok[:, base] = True
    pol = {}
    for name, mu in POLICIES.items():
        if mu is None:
            choice = np.full(M, base)
        else:
            choice = np.where(ok, h1 / h1[:, base].sum() + mu * ent / e0 * M / 100, np.inf).argmin(1)
        shape = np.array([S.index(opts[c][0]) for c in choice])
        thr = np.array([ts(opts[c][0])[opts[c][1]] for c in choice])
        pol[name] = (shape, thr)
        print(name, "entries", f"{ent[np.arange(M), choice].sum() / e0:.2f}", np.bincount(shape, minlength=len(S)))
    pickle.dump(pol, open(out / "promisc_policies.pkl", "wb"))


THR = None  # {policy: {shape index: per-column T, inf where the column's model uses another shape}}


def init(prep, out):
    global THR
    per_k.init(prep)
    pol = pickle.load(open(Path(out) / "promisc_policies.pkl", "rb"))
    THR = {}
    for name, (shape, thr) in pol.items():
        THR[name] = {}
        for si in np.unique(shape):
            t = np.where(shape == si, thr, np.inf).astype(np.float32)
            THR[name][si] = t[per_k.MODEL]


def scan_orf(seq):
    """{policy: (unique (shape, word) keys, counts)} of hits, and word occurrences by shape."""
    x = gor.encode(seq).astype(np.int64)
    by_res, to_end = per_k.BY_RES, per_k.TO_END
    n = by_res.shape[1]
    pad = max(span(s) for s in S)
    found = {p: [] for p in THR}
    wd = {si: words(x, s) for si, s in enumerate(S)}
    for c0 in range(0, n, per_k.CHUNK):
        c1 = min(c0 + per_k.CHUNK, n)
        blk = blocks(by_res[:, c0:min(c1 + pad, n)][x])
        for si, s in enumerate(S):
            users = [p for p in THR if si in THR[p]]
            if not users:
                continue
            a, g, b_ = SHAPES[s]
            v = pattern_scores(blk, (a, g, b_, 0))[:, :c1 - c0]
            for p in users:
                i, c = np.nonzero(v >= THR[p][si][c0:c1])
                c = c + c0
                keep = to_end[c] >= span(s)
                w = wd[si][i[keep]]
                w = w[w >= 0]
                found[p].append(si * 20 ** 7 + w)
    res = {}
    for p, parts in found.items():
        k = np.concatenate(parts) if parts else np.zeros(0, dtype=np.int64)
        res[p] = np.unique(k, return_counts=True)
    occ = np.concatenate([si * 20 ** 7 + w[w >= 0] for si, w in wd.items()])
    return res, np.unique(occ, return_counts=True)


def merge(acc, u, c):
    acc.append((u, c))
    if len(acc) > 50:
        u = np.concatenate([a for a, _ in acc])
        c = np.concatenate([b for _, b in acc])
        uu, inv = np.unique(u, return_inverse=True)
        acc[:] = [(uu, np.bincount(inv, weights=c).astype(np.int64))]


def collapse(acc):
    if not acc:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    merge(acc, *acc.pop())
    u = np.concatenate([a for a, _ in acc])
    c = np.concatenate([b for _, b in acc])
    uu, inv = np.unique(u, return_inverse=True)
    return uu, np.bincount(inv, weights=c).astype(np.int64)


def scan(prep, orfs, reads, workers, out):
    named = fasta(orfs, int(reads))
    seqs = list(named.values())
    half = [int(nm[1:].split("_")[0]) % 2 for nm in named]
    init(prep, out)
    hits = {(p, h): [] for p in THR for h in (0, 1)}
    occ = {h: [] for h in (0, 1)}
    with mp.get_context("fork").Pool(int(workers)) as pool:
        for n, (res, (ou, oc)) in enumerate(pool.imap(scan_orf, seqs, chunksize=2)):
            for p, (u, c) in res.items():
                merge(hits[p, half[n]], u, c)
            merge(occ[half[n]], ou, oc)
            if n % 50 == 0:
                print(n, len(seqs), flush=True)
    reads_by_half = [len({nm.split("_")[0] for nm, q in zip(named, half) if q == h}) for h in (0, 1)]
    data = dict(hits={k: collapse(v) for k, v in hits.items()}, occ={h: collapse(v) for h, v in occ.items()},
                reads=reads_by_half)
    pickle.dump(data, open(Path(out) / "promisc_scan.pkl", "wb"))
    print("scan done", {p: data["hits"][p, 0][1].sum() / reads_by_half[0] for p in THR})


def pairs(prep, orfs, table, per_bin, out):
    lo, offs = load_prep(prep)
    pol = pickle.load(open(Path(out) / "promisc_policies.pkl", "rb"))
    rows = [line.rstrip("\n").split("\t") for line in open(Path(prep) / "models.tsv")][1:]
    index = {r[0]: i for i, r in enumerate(rows)}
    ga = {r[0]: float(r[4]) for r in rows}
    seqs = fasta(orfs)
    best = {}
    for orf, acc, score, *_ in domtbl(table):
        if acc in index and orf in seqs:
            best[orf, acc] = max(score, best.get((orf, acc), -1e9))
    bins = defaultdict(list)
    for key, sc in best.items():
        bins[next(i for i, low in enumerate(BIN_LOW) if sc - ga[key[1]] >= low)].append(key)
    rng = random.Random(23)
    picks = [(key, b) for b in sorted(bins) for key in rng.sample(sorted(bins[b]), min(int(per_bin), len(bins[b])))]
    seeds = {p: [] for p in pol}
    bin_of = np.array([b for _, b in picks])
    for r, ((orf, acc), b) in enumerate(picks):
        m = index[acc]
        x = gor.encode(seqs[orf]).astype(np.int64)
        a_mat, diag, _ = diagonals(np.asarray(lo[offs[m]:offs[m + 1]]), x)
        good = np.flatnonzero(diag >= BAR)
        blk = diag_blocks(a_mat) if len(good) else None
        for p, (shape, thr) in pol.items():
            if blk is None:
                seeds[p].append(np.zeros(0, dtype=np.int64))
                continue
            si = shape[m]
            a, g, b_ = SHAPES[S[si]]
            v = pair_scores(blk, (a, g, b_, 0), good)
            pos = np.flatnonzero((v >= thr[m]).any(0))
            w = words(x, S[si])[pos]
            seeds[p].append(np.unique(si * 20 ** 7 + w[w >= 0]))
        if r % 2000 == 0:
            print(r, len(picks), flush=True)
    pickle.dump(dict(seeds=seeds, bin=bin_of), open(Path(out) / "promisc_pairs.pkl", "wb"))
    print("pairs done", len(picks))


def posting_one(job):
    p, key = job
    by_res, to_end = per_k.BY_RES, per_k.TO_END
    n = by_res.shape[1]
    si, w = divmod(int(key), 20 ** 7)
    off = offsets(S[si])
    v = np.zeros(n, dtype=np.float32)
    for k, o in enumerate(off):
        v[:n - o] += by_res[(w // 20 ** k) % 20, o:]
    return int(((v >= THR[p][si]) & (to_end >= span(S[si]))).sum())


def postings(prep, out, workers):
    """Posting length under each policy of every seed word in the pairs, by scoring it on every column."""
    out = Path(out)
    init(prep, out)
    pa = pickle.load(open(out / "promisc_pairs.pkl", "rb"))
    post = {}
    with mp.get_context("fork").Pool(int(workers)) as pool:
        for p, lists in pa["seeds"].items():
            allw = np.unique(np.concatenate(lists))
            lens = np.array(pool.map(posting_one, [(p, k) for k in allw], chunksize=64), dtype=np.int64)
            post[p] = (allw, lens)
            print(p, len(allw), "words", flush=True)
    pickle.dump(post, open(out / "promisc_postings.pkl", "wb"))


def report(out):
    out = Path(out)
    sc = pickle.load(open(out / "promisc_scan.pkl", "rb"))
    pa = pickle.load(open(out / "promisc_pairs.pkl", "rb"))
    post = pickle.load(open(out / "promisc_postings.pkl", "rb"))
    b = pa["bin"]
    logq = np.log(Q)

    def bg(keys):
        si, w = np.divmod(keys, 20 ** 7)
        out_ = np.zeros(len(keys))
        for k in range(7):
            aa = (w // 20 ** k) % 20
            use = k < np.array([len(offsets(s)) for s in S])[si]
            out_ += np.where(use, logq[aa], 0)
        return out_

    print("recall by bin " + " ".join(BIN_NAMES) + "; hits removed on even reads, by a list made as described")
    for p in pa["seeds"]:
        u0, c0 = sc["hits"][p, 1]  # odd reads: builds the sample list
        u1, c1 = sc["hits"][p, 0]  # even reads: scored
        ou, oc = (np.concatenate(z) for z in zip(sc["occ"][0], sc["occ"][1]))
        ou, inv = np.unique(ou, return_inverse=True)
        oc = np.bincount(inv, weights=oc)
        # posting length of a hit word: its hits per occurrence over both halves
        hu = np.concatenate([u0, u1])
        hc = np.concatenate([c0, c1])
        hu, inv = np.unique(hu, return_inverse=True)
        hc = np.bincount(inv, weights=hc)
        plen = hc / oc[np.searchsorted(ou, hu)]
        pw, pl = post[p]
        seeds = pa["seeds"][p]
        total = c1.sum()
        print(f"\n{p}: {total / sc['reads'][0]:,.0f} hits per read (even reads)")

        def show(label, banned):
            removed = c1[np.isin(u1, banned)].sum() / total
            ok = np.array([len(s) and not np.isin(s, banned).all() for s in seeds])
            r = [ok[b == i].mean() for i in range(5)]
            print(f"  {label:22} removed {removed:.3f}  " + " ".join(f"{v:.3f}" for v in r))

        show("none", np.zeros(0, dtype=np.int64))
        order_s = u0[np.argsort(-c0)]
        cum_s = np.sort(c0)[::-1].cumsum() / c0.sum()
        # index and expected: rank every word with a known posting length (hit words and seed words)
        known = np.unique(np.concatenate([hu, pw]))
        lens = np.zeros(len(known))
        lens[np.searchsorted(known, hu)] = plen
        lens[np.searchsorted(known, pw)] = pl
        hits_known = np.zeros(len(known))
        hits_known[np.searchsorted(known, u0)] = c0
        for rule, score in (("index", lens), ("expected", np.log(np.maximum(lens, 1e-9)) + bg(known))):
            order = np.argsort(-score)
            cum = hits_known[order].cumsum() / c0.sum()
            for share in (0.05, 0.1, 0.2, 0.3, 0.4):
                show(f"{rule} {share}", known[order[:np.searchsorted(cum, share) + 1]])
        for share in (0.05, 0.1, 0.2, 0.3, 0.4):
            show(f"sample {share}", order_s[:np.searchsorted(cum_s, share) + 1])


if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    {"policies": policies, "scan": scan, "pairs": pairs, "postings": postings, "report": report}[cmd](*args)
