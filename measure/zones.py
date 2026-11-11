"""Can a read's own sequence say where its true seeds are, so the other positions skip their lookups?

Labels come from PDB chains aligned to Pfam at GA: a residue is positive if the 4-mer around it scores
>= T at its true column. Two predictors, each trained without the clan of the model it is tested on:
    gor    GOR scoring on a 17-residue window, states: seed, in a domain without a seed, outside
    kmer4  a table of 160,000 4-mers, log-odds of being a seed
Tested on real reads: skipping positions below a score removes that share of the k4 T7 hits (from a
scan against every column) and loses the hmmsearch hits with no seed left on a diagonal >= 22 bits.

The value rule scores a 4-mer by its seed log-odds minus the log of its hits (its posting length in
the k4 T7 index), so a k-mer pays for its noise. Rarity alone is the control.

Usage: zones.py SS_DIR SS_OUT PREP ORFS.faa DOMTBL OUT N_READS WORKERS MAX_PER_BIN
       zones.py value SS_DIR PREP ORFS.faa DOMTBL k4T7.idx MAX_PER_BIN
"""

import multiprocessing as mp
import pickle
import random
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

import gor
from ss_accuracy import FOLDS, alignments, chain_folds, model_folds, read_chains
from ss_kmers import BAR, BIN_LOW, BIN_NAMES, diagonals, domtbl, fasta, load_prep

K, T = 4, 7.0
STARTS, STARTS_BYTES = np.uint64, 8  # filter's index: k and total as u64, then 20^k + 1 starts
QUANTILES = (0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.98)


def codes(x):
    """Each position's 4-mer as an integer; -1 where it runs off the end or holds an unknown residue."""
    if len(x) < K:
        return np.full(len(x), -1, dtype=np.int64)
    w = np.lib.stride_tricks.sliding_window_view(x.astype(np.int64), K)
    c = w @ (20 ** np.arange(K - 1, -1, -1))
    c[(w == 20).any(1)] = -1
    return np.concatenate([c, np.full(K - 1, -1, dtype=np.int64)])


def labels(ss, prep, chains, index, offs):
    """{chain: L x state}: 0 seed at this position's 4-mer, 1 aligned without one, 2 outside a domain."""
    lo = np.memmap(prep / "profile.f32", dtype=np.float32, mode="r").reshape(-1, 21)
    out = {n: np.full(len(x), 2, dtype=np.int64) for n, (x, _) in chains.items()}
    for acc, chain, pairs in alignments(ss / "reps.hs.txt"):
        if chain not in chains or acc not in index:
            continue
        x = chains[chain][0]
        base = offs[index[acc]] - 1
        y = out[chain]
        for j, (col, res) in enumerate(pairs):
            if not col or res < 0 or res >= len(x):
                continue
            if y[res] == 2:
                y[res] = 1
            run = pairs[j:j + K]
            if len(run) == K and all(c == col + d and r == res + d for d, (c, r) in enumerate(run)):
                if res + K <= len(x) and x[res:res + K].max() < 20:
                    if sum(lo[base + col + d, x[res + d]] for d in range(K)) >= T:
                        y[res + 1] = 0  # the 4-mer's centre, as GOR scores a centre
    return out


def train(chains, y, homes):
    """Per fold: GOR tables, and the 4-mer log-odds table."""
    models = []
    for f in range(FOLDS):
        names = [n for n in chains if f not in homes[n]]
        m = gor.Gor([(chains[n][0], y[n]) for n in names])
        seed, every = np.zeros(20 ** K), np.zeros(20 ** K)
        for n in names:
            c = codes(chains[n][0])
            ok = c >= 0
            every += np.bincount(c[ok], minlength=20 ** K)
            # The label sits on the 4-mer's centre, one past its start
            pos = np.flatnonzero(ok[:-1] & (y[n][1:] == 0)) if len(c) > 1 else []
            if len(pos):
                seed += np.bincount(c[pos], minlength=20 ** K)
        rate = seed.sum() / every.sum()
        table = np.log((seed + 10 * rate) / (every - seed + 10 * (1 - rate)))
        models.append((m, table))
        print(f"fold {f}: {len(names)} chains, seed share {rate:.3f}", flush=True)
    return models


def position_scores(x, model):
    """Both predictors' scores at each 4-mer start. Positions with no 4-mer get -inf."""
    m, table = model
    g = m.score_seq(x.astype(np.int8))[:, 0]
    c = codes(x)
    start_gor = np.full(len(x), -np.inf)
    start_tab = np.full(len(x), -np.inf)
    if len(x) >= K:
        start_gor[:len(x) - K + 1] = g[1:len(x) - K + 2]
    ok = c >= 0
    start_tab[ok] = table[c[ok]]
    return start_gor, start_tab


BY_RES = OFFS = TO_END = MODEL0 = None


def init(prep, model):
    global BY_RES, OFFS, TO_END, MODEL0
    lo, OFFS = load_prep(prep)
    BY_RES = np.array(lo).T.copy()
    TO_END = np.zeros(len(lo), dtype=np.int32)
    for a, b in zip(OFFS[:-1], OFFS[1:]):
        TO_END[a:b] = np.arange(b - a, 0, -1)
    MODEL0 = model


def scan(seq):
    """Hits of each 4-mer at T, against every column, with both predictors' scores."""
    x = gor.encode(seq).astype(np.int64)
    hits = np.zeros(len(x), dtype=np.int64)
    n = BY_RES.shape[1]
    for i in range(len(x) - K + 1):
        cur = BY_RES[x[i]].copy()
        for j in range(1, K):
            cur[:n - j] += BY_RES[x[i + j]][j:]
        idx = np.flatnonzero(cur >= T)
        hits[i] = (TO_END[idx] >= K).sum()
    return hits, *position_scores(x, MODEL0)


def value():
    ss, prep = Path(sys.argv[2]), Path(sys.argv[3])
    orfs, table, idx = sys.argv[4:7]
    per_bin = int(sys.argv[7])
    accs, fold = model_folds(prep)
    index = {a: i for i, a in enumerate(accs)}
    _, offs = load_prep(prep)
    chains = read_chains(ss)
    homes = chain_folds(ss, chains, fold)
    models = train(chains, labels(ss, prep, chains, index, offs), homes)
    with open(idx, "rb") as f:
        k, total = np.frombuffer(f.read(16), dtype=np.uint64)
        starts = np.frombuffer(f.read(STARTS_BYTES * (20 ** int(k) + 1)), dtype=STARTS)
    assert int(k) == K and int(starts[-1]) == int(total)
    post = np.diff(starts).astype(np.int64)
    print("entries", total, "4-mers with none", (post == 0).sum(), flush=True)
    all_seqs = fasta(orfs)
    codes_all = np.concatenate([codes(gor.encode(s)) for s in all_seqs.values()])
    codes_all = codes_all[codes_all >= 0]
    weight = np.bincount(codes_all, minlength=20 ** K) * post  # hits by 4-mer over the whole read set
    reads = len({n.split("_")[0] for n in all_seqs})
    print(f"{weight.sum() / reads:,.0f} hits per read from postings, {reads} reads", flush=True)
    rows = [line.rstrip("\n").split("\t") for line in open(prep / "models.tsv")][1:]
    ga = {r[0]: float(r[4]) for r in rows}
    best = {}
    for orf, acc, score, *_ in domtbl(table):
        if acc in index and orf in all_seqs:
            best[orf, acc] = max(score, best.get((orf, acc), -1e9))
    bins = defaultdict(list)
    for key, sc in best.items():
        bins[next(i for i, low in enumerate(BIN_LOW) if sc - ga[key[1]] >= low)].append(key)
    rng = random.Random(23)
    picks = [(key, b) for b in sorted(bins) for key in rng.sample(sorted(bins[b]), min(per_bin, len(bins[b])))]
    half = len(picks)
    lo = np.memmap(prep / "profile.f32", dtype=np.float32, mode="r").reshape(-1, 21)
    logpost = np.log(np.maximum(post, 1))
    rules = {"value": lambda t: t - logpost, "rarity": lambda t: -logpost, "seed rate": lambda t: t}
    seen = {n: np.full(half, -np.inf) for n in rules}
    bin_of = np.zeros(half, dtype=np.int8)
    for r, ((orf, acc), b) in enumerate(picks):
        bin_of[r] = b
        m = index[acc]
        x = gor.encode(all_seqs[orf]).astype(np.int64)
        a, diag, _ = diagonals(np.asarray(lo[offs[m]:offs[m + 1]]), x)
        good = np.flatnonzero(diag >= BAR)
        if not len(good) or len(x) < K:
            continue
        cum = sum(a[good][:, j:len(x) - K + 1 + j] for j in range(K))
        positions = np.flatnonzero((cum >= T).any(0))
        c = codes(x)[positions]
        c = c[c >= 0]
        if not len(c):
            continue
        for n, fn in rules.items():
            seen[n][r] = fn(models[fold[acc]][1])[c].max()
    print("\nrule hits-kept 4-mers-kept recall " + " ".join(BIN_NAMES))
    table0 = models[0][1]
    for n, fn in rules.items():
        score = fn(table0)
        order = np.argsort(-score)
        cum = weight[order].cumsum() / weight.sum()
        for share in (0.8, 0.6, 0.4, 0.2, 0.1, 0.05, 0.04, 0.02):
            i = int(np.searchsorted(cum, share))
            cut = score[order[i]]
            print(f"{n:9} {share:.2f} {(i + 1) / 20 ** K:.3f}  " + " ".join(f"{(seen[n][bin_of == b] >= cut).mean():.3f}" for b in range(5)), flush=True)


def main():
    ss, out, prep = (Path(a) for a in sys.argv[1:4])
    orfs, table = sys.argv[4:6]
    dest = Path(sys.argv[6])
    reads, workers, per_bin = (int(a) for a in sys.argv[7:10])
    accs, fold = model_folds(prep)
    index = {a: i for i, a in enumerate(accs)}
    _, offs = load_prep(prep)
    chains = read_chains(ss)
    homes = chain_folds(ss, chains, fold)
    y = labels(ss, prep, chains, index, offs)
    print("labels:", np.bincount(np.concatenate(list(y.values())), minlength=3), flush=True)
    models = train(chains, y, homes)

    # Noise: every k4 T7 hit of a sample of reads, by position
    seqs = list(fasta(orfs, reads).values())
    random.Random(1).shuffle(seqs)
    init(prep, models[0])
    hit, sg, st = [], [], []
    with mp.get_context("fork").Pool(workers) as pool:
        for n, (h, a, b) in enumerate(pool.imap_unordered(scan, seqs, chunksize=4)):
            hit.append(h)
            sg.append(a)
            st.append(b)
            if n % 200 == 0:
                print(n, len(seqs), flush=True)
    hit, sg, st = np.concatenate(hit), np.concatenate(sg), np.concatenate(st)
    print(f"{hit.sum() / reads:,.0f} hits per read, {len(hit)} positions", flush=True)

    # Hits of hmmsearch: the best predictor score among positions holding a seed on a good diagonal
    rows = [line.rstrip("\n").split("\t") for line in open(prep / "models.tsv")][1:]
    ga = {r[0]: float(r[4]) for r in rows}
    all_seqs = fasta(orfs)
    best = {}
    for orf, acc, score, *_ in domtbl(table):
        if acc in index and orf in all_seqs:
            best[orf, acc] = max(score, best.get((orf, acc), -1e9))
    bins = defaultdict(list)
    for key, s in best.items():
        bins[next(i for i, low in enumerate(BIN_LOW) if s - ga[key[1]] >= low)].append(key)
    rng = random.Random(23)
    picks = [(key, b) for b in sorted(bins) for key in rng.sample(sorted(bins[b]), min(per_bin, len(bins[b])))]
    lo = np.memmap(prep / "profile.f32", dtype=np.float32, mode="r").reshape(-1, 21)
    pg = np.full(len(picks), -np.inf)
    pt = np.full(len(picks), -np.inf)
    bin_of = np.zeros(len(picks), dtype=np.int8)
    for r, ((orf, acc), b) in enumerate(picks):
        bin_of[r] = b
        m = index[acc]
        x = gor.encode(all_seqs[orf]).astype(np.int64)
        a, diag, _ = diagonals(np.asarray(lo[offs[m]:offs[m + 1]]), x)
        good = np.flatnonzero(diag >= BAR)
        if not len(good) or len(x) < K:
            continue
        cum = sum(a[good][:, j:len(x) - K + 1 + j] for j in range(K))
        positions = np.flatnonzero((cum >= T).any(0))
        if not len(positions):
            continue
        g, t = position_scores(x, models[fold[acc]])
        pg[r], pt[r] = g[positions].max(), t[positions].max()
    print("\npredictor keep-share-of-positions hits-share recall " + " ".join(BIN_NAMES))
    results = {}
    for name, noise_score, pair_score in (("gor", sg, pg), ("kmer4", st, pt)):
        valid = np.isfinite(noise_score)
        print(f"{name:6} all      1.000  " + " ".join(f"{np.isfinite(pair_score[bin_of == i]).mean():.3f}" for i in range(5)))
        for q in QUANTILES:
            cut = np.quantile(noise_score[valid], q)
            keep = valid & (noise_score >= cut)
            share = hit[keep].sum() / hit.sum()
            rec = [(pair_score[bin_of == i] >= cut).mean() for i in range(5)]
            results[name, q] = (share, rec)
            print(f"{name:6} {1 - q:.2f}     {share:.3f}  " + " ".join(f"{v:.3f}" for v in rec), flush=True)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "wb") as f:
        pickle.dump(dict(results=results, hit=hit, gor=sg, kmer4=st), f)


if __name__ == "__main__":
    if sys.argv[1] == "value":
        value()
    else:
        main()
