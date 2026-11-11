"""Per-profile T for k4: how skewed are hits across models, and what could a T per model save?

per_model.py noise PREP ORFS.faa N_READS WORKERS OUT    k4 hits per model by score, every column
per_model.py pairs PREP ORFS.faa DOMTBL OUT             best k4 seed on a diagonal >= 22, every pair
per_model.py report PREP OUT

Rules, each against one global T at the same hits:
    cap     T_m: the least T whose hits from model m stay under a budget (per model, or per column)
    oracle  T_m maximising recall - lambda x hits on the pairs themselves; models with no pairs keep T7.
            In-sample, so an upper bound; 'held out' fits on odd reads and scores on even ones.
Recall is the mean over bins, each bin weighted equally, as the picks in ss_kmers.py are.
"""

import multiprocessing as mp
import pickle
import sys
from pathlib import Path

import numpy as np

import gor
from ss_kmers import BAR, BIN_LOW, BIN_NAMES, diagonals, domtbl, fasta, load_prep

K, T0, NB = 4, 7.0, 30  # bins of 0.5 bit from 7
TS = T0 + np.arange(NB) / 2
OFF = NB  # option index meaning the model is never seeded

BY_RES = TO_END = MODEL = None
M = 0


def init(prep):
    global BY_RES, TO_END, MODEL, M
    lo, offs = load_prep(prep)
    BY_RES = np.array(lo).T.copy()
    M = len(offs) - 1
    TO_END = np.zeros(len(lo), dtype=np.int32)
    MODEL = np.zeros(len(lo), dtype=np.int64)
    for m, (a, b) in enumerate(zip(offs[:-1], offs[1:])):
        TO_END[a:b] = np.arange(b - a, 0, -1)
        MODEL[a:b] = m


def noise_orf(seq):
    x = gor.encode(seq).astype(np.int64)
    n = BY_RES.shape[1]
    keys = []
    for i in range(len(x) - K + 1):
        cur = BY_RES[x[i]].copy()
        for j in range(1, K):
            cur[:n - j] += BY_RES[x[i + j]][j:]
        idx = np.flatnonzero(cur >= T0)
        idx = idx[TO_END[idx] >= K]
        keys.append(MODEL[idx] * NB + np.minimum(((cur[idx] - T0) * 2).astype(np.int64), NB - 1))
    keys = np.concatenate(keys) if keys else np.zeros(0, dtype=np.int64)
    u, c = np.unique(keys, return_counts=True)
    return u, c


def noise(prep, orfs, reads, workers, out):
    seqs = list(fasta(orfs, int(reads)).values())
    init(prep)
    hist = np.zeros(M * NB, dtype=np.int64)
    with mp.get_context("fork").Pool(int(workers)) as pool:
        for n, (u, c) in enumerate(pool.imap_unordered(noise_orf, seqs, chunksize=4)):
            hist[u] += c
            if n % 100 == 0:
                print(n, len(seqs), flush=True)
    with open(Path(out) / "per_model_noise.pkl", "wb") as f:
        pickle.dump(dict(hist=hist.reshape(M, NB), reads=int(reads)), f)
    print("noise done", hist.sum() / int(reads))


def pairs(prep, orfs, table, out):
    lo, offs = load_prep(prep)
    rows = [line.rstrip("\n").split("\t") for line in open(Path(prep) / "models.tsv")][1:]
    index = {r[0]: i for i, r in enumerate(rows)}
    ga = {r[0]: float(r[4]) for r in rows}
    seqs = fasta(orfs)
    best = {}
    for orf, acc, score, *_ in domtbl(table):
        if acc in index and orf in seqs:
            best[orf, acc] = max(score, best.get((orf, acc), -1e9))
    keys = sorted(best)
    model = np.array([index[a] for _, a in keys])
    read = np.array([int(o[1:].split("_")[0]) for o, _ in keys])
    bin_of = np.array([next(i for i, low in enumerate(BIN_LOW) if best[k] - ga[k[1]] >= low) for k in keys])
    seed = np.full(len(keys), -np.inf, dtype=np.float32)
    for r, (orf, acc) in enumerate(keys):
        m = index[acc]
        x = gor.encode(seqs[orf]).astype(np.int64)
        if len(x) < K:
            continue
        a, diag, _ = diagonals(np.asarray(lo[offs[m]:offs[m + 1]]), x)
        good = a[diag >= BAR]
        if len(good):
            seed[r] = sum(good[:, j:len(x) - K + 1 + j] for j in range(K)).max()
        if r % 5000 == 0:
            print(r, len(keys), flush=True)
    with open(Path(out) / "per_model_pairs.pkl", "wb") as f:
        pickle.dump(dict(model=model, read=read, bin=bin_of, seed=seed), f)
    print("pairs done", len(keys))


def report(prep, out):
    out = Path(out)
    no = pickle.load(open(out / "per_model_noise.pkl", "rb"))
    pa = pickle.load(open(out / "per_model_pairs.pkl", "rb"))
    rows = [line.rstrip("\n").split("\t") for line in open(Path(prep) / "models.tsv")][1:]
    _, offs = load_prep(prep)
    length = np.diff(offs)
    # hits[m, t]: hits per read from model m at T = TS[t]; column OFF is zero
    hits = np.zeros((len(length), NB + 1))
    hits[:, :NB] = no["hist"][:, ::-1].cumsum(1)[:, ::-1] / no["reads"]
    total0 = hits[:, 0].sum()
    model, b, seed, read = pa["model"], pa["bin"], pa["seed"], pa["read"]
    weight = np.array([1 / (5 * (b == i).sum()) for i in range(5)])[b]
    # option index each pair needs at most: pair found at option t if seed >= TS[t]
    need = np.where(seed >= T0, np.minimum(((seed - T0) * 2).astype(np.int64), NB - 1), -1)

    def recall(t_of_model, mask=None):
        mask = np.ones(len(b), dtype=bool) if mask is None else mask
        t = t_of_model[model]
        ok = (t < OFF) & (need >= t)
        return np.array([ok[mask & (b == i)].mean() for i in range(5)])

    def cost(t_of_model):
        return hits[np.arange(len(hits)), t_of_model].sum() / total0

    def show(label, t_of_model, mask=None):
        r = recall(t_of_model, mask)
        print(f"{label:28} hits {cost(t_of_model):.3f}  " + " ".join(f"{v:.3f}" for v in r) + f"  mean {r.mean():.3f}")

    print(f"k4 T7: {total0:,.0f} hits per read from {no['reads']} reads; {len(b)} pairs; per bin "
          + " ".join(f"{BIN_NAMES[i]} {(b == i).sum()}" for i in range(5)))

    print("\nSkew of hits across models at T7")
    h = hits[:, 0]
    order = np.argsort(-h)
    cum = h[order].cumsum() / h.sum()
    for top in (0.001, 0.01, 0.05, 0.1, 0.25, 0.5):
        print(f"  top {top:6.1%} of models ({int(top * len(h)):5d}): {cum[int(top * len(h)) - 1]:.3f} of hits")
    per_col = h / length
    oc = np.argsort(-per_col)
    print(f"  hits per column: median {np.median(per_col):.3f}, p90 {np.quantile(per_col, .9):.3f}, p99 {np.quantile(per_col, .99):.3f}, "
          f"max {per_col.max():.2f}; top 1% of models by it hold {h[oc[:len(h) // 100]].sum() / h.sum():.3f} of hits")
    print(f"  correlation of a model's hits with its length: {np.corrcoef(h, length)[0, 1]:.2f}")
    has = np.bincount(model, minlength=len(h)) > 0
    has_ga = np.bincount(model[b == 0], minlength=len(h)) > 0
    print(f"  hits from models with no pair in this sample: {h[~has].sum() / h.sum():.3f} ({(~has).sum()} models); "
          f"with no GA pair: {h[~has_ga].sum() / h.sum():.3f}")
    npairs = np.bincount(model, minlength=len(h))
    print("  top 15 models by hits: acc name length hits/read hits/column pairs (>=GA)")
    for m in order[:15]:
        print(f"    {rows[m][0]} {rows[m][1]:16} {length[m]:5d} {h[m]:8.0f} {per_col[m]:6.2f} {npairs[m]:4d} ({(model[(b == 0)] == m).sum()})")

    print("\nGlobal T (recall by bin " + " ".join(BIN_NAMES) + ")")
    for t in range(0, 9):
        show(f"T{TS[t]:g}", np.full(len(h), t))

    def at_budget(budget):
        """Least T per model with hits under budget (an array, per model)."""
        ok = hits[:, :NB] <= budget[:, None]
        t = np.where(ok.any(1), ok.argmax(1), OFF)
        return t

    print("\nCap per model (hits per read), and per column (x model length), swept to a share of hits")
    for name, scale in (("cap per model", np.ones(len(h))), ("cap per column", length.astype(float))):
        for share in (0.8, 0.6, 0.4, 0.33, 0.2, 0.1):
            lo_, hi_ = 1e-6, h.max() / scale.min()
            for _ in range(60):
                mid = np.sqrt(lo_ * hi_)
                if cost(at_budget(mid * scale)) > share:
                    hi_ = mid
                else:
                    lo_ = mid
            show(f"{name} {share}", at_budget(lo_ * scale))

    def oracle(lam, fit):
        """Per model, the option maximising its pairs' weight found minus lam x hits; no pairs: T7."""
        t = np.zeros(len(h), dtype=np.int64)
        gain = np.zeros((len(h), NB + 1))
        sel = fit & (need >= 0)
        np.add.at(gain, (model[sel], need[sel]), weight[sel])
        gain = gain[:, ::-1].cumsum(1)[:, ::-1]  # found at option t if need >= t
        gain[:, OFF] = 0
        fitted = np.bincount(model[fit], minlength=len(h)) > 0
        obj = gain - lam * hits / total0
        t[fitted] = obj[fitted].argmax(1)
        return t

    print("\nOracle per-model T, swept by lambda; models without pairs keep T7")
    odd = read % 2 == 1
    for lam in (0.003, 0.01, 0.03, 0.1, 0.3, 1, 3):
        show(f"in-sample lambda {lam}", oracle(lam, np.ones(len(b), dtype=bool)))
    print("held out: fit on odd reads' pairs, scored on even reads' (global T on even reads for comparison)")
    for t in (0, 1, 2, 3):
        show(f"  global T{TS[t]:g}", np.full(len(h), t), ~odd)
    for lam in (0.003, 0.01, 0.03, 0.1, 0.3, 1):
        show(f"  held-out lambda {lam}", oracle(lam, odd), ~odd)


if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    {"noise": noise, "pairs": pairs, "report": report}[cmd](*args)
