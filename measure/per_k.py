"""Per-profile seed shape and T: each model takes the cheapest (shape, T) its own profile says still finds
its homologs, so a model where k7 works never pays for k4.

per_k.py noise PREP ORFS.faa N_READS WORKERS OUT   hits per model by shape and score, every column
per_k.py pairs PREP ORFS.faa DOMTBL OUT            best seed per shape on a diagonal >= 22, every pair
per_k.py profile PREP OUT WINDOWS                  per model: predicted detection, index entries
per_k.py report PREP OUT

Detection is predicted from the profile alone: a window's score under emissions mixed with background
(alpha emissions, 1 - alpha background, for homologs more distant than the model's own), as a normal
approximation, then P(a seed in a 30-residue fragment) = 1 - exp(-27 x mean window rate). Index entries
are exact counts of words >= T, by convolution on sampled windows, scaled to the model's windows.
"""

import multiprocessing as mp
import pickle
import sys
from pathlib import Path

import numpy as np
from scipy.special import ndtr

import gor
from spaced import blocks, diag_blocks, pair_scores, pattern_scores
from ss_kmers import BAR, BIN_LOW, BIN_NAMES, Q, diagonals, domtbl, fasta, load_prep

SHAPES = {"k4": (2, 0, 2), "k5": (2, 0, 3), "3-1-3": (3, 1, 3), "3-2-3": (3, 2, 3), "k6": (3, 0, 3), "k7": (3, 0, 4)}
T0 = {"k4": 7.0, "k5": 7.5, "3-1-3": 8.0, "3-2-3": 8.0, "k6": 8.0, "k7": 9.0}
NB = 32  # bins of 0.5 bit from T0
CHUNK = 400_000
FRAG = 30
ALPHAS = (1.0, 0.5)


def span(s):
    a, g, b = SHAPES[s]
    return a + g + b


def offsets(s):
    a, g, b = SHAPES[s]
    return list(range(a)) + list(range(a + g, a + g + b))


def ts(s):
    return T0[s] + np.arange(NB) / 2


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
    keys = {s: [] for s in SHAPES}
    pad = max(span(s) for s in SHAPES)
    for c0 in range(0, n, CHUNK):
        c1 = min(c0 + CHUNK, n)
        blk = blocks(BY_RES[:, c0:min(c1 + pad, n)][x])
        for s, (a, g, b) in SHAPES.items():
            v = pattern_scores(blk, (a, g, b, 0))[:, :c1 - c0]
            i, c = np.nonzero(v >= T0[s])
            c = c + c0
            keep = TO_END[c] >= span(s)
            i, c = i[keep], c[keep]
            keys[s].append(MODEL[c] * NB + np.minimum(((v[i, c - c0] - T0[s]) * 2).astype(np.int64), NB - 1))
    out = {}
    for s in SHAPES:
        k = np.concatenate(keys[s])
        out[s] = np.unique(k, return_counts=True)
    return out


def noise(prep, orfs, reads, workers, out):
    named = fasta(orfs, int(reads))
    seqs = list(named.values())
    half = [int(n[1:].split("_")[0]) % 2 for n in named]  # odd and even reads, to choose on one and score on the other
    init(prep)
    hist = {s: np.zeros((2, M * NB), dtype=np.int64) for s in SHAPES}
    with mp.get_context("fork").Pool(int(workers)) as pool:
        for n, h in enumerate(pool.imap(noise_orf, seqs, chunksize=2)):
            for s, (u, c) in h.items():
                hist[s][half[n], u] += c
            if n % 50 == 0:
                print(n, len(seqs), flush=True)
    with open(Path(out) / "per_k_noise.pkl", "wb") as f:
        pickle.dump(dict(hist={s: h.sum(0).reshape(M, NB) for s, h in hist.items()},
                         halves={s: h.reshape(2, M, NB) for s, h in hist.items()},
                         half_reads=[len({n.split("_")[0] for n, q in zip(named, half) if q == i}) for i in (0, 1)],
                         reads=int(reads)), f)
    print("noise done", {s: hist[s].sum() / int(reads) for s in SHAPES})


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
    bin_of = np.array([next(i for i, low in enumerate(BIN_LOW) if best[k] - ga[k[1]] >= low) for k in keys])
    seed = {s: np.full(len(keys), -np.inf, dtype=np.float32) for s in SHAPES}
    for r, (orf, acc) in enumerate(keys):
        m = index[acc]
        x = gor.encode(seqs[orf]).astype(np.int64)
        a_mat, diag, _ = diagonals(np.asarray(lo[offs[m]:offs[m + 1]]), x)
        good = np.flatnonzero(diag >= BAR)
        if not len(good):
            continue
        blk = diag_blocks(a_mat)
        for s, (a, g, b) in SHAPES.items():
            v = pair_scores(blk, (a, g, b, 0), good)
            if v.size:
                seed[s][r] = v.max()
        if r % 5000 == 0:
            print(r, len(keys), flush=True)
    with open(Path(out) / "per_k_pairs.pkl", "wb") as f:
        pickle.dump(dict(model=model, bin=bin_of, seed=seed), f)
    print("pairs done", len(keys))


def count_at_least(cols, weights, thresholds, res=4):
    """Weighted count of words with sum >= each threshold; cols: n x 20 scores, weights: n x 20."""
    ints = np.floor(cols * res).astype(np.int64)
    low = ints.min(1)
    dist = np.ones(1)
    for c in range(len(cols)):
        h = np.bincount(ints[c] - low[c], weights=weights[c])
        dist = np.convolve(dist, h)
    base = low.sum()
    tail = dist[::-1].cumsum()[::-1]  # tail[i]: words with sum >= (base + i) / res
    i = np.ceil(np.asarray(thresholds) * res).astype(np.int64) - base
    return np.where(i < 0, tail[0], np.where(i >= len(tail), 0, tail[np.clip(i, 0, len(tail) - 1)]))


def profile(prep, out, windows):
    lo, offs = load_prep(prep)
    lo = np.asarray(lo)[:, :20].astype(np.float64)
    n, M_ = len(lo), len(offs) - 1
    to_end = np.zeros(n, dtype=np.int64)
    model = np.zeros(n, dtype=np.int64)
    for m, (a, b) in enumerate(zip(offs[:-1], offs[1:])):
        to_end[a:b] = np.arange(b - a, 0, -1)
        model[a:b] = m
    emit = Q * 2.0 ** lo
    emit /= emit.sum(1, keepdims=True)
    det, entries, nwin = {}, {}, {}
    for s in SHAPES:
        sp, off = span(s), offsets(s)
        ok = np.flatnonzero(to_end >= sp)
        nwin[s] = np.bincount(model[ok], minlength=M_)
        for alpha in ALPHAS:
            e = alpha * emit + (1 - alpha) * Q
            mu = (e * lo).sum(1)
            var = (e * lo ** 2).sum(1) - mu ** 2
            wm = sum(mu[ok + o] for o in off)
            ws = np.sqrt(sum(var[ok + o] for o in off))
            rate = np.zeros((M_, NB))
            for t_i, t in enumerate(ts(s)):
                rate[:, t_i] = np.bincount(model[ok], weights=ndtr((wm - t) / ws), minlength=M_)
            rate /= np.maximum(nwin[s], 1)[:, None]
            det[s, alpha] = 1 - np.exp(-(FRAG - sp + 1) * rate)
        # Entries: exact counts on sampled windows
        rng = np.random.default_rng(11)
        ent = np.zeros((M_, NB))
        seen = np.zeros(M_)
        starts = ok[np.argsort(rng.random(len(ok)))]
        taken = np.zeros(M_, dtype=np.int64)
        for w in starts:
            m = model[w]
            if taken[m] >= int(windows):
                continue
            taken[m] += 1
            ent[m] += count_at_least(lo[[w + o for o in off]], np.ones((len(off), 20)), ts(s))
            seen[m] += 1
        entries[s] = ent / np.maximum(seen, 1)[:, None] * nwin[s][:, None]
        print(s, "entries at T0", f"{entries[s][:, 0].sum():.3g}", flush=True)
    with open(Path(out) / "per_k_profile.pkl", "wb") as f:
        pickle.dump(dict(det=det, entries=entries), f)


def report(prep, out):
    out = Path(out)
    no = pickle.load(open(out / "per_k_noise.pkl", "rb"))
    pa = pickle.load(open(out / "per_k_pairs.pkl", "rb"))
    pr = pickle.load(open(out / "per_k_profile.pkl", "rb"))
    model, b = pa["model"], pa["bin"]
    M_ = no["hist"]["k4"].shape[0]
    opts = [(s, t) for s in SHAPES for t in range(NB)]
    # hits[m, o], entries[m, o], found[pair, o]
    hits = np.stack([no["hist"][s][:, ::-1].cumsum(1)[:, ::-1][:, t] / no["reads"] for s, t in opts], 1)
    ent = np.stack([pr["entries"][s][:, t] for s, t in opts], 1)
    found = np.stack([pa["seed"][s] >= ts(s)[t] for s, t in opts], 1)
    base = opts.index(("k4", 0))
    h0, e0 = hits[:, base].sum(), ent[:, base].sum()

    def show(label, choice):
        """choice: option index per model, -1 for off."""
        mi = np.arange(M_)
        on = choice >= 0
        h = hits[mi[on], choice[on]].sum() / h0
        e = ent[mi[on], choice[on]].sum() / e0
        c = choice[model]
        ok = (c >= 0) & found[np.arange(len(model)), np.maximum(c, 0)]
        r = np.array([ok[b == i].mean() for i in range(5)])
        print(f"{label:34} hits {h:.3f} entries {e:6.2f}  " + " ".join(f"{v:.3f}" for v in r) + f"  mean {r.mean():.3f}", flush=True)
        return h, r

    print(f"k4 T7: {h0:,.0f} hits per read, {e0:.3g} index entries; {len(b)} pairs")
    print("\nOne shape everywhere (recall by bin " + " ".join(BIN_NAMES) + ")")
    for s in SHAPES:
        for t in range(0, NB, 2):
            h = hits[:, opts.index((s, t))].sum() / h0
            if 0.01 < h < 1.3:
                show(f"{s} T{ts(s)[t]:g}", np.full(M_, opts.index((s, t))))

    print("\nPer model, under an index budget: the option with least hits + mu x entries among those whose"
          " predicted detection is >= rho x k4 T7's. Hits to choose by: odd reads, scored on even reads;"
          " the same reads (optimistic); or entries x the option's global hits per entry (profile only)")
    mi = np.arange(M_)
    h_half = [np.stack([no["halves"][s][i][:, ::-1].cumsum(1)[:, ::-1][:, t] / no["half_reads"][i] for s, t in opts], 1)
              for i in (0, 1)]
    per_entry = hits.sum(0) / np.maximum(ent.sum(0), 1e-9)
    sources = {"held out": (h_half[1], h_half[0]), "same reads": (hits, hits), "profile only": (ent * per_entry, hits)}
    for alpha in ALPHAS:
        det = np.stack([pr["det"][s, alpha][:, t] for s, t in opts], 1)
        for rho in (1.0, 0.97):
            ok = det >= rho * det[:, [base]] - 1e-12
            ok[:, base] = True
            for src, (choose_by, score_by) in sources.items():
                e_base = score_by[:, base].sum()
                for mu in (1e-2, 1e-3, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 0):
                    choice = np.where(ok, choose_by / choose_by[:, base].sum() + mu * ent / e0 * M_ / 100, np.inf).argmin(1)
                    c = choice[model]
                    found_ = found[np.arange(len(model)), c]
                    r = np.array([found_[b == i].mean() for i in range(5)])
                    print(f"alpha {alpha} rho {rho} {src:12} mu {mu:<6g} hits {score_by[mi, choice].sum() / e_base:.3f} "
                          f"entries {ent[mi, choice].sum() / e0:6.1f}  " + " ".join(f"{v:.3f}" for v in r), flush=True)

if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    {"noise": noise, "pairs": pairs, "profile": profile, "report": report}[cmd](*args)
