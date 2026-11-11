"""Spaced seeds with a drift window, by secondary structure class.

A seed is two blocks: A of a residues at column c, B of b residues at column c + a + g. On the read, B
may start anywhere in i + a + g + d for |d| <= w (and never over A), so a short indel in the gap is
tolerated; nothing hashes an exact distance. Score: A's sum plus B's best sum over the window. A seed's
class is its middle column's. w = 0 with g = 0 is a contiguous (a + b)-mer, the control.

spaced.py noise PREP SS_OUT ORFS.faa N_READS WORKERS SHARE   hits per read, against a share of models
spaced.py pairs PREP SS_OUT ORFS.faa DOMTBL MAX_PER_BIN      best seed per class on hmmsearch's hits
spaced.py report SS_OUT                                      equal-hit tables, and one seed per class
spaced.py check

Recall counts a pair when A sits on a diagonal >= 22 bits, as ss_kmers.py does for contiguous k.
"""

import multiprocessing as mp
import pickle
import random
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

import gor
from ss_kmers import BAR, BIN_LOW, BIN_NAMES, diagonals, domtbl, fasta, load_prep

PATTERNS = [(a, g, b, w) for a in (2, 3) for b in (2, 3) for g in range(5) for w in range(3)]
CLASSES = {"pred": 3, "truth": 4}
BIN = 2
TBINS = 100
MIN_SCORE = 5.0


def deltas(a, g, w):
    return range(max(-w, -g), w + 1)


def name(p):
    a, g, b, w = p
    return f"{a}-{g}-{b}" + (f"±{w}" if w else "")


# ---- noise ----

SUB = TO_END = CLS = None


def init(prep, out, share):
    """A share of the models, whole, so every column has its model's neighbours."""
    global SUB, TO_END, CLS
    lo, offs = load_prep(prep)
    rng = np.random.default_rng(5)
    keep = np.flatnonzero(rng.random(len(offs) - 1) < share)
    idx = np.concatenate([np.arange(offs[m], offs[m + 1]) for m in keep])
    SUB = np.ascontiguousarray(np.asarray(lo)[idx].T)  # 21 x columns
    TO_END = np.concatenate([np.arange(offs[m + 1] - offs[m], 0, -1) for m in keep])
    z = np.load(Path(out) / "policies.npz")
    CLS = {p: z[p][idx].astype(np.int64) for p in CLASSES}
    return len(idx), len(lo)


def blocks(rows):
    """{k: L x n block sums} for k = 2, 3, 4 from L x n per-residue scores; -inf where a block runs off."""
    length, n = rows.shape
    out = {}
    for k in (2, 3, 4):
        s = np.full((length, n), -np.inf, dtype=np.float32)
        if length >= k:
            s[:length - k + 1, :n - k + 1] = sum(rows[j:length - k + 1 + j, j:n - k + 1 + j] for j in range(k))
        out[k] = s
    return out


def pattern_scores(blk, p):
    """L x n seed scores of pattern p, anchored at A's start."""
    a, g, b, w = p
    sa, sb = blk[a], blk[b]
    length, n = sa.shape
    s = a + g
    best = np.full((length, n), -np.inf, dtype=np.float32)
    for d in deltas(a, g, w):
        r = s + d
        if r < length and s < n:
            np.maximum(best[:length - r, :n - s], sb[r:, s:], out=best[:length - r, :n - s])
    return sa + best


def noise_orf(seq):
    x = gor.encode(seq).astype(np.int64)
    blk = blocks(SUB[x])
    hist = {p: {q: np.zeros((c, TBINS), dtype=np.int64) for q, c in CLASSES.items()} for p in PATTERNS}
    for p in PATTERNS:
        a, g, b, _ = p
        span = a + g + b
        v = pattern_scores(blk, p)
        i, c = np.nonzero(v >= MIN_SCORE)
        keep = TO_END[c] >= span
        i, c = i[keep], c[keep]
        tb = np.minimum((v[i, c] * BIN).astype(np.int64), TBINS - 1)
        mid = c + span // 2
        for q, nc in CLASSES.items():
            hist[p][q] += np.bincount(CLS[q][mid] * TBINS + tb, minlength=nc * TBINS).reshape(nc, TBINS)
    return hist, len(x)


def noise(prep, out, orfs, reads, workers, share):
    seqs = list(fasta(orfs, int(reads)).values())
    random.Random(1).shuffle(seqs)
    cols, total_cols = init(prep, out, float(share))
    print(f"{cols} of {total_cols} columns", flush=True)
    total = {p: {q: np.zeros((c, TBINS), dtype=np.int64) for q, c in CLASSES.items()} for p in PATTERNS}
    residues = 0
    with mp.get_context("fork").Pool(int(workers)) as pool:
        for n, (h, length) in enumerate(pool.imap_unordered(noise_orf, seqs, chunksize=2)):
            residues += length
            for p in PATTERNS:
                for q in CLASSES:
                    total[p][q] += h[p][q]
            if n % 50 == 0:
                print(n, len(seqs), flush=True)
    with open(Path(out) / "spaced_noise.pkl", "wb") as f:
        pickle.dump(dict(hist=total, reads=int(reads), scale=total_cols / cols, residues=residues,
                         columns={q: np.bincount(CLS[q], minlength=c) for q, c in CLASSES.items()}), f)
    print("noise done", residues)


# ---- pairs ----

def diag_blocks(a_mat):
    """{k: D x L block sums along each diagonal}; -inf past the read's end."""
    nd, length = a_mat.shape
    out = {}
    for k in (2, 3, 4):
        s = np.full((nd, length), -np.inf, dtype=np.float32)
        if length >= k:
            s[:, :length - k + 1] = sum(a_mat[:, j:length - k + 1 + j] for j in range(k))
        out[k] = s
    return out


def pair_scores(blk, p, good):
    """Seed scores anchored on the good diagonals: A at (d, i), B at (d - delta, i + a + g + delta)."""
    a, g, b, w = p
    sa, sb = blk[a], blk[b]
    nd, length = sa.shape
    best = np.full((len(good), length), -np.inf, dtype=np.float32)
    for d in deltas(a, g, w):
        r = a + g + d
        rows = good - d
        ok = (rows >= 0) & (rows < nd)
        if r < length and ok.any():
            best[ok, :length - r] = np.maximum(best[ok, :length - r], sb[rows[ok], r:])  # fancy index: no out=
    return sa[good] + best


def pairs(prep, out, orfs, table, per_bin):
    lo, offs = load_prep(prep)
    z = np.load(Path(out) / "policies.npz")
    cls = {q: z[q] for q in CLASSES}
    rows = [line.rstrip("\n").split("\t") for line in open(Path(prep) / "models.tsv")][1:]
    index = {r[0]: i for i, r in enumerate(rows)}
    ga = {r[0]: float(r[4]) for r in rows}
    seqs = fasta(orfs)
    best = {}
    for orf, acc, score, *_ in domtbl(table):
        if acc in index and orf in seqs:
            best[orf, acc] = max(score, best.get((orf, acc), -1e9))
    bins = defaultdict(list)
    for key, s in best.items():
        bins[next(i for i, low in enumerate(BIN_LOW) if s - ga[key[1]] >= low)].append(key)
    rng = random.Random(23)
    picks = [(key, b) for b in sorted(bins) for key in rng.sample(sorted(bins[b]), min(int(per_bin), len(bins[b])))]
    n = len(picks)
    tbl = {p: {q: np.full((n, c), -99.0, dtype=np.float32) for q, c in CLASSES.items()} for p in PATTERNS}
    bin_of = np.zeros(n, dtype=np.int8)
    for r, ((orf, acc), b) in enumerate(picks):
        bin_of[r] = b
        m = index[acc]
        start, end = offs[m], offs[m + 1]
        x = gor.encode(seqs[orf]).astype(np.int64)
        a_mat, diag, col = diagonals(np.asarray(lo[start:end]), x)
        good = np.flatnonzero(diag >= BAR)
        if not len(good):
            continue
        blk = diag_blocks(a_mat)
        for p in PATTERNS:
            v = pair_scores(blk, p, good)
            mid = col[good] + (p[0] + p[1] + p[2]) // 2
            ok = (v >= MIN_SCORE) & (mid >= 0) & (mid < end - start)
            if not ok.any():
                continue
            for q in CLASSES:
                np.maximum.at(tbl[p][q], (r, cls[q][start + mid[ok]].astype(np.int64)), v[ok])
        if r % 2000 == 0:
            print(r, n, flush=True)
    with open(Path(out) / "spaced_pairs.pkl", "wb") as f:
        pickle.dump(dict(tbl=tbl, bin=bin_of), f)
    print("pairs done", n)


# ---- report ----

def report(out):
    out = Path(out)
    no = pickle.load(open(out / "spaced_noise.pkl", "rb"))
    pa = pickle.load(open(out / "spaced_pairs.pkl", "rb"))
    reads, scale, tbl, b = no["reads"], no["scale"], pa["tbl"], pa["bin"]
    n = len(b)
    ts = np.arange(MIN_SCORE, 30.01, 0.5)

    def hits(p, q, c, t):
        return no["hist"][p][q][c, int(round(t * BIN)):].sum() * scale / reads

    def recall(ok, mask=None):
        mask = np.ones(n, dtype=bool) if mask is None else mask
        return np.array([ok[mask & (b == i)].mean() for i in range(5)])

    k4 = (2, 0, 2, 0)
    base_hits = sum(hits(k4, "pred", c, 7) for c in range(3))
    base_ok = tbl[k4]["pred"].max(1) >= 7
    print(f"k4 T7 (2-0-2): {base_hits:,.0f} hits per read; recall " + " ".join(f"{v:.3f}" for v in recall(base_ok)))
    print("pairs per bin", {BIN_NAMES[i]: int((b == i).sum()) for i in range(5)})

    print("\nOne pattern everywhere, at the T whose hits come nearest a share of k4 T7's. Recall by bin "
          + " ".join(BIN_NAMES) + ", and the mean of the three weak bins")
    for share in (1.0, 0.33, 0.1, 0.04):
        rows = []
        for p in PATTERNS:
            h = np.array([sum(hits(p, "pred", c, t) for c in range(3)) for t in ts])
            t = ts[np.argmin(np.abs(np.log(np.maximum(h, 1e-9) / (share * base_hits))))]
            rec = recall(tbl[p]["pred"].max(1) >= t)
            rows.append((rec[2:].mean(), p, t, h[ts == t][0] / base_hits, rec))
        rows.sort(reverse=True, key=lambda r: r[0])
        print(f"\nshare {share}: best 8 of {len(rows)}, then the contiguous controls")
        contiguous = [r for r in rows if r[1][1] == 0 and r[1][3] == 0]
        for weak, p, t, hs, rec in rows[:8] + contiguous:
            print(f"  {name(p):8} T{t:<5g} hits {hs:.3f}  " + " ".join(f"{v:.3f}" for v in rec) + f"  weak {weak:.3f}")

    print("\nBy class, at equal chance: the T where the class's hits per window x column match k4 T7's in"
          " that class, then the class's recall alone (pairs found by a seed in that class)")
    for q, names in (("pred", "HEC"), ("truth", "HEC?")):
        cols = no["columns"][q]
        for c, cn in enumerate(names):
            target = hits(k4, q, c, 7)
            rows = []
            for p in PATTERNS:
                h = np.array([hits(p, q, c, t) for t in ts])
                t = ts[np.argmin(np.abs(np.log(np.maximum(h, 1e-9) / target)))]
                rec = recall(tbl[p][q][:, c] >= t)
                rows.append((rec[2:].mean(), p, t, rec))
            rows.sort(reverse=True, key=lambda r: r[0])
            k4row = next(r for r in rows if r[1] == k4)
            print(f"{q:5} {cn}  columns {cols[c]}  k4 T7: " + " ".join(f"{v:.3f}" for v in k4row[3]))
            for weak, p, t, rec in rows[:5]:
                print(f"           {name(p):8} T{t:<5g} " + " ".join(f"{v:.3f}" for v in rec) + f"  weak {weak:.3f}")
            by_w = {w: max(r[0] for r in rows if r[1][3] == w) for w in range(3)}
            print("           best weak by window w: " + " ".join(f"w{w} {v:.3f}" for w, v in by_w.items()))

    # One (pattern, T) per class: least hits at recall within delta of k4 T7 in every bin, chosen on half
    # the pairs and reported on the other half
    rng = np.random.default_rng(77)
    half = rng.random(n) < 0.5
    options = [(p, t) for p in PATTERNS for t in ts]
    need_base = recall(base_ok, half)
    print("\nOne (pattern, T) per class; hits against k4 T7, held-out recall (baseline "
          + " ".join(f"{v:.3f}" for v in recall(base_ok, ~half)) + ")")
    for q, nc in CLASSES.items():
        cost = {(c, o): hits(o[0], q, c, o[1]) for c in range(nc) for o in options}
        okf = {}

        def ok_for(c, o):
            if (c, o) not in okf:
                okf[c, o] = tbl[o[0]][q][:, c] >= o[1]
            return okf[c, o]

        for delta in (0.0, 0.01, 0.02, 0.05):
            need = need_base - delta
            pol = [(k4, 7.0)] * nc
            for _ in range(6):
                changed = False
                for c in range(nc):
                    rest = np.zeros(n, dtype=bool)
                    for c2 in range(nc):
                        if c2 != c:
                            rest |= ok_for(c2, pol[c2])
                    cur = cost[c, pol[c]]
                    for o in options:
                        if cost[c, o] < cur:
                            ok = rest | ok_for(c, o)
                            if all(ok[half & (b == i)].mean() >= need[i] for i in range(5)):
                                pol[c], cur, changed = o, cost[c, o], True
                if not changed:
                    break
            ok = np.zeros(n, dtype=bool)
            for c, o in enumerate(pol):
                ok |= ok_for(c, o)
            total = sum(cost[c, o] for c, o in enumerate(pol)) / base_hits
            # Uniform control: the same search with one option shared by every class
            uni, ucost = None, np.inf
            for o in options:
                co = sum(cost[c, o] for c in range(nc))
                if co < ucost:
                    uok = np.zeros(n, dtype=bool)
                    for c in range(nc):
                        uok |= ok_for(c, o)
                    if all(uok[half & (b == i)].mean() >= need[i] for i in range(5)):
                        uni, ucost = o, co
            uok = np.zeros(n, dtype=bool)
            for c in range(nc):
                uok |= ok_for(c, uni)
            print(f"{q:5} delta {delta:.2f}: per class {total:.3f}  " + " ".join(f"{v:.3f}" for v in recall(ok, ~half))
                  + "  " + " ".join(f"{name(o[0])} T{o[1]:g}" for o in pol))
            print(f"{'':5}             uniform   {ucost / base_hits:.3f}  " + " ".join(f"{v:.3f}" for v in recall(uok, ~half))
                  + f"  {name(uni[0])} T{uni[1]:g}", flush=True)


def check():
    """The diagonal and full-matrix routes must give the same best seed of a pattern."""
    rng = np.random.default_rng(3)
    prof = rng.normal(-1, 2, (40, 21)).astype(np.float32)
    x = rng.integers(0, 20, 25)
    a_mat, _, col = diagonals(prof, x)
    good = np.arange(a_mat.shape[0])
    dblk = diag_blocks(a_mat)
    mblk = blocks(np.ascontiguousarray(prof.T)[x])  # read x column
    for p in PATTERNS:
        v = pair_scores(dblk, p, good)
        m = pattern_scores(mblk, p)
        span = p[0] + p[1] + p[2]
        cols = col[good]
        ok = (cols >= 0) & (cols + span <= len(prof)) & np.isfinite(v)
        ii = np.broadcast_to(np.arange(len(x)), cols.shape)
        assert np.allclose(v[ok], m[ii[ok], cols[ok]]), name(p)
    # A one-residue deletion in the read between blocks is found with w = 1, not w = 0
    prof = np.full((12, 21), -4.0, dtype=np.float32)
    seq = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    for c, r in enumerate(seq):
        prof[c, r] = 3.0
    read = np.array([2, 3, 5, 6])  # column 4 deleted, between the blocks of 2-1-2
    blk = diag_blocks(diagonals(prof, read)[0])
    good = np.arange(len(prof) + len(read) - 1)
    assert pair_scores(blk, (2, 1, 2, 1), good).max() == 12
    assert pair_scores(blk, (2, 1, 2, 0), good).max() < 12
    print("spaced check ok")


if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    {"noise": noise, "pairs": pairs, "report": report, "check": check}[cmd](*args)
