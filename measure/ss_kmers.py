"""Do k-mers from different secondary structure differ in value, and can k follow the structure?

A column's class is known when the index is built, so choosing k or T per class costs nothing per
read. Measured on real reads against hmmsearch, by brute force:

ss_kmers.py classes PREP SS_OUT Pfam-A.hmm     column classes of each policy, and indel rates
ss_kmers.py noise   PREP SS_OUT ORFS.faa N WORKERS   seed hits per read by k, score and class
ss_kmers.py pairs   PREP SS_OUT ORFS.faa DOMTBL MAX_PER_BIN   the same on hmmsearch's hits
ss_kmers.py entries PREP SS_OUT N              index entries by k, T and class, from N sampled windows
ss_kmers.py report  SS_OUT                     tables, and the best policies at equal recall
ss_kmers.py web     PREP SS_OUT WEB_DATA       per-column structure for the explorer's family view
ss_kmers.py table   SS_OUT                     k against class at equal chance rates
ss_kmers.py check                              self-check of the diagonal and seed scoring

Policies are ways of classing columns: pred (GOR on emissions), truth (DSSP where a chain aligns,
else class 3), shuffle (pred rotated within each model: same segments, wrong place), info and indel
(terciles of column information and of the indel probability) and predinfo (pred x info).
"""

import gzip
import multiprocessing as mp
import pickle
import random
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

import gor

KS = (4, 5, 6, 7, 8)
MIN_T = {4: 6.0, 5: 7.0, 6: 9.0, 7: 10.0, 8: 12.0}
BIN = 2  # bins per bit
TBINS = 80
POLICIES = {"pred": 3, "truth": 4, "shuffle": 3, "info": 3, "indel": 3, "predinfo": 9}
BIN_LOW = (0, -5, -10, -20, -np.inf)
BIN_NAMES = (">=GA", "-5-GA", "-10--5", "-20--10", "<-20")
BAR = 22.0
Q = np.array([.0787945, .0151600, .0535222, .0668298, .0397062, .0695071, .0229198, .0590092,
              .0594422, .0963728, .0237718, .0414386, .0482904, .0395639, .0540978, .0683364,
              .0540687, .0673417, .0114135, .0304133])
Q /= Q.sum()


def load_prep(prep):
    prep = Path(prep)
    offs = np.fromfile(prep / "offs.u32", dtype=np.uint32).astype(np.int64)
    lo = np.memmap(prep / "profile.f32", dtype=np.float32, mode="r").reshape(-1, 21)
    return lo, offs


def fasta(path, limit=None):
    """{orf name: sequence}; limit keeps reads below that number."""
    out, name = {}, None
    with open(path) as f:
        for line in f:
            if line[0] == ">":
                name = line[1:].split()[0]
                if limit is not None and int(name[1:].split("_")[0]) >= limit:
                    name = None
            elif name:
                out[name] = line.strip()
    return out


def terciles(x):
    a, b = np.quantile(x, [1 / 3, 2 / 3])
    return ((x > a).astype(np.int8) + (x > b).astype(np.int8)).astype(np.int8)


def classes(prep, out, hmm):
    lo, offs = load_prep(prep)
    out = Path(out)
    n = len(lo)
    z = np.load(out / "classes.npz")
    pred, truth = z["pred"], z["truth"]
    tmm = np.zeros(n, dtype=np.float32)
    i = 0
    with open(hmm) as f:
        inside = False
        for line in f:
            if line.startswith("HMM "):
                inside = True
                next(f), next(f), next(f), next(f)  # transition header, COMPO, its two lines
            elif line.startswith("//"):
                inside = False
            elif inside:
                next(f)
                tmm[i] = float(next(f).split()[0])
                i += 1
    assert i == n, (i, n)
    pgap = 1 - np.exp(-tmm)
    lo20 = np.asarray(lo[:, :20])
    p = Q * np.exp2(lo20)
    info = (p * lo20).sum(1).astype(np.float32)

    def smooth(x, w=5):
        s = np.zeros(n, dtype=np.float32)
        for a, b in zip(offs[:-1], offs[1:]):
            s[a:b] = np.convolve(x[a:b], np.ones(w) / w, "same")
        return s

    rng = np.random.default_rng(9)
    shuffle = pred.copy()
    for a, b in zip(offs[:-1], offs[1:]):
        shuffle[a:b] = np.roll(pred[a:b], rng.integers(0, b - a))
    info_class = terciles(smooth(info))
    np.savez_compressed(out / "policies.npz", pred=pred, truth=truth, shuffle=shuffle, info=info_class,
                        indel=terciles(smooth(pgap)), predinfo=pred * 3 + info_class, pgap=pgap, rel=info)
    print("states", n)


class Scorer:
    def __init__(self, prep, out):
        self.lo, self.offs = load_prep(prep)
        z = np.load(Path(out) / "policies.npz")
        self.cls = {p: z[p] for p in POLICIES}
        n = len(self.lo)
        self.to_end = np.zeros(n, dtype=np.int32)
        for a, b in zip(self.offs[:-1], self.offs[1:]):
            self.to_end[a:b] = np.arange(b - a, 0, -1)


S = None
TAB = None


def init_noise(prep, out):
    global S, TAB
    S = Scorer(prep, out)
    S.by_res = np.array(S.lo).T.copy()  # 21 x states
    with open(Path(out) / "gor_tables.pkl", "rb") as f:
        TAB = pickle.load(f)[0]


def noise_orf(seq):
    """{k: {policy: classes x 3 x bins}} of seeds in one ORF, against every column."""
    x = gor.encode(seq).astype(np.int64)
    rc = TAB.predict_seq(x.astype(np.int8))
    hist = {k: {p: np.zeros((c, 3, TBINS), dtype=np.int64) for p, c in POLICIES.items()} for k in KS}
    n = S.by_res.shape[1]
    for i in range(len(x) - KS[0] + 1):
        cur = S.by_res[x[i]].copy()
        for j in range(1, KS[-1]):
            if i + j >= len(x):
                break
            cur[:n - j] += S.by_res[x[i + j]][j:]
            k = j + 1
            if k not in hist:
                continue
            idx = np.flatnonzero(cur >= MIN_T[k])
            idx = idx[S.to_end[idx] >= k]
            if not len(idx):
                continue
            tb = np.minimum((cur[idx] * BIN).astype(np.int64), TBINS - 1)
            read = rc[i + (k - 1) // 2]
            mid = idx + (k - 1) // 2
            for p, c in POLICIES.items():
                hist[k][p][:, read] += np.bincount(S.cls[p][mid].astype(np.int64) * TBINS + tb, minlength=c * TBINS).reshape(c, TBINS)
    return hist, len(x)


def noise(prep, out, orfs, reads, workers):
    seqs = list(fasta(orfs, int(reads)).values())
    random.Random(1).shuffle(seqs)
    total = {k: {p: np.zeros((c, 3, TBINS), dtype=np.int64) for p, c in POLICIES.items()} for k in KS}
    residues = 0
    ctx = mp.get_context("fork")
    init_noise(prep, out)
    with ctx.Pool(int(workers)) as pool:
        for n, (h, length) in enumerate(pool.imap_unordered(noise_orf, seqs, chunksize=4)):
            residues += length
            for k in KS:
                for p in POLICIES:
                    total[k][p] += h[k][p]
            if n % 100 == 0:
                print(n, len(seqs), residues, flush=True)
    with open(Path(out) / "noise.pkl", "wb") as f:
        pickle.dump({"hist": total, "reads": int(reads), "orfs": len(seqs), "residues": residues}, f)
    print("noise done", residues)


def diagonals(prof, x):
    """Matrix (diagonals x residues) of per-residue scores, column = d - (L - 1) + i, and each diagonal's
    best ungapped local score. Off-model cells are -100."""
    m, length = len(prof), len(x)
    d = np.arange(m + length - 1)[:, None]
    col = d - (length - 1) + np.arange(length)[None, :]
    ok = (col >= 0) & (col < m)
    a = np.full(col.shape, -100.0, dtype=np.float32)
    a[ok] = prof[col[ok], np.broadcast_to(x, col.shape)[ok]]
    h = np.zeros(m + length - 1, dtype=np.float32)
    best = h.copy()
    for i in range(length):
        h = np.maximum(h, 0) + a[:, i]
        best = np.maximum(best, h)
    return a, best, col


def domtbl(path):
    with open(path) as f:
        for line in f:
            if line[0] != "#":
                c = line.split()
                yield c[0], c[4].split(".")[0], float(c[13]), int(c[15]), int(c[16]), int(c[17]), int(c[18])


def pairs(prep, out, orfs, table, per_bin):
    s = Scorer(prep, out)
    out = Path(out)
    rows = [line.rstrip("\n").split("\t") for line in open(Path(prep) / "models.tsv")][1:]
    index = {r[0]: i for i, r in enumerate(rows)}
    ga = {r[0]: float(r[4]) for r in rows}
    with open(out / "gor_tables.pkl", "rb") as f:
        tables = pickle.load(f)
    sys.path.insert(0, str(Path(__file__).parent))
    from ss_accuracy import model_folds
    _, fold = model_folds(Path(prep))
    seqs = fasta(orfs)
    best = {}
    for orf, acc, score, hf, ht, af, at in domtbl(table):
        if acc in index and orf in seqs and ((orf, acc) not in best or score > best[orf, acc][0]):
            best[orf, acc] = (score, hf, ht, af, at)
    bins = defaultdict(list)
    for key, v in best.items():
        margin = v[0] - ga[key[1]]
        bins[next(i for i, low in enumerate(BIN_LOW) if margin >= low)].append(key)
    rng = random.Random(23)
    picks = []
    print({BIN_NAMES[b]: len(v) for b, v in sorted(bins.items())}, flush=True)
    for b in sorted(bins):
        v = sorted(bins[b])
        picks += [(key, b) for key in rng.sample(v, min(int(per_bin), len(v)))]
    n = len(picks)
    tbl = {k: {p: np.full((n, c, 3), -99.0, dtype=np.float32) for p, c in POLICIES.items()} for k in KS}
    best_diag = np.zeros(n, dtype=np.float32)
    bin_of = np.zeros(n, dtype=np.int8)
    # Windows on the alignment's own diagonal, gapless hits only: true-positive rate by score
    tp = {k: {p: np.zeros((5, c, TBINS), dtype=np.int64) for p, c in POLICIES.items()} for k in KS}
    tp_n = {k: {p: np.zeros((5, c), dtype=np.int64) for p, c in POLICIES.items()} for k in KS}
    conf = np.zeros((5, 3, 3), dtype=np.int64)  # read predicted x column predicted
    conf_t = np.zeros((5, 3, 3), dtype=np.int64)  # read predicted x column truth
    conf_p = np.zeros((5, 3, 3), dtype=np.int64)  # column predicted x column truth
    per_col = defaultdict(lambda: np.zeros(3, dtype=np.int64))
    for r, (key, b) in enumerate(picks):
        orf, acc = key
        _, hf, ht, af, at = best[key]
        m = index[acc]
        start, end = s.offs[m], s.offs[m + 1]
        prof = np.asarray(s.lo[start:end])
        x = gor.encode(seqs[orf]).astype(np.int64)
        rc = tables[fold[acc]].predict_seq(x.astype(np.int8))
        a, diag, col = diagonals(prof, x)
        best_diag[r] = diag.max()
        bin_of[r] = b
        good = np.flatnonzero(diag >= BAR)
        cls_col = {p: s.cls[p][start:end] for p in POLICIES}
        gapless = ht - hf == at - af
        if gapless:
            d_al = (hf - af) + len(x) - 1  # 1-based both, so the offset carries over
            cols = np.arange(af - 1, at)
            if col[d_al, cols[0]] != hf - 1:
                raise ValueError("diagonal bookkeeping")
            read_c, c_col = rc[cols], cols + hf - af
            pr, tr = cls_col["pred"][c_col], cls_col["truth"][c_col]
            conf[b] += np.bincount(read_c * 3 + pr, minlength=9).reshape(3, 3)
            w = tr < 3
            conf_t[b] += np.bincount(read_c[w] * 3 + tr[w], minlength=9).reshape(3, 3)
            conf_p[b] += np.bincount(pr[w].astype(np.int64) * 3 + tr[w], minlength=9).reshape(3, 3)
            for c0, rcl in zip(c_col, read_c):
                per_col[start + c0][rcl] += 1
        for k in KS:
            cum = a[:, :len(x) - k + 1].copy() if len(x) >= k else None
            if cum is None:
                continue
            for j in range(1, k):
                cum += a[:, j:len(x) - k + 1 + j]
            if gapless and at - af + 1 >= k:
                w = cum[d_al, af - 1:at - k + 1]
                cc = (col[d_al, af - 1:at - k + 1] + (k - 1) // 2)
                tb = np.clip((w * BIN).astype(np.int64), 0, TBINS - 1)  # negatives land in bin 0
                for p, nc in POLICIES.items():
                    c = cls_col[p][cc]
                    np.add.at(tp[k][p][b], (c, tb), 1)
                    tp_n[k][p][b] += np.bincount(c, minlength=nc)
            if not len(good):
                continue
            gd, gi = np.nonzero(cum[good] >= MIN_T[k])
            if not len(gd):
                continue
            v = cum[good][gd, gi]
            mid = col[good[gd], gi] + (k - 1) // 2
            rcmid = rc[gi + (k - 1) // 2]
            for p in POLICIES:
                c = cls_col[p][mid]
                t = tbl[k][p]
                np.maximum.at(t, (r, c, 0), v)
                if p == "pred":
                    agree = c == rcmid
                    np.maximum.at(t, (r, c[agree], 1), v[agree])
                    no_conflict = ~(((c == 0) & (rcmid == 1)) | ((c == 1) & (rcmid == 0)))
                    np.maximum.at(t, (r, c[no_conflict], 2), v[no_conflict])
        if r % 2000 == 0:
            print(r, n, flush=True)
    agree = other = 0
    obs = np.zeros(3)
    for v in per_col.values():
        t = v.sum()
        if t > 1:
            agree += (v * (v - 1)).sum()
            other += t * (t - 1)
            obs += v
    with open(out / "pairs.pkl", "wb") as f:
        pickle.dump(dict(tbl=tbl, best=best_diag, bin=bin_of, tp=tp, tp_n=tp_n, conf=conf, conf_t=conf_t, conf_p=conf_p,
                         same=agree / max(other, 1), chance=((obs / obs.sum()) ** 2).sum(),
                         columns=sum(v.sum() > 1 for v in per_col.values())), f)
    print("pairs done", n)


def entries(prep, out, n):
    lo, offs = load_prep(prep)
    npz = np.load(Path(out) / "policies.npz")
    z = {p: npz[p] for p in POLICIES}  # npz members decompress on every access
    rng = np.random.default_rng(31)
    states = len(lo)
    pick = np.sort(rng.choice(states, int(n), replace=False))
    r = 4  # bins per bit in the convolution
    offset = 400  # bins of slack for negative sums
    size = 1200
    h = {k: {p: np.zeros((c, TBINS), dtype=np.float64) for p, c in POLICIES.items()} for k in KS}
    seen = {k: {p: np.zeros(c, dtype=np.int64) for p, c in POLICIES.items()} for k in KS}
    ends = np.searchsorted(offs, pick, side="right")
    for s, e in zip(pick, offs[ends]):
        rows = np.asarray(lo[s:min(s + KS[-1], e), :20])
        cur = np.zeros(size)
        cur[offset] = 1
        for j, row in enumerate(rows):
            b = np.rint(row * r).astype(np.int64)
            nxt = np.zeros(size)
            nz = np.flatnonzero(cur)
            for v in b:
                t = nz + v
                t = t[(t >= 0) & (t < size)]
                nxt[t] += cur[t - v]
            cur = nxt
            k = j + 1
            if k in h:
                bits = (np.arange(size) - offset) / r
                bins = np.clip((bits * BIN).astype(np.int64), 0, TBINS - 1)
                pos = cur * (bits >= 0)
                hist = np.bincount(bins, weights=pos, minlength=TBINS)
                mid = s + (k - 1) // 2
                for p in POLICIES:
                    c = z[p][mid]
                    h[k][p][c] += hist
                    seen[k][p][c] += 1
    with open(Path(out) / "entries.pkl", "wb") as f:
        pickle.dump(dict(h=h, seen=seen, n=int(n), states=states), f)
    print("entries done")


def report(out):
    out = Path(out)
    with open(out / "noise.pkl", "rb") as f:
        no = pickle.load(f)
    with open(out / "pairs.pkl", "rb") as f:
        pa = pickle.load(f)
    with open(out / "entries.pkl", "rb") as f:
        en = pickle.load(f)
    npz = np.load(out / "policies.npz")
    z = {f: npz[f] for f in npz.files}
    reads = no["reads"]
    print(f"noise from {reads} reads, {no['orfs']} ORFs, {no['residues']} residues")

    columns = {p: np.bincount(z[p], minlength=c) for p, c in POLICIES.items()}

    def hits(k, p, c, t, gate=0):
        """Seeds per read from class c at k, T >= t."""
        hb = no["hist"][k][p][c]
        lo = int(round(t * BIN))
        if gate == 0:
            return hb[:, lo:].sum() / reads
        if gate == 1:
            return hb[c, lo:].sum() / reads
        keep = [r for r in range(3) if not ((c == 0 and r == 1) or (c == 1 and r == 0))]
        return hb[keep, lo:].sum() / reads

    def entry(k, p, c, t):
        lo = int(round(t * BIN))
        seen = en["seen"][k][p][c]
        share = columns[p][c] / seen if seen else 0
        return en["h"][k][p][c, lo:].sum() * share

    b = pa["bin"]
    tbl = pa["tbl"]
    n = len(b)
    print("pairs per bin", {BIN_NAMES[i]: int((b == i).sum()) for i in range(5)})
    base = tbl[4]["pred"][:, :, 0].max(1) >= 7
    print("baseline k4 T7, diagonal >= 22, recall by bin (plan: .994 .94 .91 .87 .83):",
          " ".join(f"{base[b == i].mean():.3f}" for i in range(5)))
    print(f"hits per read k4 T7 {sum(hits(4, 'pred', c, 7) for c in range(3)):,.0f}; "
          f"entries {sum(entry(4, 'pred', c, 7) for c in range(3)):.3g}; plan 221,000 and 2.8e8")
    print(f"entries k5 T8 {sum(entry(5, 'pred', c, 8) for c in range(3)):.3g} (plan 3.2e9)")

    print("\nRead's GOR class against the column, gapless hits (rows read H E C):")
    for i in range(5):
        c = pa["conf"][i]
        ct = pa["conf_t"][i]
        cp = pa["conf_p"][i]
        print(f"{BIN_NAMES[i]:>8}: read=emission {c.trace() / max(c.sum(), 1):.3f} of {c.sum()};"
              f" read=DSSP {ct.trace() / max(ct.sum(), 1):.3f}, emission=DSSP {cp.trace() / max(cp.sum(), 1):.3f} of {ct.sum()}")
    print(f"Two reads at one column get the same class {pa['same']:.3f} (chance {pa['chance']:.3f}), {pa['columns']} columns")

    print("\nIndel probability by column class (HMM m->m; chain gap events per aligned residue):")
    pg = z["pgap"]
    cz = np.load(out / "classes.npz")
    for p, names in (("pred", "HEC"), ("truth", "HEC")):
        for c, name in enumerate(names):
            w = z[p] == c
            print(f"  {p:6} {name}: columns {w.sum():8d}  1-tmm {pg[w].mean():.4f}  chain gaps/residue {cz['gaps'][w].sum() / max(cz['aligned'][w].sum(), 1):.4f}  info {z['rel'][w].mean():.3f}")

    print("\nPer class, seed value on real hits: rate on the true diagonal against chance per window x column")
    print("policy class k T  hits/read  entries  tp >=GA  tp -10--5  lr >=GA  lr -10--5  pairs only here (>=GA, -10--5)")
    for p, names in (("pred", "HEC"), ("truth", "HEC")):
        for k, ts in ((4, (7, 8, 9)), (5, (8, 9, 10)), (6, (10, 12)), (7, (12, 14)), (8, (14, 16))):
            for t in ts:
                for c, name in enumerate(names):
                    lo = int(round(t * BIN))
                    h_ = hits(k, p, c, t)
                    chance = h_ * reads / max(no["residues"] * columns[p][c], 1)
                    tpg = pa["tp"][k][p][0, c, lo:].sum() / max(pa["tp_n"][k][p][0, c], 1)
                    tpw = pa["tp"][k][p][2, c, lo:].sum() / max(pa["tp_n"][k][p][2, c], 1)
                    print(f"{p:6} {name} k{k} T{t:<3} {h_:10,.0f} {entry(k, p, c, t):9.3g}  {tpg:.3f}  {tpw:.3f}  {tpg / chance:9,.0f} {tpw / chance:9,.0f}")

    print("\nRead gate on k4 T7 with pred: hits per read, recall by bin")
    for gate, name in enumerate(("no gate", "read = column", "not H against E")):
        h_ = sum(hits(4, "pred", c, 7, gate) for c in range(3))
        ok = tbl[4]["pred"][:, :, gate].max(1) >= 7
        print(f"{name:16} {h_:10,.0f}  " + " ".join(f"{ok[b == i].mean():.3f}" for i in range(5)))

    # Policy search: one (k, T) per class, on half the pairs; recall reported on the other half
    rng = np.random.default_rng(77)
    half = rng.random(n) < 0.5
    options = [(k, t / 2) for k in KS for t in range(int(MIN_T[k] * 2), 60)] + [None]

    def passes(policy, p):
        ok = np.zeros(n, dtype=bool)
        for c, o in enumerate(policy):
            if o:
                ok |= tbl[o[0]][p][:, c, 0] >= o[1]
        return ok

    def cost(policy, p):
        return sum(hits(o[0], p, c, o[1]) for c, o in enumerate(policy) if o)

    def size(policy, p):
        return sum(entry(o[0], p, c, o[1]) for c, o in enumerate(policy) if o)

    base_rec = np.array([base[half & (b == i)].mean() for i in range(5)])
    base_cost = sum(hits(4, "pred", c, 7) for c in range(3))
    base_size = sum(entry(4, "pred", c, 7) for c in range(3))

    def search(p, delta, cap, classes):
        policy = [(4, 7.0)] * classes
        need = base_rec - delta
        ok_for = {}
        for c in range(classes):
            for o in options:
                ok_for[c, o] = (tbl[o[0]][p][:, c, 0] >= o[1]) if o else np.zeros(n, dtype=bool)
        for _ in range(6):
            changed = False
            for c in range(classes):
                rest = np.zeros(n, dtype=bool)
                for c2, o in enumerate(policy):
                    if c2 != c:
                        rest |= ok_for[c2, o]
                best_o, best_c = policy[c], cost(policy, p)
                other_size = size(policy, p) - (entry(policy[c][0], p, c, policy[c][1]) if policy[c] else 0)
                other_cost = best_c - (hits(policy[c][0], p, c, policy[c][1]) if policy[c] else 0)
                for o in options:
                    co = other_cost + (hits(o[0], p, c, o[1]) if o else 0)
                    if co >= best_c:
                        continue
                    if cap and other_size + (entry(o[0], p, c, o[1]) if o else 0) > cap * base_size:
                        continue
                    ok = rest | ok_for[c, o]
                    if all(ok[half & (b == i)].mean() >= need[i] for i in range(5)):
                        best_o, best_c = o, co
                if best_o != policy[c]:
                    policy[c] = best_o
                    changed = True
            if not changed:
                break
        return policy

    print("\nBest (k, T) per class; recall on held-out pairs; hits and entries against k4 T7")
    print("policy   delta cap  hits  entries  recall by bin (baseline " + " ".join(f"{base[~half & (b == i)].mean():.3f}" for i in range(5)) + ")  choice")
    for delta in (0.0, 0.005, 0.01, 0.02):
        for cap in (None, 1.5, 4):
            for p, classes in (("pred", 3), ("shuffle", 3), ("info", 3), ("indel", 3), ("predinfo", 9), ("truth", 4)):
                pol = search(p, delta, cap, classes)
                ok = passes(pol, p)
                print(f"{p:8} {delta:.3f} {str(cap):4} {cost(pol, p) / base_cost:.3f}  {size(pol, p) / base_size:6.2f}  "
                      + " ".join(f"{ok[~half & (b == i)].mean():.3f}" for i in range(5)) + "  "
                      + " ".join(f"k{o[0]}T{o[1]:g}" if o else "off" for o in pol), flush=True)
            # One (k, T) everywhere: the control for all of the above
            uni = None
            best = np.inf
            need = base_rec - delta
            all_cls = lambda o: tuple([o] * 3)
            for o in options[:-1]:
                pol = all_cls(o)
                co = cost(pol, "pred")
                if co < best and (not cap or size(pol, "pred") <= cap * base_size):
                    ok = passes(pol, "pred")
                    if all(ok[half & (b == i)].mean() >= need[i] for i in range(5)):
                        uni, best = pol, co
            if uni:
                ok = passes(uni, "pred")
                print(f"{'uniform':8} {delta:.3f} {str(cap):4} {best / base_cost:.3f}  {size(uni, 'pred') / base_size:6.2f}  "
                      + " ".join(f"{ok[~half & (b == i)].mean():.3f}" for i in range(5)) + f"  k{uni[0][0]}T{uni[0][1]:g}", flush=True)


def web(prep, out, data, chunks=200):
    """Per state, 5 bytes: GOR class, DSSP class (3 = no chain), P(H), P(E), P(C) x 255."""
    _, offs = load_prep(prep)
    z = np.load(Path(out) / "classes.npz")
    odds = np.exp(np.clip(z["scores"].astype(np.float64), -30, 30))
    p = odds / (1 + odds)  # each score is the log-odds of one state against the other two
    p /= p.sum(1, keepdims=True)
    rows = np.column_stack([z["pred"], z["truth"], np.round(p * 255)]).astype(np.uint8)
    (Path(data) / "ss").mkdir(exist_ok=True)
    n = len(offs) - 1
    for c in range(chunks):
        a, b = c * n // chunks, (c + 1) * n // chunks
        with gzip.open(Path(data) / "ss" / f"{c}.bin.gz", "wb", compresslevel=9) as f:
            f.write(rows[offs[a]:offs[b]].tobytes())
    print("wrote", chunks, "pieces")


def at_noise(chance, values, target):
    """values at the T where chance falls to target, interpolated in log chance; None if out of range."""
    w = np.flatnonzero(chance > 0)
    if not len(w) or target > chance[w[0]] or target < chance[w[-1]]:
        return None
    i = np.searchsorted(-np.log(chance[w]), -np.log(target), side="right") - 1
    i = min(i, len(w) - 2)
    a, b = w[i], w[i + 1]
    f = (np.log(chance[a]) - np.log(target)) / (np.log(chance[a]) - np.log(chance[b]))
    return values[a] + f * (values[b] - values[a])


def table(out):
    """True-diagonal rate of a k-mer window, by k and class, at the same chance rate per window x column."""
    out = Path(out)
    no = pickle.load(open(out / "noise.pkl", "rb"))
    pa = pickle.load(open(out / "pairs.pkl", "rb"))
    en = pickle.load(open(out / "entries.pkl", "rb"))
    npz = np.load(out / "policies.npz")
    grid = np.arange(TBINS)

    def tail(x):
        return x[::-1].cumsum()[::-1].astype(float)

    print(f"{no['reads']} reads, {no['residues']} residues; windows on the true diagonal, gapless hits")
    for p, title in (("truth", "DSSP"), ("pred", "GOR on emissions")):
        cols = np.bincount(npz[p], minlength=POLICIES[p])
        for target in (2.6e-4, 5e-5, 5e-6):
            print(f"\n{title}, chance {target:.1e} per window x column")
            print("k  class  T     >=GA    -5-GA   -10--5  -20--10  fold over chance -10--5  entries/column  windows -10--5")
            for k in KS:
                for c, name in enumerate("HEC"):
                    positions = no["residues"] - (k - 1) * no["orfs"]
                    chance = tail(no["hist"][k][p][c].sum(0)) / (positions * cols[c])
                    t = at_noise(chance, grid / BIN, target)
                    if t is None:
                        print(f"{k}  {name}      -")
                        continue
                    rate = [at_noise(chance, tail(pa["tp"][k][p][b, c]) / max(pa["tp_n"][k][p][b, c], 1), target) for b in range(4)]
                    entries = at_noise(chance, tail(en["h"][k][p][c]) / max(en["seen"][k][p][c], 1), target)
                    print(f"{k}  {name}  {t:5.1f}  " + "  ".join(f"{r:.3f}" for r in rate)
                          + f"  {rate[2] / target:12,.0f}  {entries:12,.0f}  {pa['tp_n'][k][p][2, c]:8d}")


def check():
    prof = np.full((6, 21), -1.0, dtype=np.float32)
    x = np.array([0, 1, 2, 3], dtype=np.int64)
    for i in range(4):
        prof[i + 2, x[i]] = 5.0  # the read lies on columns 2..5
    a, best, col = diagonals(prof, x)
    d = 2 + len(x) - 1
    assert best.argmax() == d and abs(best.max() - 20) < 1e-5, best
    assert col[d, 0] == 2 and col[d, 3] == 5
    assert (a[col < 0] == -100).all() and (a[col >= 6] == -100).all()
    print("ss_kmers check ok")


if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "check":
        check()
    elif cmd == "classes":
        classes(*args)
    elif cmd == "noise":
        noise(*args)
    elif cmd == "pairs":
        pairs(*args)
    elif cmd == "entries":
        entries(*args)
    elif cmd == "table":
        table(*args)
    elif cmd == "web":
        web(*args)
    elif cmd == "report":
        report(*args)
