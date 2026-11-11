"""How well GOR predicts secondary structure on fragments and on HMM emissions, without homology.

Models fall in five folds by clan. A fold's tables are trained on PDB chains with no GA domain in
that fold, so every prediction here is of a family the tables never saw.

Usage: ss_accuracy.py SS_DIR PREP_DIR OUT_DIR
SS_DIR holds chains.tsv, reps.faa, reps.domtbl and reps.hs.txt (hmmsearch --cut_ga on reps.faa).
Writes classes.npz, per Pfam state: predicted class, true class where a chain covers the column, and
per column gap counts. Prints the tables.
"""

import multiprocessing as mp
import sys
import zlib
from collections import defaultdict
from pathlib import Path

import numpy as np

import gor

Q = np.array([.0787945, .0151600, .0535222, .0668298, .0397062, .0695071, .0229198, .0590092,
              .0594422, .0963728, .0237718, .0414386, .0482904, .0395639, .0540978, .0683364,
              .0540687, .0673417, .0114135, .0304133])
Q /= Q.sum()
FOLDS = 5


def model_folds(prep):
    rows = [line.rstrip("\n").split("\t") for line in open(prep / "models.tsv")][1:]
    accs = [r[0] for r in rows]
    fold = {r[0]: zlib.crc32((r[2] or r[0]).encode()) % FOLDS for r in rows}
    return accs, fold


def read_chains(ss):
    reps = {line[1:].strip() for line in open(ss / "reps.faa") if line[0] == ">"}
    out = {}
    for line in open(ss / "chains.tsv"):
        name, seq, st = line.rstrip("\n").split("\t")
        if name in reps:
            out[name] = (gor.encode(seq), gor.encode_ss(st))
    return out


def chain_folds(ss, chains, fold):
    home = defaultdict(set)
    for line in open(ss / "reps.domtbl"):
        if line[0] != "#":
            c = line.split()
            home[c[0]].add(fold[c[4].split(".")[0]])
    return {n: home[n] or {zlib.crc32(n.encode()) % FOLDS} for n in chains}


def alignments(path):
    """Yields (model acc, chain, [(column or 0, residue index or -1)]) from hmmsearch's text."""
    acc = name = chain = None
    q, t, qstart, tstart = [], [], 0, 0

    def pairs():
        out, c, r = [], qstart, tstart
        for a, b in zip("".join(q), "".join(t)):
            col, res = 0, -1
            if a != ".":
                col, c = c, c + 1
            if b != "-":
                res, r = r, r + 1
            out.append((col, res))
        return out

    with open(path) as f:
        for line in f:
            if line.startswith(("  == domain", ">> ", "Internal pipeline", "Query:")):
                if q and t:
                    yield acc, chain, pairs()
                q, t = [], []
            if line.startswith("Query:"):
                name = line.split()[1]
            elif line.startswith("Accession:"):
                acc = line.split()[1].split(".")[0]
            elif line.startswith(">> "):
                chain = line.split()[1]
            else:
                c = line.split()
                if len(c) == 4 and c[1].isdigit() and c[3].isdigit():
                    if c[0] == name:
                        if not q:
                            qstart = int(c[1])
                        q.append(c[2])
                    elif c[0] == chain:
                        if not t:
                            tstart = int(c[1]) - 1
                        t.append(c[2])


def fragment_table(models, chains, homes):
    """Q3 by fragment length and by distance from the fragment's end."""
    rng = np.random.default_rng(5)
    print("\nGOR on fragments of held-out chains (no homolog in training)")
    print("length  residues  Q3    Q3+const  H recall E recall C recall  H prec  E prec")
    by_end = defaultdict(lambda: [0, 0])
    for length in (None, 100, 50, 30, 20):
        ok = ok2 = tot = 0
        conf = np.zeros((3, 3), dtype=np.int64)
        for f, m in enumerate(models):
            win, truth, ends = [], [], []
            for n, (x, st) in chains.items():
                if min(homes[n]) != f:
                    continue
                if length is None:
                    win.append(gor.windows(x))
                    truth.append(st)
                    continue
                if len(x) < length:
                    continue
                start = rng.integers(0, length)
                for a in range(start, len(x) - length + 1, length):
                    win.append(gor.windows(x[a:a + length]))
                    truth.append(st[a:a + length])
                    ends.append(np.minimum(np.arange(length), length - 1 - np.arange(length)))
            win, truth = np.concatenate(win), np.concatenate(truth)
            s = m.score_windows(win)
            ok += (s.argmax(1) == truth).sum()
            p = (s + m.const).argmax(1)
            ok2 += (p == truth).sum()
            tot += len(truth)
            conf += np.bincount(truth * 3 + p, minlength=9).reshape(3, 3)
            if length == 30:
                e = np.concatenate(ends)
                for d in range(15):
                    w = e == d
                    by_end[d][0] += (p[w] == truth[w]).sum()
                    by_end[d][1] += w.sum()
        rec = conf.diagonal() / conf.sum(1)
        pre = conf.diagonal() / conf.sum(0)
        print(f"{length or 'full':>6} {tot:9d} {ok / tot:.3f}  {ok2 / tot:.3f}   "
              f"{rec[0]:.3f}    {rec[1]:.3f}    {rec[2]:.3f}    {pre[0]:.3f}  {pre[1]:.3f}")
    print("\n30 aa fragments, Q3 by distance from the nearer end:")
    print(" ".join(f"{d}:{a / b:.3f}" for d, (a, b) in sorted(by_end.items())))


def predict_fold(args):
    m, profile_path, ranges = args
    lo = np.memmap(profile_path, dtype=np.float32, mode="r").reshape(-1, 21)
    out = []
    for a, b in ranges:
        p = Q * np.exp2(lo[a:b, :20].astype(float))
        p /= p.sum(1, keepdims=True)
        out.append(m.score_profile(p))
    return out


def main():
    ss, prep, out = (Path(a) for a in sys.argv[1:4])
    out.mkdir(parents=True, exist_ok=True)
    accs, fold = model_folds(prep)
    chains = read_chains(ss)
    homes = chain_folds(ss, chains, fold)
    print(len(chains), "chains,", sum(len(x) for x, _ in chains.values()), "residues")
    models = []
    rng = np.random.default_rng(3)
    for f in range(FOLDS):
        train = [c for n, c in chains.items() if f not in homes[n]]
        m = gor.Gor(train)
        pick = rng.choice(len(train), min(1500, len(train)), replace=False)
        s = np.concatenate([m.score_seq(train[i][0]) for i in pick])
        t = np.concatenate([train[i][1] for i in pick])
        m.const = gor.fit_constants(s, t)
        models.append(m)
        print(f"fold {f}: {len(train)} training chains, constants H {m.const[0]:.1f} E {m.const[1]:.1f}", flush=True)
    fragment_table(models, chains, homes)

    offs = np.fromfile(prep / "offs.u32", dtype=np.uint32).astype(np.int64)
    index = {a: i for i, a in enumerate(accs)}
    n = int(offs[-1])
    # Truth per state, from the chains' alignments; gap events per state too
    counts = np.zeros((n, 3), dtype=np.int32)
    aligned = np.zeros(n, dtype=np.int32)
    gaps = np.zeros(n, dtype=np.int32)
    domains = bad = 0
    for acc, chain, pairs in alignments(ss / "reps.hs.txt"):
        if chain not in chains or acc not in index:
            continue
        st = chains[chain][1]
        base = offs[index[acc]] - 1
        domains += 1
        last_gap = False
        last_col = None
        for col, res in pairs:
            if res >= len(st):
                bad += 1
                break
            if col and res >= 0:
                counts[base + col, st[res]] += 1
                aligned[base + col] += 1
                last_gap = False
                last_col = col
            elif not last_gap and last_col is not None:
                gaps[base + last_col] += 1
                last_gap = True
    print(f"\n{domains} GA domains parsed from chains, {bad} failed bounds, {(aligned > 0).sum()} of {n} states covered")
    truth = np.where(aligned > 0, counts.argmax(1), 3).astype(np.int8)

    # Emissions, each model by the fold that never saw its clan
    ranges = [[] for _ in range(FOLDS)]
    which = [[] for _ in range(FOLDS)]
    for i, a in enumerate(accs):
        ranges[fold[a]].append((int(offs[i]), int(offs[i + 1])))
        which[fold[a]].append(i)
    with mp.Pool(4) as pool:
        done = pool.map(predict_fold, [(models[f], str(prep / "profile.f32"), ranges[f]) for f in range(FOLDS)])
    scores = np.zeros((n, 3), dtype=np.float32)
    for f in range(FOLDS):
        for (a, b), s in zip(ranges[f], done[f]):
            scores[a:b] = s
    # Decision constants for emissions, from the other folds' columns
    state_fold = np.zeros(n, dtype=np.int8)
    for f in range(FOLDS):
        for a, b in ranges[f]:
            state_fold[a:b] = f
    pred = np.zeros(n, dtype=np.int8)
    known = truth < 3
    print("\nEmission prediction against the chain's true state at the column")
    print("fold  columns  Q3 raw  Q3 const  constants")
    for f in range(FOLDS):
        fit = known & (state_fold != f)
        c = gor.fit_constants(scores[fit], truth[fit])
        here = state_fold == f
        pred[here] = (scores[here] + c).argmax(1)
        k = known & here
        print(f"{f}  {k.sum():8d}  {(scores[k].argmax(1) == truth[k]).mean():.3f}  {(pred[k] == truth[k]).mean():.3f}  H {c[0]:.1f} E {c[1]:.1f}", flush=True)
    k = known
    conf = np.bincount(truth[k] * 3 + pred[k], minlength=9).reshape(3, 3)
    print("all", k.sum(), f"Q3 {(pred[k] == truth[k]).mean():.3f}")
    print("truth rows H E C, predicted columns H E C:\n", conf)
    print("share of all Pfam columns predicted H E C:", np.round(np.bincount(pred, minlength=3) / n, 3))
    print("share of covered columns truly H E C:", np.round(np.bincount(truth[k], minlength=3) / k.sum(), 3))
    # Same chains, sequence alone, within the same columns: does the profile beat the sequence?
    np.savez_compressed(out / "classes.npz", pred=pred, truth=truth, scores=scores, aligned=aligned, gaps=gaps, fold=state_fold)
    import pickle
    with open(out / "gor_tables.pkl", "wb") as o:
        pickle.dump(models, o)
    print("saved", out)


if __name__ == "__main__":
    main()
