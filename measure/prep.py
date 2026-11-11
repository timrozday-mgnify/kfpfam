"""Pfam-A.hmm.gz + clans to the arrays and model table every later step reads.

Usage: prep.py Pfam-A.hmm.gz Pfam-A.clans.tsv.gz outdir
Writes profile.f32 (states x 21, log-odds bits, last column for unknown residues), offs.u32,
accs.txt and models.tsv.
"""

import gzip
import sys
from pathlib import Path

import numpy as np

Q = np.array([.0787945, .0151600, .0535222, .0668298, .0397062, .0695071, .0229198, .0590092,
              .0594422, .0963728, .0237718, .0414386, .0482904, .0395639, .0540978, .0683364,
              .0540687, .0673417, .0114135, .0304133])
Q /= Q.sum()


def parse(path):
    meta, rows, offs, m = [], [], [0], {}
    with gzip.open(path, "rt") as f:
        for line in f:
            key = line[:5].strip()
            if key in ("NAME", "ACC", "LENG", "EFFN", "NSEQ", "GA"):
                m[key] = line.split()[1]
            elif line.startswith("STATS LOCAL MSV"):
                m["MU"], m["LAMBDA"] = line.split()[3:5]
            elif key == "HMM":
                for _ in range(4):
                    next(f)  # transition header, COMPO and the begin-state lines
                for _ in range(int(m["LENG"])):
                    rows.append(np.array(next(f).split()[1:21], dtype=float))
                    next(f), next(f)
                offs.append(len(rows))
                meta.append(m)
                m = {}
    return np.array(rows), np.array(offs, dtype=np.uint32), meta


def main():
    hmm, clans, out = sys.argv[1:4]
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    emit, offs, meta = parse(hmm)
    p = np.exp(-emit)
    lo = np.log2(p / Q)
    rel = (p * lo).sum(1)
    table = np.hstack([lo, np.full((len(lo), 1), -100.0)]).astype(np.float32)
    table.tofile(out / "profile.f32")
    offs.tofile(out / "offs.u32")
    clan = {}
    with gzip.open(clans, "rt") as f:
        for line in f:
            c = line.split("\t")
            clan[c[0]] = c[1]
    with open(out / "models.tsv", "w") as o, open(out / "accs.txt", "w") as a:
        o.write("acc\tname\tclan\tlength\tga\tnseq\teffn\tmu\tlambda\tmean_info\ttop_info\t"
                "conserved\n")
        for i, m in enumerate(meta):
            acc = m["ACC"].split(".")[0]
            r = rel[offs[i]:offs[i + 1]]
            top = np.convolve(r, np.ones(5), "valid").max() if len(r) >= 5 else r.sum()
            a.write(acc + "\n")
            o.write(f"{acc}\t{m['NAME']}\t{clan.get(acc, '')}\t{m['LENG']}\t{m['GA']}\t{m['NSEQ']}"
                    f"\t{m['EFFN']}\t{m['MU']}\t{m['LAMBDA']}\t{r.mean():.3f}\t{top:.2f}"
                    f"\t{(r >= 1).mean():.3f}\n")
    print(len(meta), "models,", len(table), "states")


if __name__ == "__main__":
    main()
