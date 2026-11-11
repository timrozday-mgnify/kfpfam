"""Six-frame ORFs (split at stops, >= 15 aa) of real reads, in a seeded random order, so the first
N reads of any run are a random subset and their names match across subset and whole-run FASTAs.

Usage: real_orfs.py OUT_PREFIX SUBSET R1.fastq.gz [R2.fastq.gz ...]
Writes OUT_PREFIX.subset.faa, OUT_PREFIX.all.faa, OUT_PREFIX.coords.tsv (each ORF's strand, first base on
its strand and the read length: recruit.py places alignments in the read with it) and prints per-read figures. Names are R<read>_<orf>.
"""

import gzip
import itertools
import random
import sys

BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {"".join(c): AA[i] for i, c in enumerate(itertools.product(BASES, repeat=3))}
RC = str.maketrans("ACGTN", "TGCAN")


def orfs_at(seq, min_len=15):
    """[(protein, strand, start)]: start is the ORF's first base on its own strand (0-based), so on the
    forward read the ORF covers start..start + 3 len for '+' and len(seq) - start - 3 len..len(seq) - start
    for '-'."""
    out = []
    for sign, strand in (("+", seq), ("-", seq[::-1].translate(RC))):
        for frame in range(3):
            prot = "".join(CODON.get(strand[i:i + 3], "X") for i in range(frame, len(strand) - 2, 3))
            at = 0
            for p in prot.split("*"):
                if len(p) >= min_len:
                    out.append((p, sign, frame + 3 * at))
                at += len(p) + 1
    return out


def orfs(seq, min_len=15):
    return [p for p, _, _ in orfs_at(seq, min_len)]


def write_coords(f, name, sign, start, read_len):
    f.write(f"{name}\t{sign}\t{start}\t{read_len}\n")


def main():
    prefix, subset = sys.argv[1], int(sys.argv[2])
    seqs = []
    for path in sys.argv[3:]:
        with gzip.open(path, "rt") as f:
            seqs += [line.strip() for i, line in enumerate(f) if i % 4 == 1]
    random.Random(5).shuffle(seqs)
    n_orfs = residues = 0
    with open(f"{prefix}.subset.faa", "w") as a, open(f"{prefix}.all.faa", "w") as b, \
            open(f"{prefix}.coords.tsv", "w") as c:
        c.write("orf\tstrand\tstart\tread_len\n")
        for r, s in enumerate(seqs):
            for j, (p, sign, start) in enumerate(orfs_at(s)):
                line = f">R{r}_{j}\n{p}\n"
                b.write(line)
                write_coords(c, f"R{r}_{j}", sign, start, len(s))
                if r < subset:
                    a.write(line)
                    n_orfs += 1
                    residues += len(p)
    print(len(seqs), "reads;", n_orfs / subset, "ORFs and", residues / subset, "residues per read")


if __name__ == "__main__":
    main()
