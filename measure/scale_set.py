"""ORFs shaped like translated reads, for timing at scale: per read, one true-frame window of a
protein and the rest shuffled residues, with the real read set's ORF count and lengths.

Usage: scale_set.py PROTEINS.faa READS_ORFS.faa NREADS OUT.faa
"""

import random
import sys
from collections import defaultdict

from sets import fasta

proteins, orfs, n, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
rng = random.Random(31)
long_enough = [s for _, s in fasta(proteins) if len(s) >= 50]
shapes = defaultdict(list)
for name, s in fasta(orfs):
    shapes[name.split("_")[0]].append(len(s))
shapes = list(shapes.values())
with open(out, "w") as o:
    for r in range(n):
        lengths = rng.choice(shapes)
        for j, length in enumerate(lengths):
            s = rng.choice(long_enough)
            st = rng.randrange(max(1, len(s) - length))
            w = s[st:st + length]
            if j:
                w = "".join(rng.sample(w, len(w)))
            o.write(f">R{r}_{j}\n{w}\n")
