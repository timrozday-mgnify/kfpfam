"""The pipeline's chosen ORF rule, as a FASTA for the filter: plain ORFs in FGS's frame, plus other
frames' ORFs with hexamer coding score >= MIN (tables trained on one half of the reads, applied to the
other, as orf_filters.py frames does). The same names are kept from the reversed ORFs, the decoys.

Usage: orf_rule.py READS.fasta CALLS.faa ORFS.faa DECOY.faa MIN OUT_PREFIX
"""

import sys
from collections import Counter, defaultdict

from orf_filters import fasta, hexamers, log_odds, placed
from real_orfs import RC

reads_fa, calls_faa, orfs_faa, decoy_faa, low, out = sys.argv[1:]
low = float(low)
orfs = {}
for r, s in fasta(reads_fa):
    orfs[r] = [(f"{r}_{j}", a, strand, s[a:b] if not strand else s[a:b][::-1].translate(RC))
               for j, (_, a, b, strand) in enumerate(placed(s))]
called = defaultdict(set)
for name, _ in fasta(calls_faa):
    r, a, b, strand = name.rsplit("_", 3)
    called[r] |= {(strand == "-", (int(a) - 1) % 3), (strand == "-", int(b) % 3)}
fold = lambda r: int(r[1:]) % 2  # noqa: E731
coding, other = [Counter(), Counter()], [Counter(), Counter()]
for r, frames in called.items():
    for _, a, st, nt in orfs[r]:
        (coding if (bool(st), a % 3) in frames else other)[fold(r)].update(hexamers(nt))
tables = [log_odds(coding[1 - i], other[1 - i]) for i in (0, 1)]
keep = set()
for r, rows in orfs.items():
    t = tables[fold(r)]
    for name, a, st, nt in rows:
        if (bool(st), a % 3) in called.get(r, ()) or sum(t.get(h, 0.0) for h in hexamers(nt)) >= low:
            keep.add(name)
for src, dst in ((orfs_faa, out + ".faa"), (decoy_faa, out + ".decoy.faa")):
    n = res = 0
    with open(dst, "w") as f:
        for name, seq in fasta(src):
            if name in keep:
                f.write(f">{name}\n{seq}\n")
                n += 1
                res += len(seq)
    print(f"{dst}: {n} ORFs, {res / 50000:.1f} residues per read")
