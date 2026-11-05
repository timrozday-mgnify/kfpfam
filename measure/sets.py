"""Evaluation sets from MGnify proteins and real reads. Nothing here feeds the k-mer choice.

sets.py fragments DOMTBL[,DOMTBL] PROTEINS.faa OUT PER_FAMILY
    Fragments of GA-passing domains: 50 aa (150 bp) and 33 aa (100 bp), up to PER_FAMILY per
    family and length, so rare families count. Each has a weight restoring the family's abundance.
sets.py pairs FRAGMENTS_META FRAG_DOMTBL READS_ORFS.faa READS_DOMTBL MODELS.tsv OUT
    The (sequence, model) pairs hmmsearch reports at E <= 1e-3 or GA: what a pre-filter must keep.
"""

import csv
import random
import sys
from collections import defaultdict


def fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            if line[0] == ">":
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def domtbl(path):
    with open(path) as f:
        for line in f:
            if line[0] != "#":
                c = line.split()
                yield c[0], c[4].split(".")[0], float(c[6]), float(c[13]), int(c[17]), int(c[18])


def fragments(tables, proteins, out, per_family):
    rng = random.Random(17)
    seqs = dict(fasta(proteins))
    fam = defaultdict(list)
    for table in tables.split(","):
        for protein, acc, _, score, af, at in domtbl(table):
            if protein in seqs:
                fam[acc].append((protein, af, at, score))
    print(len(fam), "families,", sum(map(len, fam.values())), "domains")
    with open(out + ".faa", "w") as fa, open(out + ".tsv", "w") as meta:
        meta.write("id\tacc\tlength\tweight\n")
        n = 0
        for acc in sorted(fam):
            for length, min_overlap in ((50, 30), (33, 25)):
                doms = fam[acc]
                picks = rng.sample(doms, min(per_family, len(doms)))
                for protein, af, at, _ in picks:
                    s = seqs[protein]
                    if len(s) < length:
                        continue
                    lo = max(0, af - 1 - (length - min_overlap))
                    hi = min(len(s) - length, at - min_overlap)
                    if lo > hi:
                        lo, hi = max(0, af - 1), min(len(s) - length, at - 1)
                        if lo > hi:
                            continue
                    st = rng.randint(lo, hi)
                    fa.write(f">F{n}\n{s[st:st + length]}\n")
                    meta.write(f"F{n}\t{acc}\t{length}\t{len(doms) / len(picks):.2f}\n")
                    n += 1
    print(n, "fragments")


def pairs(meta_path, frag_table, frags, orfs, reads_table, models, out, max_reads=3000, per_family=3):
    rng = random.Random(23)
    with open(models) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    index = {r["acc"]: i for i, r in enumerate(rows)}
    ga = [float(r["ga"]) for r in rows]
    with open(meta_path) as f:
        meta = {r["id"]: r for r in csv.DictReader(f, delimiter="\t")}
    seqs = dict(fasta(frags))
    seqs.update(fasta(orfs))

    def hits(table):
        best = {}
        for q, acc, evalue, score, *_ in domtbl(table):
            key = (q, index[acc])
            if key not in best or score > best[key][1]:
                best[key] = (evalue, score)
        return best

    with open(out + ".seqs.tsv", "w") as s, open(out + ".meta.tsv", "w") as o:
        o.write("set\tquery\tmodel\tweight\tevalue\tscore\tga_pass\n")
        by = defaultdict(list)
        for (q, m), (e, sc) in hits(frag_table).items():
            r = meta[q]
            if index[r["acc"]] == m and (e <= 1e-3 or sc >= ga[m]):
                by[(r["length"], m)].append((q, float(r["weight"]), e, sc))
        n = 0
        for (length, m), items in sorted(by.items()):
            rng.shuffle(items)
            kept = items[:per_family]
            scale = len(items) / len(kept)
            for q, w, e, sc in kept:
                s.write(f"{seqs[q]}\t{m}\n")
                o.write(f"f{length}\t{q}\t{m}\t{w * scale:.2f}\t{e:.3g}\t{sc}\t{int(sc >= ga[m])}\n")
                n += 1
        print(n, "fragment pairs,", len({k[1] for k in by}), "families")
        rh = [(k, v) for k, v in hits(reads_table).items() if v[0] <= 1e-3 or v[1] >= ga[k[1]]]
        rng.shuffle(rh)
        print(len(rh), "read pairs available")
        for (q, m), (e, sc) in rh[:max_reads]:
            s.write(f"{seqs[q]}\t{m}\n")
            o.write(f"reads\t{q}\t{m}\t1\t{e:.3g}\t{sc}\t{int(sc >= ga[m])}\n")


if __name__ == "__main__":
    if sys.argv[1] == "fragments":
        fragments(sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]))
    else:
        pairs(*sys.argv[2:9])
