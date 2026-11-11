"""Known-origin weak hits from real reads: windows of w aa cut from ORFs called at GA, centred at a random
point of the aligned region, plus as many windows from the reversed ORFs; pyhmmer against the confirmed
targets only (Z = 20 M, E <= 1000). Writes window, origin models (the ORF's GA calls), model, score, E.

Usage (pyhmmer container): recruit_trunc.py PFAM.h3m REAL.domtbl ORFS.faa DECOY.faa MODELS.tsv CPUS OUT.tsv
"""

import csv
import random
import sys

import pyhmmer

hmm_path, domtbl, orfs, decoys, models, cpus, out = sys.argv[1:]
ga = {r["acc"]: float(r["ga"]) for r in csv.DictReader(open(models), delimiter="\t")}
calls = {}
for line in open(domtbl):
    if line[0] != "#":
        c = line.split()
        a = c[4].split(".")[0]
        if a in ga and float(c[13]) >= ga[a]:
            calls.setdefault(c[0], {}).setdefault(a, (int(c[17]), int(c[18])))


def fasta(path, keep):
    seqs, name = {}, None
    for line in open(path):
        if line[0] == ">":
            name = line[1:].split()[0]
        elif name in keep:
            seqs[name] = seqs.get(name, "") + line.strip()
    return seqs


real, rev = fasta(orfs, calls), fasta(decoys, calls)
rng = random.Random(1)
alphabet = pyhmmer.easel.Alphabet.amino()
seqs, origin = [], {}
for w in (10, 15, 20, 25):
    for orf, hits in sorted(calls.items()):
        s = real[orf]
        for a, (f, t) in sorted(hits.items()):
            c = rng.randint(f, t) - 1
            st = max(0, min(c - w // 2, len(s) - w))
            name = f"T{w}_{orf}_{a}"
            seqs.append(pyhmmer.easel.TextSequence(name=name, sequence=s[st:st + w]))
            origin[name] = ",".join(sorted(hits))
        if orf in rev and len(rev[orf]) >= w:
            st = rng.randint(0, len(rev[orf]) - w)
            name = f"D{w}_{orf}"
            seqs.append(pyhmmer.easel.TextSequence(name=name, sequence=rev[orf][st:st + w]))
            origin[name] = "-"
block = pyhmmer.easel.TextSequenceBlock(seqs).digitize(alphabet)
conf = {a for h in calls.values() for a in h}
with pyhmmer.plan7.HMMFile(hmm_path) as f:
    hmms = [h for h in f if str(h.accession).split(".")[0] in conf]
print(len(seqs), "windows,", len(hmms), "models", flush=True)
with open(out, "w") as o:
    for hits in pyhmmer.hmmsearch(hmms, block, cpus=int(cpus), Z=20_000_000, E=1000, domE=1000):
        acc = str(hits.query.accession).split(".")[0]
        for h in hits:
            n = h.name.decode() if isinstance(h.name, bytes) else h.name
            o.write(f"{n}\t{origin[n]}\t{acc}\t{max(d.score for d in h.domains):.2f}\t{h.evalue:.3g}\n")
