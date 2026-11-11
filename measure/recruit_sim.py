"""Simulated metagenome with known read origin, for recruitment (docs/recruit-plan.md, evaluation b).

annotate (pyhmmer container): N random MGnify proteins, their Pfam domains at E <= 1 (score, GA flag,
  envelope).  recruit_sim.py annotate PFAM.h3m MODELS.tsv PROTEINS.faa N CPUS OUT_PREFIX
reads (local): proteins chained into contigs of 2-6 genes on random strands with random spacers; each
  contig gets an identity (1.0 to 0.5, BLOSUM62 substitutions) and a coverage (log-uniform 0.3-30), is
  back-translated with random codons, and sampled as 100 or 150 bp reads from either strand; six-frame
  ORFs as measure/real_orfs.py.  recruit_sim.py reads OUT_PREFIX
Writes OUT_PREFIX.orfs.faa (R<read>_<orf>), OUT_PREFIX.coords.tsv and OUT_PREFIX.truth.tsv, one row per
gene segment of a read (>= 30 bp): read, contig, identity, coverage, protein, the segment on the forward
read (bp), domains at GA overlapping it by >= 10 aa, other domains overlapping.
"""

import math
import random
import sys

from real_orfs import CODON, RC, orfs_at

IDENTITIES = (1.0, 0.9, 0.8, 0.7, 0.6, 0.5)
AAS = "ACDEFGHIKLMNPQRSTVWY"


def fasta_sample(path, n, seed=11):
    names = [line[1:].split()[0] for line in open(path) if line[0] == ">"]
    keep = set(random.Random(seed).sample(names, n))
    out, name = {}, None
    for line in open(path):
        if line[0] == ">":
            name = line[1:].split()[0]
        elif name in keep:
            out[name] = out.get(name, "") + line.strip()
    return out


def annotate(hmm_path, models, proteins, n, cpus, prefix):
    import pyhmmer

    ga = {line.split("\t")[0]: float(line.split("\t")[4]) for line in open(models) if not line.startswith("acc\t")}
    seqs = fasta_sample(proteins, int(n))
    with open(f"{prefix}.proteins.faa", "w") as o:
        for k, s in seqs.items():
            o.write(f">{k}\n{s}\n")
    alphabet = pyhmmer.easel.Alphabet.amino()
    block = pyhmmer.easel.TextSequenceBlock(
        pyhmmer.easel.TextSequence(name=k, sequence=s) for k, s in seqs.items()).digitize(alphabet)
    with pyhmmer.plan7.HMMFile(hmm_path) as f, open(f"{prefix}.domains.tsv", "w") as o:
        o.write("protein\tacc\tscore\tga\tenv_from\tenv_to\n")
        hmms = list(f)
        for hits in pyhmmer.hmmsearch(hmms, block, cpus=int(cpus), Z=20_000_000, E=1, domE=1):
            acc = hits.query.accession
            acc = (acc.decode() if isinstance(acc, bytes) else acc).split(".")[0]
            for h in hits:
                name = h.name.decode() if isinstance(h.name, bytes) else h.name
                for d in h.domains:
                    o.write(f"{name}\t{acc}\t{d.score:.2f}\t{int(d.score >= ga.get(acc, 1e9))}\t{d.env_from}\t{d.env_to}\n")


def blosum_conditional():
    """P(b | a, b != a) from BLOSUM62 (half-bit scores, so p(a, b) = q_a q_b 2^(s / 2))."""
    from Bio.Align import substitution_matrices

    m = substitution_matrices.load("BLOSUM62")
    out = {}
    for a in AAS:
        w = [(b, 2 ** (m[a][b] / 2) * 1.0) for b in AAS if b != a]
        tot = sum(x for _, x in w)
        out[a] = ([b for b, _ in w], [x / tot for _, x in w])
    return out


def mutate(seq, ident, cond, rng):
    return "".join(c if c in cond and rng.random() < ident else
                   (rng.choices(*cond[c])[0] if c in cond else rng.choice(AAS)) for c in seq)


def reads(prefix, seed=7):
    """Contigs of 2-6 neighbouring genes, each on a random strand with 0-150 random bases between them.
    A contig has one identity and a coverage (log-uniform 0.3-30); reads of 100 or 150 bp from either
    strand. Truth is per gene segment of a read: the read stretch it covers and the domains there."""
    rng = random.Random(seed)
    cond = blosum_conditional()
    back = {}
    for codon, aa in CODON.items():
        if aa != "*":
            back.setdefault(aa, []).append(codon)
    seqs, name = {}, None
    for line in open(f"{prefix}.proteins.faa"):
        if line[0] == ">":
            name = line[1:].strip()
        else:
            seqs[name] = line.strip()
    doms = {}
    for line in list(open(f"{prefix}.domains.tsv"))[1:]:
        p, a, s, g, f, t = line.split()
        doms.setdefault(p, []).append((a, int(g), int(f) - 1, int(t)))
    names = sorted(seqs)
    rng.shuffle(names)
    contigs = []
    while names:
        k = rng.randint(2, 6)
        contigs.append(names[:k])
        names = names[k:]
    r = 0
    with open(f"{prefix}.orfs.faa", "w") as fo, open(f"{prefix}.truth.tsv", "w") as ft, \
            open(f"{prefix}.coords.tsv", "w") as fc:
        ft.write("read\tcontig\tidentity\tcoverage\tprotein\tseg_from\tseg_to\tga\tother\n")
        fc.write("orf\tstrand\tstart\tread_len\n")
        for ci, genes in enumerate(contigs):
            ident = rng.choice(IDENTITIES)
            cov = math.exp(rng.uniform(math.log(0.3), math.log(30)))
            dna, at = "", []
            for g in genes:
                dna += "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 150)))
                cds = "".join(rng.choice(back[c]) for c in mutate(seqs[g], ident, cond, rng)) + "TAA"
                strand = rng.choice("+-")
                at.append((g, strand, len(dna), len(dna) + len(cds) - 3))  # the stop codon left out
                dna += cds if strand == "+" else cds[::-1].translate(RC)
            dna += "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 150)))
            n_reads = max(1, round(cov * len(dna) / 125))
            for _ in range(n_reads):
                n = rng.choice((100, 150))
                if len(dna) < n:
                    continue
                r0 = rng.randrange(0, len(dna) - n + 1)
                r1 = r0 + n
                read = dna[r0:r1]
                minus = rng.random() < 0.5
                if minus:
                    read = read[::-1].translate(RC)
                rows = []
                for g, strand, g0, g1 in at:
                    lo, hi = max(r0, g0), min(r1, g1)
                    if hi - lo < 30:
                        continue
                    if strand == "+":
                        a0, a1 = (lo - g0) // 3, (hi - g0) // 3
                    else:  # the protein starts at the gene's right end
                        a0, a1 = (g1 - hi) // 3, (g1 - lo) // 3
                    ga = sorted({a for a, gg, f, t in doms.get(g, []) if gg and min(t, a1) - max(f, a0) >= 10})
                    other = sorted({a for a, gg, f, t in doms.get(g, []) if min(t, a1) - max(f, a0) >= 10} - set(ga))
                    s0, s1 = (r1 - hi, r1 - lo) if minus else (lo - r0, hi - r0)
                    rows.append((g, s0, s1, ",".join(ga), ",".join(other)))
                for g, s0, s1, ga, other in rows or [("-", 0, n, "", "")]:
                    ft.write(f"R{r}\tC{ci}\t{ident}\t{cov:.3f}\t{g}\t{s0}\t{s1}\t{ga}\t{other}\n")
                for j, (o, sign, start) in enumerate(orfs_at(read)):
                    fo.write(f">R{r}_{j}\n{o}\n")
                    fc.write(f"R{r}_{j}\t{sign}\t{start}\t{n}\n")
                r += 1
    print(len(contigs), "contigs,", r, "reads")


if __name__ == "__main__":
    if sys.argv[1] == "annotate":
        annotate(*sys.argv[2:])
    else:
        reads(sys.argv[2])
