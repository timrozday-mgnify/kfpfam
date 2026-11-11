"""HMMER on the filter's candidates, writing every pair recruitment needs (runs in the pyhmmer container).

Per (ORF, model) the strongest domain: ORF, accession, model length M, ORF length L, domain score,
sequence E-value, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to, c_evalue, i_evalue, cigar
(pfam-reannotate's convention: M match, I sequence only, D model only). Z = 20 M and E <= 1000, as the
truth was run. domZ is HMMER's default, the hits a search reports, so c_evalue depends on which ORFs a
model was searched with, including the chunk (scores, evalue and i_evalue do not; tested by chunking).

Candidates are read in chunks of whole filter blocks (CHUNK pairs or more), with their ORFs read forward
from ORFS.faa, so memory does not grow with the sample: all at once held about 3.5 GB per million reads.

Usage: pairs.py PFAM.h3m MODELS.tsv ORFS.faa CANDIDATES BAR CPUS OUT.tsv
  CANDIDATES: the filter's TSV (ORF, model index, diagonal, gapped), pairs with gapped >= BAR searched,
              in `filter join`'s order against the same ORFS.faa (one line per sequence);
              or 'all:ACCS.txt' to search every ORF against the models listed (one accession per line).
"""

import gc
import multiprocessing as mp
import re
import sys
import threading
from collections import defaultdict

import pyhmmer

Z = 20_000_000
CHUNK = 8_000_000  # pairs per chunk: ~1.5 GB; each chunk rebuilds every model's profile (~2 core-s)
OP = re.compile(r"(\d+)([MID])")
COLS = ("orf", "acc", "M", "L", "score", "evalue", "hmm_from", "hmm_to", "ali_from", "ali_to",
        "env_from", "env_to", "c_evalue", "i_evalue", "cigar")


def fasta(path, keep=None):
    out, name = {}, None
    for line in open(path):
        if line[0] == ">":
            name = line[1:].split()[0]
        elif keep is None or name in keep:
            out[name] = out.get(name, "") + line.strip()
    return out


def chunks(path, accs, bar, size):
    """The filter's candidates, {acc: ORF names} with gapped >= bar, in runs of whole blocks of at least size
    pairs. `filter join` writes a block's pairs in model order, so a new block starts where the model index
    falls; a fall missed only merges two blocks."""
    cands, n, last = defaultdict(set), 0, -1
    for line in open(path):
        q, m, _, g = line.split()
        m = int(m)
        if m < last and n >= size:
            yield cands
            cands, n = defaultdict(set), 0
        last = m
        if float(g) >= bar:
            cands[accs[m]].add(q)
            n += 1
    if n:
        yield cands


class Ahead:
    """A one-line-per-sequence FASTA read forward: take(names) gives {name: sequence}, reading on from where
    the last call stopped (candidates come in the FASTA's order, block by block)."""

    def __init__(self, path):
        self.path, self.f, self.name = path, open(path), None

    def take(self, names):
        out = {}
        for line in self.f:
            if line[0] == ">":
                self.name = line[1:].split()[0]
            elif self.name in names:
                out[self.name] = line.strip()
                if len(out) == len(names):
                    return out
        missing = names - out.keys()
        raise SystemExit(f"{len(missing)} candidate ORFs (e.g. {next(iter(missing))}) not found ahead in {self.path}: "
                         "candidates must be in filter join's order, against the FASTA it read")


def acc_of(hmm):
    a = hmm.accession
    return (a.decode() if isinstance(a, bytes) else a).split(".")[0]


def cigar(model, target):
    """Run-length operations over the alignment, as pfam-reannotate's reannotate.py writes them: M match,
    I residue in the sequence only (model line '.'), D model position only (sequence line '-')."""
    ops = []
    for m, t in zip(model, target):
        op = "I" if m == "." else "D" if t == "-" else "M"
        if ops and ops[-1][0] == op:
            ops[-1][1] += 1
        else:
            ops.append([op, 1])
    return "".join(f"{n}{op}" for op, n in ops)


def rows(hmm, hits, lengths):
    acc = acc_of(hmm)
    for h in hits:
        if not h.domains:
            continue
        d = max(h.domains, key=lambda d: d.score)
        n = h.name.decode() if isinstance(h.name, bytes) else h.name
        a = d.alignment
        ops = cigar(a.hmm_sequence, a.target_sequence)
        used = {op: sum(int(x) for x, o in OP.findall(ops) if o in op) for op in ("MD", "MI")}
        assert (used["MD"], used["MI"]) == (a.hmm_to - a.hmm_from + 1, a.target_to - a.target_from + 1), (n, acc, ops)
        yield (n, acc, hmm.M, lengths[n], f"{d.score:.2f}", f"{h.evalue:.3g}",
               a.hmm_from, a.hmm_to, a.target_from, a.target_to, d.env_from, d.env_to,
               f"{d.c_evalue:.3g}", f"{d.i_evalue:.3g}", ops)


_JOB = {}


def _work(part):
    hmms, keep, pipe_args = _JOB["hmms"], _JOB["keep"], _JOB["pipe_args"]
    alphabet = pyhmmer.easel.Alphabet.amino()
    pipe = pyhmmer.plan7.Pipeline(alphabet, Z=Z, **pipe_args)
    bg = pyhmmer.plan7.Background(alphabet)
    out = []
    for key, orfs in part:
        hmm = hmms[key]
        prof = pyhmmer.plan7.Profile(hmm.M, alphabet)
        prof.configure(hmm, bg)
        block = pyhmmer.easel.TextSequenceBlock(
            pyhmmer.easel.TextSequence(name=q.encode(), sequence=s) for q, s in orfs).digitize(alphabet)
        out.extend(keep(hmm, pipe.search_hmm(prof.to_optimized(), block)))
        pipe.clear()
    return out


def search_each(hmms, cands, seqs, cpus, keep, **pipe_args):
    """Each model of hmms ({key: HMM}) on its own candidates (cands, {key: ORF names}), yielding what
    keep(hmm, hits) returns. Fork workers, not threads: building blocks in Python holds the interpreter
    lock, which tripled the core time at 6 threads (docs/plan.md, 2026-11-11). Biggest jobs are dealt out
    first and interleaved, so no worker is left with the tail."""
    _JOB.update(hmms=hmms, keep=keep, pipe_args=pipe_args or {"E": 1000, "domE": 1000})
    order = sorted(cands, key=lambda k: -len(cands[k]) * hmms[k].M)
    n = cpus * 32
    # each part's ORFs travel with it: a worker reading the inherited dict would copy most of its pages
    # (reference counts) while using a sixth of it, about 2 GB per worker at 1 M reads
    # The pool reads its task iterator ahead without limit, so at most 4 parts per worker are built and
    # queued at once (without the cap, 10 GB at 1 M reads and 19 candidates per read)
    room = threading.BoundedSemaphore(4 * cpus)

    def parts():
        for i in range(n):
            room.acquire()
            yield [(k, [(q, seqs[q]) for q in sorted(cands[k])]) for k in order[i::n]]

    gc.freeze()  # nor should workers' collector walk, and so copy, inherited objects
    with mp.get_context("fork").Pool(cpus) as pool:
        for part in pool.imap(_work, parts()):
            room.release()
            yield from part
    _JOB.clear()


def main():
    hmm_path, models, faa, cands_arg, bar, cpus, out = sys.argv[1:]
    cpus = int(cpus)
    accs = [line.split("\t")[0] for line in open(models)][1:]
    alphabet = pyhmmer.easel.Alphabet.amino()
    if cands_arg.startswith("all:"):
        need = {line.strip() for line in open(cands_arg[4:]) if line.strip()}
        seqs = fasta(faa)
    else:
        need = set(accs)
    with pyhmmer.plan7.HMMFile(hmm_path) as f:
        hmms = {acc_of(h): h for h in f if acc_of(h) in need}
    with open(out, "w") as o:
        o.write("\t".join(COLS) + "\n")
        if cands_arg.startswith("all:"):
            lengths = {n: len(s) for n, s in seqs.items()}
            print(f"{len(seqs)} ORFs, {len(hmms)} models", flush=True)
            block = pyhmmer.easel.TextSequenceBlock(
                pyhmmer.easel.TextSequence(name=n, sequence=s) for n, s in seqs.items()).digitize(alphabet)
            for hits in pyhmmer.hmmsearch(list(hmms.values()), block, cpus=cpus, Z=Z, E=1000, domE=1000):
                for r in rows(hits.query, hits, lengths):
                    o.write("\t".join(map(str, r)) + "\n")
            return
        ahead, n_orfs, used = Ahead(faa), 0, set()
        for cands in chunks(cands_arg, accs, float(bar), CHUNK):
            seqs = ahead.take({q for v in cands.values() for q in v})
            lengths = {n: len(s) for n, s in seqs.items()}
            n_orfs += len(seqs)
            used |= cands.keys()
            for r in search_each(hmms, cands, seqs, cpus, lambda hmm, hits, L=lengths: list(rows(hmm, hits, L))):
                o.write("\t".join(map(str, r)) + "\n")
            del cands, seqs, lengths  # before the next chunk is read, not after
        print(f"{n_orfs} ORFs, {len(used)} models", flush=True)


if __name__ == "__main__":
    main()
