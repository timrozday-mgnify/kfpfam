"""V1 of the pipeline plan: hmmer_pairs.py on every (model, ORF) pair against hmmsearch --cut_ga.

Usage: pairs_v1.py ORFS.faa NREADS MODELS.tsv PIPELINE_DIR OUT_DIR
Compares every domtbl column, the coverage profile and the alignments in Parquet.
Run hmmsearch first if OUT_DIR/hmmsearch.domtbl is missing, with the printed command.
"""

import csv
import gzip
import random
import re
import subprocess
import sys
import time
from array import array
from collections import Counter
from pathlib import Path

import pyarrow.parquet as pq

orfs, nreads, models, pipe, out = sys.argv[1:6]
nreads, out, pipe = int(nreads), Path(out), Path(pipe)
out.mkdir(exist_ok=True, parents=True)
sys.path.insert(0, str(pipe / "modules/local/hmmer_pairs/bin"))
import hmmer_pairs as hp  # noqa: E402

HMM = str(Path(orfs).parents[3] / "pfam/Pfam-A.hmm")
Z = 30134
sample = out / "orfs.faa"
if not sample.exists():
    name = lambda line: line[1:].split("_frame=")[0]  # noqa: E731
    with open(orfs) as f:
        reads = sorted({name(line) for line in f if line[0] == ">"})
    keep = set(random.Random(31).sample(reads, nreads))
    print("reads", len(reads), flush=True)
    with open(orfs) as f, open(sample, "w") as o:
        on = False
        for line in f:
            if line[0] == ">":
                on = name(line) in keep
            if on:
                o.write(line)
print("ORFs in sample", sum(1 for line in open(sample) if line[0] == ">"), flush=True)
if not (out / "hmmsearch.domtbl").exists():
    sys.exit(f"run hmmsearch -Z {Z} --cut_ga --notextw -o hmmsearch.txt --domtblout hmmsearch.domtbl Pfam-A.hmm orfs.faa")

ids = {}
with open(sample) as f:
    for line in f:
        if line[0] == ">":
            ids[line[1:].split()[0]] = len(ids)
seqs = hp.read_orfs(sample, ids)
every = array("I", range(len(seqs)))
with open(models) as f:
    pairs = {r["name"]: every for r in csv.DictReader(f, delimiter="\t")}
t = time.time()
hp.run(pairs, seqs, HMM, Z, 6, str(out / "pairs.domtbl.gz"), str(out / "pairs.parquet"))
print(f"hmmer_pairs on {len(pairs)} models x {len(seqs)} ORFs: {time.time() - t:.0f} s wall", flush=True)


def rows(lines):
    return {(c[0], c[3], c[9]): c for c in (line.split(None, 22) for line in lines if line[0] != "#")}


a = rows(open(out / "hmmsearch.domtbl"))
b = rows(gzip.open(out / "pairs.domtbl.gz", "rt"))
print("domains: hmmsearch", len(a), "pairs", len(b), "only hmmsearch", len(a.keys() - b.keys()), "only pairs", len(b.keys() - a.keys()))
names = "target acc tlen query qacc qlen E score bias n of cE iE dscore dbias hfrom hto afrom ato efrom eto accu desc".split()
diff = Counter()
examples = {}
for k in a.keys() & b.keys():
    for n, x, y in zip(names, a[k], b[k]):
        if x.strip() != y.strip():
            diff[n] += 1
            examples.setdefault(n, (k, x.strip(), y.strip()))
print("columns differing:", dict(diff) or "none")
for n, e in examples.items():
    print("  example", n, e)

for label, path in (("hmmsearch", out / "hmmsearch.domtbl"), ("pairs", out / "pairs.domtbl.gz")):
    cat = "gunzip -c" if str(path).endswith(".gz") else "cat"
    subprocess.run(f"{cat} {path} | python3 {pipe}/modules/local/parsehmmsearchcoverage/bin/hmmer_domtbl_parse_coverage.py -i - -o {out}/{label}.tsv -s {out}/{label}.json", shell=True, check=True)
same = lambda x: (out / f"hmmsearch.{x}").read_text() == (out / f"pairs.{x}").read_text()  # noqa: E731
print("profile identical:", same("tsv"), " stats identical:", same("json"), (out / "pairs.json").read_text())

# Alignments as hmmsearch printed them
ali = {}
query = acc = target = dom = None
buf = []
line_re = re.compile(r"^\s*(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s*$")


def flush():
    if dom is None or not buf:
        return
    parsed = [m.groups() for m in map(line_re.match, buf) if m]
    model = next(p for p in parsed if p[0] == query)
    tgt = next(p for p in parsed if p[0] == target)
    pp = next(line.split()[0] for line in buf if line.rstrip().endswith(" PP"))
    ali[(target, acc, dom)] = (hp.cigar(model[2], tgt[2]), tgt[2].replace("-", "").upper(), pp, int(model[1]), int(model[3]), int(tgt[1]), int(tgt[3]), tgt[2])


with open(out / "hmmsearch.txt") as f:
    for line in f:
        if line.startswith("Query:"):
            flush(); dom = None; buf = []
            query = line.split()[1]
        elif line.startswith("Accession:"):
            acc = line.split()[1].split(".")[0]
        elif line.startswith(">> "):
            flush(); dom = None; buf = []
            target = line.split()[1]
        elif line.lstrip().startswith("== domain"):
            flush(); buf = []
            dom = int(line.split()[2])
        elif dom is not None and line.strip():
            buf.append(line)
        elif dom is not None and buf and not line.strip() and any(x.rstrip().endswith(" PP") for x in buf):
            flush(); dom = None; buf = []
flush()

t = pq.read_table(out / "pairs.parquet").to_pylist()
bad = Counter()
seen = 0
for r in t:
    k = (f"{r['read']}_frame={r['frame']}_begin={r['orf_begin']}_end={r['orf_end']}", r["pfam"], r["domain"])
    h = ali.get(k)
    if h is None:
        bad["not in hmmsearch text"] += 1
        continue
    seen += 1
    for name, x, y in (("cigar", r["cigar"], h[0]), ("residues", r["residues"], h[1]), ("pp", r["pp"], h[2]), ("hmm_from", r["hmm_from"], h[3]), ("hmm_to", r["hmm_to"], h[4]), ("ali_from", r["ali_from"], h[5]), ("ali_to", r["ali_to"], h[6])):
        bad[name] += x != y
    # The round trip the plan promises: cigar and residues back to the aligned target
    rebuilt, i = [], 0
    for n, op in re.findall(r"(\d+)([MID])", r["cigar"]):
        n = int(n)
        if op == "D":
            rebuilt.append("-" * n)
        else:
            rebuilt.append(r["residues"][i:i + n])
            i += n
    bad["round trip"] += "".join(rebuilt) != h[7].upper()
print("parquet rows", len(t), "parsed from hmmsearch text", len(ali), "compared", seen, "mismatches:", {k: v for k, v in bad.items() if v} or "none")
print("parquet bytes", (out / "pairs.parquet").stat().st_size, "per domain", round((out / "pairs.parquet").stat().st_size / max(1, len(t)), 1))
