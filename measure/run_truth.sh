#!/bin/bash
# Build the evaluation sets and the dashboard's measurements. One heavy job at a time.
# Run from the repo root: measure/run_truth.sh  (needs data/explore, data/prep, the kfpfam-tools image)
set -eu
PY=${PY:-python3}
THREADS=${THREADS:-6}
D=$PWD/data
T=$D/truth
S=$D/scan
mkdir -p "$T" "$S"/{tables,truth,bg}
tools () { docker run --rm --platform linux/arm64 --cpus "$THREADS" --memory 14g -v "$D":/data kfpfam-tools "$@"; }
FAA=../kmer-functional-profiler/data/mgnify/mgy_clusters_head.faa
[ -s "$T/proteins.faa" ] || cp "${MGY_FASTA:-$FAA}" "$T/proteins.faa"
# 1. Pfam on every protein not yet annotated
if [ ! -s "$T/rest.domtbl" ]; then
  $PY - "$T/proteins.faa" "$D/explore/truth/proteins40k.faa" "$T/rest.faa" <<'PY'
import sys
done = {l[1:].split()[0] for l in open(sys.argv[2]) if l[0] == ">"}
keep = False
with open(sys.argv[3], "w") as o:
    for l in open(sys.argv[1]):
        if l[0] == ">":
            keep = l[1:].split()[0] not in done
        if keep:
            o.write(l)
PY
  tools hmmsearch --cut_ga --cpu "$THREADS" --noali -o /dev/null --domtblout /data/truth/rest.domtbl /data/pfam/Pfam-A.hmm /data/truth/rest.faa
fi
TABLES=$D/explore/truth/proteins40k.domtbl,$T/rest.domtbl
# 2. Family-balanced fragments, and hmmsearch on them and on the reads, at a realistic Z
[ -s "$T/frag.faa" ] || $PY measure/sets.py fragments "$TABLES" "$T/proteins.faa" "$T/frag" 4
Z=20000000
[ -s "$T/frag.domtbl" ] || tools hmmsearch -E 10 --domE 10 -Z $Z --domZ $Z --cpu "$THREADS" --noali -o /dev/null --domtblout /data/truth/frag.domtbl /data/pfam/Pfam-A.hmm /data/truth/frag.faa
[ -s "$T/reads.domtbl" ] || tools hmmsearch -E 10 --domE 10 -Z $Z --domZ $Z --cpu "$THREADS" --noali -o /dev/null --domtblout /data/truth/reads.domtbl /data/pfam/Pfam-A.hmm /data/explore/truth/reads10k.faa
# 3. Pairs and background
[ -s "$S/pairs.meta.tsv" ] || $PY measure/sets.py pairs "$T/frag.tsv" "$T/frag.domtbl" "$T/frag.faa" "$D/explore/truth/reads10k.faa" "$T/reads.domtbl" "$D/prep/models.tsv" "$S/pairs"
[ -s "$S/bg.orfs.tsv" ] || $PY - "$D/explore/truth/reads10k.json" "$S/bg.orfs.tsv" 300 <<'PY'
import json, sys
reads = json.load(open(sys.argv[1]))[:int(sys.argv[3])]
with open(sys.argv[2], "w") as o:
    for i, r in enumerate(reads):
        for orf in r["orfs"]:
            o.write(f"{i}\t{orf}\n")
PY
$PY - "$D/explore/truth/proteins40k.domtbl" "$T/rest.domtbl" "$S/abundance.tsv" <<'PY'
import sys
from collections import Counter
c = Counter()
for p in sys.argv[1:3]:
    for l in open(p):
        if l[0] != "#":
            c[l.split()[4].split(".")[0]] += 1
with open(sys.argv[3], "w") as o:
    for k, v in c.items():
        o.write(f"{k}\t{v}\n")
print(len(c), "families in the proteins")
PY
# 4. The scans
SCAN=measure/scan/target/release/scan
export RAYON_NUM_THREADS=$THREADS
[ -s "$S/tables/p33.strata.f64" ] || $SCAN tables "$D/prep" "$S/tables"
[ -s "$S/truth/p33.rows.gz" ] || $SCAN truth "$D/prep" "$S/pairs.seqs.tsv" "$S/truth"
[ -s "$S/bg/p33.bg.gz" ] || $SCAN bg "$D/prep" "$S/bg.orfs.tsv" 16 "$S/bg"
[ -s "$S/bg.best.tsv" ] || $SCAN best "$D/prep" "$S/bg.orfs.tsv" "$S/bg.best.tsv"
if [ ! -s "$S/pairs.best.tsv" ]; then
  awk -F'\t' '{print NR-1"\t"$1}' "$S/pairs.seqs.tsv" > "$S/pairs.ids.tsv"
  $SCAN best "$D/prep" "$S/pairs.ids.tsv" "$S/pairs.best.tsv"
fi
echo all-measured
