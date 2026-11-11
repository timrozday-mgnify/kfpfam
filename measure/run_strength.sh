#!/bin/bash
# Truth down to E <= 1000, reversed-read decoys, and the scans to break recall and precision down by
# match strength. Run from the repo root once the tools image and data/scan exist.
set -eu
PY=${PY:-python3}
D=$PWD/data
T=$D/truth
S=$D/scan2
Z=20000000
tools () { docker run --rm --platform linux/arm64 --cpus 6 --memory 14g -v "$D":/data kfpfam-tools "$@"; }
search () { # query fasta under /data, output
  [ -s "$D/$2" ] || tools hmmsearch -E 1000 --domE 1000 -Z $Z --domZ $Z --cpu 6 --noali -o /dev/null --domtblout "/data/$2" /data/pfam/Pfam-A.hmm "/data/$1"
}
mkdir -p "$S/truth" "$S/decoy" "$S/hmm"
ln -sfn ../scan/tables "$S/tables"
ln -sfn ../scan/bg "$S/bg"
# Reversed, not reverse-complemented: same composition, no homology
[ -s "$T/decoy.faa" ] || $PY -c "
import sys
for l in open('$D/explore/truth/reads10k.faa'):
    sys.stdout.write(l if l[0]=='>' else l.strip()[::-1]+'\n')" > "$T/decoy.faa"
search truth/frag.faa truth/frag.e1000.domtbl
search explore/truth/reads10k.faa truth/reads.e1000.domtbl
search truth/decoy.faa truth/decoy.e1000.domtbl
search real/ERR10889147/orfs.subset.faa real/ERR10889147/subset.e1000.domtbl
[ -s "$D/baseline/nail_decoy.tbl" ] || tools sh -c "nail search -t 6 -X --mmseqs-path /opt/mmseqs-15-6f452/bin/mmseqs --tmp-dir /data/baseline/tmp_decoy --tbl-out /data/baseline/nail_decoy.tbl --seeds-out /data/baseline/nail_decoy.seeds.tsv --use-accession /data/pfam/Pfam-A.hmm /data/truth/decoy.faa > /data/baseline/nail_decoy.log 2>&1; rm -rf /data/baseline/tmp_decoy"
ln -sf ../../truth/reads.e1000.domtbl "$S/hmm/reads.domtbl"
ln -sf ../../truth/decoy.e1000.domtbl "$S/hmm/decoy.domtbl"
cp -n "$D/scan/abundance.tsv" "$D/scan/bg.orfs.tsv" "$D/scan/bg.best.tsv" "$S/" 2>/dev/null || true
[ -s "$S/pairs.meta.tsv" ] || $PY measure/sets.py pairs "$T/frag.tsv" "$T/frag.e1000.domtbl" "$T/frag.faa" "$D/explore/truth/reads10k.faa" "$T/reads.e1000.domtbl" "$D/prep/models.tsv" "$S/pairs"
SCAN=measure/scan/target/release/scan
export RAYON_NUM_THREADS=6
[ -s "$S/truth/p33.rows.gz" ] || $SCAN truth "$D/prep" "$S/pairs.seqs.tsv" "$S/truth"
# Decoys of the same 300 reads
if [ ! -s "$S/decoy.orfs.tsv" ]; then
  awk -F'\t' 'BEGIN{OFS="\t"}{n=split($2,a,""); s=""; for(i=n;i>0;i--) s=s a[i]; print $1,s}' "$S/bg.orfs.tsv" > "$S/decoy.orfs.tsv"
fi
[ -s "$S/decoy/p33.bg.gz" ] || $SCAN bg "$D/prep" "$S/decoy.orfs.tsv" 16 "$S/decoy"
[ -s "$S/decoy.best.tsv" ] || $SCAN best "$D/prep" "$S/decoy.orfs.tsv" "$S/decoy.best.tsv"
# Best diagonals: only sequences not scored before
if [ ! -s "$S/pairs.best.tsv" ]; then
  $PY - "$D/scan/pairs.seqs.tsv" "$D/scan/pairs.best.tsv" "$S/pairs.seqs.tsv" "$S/new.ids.tsv" "$S/old.best.tsv" <<'PY'
import sys
old_seq, old_best, new_seq, new_ids, carried = sys.argv[1:6]
seqs = [l.split("\t")[0] for l in open(old_seq)]
best = {s: l.rstrip("\n").split("\t")[1:] for s, l in zip(seqs, open(old_best))}
with open(new_ids, "w") as n, open(carried, "w") as c:
    for i, l in enumerate(open(new_seq)):
        s = l.split("\t")[0]
        if s in best:
            c.write(f"{i}\t{best[s][0]}\t{best[s][1]}\n")
        else:
            n.write(f"{i}\t{s}\n")
PY
  $SCAN best "$D/prep" "$S/new.ids.tsv" "$S/new.best.tsv"
  cat "$S/old.best.tsv" "$S/new.best.tsv" | sort -n -k1,1 > "$S/pairs.best.tsv"
fi
$PY measure/baseline_bins.py "$T" "$D/baseline" "$S/pairs.meta.tsv" "$D/prep/models.tsv" 10000 "$S/byBin.json"
$PY - "$S/byBin.json" web/data/baselines.json <<'PY'
import json, sys
d = json.load(open(sys.argv[2]))
d["byBin"] = json.load(open(sys.argv[1]))
d["provenance"] = "Truth: hmmsearch on fragments of MGnify proteins and reads simulated from reference genomes, at Z = 20 million, down to E 1000; decoys are the reads reversed."
json.dump(d, open(sys.argv[2], "w"))
PY
$PY measure/export.py "$S" "$D/prep" "$S/abundance.tsv" web/data
echo strength-done
