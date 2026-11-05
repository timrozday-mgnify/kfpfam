#!/bin/bash
# Score the policies on a real metagenome: run_real.sh's reads, against hmmsearch.
# Recall and precision of reported hits come from every hit of the subset; candidates hmmsearch does
# not report need a scan over all models, so those use BG reads. Usage: run_real_policy.sh NAME BG
set -eu
NAME=$1; BG=$2
PY=${PY:-python3}
D=$PWD/data
R=$D/real/$NAME
S=$D/scan3
Z=20000000
OLD=$D/scan2
tools () { docker run --rm --platform linux/arm64 --cpus 6 --memory 12g -v "$D":/data kfpfam-tools "$@"; }
SCAN=measure/scan/target/release/scan
export RAYON_NUM_THREADS=6
# Reversed, not reverse-complemented: same composition, no homology
[ -s "$R/decoy.faa" ] || $PY -c "
import sys
for l in open('$R/orfs.subset.faa'):
    sys.stdout.write(l if l[0]=='>' else l.strip()[::-1]+'\n')" > "$R/decoy.faa"
[ -s "$R/decoy.e1000.domtbl" ] || tools hmmsearch -E 1000 --domE 1000 -Z $Z --domZ $Z --cpu 6 --noali -o /dev/null --domtblout "/data/real/$NAME/decoy.e1000.domtbl" /data/pfam/Pfam-A.hmm "/data/real/$NAME/decoy.faa"
[ -s "$R/nail_decoy.tbl" ] || tools sh -c "nail search -t 6 -X --mmseqs-path /opt/mmseqs-15-6f452/bin/mmseqs --tmp-dir /data/real/$NAME/tmp_decoy --tbl-out /data/real/$NAME/nail_decoy.tbl --seeds-out /data/real/$NAME/nail_decoy.seeds.tsv --use-accession /data/pfam/Pfam-A.hmm /data/real/$NAME/decoy.faa > /data/real/$NAME/nail_decoy.log 2>&1; rm -rf /data/real/$NAME/tmp_decoy"
mkdir -p "$S/truth" "$S/real/bg" "$S/real/decoy" "$S/real/hmm"
for x in tables bg decoy hmm; do ln -sfn "../scan2/$x" "$S/$x"; done
for x in abundance.tsv bg.orfs.tsv bg.best.tsv decoy.orfs.tsv decoy.best.tsv; do ln -sf "../scan2/$x" "$S/$x"; done
# Pairs: everything scored before, then the real reads and their decoys
if [ ! -s "$S/pairs.best.tsv" ]; then
  cp "$OLD/pairs.seqs.tsv" "$OLD/pairs.meta.tsv" "$S/"
  OLDN=$(wc -l < "$OLD/pairs.seqs.tsv")
  $PY measure/sets.py extend real "$R/orfs.subset.faa" "$R/subset.e1000.domtbl" "$D/prep/models.tsv" "$S/pairs" 3000
  $PY measure/sets.py extend realdecoy "$R/decoy.faa" "$R/decoy.e1000.domtbl" "$D/prep/models.tsv" "$S/pairs" 3000
  $SCAN truth "$D/prep" "$S/pairs.seqs.tsv" "$S/truth"
  awk -F'\t' -v n="$OLDN" 'NR>n{print NR-1"\t"$1}' "$S/pairs.seqs.tsv" > "$S/new.ids.tsv"
  $SCAN best "$D/prep" "$S/new.ids.tsv" "$S/new.best.tsv"
  cat "$OLD/pairs.best.tsv" "$S/new.best.tsv" | sort -n -k1,1 > "$S/pairs.best.tsv"
fi
# Background and decoys: the first BG reads, a random sample as real_orfs.py shuffles
Q=$S/real
grep -c '>' "$R/orfs.subset.faa" > /dev/null
$PY - "$R/orfs.subset.faa" "$BG" "$Q/bg.orfs.tsv" "$Q/decoy.orfs.tsv" "$Q/truth_reads" <<'PY'
import sys
path, bg, a, b, n = sys.argv[1:6]
bg, top = int(bg), 0
with open(path) as f, open(a, "w") as x, open(b, "w") as y:
    for line in f:
        if line[0] == ">":
            read = int(line[2:].split("_")[0])
            top = max(top, read)
        elif read < bg:
            s = line.strip()
            x.write(f"{read}\t{s}\n")
            y.write(f"{read}\t{s[::-1]}\n")
open(n, "w").write(str(top + 1))
PY
ln -sf "$R/subset.e1000.domtbl" "$Q/hmm/reads.domtbl"
ln -sf "$R/decoy.e1000.domtbl" "$Q/hmm/decoy.domtbl"
[ -s "$Q/bg/p33.bg.gz" ] || $SCAN bg "$D/prep" "$Q/bg.orfs.tsv" 16 "$Q/bg"
[ -s "$Q/decoy/p33.bg.gz" ] || $SCAN bg "$D/prep" "$Q/decoy.orfs.tsv" 16 "$Q/decoy"
[ -s "$Q/bg.best.tsv" ] || $SCAN best "$D/prep" "$Q/bg.orfs.tsv" "$Q/bg.best.tsv"
[ -s "$Q/decoy.best.tsv" ] || $SCAN best "$D/prep" "$Q/decoy.orfs.tsv" "$Q/decoy.best.tsv"
$PY measure/export.py "$S" "$D/prep" "$S/abundance.tsv" web/data
echo real-policy-done
