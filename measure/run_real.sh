#!/bin/bash
# hmmsearch truth on a random subset of a real metagenome, then nail on the subset and on all reads.
# Nail waits for other scans so its timing is clean. Usage: measure/run_real.sh NAME SUBSET R1 [R2 ...]
set -eu
NAME=$1; SUBSET=$2; shift 2
D=$PWD/data
R=$D/real/$NAME
mkdir -p "$R"
tools () { docker run --rm --platform linux/arm64 --cpus 6 --memory 14g -v "$D":/data kfpfam-tools "$@"; }
[ -s "$R/orfs.all.faa" ] || ${PY:-python3} measure/real_orfs.py "$R/orfs" "$SUBSET" "$@"
Z=20000000
if [ ! -s "$R/subset.domtbl" ]; then
  tools sh -c "/usr/bin/time -v hmmsearch -E 10 --domE 10 -Z $Z --domZ $Z --cpu 6 --noali -o /dev/null --domtblout /data/real/$NAME/subset.domtbl /data/pfam/Pfam-A.hmm /data/real/$NAME/orfs.subset.faa" 2> "$R/hmmsearch.time"
fi
while pgrep -f 'scan bg' > /dev/null; do sleep 20; done
nail () { # label, fasta
  tools sh -c "/usr/bin/time -v nail search -t 6 -X --mmseqs-path /opt/mmseqs-15-6f452/bin/mmseqs --tmp-dir /data/real/$NAME/tmp_$1 --tbl-out /data/real/$NAME/nail_$1.tbl --seeds-out /data/real/$NAME/nail_$1.seeds.tsv --use-accession /data/pfam/Pfam-A.hmm /data/real/$NAME/$2 > /data/real/$NAME/nail_$1.log 2>&1; rm -rf /data/real/$NAME/tmp_$1"
}
[ -s "$R/nail_subset.tbl" ] || nail subset orfs.subset.faa
[ -s "$R/nail_all.tbl" ] || nail all orfs.all.faa
echo real-finished
