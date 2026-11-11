#!/bin/bash
# The sketch sweep on ERR10889147, one job at a time: 50,000 reads and all 1.07 M reads.
set -u
cd "$(dirname "$0")/.."
R=data/real/ERR10889147
pair() { # tag k T scaled orfs
  local o=data/sketch/$1.k$2T$3.s$4
  [ -s $o.fwd.tsv ] || measure/sketch_run.sh $5 $2 $3 $4 fwd $o.fwd.tsv
  [ -s $o.rev.tsv ] || measure/sketch_run.sh $5 $2 $3 $4 rev $o.rev.tsv
  echo "done $o"
}
for c in "5 10" "6 12" "7 14" "8 18" "12 24" "14 26"; do read -r k t <<<"$c"; pair sub $k $t 1 $R/orfs.subset.faa; done
for s in 8 64; do pair sub 10 20 $s $R/orfs.subset.faa; done
for c in "7 14" "10 20" "12 24"; do read -r k t <<<"$c"; pair all $k $t 1 $R/orfs.all.faa; done
