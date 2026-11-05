#!/bin/bash
# Phase 1: MMseqs2 and nail as pre-filters on the real read ORFs. Run inside kfpfam-tools, data at /data.
set -u
T=${THREADS:-6}
B=/data/baseline
Q=/data/explore/truth/reads10k.faa
M18=/opt/mmseqs-18-8cc5c/bin/mmseqs
M15=/opt/mmseqs-15-6f452/bin/mmseqs
mkdir -p $B/tmp
cd $B
if [ ! -s profDb.dbtype ]; then
  $M18 convertmsa /data/pfam/Pfam-A.seed.gz msaDb --identifier-field 1 -v 1
  $M18 msa2profile msaDb profDb --match-mode 1 --match-ratio 0.5 --threads $T -v 1
fi
[ -s orfDb.dbtype ] || $M18 createdb $Q orfDb -v 1
run_pref () { # name, extra flags
  name=$1; shift
  [ -s pref_$name.tsv ] && return
  /usr/bin/time -v $M18 prefilter orfDb profDb pref_$name --threads $T --max-seqs 2000 "$@" > pref_$name.log 2>&1
  $M18 createtsv orfDb profDb pref_$name pref_$name.tsv -v 1
  rm -f pref_$name pref_$name.[0-9]* pref_$name.index pref_$name.dbtype
}
for s in 1 4 5.7 7.5; do run_pref s$s -s $s; done
run_pref s7.5u0 -s 7.5 --min-ungapped-score 0
run_pref s5.7k5 -s 5.7 -k 5
run_pref s7.5k5 -s 7.5 -k 5
# nail, with its own pipeline on the released HMMs
if [ ! -s nail.tbl ]; then
  /usr/bin/time -v nail search -t $T -X --mmseqs-path $M15 --tmp-dir $B/tmp --tbl-out nail.tbl --seeds-out nail.seeds.tsv --use-accession /data/pfam/Pfam-A.hmm $Q > nail.log 2>&1
fi
echo finished
