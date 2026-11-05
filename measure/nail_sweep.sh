#!/bin/bash
# nail with its filters relaxed, on reads, decoys (names D...) and fragments together. Inside kfpfam-tools.
M15=/opt/mmseqs-15-6f452/bin/mmseqs
B=/data/baseline/sweep
cd $B
run () { # tag, nail flags
  tag=$1; shift
  [ -s $tag.tbl ] && return
  /usr/bin/time -v nail search -t ${THREADS:-6} -X --mmseqs-path $M15 --tmp-dir $B/tmp_$tag --tbl-out $tag.tbl --seeds-out $tag.seeds.tsv --use-accession "$@" /data/pfam/Pfam-A.hmm $B/all.faa > $tag.log 2>&1
  rm -rf $B/tmp_$tag
}
run s2 -S 1e-2 -E 1000
run s0 -S 1 -E 1000
run all1 -S 1 -C 1 -F 1 -E 1000
run static -S 1 -C 1 -F 1 -E 1000 --seed-mode static
run k5 -S 1 -C 1 -F 1 -E 1000 --mmseqs-k 5
echo finished
