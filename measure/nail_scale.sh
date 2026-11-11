#!/bin/bash
# nail on fragments (recall) and at growing read counts (per-read cost). Inside kfpfam-tools.
M15=/opt/mmseqs-15-6f452/bin/mmseqs
B=/data/baseline
cd $B
run () { # name fasta
  [ -s nail_$1.tbl ] && return
  /usr/bin/time -v nail search -t ${THREADS:-6} -X --mmseqs-path $M15 --tmp-dir $B/tmp_$1 --tbl-out nail_$1.tbl --seeds-out nail_$1.seeds.tsv --use-accession /data/pfam/Pfam-A.hmm "$2" > nail_$1.log 2>&1
  rm -rf $B/tmp_$1
}
run frag /data/truth/frag.faa
run s100k $B/scale100k.faa
run s1m $B/scale1m.faa
echo finished
