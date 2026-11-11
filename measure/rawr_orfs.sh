#!/bin/bash
# Amino-acid sequences as raw-reads-analysis-pipeline makes them for Pfam profiling, without its
# decontamination: standardise, repair, fastp QC, fastp merge (unmerged kept), six-frame translation.
# Same containers and arguments as v6.0.0's conf/modules.config. Usage: rawr_orfs.sh OUT R1.fastq.gz R2.fastq.gz
set -eu
OUT=$1; R1=$(realpath "$2"); R2=$(realpath "$3"); mkdir -p "$OUT"
IN=$(dirname "$R1"); A=/in/$(basename "$R1"); B=/in/$(basename "$R2")
run () { img=$1; shift; docker run --rm --platform linux/amd64 --cpus 6 --memory 10g -v "$IN":/in:ro -v "$(realpath "$OUT")":/w -w /w "$img" "$@"; }
BB=community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9
FP=community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690
SK=quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0
TR="--threads 4 --append-frame --allow-unknown-codon --frame 6 --out-subseqs --min-len 15"
run $BB reformat.sh -Xmx8g in=$A in2=$B out=std_1.fastq.gz out2=std_2.fastq.gz outs=std_s.fastq.gz threads=1 qin=33 qout=33 allowidenticalnames=t trimreaddescription=t addslash=t spaceslash=f verifypaired=f 2> "$OUT/std.log"
run $BB repair.sh -Xmx8g in=std_1.fastq.gz in2=std_2.fastq.gz out=rep_1.fastq.gz out2=rep_2.fastq.gz outs=rep_s.fastq.gz threads=1 2> "$OUT/rep.log"
run $FP fastp --in1 rep_1.fastq.gz --in2 rep_2.fastq.gz --out1 qc_1.fastq.gz --out2 qc_2.fastq.gz --json qc.json --html qc.html --thread 4 --detect_adapter_for_pe 2> "$OUT/qc.log"
run $FP fastp --in1 qc_1.fastq.gz --in2 qc_2.fastq.gz --out1 un_1.fastq.gz --out2 un_2.fastq.gz --json merge.json --html merge.html -m --merged_out merged.fastq.gz --thread 2 --detect_adapter_for_pe --include_unmerged -A -G -Q -L 2> "$OUT/merge.log"
run $BB reformat.sh -Xmx8g in=merged.fastq.gz out=merged.fasta threads=1 qin=33 qout=33 uniquenames=t trimreaddescription=t 2> "$OUT/fasta.log"
run $SK seqkit translate $TR merged.fasta > "$OUT/orfs.faa"
# Control: the same translator on the untouched reads, to separate it from QC and merging
run $SK sh -c "zcat $A $B | seqkit translate $TR" > "$OUT/raw_orfs.faa"
rm -f "$OUT"/std_* "$OUT"/rep_* "$OUT"/un_*
echo rawr-finished
