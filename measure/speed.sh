#!/bin/bash
# Speed plan step 0: one engine run on the 50,000 ERR10889147 reads and their reversals, timed, scored
# against hmmsearch, and checked against the reference index search (data/index, today's filter search).
# Usage: measure/speed.sh NAME THREADS SUBCOMMAND [ARGS...]
#   runs: filter SUBCOMMAND data/prep ARGS... ORFS OUT.tsv   (ORFS and OUT appended), writing data/speed/NAME/
# Output TSV: query, model, diagonal, gapped, as filter search writes. One job at a time, as the laptop needs.
set -euo pipefail
cd "$(dirname "$0")/.."
NAME=$1 THREADS=$2 SUB=$3
shift 3
R=data/real/ERR10889147
OUT=data/speed/$NAME
mkdir -p "$OUT"
# ORFS and DECOY override the ORF sets (e.g. a frame rule's output); the truth stays hmmsearch on all ORFs
for side in real decoy; do
  orfs=${ORFS:-$R/orfs.subset.faa}
  [ $side = decoy ] && orfs=${DECOY:-$R/decoy.faa}
  RAYON_NUM_THREADS=$THREADS /usr/bin/time -l filter/target/release/filter "$SUB" data/prep "$@" "$orfs" "$OUT/$side.tsv" \
    > "$OUT/$side.log" 2>&1
  grep -E "ORFs|wall|real|maximum resident" "$OUT/$side.log" | head -4
done
python3 measure/index_score.py "$OUT/real.tsv" "$OUT/decoy.tsv" $R/subset.e1000.domtbl $R/decoy.e1000.domtbl \
  data/prep/models.tsv 50000 > "$OUT/score.txt"
cat "$OUT/score.txt"
python3 - "$OUT/real.tsv" data/index/real.tsv <<'EOF'
import sys
def pairs(p):
    out = {}
    for line in open(p):
        q, m, d, g = line.split()
        if float(d) >= 22:
            out[q, m] = max(out.get((q, m), 0), float(d))
    return out
new, ref = pairs(sys.argv[1]), pairs(sys.argv[2])
shared = new.keys() & ref.keys()
worst = max((abs(new[k] - ref[k]) for k in shared), default=0)
print(f"against the reference at diagonal >= 22: {len(shared)} shared, {len(new.keys() - ref.keys())} new only, "
      f"{len(ref.keys() - new.keys())} reference only, largest score difference {worst:.2f}")
EOF
