#!/bin/bash
# One sketch job at a time under a memory cap (the laptop has crashed on local runs).
# Usage: sketch_run.sh ORFS K T SCALED fwd|rev OUT
set -u
CAP_KB=${CAP_KB:-10000000}
export RAYON_NUM_THREADS=${RAYON_NUM_THREADS:-4}
cd "$(dirname "$0")/.."
filter/target/release/filter sketch data/prep "$1" "$2" "$3" "$4" "$5" "$6" 2>>data/sketch/runs.log &
pid=$!
while kill -0 $pid 2>/dev/null; do
  rss=$(ps -o rss= -p $pid | tr -d ' ')
  if [ "${rss:-0}" -gt "$CAP_KB" ]; then
    kill -9 $pid
    echo "KILLED over cap: $*" | tee -a data/sketch/runs.log
    exit 1
  fi
  sleep 2
done
wait $pid
