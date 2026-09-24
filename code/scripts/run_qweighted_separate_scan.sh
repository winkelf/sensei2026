#!/usr/bin/env bash
# Launches the 6-point (w2, w3) weight scan for qWeightedSeparate.C at
# radius=15, one job per weight pair, in parallel. Each job writes to its
# own {txts,results,root}/qweighted_separate/r15_w2-<w2>_w3-<w3>/ directories and its own
# log file under code/logs/.
#
# Precompiles once first (ACLiC would otherwise race if 6 processes tried
# to compile the same qWeightedSeparate_C.so at once), then launches all 6
# execution-only runs together.
set -e
# ROOT macros run from cxx/ (they use ../configs, ../txts, ../results, ../root)
cd "$(dirname "$0")/../cxx"

mkdir -p ../logs

echo "Precompiling qWeightedSeparate.C..."
root -l -b -q -e 'gROOT->ProcessLine(".L qWeightedSeparate.C+")' > ../logs/precompile.log 2>&1

WEIGHTS=("1 1" "10 10" "10 100" "100 10" "100 1000" "1000 100")

for wpair in "${WEIGHTS[@]}"; do
  set -- $wpair
  w2=$1; w3=$2
  tag="qweighted_separate/r15_w2-${w2}_w3-${w3}"
  mkdir -p "../txts/${tag}" "../results/${tag}" "../root/${tag}"
  logfile="../logs/w2_${w2}_w3_${w3}.log"
  echo "Launching w2=$w2 w3=$w3 -> $logfile"
  ( root -l -b -q "qWeightedSeparate.C+(15, ${w2}.0, ${w3}.0)" > "$logfile" 2>&1; echo "EXIT:$?" >> "$logfile" ) &
done

echo "All 6 jobs launched. Waiting for completion..."
wait
echo "All jobs finished. Check logs/, txts/qweighted_separate/ and results/qweighted_separate/ for output."
