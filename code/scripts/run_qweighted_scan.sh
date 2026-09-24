#!/usr/bin/env bash
# Launches the 6-point (w2, w3) weight scan for qWeighted.C (single
# combined-q statistic, one shared mask) at radius=15, one job per weight
# pair, in parallel. Each job writes to its own
# {txts,results,root}/qweighted/r15_w2-<w2>_w3-<w3>/ directories and its own log file
# under code/logs/.
#
# Precompiles once first (ACLiC would otherwise race if 6 processes tried
# to compile the same qWeighted_C.so at once), then launches all 6
# execution-only runs together.
#
# Output dirs are pre-created here (rather than left to gSystem->mkdir
# inside the macro) because 6 processes calling a recursive mkdir at the
# same instant previously lost a race and silently dropped one job's txt
# output (see qWeightedSeparate's r15_w2-1000_w3-100 rerun).
set -e
# ROOT macros run from cxx/ (they use ../configs, ../txts, ../results, ../root)
cd "$(dirname "$0")/../cxx"

mkdir -p ../logs

echo "Precompiling qWeighted.C..."
root -l -b -q -e 'gROOT->ProcessLine(".L qWeighted.C+")' > ../logs/precompile_qweighted.log 2>&1

WEIGHTS=("1 1" "10 10" "10 100" "100 10" "100 1000" "1000 100" "0.1 1" "0.1 0.1" "0.01 1" "0.01 0.1")

for wpair in "${WEIGHTS[@]}"; do
  set -- $wpair
  w2=$1; w3=$2
  tag="qweighted/r15_w2-${w2}_w3-${w3}"
  mkdir -p "../txts/${tag}" "../results/${tag}" "../root/${tag}"
  logfile="../logs/qweighted_w2_${w2}_w3_${w3}.log"
  echo "Launching w2=$w2 w3=$w3 -> $logfile"
  ( root -l -b -q "qWeighted.C+(15, ${w2}, ${w3})" > "$logfile" 2>&1; echo "EXIT:$?" >> "$logfile" ) &
done

echo "All 6 jobs launched. Waiting for completion..."
wait
echo "All jobs finished. Check logs/, txts/qweighted/ and results/qweighted/ for output."
