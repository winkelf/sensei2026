#!/bin/bash
# Runs the window studies for one dataset folder that has its own config.yaml
# (e.g. results/hits_skp_2025-09-28/). All outputs land inside that folder,
# so datasets never overwrite each other or the original results.
#
# usage: run_dataset.sh <dataset_dir> [stage ...]
#   stages: scan  -> MaskScanManye.C      (window trees, standard mask)
#           lec   -> MaskScanManyeLEC.C   (window trees, LEC mask)
#           rates -> ImageRates.C         (per-image n2 and Npix)
#           plots -> PlotK2EqualsTwo.C, PlotQScenarios.C, PlotK1Independence.C
#   default: all four, in that order. The generation stages skip any output
#   file that already exists, so an interrupted run can just be restarted.

set -e

SRC=$(cd "$(dirname "$0")" && pwd)
DIR=$1
shift || true

if [ -z "$DIR" ] || [ ! -f "$DIR/config.yaml" ]; then
    echo "usage: $0 <dataset_dir with config.yaml> [scan|lec|rates|plots ...]" >&2
    exit 1
fi

cd "$DIR"
STAGES=${*:-scan lec rates plots}

for s in $STAGES; do
    echo "=== $(date '+%F %T')  $DIR: $s ==="
    case $s in
        scan)  root -l -b -q "$SRC/MaskScanManye.C" ;;
        lec)   root -l -b -q "$SRC/MaskScanManyeLEC.C" ;;
        rates) root -l -b -q "$SRC/ImageRates.C" ;;
        plots) root -l -b -q "$SRC/PlotK2EqualsTwo.C"
               root -l -b -q "$SRC/PlotQScenarios.C"
               root -l -b -q "$SRC/PlotK1Independence.C" ;;
        *)     echo "unknown stage: $s" >&2; exit 1 ;;
    esac
done

echo "=== $(date '+%F %T')  $DIR: done ==="
