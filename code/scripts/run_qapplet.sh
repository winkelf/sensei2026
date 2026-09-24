#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

(cd "$SCRIPT_DIR/../cxx"    && root -l -b -q qAppletExport.C)
(cd "$SCRIPT_DIR/../python" && python3 qApplet.py "$@")
