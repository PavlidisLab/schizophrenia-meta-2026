#!/usr/bin/env bash
# Rebuild every figure from the figure data in this repository.
#   ./reproduce_figures.sh          check the data against the manifests, then rebuild
# Figures are written to figures/results/ (SVG, PNG and a *_meta.json with input and output sha256).
set -euo pipefail
cd "$(dirname "$0")"
[ $# -eq 0 ] || { echo "unknown option $1" >&2; exit 2; }
python3 figures/scripts/verify_data.py
python3 figures/scripts/rebuild_all.py
