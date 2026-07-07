#!/usr/bin/env bash
# ============================================================================
# run_apical_nearfield_pipeline.sh - reproduce the apical orientation ->
#   near-field citrate emission analysis + LOCKED publication figure.
#   See docs/apical_nearfield.md. Reuses the two-annotator CONSENSUS apical map
#   (results/annotation/apical_map_consensus.csv) and the finalized cache
#   (curated instances/zones); no PDF re-parse / python needed.
# Usage:  ./run_apical_nearfield_pipeline.sh
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
RS="/c/Program Files/R/R-4.4.2/bin/Rscript.exe"

echo "== [1/2] consensus apical report + apical_gradient_per_organoid.csv =="
"$RS" R/11_per_organoid_final/03_apical_report.R

echo "== [2/2] LOCKED four-view near-field publication figure =="
"$RS" R/11_per_organoid_final/04_nearfield_figure.R

echo "== DONE =="
echo "  report  : figures/annotation/apical_citrate_dha_report.pdf"
echo "  figure  : figures/annotation/apical_nearfield_emission_figure.pdf"
echo "  panels  : figures/annotation/nearfield_panels/"
