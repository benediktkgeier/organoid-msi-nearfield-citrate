#!/usr/bin/env bash
# ============================================================================
# run_wholeregion_figures.sh - regenerate the LOCKED whole-region figures built
#   on the consensus-curated segmentation:
#     1) per-dataset citrate/DHA gradient report, ALL 20 datasets (v3)
#     2) whole-region near-field citrate emission figure, ALL 20 datasets
#        - default heatmap (per-section auto-scale)
#        - globalheat heatmap (one scale pooled across all 20 datasets)
#
#   Reads ONLY finalized cache/*.rds (no preprocessing / clustering / annotation
#   gate is recomputed). See docs/citrate_gradient_perdataset.md and
#   docs/nearfield_wholeregion.md.
#
#   Environment: R 4.4.2 ONLY.
#   Usage:  ./run_wholeregion_figures.sh
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
RS="/c/Program Files/R/R-4.4.2/bin/Rscript.exe"

echo "== [1/3] per-dataset citrate/DHA gradient report - ALL 20 datasets (v3) =="
"$RS" R/11_per_organoid_final/09_citrate_gradient_perdataset_v3.R all

echo "== [2/3] whole-region near-field figure - ALL 20 (per-section heatmap) =="
"$RS" R/11_per_organoid_final/08_nearfield_wholeregion.R all

echo "== [3/3] whole-region near-field figure - ALL 20 (GLOBAL heatmap scale) =="
"$RS" R/11_per_organoid_final/08_nearfield_wholeregion.R all globalheat

echo "== DONE =="
echo "  perdataset v3 : figures/gradient/citrate_gradient_perdataset_v3.pdf            (21 pages)"
echo "  nearfield all : figures/annotation/apical_nearfield_emission_figure_all.pdf     (20 pages)"
echo "  nearfield gh  : figures/annotation/apical_nearfield_emission_figure_all_globalheat.pdf (20 pages)"
echo ""
echo "  Single section (optional):"
echo "    Rscript R/11_per_organoid_final/08_nearfield_wholeregion.R AO_0h_sl6A_sec2a [globalheat]"
echo "    Rscript R/11_per_organoid_final/09_citrate_gradient_perdataset_v3.R AO_0h_sl6A_sec2a"
