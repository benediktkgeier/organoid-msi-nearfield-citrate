#!/usr/bin/env bash
# ============================================================================
# regen_reports.sh — regenerate report/figure PDFs from the EXISTING cache/*.rds.
#   Non-interactive only: this runs the terminal PDF/figure generators that read
#   already-finalized cache. It does NOT recompute preprocessing, SSC clustering,
#   segmentation, or any manual-annotation gate.
#
#   Each script's stdout+stderr is captured to results/_regen_logs/<name>.log.
#   A failing script is logged and SKIPPED (its prior PDF is left in place); the
#   run continues. A pass/fail summary is written to results/_regen_logs/summary.txt.
#
#   Environment: R 4.4.2 ONLY.
#   Usage:  ./regen_reports.sh
# ============================================================================
set -uo pipefail
cd "$(dirname "$0")"
RS="/c/Program Files/R/R-4.4.2/bin/Rscript.exe"
LOGDIR="results/_regen_logs"
mkdir -p "$LOGDIR"
SUMMARY="$LOGDIR/summary.txt"
: > "$SUMMARY"

# Report/figure generators that read finalized cache. Order is informative only.
SCRIPTS=(
  "R/01_preprocess/qc_render.R"
  "R/04_ssc_ontissue/05_ssc_report.R"
  "R/07_metabolite_id/02_metabolite_report.R"
  "R/07_metabolite_id/03_citrate_resolution.R"
  "R/07_metabolite_id/04_citrate_window_images.R"
  "R/07_metabolite_id/05_citrate_isotopes_adducts.R"
  "R/05_registration_refine/04_overlay_report.R"
  "R/06_if_registration/12_dataset_pairs_hq.R"
  "R/08_organoid_gradient_survey/04_report_pdf.R"
  "R/11_per_organoid_final/03_apical_report.R"
  "R/11_per_organoid_final/04_nearfield_figure.R"
  "R/11_per_organoid_final/07_gradient_test_perorganoid.R"
  "R/11_per_organoid_final/08_nearfield_wholeregion.R"
  "R/11_per_organoid_final/09_citrate_gradient_perdataset_v3.R"
  "R/11_per_organoid_final/10_citrate_gradient_report_final.R"
  "R/11_per_organoid_final/13_citrate_gradient_report_3class.R"
)

START=$(date +%s 2>/dev/null || echo 0)
pass=0; fail=0
for s in "${SCRIPTS[@]}"; do
  name=$(echo "$s" | sed 's#/#__#g; s#\.R$##')
  log="$LOGDIR/$name.log"
  echo "=== running $s ==="
  if "$RS" "$s" > "$log" 2>&1; then
    echo "PASS  $s" | tee -a "$SUMMARY"; pass=$((pass+1))
  else
    echo "FAIL  $s   (see $log)" | tee -a "$SUMMARY"; fail=$((fail+1))
  fi
done

# WITHOUT-mixed 2-class variant of the apical-class report (13 with the nomixed flag)
echo "=== running R/11_per_organoid_final/13_citrate_gradient_report_3class.R all nomixed ==="
nm_log="$LOGDIR/R__11_per_organoid_final__13_citrate_gradient_report_3class_nomixed.log"
if "$RS" R/11_per_organoid_final/13_citrate_gradient_report_3class.R all nomixed > "$nm_log" 2>&1; then
  echo "PASS  13_citrate_gradient_report_3class.R all nomixed" | tee -a "$SUMMARY"; pass=$((pass+1))
else
  echo "FAIL  13_citrate_gradient_report_3class.R all nomixed   (see $nm_log)" | tee -a "$SUMMARY"; fail=$((fail+1))
fi

echo "" | tee -a "$SUMMARY"
echo "TOTAL: $pass passed, $fail failed" | tee -a "$SUMMARY"
echo "Logs in $LOGDIR/" | tee -a "$SUMMARY"
