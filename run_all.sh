#!/usr/bin/env bash
# ============================================================================
# run_all.sh — regenerate the publication figure/result set for the fork.
#   Organoid MSI near-field citrate project (apical-out vs basolateral-out;
#   single condition: all sections incubated in CMC >=5 min, treated identically).
#
#   THIS FORK REUSES the validated, timepoint-agnostic intermediate caches from
#   the working analysis (Analysis_R/cache) READ-ONLY, via CACHE_SRC + cache_in()
#   in R/00_lib/lib_paths.R. It does NOT rebuild preprocessing / SSC / segmentation
#   and does NOT re-run the manual annotation gates (their marked artifacts are
#   reused). It regenerates the analysis/figure layer only, in dependency order.
#
#   Each script's stdout+stderr -> results/_regen_logs/<name>.log; a failing
#   script is logged and SKIPPED (run continues). Summary at the end.
#
#   The upstream QC figures (citrate-standard gate = Phase 02, metabolite ID =
#   Phase 07) read heavier raw data and are OFF by default; enable with
#   `./run_all.sh qc`.
#
#   SSC tuning: this regenerates the SSC report (04/05) from the cached mask. To
#   change on-tissue clustering, first re-run `R/04_ssc_ontissue/04_ssc_tissue_mask.R`
#   (rebuilds cache/peaks_tissue_combined.rds into the fork-local cache/), THEN
#   `./run_all.sh` to refresh the report and every downstream figure.
#
#   Environment: R 4.4.2 ONLY (Cardinal v3). Never setCardinalBPPARAM().
#   Usage:  ./run_all.sh        # regenerate the publication figure set
#           ./run_all.sh qc     # also regenerate upstream QC (citrate/metabolite)
# ============================================================================
set -uo pipefail
cd "$(dirname "$0")"
RS="/c/Program Files/R/R-4.4.2/bin/Rscript.exe"
LOGDIR="results/_regen_logs"; mkdir -p "$LOGDIR"
SUMMARY="$LOGDIR/summary.txt"; : > "$SUMMARY"

# Dependency-ordered figure/result generators (each reads reused cache).
SCRIPTS=(
  # SSC on-tissue clustering + organoid segmentation report (re-run after tuning SSC
  # params in R/04_ssc_ontissue/04_ssc_tissue_mask.R to see the effect)
  "R/04_ssc_ontissue/05_ssc_report.R"
  # pooled gradient survey (supporting context)
  "R/08_organoid_gradient_survey/03_gradient_stats.R"
  "R/08_organoid_gradient_survey/04_report_pdf.R"
  # apical near-field chain (headline) — 03 reads the CONSENSUS map by default
  "R/11_per_organoid_final/03_apical_report.R"
  "R/11_per_organoid_final/04_nearfield_figure.R"
  "R/11_per_organoid_final/07_gradient_test_perorganoid.R"
  "R/11_per_organoid_final/08_nearfield_wholeregion.R"
  "R/11_per_organoid_final/09_citrate_gradient_perdataset_v3.R"
  # final combined per-dataset report (MSI / BF / SSC / IF / overlay / gradient map)
  "R/11_per_organoid_final/10_citrate_gradient_report_final.R"
  # apical-class variant of the final report: default = WITH mixed (mixed = own grey
  # trend line); the WITHOUT-mixed 2-class variant is a 2nd invocation below
  "R/11_per_organoid_final/13_citrate_gradient_report_3class.R"
  # IF pairs report (ZO-1 red / DAPI cyan)
  "R/06_if_registration/12_dataset_pairs_hq.R"
  # methods & pipeline report
  "figures/methods_report/generate_methods_report.R"
)

if [ "${1:-}" = "qc" ]; then
  SCRIPTS=(
    "R/02_CitrateStandard/01_standard_spectra_mz.R"
    "R/02_CitrateStandard/02_calibration_curve.R"
    "R/02_CitrateStandard/03_id_fingerprint.R"
    "R/02_CitrateStandard/04_standard_anchored_citrate.R"
    "R/07_metabolite_id/02_metabolite_report.R"
    "R/07_metabolite_id/03_citrate_resolution.R"
    "R/07_metabolite_id/04_citrate_window_images.R"
    "R/07_metabolite_id/05_citrate_isotopes_adducts.R"
    "${SCRIPTS[@]}"
  )
fi

pass=0; fail=0
for s in "${SCRIPTS[@]}"; do
  name=$(echo "$s" | sed 's#/#__#g; s#\.R$##')
  log="$LOGDIR/$name.log"
  echo; echo ">>> $s"
  if "$RS" "$s" > "$log" 2>&1; then
    echo "PASS  $s" | tee -a "$SUMMARY"; pass=$((pass+1))
  else
    echo "FAIL  $s   (see $log)" | tee -a "$SUMMARY"; fail=$((fail+1))
  fi
done

# whole-region GLOBAL-heatmap variant (2nd invocation of 08 with the globalheat flag)
echo; echo ">>> R/11_per_organoid_final/08_nearfield_wholeregion.R all globalheat"
gh_log="$LOGDIR/R__11_per_organoid_final__08_nearfield_wholeregion_globalheat.log"
if "$RS" R/11_per_organoid_final/08_nearfield_wholeregion.R all globalheat > "$gh_log" 2>&1; then
  echo "PASS  08_nearfield_wholeregion.R all globalheat" | tee -a "$SUMMARY"; pass=$((pass+1))
else
  echo "FAIL  08_nearfield_wholeregion.R all globalheat   (see $gh_log)" | tee -a "$SUMMARY"; fail=$((fail+1))
fi

# WITHOUT-mixed 2-class variant of the apical-class report (2nd invocation of 13 with nomixed)
echo; echo ">>> R/11_per_organoid_final/13_citrate_gradient_report_3class.R all nomixed"
nm_log="$LOGDIR/R__11_per_organoid_final__13_citrate_gradient_report_3class_nomixed.log"
if "$RS" R/11_per_organoid_final/13_citrate_gradient_report_3class.R all nomixed > "$nm_log" 2>&1; then
  echo "PASS  13_citrate_gradient_report_3class.R all nomixed" | tee -a "$SUMMARY"; pass=$((pass+1))
else
  echo "FAIL  13_citrate_gradient_report_3class.R all nomixed   (see $nm_log)" | tee -a "$SUMMARY"; fail=$((fail+1))
fi

echo "" | tee -a "$SUMMARY"
echo "TOTAL: $pass passed, $fail failed" | tee -a "$SUMMARY"
echo "Logs in $LOGDIR/  |  figures in figures/  |  results in results/"
echo "############ run_all.sh COMPLETE ############"
