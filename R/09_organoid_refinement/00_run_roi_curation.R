#!/usr/bin/env Rscript
# ============================================================================
# 00_run_roi_curation.R - one-shot runner for the organoid ROI-curation phase.
# ----------------------------------------------------------------------------
# Runs the APPLY half of the curation loop in the correct order, after you have
# marked the canvas in Adobe. Run this RIGHT AFTER organoid segmentation
# (R/08_organoid_gradient_survey/01_segment_organoids.R) and re-run it each
# curation iteration. See docs/roi_curation.md for the full workflow and
# marking conventions.
#
#   01_organoid_split_apply  (apply green-line cuts -> instances_split)
#   03_island_cleanup_apply  (merge/delete from organoid_island_actions.csv -> instances_clean)
#   04_finalize_instances    (one connected ROI = one id; drop specks; honour removals -> instances_final)
#   02_island_cleanup_canvas (re-render organoid_island_cleanup.pdf for the next iteration)
#
# Usage:
#   Rscript R/09_organoid_refinement/00_run_roi_curation.R [pdf1 pdf2 ...]   # pdfs passed to 01_organoid_split_apply (new cuts)
#   Rscript R/09_organoid_refinement/00_run_roi_curation.R --apply-only      # re-derive from stored strokes only
# Any extra flags (e.g. --any-color, --apply-only) are forwarded to 01_organoid_split_apply.
# ============================================================================

R_DIR   <- file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R")
RSCRIPT <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
args <- commandArgs(trailingOnly = TRUE)

run <- function(script, pass = character(0)) {
  cat(sprintf("\n=== %s %s ===\n", script, paste(pass, collapse = " ")))
  code <- system2(RSCRIPT, c(shQuote(file.path(R_DIR, script)), shQuote(pass)))
  if (code != 0) stop(sprintf("%s exited with code %d", script, code))
}

run("41_organoid_split_apply.R", args)   # forwards pdfs / --any-color / --apply-only
run("43_island_cleanup_apply.R")
run("44_finalize_instances.R")
run("42_island_cleanup_canvas.R")
cat("\n[45] ROI curation applied. Inspect figures/annotation/organoid_island_cleanup.pdf;\n")
cat("[45] mark another round and re-run, or proceed to apical scoring (R/10_apical_annotation/01_apical_annotate.R).\n")
