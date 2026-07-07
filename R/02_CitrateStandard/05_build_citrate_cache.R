#!/usr/bin/env Rscript
# 02_CitrateStandard / 05_build_citrate_cache.R   (GATE precompute step)
# Precompute the anchored citrate per pixel for EVERY inventory sample, once, so
# downstream citrate consumers (phases 05/06/08/11) just load it instead of each
# re-reading the raw imzML. One pass; cache keyed by (sample_id, x, y).
#
#   cache/citrate_anchored_<sample_id>.rds = data.frame(x, y, cit_raw, tic)
#   where cit_raw = sum of centroid intensities in CITRATE_ANCHOR_MZ +- CITRATE_WIN_PPM,
#   tic = per-pixel total ion current (for TIC-normalisation downstream).
#
# Run this AFTER reviewing the standard-anchored report (04) and confirming the
# citrate spectral signal makes sense -- this is the gate before image analysis.

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_citrate.R"))
TAG <- "02.05"

inv <- load_inventory()
log_msg <- function(...) message(sprintf("[%s] %s", TAG, sprintf(...)))
log_msg("anchor %.5f +-%d ppm | building cache for %d samples", CITRATE_ANCHOR_MZ, CITRATE_WIN_PPM, nrow(inv))

ok <- 0L; skipped <- character(0)
for (i in seq_len(nrow(inv))) {
  sid <- inv$sample_id[i]; imz <- inv$imzml_path[i]
  if (!file.exists(imz)) { skipped <- c(skipped, sid); log_msg("  SKIP %s (imzML missing)", sid); next }
  df <- citrate_raw_pixels(sid, imz)
  attr(df, "anchor_mz")  <- CITRATE_ANCHOR_MZ
  attr(df, "win_ppm")    <- CITRATE_WIN_PPM
  attr(df, "imzml_path") <- imz
  attr(df, "built")      <- "02_CitrateStandard/05_build_citrate_cache.R"
  saveRDS(df, CITRATE_CACHE(sid))
  ok <- ok + 1L
  log_msg("  %-22s %5d px, %5d cit>0 -> %s", sid, nrow(df), sum(df$cit_raw>0), basename(CITRATE_CACHE(sid)))
}
log_msg("DONE: %d caches written%s", ok,
        if (length(skipped)) sprintf(" ; skipped: %s", paste(skipped, collapse=", ")) else "")
