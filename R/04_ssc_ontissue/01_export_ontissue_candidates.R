#!/usr/bin/env Rscript
# 01_export_ontissue_candidates.R
# Fresh on-tissue selection Mode-B (global cross-sample p99.5 clip) export of the SCiLS
# 20-section feature set, restricted to the 3 sections that already have a
# prior (TIMSCONVERT-era) on-tissue selection annotation output. The prior annotations were
# made on a DIFFERENT picking (TIMSCONVERT, 5-section) so they don't transfer to
# the SCiLS grid -- the user re-filters these 3 from scratch before SSC.
#
# Differs from the (now deprecated) all-20-section global-clip export only in:
#   - renders ONLY the 3 target sections (not all 20)
#   - global clip = p99.5 across the 3 EXPORTED sections (self-consistent set)
#   - writes to results/peakme/ (the invalid TIMSCONVERT zips were deleted)
#   - materialises only the 3 sections' pixels (~9.8k px) -> light on RAM
#
# Output (per-sample, FLAT zips per locked on-tissue selection rules):
#   results/peakme/peakme_upload_globalclip_<sid>/   PNG dir
#   results/peakme/peakme_upload_globalclip_<sid>.zip
#   results/peakme/global_clip_per_feature_3section.rds
# Usage: Rscript R/04_ssc_ontissue/01_export_ontissue_candidates.R

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_paths.R"))
suppressPackageStartupMessages({
  library(Cardinal)
  library(viridisLite)
  library(parallel)
})

TARGET_SIDS <- c("AO_0h_sl6A_sec1a", "AO_0h_sl6A_sec4b", "AO_20h_sl4A_sec5a")
OUT_DIR       <- file.path(RES_DIR, "peakme")
PNG_W <- 720L; PNG_H <- 720L
CLIP_QUANTILE <- IMG_CLIP_HI   # 0.995
GAMMA         <- IMG_GAMMA     # 1.0

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
t_all <- Sys.time()

# ---- 1. Load 20-section MSE ------------------------------------------------
cat("[05] Loading peaks_combined.rds (20-section SCiLS)...\n")
mse <- readRDS(file.path(CACHE_DIR, "peaks_combined.rds"))
fd  <- as.data.frame(featureData(mse))
mz_vec <- mz(mse)
pd  <- as.data.frame(pixelData(mse))
all_sids <- levels(pd$sample_id); if (is.null(all_sids)) all_sids <- unique(as.character(pd$sample_id))
miss <- setdiff(TARGET_SIDS, all_sids)
if (length(miss) > 0) stop("Target sample(s) not in MSE: ", paste(miss, collapse = ", "))
n_feat <- length(mz_vec)
cat(sprintf("[05] MSE: %d features x %d pixels; exporting %d sections: %s\n",
            n_feat, ncol(mse), length(TARGET_SIDS), paste(TARGET_SIDS, collapse = ", ")))

# ---- 2. Materialise ONLY the 3 target sections' pixels ---------------------
cols <- which(as.character(pd$sample_id) %in% TARGET_SIDS)
cat(sprintf("[05] Materialising spectra for %d pixels (3 sections)...\n", length(cols)))
t0 <- Sys.time()
sp3 <- as.matrix(spectra(mse)[, cols, drop = FALSE])   # features x pixels(3 sections)
cat(sprintf("  dense slice %.2f GB (%d x %d) in %.0fs\n",
            object.size(sp3) / 1024^3, nrow(sp3), ncol(sp3),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
pd3 <- pd[cols, , drop = FALSE]          # pixelData for the 3 sections, col-aligned to sp3

# ---- 3. Global p99.5 clip per feature across the 3 EXPORTED sections --------
cat("[05] Computing global p99.5 per feature across the 3 sections...\n")
tmp <- sp3; tmp[tmp == 0] <- NA_real_
global_clip <- apply(tmp, 1, function(v) {
  if (all(is.na(v))) return(1)
  quantile(v, CLIP_QUANTILE, na.rm = TRUE)
})
rm(tmp); gc()
saveRDS(global_clip, file.path(OUT_DIR, "global_clip_per_feature_3section.rds"))
cat(sprintf("  clip head: %.3g %.3g %.3g\n", global_clip[1], global_clip[2], global_clip[3]))

# ---- 4. Per-sample layout + dense slice ------------------------------------
per_sample <- list()
for (sid in TARGET_SIDS) {
  msk <- as.character(pd3$sample_id) == sid
  x <- pd3$x[msk]; y <- pd3$y[msk]
  xs <- sort(unique(x)); ys <- sort(unique(y))
  per_sample[[sid]] <- list(
    xs = xs, ys = ys, ix = match(x, xs), iy = match(y, ys),
    nx = length(xs), ny = length(ys),
    sp_dense = sp3[, msk, drop = FALSE]
  )
  cat(sprintf("  %s: %d px -> %d x %d grid\n", sid, sum(msk), length(xs), length(ys)))
}
rm(sp3); gc()

# ---- 5. metadata template: rank by combined (all-sample) mean intensity desc
md_template <- data.frame(filename = sprintf("%.4f.png", mz_vec),
                          mz_value = mz_vec, stringsAsFactors = FALSE)
md_template$rank_mean <- fd$mean
md_template <- md_template[order(-md_template$rank_mean), , drop = FALSE]
md_template$rank <- seq_len(nrow(md_template))
md_template$rank_mean <- NULL
rownames(md_template) <- NULL

# ---- 6. Per-sample renderer (flat zip) -------------------------------------
render_sample <- function(sid, per_sample, global_clip, mz_vec,
                          md_template, out_dir, png_w, png_h, gamma) {
  suppressPackageStartupMessages(library(viridisLite))
  ps <- per_sample[[sid]]; sp_mat <- ps$sp_dense
  output_dir <- file.path(out_dir, sprintf("peakme_upload_globalclip_%s", sid))
  zip_path   <- paste0(output_dir, ".zip")
  if (dir.exists(output_dir)) {
    old <- list.files(output_dir, pattern = "\\.(png|csv)$", full.names = TRUE)
    if (length(old) > 0) file.remove(old)
  } else dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pal <- viridis(256); t0 <- Sys.time()
  for (i in seq_along(mz_vec)) {
    mat <- matrix(0, nrow = ps$nx, ncol = ps$ny)
    mat[cbind(ps$ix, ps$iy)] <- sp_mat[i, ]
    hi <- global_clip[i]; if (!is.finite(hi) || hi <= 0) hi <- 1
    scaled <- pmin(mat / hi, 1); if (gamma != 1) scaled <- scaled ^ gamma
    fn <- file.path(output_dir, sprintf("%.4f.png", mz_vec[i]))
    png(fn, width = png_w, height = png_h, units = "px", bg = "black")
    par(mar = c(0, 0, 0, 0))
    image(ps$xs, ps$ys, scaled, col = pal, asp = 1, useRaster = TRUE,
          zlim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
    dev.off()
  }
  write.csv(md_template, file.path(output_dir, "metadata.csv"), row.names = FALSE)
  old_wd <- getwd(); setwd(output_dir)
  files <- list.files(".", pattern = "\\.(png|csv)$")
  if (file.exists(zip_path)) file.remove(zip_path)
  utils::zip(zip_path, files = files, flags = "-r9X")
  setwd(old_wd)
  list(sid = sid, zip = zip_path, n_pngs = length(mz_vec),
       elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}

# ---- 7. Parallel: one worker per section -----------------------------------
cat(sprintf("\n[05] Rendering %d sections in parallel (%d ions each)...\n",
            length(TARGET_SIDS), n_feat))
cl <- makeCluster(length(TARGET_SIDS), type = "PSOCK")
clusterExport(cl, c("render_sample", "per_sample", "global_clip", "mz_vec",
                    "md_template", "OUT_DIR", "PNG_W", "PNG_H", "GAMMA"),
              envir = environment())
results <- parLapply(cl, TARGET_SIDS, function(sid)
  render_sample(sid, per_sample, global_clip, mz_vec, md_template,
                OUT_DIR, PNG_W, PNG_H, GAMMA))
stopCluster(cl)

cat(sprintf("\n[05] DONE in %.1f min - SCiLS 3-section on-tissue selection zips (Mode B globalclip):\n",
            as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
for (res in results) {
  if (file.exists(res$zip))
    cat(sprintf("  %s  (%.1f MB, %d ions, %.1f min)\n", res$zip,
                file.info(res$zip)$size / 1024^2, res$n_pngs, res$elapsed_min))
  else cat(sprintf("  %s  [MISSING]\n", res$zip))
}
