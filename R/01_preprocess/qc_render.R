#!/usr/bin/env Rscript
# qc_render.R
# Standalone re-render of the Phase-1 QC PDF from the SAVED combined MSE.
# Use this when 01_preprocess.R built peaks_combined.rds successfully but the
# embedded QC render failed (e.g. the n=20 "figure margins too large" layout
# bug) -- it avoids repeating the ~90 min preprocess. Logic mirrors section 9
# of 01_preprocess.R but lays panels out in a square-ish grid (.grid_dims).
#
# Reads:  cache/peaks_combined.rds, cache/feature_filter_thresholds.rds,
#         cache/meanspec_<sid>.rds, results/offset_<sid>.csv
# Writes: figures/preprocess/qc_combined.pdf
# Usage: Rscript R/01_preprocess/qc_render.R

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_paths.R"))
suppressPackageStartupMessages({
  library(Cardinal)
  library(viridisLite)
})
ION_LIB <- "D:/R/Projects/BIG_MSI/HP_Comparative_2026/R/lib_ion_image.R"
if (file.exists(ION_LIB)) source(ION_LIB)

log_msg <- function(...) message(sprintf("[%s] %s",
                                         format(Sys.time(), "%H:%M:%S"),
                                         sprintf(...)))

inv <- load_inventory()
peaks_rds <- file.path(CACHE_DIR, "peaks_combined.rds")
log_msg("Loading %s ...", basename(peaks_rds))
mse_all <- readRDS(peaks_rds)
log_msg("Loaded MSE: %d features x %d pixels", length(mz(mse_all)), ncol(mse_all))

th <- readRDS(file.path(CACHE_DIR, "feature_filter_thresholds.rds"))
INTENSITY_CUTOFF <- th$intensity_cutoff
log_msg("Tier-2 floor (from thresholds rds): %.3g", INTENSITY_CUTOFF)

pd_all <- as.data.frame(pixelData(mse_all))
sp_mat <- spectra(mse_all)
mz_vec <- mz(mse_all)

# Per-sample lock-mass residual medians for the summary page (from offset CSVs)
all_offset_rows <- list()
for (sid in inv$sample_id) {
  fp <- file.path(RES_DIR, sprintf("offset_%s.csv", sid))
  if (file.exists(fp)) all_offset_rows[[sid]] <- read.csv(fp, stringsAsFactors = FALSE)
}

# ---- helpers (verbatim from 01_preprocess.R section 9, grid layout) ----------
.grid_dims <- function(n) {
  nc <- ceiling(sqrt(n)); nr <- ceiling(n / nc); c(nr, nc)
}
.render_xy_panel <- function(x, y, vals, main, hi_abs = NULL) {
  if (length(x) == 0) { plot.new(); title(main); return(invisible(NULL)) }
  xs <- sort(unique(x)); ys <- sort(unique(y))
  mat <- matrix(NA_real_, nrow = length(xs), ncol = length(ys))
  rownames(mat) <- xs; colnames(mat) <- ys
  mat[cbind(match(x, xs), match(y, ys))] <- as.numeric(vals)
  mat[is.na(mat)] <- 0
  if (!is.null(hi_abs) && is.finite(hi_abs) && hi_abs > 0) {
    scaled <- pmin(mat / hi_abs, 1)
    if (IMG_GAMMA != 1) scaled <- scaled ^ IMG_GAMMA
    par(mar = c(2, 1.5, 3, 4))   # compact margins for the grid
    image(xs, ys, scaled, col = viridis(256), asp = 1, main = main,
          xlab = "", ylab = "", useRaster = TRUE, zlim = c(0, 1),
          cex.main = 0.6)
    cb_x0 <- grconvertX(1.04, "npc", "user"); cb_x1 <- grconvertX(1.10, "npc", "user")
    cb_y0 <- grconvertY(0.05, "npc", "user"); cb_y1 <- grconvertY(0.95, "npc", "user")
    n_col <- 256
    y_breaks <- seq(cb_y0, cb_y1, length.out = n_col + 1)
    rect(cb_x0, head(y_breaks, -1), cb_x1, tail(y_breaks, -1),
         col = viridis(n_col), border = NA, xpd = NA)
    rect(cb_x0, cb_y0, cb_x1, cb_y1, col = NA, border = "black", lwd = 0.4, xpd = NA)
    tick_fracs <- c(0, 0.5, 1)
    tick_y <- grconvertY(0.05 + 0.9 * tick_fracs, "npc", "user")
    text(cb_x1, tick_y, sprintf("%.2g", tick_fracs * hi_abs),
         pos = 4, cex = 0.5, xpd = NA, offset = 0.15)
  } else {
    # TIC page: self-contained image + colorbar at this panel's own p99.5
    pos_vals <- mat[mat > 0]
    hi <- if (length(pos_vals)) as.numeric(quantile(pos_vals, IMG_CLIP_HI)) else max(mat)
    scaled <- if (hi > 0) pmin(mat / hi, 1) else mat
    par(mar = c(2, 1.5, 3, 4))
    image(xs, ys, scaled, col = viridis(256), asp = 1, main = main,
          xlab = "", ylab = "", useRaster = TRUE, zlim = c(0, 1), cex.main = 0.6)
    cb_x0 <- grconvertX(1.04, "npc", "user"); cb_x1 <- grconvertX(1.10, "npc", "user")
    cb_y0 <- grconvertY(0.05, "npc", "user"); cb_y1 <- grconvertY(0.95, "npc", "user")
    y_breaks <- seq(cb_y0, cb_y1, length.out = 257)
    rect(cb_x0, head(y_breaks, -1), cb_x1, tail(y_breaks, -1),
         col = viridis(256), border = NA, xpd = NA)
    text(cb_x1, grconvertY(c(0.05, 0.95), "npc", "user"),
         sprintf("%.2g", c(0, hi)), pos = 4, cex = 0.5, xpd = NA, offset = 0.15)
  }
}

qc_pdf <- file.path(FIG_DIR, "preprocess", "qc_combined.pdf")
dir.create(dirname(qc_pdf), showWarnings = FALSE, recursive = TRUE)
log_msg("Rendering grid-layout QC PDF: %s", qc_pdf)

pdf(qc_pdf, width = 14, height = 9)   # wider page for the 4x5 grid

# Page 1: per-sample TIC maps (grid)
par(mfrow = .grid_dims(length(inv$sample_id)))
for (sid in inv$sample_id) {
  msk <- pd_all$sample_id == sid
  .render_xy_panel(pd_all$x[msk], pd_all$y[msk], pd_all$TIC[msk],
                   main = sprintf("%s TIC", sid))
}

# Page 2: mean spectra overlaid
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
n_t0  <- sum(inv$group == "t0_instant")
n_t20 <- sum(inv$group == "t20_incubated")
blues <- colorRampPalette(c("#1f77b4", "#0a4f7a"))(max(n_t0, 1))
reds  <- colorRampPalette(c("#d62728", "#7f0a0a"))(max(n_t20, 1))
sample_colors <- setNames(rep(NA_character_, nrow(inv)), inv$sample_id)
sample_colors[inv$group == "t0_instant"]    <- blues[seq_len(n_t0)]
sample_colors[inv$group == "t20_incubated"] <- reds[seq_len(n_t20)]
ms_list <- lapply(inv$sample_id, function(sid)
  readRDS(file.path(CACHE_DIR, sprintf("meanspec_%s.rds", sid))))
names(ms_list) <- inv$sample_id
ymax <- max(sapply(ms_list, function(d) max(log10(pmax(d$mean, 1e-12)))))
ymin <- log10(INTENSITY_CUTOFF) - 0.3
plot(NULL, xlim = c(100, 900), ylim = c(ymin, ymax + 0.5),
     xlab = "m/z", ylab = "log10(mean intensity, TIC-norm)",
     main = sprintf("Mean spectra by sample (Tier-2 floor at log10 %.2f = %.3g)",
                    log10(INTENSITY_CUTOFF), INTENSITY_CUTOFF))
abline(h = log10(INTENSITY_CUTOFF), lty = 3, col = "gray40")
for (sid in inv$sample_id) {
  d <- ms_list[[sid]]
  lines(d$mz, log10(pmax(d$mean, 1e-12)), type = "h",
        col = sample_colors[sid], lwd = 0.5)
}
abline(v = ALIGN_ANCHORS, col = "red", lty = 2, lwd = 0.7)
legend("topright", legend = names(sample_colors), col = sample_colors,
       lwd = 2, bty = "n", cex = 0.7, ncol = 2)

# Pages 3+: anchor ion images (grid), global p99.5 clip across all samples
for (k in seq_along(ALIGN_ANCHORS)) {
  tgt <- ALIGN_ANCHORS[k]; lbl <- names(ALIGN_ANCHORS)[k]
  d_ppm_all <- abs(mz_vec - tgt) / tgt * 1e6
  best <- which.min(d_ppm_all)
  par(mfrow = .grid_dims(length(inv$sample_id)))
  if (length(best) == 0 || d_ppm_all[best] > 50) {
    for (sid in inv$sample_id) { plot.new(); title(sprintf("%s\n%s NO MATCH", sid, lbl)) }
    next
  }
  mz_obs <- mz_vec[best]
  feat_vec <- as.numeric(sp_mat[best, ])
  pos_vals <- feat_vec[feat_vec > 0]
  global_hi <- if (length(pos_vals) > 0)
    as.numeric(quantile(pos_vals, IMG_CLIP_HI, na.rm = TRUE)) else NA_real_
  for (sid in inv$sample_id) {
    msk <- pd_all$sample_id == sid
    ttl <- sprintf("%s | %s\nobs %.4f (%+.2f ppm)", sid, lbl, mz_obs,
                   (mz_obs - tgt) / tgt * 1e6)
    .render_xy_panel(pd_all$x[msk], pd_all$y[msk], feat_vec[msk],
                     main = ttl, hi_abs = global_hi)
  }
}

# Final page: summary
par(mfrow = c(1, 1), mar = c(2, 2, 2, 2)); plot.new()
summary_lines <- c(
  "Phase 1 combined preprocess (20 sections) -- QC re-render",
  sprintf("Combined MSE:   %d features x %d pixels",
          length(mz(mse_all)), ncol(mse_all)),
  sprintf("Tier-2 floor:   mean >= %.3g (Kneedle auto)", INTENSITY_CUTOFF),
  "",
  "Lock-mass residual medians (QC only, no shift applied):")
for (sid in inv$sample_id) {
  o <- all_offset_rows[[sid]]
  if (is.null(o)) next
  m  <- unique(o$median_residual_ppm)
  nv <- sum(!is.na(o$ppm_residual))
  summary_lines <- c(summary_lines,
    sprintf("  %-20s median %+.2f ppm  (%d/%d anchors)",
            sid, m, nv, length(ALIGN_ANCHORS)))
}
text(0, 1, paste(summary_lines, collapse = "\n"),
     adj = c(0, 1), family = "mono", cex = 0.8)

dev.off()
log_msg("QC PDF written: %s", qc_pdf)
