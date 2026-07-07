#!/usr/bin/env Rscript
# 01_preprocess.R (COMBINED workflow)
# Cardinal v3 multi-sample pattern:
#   readMSIData() -> MSImagingArrays per sample
#   c(arr1, arr2, ...) -> combined arrays
#   normalize(method="tic")  ::  lazy
#   process()                ::  materializes a unified-m/z MSImagingExperiment
#
# Lock-mass QC is a DIAGNOSTIC only (per-sample residual measurement, no shift
# applied) -- the Bruker instrument lock-mass + TIMSCONVERT "use recalibrated
# data" already give us calibrated m/z. We skip Cardinal's recalibrate() for
# the first pass (audit cross-sample alignment in step 12; add if needed).
#
# Output:
#   cache/peaks_combined.rds            MSImagingExperiment, unified m/z grid,
#                                       pixelData$sample_id tags each pixel
#   cache/meanspec_<sample_id>.rds      per-sample mean spectrum data.frame
#   results/offset_<sample_id>.csv      per-sample lock-mass residuals
#   results/preprocess_log.csv              one-row log (combined run)
#   figures/preprocess/qc_combined.pdf

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_paths.R"))
suppressPackageStartupMessages({
  library(Cardinal)
  library(viridisLite)
})

ION_LIB <- "D:/R/Projects/BIG_MSI/HP_Comparative_2026/R/lib_ion_image.R"
if (file.exists(ION_LIB)) source(ION_LIB)

set.seed(42)
# Register SnowParam (Windows-only path). Cardinal v3 uses the BiocParallel
# registered default. Do NOT call setCardinalBPPARAM() — it triggers an
# internal BPPARAM-injection bug in peakAlign.
register_parallel()
cat(sprintf("[parallel] Registered SnowParam workers=%d\n", N_WORKERS))

t0 <- Sys.time()
log_msg <- function(...) message(sprintf("[%s] %s",
                                         format(Sys.time(), "%H:%M:%S"),
                                         sprintf(...)))

inv <- load_inventory()
log_msg("=== Phase 1 combined preprocess (%d samples) ===", nrow(inv))

# ---------- 1. Read each imzML as MSImagingArrays ----------
arr_list <- list()
for (i in seq_len(nrow(inv))) {
  sid <- inv$sample_id[i]
  fp  <- inv$imzml_path[i]
  log_msg("Reading %s -> %s", sid, basename(fp))
  arr_list[[sid]] <- readMSIData(fp)
}

# ---------- 2. Combine via c() ----------
log_msg("Combining %d MSImagingArrays via c()...", length(arr_list))
arr_all <- do.call(c, unname(arr_list))
log_msg("Combined arrays: %d pixels total across %d runs",
        length(arr_all), length(unique(run(arr_all))))

# ---------- 3. Normalize (TIC), lazy ----------
log_msg("Queueing TIC normalize (lazy)...")
arr_all <- normalize(arr_all, method = "tic")

# ---------- 4. process() executes the lazy ops (still MSImagingArrays) ----------
log_msg("Running process() on combined arrays...")
arr_all <- process(arr_all)
log_msg("After process(): class=%s, length=%d",
        class(arr_all)[1], length(arr_all))

# ---------- 5. Build reference m/z grid + convert Arrays -> MSImagingExperiment ----------
# Cardinal v3's estimateReferencePeaks() is tuned for profile spectra and is too
# strict on centroid input (returned 2 peaks of ~66k expected). Build the ref
# grid manually: pool every pixel's centroid m/z, single-linkage at 15 ppm.
# (This is the same logic that 10_build_reference_peaks.R used previously.)
log_msg("Building reference m/z grid: pooling per-pixel centroid m/z arrays...")
mz_per_pixel <- mz(arr_all)
n_pix_total  <- length(mz_per_pixel)
all_lengths  <- lengths(mz_per_pixel)

# --- Density-invariant reference seed (FIX 2026-06-15) -----------------------
# Single-linkage at PPM_TOL CHAINS as pixel count grows. Pooling all 90,760
# pixels' ~179M centroids leaves almost no >25 ppm gaps on the m/z axis, so the
# seed grid collapsed to ~31 mega-clusters (vs the validated 269 at the
# 5-section / 15,620-px density). The reference grid is only a COARSE SEED --
# every pixel is still binned in convert/peakAlign below -- so build it from a
# fixed-size random pixel subsample, holding the pooling density at the
# validated regime regardless of total pixel count. Rare m/z are unaffected:
# the Tier-1 freq>0.001 (~90 px) filter drops anything a 15k-px sample misses.
# set.seed(42) above makes the subsample reproducible.
N_REF_PIX <- 15000L
if (n_pix_total > N_REF_PIX) {
  ref_pix_idx <- sort(sample.int(n_pix_total, N_REF_PIX))
  log_msg("  reference seed: %d/%d pixel subsample (density-invariant)",
          N_REF_PIX, n_pix_total)
} else {
  ref_pix_idx <- seq_len(n_pix_total)
  log_msg("  reference seed: all %d pixels", n_pix_total)
}
n_peaks_total <- sum(all_lengths[ref_pix_idx])
log_msg("  total per-pixel centroid peaks pooled (seed): %d", n_peaks_total)
log_msg("  per-pixel peaks (all %d px): min=%d median=%d max=%d",
        n_pix_total, min(all_lengths),
        as.integer(median(all_lengths)), max(all_lengths))

# Pool + single-linkage at PPM_TOL (sort, take intensity-blind median per cluster)
# Cardinal v3 returns an S4 list-like object for mz(arr); plain unlist() fails.
# Pre-allocate and fill explicitly. Iterate only the subsampled seed pixels
# (mz_per_pixel[[i]] double-bracket indexing is the verified-safe accessor).
all_mz <- numeric(n_peaks_total)
pos <- 1L
for (i in ref_pix_idx) {
  v <- as.numeric(mz_per_pixel[[i]])
  n <- length(v)
  if (n > 0) {
    all_mz[pos:(pos + n - 1L)] <- v
    pos <- pos + n
  }
}
all_mz <- sort(all_mz)
log_msg("  pooled+sorted unique m/z floats: %d", length(all_mz))
gaps_ppm <- c(Inf, diff(all_mz) / head(all_mz, -1) * 1e6)
cluster_id <- cumsum(gaps_ppm >= PPM_TOL)
ref_peaks <- as.numeric(tapply(all_mz, cluster_id, median))
n_ref <- length(ref_peaks)
log_msg("Reference m/z grid: %d peaks (single-linkage at %d ppm)",
        n_ref, PPM_TOL)

log_msg("convertMSImagingArrays2Experiment(tolerance=%d ppm)...", PPM_TOL)
mse_all <- convertMSImagingArrays2Experiment(
  arr_all, mz = ref_peaks, tolerance = PPM_TOL, units = "ppm")
n_after_convert <- length(mz(mse_all))
n_pix <- ncol(mse_all)
log_msg("After convertMSImagingArrays2Experiment: %d features x %d pixels",
        n_after_convert, n_pix)

# ---------- 5b. peakAlign at 15 ppm (deduplicates the expanded grid) ----------
log_msg("peakAlign(tolerance=%d ppm) on combined MSE...", PPM_TOL)
mse_all <- peakAlign(mse_all, tolerance = PPM_TOL, units = "ppm")
n_after_align <- length(mz(mse_all))
log_msg("After peakAlign(%d ppm): %d features (reduced %.1fx from convert)",
        PPM_TOL, n_after_align, n_after_convert / n_after_align)

# ---------- 5c. Compute per-feature stats for filtering ----------
log_msg("summarizeFeatures(mean, max, nz) on aligned MSE...")
mse_all <- summarizeFeatures(
  mse_all, stat = c(mean = "mean", max = "max", nz = "nnzero"))
featureData(mse_all)$freq <- featureData(mse_all)$nz / ncol(mse_all)

# ---------- 5d. Tier-1 freq filter (drop singletons) ----------
FREQ_MIN <- 0.001  # > 0.1% of pixels (~10/10725)
keep_freq <- featureData(mse_all)$freq > FREQ_MIN
n_after_freq <- sum(keep_freq)
log_msg("Tier-1 freq filter (freq > %.4f): %d features kept (%d dropped)",
        FREQ_MIN, n_after_freq, sum(!keep_freq))
mse_all <- subsetFeatures(mse_all, freq > FREQ_MIN)

# Checkpoint: post-freq, pre-intensity. Enables future cutoff iteration without
# re-running the heavy peakAlign + summarizeFeatures steps.
freq_rds <- file.path(CACHE_DIR, "peaks_after_freq.rds")
saveRDS(mse_all, freq_rds)
log_msg("Checkpoint saved: %s (%d features)", basename(freq_rds), length(mz(mse_all)))

# ---------- 5e. Noise-floor diagnostic + Kneedle elbow on log10(mean) ----------
mean_vec <- featureData(mse_all)$mean
freq_vec <- featureData(mse_all)$freq
log_msg("Detecting noise-floor cutoff via Kneedle on log10(mean)...")

# Kneedle algorithm (Satopaa 2011): on a sorted-descending curve, the knee is
# the point of maximum distance from the line connecting endpoints. We apply
# it to log10(mean) sorted descending vs normalized rank.
.kneedle_knee <- function(y_desc) {
  n <- length(y_desc)
  if (n < 10) return(list(idx = max(1L, n %/% 2L), cutoff = y_desc[max(1L, n %/% 2L)]))
  x <- seq(0, 1, length.out = n)
  y <- (y_desc - min(y_desc)) / (max(y_desc) - min(y_desc))
  # distance from each (x,y) to the line from (0,1) to (1,0)
  # line equation: x + y - 1 = 0 ; distance = |x + y - 1| / sqrt(2)
  d <- abs(x + y - 1) / sqrt(2)
  idx <- which.max(d)
  list(idx = idx, cutoff = y_desc[idx])
}
log10_mean_desc <- sort(log10(pmax(mean_vec, 1e-12)), decreasing = TRUE)
knee <- .kneedle_knee(log10_mean_desc)
auto_intensity_cutoff <- 10 ^ knee$cutoff
log_msg("Kneedle knee at rank %d / %d -> intensity cutoff = %.3g (log10 = %.2f)",
        knee$idx, length(log10_mean_desc), auto_intensity_cutoff, knee$cutoff)

# User override: if cache/feature_filter_thresholds.rds exists, use it.
thresholds_rds <- file.path(CACHE_DIR, "feature_filter_thresholds.rds")
if (file.exists(thresholds_rds)) {
  th <- readRDS(thresholds_rds)
  log_msg("Using OVERRIDE thresholds from %s: freq_min=%.4f intensity_cutoff=%.3g",
          basename(thresholds_rds), th$freq_min, th$intensity_cutoff)
  INTENSITY_CUTOFF <- th$intensity_cutoff
} else {
  INTENSITY_CUTOFF <- auto_intensity_cutoff
  saveRDS(list(freq_min = FREQ_MIN, intensity_cutoff = INTENSITY_CUTOFF),
          thresholds_rds)
  log_msg("Saved auto thresholds: %s", thresholds_rds)
}

# Diagnostic PDF (4 panels)
fd_pdf <- file.path(FIG_DIR, "preprocess",
                    "feature_filter_diagnostic.pdf")
fd_err <- tryCatch({
pdf(fd_pdf, width = 11, height = 8.5)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
# Panel A: rank vs log10(mean) with elbow
plot(seq_along(log10_mean_desc), log10_mean_desc,
     type = "l", col = "steelblue", lwd = 1.2,
     xlab = "Feature rank (sorted desc by mean)",
     ylab = "log10(mean intensity)",
     main = sprintf("A. Rank vs log10(mean) (Kneedle elbow at rank %d)",
                    knee$idx))
abline(v = knee$idx, lty = 2, col = "red")
abline(h = knee$cutoff, lty = 2, col = "red")
points(knee$idx, knee$cutoff, pch = 19, col = "red", cex = 1.5)
text(knee$idx, knee$cutoff,
     sprintf("  cutoff=%.3g", auto_intensity_cutoff),
     pos = 4, col = "red", cex = 0.9)
# Panel B: rank vs freq
freq_desc <- sort(freq_vec, decreasing = TRUE)
plot(seq_along(freq_desc), freq_desc,
     type = "l", col = "darkorange", lwd = 1.2,
     xlab = "Feature rank (sorted desc by freq)",
     ylab = "freq (fraction of pixels nonzero)",
     log = "y",
     main = "B. Rank vs freq (log y)")
abline(h = FREQ_MIN, lty = 2, col = "red")
text(length(freq_desc) * 0.8, FREQ_MIN,
     sprintf("  freq_min=%.4f", FREQ_MIN),
     pos = 3, col = "red", cex = 0.9)
# Panel C: scatter log10(mean) vs freq
plot(freq_vec, log10(pmax(mean_vec, 1e-12)),
     pch = 16, cex = 0.3, col = adjustcolor("black", 0.4),
     xlab = "freq", ylab = "log10(mean)",
     main = "C. log10(mean) vs freq (real signal upper-right)")
abline(v = FREQ_MIN, lty = 2, col = "red")
abline(h = knee$cutoff, lty = 2, col = "red")
# Panel D: cumulative intensity fraction
cum_int <- cumsum(10 ^ log10_mean_desc) / sum(10 ^ log10_mean_desc)
plot(seq_along(cum_int), cum_int,
     type = "l", col = "purple", lwd = 1.2,
     xlab = "Top-N features (rank by mean intensity)",
     ylab = "Cumulative intensity fraction",
     main = "D. Cumulative captured intensity")
abline(v = knee$idx, lty = 2, col = "red")
abline(h = cum_int[knee$idx], lty = 2, col = "red")
text(knee$idx, cum_int[knee$idx],
     sprintf("  %.1f%% captured at knee",
             100 * cum_int[knee$idx]),
     pos = 4, col = "red", cex = 0.9)
dev.off()
NULL
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  log_msg("Diagnostic PDF FAILED: %s", conditionMessage(e))
  conditionMessage(e)
})
if (is.null(fd_err)) log_msg("Filter diagnostic: %s", fd_pdf)

# ---------- 5f. Tier-2 intensity filter ----------
keep_int <- featureData(mse_all)$mean >= INTENSITY_CUTOFF
n_after_int <- sum(keep_int)
log_msg("Tier-2 intensity filter (mean >= %.3g): %d features kept (%d dropped)",
        INTENSITY_CUTOFF, n_after_int, sum(!keep_int))
mse_all <- subsetFeatures(mse_all, mean >= INTENSITY_CUTOFF)

n_feat <- length(mz(mse_all))
log_msg(paste0("FEATURE COUNT BUDGET",
               "  per-pixel-peaks=%d  pooled=%d  ref-grid=%d",
               "  after-convert=%d  after-peakAlign=%d",
               "  after-freq=%d  after-intensity=%d  FINAL=%d"),
        n_peaks_total, length(all_mz), n_ref,
        n_after_convert, n_after_align, n_after_freq, n_after_int, n_feat)

# ---------- 5. Annotate pixelData$sample_id from run ----------
# Map basename(imzml_path) without extension -> sample_id
inv$run_name <- sub("\\.imzML$", "", basename(inv$imzml_path),
                    ignore.case = TRUE)
run_to_sid <- setNames(inv$sample_id, inv$run_name)
px_run <- as.character(pixelData(mse_all)$run)
px_sid <- run_to_sid[px_run]
if (any(is.na(px_sid))) {
  stop(sprintf("Could not map %d run values to sample_id. Unique runs: %s",
               sum(is.na(px_sid)),
               paste(unique(px_run), collapse = ", ")))
}
pixelData(mse_all)$sample_id <- factor(px_sid, levels = inv$sample_id)
log_msg("Per-sample pixel counts:")
print(table(pixelData(mse_all)$sample_id))

# ---------- 6. Per-pixel TIC (diagnostic) ----------
log_msg("Computing per-pixel TIC...")
mse_all <- summarizePixels(mse_all, stat = c(TIC = "sum"))
tic <- pixelData(mse_all)$TIC
log_msg("TIC: min=%.3g median=%.3g max=%.3g (zero=%d)",
        min(tic, na.rm = TRUE), median(tic, na.rm = TRUE),
        max(tic, na.rm = TRUE), sum(tic == 0, na.rm = TRUE))

# ---------- 7. Save combined peaks RDS ----------
peaks_rds <- file.path(CACHE_DIR, "peaks_combined.rds")
saveRDS(mse_all, peaks_rds)
log_msg("Saved combined peaks RDS: %s", peaks_rds)

# ---------- 8. Per-sample mean specs + lock-mass diagnostic ----------
# Avoid subsetting the MSE (Cardinal v3 [-method is fiddly). Extract spectra
# matrix once and slice columns via matter's column-indexing.
all_offset_rows <- list()
log_msg("Per-sample mean spec + lock-mass QC...")
sp_mat <- spectra(mse_all)            # features x pixels (matter sparse)
mz_vec <- mz(mse_all)
for (sid in inv$sample_id) {
  cols <- which(pixelData(mse_all)$sample_id == sid)
  log_msg("  %s: %d pixels", sid, length(cols))
  sub <- sp_mat[, cols, drop = FALSE]
  mean_vec <- as.numeric(rowMeans(sub, na.rm = TRUE))
  nnz_vec  <- as.integer(rowSums(sub > 0, na.rm = TRUE))
  ms_df <- data.frame(
    mz   = mz_vec,
    mean = mean_vec,
    max  = NA_real_,
    nnz  = nnz_vec,
    freq = nnz_vec / length(cols)
  )
  mean_rds <- file.path(CACHE_DIR, sprintf("meanspec_%s.rds", sid))
  saveRDS(ms_df, mean_rds)
  log_msg("    saved %s", basename(mean_rds))

  # Lock-mass residual diagnostic (QC only -- no shift applied)
  mz_grid <- ms_df$mz
  mean_int <- ms_df$mean
  ot <- data.frame(
    sample_id     = sid,
    anchor_name   = names(ALIGN_ANCHORS),
    anchor_target = unname(ALIGN_ANCHORS),
    obs_mz        = NA_real_,
    ppm_residual  = NA_real_,
    obs_intensity = NA_real_,
    stringsAsFactors = FALSE
  )
  for (k in seq_along(ALIGN_ANCHORS)) {
    tgt <- ALIGN_ANCHORS[k]
    win_idx <- which(abs(mz_grid - tgt) / tgt * 1e6 < ALIGN_TOL_PPM)
    if (length(win_idx) == 0) next
    imax <- win_idx[which.max(mean_int[win_idx])]
    obs <- mz_grid[imax]
    ot$obs_mz[k]        <- obs
    ot$ppm_residual[k]  <- (obs - tgt) / tgt * 1e6
    ot$obs_intensity[k] <- mean_int[imax]
  }
  median_res <- median(ot$ppm_residual, na.rm = TRUE)
  n_valid    <- sum(!is.na(ot$ppm_residual))
  log_msg("    lock-mass residuals (%d/%d valid, median %+.2f ppm):",
          n_valid, length(ALIGN_ANCHORS), median_res)
  for (k in seq_len(nrow(ot))) {
    if (is.na(ot$obs_mz[k])) {
      log_msg("      %-22s @ %.4f  NO MATCH within %d ppm",
              ot$anchor_name[k], ot$anchor_target[k], ALIGN_TOL_PPM)
    } else {
      log_msg("      %-22s @ %.4f -> obs %.4f, %+.2f ppm, int %.3g",
              ot$anchor_name[k], ot$anchor_target[k], ot$obs_mz[k],
              ot$ppm_residual[k], ot$obs_intensity[k])
    }
  }
  ot$median_residual_ppm <- median_res
  write.csv(ot, file.path(RES_DIR, sprintf("offset_%s.csv", sid)),
            row.names = FALSE)
  all_offset_rows[[sid]] <- ot
}

# ---------- 9. QC PDF ----------
qc_pdf <- file.path(FIG_DIR, "preprocess", "qc_combined.pdf")
dir.create(dirname(qc_pdf), showWarnings = FALSE, recursive = TRUE)
log_msg("Rendering combined QC PDF: %s", qc_pdf)

.get_xy <- function(mse_obj) {
  src <- tryCatch(as.data.frame(coord(mse_obj)), error = function(e) NULL)
  if (is.null(src) || !"x" %in% names(src) || !"y" %in% names(src)) {
    src <- as.data.frame(pixelData(mse_obj))
  }
  list(x = src$x, y = src$y)
}

.render_xy_panel <- function(x, y, vals, main, hi_abs = NULL) {
  # If hi_abs is supplied, use v3-locked GLOBAL clip (cross-sample comparison)
  # with a panel-local colorbar. Otherwise fall back to BIG_MSI render_ion_image
  # (per-panel p99.5). NA cells are filled with 0 so they render as the darkest
  # viridis color (NOT the device background).
  if (length(x) == 0) { plot.new(); title(main); return(invisible(NULL)) }
  xs <- sort(unique(x)); ys <- sort(unique(y))
  mat <- matrix(NA_real_, nrow = length(xs), ncol = length(ys))
  rownames(mat) <- xs; colnames(mat) <- ys
  mat[cbind(match(x, xs), match(y, ys))] <- as.numeric(vals)
  # Fill NA cells with 0 so they render as viridis[1] (dark purple) not bg
  mat[is.na(mat)] <- 0

  if (!is.null(hi_abs) && is.finite(hi_abs) && hi_abs > 0) {
    # v3-locked: linear, GLOBAL p99.5 clip, gamma 1.0, viridis
    scaled <- pmin(mat / hi_abs, 1)
    if (IMG_GAMMA != 1) scaled <- scaled ^ IMG_GAMMA
    par(mar = c(3, 2, 4, 6))   # wider right margin for the colorbar
    image(xs, ys, scaled, col = viridis(256), asp = 1, main = main,
          xlab = "", ylab = "", useRaster = TRUE, zlim = c(0, 1))

    # Right-side colorbar using NPC (panel-relative) coordinates
    cb_x0 <- grconvertX(1.04, "npc", "user")
    cb_x1 <- grconvertX(1.09, "npc", "user")
    cb_y0 <- grconvertY(0.05, "npc", "user")
    cb_y1 <- grconvertY(0.95, "npc", "user")
    n_col <- 256
    y_breaks <- seq(cb_y0, cb_y1, length.out = n_col + 1)
    rect(cb_x0, head(y_breaks, -1), cb_x1, tail(y_breaks, -1),
         col = viridis(n_col), border = NA, xpd = NA)
    rect(cb_x0, cb_y0, cb_x1, cb_y1, col = NA, border = "black",
         lwd = 0.4, xpd = NA)
    # tick labels at 0, 50%, 100%
    tick_fracs <- c(0, 0.5, 1)
    tick_vals  <- tick_fracs * hi_abs
    tick_y     <- grconvertY(0.05 + 0.9 * tick_fracs, "npc", "user")
    text(cb_x1, tick_y, sprintf("%.2g", tick_vals),
         pos = 4, cex = 0.65, xpd = NA, offset = 0.2)
    # label above colorbar
    text((cb_x0 + cb_x1) / 2, grconvertY(1.00, "npc", "user"),
         "linear", cex = 0.55, xpd = NA)
    text((cb_x0 + cb_x1) / 2, grconvertY(0.97, "npc", "user"),
         "p99.5", cex = 0.55, xpd = NA)
  } else if (exists("render_ion_image", mode = "function")) {
    img <- list(matrix = mat, x_coord = xs, y_coord = ys,
                mz_matched = NA, ppm_err = NA, n_nonzero = NA, n_total = NA)
    render_ion_image(img, clip_hi = IMG_CLIP_HI, gamma = IMG_GAMMA,
                     palette = viridis(256), main = main,
                     add_colorbar = TRUE)
  } else {
    image(xs, ys, mat, col = viridis(256), asp = 1, main = main,
          xlab = "x", ylab = "y", useRaster = TRUE)
  }
}

qc_err <- tryCatch({
pdf(qc_pdf, width = 11, height = 8.5)

# Square-ish panel grid: c(1, n) collapses to margins-too-large beyond ~6
# samples (each panel narrower than its margins). Works at n=5 and n=20.
.grid_dims <- function(n) {
  nc <- ceiling(sqrt(n)); nr <- ceiling(n / nc); c(nr, nc)
}
# Page 1: per-sample TIC maps in a grid (no MSE subset; use pixelData directly)
par(mfrow = .grid_dims(length(inv$sample_id)), mar = c(2, 2, 2.5, 4))
pd_all <- as.data.frame(pixelData(mse_all))
for (sid in inv$sample_id) {
  msk <- pd_all$sample_id == sid
  .render_xy_panel(pd_all$x[msk], pd_all$y[msk], pd_all$TIC[msk],
                   main = sprintf("%s TIC", sid))
}

# Page 2: mean spectra overlaid
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
# Color per sample group (0h = cool blue tones, 20h = warm red tones), scalable
n_t0  <- sum(inv$group == "t0_instant")
n_t20 <- sum(inv$group == "t20_incubated")
blues  <- colorRampPalette(c("#1f77b4", "#0a4f7a"))(max(n_t0,  1))
reds   <- colorRampPalette(c("#d62728", "#7f0a0a"))(max(n_t20, 1))
sample_colors <- setNames(rep(NA_character_, nrow(inv)), inv$sample_id)
sample_colors[inv$group == "t0_instant"]    <- blues[seq_len(n_t0)]
sample_colors[inv$group == "t20_incubated"] <- reds[seq_len(n_t20)]
ms_list <- lapply(inv$sample_id, function(sid) {
  readRDS(file.path(CACHE_DIR, sprintf("meanspec_%s.rds", sid)))
})
names(ms_list) <- inv$sample_id
ymax <- max(sapply(ms_list, function(d) max(log10(pmax(d$mean, 1e-12)))))
# Tightened y-axis: bottom = filter floor with 0.3 padding (no empty padding below)
ymin <- log10(INTENSITY_CUTOFF) - 0.3
plot(NULL, xlim = c(100, 900), ylim = c(ymin, ymax + 0.5),
     xlab = "m/z", ylab = "log10(mean intensity, TIC-norm)",
     main = sprintf("Mean spectra by sample (Tier-2 floor at log10 %.2f = %.3g)",
                    log10(INTENSITY_CUTOFF), INTENSITY_CUTOFF))
abline(h = log10(INTENSITY_CUTOFF), lty = 3, col = "gray40")
for (sid in inv$sample_id) {
  d <- ms_list[[sid]]
  lines(d$mz, log10(pmax(d$mean, 1e-12)),
        type = "h", col = sample_colors[sid], lwd = 0.5)
}
abline(v = ALIGN_ANCHORS, col = "red", lty = 2, lwd = 0.7)
legend("topright", legend = names(sample_colors),
       col = sample_colors, lwd = 2, bty = "n", cex = 0.9)

# Pages 3+: anchor ion images, per sample side by side.
# v3-locked: GLOBAL p99.5 clip across both samples per ion, viridis, linear,
# gamma 1.0. Both panels share the same colorbar -> cross-sample comparable.
for (k in seq_along(ALIGN_ANCHORS)) {
  tgt <- ALIGN_ANCHORS[k]
  lbl <- names(ALIGN_ANCHORS)[k]
  d_ppm_all <- abs(mz_vec - tgt) / tgt * 1e6
  best <- which.min(d_ppm_all)
  par(mfrow = .grid_dims(length(inv$sample_id)), mar = c(2, 2, 2.5, 4))
  if (length(best) == 0 || d_ppm_all[best] > 50) {
    for (sid in inv$sample_id) {
      plot.new(); title(sprintf("%s\n%s NO MATCH", sid, lbl))
    }
    next
  }
  mz_obs <- mz_vec[best]
  feat_vec <- as.numeric(sp_mat[best, ])

  # GLOBAL p99.5 clip across ALL samples for this ion (v3 locked rule)
  pos_vals <- feat_vec[feat_vec > 0]
  global_hi <- if (length(pos_vals) > 0) {
    as.numeric(quantile(pos_vals, IMG_CLIP_HI, na.rm = TRUE))
  } else NA_real_

  for (sid in inv$sample_id) {
    msk <- pd_all$sample_id == sid
    ttl <- sprintf("%s | %s\nobs %.4f (%+.2f ppm) | global p99.5 = %.3g",
                   sid, lbl, mz_obs, (mz_obs - tgt) / tgt * 1e6,
                   ifelse(is.finite(global_hi), global_hi, NA))
    .render_xy_panel(pd_all$x[msk], pd_all$y[msk], feat_vec[msk],
                     main = ttl, hi_abs = global_hi)
  }
}

# Final page: summary
par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))
plot.new()
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
offset_summary <- do.call(rbind, all_offset_rows)
summary_lines <- c(
  "Phase 1 combined preprocess",
  sprintf("Samples:        %d (%s)", nrow(inv),
          paste(inv$sample_id, collapse = ", ")),
  sprintf("Combined MSE:   %d features x %d pixels",
          length(mz(mse_all)), ncol(mse_all)),
  sprintf("Processing:     readMSIData -> c() -> normalize(TIC) -> process()"),
  sprintf("                NO recalibrate (lock-mass trusted; audit in step 12)"),
  "",
  "Lock-mass residual medians (QC only, no shift applied):"
)
for (sid in inv$sample_id) {
  m <- unique(all_offset_rows[[sid]]$median_residual_ppm)
  nv <- sum(!is.na(all_offset_rows[[sid]]$ppm_residual))
  summary_lines <- c(summary_lines,
    sprintf("  %-12s median %+.2f ppm  (%d/%d anchors)",
            sid, m, nv, length(ALIGN_ANCHORS)))
}
summary_lines <- c(summary_lines, "",
  "RDS outputs:",
  sprintf("  %s", peaks_rds),
  paste0("  ", file.path(CACHE_DIR,
                        sprintf("meanspec_%s.rds", inv$sample_id))),
  "",
  sprintf("Elapsed: %.1f min", elapsed))
text(0, 1, paste(summary_lines, collapse = "\n"),
     adj = c(0, 1), family = "mono", cex = 0.85)

dev.off()
NULL
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  log_msg("QC PDF FAILED: %s (continuing; RDS already saved)",
          conditionMessage(e))
  conditionMessage(e)
})
if (is.null(qc_err)) log_msg("QC PDF: %s", qc_pdf)

# ---------- 10. Log row ----------
log_csv <- file.path(RES_DIR, "phase1_log.csv")
log_entry <- data.frame(
  finished_at      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  n_samples        = nrow(inv),
  sample_ids       = paste(inv$sample_id, collapse = ";"),
  combined_features = length(mz(mse_all)),
  combined_pixels  = ncol(mse_all),
  per_sample_residual_ppm = paste(sprintf("%s=%+.2f",
                                          inv$sample_id,
                                          sapply(inv$sample_id, function(s)
                                            unique(all_offset_rows[[s]]$median_residual_ppm))),
                                  collapse = ";"),
  elapsed_min      = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
  peaks_rds        = peaks_rds,
  qc_pdf           = qc_pdf,
  stringsAsFactors = FALSE
)
write.table(log_entry, log_csv,
            append    = file.exists(log_csv),
            sep       = ",",
            col.names = !file.exists(log_csv),
            row.names = FALSE,
            quote     = TRUE)

log_msg("DONE in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
