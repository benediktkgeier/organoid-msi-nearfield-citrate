#!/usr/bin/env Rscript
# 03_refilter_freq_floor.R
# Apply a spatial-coverage (frequency) floor on top of the Tier-2 Kneedle cut so
# the base MSE only carries features with enough pixels to support per-zone
# gradient statistics. User decision 2026-06-16: freq >= 0.05 (present in >=5%
# of pixels, ~4,500 of 90,760) on top of the existing mean >= 0.0138 Kneedle cut
# -> 11,230 -> 8,772 features (keeps ~96% of total ion signal).
#
# The current peaks_combined.rds (self-contained, in-memory sparse, 11,230 feat)
# already has the intensity cut, and the target is exactly its freq>=0.05 subset,
# so we subset in place (stays self-contained -- no re-realize) and re-save.
# Usage: Rscript R/01_preprocess/03_refilter_freq_floor.R

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_paths.R"))
suppressPackageStartupMessages(library(Cardinal))

PEAKS  <- file.path(CACHE_DIR, "peaks_combined.rds")
TMP    <- file.path(CACHE_DIR, "peaks_combined.rds.tmp")
THRESH <- file.path(CACHE_DIR, "feature_filter_thresholds.rds")
FREQ_FLOOR <- 0.05
log_msg <- function(...) message(sprintf("[07 %s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

log_msg("Loading %s ...", basename(PEAKS))
mse <- readRDS(PEAKS)
fd  <- as.data.frame(featureData(mse))
stopifnot("freq" %in% names(fd))
log_msg("Before: %d features x %d pixels; spectra class=%s (file-backed=%s)",
        nrow(mse), ncol(mse), class(spectra(mse))[1], !is.null(attr(spectra(mse), "path")))
log_msg("  freq: min=%.4f median=%.4f max=%.4f", min(fd$freq), median(fd$freq), max(fd$freq))

keep <- fd$freq >= FREQ_FLOOR
log_msg("Applying freq >= %.2f floor: keep %d, drop %d", FREQ_FLOOR, sum(keep), sum(!keep))

mse <- mse[keep, ]
log_msg("After: %d features x %d pixels", nrow(mse), ncol(mse))
stopifnot(is.null(attr(spectra(mse), "path")))   # still in-memory / self-contained

# self-containment proof on the subset
probe <- as.matrix(spectra(mse)[seq_len(min(3, nrow(mse))),
                                round(ncol(mse) * 0.6) + 0:4, drop = FALSE])
log_msg("Mid-matrix read OK (%s)", paste(dim(probe), collapse = "x"))

# Save (temp) -> reload verify -> atomic replace
saveRDS(mse, TMP)
mse2 <- readRDS(TMP)
stopifnot(nrow(mse2) == sum(keep), is.null(attr(spectra(mse2), "path")))
r2 <- as.matrix(spectra(mse2)[seq_len(min(3, nrow(mse2))),
                              round(ncol(mse2) * 0.6) + 0:4, drop = FALSE])
stopifnot(isTRUE(all.equal(r2, probe)))
file.rename(TMP, PEAKS)
log_msg("Replaced %s (%d features, self-contained).", basename(PEAKS), nrow(mse2))

# Record the effective filter (freq floor now 0.05; intensity Kneedle unchanged)
th <- if (file.exists(THRESH)) readRDS(THRESH) else list()
th$freq_min          <- FREQ_FLOOR
th$freq_floor_note   <- "spatial-coverage floor for gradient stats (2026-06-16)"
# keep existing intensity_cutoff (Kneedle 0.0138)
saveRDS(th, THRESH)
log_msg("Updated thresholds: freq_min=%.2f, intensity_cutoff=%.4g",
        th$freq_min, ifelse(is.null(th$intensity_cutoff), NA, th$intensity_cutoff))
log_msg("DONE. peaks_combined.rds = %.0f MB", file.info(PEAKS)$size / 1024^2)
