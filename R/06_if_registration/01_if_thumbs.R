#!/usr/bin/env Rscript
# ============================================================================
# 01_if_thumbs.R - STEP 0: cache DAPI thumbnails for IF registration + a
#   thumbnail-preview report (so each step has a report).
#
#   Per slide  : overview DAPI block-mean thumbs  F=16 (29.3 um/px, matches BF
#                nd2thumb) and F=4 (7.3 um/px, coarse locate scale).
#   Per section: hi-res DAPI block-mean thumbs     F=5 (1.81 um/px, ~ov native)
#                and F=20 (7.24 um/px, coarse locate scale).
#
# Input : LM_Bsections/*_over.nd2, *_sec*.nd2
# Output: cache/register_if/ovthumb_<slide>.rds, hrthumb_<sid_if>.rds
#         figures/if_registration/step0_thumbnails.pdf
# Usage : Rscript R/06_if_registration/01_if_thumbs.R [all | sid_if ...]   (default = pilot)
# ============================================================================

if (Sys.getenv("JAVA_HOME") == "") Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jre1.8.0_491")
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/if_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register_if.R"))
suppressPackageStartupMessages({ library(RBioFormats) })

args <- commandArgs(trailingOnly = TRUE)
SEC  <- if_sections()
sel  <- if (length(args) == 0) PILOT_SIDS else if (identical(args, "all")) SEC$sid_if else args
SEC  <- SEC[SEC$sid_if %in% sel, , drop = FALSE]
slides_needed <- unique(SEC$slide)

# ---- overview thumbs (per slide) ------------------------------------------
for (sl in slides_needed) {
  s <- IF_SLIDES[IF_SLIDES$slide == sl, ]
  of <- file.path(IF_CACHE, sprintf("ovthumb_%s.rds", sl))
  if (file.exists(of)) { cat(sprintf("[90] ov %s cached\n", sl)); next }
  cat(sprintf("[90] overview %s: block-mean DAPI (F=%d, F=%d)...\n", sl, F_OV, F_OV_C))
  m16 <- nd2_channel_block_mean(s$over, ch = DAPI_CH_DEFAULT, F = F_OV)
  m4  <- nd2_channel_block_mean(s$over, ch = DAPI_CH_DEFAULT, F = F_OV_C)
  saveRDS(list(m16 = m16$m, m4 = m4$m, F16 = F_OV, F4 = F_OV_C, SX = m16$SX, SY = m16$SY, slide = sl), of)
  cat(sprintf("     -> ovthumb_%s.rds  m16=%dx%d  m4=%dx%d\n", sl, nrow(m16$m), ncol(m16$m), nrow(m4$m), ncol(m4$m)))
}

# ---- hi-res DAPI thumbs (per section) -------------------------------------
for (i in seq_len(nrow(SEC))) {
  r <- SEC[i, ]; hf <- file.path(IF_CACHE, sprintf("hrthumb_%s.rds", r$sid_if))
  if (file.exists(hf)) { cat(sprintf("[90] hr %s cached\n", r$sid_if)); next }
  if (!file.exists(r$hr_path)) { cat(sprintf("[90] MISSING %s\n", r$hr_path)); next }
  cat(sprintf("[90] hi-res %s: block-mean DAPI (F=%d, F=%d)...\n", r$sid_if, F_HR, F_HR_C))
  m5  <- nd2_channel_block_mean(r$hr_path, ch = DAPI_CH_DEFAULT, F = F_HR)
  m20 <- nd2_channel_block_mean(r$hr_path, ch = DAPI_CH_DEFAULT, F = F_HR_C)
  saveRDS(list(m5 = m5$m, m20 = m20$m, F5 = F_HR, F20 = F_HR_C, SX = m5$SX, SY = m5$SY,
               dapi_ch = DAPI_CH_DEFAULT, nC = m5$nC, sid_if = r$sid_if), hf)
  cat(sprintf("     -> hrthumb_%s.rds  m5=%dx%d  m20=%dx%d\n", r$sid_if, nrow(m5$m), ncol(m5$m), nrow(m20$m), ncol(m20$m)))
}

# ---- report: thumbnail preview --------------------------------------------
gray_raster <- function(m) { v <- norm01(m); as.raster(v) }
pdf(file.path(IF_FIG, "step0_thumbnails.pdf"), width = 11, height = 7)
for (sl in slides_needed) {
  ov <- readRDS(file.path(IF_CACHE, sprintf("ovthumb_%s.rds", sl)))
  secs <- SEC[SEC$slide == sl, ]
  par(mfrow = c(1, 1), mar = c(1, 1, 3, 1))
  plot.new(); plot.window(c(1, ncol(ov$m16)), c(nrow(ov$m16), 1), asp = 1)
  rasterImage(gray_raster(ov$m16), 1, nrow(ov$m16), ncol(ov$m16), 1)
  title(sprintf("STEP 0  overview DAPI thumb  %s  (F=%d, %.1f um/px)  %d x %d native", sl, ov$F16, OV_UMPX*ov$F16, ov$SX, ov$SY), cex.main = 1)
  # per-section hi-res DAPI thumbs grid
  n <- nrow(secs); if (n == 0) next
  nc <- min(3, n); nr <- ceiling(n / nc)
  par(mfrow = c(nr, nc), mar = c(1, 1, 2.5, 1))
  for (j in seq_len(n)) {
    h <- readRDS(file.path(IF_CACHE, sprintf("hrthumb_%s.rds", secs$sid_if[j])))
    plot.new(); plot.window(c(1, ncol(h$m5)), c(nrow(h$m5), 1), asp = 1)
    rasterImage(gray_raster(h$m5), 1, nrow(h$m5), ncol(h$m5), 1)
    title(sprintf("%s  hi-res DAPI (F=%d, %.2f um/px)", secs$sid_if[j], h$F5, HR_UMPX*h$F5), cex.main = 0.85)
  }
}
dev.off()
cat(sprintf("[90] DONE -> step0_thumbnails.pdf  (%d slides, %d sections)\n", length(slides_needed), nrow(SEC)))
