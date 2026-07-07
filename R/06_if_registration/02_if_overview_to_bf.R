#!/usr/bin/env Rscript
# ============================================================================
# 02_if_overview_to_bf.R - STEP 1: locate each hi-res IF section within its
#   whole-slide overview by FFT normalized cross-correlation on DAPI texture
#   maps (same B-slide, same stain, known x5.07 scale -> reliable). Produces the
#   per-section similarity B_hr_ov (hi-res native px -> overview native px) and a
#   report showing every section boxed in the overview.
#
#   [Method note] The originally-planned overview-DAPI <-> BF-image IoU match is
#   not usable: the MSI brightfield .nd2 thumb is low-contrast / grid-dominated
#   and the overview has gel bubbles. Instead we locate sections in the overview
#   here (STEP 1) and, in STEP 2 (R/06_if_registration/03_if_overview_to_msi.R), fit overview->BF from the correspondence
#   {section center in overview} <-> {MSI section-N footprint centroid in BF}.
#
# Input : cache/register_if/ovthumb_<slide>.rds, hrthumb_<sid_if>.rds
# Output: cache/register_if/locate_<sid_if>.rds
#         figures/if_registration/step1_locate_overview.pdf, results/if_registration/locate_summary.csv
# Usage : Rscript R/06_if_registration/02_if_overview_to_bf.R [all | sid_if ...]   (default = all sections)
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/if_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register_if.R"))
suppressPackageStartupMessages({ library(EBImage) })

args <- commandArgs(trailingOnly = TRUE)
SEC  <- if_sections()
sel  <- if (length(args) == 0 || identical(args, "all")) SEC$sid_if else args
SEC  <- SEC[SEC$sid_if %in% sel, , drop = FALSE]

ovcache <- list(); getov <- function(sl) { if (is.null(ovcache[[sl]])) ovcache[[sl]] <<- readRDS(file.path(IF_CACHE, sprintf("ovthumb_%s.rds", sl))); ovcache[[sl]] }

summ <- list(); loc <- list()
for (i in seq_len(nrow(SEC))) {
  r <- SEC[i, ]; hf <- file.path(IF_CACHE, sprintf("hrthumb_%s.rds", r$sid_if))
  if (!file.exists(hf)) { cat(sprintf("[91] no thumb for %s (run R/90 all)\n", r$sid_if)); next }
  ov <- getov(r$slide); h <- readRDS(hf)
  L  <- locate_hires_in_overview(ov, h)
  saveRDS(c(L, list(sid_if = r$sid_if, slide = r$slide, secn = r$secn, msi_slide = r$msi_slide,
                    SX_hr = h$SX, SY_hr = h$SY)),
          file.path(IF_CACHE, sprintf("locate_%s.rds", r$sid_if)))
  summ[[r$sid_if]] <- data.frame(sid_if = r$sid_if, slide = r$slide, secn = r$secn,
                                 ncc = round(L$ncc, 3), ang = L$ang,
                                 ov_x = round(L$center_ov[1]), ov_y = round(L$center_ov[2]))
  loc[[r$sid_if]] <- L
  cat(sprintf("[91] %s: NCC=%.3f ang=%d center_ov=(%.0f,%.0f)\n", r$sid_if, L$ncc, L$ang, L$center_ov[1], L$center_ov[2]))
}
summ_df <- do.call(rbind, summ); write.csv(summ_df, file.path(IF_RES, "locate_summary.csv"), row.names = FALSE)

# ---- report: each section boxed in its overview ---------------------------
pdf(file.path(IF_FIG, "step1_locate_overview.pdf"), width = 13, height = 5.5)
for (sl in unique(SEC$slide)) {
  ov <- getov(sl); secs <- SEC[SEC$slide == sl, ]
  par(mar = c(1, 1, 3, 1))
  plot.new(); plot.window(c(1, ncol(ov$m16)), c(nrow(ov$m16), 1), asp = 1)
  rasterImage(as.raster(norm01(ov$m16)), 1, nrow(ov$m16), ncol(ov$m16), 1)
  pal <- rainbow(nrow(secs))
  for (j in seq_len(nrow(secs))) {
    L <- loc[[secs$sid_if[j]]]; if (is.null(L)) next
    b <- L$box; rect(b["c0"], b["r0"], b["c0"]+b["tw"]-1, b["r0"]+b["th"]-1, border = pal[j], lwd = 2.5)
    text(b["c0"]+b["tw"]/2, b["r0"]+b["th"]/2, sprintf("sec%d\n%.2f", secs$secn[j], L$ncc), col = pal[j], cex = 0.8, font = 2)
  }
  title(sprintf("STEP1  hi-res sections located in overview %s (NCC on DAPI texture)", sl), cex.main = 1)
}
dev.off()
cat(sprintf("[91] DONE -> step1_locate_overview.pdf (%d sections)\n", nrow(summ_df)))
print(summ_df, row.names = FALSE)
