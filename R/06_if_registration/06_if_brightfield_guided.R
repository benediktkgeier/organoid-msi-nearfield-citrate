#!/usr/bin/env Rscript
# ============================================================================
# 06_if_brightfield_guided.R - brightfield A-slide annotate sheet WITH a faint
#   guide of where the MSI measured the 6 organoid sections (is_tissue footprints
#   projected via the existing nd2final transforms). MSI-derived (not DAPI). Use
#   this if the raw post-MALDI brightfield is too faint to find the organoids by
#   eye; the guide marks the MSI sections so you can draw/adjust the 6 rectangles
#   (and add any organoids MSI may have missed).
#
# Input : cache/register_if/bfmin_<msi_slide>_F4.rds, nd2final_*.rds, TISSUE_MSE
# Output: figures/if_registration/annotate_brightfield_<msi_slide>_guided.pdf
# Usage : Rscript R/06_if_registration/06_if_brightfield_guided.R [all | slide ...]
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/if_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register_if.R"))
suppressPackageStartupMessages({ library(Cardinal); library(EBImage) })

args   <- commandArgs(trailingOnly = TRUE)
slides <- if (length(args) == 0 || identical(args, "all")) IF_SLIDES$slide else args
ANNOT_BF_F <- 4L
mse <- readRDS(TISSUE_MSE); pd <- as.data.frame(pixelData(mse))

for (sl in slides) {
  msi_slide <- IF_SLIDES$msi_slide[IF_SLIDES$slide == sl]
  bmf <- file.path(IF_CACHE, sprintf("bfmin_%s_F%d.rds", msi_slide, ANNOT_BF_F))
  if (!file.exists(bmf)) { cat(sprintf("[95] missing %s (run R/94 first)\n", bmf)); next }
  bm <- readRDS(bmf); v <- norm01(bm); kb <- makeBrush(51, "disc")
  bfd <- norm01(v - filter2(v, kb/sum(kb), boundary = "replicate"), 0.01, 0.99)
  F <- ANNOT_BF_F

  # MSI section footprints -> brightfield-thumb px, colored per section number
  sids <- grep(sprintf("_%s_", msi_slide), unique(as.character(pd$sample_id)), value = TRUE)
  secn_of <- function(s) as.integer(sub(".*_sec(\\d+)[ab]?$", "\\1", s))
  pal <- c("#e41a1c","#ff7f00","#ffd400","#4daf4a","#377eb8","#984ea3")
  pdf(file.path(IF_FIG, sprintf("annotate_brightfield_%s_guided.pdf", msi_slide)),
      width = 16, height = 16*nrow(bfd)/ncol(bfd))
  par(mar = c(0,0,0,0)); plot.new(); plot.window(c(1, ncol(bfd)), c(nrow(bfd), 1), asp = 1)
  rasterImage(as.raster(bfd), 1, nrow(bfd), ncol(bfd), 1)
  for (s in sids) {
    xf <- file.path(REG_CACHE, sprintf("nd2final_%s.rds", s)); if (!file.exists(xf)) next
    B <- readRDS(xf)$B_msi_nd2; sub <- pd[as.character(pd$sample_id)==s & pd$is_tissue, c("x","y")]
    if (!nrow(sub)) next
    nd <- apply_affine(B, cbind(sub$x, sub$y)); n <- secn_of(s)
    points(nd[,1]/F, nd[,2]/F, pch = 15, cex = 0.18, col = adjustcolor(pal[((n-1)%%6)+1], 0.30))
    ctr <- colMeans(nd)/F
    text(ctr[1], ctr[2], n, col = pal[((n-1)%%6)+1], font = 2, cex = 2.2)
  }
  dev.off()
  cat(sprintf("[95] %s -> annotate_brightfield_%s_guided.pdf\n", msi_slide, msi_slide))
}
cat("[95] DONE\n")
