#!/usr/bin/env Rscript
# ============================================================================
# 09b_citrate_ion_hires.R - HIGH-RES export of the "Citrate [M-H]- 191.02" ion
#   image panel ONLY, one PNG per dataset, at 600 dpi.
#   Companion to 09_citrate_gradient_perdataset_v3.R. That LOCKED driver renders a
#   21-page PDF (1 summary + 20 per-dataset 2x3 pages); panel 2 of each per-dataset
#   page is the citrate single-ion image. Here we extract ONLY that panel for all 20
#   datasets, each as a standalone high-resolution PNG - WITH its scale bar and
#   colorbar, exactly as in the report.
#
#   The panel helpers live INSIDE the locked driver (not in R/00_lib), so to preserve
#   the FROZEN visual spec byte-for-byte they are COPIED VERBATIM below:
#     sec_data(), sdist_mat(), scale_bar_msi(), ion_panel()  (+ constants pal / SCALE_UM
#     / PX_UM / TRI_MAR / frame colour). Do NOT edit them here; edit the locked driver
#     first if the spec ever changes. Global p99.5 clip HI is computed identically, so
#     every PNG shares the SAME viridis scale as the report.
#   ONE INTENTIONAL DEVIATION from the locked report: the colorbar is labelled with a
#     simplified RELATIVE 0..1 scale ("citrate (rel.)") instead of the raw p99.5 clip
#     value, and the right margin is widened so the label never clips. The absolute
#     clip value is still printed once in the text subtitle for the record.
#
# In : cache/peaks_tissue_combined.rds (TISSUE_MSE), cache/zones_<sid>.rds  (dotted outline)
# Out: figures/gradient/citrate_ion_image_hires/citrate_191_<sid>.png  (x20)
# Usage: Rscript R/11_per_organoid_final/09b_citrate_ion_hires.R [all|<sid>]
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
suppressPackageStartupMessages({ library(Cardinal); library(png); library(viridisLite) })

# ---- constants COPIED VERBATIM from the locked 09 driver --------------------
SCALE_UM <- 100; PX_UM <- MSI_PIXEL_UM
CIT_MZ   <- 191.0217
pal      <- viridis(256)
TRI_MAR  <- c(3.0, 1.0, 2.4, 4.4)   # wider right margin than locked 09 (3.6) so the colorbar label never clips
DPI      <- 600L

# ---- backbone MSE + global citrate clip (identical to locked 09, lines 51-65)
mse <- readRDS(TISSUE_MSE)
pd  <- as.data.frame(pixelData(mse)); pd$gidx <- seq_len(nrow(pd))
mzs <- mz(mse); cit <- which.min(abs(mzs - CIT_MZ))
cat(sprintf("[09b] citrate feat %d (%.4f)\n", cit, mzs[cit]))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_citrate.R"))
val_cit <- citrate_onto_pd(pd)                                   # anchored citrate (raw imzML, TIC-norm)
HI <- as.numeric(quantile(val_cit[val_cit > 0], IMG_CLIP_HI, na.rm = TRUE))
if (!is.finite(HI) || HI <= 0) HI <- max(val_cit, na.rm = TRUE)
cat(sprintf("[09b] global p99.5 citrate clip: %.3g\n", HI))

SIDS20 <- levels(pixelData(mse)$sample_id); if (is.null(SIDS20)) SIDS20 <- sort(unique(as.character(pd$sample_id)))
ord20  <- sort(SIDS20)
FC     <- setNames(rep("#444444", length(SIDS20)), SIDS20)      # neutral per-dataset frame colour (verbatim)

a <- commandArgs(trailingOnly = TRUE)
RENDER_SIDS <- if (length(a) == 0 || tolower(a[1]) == "all") ord20 else a[1]
RENDER_SIDS <- RENDER_SIDS[RENDER_SIDS %in% SIDS20]
stopifnot(length(RENDER_SIDS) >= 1)

# ---- helpers COPIED VERBATIM from the locked 09 driver (lines 83-115) -------
sec_data <- function(sid) {
  m <- which(as.character(pd$sample_id) == sid); x <- pd$x[m]; y <- pd$y[m]
  xs <- sort(unique(x)); ys <- sort(unique(y)); ix <- match(x, xs); iy <- match(y, ys)
  fill <- function(v){ mm <- matrix(0, length(xs), length(ys)); mm[cbind(ix, iy)] <- v; mm }
  list(xs = xs, ys = ys, cit = fill(val_cit[m]))
}
sdist_mat <- function(z, xs, ys) {
  if (is.null(z)) return(NULL)
  m <- matrix(NA_real_, length(xs), length(ys), dimnames = list(xs, ys))
  m[cbind(match(z$x, xs), match(z$y, ys))] <- z$signed_dist_um; m
}
scale_bar_msi <- function(xs, ys) {
  blen <- SCALE_UM / PX_UM; xr <- max(xs); yb <- min(ys) - 0.06*diff(range(ys))
  segments(xr - blen, yb, xr, yb, col = "black", lwd = 3, xpd = NA)
  text(xr - blen/2, yb - 0.09*diff(range(ys)), sprintf("%d um", SCALE_UM), col = "black", cex = 0.55, adj = c(0.5, 1), xpd = NA)
}
ion_panel <- function(mat, xs, ys, hi, main, fc, mar, colorbar, outline = NULL) {
  sc <- pmin(mat / hi, 1); par(mar = mar)
  image(xs, ys, sc, col = pal, asp = 1, useRaster = TRUE, zlim = c(0,1),
        axes = FALSE, xlab = "", ylab = "", main = main, cex.main = 0.95, col.main = fc)
  box(col = fc, lwd = 1.6)
  if (!is.null(outline) && any(is.finite(outline)))
    contour(xs, ys, outline, levels = 0, add = TRUE, drawlabels = FALSE, col = "white", lty = 3, lwd = 0.8)
  if (colorbar) {                                 # SIMPLIFIED vs locked 09: relative 0..1 label (not the raw p99.5 clip)
    cx0<-grconvertX(1.02,"npc","user"); cx1<-grconvertX(1.05,"npc","user")
    cy0<-grconvertY(0.06,"npc","user"); cy1<-grconvertY(0.94,"npc","user")
    yb<-seq(cy0,cy1,length.out=257); rect(cx0,head(yb,-1),cx1,tail(yb,-1),col=pal,border=NA,xpd=NA)
    text(cx1, c(cy0,cy1), c("0","1"), pos=4, cex=0.72, xpd=NA, offset=0.25)
    text((cx0+cx1)/2, grconvertY(1.0,"npc","user"), "citrate\n(rel.)", cex=0.62, xpd=NA, adj=c(0.5,0.5))
  }
  scale_bar_msi(xs, ys)
}

# ---- output folder ---------------------------------------------------------
OUT_DIR <- file.path(GRAD_FIG, "citrate_ion_image_hires")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- one high-res citrate ion-image PNG for a section ----------------------
render_cit <- function(sid) {
  sd <- sec_data(sid); fc <- FC[sid]
  fz <- cache_in(sprintf("zones_%s.rds", sid))
  z  <- if (file.exists(fz)) readRDS(fz) else NULL
  sdm <- sdist_mat(z, sd$xs, sd$ys)                              # white-dotted organoid outline
  n_org <- if (is.null(z)) 0L else length(unique(z$instance[z$instance > 0]))
  cat(sprintf("[09b] %s: %d px grid %dx%d, %d organoid%s\n", sid, length(sd$xs)*length(sd$ys),
              length(sd$xs), length(sd$ys), n_org, if (n_org==1) "" else "s"))

  # size to the section grid aspect (asp=1 => no distortion; this only trims letterboxing).
  asp_img <- length(sd$xs) / length(sd$ys)
  if (asp_img >= 1) { w_core <- 8.5; h_core <- 8.5 / asp_img } else { h_core <- 8.5; w_core <- 8.5 * asp_img }
  w_in <- w_core + 2.0                                           # left + right (colorbar) margin room
  h_in <- h_core + 1.9                                           # title + scale-bar room

  out_png <- file.path(OUT_DIR, sprintf("citrate_191_%s.png", sid))
  png(out_png, width = w_in, height = h_in, units = "in", res = DPI, type = "cairo")
  par(oma = c(1.4, 0.4, 3.2, 0.4))
  ion_panel(sd$cit, sd$xs, sd$ys, HI,
            "Citrate [M-H]- 191.02  (white dotted = organoid outline)", fc, TRI_MAR, TRUE, outline = sdm)
  mtext(sprintf("%s - Citrate [M-H]- 191.02   (%d organoid%s)",
                disp_id(sid), n_org, if (n_org==1) "" else "s"),
        outer = TRUE, font = 2, cex = 1.25, line = 1.4, col = fc)
  mtext(sprintf("viridis, GLOBAL p99.5 clip 0..%.3g (across all 20 datasets); 100 um scale bar; white-dotted organoid segmentation outline (consensus-curated)",
                HI),
        outer = TRUE, cex = 0.72, line = 0.1, col = "grey30")
  dev.off()
  cat(sprintf("[09b]   -> %s (%.1f x %.1f in @ %d dpi)\n", out_png, w_in, h_in, DPI))
}

for (sid in RENDER_SIDS) render_cit(sid)
cat(sprintf("[09b] DONE -> %s (%d PNG%s @ %d dpi, GLOBAL p99.5 clip %.3g)\n",
            OUT_DIR, length(RENDER_SIDS), if (length(RENDER_SIDS)==1) "" else "s", DPI, HI))
