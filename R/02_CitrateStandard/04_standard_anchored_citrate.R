#!/usr/bin/env Rscript
# 02_CitrateStandard / 04_standard_anchored_citrate.R
# Extract tissue citrate at the AUTHENTIC citrate mass measured in the standard,
# sweeping the window half-width over 5 / 7 / 10 ppm so the citrate-vs-co-isobar
# trade-off is visible. All acquisition parameters were identical between the
# standard and tissue runs (user), so the standard mass scale transfers.
#
#   anchor = intensity-weighted [M-H]- centroid of the top-conc standard spot
#   windows = anchor +- {5,7,10} ppm
# The co-isobar shoulder is at ~191.0214 (+10 ppm): a +-5 ppm window excludes it,
# +-7 ppm just reaches it, +-10 ppm INCLUDES it (high edge 191.0216).
#
# NOTE: data are CENTROIDED -> each pixel has ONE merged 191 centroid (citrate +
# co-isobar unresolved, R~2400). A narrow window therefore SELECTS the pixels
# whose merged centroid is citrate-side; it does not integrate citrate area. As
# the window widens past the shoulder it starts admitting shoulder-side pixels.
#
# Out: figures/citrate_standard/04_standard_anchored_citrate.pdf
#      results/citrate_standard/citrate_anchor.csv

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/02_CitrateStandard/00_config.R"))
suppressPackageStartupMessages(library(viridisLite))
TAG <- "02.04"
OUT_PDF <- file.path(OUT_FIG, "04_standard_anchored_citrate_5-7-10ppm.pdf")
OUT_CSV <- file.path(OUT_RES, "citrate_anchor.csv")

WIDTHS     <- c(5, 7, 10)               # window half-widths (ppm) to compare
OLD_CENTER <- 191.0217                  # blended SCiLS "citrate" feature (downstream)
SHOULDER   <- 191.0214                  # measured tissue co-isobar shoulder (+10 ppm)
ISOBAR     <- 191.0461                  # dominant +138 ppm isobar
HZONE      <- c(191.0210, 191.0220)     # "shoulder-side" zone (contamination flag)
HSPAN      <- c(190.985, 191.060)
wcol       <- c("#1a7a3a", "#d08a1f", "#b3331f")   # 5 / 7 / 10 ppm colours

# ---- 1. anchor = standard citrate centroid (top-conc spot) ------------------
topf <- SPOTS$file[which.max(SPOTS$conc_M)]
sd0  <- read_file(STD_DIR, topf, wins=list(mz_win(C12,30)), hrange=HSPAN)
sel  <- sd0$hctr >= 191.013 & sd0$hctr <= 191.026
ANCHOR <- sum(sd0$hctr[sel]*sd0$hw[sel]) / sum(sd0$hw[sel])
WIN <- lapply(WIDTHS, function(p) c(ANCHOR*(1-p/1e6), ANCHOR*(1+p/1e6)))   # nested
log_cfg(TAG, "anchor = %.5f (%+.2f ppm vs theo)", ANCHOR, ppm_of(ANCHOR, C12))

# ---- GATE self-check: measured standard anchor must match the locked constant -
# CITRATE_ANCHOR_MZ / CITRATE_WIN_PPM (lib_paths.R) drive ALL pipeline citrate
# extraction. The freshly-measured standard centroid must agree; drift => re-lock.
.dppm <- ppm_of(ANCHOR, CITRATE_ANCHOR_MZ)
log_cfg(TAG, "GATE anchor vs locked CITRATE_ANCHOR_MZ %.5f: %+.2f ppm | locked window +-%d ppm",
        CITRATE_ANCHOR_MZ, .dppm, CITRATE_WIN_PPM)
if (abs(.dppm) > 2)
  warning(sprintf("[%s] GATE FAIL: standard anchor %.5f drifts %+.2f ppm from locked CITRATE_ANCHOR_MZ (>2 ppm) -- re-lock lib_paths.R / investigate before downstream image analysis.",
                  TAG, ANCHOR, .dppm))
# The pipeline-locked window is +-CITRATE_WIN_PPM (=7 ppm); the 5/7/10 sweep below
# is the supporting evidence for that choice (max citrate capture, 0% shoulder).
for (i in seq_along(WIDTHS))
  log_cfg(TAG, "  +-%2d ppm = %.5f - %.5f  (high edge %+.1f ppm vs shoulder)",
          WIDTHS[i], WIN[[i]][1], WIN[[i]][2], ppm_of(WIN[[i]][2], SHOULDER))

# windows passed to the reader: 3 widths, then shoulder-zone, then OLD feature
WINS <- c(WIN, list(HZONE, c(OLD_CENTER*(1-5/1e6), OLD_CENTER*(1+5/1e6))))
iHZ <- 4; iOLD <- 5

# ---- 2. read tissue sections ------------------------------------------------
TIS <- active_tissue()
log_cfg(TAG, "tissue side: %s (%d sections)", TIS$mode, length(TIS$sections))
SEC <- lapply(TIS$sections, function(s) {
  d <- read_file(TIS$dir, s$imz, wins=WINS, hrange=HSPAN)
  c(s, list(x=d$x, y=d$y, M=d$M, hctr=d$hctr, hw=d$hw))
})
tis_hw <- Reduce(`+`, lapply(SEC, `[[`, "hw")); tis_ctr <- SEC[[1]]$hctr
std_hw <- sd0$hw; std_ctr <- sd0$hctr

# per-width stats pooled over sections
stat <- lapply(seq_along(WIDTHS), function(i) {
  npx <- sum(sapply(SEC, function(s) sum(s$M[,i]>0)))
  tot <- sum(sapply(SEC, function(s) sum(s$M[,i])))
  shoulderside <- sum(sapply(SEC, function(s) sum(s$M[,i]>0 & s$M[,iHZ]>0)))
  cooccur_old  <- sum(sapply(SEC, function(s) sum(s$M[,i]>0 & s$M[,iOLD]>0)))
  list(ppm=WIDTHS[i], n_px=npx, total=tot, shoulderside=shoulderside, cooccur_old=cooccur_old)
})
for (st in stat)
  log_cfg(TAG, "  +-%2d ppm: pixels %d, shoulder-side %d (%.0f%%), co-occur w/OLD %d",
          st$ppm, st$n_px, st$shoulderside, 100*st$shoulderside/max(st$n_px,1), st$cooccur_old)

# ============================ RENDER =======================================
pal <- viridis(256)
pdf(OUT_PDF, width=11, height=8.5)

## ---- PAGE 1: window placement (3 widths) on standard + tissue spectra ------
par(mfrow=c(1,1), mar=c(4.8,4.8,4.4,1.2))
OX <- c(191.004, 191.050)
ssel <- std_ctr>=OX[1]&std_ctr<=OX[2]; tsel <- tis_ctr>=OX[1]&tis_ctr<=OX[2]
cmaxS <- max(std_hw[std_ctr>=191.010&std_ctr<=191.032]); cmaxT <- max(tis_hw[tis_ctr>=191.010&tis_ctr<=191.032])
plot(NA, xlim=OX, ylim=c(0,1.20), xaxs="i", yaxs="i", axes=FALSE, xlab="",
     ylab="normalised intensity (each to its 191-cluster max)")
# nested windows as brackets at descending heights (widest lightest at back)
for (i in rev(seq_along(WIDTHS))) rect(WIN[[i]][1], 0, WIN[[i]][2], 1.20, col=paste0(wcol[i],"14"), border=NA)
lines(std_ctr[ssel], (std_hw/cmaxS)[ssel], col="#2c7a2c", lwd=2.2)
lines(tis_ctr[tsel], (tis_hw/cmaxT)[tsel], col="#555555", lwd=1.8)
abline(v=ANCHOR, col="#1a7a3a", lwd=1.8)
abline(v=SHOULDER, col="#b3331f", lwd=1.4, lty=2)
abline(v=OLD_CENTER, col="#7d3c98", lwd=1.2, lty=3)
for (i in seq_along(WIDTHS)) { yb <- 1.13 - 0.05*(i-1)
  segments(WIN[[i]][1], yb, WIN[[i]][2], yb, col=wcol[i], lwd=2.6)
  segments(WIN[[i]][c(1,2)], yb-0.012, WIN[[i]][c(1,2)], yb+0.012, col=wcol[i], lwd=2.6)
  text(WIN[[i]][2], yb, sprintf(" +-%d", WIDTHS[i]), col=wcol[i], adj=c(0,0.5), cex=0.72, font=2, xpd=NA) }
axis(1, at=seq(191.00,191.05,0.01)); axis(2, las=1, at=seq(0,1,0.2))
pp <- seq(0,150,20); axis(3, at=C12*(1+pp/1e6), labels=sprintf("%+d",pp), col.axis="grey35", col="grey55")
mtext("ppm relative to citrate", side=3, line=1.9, cex=0.8, col="grey35"); mtext("m/z", side=1, line=2.4, cex=0.92)
title("Standard-anchored citrate window: 5 / 7 / 10 ppm half-widths", cex.main=1.2, line=3.3, adj=0)
legend("right", bty="n", cex=0.8,
       legend=c(sprintf("standard citrate (anchor) %.4f", ANCHOR), "tissue 191 (pooled)",
                "shoulder 191.0214 (+10 ppm)", "old feature 191.0217"),
       lwd=c(2.2,1.8,1.4,1.2), lty=c(1,1,2,3), col=c("#2c7a2c","#555555","#b3331f","#7d3c98"))
mtext(sprintf("+-5 excludes the +10 ppm shoulder (high edge %.4f); +-7 reaches it (%.4f); +-10 INCLUDES it (%.4f). Same params both runs -> standard mass transfers.",
        WIN[[1]][2], WIN[[2]][2], WIN[[3]][2]), side=1, line=3.7, cex=0.64, col="grey30")

## ---- PAGE 2: tissue ion images, one row per width -------------------------
render <- function(s, v, hi) {
  xs<-sort(unique(s$x)); ys<-sort(unique(s$y)); m<-matrix(0,length(xs),length(ys))
  m[cbind(match(s$x,xs),match(s$y,ys))] <- v
  par(mar=c(0.5,0.4,1.6,0.2))
  image(xs,ys,pmin(m/hi,1),col=pal,asp=1,useRaster=TRUE,zlim=c(0,1),axes=FALSE,xlab="",ylab="",main="")
  box(col="#444444",lwd=1.5)
  segments(xs[1]+0.05*diff(range(xs)), ys[1]+0.06*diff(range(ys)),
           xs[1]+0.05*diff(range(xs))+SCALE_UM/PX_UM, ys[1]+0.06*diff(range(ys)), col="white", lwd=2, xpd=NA)
}
hiW <- sapply(seq_along(WIDTHS), function(i)
  as.numeric(quantile(unlist(lapply(SEC, function(s) s$M[s$M[,i]>0, i])), IMG_CLIP_HI, na.rm=TRUE)))
layout(matrix(1:18, nrow=3, byrow=TRUE), widths=c(0.6,1,1,1,1,1))
par(oma=c(2.4,0.5,4.4,0.5))
for (i in seq_along(WIDTHS)) {
  par(mar=c(0.5,0.3,1.6,0.1)); plot.new()
  text(0.5,0.5,sprintf("+-%d ppm\n%.0f%% shoulder", WIDTHS[i], 100*stat[[i]]$shoulderside/max(stat[[i]]$n_px,1)),
       font=2, cex=0.92, col=wcol[i])
  for (si in seq_along(SEC)) { s<-SEC[[si]]; render(s, s$M[,i], hiW[i])
    if (i==1) mtext(sprintf("%s (%s)", s$sid, s$grp), side=3, line=0.3, cex=0.58,
                    col="#444444") }
}
mtext("Tissue citrate ion images at the standard-anchored mass: +-5 / +-7 / +-10 ppm windows",
      outer=TRUE, line=2.6, cex=1.08, font=2)
mtext(sprintf("viridis, linear, gamma 1.0, 500 um bar, per-row global p99.5 clip. %s tissue. CENTROIDED data: window SELECTS citrate-side pixels; '%% shoulder' = captured pixels that are also shoulder-side (191.0210-191.0220) -> contamination as the window widens past +10 ppm.",
      TIS$mode), outer=TRUE, line=1.0, cex=0.58, col="grey30")
dev.off()

# ---- CSV: anchor + per-width capture / contamination -----------------------
out <- data.frame(
  anchor_mz=round(ANCHOR,5), anchor_dppm_vs_theo=round(ppm_of(ANCHOR,C12),2),
  window_ppm=WIDTHS,
  win_lo=round(sapply(WIN,`[`,1),5), win_hi=round(sapply(WIN,`[`,2),5),
  high_edge_vs_shoulder_ppm=round(sapply(WIN, function(w) ppm_of(w[2], SHOULDER)),1),
  n_pixels=sapply(stat,`[[`,"n_px"),
  total_intensity=round(sapply(stat,`[[`,"total")),
  shoulderside_px=sapply(stat,`[[`,"shoulderside"),
  shoulderside_pct=round(100*sapply(stat,`[[`,"shoulderside")/pmax(sapply(stat,`[[`,"n_px"),1),1),
  cooccur_with_old_px=sapply(stat,`[[`,"cooccur_old"),
  tissue_mode=TIS$mode,
  note="centroided data: window selects citrate-side pixels (not area integration); shoulderside_pct = contamination from the +10 ppm co-isobar; profile spectra needed for true integration",
  stringsAsFactors=FALSE)
write.csv(out, OUT_CSV, row.names=FALSE)
log_cfg(TAG, "DONE -> %s ; %s", OUT_PDF, OUT_CSV)
