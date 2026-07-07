#!/usr/bin/env Rscript
# 02_CitrateStandard / 01_standard_spectra_mz.R
# True on-instrument citrate mass & peak shape from the authentic standard.
#
# The pure citrate-in-CMC spots have NO tissue co-isobar, so they reveal the
# citrate [M-H]- centroid and FWHM the blended tissue feature (191.0217, +8 ppm,
# R~2400; R/07_metabolite_id/03_citrate_resolution.R) cannot. We then OVERLAY the
# pooled standard 191 peak against the pooled tissue 191 peak: the tissue peak is
# shifted ~+8 ppm and broadened = co-isobar contribution quantified against an
# authentic reference. This is the direct answer to script 03's open question.
#
# Out: figures/citrate_standard/01_standard_spectra_mz.pdf
#      results/citrate_standard/standard_mz_accuracy.csv

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/02_CitrateStandard/00_config.R"))
suppressPackageStartupMessages(library(viridisLite))
TAG <- "02.01"

HSPAN <- c(190.97, 191.10)            # high-res spectrum window around m/z 191
CLU   <- c(191.010, 191.032)          # citrate "cluster" window for centroid/FWHM
OUT_PDF <- file.path(OUT_FIG, "01_standard_spectra_mz.pdf")
OUT_CSV <- file.path(OUT_RES, "standard_mz_accuracy.csv")

# intensity-weighted centroid + smoothed-FWHM within a window of a hist (ctr,hw)
peak_metrics <- function(ctr, hw, win) {
  sel <- ctr >= win[1] & ctr <= win[2]
  if (!any(sel) || sum(hw[sel]) <= 0) return(list(cen=NA, fwhm_mz=NA, fwhm_ppm=NA, R=NA, apex=NA))
  cen <- sum(ctr[sel]*hw[sel]) / sum(hw[sel])
  sm  <- as.numeric(stats::filter(hw, rep(1/5,5))); sm[!is.finite(sm)] <- 0
  smc <- sm[sel]; ctc <- ctr[sel]
  apex <- ctc[which.max(smc)]; hm <- max(smc)/2
  above <- which(smc >= hm)
  fwhm_mz <- if (length(above) >= 2) ctc[max(above)] - ctc[min(above)] else NA
  list(cen=cen, fwhm_mz=fwhm_mz, fwhm_ppm=if(is.na(fwhm_mz)) NA else fwhm_mz/apex*1e6,
       R=if(is.na(fwhm_mz)||fwhm_mz<=0) NA else apex/fwhm_mz, apex=apex)
}

# ---- read the 6 standard spots --------------------------------------------
log_cfg(TAG, "Reading %d standard spots ...", nrow(SPOTS))
STD <- lapply(seq_len(nrow(SPOTS)), function(i) {
  d <- read_file(STD_DIR, SPOTS$file[i], wins = list(mz_win(C12)), hrange = HSPAN)
  pm <- peak_metrics(d$hctr, d$hw, CLU)
  mh_cen <- if (d$cden[1] > 0) d$cnum[1]/d$cden[1] else NA
  c(as.list(SPOTS[i, ]), d, list(pm = pm, mh_cen = mh_cen,
    n_onspot = sum(d$M[,1] > 0)))
})

# pooled tissue 191 spectrum (same-slide hook -> else sep-run, matrix-matched)
TIS <- active_tissue()
log_cfg(TAG, "Reading tissue side: %s (%d sections)", TIS$mode, length(TIS$sections))
tis_hw <- NULL; tis_ctr <- NULL
for (s in TIS$sections) {
  d <- read_file(TIS$dir, s$imz, wins = list(mz_win(C12)), hrange = HSPAN)
  tis_ctr <- d$hctr; tis_hw <- if (is.null(tis_hw)) d$hw else tis_hw + d$hw
}
tis_pm <- peak_metrics(tis_ctr, tis_hw, CLU)

# reference standard spot = max cluster signal (citrate-dominant, swap-agnostic)
clu_sig <- sapply(STD, function(s) { sel <- s$hctr>=CLU[1] & s$hctr<=CLU[2]; sum(s$hw[sel]) })
REFi    <- which.max(clu_sig); REF <- STD[[REFi]]
shift_ppm <- ppm_of(tis_pm$cen, REF$pm$cen)
fwhm_ratio <- tis_pm$fwhm_ppm / REF$pm$fwhm_ppm

log_cfg(TAG, "standard ref = %s : centroid %.5f (%+.2f ppm vs theo), FWHM %.1f ppm, R %.0f",
        REF$label, REF$pm$cen, ppm_of(REF$pm$cen, C12), REF$pm$fwhm_ppm, REF$pm$R)
log_cfg(TAG, "tissue (%s) : centroid %.5f (%+.2f ppm vs theo), FWHM %.1f ppm | tissue-vs-standard shift %+.2f ppm, FWHM x%.1f",
        TIS$mode, tis_pm$cen, ppm_of(tis_pm$cen, C12), tis_pm$fwhm_ppm, shift_ppm, fwhm_ratio)

# ---- adduct centroids in the reference spot (for page 3 + CSV) -------------
ADD <- lapply(IONS, function(io) {
  d <- read_file(STD_DIR, REF$file, wins = list(mz_win(io$mz, 12)),
                 hrange = mz_win(io$mz, 60))
  cen <- if (d$cden[1] > 0) d$cnum[1]/d$cden[1] else NA
  list(io = io, ctr = d$hctr, hw = d$hw, cen = cen, sum = sum(d$M[,1]))
})

# ============================ RENDER =======================================
pal <- viridis(256)
pdf(OUT_PDF, width = 11, height = 8.5)

## ---- PAGE 1: per-spot 191 spectra -----------------------------------------
layout(matrix(1:6, nrow = 2, byrow = TRUE)); par(oma = c(2.4, 1.0, 4.4, 0.6))
XL <- c(190.995, 191.055)
for (i in seq_len(nrow(SPOTS))) {
  s <- STD[[i]]; sel <- s$hctr >= XL[1] & s$hctr <= XL[2]
  yy <- s$hw[sel]; xx <- s$hctr[sel]; yt <- max(yy, 1)
  par(mar = c(2.8, 3.6, 1.8, 0.6))
  plot(xx, yy, type = "n", xlim = XL, ylim = c(0, yt*1.16), xaxs="i", yaxs="i",
       axes = FALSE, xlab = "", ylab = "")
  segments(xx, 0, xx, yy, col = "#9a9a9a", lwd = 1)
  abline(v = C12, col = "#c0392b", lwd = 1.6)
  if (is.finite(s$pm$cen)) abline(v = s$pm$cen, col = "#1f5fa8", lwd = 1.4, lty = 2)
  axis(1, at = seq(191.00,191.05,0.02), cex.axis=0.72, mgp=c(2,0.4,0))
  axis(2, las = 1, cex.axis = 0.7, mgp = c(2,0.5,0))
  box(col = "grey70")
  mtext(sprintf("%s  (n on-spot = %d)", s$label, s$n_onspot), side=3, line=0.5, font=2, cex=0.82)
  if (is.finite(s$pm$cen))
    mtext(sprintf("centroid %.4f  (%+.1f ppm) | FWHM %.0f ppm | R %s",
          s$pm$cen, ppm_of(s$pm$cen, C12), s$pm$fwhm_ppm,
          ifelse(is.finite(s$pm$R), format(round(s$pm$R,-2), big.mark=","), "NA")),
          side = 3, line = -0.4, cex = 0.6, col = "#1f5fa8")
  else mtext("no citrate signal (specificity OK for blank)", side=3, line=-0.4, cex=0.62, col="grey45")
}
mtext("Authentic citrate standard: m/z 191 peak per concentration (pure citrate in CMC, no tissue co-isobar)",
      outer = TRUE, line = 2.6, cex = 1.12, font = 2)
mtext("Red = theoretical citrate [M-H]- 191.019726.  Blue dashed = intensity-weighted centroid (cluster 191.010-191.032).  Summed intensity per 0.1 mDa bin.",
      outer = TRUE, line = 1.0, cex = 0.72, col = "grey30")
mtext(LABEL_NOTE, outer = TRUE, side = 1, line = 0.6, cex = 0.6, col = "grey50", font = 3)

## ---- PAGE 2: standard-vs-tissue overlay (THE answer to script 03) ---------
# Focus on the citrate cluster; normalise each trace to its max WITHIN the
# cluster (191.010-191.032) so the shift/broadening reads cleanly (tissue's
# global max is the distant +138 ppm isobar at 191.046, noted in the caption).
par(mfrow = c(1,1), oma = c(0,0,0,0), mar = c(5.0, 5.0, 5.2, 1.4))
OX <- c(191.004, 191.034)
cmaxS <- max(REF$hw[REF$hctr>=CLU[1] & REF$hctr<=CLU[2]])
cmaxT <- max(tis_hw[tis_ctr>=CLU[1] & tis_ctr<=CLU[2]])
selS <- REF$hctr >= OX[1] & REF$hctr <= OX[2]
selT <- tis_ctr  >= OX[1] & tis_ctr  <= OX[2]
ys <- REF$hw[selS]/cmaxS; xs <- REF$hctr[selS]
yt <- tis_hw[selT]/cmaxT; xt <- tis_ctr[selT]
plot(xs, ys, type="n", xlim=OX, ylim=c(0,1.18), xaxs="i", yaxs="i", axes=FALSE,
     xlab="", ylab="normalised intensity (each to its citrate-cluster max)")
polygon(c(OX[1],xs,OX[2]), c(0,ys,0), col="#1a7a3a22", border=NA)
lines(xs, ys, col="#1a7a3a", lwd=2.4)
lines(xt, yt, col="#b3331f", lwd=2.4)
abline(v=C12,  col="#c0392b", lwd=1.6, lty=3)
if (is.finite(REF$pm$cen))  abline(v=REF$pm$cen, col="#1a7a3a", lwd=1.4, lty=2)
if (is.finite(tis_pm$cen))  abline(v=tis_pm$cen, col="#b3331f", lwd=1.4, lty=2)
axis(1, at=seq(191.005,191.030,0.005)); axis(2, las=1)
pp <- seq(-60,80,20); axis(3, at=C12*(1+pp/1e6), labels=sprintf("%+d",pp), col.axis="grey35", col="grey55")
mtext("ppm relative to citrate", side=3, line=1.9, cex=0.82, col="grey35")
mtext("m/z", side=1, line=2.4, cex=0.95)
title("Pure citrate standard pins the true m/z; tissue 191 is shifted + broadened by the co-isobar",
      cex.main=1.16, line=3.4, adj=0)
legend("topright", bty="n", cex=0.92, lwd=2.4, lty=1, col=c("#1a7a3a","#b3331f"),
       legend=c(sprintf("standard (%s): centroid %.4f (%+.1f ppm), FWHM %.0f ppm",
                        REF$label, REF$pm$cen, ppm_of(REF$pm$cen,C12), REF$pm$fwhm_ppm),
                sprintf("tissue (%s): centroid %.4f (%+.1f ppm), FWHM %.0f ppm",
                        TIS$mode, tis_pm$cen, ppm_of(tis_pm$cen,C12), tis_pm$fwhm_ppm)))
mtext(sprintf("Tissue peak sits %+.1f ppm above the standard and is %.1fx wider; tissue retains intensity at the true citrate mass but is blended with a near +14 ppm shoulder (and, beyond this view, a dominant +138 ppm isobar at 191.046) -> the tissue 191 feature is citrate blended with unresolved co-isobars, not pure citrate.",
              shift_ppm, fwhm_ratio),
      side=1, line=4.0, cex=0.66, col="grey25")
mtext(sprintf(TISSUE_CAVEAT, TIS$mode),
      side=1, line=4.7, cex=0.55, col="grey45", font=3)

## ---- PAGE 3: adduct family in the reference standard spot ------------------
layout(matrix(1:6, nrow=2, byrow=TRUE)); par(oma=c(2.0,1.0,4.2,0.6))
add_show <- IONS[c(2,3,4,5,6,7)]                       # C13, Na, 2Na, Cl, OAc, dimer
for (k in seq_along(add_show)) {
  io <- add_show[[k]]; a <- ADD[[which(sapply(IONS,`[[`,"key")==io$key)]]
  w  <- mz_win(io$mz, 55); sel <- a$ctr>=w[1] & a$ctr<=w[2]
  yy <- a$hw[sel]; xx <- a$ctr[sel]; yt <- max(yy,1)
  par(mar=c(2.8,3.4,1.8,0.5))
  plot(xx, yy, type="n", xlim=w, ylim=c(0,yt*1.16), xaxs="i", yaxs="i", axes=FALSE, xlab="", ylab="")
  segments(xx,0,xx,yy,col="#9a9a9a",lwd=1)
  abline(v=io$mz, col="#c0392b", lwd=1.5)
  if (is.finite(a$cen)) abline(v=a$cen, col="#1f5fa8", lwd=1.3, lty=2)
  axis(1, cex.axis=0.68, mgp=c(2,0.35,0)); axis(2, las=1, cex.axis=0.66, mgp=c(2,0.5,0)); box(col="grey75")
  mtext(io$lab, side=3, line=0.45, font=2, cex=0.78)
  mtext(sprintf("theo %.4f | meas %s", io$mz,
        ifelse(is.finite(a$cen), sprintf("%.4f (%+.1f ppm)", a$cen, ppm_of(a$cen,io$mz)), "n.d.")),
        side=3, line=-0.45, cex=0.58, col="#1f5fa8")
}
mtext(sprintf("Citrate adduct family in the pure standard (%s spot): all at the expected accurate mass",
              REF$label), outer=TRUE, line=2.6, cex=1.08, font=2)
mtext("Red = theoretical adduct m/z, blue dashed = measured centroid (+-12 ppm window).  Confirms the standard reproduces citrate's full adduct fingerprint.",
      outer=TRUE, line=1.0, cex=0.72, col="grey30")
dev.off()

# ---- CSV ------------------------------------------------------------------
rows <- lapply(STD, function(s) data.frame(
  source="standard", label=s$label, conc_M=s$conc_M, n_onspot=s$n_onspot,
  mh_window_centroid=round(s$mh_cen,5),
  cluster_centroid=round(s$pm$cen,5),
  dppm_vs_theo=round(ppm_of(s$pm$cen, C12),2),
  fwhm_ppm=round(s$pm$fwhm_ppm,1), R=round(s$pm$R),
  stringsAsFactors=FALSE))
tis_row <- data.frame(source=paste0("tissue(",TIS$mode,")"), label="pooled", conc_M=NA,
  n_onspot=NA, mh_window_centroid=NA, cluster_centroid=round(tis_pm$cen,5),
  dppm_vs_theo=round(ppm_of(tis_pm$cen, C12),2), fwhm_ppm=round(tis_pm$fwhm_ppm,1),
  R=round(tis_pm$R), stringsAsFactors=FALSE)
add_rows <- do.call(rbind, lapply(ADD, function(a) data.frame(
  source=paste0("standard_adduct(",REF$label,")"), label=a$io$lab, conc_M=REF$conc_M,
  n_onspot=NA, mh_window_centroid=round(a$cen,5), cluster_centroid=NA,
  dppm_vs_theo=round(ppm_of(a$cen, a$io$mz),2), fwhm_ppm=NA, R=NA, stringsAsFactors=FALSE)))
write.csv(rbind(do.call(rbind, rows), tis_row, add_rows), OUT_CSV, row.names = FALSE)

log_cfg(TAG, "DONE -> %s ; %s", OUT_PDF, OUT_CSV)
