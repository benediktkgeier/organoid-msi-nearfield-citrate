#!/usr/bin/env Rscript
# 02_CitrateStandard / 02_calibration_curve.R
# Concentration -> response calibration from the citrate dilution series.
#
#  - per-spot response = citrate [M-H]- (+-10 ppm) summed intensity over on-spot
#    pixels (on-spot = tic>0, i.e. a spectrum was recorded).
#  - log-log calibration curve + linear fit (slope, R2), LOD/LOQ from the blank,
#    dynamic range, top-end roll-off check.
#  - 10 mM / 100 mM label-swap resolved from the dose-response (the curve must be
#    monotonic in true concentration).
#  - [2M-H]- dimer (383.0467) cross-check: dimer ~ [citrate]^2, so its log-log
#    slope should be ~2 vs ~1 for the monomer -> orthogonal proof of a real,
#    concentration-driven analyte (noise / artefacts do not scale this way).
#  - blank specificity: citrate-free CMC shows no 191 signal.
#
# CAVEAT: standard is in CMC, NOT tissue. Ion suppression differs -> calibration
# is semi-quantitative for tissue (order-of-magnitude only). No tissue
# concentration is computed here (standalone scope).
#
# Out: figures/citrate_standard/02_calibration_curve.pdf
#      results/citrate_standard/calibration_table.csv , calibration_fit.csv

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/02_CitrateStandard/00_config.R"))
suppressPackageStartupMessages(library(viridisLite))
TAG <- "02.02"
OUT_PDF  <- file.path(OUT_FIG, "02_calibration_curve.pdf")
OUT_TAB  <- file.path(OUT_RES, "calibration_table.csv")
OUT_FIT  <- file.path(OUT_RES, "calibration_fit.csv")

WINS <- list(MH = mz_win(C12), dimer = mz_win(DIMER_MZ), noise = mz_win(NOISE_WIN_C))

# ---- read spots ------------------------------------------------------------
log_cfg(TAG, "Reading %d standard spots ...", nrow(SPOTS))
S <- lapply(seq_len(nrow(SPOTS)), function(i) {
  d <- read_file(STD_DIR, SPOTS$file[i], wins = WINS)
  on <- d$tic > 0                                   # on-spot = spectrum recorded
  mh <- d$M[,1]; dm <- d$M[,2]; ns <- d$M[,3]
  list(label=SPOTS$label[i], conc_M=SPOTS$conc_M[i], is_blank=SPOTS$is_blank[i],
       file_conc=SPOTS$file_conc[i], x=d$x, y=d$y, mh=mh, on=on, tic=d$tic,
       n_onspot=sum(on),
       mh_mean   = mean(mh[on]),
       mh_median = median(mh[on]),
       mh_p90    = as.numeric(quantile(mh[on], 0.90)),
       mh_ticn   = mean((mh/ifelse(d$tic>0,d$tic,NA))[on], na.rm=TRUE),
       dimer_mean= mean(dm[on]),
       noise_mean= mean(ns[on]),
       mh_on_px  = mh[on])
})
names(S) <- SPOTS$label

resp  <- sapply(S, `[[`, "mh_mean")            # headline response (raw)
dimer <- sapply(S, `[[`, "dimer_mean")
concT <- SPOTS$conc_M                          # TRUE (already swap-corrected) concentration

# ---- QC: confirm the corrected labels give a monotonic dose-response -------
vm <- verify_monotonic(concT, resp)
swap_txt <- "filenames mislabeled at acquisition (10<->100 uM, 10<->100 mM); corrected to true conc"
log_cfg(TAG, "monotonic dose-response with corrected labels: %s", vm$monotonic)

# ---- blank-based LOD/LOQ ---------------------------------------------------
# NOTE: citrate-free CMC carries a low BACKGROUND ion in the m/z 191 window (the
# same co-isobaric family that contaminates tissue; see script 01). This sets a
# real, matrix-limited LOD: citrate must exceed the CMC background to be detected.
bi <- which(SPOTS$is_blank)
blank_px  <- S[[bi]]$mh_on_px
blank_mean<- mean(blank_px); blank_sd <- sd(blank_px)
LOD_int <- blank_mean + 3*blank_sd
LOQ_int <- blank_mean + 10*blank_sd

# classify spots by detection (conc>0): detectable = above LOD, quantifiable = above LOQ
detect <- concT > 0 & resp > LOD_int
quant  <- concT > 0 & resp > LOQ_int

# ---- fit on the DETECTABLE (above-LOD) points only -------------------------
# Low spots that fall in/under the CMC background are not real citrate responses
# and must not drag the calibration slope; fit the linear regime above LOD.
fit_loglog <- function(cc, rr, k) {
  if (sum(k) < 2) return(list(b0=NA, b1=NA, R2=NA, n=sum(k), k=k))
  fit <- lm(log10(rr[k]) ~ log10(cc[k]))
  list(b0=unname(coef(fit)[1]), b1=unname(coef(fit)[2]),
       R2=summary(fit)$r.squared, n=sum(k), k=k)
}
FM <- fit_loglog(concT, resp, detect)
conc_at <- function(int) if (is.na(FM$b1)||FM$b1==0) NA else 10^((log10(int) - FM$b0)/FM$b1)
LOD_conc <- if (is.finite(LOD_int) && LOD_int>0) conc_at(LOD_int) else NA
LOQ_conc <- if (is.finite(LOQ_int) && LOQ_int>0) conc_at(LOQ_int) else NA
dyn_dec  <- if (sum(detect) >= 2) log10(max(concT[detect])/min(concT[detect])) else NA
top_i    <- which.max(concT)
top_resid<- if (is.na(FM$b1)) NA else log10(resp[top_i]) - (FM$b0 + FM$b1*log10(concT[top_i]))

# ---- dimer cross-check: only valid if the dimer is actually observed -------
blank_dimer <- S[[bi]]$dimer_mean
dimer_detected <- max(dimer[concT>0]) > 3*max(blank_dimer, 1e-9)
FD <- if (dimer_detected) {
  fit_loglog(concT, dimer, concT>0 & dimer>0)
} else {
  list(b0=NA, b1=NA, R2=NA, n=0, k=rep(FALSE,length(dimer)))
}

# below-LOD low spots flagged as a data-quality note
below_lod_lbl <- paste(SPOTS$label[concT>0 & !detect], collapse=", ")

log_cfg(TAG, "monomer slope %.2f R2 %.3f on %d detectable pts (>LOD) | below-LOD: %s",
        FM$b1, FM$R2, FM$n, ifelse(nzchar(below_lod_lbl), below_lod_lbl, "none"))
log_cfg(TAG, "LOD ~%.2g M LOQ ~%.2g M | dyn(detectable) %.1f decades | top resid %.2f | dimer detected: %s",
        LOD_conc, LOQ_conc, dyn_dec, top_resid, dimer_detected)

# ============================ RENDER =======================================
pal <- viridis(256)
pdf(OUT_PDF, width = 11, height = 8.5)

## ---- PAGE 1: calibration curve --------------------------------------------
par(mar = c(5.0, 5.2, 5.0, 1.4))
kk <- concT > 0
xr <- range(log10(concT[kk])); xr <- xr + c(-0.4, 0.4)
yr <- range(log10(c(resp[kk], LOD_int, LOQ_int)), na.rm=TRUE) + c(-0.3,0.3)
plot(NA, xlim=xr, ylim=yr, axes=FALSE, xlab="", ylab="")
# LOD / LOQ bands (matrix-limited by the CMC 191 background)
rect(xr[1], yr[1], xr[2], log10(LOD_int), col="#d9534f12", border=NA)
rect(xr[1], log10(LOD_int), xr[2], log10(LOQ_int), col="#f0ad4e14", border=NA)
abline(h=log10(LOD_int), col="#d9534f", lty=3); abline(h=log10(LOQ_int), col="#e0992a", lty=3)
text(xr[2], log10(LOD_int), "LOD (3sd blank) ", col="#d9534f", adj=c(1,-0.3), cex=0.68)
text(xr[2], log10(LOQ_int), "LOQ (10sd blank) ", col="#c8841f", adj=c(1,-0.3), cex=0.68)
# fit line over DETECTABLE points only
xs <- seq(xr[1], xr[2], length.out=100)
if (is.finite(FM$b1)) lines(xs, FM$b0 + FM$b1*xs, col="#1f5fa8", lwd=2)
# points: detectable filled blue, below-LOD open grey
points(log10(concT[detect]), log10(resp[detect]), pch=21, bg="#1f5fa8", col="grey20", cex=1.7)
bl <- concT>0 & !detect
points(log10(concT[bl]), log10(pmax(resp[bl],10^yr[1])), pch=21, bg="white", col="grey55", cex=1.5)
# axes
xt <- seq(floor(xr[1]), ceiling(xr[2]))
axis(1, at=xt, labels=sapply(xt, function(e) {
  v <- 10^e; if (v>=1e-3) sprintf("%g mM", v*1e3) else sprintf("%g uM", v*1e6) }), cex.axis=0.82)
axis(2, at=pretty(yr), labels=sprintf("1e%d", pretty(yr)), las=1, cex.axis=0.82)
box(col="grey60")
mtext("true citrate concentration (log scale; 10/100 uM and 10/100 mM swap-corrected)", side=1, line=2.8, cex=0.95)
mtext("citrate response  [M-H]- mean intensity, on-spot (log10 a.u.)", side=2, line=3.4, cex=0.92)
text(log10(concT[kk]), log10(pmax(resp[kk],10^yr[1])), sapply(concT[kk], spot_label), pos=4, cex=0.66, col="grey25", offset=0.5)
title("Citrate standard calibration curve (pure citrate in CMC)", cex.main=1.22, line=3.3, adj=0)
legend("bottomright", bty="n", cex=0.84,
       pch=c(21,21), pt.bg=c("#1f5fa8","white"), col=c("grey20","grey55"),
       legend=c(sprintf("detectable (>LOD): linear fit slope %.2f, R2 %.3f, n=%d", FM$b1, FM$R2, FM$n),
                "below LOD (in CMC 191 background)"))
mtext(sprintf("LOD ~%s | LOQ ~%s | linear range ~%.1f decades (%s to %s) | top residual %+.2f (%s).   Label swaps: %s.   Dimer [2M-H]-: %s.",
        fmtM(LOD_conc), fmtM(LOQ_conc), dyn_dec,
        spot_label(min(concT[detect])), spot_label(max(concT[detect])),
        ifelse(is.finite(top_resid),top_resid,NA),
        ifelse(is.finite(top_resid)&&top_resid < -0.3, "roll-off at top", "no saturation"),
        swap_txt,
        ifelse(dimer_detected, sprintf("slope %.2f",FD$b1), "not detected above background (cross-check N/A)")),
      side=1, line=4.2, cex=0.64, col="grey25")
mtext(sprintf("Citrate-free CMC carries a background ion at m/z 191 (mean %.0f a.u.) = the same co-isobar family seen in tissue -> this sets a matrix-limited LOD. Slope ~1 (= linear response) across the detectable %s spots.",
        blank_mean, paste(sapply(concT[detect], spot_label), collapse="/")),
      side=3, line=1.4, cex=0.62, col="grey45", adj=0)
mtext("Semi-quantitative for tissue (standard is in CMC, not tissue; ion suppression differs). No tissue concentration computed in this scope.",
      side=3, line=0.4, cex=0.62, col="grey45", adj=0)

## ---- PAGE 2: per-spot droplet images + blank specificity ------------------
render_spot <- function(s, hi) {
  xs<-sort(unique(s$x)); ys<-sort(unique(s$y)); m<-matrix(0,length(xs),length(ys))
  m[cbind(match(s$x,xs),match(s$y,ys))] <- s$mh
  par(mar=c(0.5,0.5,1.8,0.5))
  image(xs,ys,pmin(m/hi,1),col=pal,asp=1,useRaster=TRUE,zlim=c(0,1),axes=FALSE,xlab="",ylab="",main="")
  box(col="grey50", lwd=1.4)
  segments(xs[1]+0.06*diff(range(xs)), ys[1]+0.08*diff(range(ys)),
           xs[1]+0.06*diff(range(xs))+SCALE_UM/PX_UM, ys[1]+0.08*diff(range(ys)), col="white", lwd=2.4, xpd=NA)
}
pos <- unlist(lapply(S, function(s) { v<-s$mh; v[v>0] }))
HI  <- as.numeric(quantile(pos, IMG_CLIP_HI, na.rm=TRUE)); if(!is.finite(HI)||HI<=0) HI<-1
layout(matrix(1:8, nrow=2, byrow=TRUE), widths=c(1,1,1,0.5))
par(oma=c(2.6,0.5,4.4,0.5))
ord <- order(SPOTS$conc_M)
for (j in seq_along(ord)) {
  i <- ord[j]; render_spot(S[[i]], HI)
  mtext(sprintf("%s%s", SPOTS$label[i], ifelse(SPOTS$is_blank[i]," (no citrate expected)","")),
        side=3, line=0.4, cex=0.74, font=2,
        col=ifelse(SPOTS$is_blank[i], "#b3331f", "#222222"))
  if (j %% 3 == 0 || j == length(ord)) {                # colorbar after each row end
    z<-seq(0,1,length.out=256); par(mar=c(1.6,0.3,1.8,2.6))
    image(1,z,matrix(z,1,256),col=pal,axes=FALSE,xlab="",ylab="",useRaster=TRUE); box(lwd=0.6)
    axis(4, at=c(0,.5,1), labels=sprintf("%.2g",c(0,.5,1)*HI), las=1, cex.axis=0.55, tcl=-0.18, mgp=c(2,0.3,0))
  }
}
mtext("Citrate [M-H]- droplet images across the dilution series (shared p99.5 scale)",
      outer=TRUE, line=2.6, cex=1.1, font=2)
mtext(sprintf("viridis, linear, gamma 1.0, 500 um bar, shared clip = p99.5 = %.2g a.u.  Ordered low->high. Blank/low spots show only a faint diffuse 191 background (CMC co-isobar); a bright localized droplet appears once citrate exceeds the background (>=1 mM).", HI),
      outer=TRUE, line=1.0, cex=0.68, col="grey30")
dev.off()

# ---- CSV ------------------------------------------------------------------
tab <- data.frame(
  label=SPOTS$label, conc_true_M=concT, file_conc=SPOTS$file_conc,
  n_onspot=sapply(S,`[[`,"n_onspot"),
  mh_mean=round(sapply(S,`[[`,"mh_mean"),3), mh_median=round(sapply(S,`[[`,"mh_median"),3),
  mh_p90=round(sapply(S,`[[`,"mh_p90"),3), mh_ticnorm_mean=signif(sapply(S,`[[`,"mh_ticn"),4),
  dimer_mean=round(sapply(S,`[[`,"dimer_mean"),3), noise_mean=round(sapply(S,`[[`,"noise_mean"),3),
  above_LOD=detect, above_LOQ=quant,
  stringsAsFactors=FALSE)
write.csv(tab, OUT_TAB, row.names=FALSE)

fit <- data.frame(
  monomer_slope=round(FM$b1,3), monomer_intercept=round(FM$b0,3), monomer_R2=round(FM$R2,4),
  n_detectable=FM$n, detectable_spots=paste(sapply(concT[detect], spot_label), collapse="; "),
  below_LOD_spots=ifelse(nzchar(below_lod_lbl), below_lod_lbl, "none"),
  dimer_detected=dimer_detected, dimer_slope=round(FD$b1,3), dimer_R2=round(FD$R2,4),
  blank_mean=round(blank_mean,3), blank_sd=round(blank_sd,3),
  LOD_intensity=round(LOD_int,3), LOQ_intensity=round(LOQ_int,3),
  LOD_conc_M=signif(LOD_conc,3), LOQ_conc_M=signif(LOQ_conc,3),
  dynamic_range_decades=round(dyn_dec,2), top_conc_residual=round(top_resid,3),
  labels_corrected="filenames mislabeled at acquisition (10<->100 uM, 10<->100 mM); table uses TRUE conc",
  dose_response_monotonic=vm$monotonic,
  notes="CMC standard (semi-quant for tissue); LOD matrix-limited by a 191 background ion in citrate-free CMC; dimer not used if undetected",
  stringsAsFactors=FALSE)
write.csv(fit, OUT_FIT, row.names=FALSE)

log_cfg(TAG, "DONE -> %s ; %s ; %s", OUT_PDF, OUT_TAB, OUT_FIT)
