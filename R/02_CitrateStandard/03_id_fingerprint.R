#!/usr/bin/env Rscript
# 02_CitrateStandard / 03_id_fingerprint.R
# Orthogonal ID confidence: does tissue 191 carry citrate's isotope + adduct
# fingerprint measured on an AUTHENTIC standard?
#
#  - 13C/12C isotope ratio: theoretical ~6.5-7% for 6-carbon citrate; measured in
#    the pure standard and in tissue.
#  - adduct ratios ([M-2H+Na]-, [M-3H+2Na]-, [M+Cl]-, [M+CH3COO]-, [2M-H]- dimer)
#    each relative to [M-H]-, standard vs tissue.
#  Matching isotope ratio + matching adduct pattern between authentic standard and
#  tissue = strong orthogonal evidence the tissue citrate component is citrate.
#
# Out: figures/citrate_standard/03_id_fingerprint.pdf
#      results/citrate_standard/fingerprint_compare.csv

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/02_CitrateStandard/00_config.R"))
suppressPackageStartupMessages(library(viridisLite))
TAG <- "02.03"
OUT_PDF <- file.path(OUT_FIG, "03_id_fingerprint.pdf")
OUT_CSV <- file.path(OUT_RES, "fingerprint_compare.csv")

ENV  <- c(190.990, 192.100)                       # C12 + C13 envelope window
WINS <- c(lapply(IONS, function(io) mz_win(io$mz)), list(noise = mz_win(NOISE_WIN_C)))
KEYS <- sapply(IONS, `[[`, "key")
LABS <- sapply(IONS, `[[`, "lab")

# pooled on-spot ion sums for one file
ion_sums <- function(dir, file, with_env = FALSE) {
  d <- read_file(dir, file, wins = WINS, hrange = if (with_env) ENV else NULL)
  on <- d$tic > 0
  s  <- colSums(d$M[on, , drop = FALSE])
  out <- list(sum = setNames(s[seq_along(IONS)], KEYS), noise = s[length(WINS)], n = sum(on))
  if (with_env) { out$hctr <- d$hctr; out$hw <- d$hw }
  out
}

# ---- standard: pick reference spot = max [M-H]- signal ---------------------
log_cfg(TAG, "Reading %d standard spots ...", nrow(SPOTS))
std_all <- lapply(seq_len(nrow(SPOTS)), function(i) ion_sums(STD_DIR, SPOTS$file[i]))
mh_each <- sapply(std_all, function(z) z$sum[["MH"]])
REFi <- which.max(mh_each); REF_LBL <- SPOTS$label[REFi]
log_cfg(TAG, "reference standard spot = %s", REF_LBL)
std <- ion_sums(STD_DIR, SPOTS$file[REFi], with_env = TRUE)

# ---- tissue: pool sections (same-slide hook -> else sep-run, matrix-matched)
TIS <- active_tissue()
log_cfg(TAG, "Reading tissue side: %s (%d sections)", TIS$mode, length(TIS$sections))
tis_sum <- setNames(numeric(length(IONS)), KEYS); tis_hw <- NULL; tis_ctr <- NULL; tis_n <- 0
for (s in TIS$sections) {
  z <- ion_sums(TIS$dir, s$imz, with_env = TRUE)
  tis_sum <- tis_sum + z$sum; tis_n <- tis_n + z$n
  tis_ctr <- z$hctr; tis_hw <- if (is.null(tis_hw)) z$hw else tis_hw + z$hw
}

# ---- ratios ----------------------------------------------------------------
ratio_to_MH <- function(s) s / s[["MH"]]
r_std <- ratio_to_MH(std$sum)
r_tis <- ratio_to_MH(tis_sum)
iso_std <- r_std[["C13"]]; iso_tis <- r_tis[["C13"]]

# concordance: fold difference standard vs tissue (lower = better)
fold <- function(a, b) { a<-max(a,1e-12); b<-max(b,1e-12); max(a/b, b/a) }
conc_fold <- sapply(KEYS, function(k) fold(r_std[[k]], r_tis[[k]]))

log_cfg(TAG, "13C/12C : standard %.1f%% | tissue %.1f%% | theoretical ~%.1f%%",
        100*iso_std, 100*iso_tis, 100*ISO_THEO)
for (k in KEYS[-1]) log_cfg(TAG, "  %-6s ratio/[M-H]-: std %.4f  tis %.4f  (x%.1f)",
                            k, r_std[[k]], r_tis[[k]], conc_fold[[k]])

# ============================ RENDER =======================================
pdf(OUT_PDF, width = 11, height = 8.5)

## ---- PAGE 1: isotope envelope ---------------------------------------------
layout(matrix(c(1,2), nrow=1), widths=c(1.25,1)); par(oma=c(2.0,1.0,4.2,0.6))
# A: full envelope, log1p, normalized to C12 peak height
c12win <- function(ctr,hw){ s<-ctr>=mz_win(C12)[1]&ctr<=mz_win(C12)[2]; if(any(s)) max(hw[s]) else 1 }
ns <- std$hw/c12win(std$hctr,std$hw); nt <- tis_hw/c12win(tis_ctr,tis_hw)
par(mar=c(3.6,3.8,1.8,0.6))
plot(std$hctr, log1p(ns*100), type="n", xlim=ENV, xaxs="i",
     ylim=c(0, log1p(max(c(ns,nt),na.rm=TRUE)*100)*1.08), axes=FALSE, xlab="", ylab="")
lines(std$hctr, log1p(ns*100), col="#1a7a3a", lwd=1.6)
lines(tis_ctr, log1p(nt*100), col="#b3331f", lwd=1.6)
abline(v=C12, col="#c0392b", lwd=1.3, lty=3); abline(v=C13_MZ, col="#1f5fa8", lwd=1.3, lty=2)
axis(1, at=seq(191.0,192.0,0.5), cex.axis=0.74, mgp=c(2,0.4,0))
axis(2, las=1, cex.axis=0.7, mgp=c(2,0.5,0)); box(col="grey70")
mtext("m/z", side=1, line=2.1, cex=0.8); mtext("log1p(% of C12 peak)", side=2, line=2.5, cex=0.78)
mtext("Isotope envelope (normalised to C12)", side=3, line=0.4, font=2, cex=0.86)
legend("topright", bty="n", cex=0.78, lwd=1.6, col=c("#1a7a3a","#b3331f"),
       legend=c(sprintf("standard (%s)",REF_LBL), sprintf("tissue (%s)",TIS$mode)))
# B: C13 zoom, linear normalized, with theoretical stick
zw <- mz_win(C13_MZ, 60)
par(mar=c(3.6,3.8,1.8,0.6))
plot(NA, xlim=zw, ylim=c(0, max(0.12, iso_std, iso_tis, ISO_THEO)*1.18*100),
     axes=FALSE, xlab="", ylab="")
sS<-std$hctr>=zw[1]&std$hctr<=zw[2]; sT<-tis_ctr>=zw[1]&tis_ctr<=zw[2]
lines(std$hctr[sS], ns[sS]*100, col="#1a7a3a", lwd=1.8)
lines(tis_ctr[sT], nt[sT]*100, col="#b3331f", lwd=1.8)
abline(v=C13_MZ, col="#1f5fa8", lwd=1.3, lty=2)
segments(C13_MZ, 0, C13_MZ, ISO_THEO*100, col="#7d3c98", lwd=2.4)
points(C13_MZ, ISO_THEO*100, pch=18, col="#7d3c98", cex=1.3)
axis(1, cex.axis=0.7, mgp=c(2,0.4,0)); axis(2, las=1, cex.axis=0.7, mgp=c(2,0.5,0)); box(col="grey70")
mtext("m/z", side=1, line=2.1, cex=0.8); mtext("% of C12 peak", side=2, line=2.5, cex=0.78)
mtext("13C zoom vs theoretical", side=3, line=0.4, font=2, cex=0.86)
legend("topright", bty="n", cex=0.74, lwd=c(1.8,1.8,2.4), pch=c(NA,NA,18),
       col=c("#1a7a3a","#b3331f","#7d3c98"),
       legend=c(sprintf("std  %.1f%%", 100*iso_std), sprintf("tissue %.1f%%", 100*iso_tis),
                sprintf("theo ~%.1f%%", 100*ISO_THEO)))
mtext("Citrate 13C/12C isotope ratio: authentic standard vs theoretical, and the inflated tissue ratio (co-isobar signature)",
      outer=TRUE, line=2.4, cex=1.02, font=2)
mtext(sprintf("The standard reproduces citrate's theoretical 13C ratio (%.1f%% vs ~%.1f%%) -> the standard is citrate.  Tissue is inflated (%.1f%%): the expected signature of the unresolved co-isobar (scripts 01 & 07_metabolite_id) inside the 191/192 windows.  Tissue side = %s.",
        100*iso_std, 100*ISO_THEO, 100*iso_tis, TIS$mode),
      outer=TRUE, line=0.9, cex=0.66, col="grey30")
mtext(sprintf("Reference standard spot = %s (true conc).  %s", REF_LBL, LABEL_NOTE),
      outer=TRUE, side=1, line=0.5, cex=0.58, col="grey50", font=3)

## ---- PAGE 2: adduct-ratio fingerprint bars + verdict ----------------------
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(7.5,5.0,5.0,1.2))
ak <- KEYS[!KEYS %in% c("MH","C13")]             # adducts only ([M-H]-=ref; 13C is on page 1)
al <- LABS[match(ak, KEYS)]
m  <- rbind(standard = r_std[ak], tissue = r_tis[ak])
ymin <- min(m[m>0], na.rm=TRUE)/2; ymax <- max(m, na.rm=TRUE)*1.6
bp <- barplot(m, beside=TRUE, log="y", ylim=c(ymin, ymax), col=c("#1a7a3a","#b3331f"),
              border="grey25", names.arg=rep("",length(ak)), las=1,
              ylab="abundance ratio to [M-H]-  (log)", axes=FALSE)
axis(2, las=1, cex.axis=0.8); box(col="grey60")
ylab_y <- 10^(par("usr")[3] - 0.04*diff(par("usr")[3:4]))     # just below the axis (log scale)
text(colMeans(bp), ylab_y, labels=al, srt=30, adj=c(1,0.5), xpd=NA, cex=0.82)
for (j in seq_along(ak)) text(colMeans(bp)[j], ymax, sprintf("x%.1f", conc_fold[[ak[j]]]),
                              cex=0.72, col=ifelse(conc_fold[[ak[j]]]<=3,"#1a7a3a","grey45"), font=2)
legend("topright", bty="n", cex=0.9, fill=c("#1a7a3a","#b3331f"), border="grey25",
       legend=c(sprintf("standard (%s) = citrate reference",REF_LBL),
                sprintf("tissue (%s) = window intensity",TIS$mode)))
title("Adduct profile: authentic citrate standard (reference) vs tissue window intensities", cex.main=1.06, line=2.6, adj=0)
n_conc <- sum(conc_fold[ak] <= 3)
mtext(sprintf("The standard shows a clean citrate adduct profile: [M-H]- dominant with a minor Na adduct (~%.0f%%). Only Na agrees in tissue (x%.1f); Na2/Cl/acetate/dimer are %d-%dx higher in tissue.",
              100*r_std[["Na"]], conc_fold[["Na"]], round(min(conc_fold[c("Na2","Cl","OAc","dimer")])), round(max(conc_fold[c("Na2","Cl","OAc","dimer")]))),
      side=1, line=5.4, cex=0.68, col="grey25")
mtext("Those tissue windows are dominated by UNRELATED co-isobaric tissue ions (not citrate adducts) -> tissue adduct ratios are NOT a clean citrate readout, consistent with the documented 191 blend. The standard's match to theory (mass + 13C) is the citrate evidence; tissue ID rests on the citrate-mass component + spatial correlation (07_metabolite_id).",
      side=1, line=6.3, cex=0.62, col="grey45", font=3)
dev.off()

# ---- CSV ------------------------------------------------------------------
out <- data.frame(
  ion = c("13C/12C", al),
  key = c("C13", ak),
  ratio_standard = signif(c(iso_std, r_std[ak]), 4),
  ratio_tissue   = signif(c(iso_tis, r_tis[ak]), 4),
  theoretical    = c(signif(ISO_THEO,3), rep(NA, length(ak))),
  fold_std_vs_tis= round(c(conc_fold[["C13"]], conc_fold[ak]), 2),
  within_3x       = c(conc_fold[["C13"]], conc_fold[ak]) <= 3,
  note = c("standard matches theory; tissue inflated by co-isobar",
           rep("std = clean citrate adduct; tissue window has unrelated co-isobaric ions", length(ak))),
  stringsAsFactors = FALSE)
attr(out, "tissue_mode") <- TIS$mode
write.csv(out, OUT_CSV, row.names = FALSE)

log_cfg(TAG, "DONE -> %s ; %s", OUT_PDF, OUT_CSV)
