#!/usr/bin/env Rscript
# ============================================================================
# 04_nearfield_figure.R - publication figure for the near-field citrate
#   emission result: apical-OUT organoids emit more absolute citrate [M-H]- into
#   the immediate gel (0-50/0-100 um) than basolateral-out. Thin DRIVER over the LOCKED
#   visualization toolkit in R/00_lib/lib_nearfield_viz.R (4 views: overlay / gradient
#   map / weather-rainbow heatmap / outward emission vectors). See
#   docs/apical_nearfield.md for the full spec and run order.
#
# In : cache/peaks_tissue_combined.rds, cache/zones_<sid>.rds (R/11_per_organoid_final/01_zones_curated.R),
#      cache/register/nd2final_<sid>.rds, figures/registration/crops/optical_<sid>.png,
#      results/annotation/apical_gradient_per_organoid.csv  (from 03_apical_report.R, CONSENSUS map)
# Out: figures/annotation/apical_nearfield_emission_figure.pdf  (hero/montage/anchor)
#      figures/annotation/nearfield_panels/<sid>_inst<k>_<view>.{pdf,png}
# Usage: Rscript R/11_per_organoid_final/04_nearfield_figure.R
# ============================================================================

ROOT <- file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final")
source(file.path(ROOT, "R/00_lib/gradient_config.R"))
source(file.path(ROOT, "R/00_lib/lib_register.R"))
source(file.path(ROOT, "R/00_lib/lib_nearfield_viz.R"))   # LOCKED 4-view toolkit

# Optional [1] flag drops the ambiguous "mixed" class -> two-class figure; reads the
# matching _nomixed gradient CSV (from 03_apical_report.R) and suffixes all outputs.
args <- commandArgs(trailingOnly = TRUE)
DROP_MIXED <- length(args) >= 1 && nzchar(args[1]) &&
                tolower(args[1]) %in% c("1", "true", "nomixed", "drop_mixed")
SUFFIX <- if (DROP_MIXED) "_nomixed" else ""
NF_CLASSES <- c("basolateral_out", "apical_out", "mixed")
if (DROP_MIXED) NF_CLASSES <- setdiff(NF_CLASSES, "mixed")

ANNOT_FIG <- file.path(FIG_DIR, "annotation")
PANEL_DIR <- file.path(ANNOT_FIG, sprintf("nearfield_panels%s", SUFFIX))
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PDF  <- file.path(ANNOT_FIG, sprintf("apical_nearfield_emission_figure%s.pdf", SUFFIX))
GRAD_CSV <- file.path(RES_DIR, "annotation", sprintf("apical_gradient_per_organoid%s.csv", SUFFIX))

nfviz_load_ion(CIT_MZ)                              # sets val_cit, HI_CIT, pd (citrate [M-H]-)

# ---- pick exemplars (top-3 apical-out / bottom-3 basolateral-out by near50) -------
grad <- read.csv(GRAD_CSV, stringsAsFactors = FALSE)
ao <- grad[grad$apical_class == "apical_out" & is.finite(grad$near50), ]; ao <- ao[order(-ao$near50), ]
ai <- grad[grad$apical_class == "basolateral_out"  & is.finite(grad$near50), ]; ai <- ai[order(ai$near50), ]
AO3 <- head(ao, 3); AI3 <- head(ai, 3); HERO <- ao[1, ]; EX <- rbind(AO3, AI3)
cat("[106] apical-OUT exemplars:\n"); print(AO3[, c("sid","instance","near50","near100")], row.names = FALSE)
cat("[106] basolateral-out exemplars:\n");  print(AI3[, c("sid","instance","near50","near100")], row.names = FALSE)
cliffs <- function(a,b){a<-a[is.finite(a)];b<-b[is.finite(b)];mean(outer(a,b,">"))-mean(outer(a,b,"<"))}
D50  <- cliffs(grad$near50[grad$apical_class=="apical_out"],  grad$near50[grad$apical_class=="basolateral_out"])
D100 <- cliffs(grad$near100[grad$apical_class=="apical_out"], grad$near100[grad$apical_class=="basolateral_out"])

nfviz_arrow_range(unique(EX[, c("sid","instance")]))   # sets GLO, GHI (arrow thickness scale)

# ===========================================================================
# RENDER multi-page PDF
# ===========================================================================
hero_sec <- prep_sec(HERO$sid)
pdf(OUT_PDF, width = 11, height = 8.5)

# ---- Page 1: HERO (all 4 views, 2x2) --------------------------------------
layout(matrix(1:4, nrow = 2, byrow = TRUE)); par(oma = c(0,0,4.2,0))
draw_overlay(hero_sec, HERO$instance, "1) citrate on brightfield (overlay)")
draw_gradmap(hero_sec, HERO$instance, "2) citrate ion image + 0/50/100 um rings")
draw_heatmap(hero_sec, HERO$instance, "3) interpolated citrate heatmap on BF (weather rainbow)")
draw_vectors(hero_sec, HERO$instance, "4) emission vectors (length=absolute, thickness=relative)")
mtext("Near-field citrate emission - hero apical-out organoid", outer = TRUE, font = 2, cex = 1.1, line = 1.9)
mtext(sprintf("%s inst %d   near50 = %.2f   (apical-out vs basolateral-out, 0-50 um: Cliff's d = %+.2f, ***)",
              sub("AO_","",HERO$sid), HERO$instance, HERO$near50, D50), outer = TRUE, cex = 0.82, line = 0.6)

# ---- Pages 2-3: contrast montage (3 apical-out, then 3 basolateral-out) ----------
draw_set <- function(rows, tag) for (r in seq_len(nrow(rows))) {
  sec <- prep_sec(rows$sid[r]); k <- rows$instance[r]
  lab <- sprintf("%s %s i%d  near50=%.2f", tag, sub("AO_","",rows$sid[r]), k, rows$near50[r])
  draw_overlay(sec, k, lab); draw_gradmap(sec, k, "ion image + rings"); draw_vectors(sec, k, "emission vectors")
}
layout(matrix(1:9, nrow = 3, byrow = TRUE)); par(oma = c(0,0,3,0)); draw_set(AO3, "OUT")
mtext("Apical-OUT (high near-field emission): overlay | ion image+rings | emission vectors",
      outer = TRUE, font = 2, cex = 0.92, line = 1.0)
layout(matrix(1:9, nrow = 3, byrow = TRUE)); par(oma = c(0,0,3,0)); draw_set(AI3, "IN")
mtext("basolateral-out (low near-field emission): overlay | ion image+rings | emission vectors",
      outer = TRUE, font = 2, cex = 0.92, line = 1.0)

# ---- Page 4: quantitative anchor (dot+box + radial profiles) ---------------
layout(matrix(c(1,2,3,3), nrow = 2, byrow = TRUE)); par(oma = c(0,0,3,0), mar = c(4,4.5,3,1))
dotbox_simple <- function(metric, ylab, main) {
  cls <- factor(grad$apical_class, levels = NF_CLASSES)
  yv <- grad[[metric]]; ok <- is.finite(yv)
  plot(NA, xlim = c(0.5, length(NF_CLASSES)+0.5), ylim = range(yv[ok])*c(0.95,1.08), xaxt="n", xlab="", ylab=ylab, main=main)
  cols <- APICAL_COLS   # locked: green basolateral-out / magenta apical-out / grey mixed
  for (i in seq_along(NF_CLASSES)) { lv <- levels(cls)[i]; v <- yv[ok & cls==lv]; if(!length(v)) next
    q <- quantile(v,c(.25,.5,.75)); rect(i-.28,q[1],i+.28,q[3],border="grey35"); segments(i-.28,q[2],i+.28,q[2],lwd=2.4)
    set.seed(i); points(i+runif(length(v),-.15,.15), v, pch=21, bg=adjustcolor(cols[lv],.7), col="white", cex=1) }
  axis(1, at=seq_along(NF_CLASSES), labels=APICAL_LABS[NF_CLASSES], cex.axis=0.9)
}
dotbox_simple("near50", "citrate (TIC-norm)", sprintf("0-50 um  (d=%+.2f)", D50))
dotbox_simple("near100","citrate (TIC-norm)", sprintf("0-100 um  (d=%+.2f)", D100))
rp <- radial_df(rbind(ao[,c("sid","instance","apical_class")], ai[,c("sid","instance","apical_class")]))
med <- lapply(c("basolateral_out","apical_out"), function(clz) tapply(rp$cit[rp$apical_class==clz], rp$dist[rp$apical_class==clz], median, na.rm=TRUE))
ymax <- max(unlist(med), na.rm = TRUE) * 1.08
par(mar = c(4,4.5,3,1))
plot(NA, xlim = c(-40,160), ylim = c(0, ymax),
     xlab = "signed distance from surface (um;  <0 inside organoid, >0 gel)", ylab = "citrate (TIC-norm)",
     main = "Radial citrate profile: basolateral-out vs apical-out (group median)")
abline(v = 0, lty = 3, col = "grey50"); rect(0,par("usr")[3],50,par("usr")[4],col=adjustcolor("grey80",0.25),border=NA)
for (clz in c("basolateral_out","apical_out")) { s <- rp[rp$apical_class==clz,]
  agg <- tapply(s$cit, s$dist, median, na.rm=TRUE)
  lines(as.numeric(names(agg)), as.numeric(agg), col = if(clz=="basolateral_out") CIN else COUT, lwd = 3, type="b", pch=19) }
legend("topright", legend=c("basolateral-out","apical-out","0-50 um window"), col=c(CIN,COUT,adjustcolor("grey70",0.6)),
       lwd=c(3,3,8), pch=c(19,19,NA), bty="n")
mtext("Quantitative anchor: per-organoid near-field emission + radial decay (organoid = unit, descriptive)",
      outer = TRUE, font = 2, cex = 1.05, line = 1.0)
dev.off()
cat(sprintf("[106] DONE -> %s\n", OUT_PDF))

# ===========================================================================
# standalone high-res panels (hero 4 views + each montage cell)
# ===========================================================================
VIEWS <- list(overlay = draw_overlay, gradmap = draw_gradmap, heatmap = draw_heatmap, vectors = draw_vectors)
export_panel <- function(sid, k, view, fn) {
  sec <- prep_sec(sid); base <- file.path(PANEL_DIR, sprintf("%s_inst%d_%s", sub("AO_","",sid), k, view))
  png(paste0(base, ".png"), width = 1500, height = 1500, res = 300); par(mar=c(2,1,2,1)); fn(sec, k, ""); dev.off()
  pdf(paste0(base, ".pdf"), width = 4, height = 4); par(mar=c(2,1,2,1)); fn(sec, k, ""); dev.off()
}
for (v in names(VIEWS)) export_panel(HERO$sid, HERO$instance, v, VIEWS[[v]])
for (r in seq_len(nrow(EX))) for (v in c("overlay","gradmap","vectors"))
  export_panel(EX$sid[r], EX$instance[r], v, VIEWS[[v]])
cat(sprintf("[106] panels -> %s\n", PANEL_DIR))
