#!/usr/bin/env Rscript
# ============================================================================
# 09_citrate_gradient_perdataset_v3.R - v3 of the per-dataset combined report.
#   *** LOCKED *** - frozen visual spec; see docs/citrate_gradient_perdataset.md.
#   Frozen visual spec (developed from an earlier single-TEST-section prototype
#   that has since been dropped in this fork), but:
#     - DEFAULT scope is ALL 20 datasets.
#     - Output is citrate_gradient_perdataset_v3.pdf.
#   Built on the CURRENT consensus-curated segmentation: each cache/zones_<sid>.rds
#   is rebuilt by R/11_per_organoid_final/01_zones_curated.R from the curated
#   instances_{final,clean,split}_<sid>.rds, whose instance IDs match the apical
#   consensus annotations. No class recolour - organoid outlines stay white-dotted
#   and per-organoid gradient curves keep the rainbow palette.
#
#   Page layout:
#     (a) thin WHITE DOTTED organoid outline (segmentation, signed-dist=0) on the
#         Citrate [M-H]- and DHA [M-H]- single-ion panels of each per-dataset page.
#     (b) a SUMMARY PAGE: all 20 datasets tiled as the citrate(yellow)+DHA(purple)
#         overlay, PLUS the citrate-vs-DHA outward-decay curve comparison
#         (per-organoid faint + bold median; citrate=yellow, DHA=purple).
#
# Usage:
#   Rscript R/11_per_organoid_final/09_citrate_gradient_perdataset_v3.R        # all 20 pages + summary
#   Rscript R/11_per_organoid_final/09_citrate_gradient_perdataset_v3.R all    # all 20 pages + summary
#   Rscript R/11_per_organoid_final/09_citrate_gradient_perdataset_v3.R <sid>  # one named section + summary
#
# In : cache/peaks_tissue_combined.rds, cache/zones_<sid>.rds (R/11_per_organoid_final/01_zones_curated.R),
#      cache/register/nd2final_<sid>.rds, figures/registration/crops/optical_<sid>.png
# Out: figures/gradient/citrate_gradient_perdataset_v3[_TEST_<sid>].pdf
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
suppressPackageStartupMessages({ library(Cardinal); library(png); library(viridisLite) })

REG_CACHE <- cache_in("register")
CROP_DIR  <- file.path(FIG_DIR, "registration", "crops")
SCALE_UM  <- 100; PX_UM <- MSI_PIXEL_UM
CIT_MZ <- 191.0217; DHA_MZ <- 327.2330
CIT_RGB <- c(1, 1, 0); DHA_RGB <- c(0.7, 0, 1)   # yellow ; purple
CIT_LINE <- "#d4b800"; DHA_LINE <- "#8e00ff"     # readable line colours (yellow/purple) for curves
OVERLAY_ALPHA <- 0.65
GRID_D <- 2L
MIN_ZONE_PX <- 3L
pal  <- viridis(256); ramp <- viridis(64)
ORG_PAL <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628","#f781bf",
             "#1b9e77","#d95f02","#7570b3","#66a61e","#e7298a")
org_colors <- function(ids) setNames(ORG_PAL[((seq_along(ids)-1) %% length(ORG_PAL))+1], as.character(ids))
log_msg <- function(...) message(sprintf("[v3] %s", sprintf(...)))

# ---- backbone MSE ----------------------------------------------------------
mse <- readRDS(TISSUE_MSE)
pd  <- as.data.frame(pixelData(mse)); pd$gidx <- seq_len(nrow(pd))
mzs <- mz(mse)
X   <- as.matrix(spectra(mse))
cit <- which.min(abs(mzs - CIT_MZ)); dha <- which.min(abs(mzs - DHA_MZ))
log_msg("citrate feat %d (%.4f) | DHA feat %d (%.4f)", cit, mzs[cit], dha, mzs[dha])
val_dha <- X[dha, ]
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_citrate.R"))
val_cit <- citrate_onto_pd(pd)   # anchored citrate (raw imzML, TIC-norm); cit kept for the log only
HI     <- as.numeric(quantile(val_cit[val_cit > 0], IMG_CLIP_HI, na.rm = TRUE))
DHA_HI <- as.numeric(quantile(val_dha[val_dha > 0], IMG_CLIP_HI, na.rm = TRUE))
if (!is.finite(HI)     || HI     <= 0) HI     <- max(val_cit, na.rm = TRUE)
if (!is.finite(DHA_HI) || DHA_HI <= 0) DHA_HI <- max(val_dha, na.rm = TRUE)
log_msg("global p99.5 clips: citrate %.3g | DHA %.3g", HI, DHA_HI)

SIDS20  <- levels(pixelData(mse)$sample_id); if (is.null(SIDS20)) SIDS20 <- sort(unique(as.character(pd$sample_id)))
FC      <- setNames(rep("#444444", length(SIDS20)), SIDS20)   # neutral per-dataset frame colour
ord20   <- sort(SIDS20)

a <- commandArgs(trailingOnly = TRUE)
# v3: default (no args) -> ALL 20 datasets.  `all` -> all 20.  `<sid>` -> one named section.
RENDER_SIDS <- if (length(a) == 0 || tolower(a[1]) == "all") ord20 else a[1]
RENDER_SIDS <- RENDER_SIDS[RENDER_SIDS %in% SIDS20]
stopifnot(length(RENDER_SIDS) >= 1)
OUT_PDF <- if (length(RENDER_SIDS) == 1) {
  file.path(GRAD_FIG, sprintf("citrate_gradient_perdataset_v3_TEST_%s.pdf", RENDER_SIDS))
} else file.path(GRAD_FIG, "citrate_gradient_perdataset_v3.pdf")

# ============================================================================
# RENDERERS
# ============================================================================
sec_data <- function(sid) {
  m <- which(as.character(pd$sample_id) == sid); x <- pd$x[m]; y <- pd$y[m]
  xs <- sort(unique(x)); ys <- sort(unique(y)); ix <- match(x, xs); iy <- match(y, ys)
  fill <- function(v){ mm <- matrix(0, length(xs), length(ys)); mm[cbind(ix, iy)] <- v; mm }
  list(xs = xs, ys = ys, cit = fill(val_cit[m]), dha = fill(val_dha[m]))
}
# signed-distance matrix aligned to a section's (xs,ys) grid -> for the dotted outline
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
# optional `outline` (signed-dist matrix) -> thin white dotted organoid boundary at level 0
ion_panel <- function(mat, xs, ys, hi, main, fc, mar, colorbar, outline = NULL) {
  sc <- pmin(mat / hi, 1); par(mar = mar)
  image(xs, ys, sc, col = pal, asp = 1, useRaster = TRUE, zlim = c(0,1),
        axes = FALSE, xlab = "", ylab = "", main = main, cex.main = 0.95, col.main = fc)
  box(col = fc, lwd = 1.6)
  if (!is.null(outline) && any(is.finite(outline)))
    contour(xs, ys, outline, levels = 0, add = TRUE, drawlabels = FALSE, col = "white", lty = 3, lwd = 0.8)
  if (colorbar) {
    cx0<-grconvertX(1.02,"npc","user"); cx1<-grconvertX(1.06,"npc","user")
    cy0<-grconvertY(0.05,"npc","user"); cy1<-grconvertY(0.95,"npc","user")
    yb<-seq(cy0,cy1,length.out=257); rect(cx0,head(yb,-1),cx1,tail(yb,-1),col=pal,border=NA,xpd=NA)
    text(cx1,grconvertY(c(0.05,0.95),"npc","user"),sprintf("%.2g",c(0,hi)),pos=4,cex=0.6,xpd=NA,offset=0.15)
  }
  scale_bar_msi(xs, ys)
}
overlay_panel <- function(cit, dha, xs, ys, main, fc, mar, scalebar = TRUE) {
  cn <- pmin(cit/HI, 1); dn <- pmin(dha/DHA_HI, 1)
  R <- pmin(cn*CIT_RGB[1] + dn*DHA_RGB[1], 1); G <- pmin(cn*CIT_RGB[2] + dn*DHA_RGB[2], 1); B <- pmin(cn*CIT_RGB[3] + dn*DHA_RGB[3], 1)
  fl <- function(ch) t(ch[, rev(seq_len(ncol(ch)))])
  ras <- as.raster(array(c(fl(R), fl(G), fl(B)), dim = c(length(ys), length(xs), 3)))
  par(mar = mar)
  plot(range(xs), range(ys), type="n", asp=1, axes=FALSE, xlab="", ylab="", main = main, cex.main = 0.9, col.main = fc)
  rasterImage(ras, min(xs), min(ys), max(xs), max(ys), interpolate = FALSE)
  box(col = fc, lwd = 1.4); if (scalebar) scale_bar_msi(xs, ys)
}
inv_affine <- function(B,NX,NY){ M<-rbind(B[1,],B[2,]); t0<-B[3,]; xy<-solve(M, rbind(NX-t0[1],NY-t0[2])); list(x=xy[1,],y=xy[2,]) }

bf_panels <- function(sid, z, org_col) {
  xf <- file.path(REG_CACHE, sprintf("nd2final_%s.rds", sid)); pf <- file.path(CROP_DIR, sprintf("optical_%s.png", sid))
  draw_missing <- function(t){ par(mar=c(1,1,3,1)); plot.new(); title(t, cex.main=0.95); text(0.5,0.5,"(no native crop)",col="grey50") }
  if (!file.exists(xf) || !file.exists(pf)) { draw_missing(paste(sid,"native BF")); draw_missing("+ citrate overlay"); return(invisible()) }
  Xr <- readRDS(xf); B <- Xr$B_msi_nd2; smn <- Xr$scale_msi_nd2; cx0 <- Xr$crop[1]; cy0 <- Xr$crop[2]
  om <- png::readPNG(pf); if (length(dim(om))==3) om <- om[,,1]; cw <- ncol(om); ch <- nrow(om)
  sub <- pd[as.character(pd$sample_id)==sid, c("x","y","is_tissue")]; W<-max(sub$x); H<-max(sub$y)
  secval <- val_cit[as.character(pd$sample_id)==sid]
  gx<-seq(1,cw,by=GRID_D); gy<-seq(1,ch,by=GRID_D)
  NX<-cx0+rep(gx,times=length(gy))-1; NY<-cy0+rep(gy,each=length(gx))-1
  ms<-inv_affine(B,NX,NY); inb<-ms$x>=0.5&ms$x<=W+0.5&ms$y>=0.5&ms$y<=H+0.5
  lut<-matrix(NA_real_,W,H); lut[cbind(sub$x,sub$y)]<-secval
  val<-rep(NA_real_,length(NX)); if(any(inb)){ xi<-pmin(pmax(round(ms$x[inb]),1),W); yi<-pmin(pmax(round(ms$y[inb]),1),H); val[inb]<-lut[cbind(xi,yi)] }
  vm<-matrix(val,length(gy),length(gx),byrow=TRUE)
  hi<-quantile(vm[is.finite(vm)&vm>0],IMG_CLIP_HI,na.rm=TRUE); if(!is.finite(hi)||hi<=0)hi<-1; vs<-pmin(vm/hi,1)
  drawbf<-function(t){ par(mar=c(1,1,3,1)); plot.new(); plot.window(c(1,cw),c(ch,1),asp=1); rasterImage(om/max(om),1,ch,cw,1); title(t,cex.main=0.95) }
  drawbf("native BF"); umpx<-PX_UM/smn; bp<-SCALE_UM/umpx
  segments(cw-bp-10,ch-14,cw-10,ch-14,lwd=3,col="black"); text(cw-bp/2-10,ch-28,sprintf("%d um",SCALE_UM),cex=0.7)
  drawbf(sprintf("+ citrate [M-H]- (%.0f%% opacity) + zones", OVERLAY_ALPHA*100))
  fin <- is.finite(vs)
  if (any(fin)) { cidx <- pmin(pmax(round(vs*63)+1,1),64); cmat <- matrix("#00000000", length(gy), length(gx))
    cmat[fin] <- adjustcolor(ramp[cidx[fin]], alpha.f = OVERLAY_ALPHA); rasterImage(as.raster(cmat), 1, ch, cw, 1, interpolate=FALSE) }
  if (!is.null(z)) {
    lutd <- matrix(NA_real_, W, H); lutd[cbind(z$x, z$y)] <- z$signed_dist_um
    dnat <- rep(NA_real_, length(NX)); if (any(inb)) { xi<-pmin(pmax(round(ms$x[inb]),1),W); yi<-pmin(pmax(round(ms$y[inb]),1),H); dnat[inb]<-lutd[cbind(xi,yi)] }
    dmat <- matrix(dnat, length(gy), length(gx), byrow=TRUE)
    if (any(is.finite(dmat)))
      contour(gx, gy, t(dmat), levels=OUT_ZONE_UM, add=TRUE, drawlabels=FALSE, col=adjustcolor("white",0.7), lwd=0.6)
    for (k in names(org_col)) {
      sf <- z[z$is_surface & z$instance == as.integer(k), c("x","y")]
      if (nrow(sf)) { tn <- apply_affine(B, as.matrix(sf)); points(tn[,1]-cx0+1, tn[,2]-cy0+1, pch=15, cex=0.5, col=org_col[k]) }
    }
  }
  invisible()
}

spr <- function(v) { idx <- which(!is.na(v)); if (length(idx) < 3 || sd(v[idx]) == 0) return(NA_real_)
  suppressWarnings(cor(idx, v[idx], method = "spearman")) }

gradient_panel <- function(sid, z, org_col, fc) {
  par(mar = c(4.2, 4.2, 3, 1))
  if (is.null(z)) { plot.new(); title("outward gradient", cex.main=0.95); text(0.5,0.5,"(no zones cache)",col="grey50"); return(invisible()) }
  nz <- length(OUT_ZONE_UM); zk <- seq_len(nz)
  zmean_k <- function(v, k) vapply(zk, function(zz){ g <- z$gidx[which(z$instance_catch==k & z$zone_out==zz)]; if(!length(g)) NA_real_ else mean(v[g]) }, numeric(1))
  zmean_all <- function(v) vapply(zk, function(zz){ g <- z$gidx[which(z$zone_out==zz)]; if(!length(g)) NA_real_ else mean(v[g]) }, numeric(1))
  inside_k  <- function(v, k) { g <- z$gidx[which(z$instance == k)];   if(!length(g)) NA_real_ else mean(v[g]) }
  inside_all<- function(v)    { g <- z$gidx[which(z$instance > 0)];    if(!length(g)) NA_real_ else mean(v[g]) }
  ids <- as.integer(names(org_col))
  citM <- sapply(ids, function(k) zmean_k(val_cit, k)); if (is.null(dim(citM))) citM <- matrix(citM, nrow=nz)
  cit_in  <- vapply(ids, function(k) inside_k(val_cit, k), numeric(1))
  rho_org <- apply(citM, 2, spr)
  md <- zmean_all(val_dha); dha_in <- inside_all(val_dha)
  cmax <- max(c(citM, cit_in), na.rm=TRUE); if(!is.finite(cmax)||cmax<=0) cmax <- 1
  dmax <- max(c(md, dha_in),   na.rm=TRUE); if(!is.finite(dmax)||dmax<=0) dmax <- 1
  xp <- 0:nz
  plot(NA, xlim=c(0,nz), ylim=c(0,1.05), xaxt="n", xlab="distance from organoid surface (um)",
       ylab="relative mean intensity (per ion)", main="outward gradient per organoid", cex.main=0.95, col.main=fc)
  axis(1, at=xp, labels=c("in", OUT_ZONE_UM), cex.axis=0.75)
  rect(par("usr")[1], par("usr")[3], 0.5, par("usr")[4], col="#00000010", border=NA)
  abline(v=0.5, lty=3, col="grey60")
  lines(xp, c(dha_in, md)/dmax, col="grey50", lwd=2, lty=2, type="b", pch=17)
  for (j in seq_along(ids)) lines(xp, c(cit_in[j], citM[,j])/cmax, col=org_col[as.character(ids[j])], lwd=2.4, type="b", pch=19)
  leg <- c(sprintf("organoid %d (rho_out=%s)", ids, ifelse(is.na(rho_org),"NA",sprintf("%+.2f",rho_org))), "DHA C22:6 control")
  legend("topright", bty="n", cex=0.78, legend=leg, col=c(unname(org_col[as.character(ids)]), "grey50"),
         lwd=c(rep(2.4,length(ids)),2), lty=c(rep(1,length(ids)),2), pch=c(rep(19,length(ids)),17))
  invisible()
}

# per-organoid outward profile normalized to its 0-10um surface band
collect_profiles <- function() {
  nz <- length(OUT_ZONE_UM); rows <- list()
  for (sid in ord20) {
    fz <- cache_in(sprintf("zones_%s.rds", sid)); if (!file.exists(fz)) next
    z <- readRDS(fz)
    for (k in sort(unique(z$instance[z$instance > 0]))) {
      surf_g <- z$gidx[z$instance == k & !is.na(z$zone_in) & z$zone_in == 1]
      rc <- mean(val_cit[surf_g]); rd <- mean(val_dha[surf_g])
      mc <- md <- rep(NA_real_, nz)
      for (zz in seq_len(nz)) { gg <- z$gidx[z$instance_catch == k & !is.na(z$zone_out) & z$zone_out == zz]
        if (length(gg) >= MIN_ZONE_PX) { mc[zz] <- mean(val_cit[gg]); md[zz] <- mean(val_dha[gg]) } }
      rows[[length(rows)+1]] <- data.frame(sample_id = sid, key = paste(sid,k), zone_um = OUT_ZONE_UM,
        cit = if (is.finite(rc) && rc>0) mc/rc else mc*NA, dha = if (is.finite(rd) && rd>0) md/rd else md*NA,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

# small per-dataset decay plot (citrate=yellow vs DHA=purple), drawn beside its overlay
decay_mini <- function(sid, prof) {
  par(mar = c(2.1, 2.2, 1.6, 0.6))
  s <- prof[prof$sample_id == sid, ]
  if (!nrow(s) || all(is.na(c(s$cit, s$dha)))) { plot.new(); title("decay (no data)", cex.main = 0.7); return(invisible()) }
  ymax <- max(1.2, quantile(c(s$cit, s$dha), 0.97, na.rm = TRUE), na.rm = TRUE)
  plot(NA, xlim = range(OUT_ZONE_UM), ylim = c(0, ymax), log = "x", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", main = sprintf("decay %s", sub("AO_","",sid)), cex.main = 0.72)
  axis(1, at = c(10,100,500), labels = c("10","100","500"), cex.axis = 0.5, tcl = -0.2, mgp = c(0,0.05,0))
  axis(2, at = c(0,1), labels = c("0","1"), cex.axis = 0.5, las = 1, tcl = -0.2, mgp = c(0,0.3,0))
  abline(h = 1, lty = 3, col = "grey75")
  for (kk in unique(s$key)) { ss <- s[s$key == kk, ]; ss <- ss[order(ss$zone_um), ]
    lines(ss$zone_um, ss$cit, col = adjustcolor(CIT_LINE, 0.30), lwd = 0.6)
    lines(ss$zone_um, ss$dha, col = adjustcolor(DHA_LINE, 0.30), lwd = 0.6) }
  mc <- tapply(s$cit, s$zone_um, median, na.rm = TRUE); md <- tapply(s$dha, s$zone_um, median, na.rm = TRUE)
  lines(as.numeric(names(mc)), as.numeric(mc), col = CIT_LINE, lwd = 2.2, type = "b", pch = 19, cex = 0.5)
  lines(as.numeric(names(md)), as.numeric(md), col = DHA_LINE, lwd = 2.2, type = "b", pch = 17, cex = 0.5)
}

# ============================================================================
# REPORT
# ============================================================================
dir.create(GRAD_FIG, showWarnings = FALSE, recursive = TRUE)
pdf(OUT_PDF, width = 16, height = 10)
TRI_MAR <- c(3.0, 1.0, 2.4, 3.6)

# ---- PAGE 1 = SUMMARY: 20 (overlay + its own decay) PAIRS tiled across the page
prof <- collect_profiles()
layout(matrix(1:40, nrow = 5, byrow = TRUE))      # 4 pairs/row x 5 rows = 8 cols x 5 rows
par(oma = c(0.5, 1, 4.5, 1))
OV_MAR <- c(0.5, 0.4, 1.5, 0.2)
for (sid in ord20) {
  sd <- sec_data(sid)
  overlay_panel(sd$cit, sd$dha, sd$xs, sd$ys, disp_id(sid), FC[sid], OV_MAR, scalebar = FALSE)
  decay_mini(sid, prof)
}
mtext("SUMMARY - all 20 datasets: each overlay (citrate=yellow, DHA=purple) paired with its outward decay",
      outer = TRUE, line = 2.7, cex = 1.3, font = 2)
mtext(sprintf("decay = per-organoid outward profile normalized to organoid surface (1.0); CITRATE=yellow, DHA=purple; bold=section median, faint=organoids; x=10-500 um log; n=%d organoids",
              length(unique(prof$key))), outer = TRUE, line = 1.1, cex = 0.9, col = "grey25")

# ---- per-dataset pages ------------------------------------------------------
for (sid in RENDER_SIDS) {
  sd <- sec_data(sid); fc <- FC[sid]
  fz <- cache_in(sprintf("zones_%s.rds", sid))
  z <- if (file.exists(fz)) readRDS(fz) else NULL
  inst_ids <- if (is.null(z)) integer(0) else sort(unique(z$instance[z$instance > 0]))
  org_col  <- org_colors(inst_ids)
  sdm <- sdist_mat(z, sd$xs, sd$ys)                                   # for the dotted outline
  layout(rbind(c(1,2,3), c(4,5,6))); par(oma = c(2.2, 1, 5, 1))
  overlay_panel(sd$cit, sd$dha, sd$xs, sd$ys, "Overlay (citrate=yellow, DHA=purple)", fc, TRI_MAR)
  ion_panel(sd$cit, sd$xs, sd$ys, HI,     "Citrate [M-H]- 191.02  (white dotted = organoid outline)",  fc, TRI_MAR, TRUE, outline = sdm)
  ion_panel(sd$dha, sd$xs, sd$ys, DHA_HI, "DHA C22:6 [M-H]- 327.23  (white dotted = organoid outline)", fc, TRI_MAR, TRUE, outline = sdm)
  bf_panels(sid, z, org_col)
  gradient_panel(sid, z, org_col, fc)
  mtext(sprintf("%s   -   %d organoid%s", disp_id(sid), length(inst_ids), if(length(inst_ids)==1)"" else "s"),
        outer=TRUE, line=3.1, cex=1.6, font=2, col=fc)
  mtext(sprintf("global p99.5 clips:  citrate %.3g  |  DHA %.3g        BF overlay = citrate, viridis, %.0f%% constant opacity; white contours = outward zones (%s um)",
                HI, DHA_HI, OVERLAY_ALPHA*100, paste(OUT_ZONE_UM, collapse="/")), outer=TRUE, line=1.5, cex=0.86)
  mtext("white dotted = organoid segmentation outline (consensus-curated)        'in' = mean citrate within organoid body (0 um source); gel pixels grouped by nearest-organoid catchment; line colour = organoid outline colour",
        outer=TRUE, line=0.3, cex=0.70, col="grey30")
}

dev.off()
log_msg("DONE -> %s (%d dataset page%s + 1 summary)", OUT_PDF, length(RENDER_SIDS), if (length(RENDER_SIDS)==1) "" else "s")
