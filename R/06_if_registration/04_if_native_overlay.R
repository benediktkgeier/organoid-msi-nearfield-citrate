#!/usr/bin/env Rscript
# ============================================================================
# 04_if_native_overlay.R - STEP 3: native-resolution IF overlay report. Per
#   section, 4 panels:
#     A  hi-res DAPI located in the overview            (validates B_hr_ov, R/06_if_registration/02_if_overview_to_bf.R)
#     B  native BF crop + MSI is_tissue (green) + IF DAPI footprint (red)
#                                                       (validates B_hr_bf, R/06_if_registration/03_if_overview_to_msi.R)
#     C  native BF crop + 4-channel IF composite (homogeneous opacity)
#        DAPI=cyan, b-catenin(GFP)=green, ZO-1(mCherry)=red, F-actin(Cy5)=magenta
#     D  MSI citrate [M-H]- 191.0217 (viridis, 65% opacity, locked style) with a
#        combined DAPI(cyan)+ZO-1(red) IF overlay in the same frame  (user request)
#
#   Histograms shown as-acquired (block-mean, min-max scaled) - no extra tuning.
#
# Input : cache/register_if/hr_to_bf_<sid_if>.rds, locate_<sid_if>.rds, ovthumb,
#         cache/register/nd2final_*.rds, TISSUE_MSE, hi-res + BF .nd2
# Output: figures/if_registration/step3_native_overlay.pdf, crops/if_optical_<sid_if>.png
# Usage : Rscript R/06_if_registration/04_if_native_overlay.R [all | sid_if ...]   (default = pilot)
# ============================================================================

Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jre1.8.0_491")
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/if_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register_if.R"))
suppressPackageStartupMessages({ library(Cardinal); library(RBioFormats); library(png); library(viridisLite) })

ION_MZ <- 191.0217; OVERLAY_ALPHA <- 0.65; GRID_D <- 2L
args <- commandArgs(trailingOnly = TRUE)
SEC  <- if_sections()
sel  <- if (length(args) == 0) PILOT_SIDS else if (identical(args, "all")) SEC$sid_if else args
SEC  <- SEC[SEC$sid_if %in% sel, , drop = FALSE]

mse <- readRDS(TISSUE_MSE); pd <- as.data.frame(pixelData(mse)); mzs <- mz(mse)
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_citrate.R"))
val_ion <- citrate_onto_pd(pd)   # anchored citrate (raw imzML, TIC-norm); was grid feat ION_MZ
inv_affine_pts <- function(B, q) sweep(q, 2, B[3, ]) %*% solve(B[1:2, ])
clip01 <- function(m, p = 0.995, g = 1) { hi <- quantile(m[is.finite(m) & m > 0], p, na.rm = TRUE); if (!is.finite(hi) || hi <= 0) hi <- 1; pmin(m/hi, 1)^g }
IF_GAMMA <- 0.55   # display brightening for faint IF channels (data unchanged)

pdf(file.path(IF_FIG, "step3_native_overlay.pdf"), width = 15, height = 4.6)
for (i in seq_len(nrow(SEC))) {
  r <- SEC[i, ]; sid <- r$sid_if
  hf <- file.path(IF_CACHE, sprintf("hr_to_bf_%s.rds", sid)); lf <- file.path(IF_CACHE, sprintf("locate_%s.rds", sid))
  if (!file.exists(hf)) { cat(sprintf("[93] no transform for %s\n", sid)); next }
  X <- readRDS(hf); L <- readRDS(lf); B <- X$B_hr_bf; SXh <- X$SX_hr; SYh <- X$SY_hr
  ov <- readRDS(file.path(IF_CACHE, sprintf("ovthumb_%s.rds", r$slide)))

  # BF crop bbox = MSI organoid footprint bbox (zoom to organoid), 1.6x margin
  msi_sids <- msi_sids_for_section(r$msi_slide, r$secn)
  foot <- do.call(rbind, lapply(msi_sids, function(s){ Bm<-readRDS(file.path(REG_CACHE,sprintf("nd2final_%s.rds",s)))$B_msi_nd2; sub<-pd[as.character(pd$sample_id)==s & pd$is_tissue,c("x","y")]; apply_affine(Bm,cbind(sub$x,sub$y)) }))
  ext <- apply(foot, 2, range); ctr <- colMeans(ext); half <- max(diff(ext[,1]), diff(ext[,2]))/2 * 1.6 + 60
  bfdim <- coreMetadata(read.metadata(r$msi_bf), series=1)
  cx0 <- max(1, floor(ctr[1]-half)); cy0 <- max(1, floor(ctr[2]-half))
  cx1 <- min(bfdim$sizeX, ceiling(ctr[1]+half)); cy1 <- min(bfdim$sizeY, ceiling(ctr[2]+half))
  bf <- read.image(r$msi_bf, series=1, subset=list(x=cx0:cx1, y=cy0:cy1), normalize=TRUE)
  a <- as.array(bf); bfm <- if (length(dim(a))==2) a else a[,,1]; bfm <- t(bfm)   # [row=y,col=x]
  cw <- ncol(bfm); ch <- nrow(bfm)

  # 4-channel hi-res thumbs at F=5 (~1.81 um/px), single-pass block-mean
  ALL <- nd2_block_mean_allch(r$hr_path, F=F_HR); CH <- lapply(1:dim(ALL)[3], function(k) ALL[,,k])
  m5dim <- dim(CH[[1]])
  # map a BF-crop grid -> hi-res m5 px, sample channels
  gx <- seq(1, cw, by=GRID_D); gy <- seq(1, ch, by=GRID_D)
  NX <- cx0 + rep(gx, times=length(gy)) - 1; NY <- cy0 + rep(gy, each=length(gx)) - 1
  hn <- inv_affine_pts(B, cbind(NX, NY)); m5x <- hn[,1]/F_HR + 0.5; m5y <- hn[,2]/F_HR + 0.5
  inb <- m5x>=1 & m5x<=m5dim[2] & m5y>=1 & m5y<=m5dim[1]
  samp <- function(M) { v <- rep(NA_real_, length(NX)); if (any(inb)) v[inb] <- M[cbind(pmin(pmax(round(m5y[inb]),1),m5dim[1]), pmin(pmax(round(m5x[inb]),1),m5dim[2]))]; matrix(v, length(gy), length(gx), byrow=TRUE) }
  DAPI <- clip01(samp(CH[[1]]), g=IF_GAMMA); GFP <- clip01(samp(CH[[2]]), g=IF_GAMMA); ZO1 <- clip01(samp(CH[[3]]), g=IF_GAMMA); FAC <- clip01(samp(CH[[4]]), g=IF_GAMMA)

  # MSI citrate projected into BF crop (per overlapping MSI dataset)
  ionm <- matrix(NA_real_, length(gy), length(gx))
  for (sid_msi in msi_sids_for_section(r$msi_slide, r$secn)) {
    Bm <- readRDS(file.path(REG_CACHE, sprintf("nd2final_%s.rds", sid_msi)))$B_msi_nd2
    sub <- pd[as.character(pd$sample_id)==sid_msi, c("x","y")]; W<-max(sub$x); H<-max(sub$y)
    sv <- val_ion[as.character(pd$sample_id)==sid_msi]
    lut <- matrix(NA_real_, W, H); lut[cbind(sub$x, sub$y)] <- sv
    ms <- inv_affine_pts(Bm, cbind(NX, NY)); ok <- ms[,1]>=0.5 & ms[,1]<=W+0.5 & ms[,2]>=0.5 & ms[,2]<=H+0.5
    vv <- rep(NA_real_, length(NX)); if (any(ok)) vv[ok] <- lut[cbind(pmin(pmax(round(ms[ok,1]),1),W), pmin(pmax(round(ms[ok,2]),1),H))]
    vmm <- matrix(vv, length(gy), length(gx), byrow=TRUE); ionm <- pmax(ionm, vmm, na.rm=TRUE)
  }
  ions <- clip01(ionm)

  # MSI is_tissue + IF DAPI footprint outlines in crop px
  tis <- pd[as.character(pd$sample_id) %in% msi_sids_for_section(r$msi_slide, r$secn) & pd$is_tissue, ]
  tnf <- NULL
  if (nrow(tis)) { for (sid_msi in unique(as.character(tis$sample_id))) { Bm <- readRDS(file.path(REG_CACHE, sprintf("nd2final_%s.rds", sid_msi)))$B_msi_nd2; tt <- tis[as.character(tis$sample_id)==sid_msi,]; nd <- apply_affine(Bm, cbind(tt$x,tt$y)); tnf <- rbind(tnf, nd) } }
  fgf <- which(DAPI > 0.45, arr.ind=TRUE)   # IF DAPI footprint in crop grid

  # ---------- draw ----------
  layout(matrix(1:4, nrow=1)); par(mar=c(1,1,3,1))
  # A: overview context
  b <- L$box; ox0<-max(1,b["c0"]-30); oy0<-max(1,b["r0"]-30); ox1<-min(ncol(ov$m16),b["c0"]+b["tw"]+30); oy1<-min(nrow(ov$m16),b["r0"]+b["th"]+30)
  sub16 <- ov$m16[oy0:oy1, ox0:ox1]
  plot.new(); plot.window(c(1,ncol(sub16)), c(nrow(sub16),1), asp=1); rasterImage(as.raster(norm01(sub16)),1,nrow(sub16),ncol(sub16),1)
  rect(b["c0"]-ox0+1, b["r0"]-oy0+1, b["c0"]+b["tw"]-ox0, b["r0"]+b["th"]-oy0, border="red", lwd=2)
  title(sprintf("A %s in overview (NCC=%.2f)", sid, L$ncc), cex.main=0.95)
  # B: BF crop + footprints
  drawbf <- function(t){ plot.new(); plot.window(c(1,cw),c(ch,1),asp=1); rasterImage(bfm/max(bfm),1,ch,cw,1); title(t,cex.main=0.95) }
  drawbf("B native BF + MSI tissue(grn)/IF DAPI(red)")
  if (!is.null(tnf)) points(tnf[,1]-cx0+1, tnf[,2]-cy0+1, pch=15, cex=0.3, col=adjustcolor("green3",0.5))
  if (nrow(fgf)) points((fgf[,2]-1)*GRID_D+1, (fgf[,1]-1)*GRID_D+1, pch=15, cex=0.25, col=adjustcolor("red",0.4))
  umpx <- BF_UMPX; bp <- 100/umpx; segments(cw-bp-10,ch-14,cw-10,ch-14,lwd=3,col="black"); text(cw-bp/2-10,ch-26,"100 um",cex=0.7)
  # C: 4-ch IF composite over BF
  drawbf("C IF composite (DAPI/b-cat/ZO-1/F-actin)")
  R <- pmin(ZO1 + FAC, 1); Gc <- GFP; Bc <- pmin(DAPI + FAC, 1)
  fin <- is.finite(DAPI); col <- matrix("#00000000", length(gy), length(gx))
  col[fin] <- rgb(R[fin], Gc[fin], Bc[fin], alpha=0.85)
  rasterImage(as.raster(col), 1, ch, cw, 1, interpolate=FALSE)
  # D: MSI citrate (viridis 0.65) + DAPI+ZO-1 IF
  drawbf(sprintf("D MSI citrate [M-H]- (%.0f%%) + DAPI(blu)/ZO-1(red)", OVERLAY_ALPHA*100))
  ramp <- viridisLite::viridis(64); finI <- is.finite(ions)
  if (any(finI)) { cidx <- pmin(pmax(round(ions*63)+1,1),64); cm <- matrix("#00000000",length(gy),length(gx)); cm[finI] <- adjustcolor(ramp[cidx[finI]], alpha.f=OVERLAY_ALPHA); rasterImage(as.raster(cm),1,ch,cw,1,interpolate=FALSE) }
  cif <- matrix("#00000000", length(gy), length(gx)); ifm <- is.finite(DAPI) & (DAPI>0.3 | ZO1>0.3)
  cif[ifm] <- rgb(pmin(ZO1[ifm],1), pmin(DAPI[ifm],1), pmin(DAPI[ifm],1), alpha=0.5)   # DAPI cyan (G+B)
  rasterImage(as.raster(cif), 1, ch, cw, 1, interpolate=FALSE)
  mtext(sprintf("%s  (sec%d <-> MSI %s)", sid, r$secn, paste(msi_sids_for_section(r$msi_slide,r$secn),collapse="+")), line=-1.0, cex=0.8, font=2)
  cat(sprintf("[93] %s done (crop %dx%d)\n", sid, cw, ch))
}
dev.off()
cat("[93] DONE -> step3_native_overlay.pdf\n")
