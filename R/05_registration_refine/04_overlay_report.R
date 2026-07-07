#!/usr/bin/env Rscript
# ============================================================================
# 04_overlay_report.R - registration overlay report with HOMOGENEOUS-opacity
#   MSI/brightfield blending (user request). The MSI ion image is drawn over the
#   native .nd2 brightfield at ONE constant opacity everywhere (no intensity
#   threshold, no per-pixel alpha variation); colour still encodes intensity.
#   Report-only: reads saved native crops (R/05_registration_refine/03_native_crops.R PNGs) + transforms - no .nd2 read.
#
# Tunables: OVERLAY_ALPHA (constant opacity), ION_MZ, colormap.
# Input : cache/register/nd2final_<sid>.rds, figures/registration/crops/optical_<sid>.png
# Output: figures/registration/registration_native.pdf
# Usage : Rscript R/05_registration_refine/04_overlay_report.R
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
suppressPackageStartupMessages({ library(Cardinal); library(png); library(viridisLite) })
REG_CACHE <- file.path(CACHE_DIR,"register"); REG_FIG <- file.path(FIG_DIR,"registration"); CROP_DIR <- file.path(REG_FIG,"crops")

OVERLAY_ALPHA <- 0.65         # <-- constant opacity of the MSI layer, same all across (user-set)
ION_MZ <- 191.0217            # citrate [M-H]-
GRID_D <- 2L                  # native-pixel step for the overlay raster

mse <- readRDS(TISSUE_MSE); pd <- as.data.frame(pixelData(mse)); mzs <- mz(mse)
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_citrate.R"))
val_ion <- citrate_onto_pd(pd)   # anchored citrate (raw imzML, TIC-norm); ION_MZ kept for the title
SIDS <- levels(pixelData(mse)$sample_id); if (is.null(SIDS)) SIDS <- sort(unique(as.character(pd$sample_id)))
inv_affine <- function(B,NX,NY){ M<-rbind(B[1,],B[2,]); t0<-B[3,]; xy<-solve(M, rbind(NX-t0[1],NY-t0[2])); list(x=xy[1,],y=xy[2,]) }
ramp <- viridisLite::viridis(64)

pdf(file.path(REG_FIG,"registration_native.pdf"), width=13, height=5.2)
for (sid in SIDS) {
  xf <- file.path(REG_CACHE,sprintf("nd2final_%s.rds",sid)); pf <- file.path(CROP_DIR,sprintf("optical_%s.png",sid))
  if (!file.exists(xf)||!file.exists(pf)) next
  X <- readRDS(xf); B <- X$B_msi_nd2; smn <- X$scale_msi_nd2; cx0 <- X$crop[1]; cy0 <- X$crop[2]
  om <- png::readPNG(pf); if (length(dim(om))==3) om <- om[,,1]; cw <- ncol(om); ch <- nrow(om)
  sub <- pd[as.character(pd$sample_id)==sid, c("x","y","is_tissue")]; W<-max(sub$x); H<-max(sub$y); tis<-sub[sub$is_tissue,]
  secval <- val_ion[as.character(pd$sample_id)==sid]

  # project ion (ON-TISSUE) onto the native crop grid
  gx<-seq(1,cw,by=GRID_D); gy<-seq(1,ch,by=GRID_D)
  NX<-cx0+rep(gx,times=length(gy))-1; NY<-cy0+rep(gy,each=length(gx))-1
  ms<-inv_affine(B,NX,NY); inb<-ms$x>=0.5&ms$x<=W+0.5&ms$y>=0.5&ms$y<=H+0.5
  lut<-matrix(NA_real_,W,H); lut[cbind(sub$x,sub$y)]<-secval   # WHOLE MSI image (full measured area)
  val<-rep(NA_real_,length(NX)); if(any(inb)){ xi<-pmin(pmax(round(ms$x[inb]),1),W); yi<-pmin(pmax(round(ms$y[inb]),1),H); val[inb]<-lut[cbind(xi,yi)] }
  vm<-matrix(val,length(gy),length(gx),byrow=TRUE)
  hi<-quantile(vm[is.finite(vm)&vm>0],0.995,na.rm=TRUE); if(!is.finite(hi)||hi<=0)hi<-1; vs<-pmin(vm/hi,1)

  layout(matrix(1:3,nrow=1)); par(mar=c(1,1,3,1))
  drawbf<-function(t){ plot.new(); plot.window(c(1,cw),c(ch,1),asp=1); rasterImage(om/max(om),1,ch,cw,1); title(t,cex.main=0.95) }
  # (1) native BF + scale bar
  drawbf(sprintf("%s native BF", sid)); umpx<-10/smn; bp<-100/umpx
  segments(cw-bp-10,ch-14,cw-10,ch-14,lwd=3,col="black"); text(cw-bp/2-10,ch-28,"100 um",cex=0.7)
  # (2) BF + refined is_tissue outline
  drawbf("+ is_tissue (refined)"); tn<-apply_affine(B,cbind(tis$x,tis$y)); points(tn[,1]-cx0+1,tn[,2]-cy0+1,pch=15,cex=0.35,col=adjustcolor("green3",0.5))
  # (3) BF + MSI ion at HOMOGENEOUS opacity (constant alpha for every MSI pixel; colour = intensity)
  drawbf(sprintf("+ MSI citrate [M-H]- (%.0f%% opacity, whole MSI image)", OVERLAY_ALPHA*100))
  fin <- is.finite(vs)                                   # whole measured MSI area (tissue + gel), no threshold
  if (any(fin)) {
    cidx <- pmin(pmax(round(vs*63)+1,1),64)
    cmat <- matrix("#00000000", length(gy), length(gx))
    cmat[fin] <- adjustcolor(ramp[cidx[fin]], alpha.f = OVERLAY_ALPHA)
    rasterImage(as.raster(cmat), 1, ch, cw, 1, interpolate=FALSE)
  }
  mtext(sid, line=-1.2, cex=0.9, font=2)
}
dev.off()
cat(sprintf("[35] DONE -> registration_native.pdf (uniform MSI opacity = %.0f%%)\n", OVERLAY_ALPHA*100))
