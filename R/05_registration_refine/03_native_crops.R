#!/usr/bin/env Rscript
# ============================================================================
# 03_native_crops.R - STEP 4: native-res .nd2 optical crops + MSI overlay.
#   Compose MSI -> JPG (refined, R/05_registration_refine/01_refine_jpg.R) with JPG -> .nd2 (V-flip + scale + offset,
#   R/05_registration_refine/02_jpg_to_nd2_offset.R) -> MSI -> native .nd2. Per section: read the native-res .nd2 crop, do a
#   TINY (+/-3 MSI px) polish on the sharp .nd2 organoid darkness, save the crop +
#   a 3-panel page (native BF | +is_tissue outline | +MSI citrate projected on BF).
#
# Input : cache/register/jpgxform_<sid>.rds, cache/register/jpg2nd2off_<slide>.rds
# Output: cache/register/nd2final_<sid>.rds, figures/registration/crops/optical_<sid>.png
#         figures/registration/registration_native.pdf, results/registration/native_summary.csv
# Usage : Rscript R/05_registration_refine/03_native_crops.R [sid ...]
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
suppressPackageStartupMessages({ library(Cardinal); library(RBioFormats); library(png); library(viridisLite) })
REG_CACHE <- file.path(CACHE_DIR,"register"); REG_RES <- file.path(RES_DIR,"registration"); REG_FIG <- file.path(FIG_DIR,"registration")
CROP_DIR <- file.path(REG_FIG,"crops"); dir.create(CROP_DIR, showWarnings=FALSE, recursive=TRUE)
slide_of <- function(sid) if (grepl("sl6A",sid)) "sl6A" else "sl4A"
nd2_path <- function(sl) if (sl=="sl6A") file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"06102026_AO_0h_sl6A.nd2") else file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"06102026_AO_20h_sl4A.nd2")
POLISH_MSI <- 3

args <- commandArgs(trailingOnly=TRUE)
mse <- readRDS(TISSUE_MSE); pd <- as.data.frame(pixelData(mse)); mzs <- mz(mse)
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_citrate.R"))
val_cit <- citrate_onto_pd(pd)   # anchored citrate: raw imzML +-CITRATE_WIN_PPM, TIC-norm (was grid feat 191.0217)
SIDS <- if (length(args)) args else { s<-levels(pixelData(mse)$sample_id); if(is.null(s)) sort(unique(as.character(pd$sample_id))) else s }
offc <- list(); getoff <- function(sl){ if(is.null(offc[[sl]])) offc[[sl]]<<-readRDS(file.path(REG_CACHE,sprintf("jpg2nd2off_%s.rds",sl))); offc[[sl]] }
ndd <- list(); ndim <- function(sl){ if(is.null(ndd[[sl]])){ cm<-coreMetadata(read.metadata(nd2_path(sl)),series=1); ndd[[sl]]<<-c(cm$sizeX,cm$sizeY)}; ndd[[sl]] }

msi_to_nd2 <- function(msi, Bj, O){ jp<-apply_affine(Bj,msi); cbind(O$scale*jp[,1]+O$tx, O$scale*(O$jpg_H-jp[,2])+O$ty) }
inv_affine <- function(B,NX,NY){ M<-rbind(B[1,],B[2,]); t0<-B[3,]; xy<-solve(M, rbind(NX-t0[1],NY-t0[2])); list(x=xy[1,],y=xy[2,]) }

summ <- list(); pages <- list()
for (sid in SIDS) {
  sl<-slide_of(sid); R<-readRDS(file.path(REG_CACHE,sprintf("jpgxform_%s.rds",sid))); O<-getoff(sl); SD<-ndim(sl); SX<-SD[1]; SY<-SD[2]
  sub<-pd[as.character(pd$sample_id)==sid,c("x","y","is_tissue")]; W<-max(sub$x); H<-max(sub$y)
  nd0 <- msi_to_nd2(cbind(sub$x,sub$y), R$B_msi_jpg, O)
  B0 <- fit_affine(cbind(sub$x,sub$y), nd0); smn <- sqrt(abs(B0[1,1]*B0[2,2]-B0[1,2]*B0[2,1])); f<-max(1L,round(smn))

  # crop bbox + margin, read native nd2
  mg <- round(10*smn)
  cx0<-max(1,floor(min(nd0[,1])-mg)); cy0<-max(1,floor(min(nd0[,2])-mg)); cx1<-min(SX,ceiling(max(nd0[,1])+mg)); cy1<-min(SY,ceiling(max(nd0[,2])+mg))
  img<-read.image(nd2_path(sl),series=1,subset=list(x=cx0:cx1,y=cy0:cy1),normalize=TRUE); a<-as.array(img); om<-if(length(dim(a))==2)a else a[,,1]; om<-t(om)  # [row=y,col=x]
  Dn<-ds_mean(1-om/max(om), f); Hn<-nrow(Dn); Wn<-ncol(Dn)
  nd2dn<-function(NX,NY){ p<-cbind((NX-cx0)/f+1,(NY-cy0)/f+1); c<-pmin(pmax(round(p[,1]),1),Wn); r<-pmin(pmax(round(p[,2]),1),Hn); Dn[cbind(r,c)] }
  tis<-sub[sub$is_tissue,]; bg<-sub[!sub$is_tissue,]; if(nrow(bg)>600)bg<-bg[sample(nrow(bg),600),]
  ndt<-msi_to_nd2(cbind(tis$x,tis$y),R$B_msi_jpg,O); ndb<-msi_to_nd2(cbind(bg$x,bg$y),R$B_msi_jpg,O); cX<-mean(ndt[,1]); cY<-mean(ndt[,2])
  sc0<-mean(nd2dn(ndt[,1],ndt[,2]))-mean(nd2dn(ndb[,1],ndb[,2]))
  # tiny polish (translation only)
  best<-list(sc=sc0,tx=0,ty=0)
  for(tx in seq(-POLISH_MSI,POLISH_MSI,1)*smn) for(ty in seq(-POLISH_MSI,POLISH_MSI,1)*smn){ s<-mean(nd2dn(ndt[,1]+tx,ndt[,2]+ty))-mean(nd2dn(ndb[,1]+tx,ndb[,2]+ty)); if(s>best$sc) best<-list(sc=s,tx=tx,ty=ty) }
  Bf <- fit_affine(cbind(sub$x,sub$y), cbind(nd0[,1]+best$tx, nd0[,2]+best$ty))
  saveRDS(list(sid=sid,slide=sl,B_msi_nd2=Bf,scale_msi_nd2=smn,polish=c(best$tx,best$ty)/smn,contrast=best$sc,crop=c(cx0,cy0,cx1,cy1)),
          file.path(REG_CACHE,sprintf("nd2final_%s.rds",sid)))
  png::writePNG(om/max(om), file.path(CROP_DIR,sprintf("optical_%s.png",sid)))
  summ[[sid]]<-data.frame(sample_id=sid,slide=sl,polish_x=round(best$tx/smn,1),polish_y=round(best$ty/smn,1),contrast=round(best$sc,4),stringsAsFactors=FALSE)
  pages[[sid]]<-list(om=om,Bf=Bf,cx0=cx0,cy0=cy0,smn=smn,W=W,H=H)
  cat(sprintf("[34] %s: crop %dx%d polish(%.1f,%.1f)MSIpx contrast=%.3f\n", sid, ncol(om),nrow(om), best$tx/smn,best$ty/smn, best$sc))
}
summ_df<-do.call(rbind,summ); write.csv(summ_df, file.path(REG_RES,"native_summary.csv"), row.names=FALSE)

# ---- report: per section native BF | +mask | +citrate ----
pdf(file.path(REG_FIG,"registration_native.pdf"), width=13, height=5.2)
for (sid in names(pages)) {
  G<-pages[[sid]]; om<-G$om; B<-G$Bf; cw<-ncol(om); ch<-nrow(om); W<-G$W; H<-G$H; smn<-G$smn
  sub<-pd[as.character(pd$sample_id)==sid,c("x","y","is_tissue")]; tis<-sub[sub$is_tissue,]
  secval<-val_cit[as.character(pd$sample_id)==sid]
  D<-3L; gx<-seq(1,cw,by=D); gy<-seq(1,ch,by=D); NX<-G$cx0+rep(gx,times=length(gy))-1; NY<-G$cy0+rep(gy,each=length(gx))-1
  ms<-inv_affine(B,NX,NY); inb<-ms$x>=0.5&ms$x<=W+0.5&ms$y>=0.5&ms$y<=H+0.5
  lut<-matrix(NA_real_,W,H); tt<-sub$is_tissue; lut[cbind(sub$x[tt],sub$y[tt])]<-secval[tt]
  val<-rep(NA_real_,length(NX)); if(any(inb)){ xi<-pmin(pmax(round(ms$x[inb]),1),W); yi<-pmin(pmax(round(ms$y[inb]),1),H); val[inb]<-lut[cbind(xi,yi)] }
  vm<-matrix(val,length(gy),length(gx),byrow=TRUE); hi<-quantile(vm[is.finite(vm)&vm>0],0.995,na.rm=TRUE); if(!is.finite(hi)||hi<=0)hi<-1; vs<-pmin(vm/hi,1)
  layout(matrix(1:3,nrow=1)); par(mar=c(1,1,3,1))
  drawbf<-function(t){ plot.new(); plot.window(c(1,cw),c(ch,1),asp=1); rasterImage(om/max(om),1,ch,cw,1); title(t,cex.main=0.95) }
  drawbf(sprintf("%s native BF", sid)); umpx<-10/smn; bp<-100/umpx; segments(cw-bp-10,ch-14,cw-10,ch-14,lwd=3,col="black"); text(cw-bp/2-10,ch-28,"100 um",cex=0.7)
  drawbf("+ is_tissue (refined)"); tn<-apply_affine(B,cbind(tis$x,tis$y)); points(tn[,1]-G$cx0+1,tn[,2]-G$cy0+1,pch=15,cex=0.35,col=adjustcolor("green3",0.5))
  drawbf("+ MSI citrate [M-H]- (native)"); ramp<-viridisLite::viridis(64); fin<-is.finite(vs)&vs>0.05
  if(any(fin)){ cidx<-pmin(pmax(round(vs*63)+1,1),64); cmat<-matrix("#00000000",length(gy),length(gx)); cmat[fin]<-adjustcolor(ramp[cidx[fin]],alpha.f=0.7); rasterImage(as.raster(cmat),1,ch,cw,1,interpolate=FALSE) }
  mtext(sid, outer=FALSE, line=-1.2, cex=0.9, font=2)
}
dev.off(); cat(sprintf("[34] DONE -> registration_native.pdf (%d sections)\n", length(pages))); print(summ_df, row.names=FALSE)
