#!/usr/bin/env Rscript
# ============================================================================
# 09_if_to_bf_register.R - register each consecutive IF section (B-slide) to its
#   brightfield rectangle (A-slide) by AUTO organoid-shape fit.
#   Source = IF DAPI (ch4) tissue mask.  Target = registered MSI is_tissue
#   footprint within the rectangle (robust organoid anchor in BF coords) AND the
#   brightfield organoid (texture mask) for the fine refinement.
#   Transform = moment-based similarity (centroid + principal axes + isotropic
#   scale, 4 flip/reflection candidates) -> pick best mask IoU -> small refine.
#   Serial sections differ, so this is organoid-level; overlay validates it.
#
# Output: cache/register_if/iftobf_<sid_if>.rds (B_if_bf), figures/if_registration/
#         if_to_bf_register.pdf  (IF projected on brightfield + MSI outline)
# Usage : Rscript R/06_if_registration/09_if_to_bf_register.R [pilot|all|sid_if ...]
# ============================================================================

if (Sys.getenv("JAVA_HOME") == "") Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jre1.8.0_491")
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/if_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register_if.R"))
suppressPackageStartupMessages({ library(RBioFormats); library(EBImage); library(Cardinal) })

rect <- read.csv(file.path(IF_RES,"bf_rectangles.csv"), stringsAsFactors=FALSE)
.mse <- readRDS(TISSUE_MSE); .pd <- as.data.frame(pixelData(.mse))
slide_b <- function(A) if (A=="sl6A") "sl6b" else "sl4b"; grp_of <- function(A) if (A=="sl6A") "0h" else "20h"
msi_pts <- function(A){ sids<-grep(sprintf("_%s_",A),unique(as.character(.pd$sample_id)),value=TRUE)
  do.call(rbind,lapply(sids,function(s){B<-readRDS(file.path(REG_CACHE,sprintf("nd2final_%s.rds",s)))$B_msi_nd2; sub<-.pd[as.character(.pd$sample_id)==s&.pd$is_tissue,c("x","y")]; apply_affine(B,cbind(sub$x,sub$y))})) }

args <- commandArgs(trailingOnly=TRUE)
PILOT <- c("AO_0h_sl6b_sec1","AO_0h_sl6b_sec3","AO_20h_sl4b_sec5")
SECS <- do.call(rbind, lapply(seq_len(nrow(rect)), function(i){ A<-rect$slide[i]; data.frame(A=A, slb=slide_b(A), n=rect$section[i], sid=sprintf("AO_%s_%s_sec%d",grp_of(A),slide_b(A),rect$section[i]), stringsAsFactors=FALSE) }))
sel <- if (length(args)==0||identical(args,"pilot")) PILOT else if (identical(args,"all")) SECS$sid else args
SECS <- SECS[SECS$sid %in% sel,,drop=FALSE]

# rigid + FIXED isotropic scale (s0 from known optical px sizes); fit rotation
# (PCA) + flip + translation (centroid), pick best mask IoU, small refine.
g_bin <- 12
dmask <- function(P, r){ Wn<-max(3L,ceiling((r[2,1]-r[1,1])/g_bin)); Hn<-max(3L,ceiling((r[2,2]-r[1,2])/g_bin)); M<-matrix(FALSE,Wn,Hn)
  xi<-pmin(pmax(round((P[,1]-r[1,1])/g_bin)+1,1),Wn); yi<-pmin(pmax(round((P[,2]-r[1,2])/g_bin)+1,1),Hn); M[cbind(xi,yi)]<-TRUE; M }
iou_of <- function(P, dst, rng){ rr<-rbind(pmin(apply(P,2,min),rng[1,]),pmax(apply(P,2,max),rng[2,]))
  a<-dmask(P,rr); b<-dmask(dst,rr); H<-min(nrow(a),nrow(b)); W<-min(ncol(a),ncol(b)); sum(a[1:H,1:W]&b[1:H,1:W])/max(1,sum(a[1:H,1:W]|b[1:H,1:W])) }
fit_rigid <- function(src, dst, s0){
  cs<-colMeans(src); cd<-colMeans(dst); rng<-apply(dst,2,range)
  Vs<-eigen(crossprod(sweep(src,2,cs)),symmetric=TRUE)$vectors; Vd<-eigen(crossprod(sweep(dst,2,cd)),symmetric=TRUE)$vectors
  best<-NULL
  for (f1 in c(1,-1)) for (f2 in c(1,-1)){ R<-Vd%*%diag(c(f1,f2))%*%t(Vs); M<-s0*t(R); B<-rbind(M, cd - cs%*%M)
    io<-iou_of(apply_affine(B,src),dst,rng); if(is.null(best)||io>best$io) best<-list(B=B,io=io,flip=c(f1,f2)) }
  # refine: small output rotation about target centroid + translation
  M0<-best$B[1:2,]; bestB<-best$B; bestio<-best$io
  for (th in seq(-12,12,3)*pi/180){ ct<-cos(th); st<-sin(th); Rr<-rbind(c(ct,st),c(-st,ct)); M<-M0%*%Rr; tt<-cd - cs%*%M
    for (tx in seq(-60,60,20)) for (ty in seq(-60,60,20)){ B<-rbind(M, tt+c(tx,ty))
      io<-iou_of(apply_affine(B,src),dst,rng); if(io>bestio){bestio<-io; bestB<-B} } }
  list(B=bestB, io=bestio, flip=best$flip, s=s0)
}

pages<-list(); summ<-list()
for (k in seq_len(nrow(SECS))){ A<-SECS$A[k]; slb<-SECS$slb[k]; n<-SECS$n[k]; sid<-SECS$sid[k]
  rs<-rect[rect$slide==A&rect$section==n,]; if(!nrow(rs)) next
  cm<-coreMetadata(read.metadata(IF_SLIDES$msi_bf[IF_SLIDES$msi_slide==A]),series=1); SX<-cm$sizeX; SY<-cm$sizeY
  bx0<-max(1,round(rs$fx0*SX)); bx1<-min(SX,round(rs$fx1*SX)); by0<-max(1,round(rs$fy0*SY)); by1<-min(SY,round(rs$fy1*SY))
  # source: IF DAPI tissue (ch4, F=5) -> IF native px
  ALL<-readRDS(file.path(IF_CACHE,sprintf("hr4raw_%s.rds",sid))); v<-norm01(ALL[,,IF_CH["DAPI"]])
  bl<-filter2(v, makeBrush(21,"disc")/sum(makeBrush(21,"disc")), boundary="replicate"); fg<-which(bl>quantile(bl,0.88),arr.ind=TRUE)
  if(nrow(fg)>6000) fg<-fg[sample(nrow(fg),6000),]
  src<-cbind((fg[,2]-0.5)*F_HR, (fg[,1]-0.5)*F_HR)   # hr4raw F=5 px -> IF native px
  # target: MSI footprint within rect (BF native) - organoid anchor in BF coords
  MP<-msi_pts(A); ib<-MP[,1]>=bx0&MP[,1]<=bx1&MP[,2]>=by0&MP[,2]<=by1; dst<-MP[ib,,drop=FALSE]
  if(nrow(dst)<20){ cat(sprintf("[98] %s: too few MSI target px - skip\n",sid)); next }
  fit<-fit_rigid(src, dst, HR_UMPX/BF_UMPX); B<-fit$B   # FIXED optical scale
  saveRDS(list(sid=sid, A=A, n=n, B_if_bf=B, flip=fit$flip, scale=fit$s, iou=fit$io, rect=c(bx0,by0,bx1,by1)),
          file.path(IF_CACHE, sprintf("iftobf_%s.rds",sid)))
  summ[[sid]]<-data.frame(sid=sid,A=A,n=n,scale=round(fit$s,4),iou=round(fit$io,3))
  pages[[sid]]<-list(A=A,n=n,sid=sid,B=B,bx0=bx0,by0=by0,bx1=bx1,by1=by1,dst=dst,ALL=ALL,slb=slb)
  cat(sprintf("[98] %s: scale=%.3f IoU=%.3f flip=(%d,%d)\n",sid,fit$s,fit$io,fit$flip[1],fit$flip[2]))
}
if(length(summ)) write.csv(do.call(rbind,summ), file.path(IF_RES,"if_to_bf_summary.csv"), row.names=FALSE)

# ---- overlay report: IF (DAPI cyan + ZO-1 red) projected on brightfield + MSI outline ----
lut_apply <- function(ch,L) pmin(pmax((ch-L["lo"])/(L["hi"]-L["lo"]),0),1)^(1/L["g"])
inv_aff <- function(B,q) sweep(q,2,B[3,]) %*% solve(B[1:2,])
pdf(file.path(IF_FIG,"if_to_bf_register.pdf"), width=13, height=6.5)
for(sid in names(pages)){ G<-pages[[sid]]; A<-G$A; LUT<-IF_LUT[[G$slb]]
  img<-read.image(IF_SLIDES$msi_bf[IF_SLIDES$msi_slide==A], series=1, subset=list(x=G$bx0:G$bx1, y=G$by0:G$by1), normalize=TRUE)
  a<-as.array(img); om<-if(length(dim(a))==2)a else a[,,1]; om<-t(om); bfd<-om/max(om); cw<-ncol(bfd); ch<-nrow(bfd)
  # project IF DAPI+ZO-1 onto BF crop via inverse B (BF native -> IF native -> IF F5 px)
  D<-3L; gx<-seq(1,cw,by=D); gy<-seq(1,ch,by=D); NX<-G$bx0+rep(gx,times=length(gy))-1; NY<-G$by0+rep(gy,each=length(gx))-1
  ifp<-inv_aff(G$B, cbind(NX,NY)); fx<-ifp[,1]/F_HR+0.5; fy<-ifp[,2]/F_HR+0.5; W5<-ncol(G$ALL); H5<-nrow(G$ALL)
  inb<-fx>=1&fx<=W5&fy>=1&fy<=H5
  sd<-function(M){ vv<-rep(NA_real_,length(NX)); if(any(inb)) vv[inb]<-M[cbind(pmin(pmax(round(fy[inb]),1),H5),pmin(pmax(round(fx[inb]),1),W5))]; matrix(vv,length(gy),length(gx),byrow=TRUE) }
  DAPI<-lut_apply(sd(G$ALL[,,IF_CH["DAPI"]]),LUT$DAPI); ZO1<-lut_apply(sd(G$ALL[,,IF_CH["ZO1"]]),LUT$ZO1)
  layout(matrix(1:2,nrow=1)); par(mar=c(1,1,3,1))
  plot.new(); plot.window(c(1,cw),c(ch,1),asp=1); rasterImage(bfd,1,ch,cw,1); title(sprintf("%s sec%d  brightfield + MSI outline (green)",A,G$n),cex.main=0.95)
  ib<-G$dst[,1]>=G$bx0&G$dst[,1]<=G$bx1&G$dst[,2]>=G$by0&G$dst[,2]<=G$by1; points(G$dst[ib,1]-G$bx0+1,G$dst[ib,2]-G$by0+1,pch=15,cex=0.3,col=adjustcolor("#39ff14",0.5))
  plot.new(); plot.window(c(1,cw),c(ch,1),asp=1); rasterImage(bfd,1,ch,cw,1)
  fin<-is.finite(DAPI); cmat<-matrix("#00000000",length(gy),length(gx)); cmat[fin]<-rgb(pmin(ZO1[fin],1),pmin(DAPI[fin],1),pmin(DAPI[fin],1),alpha=0.55); rasterImage(as.raster(cmat),1,ch,cw,1,interpolate=FALSE)
  points(G$dst[ib,1]-G$bx0+1,G$dst[ib,2]-G$by0+1,pch=15,cex=0.25,col=adjustcolor("#39ff14",0.35))
  title(sprintf("+ IF registered: DAPI(cyan)+ZO-1(red)   [IoU %.2f]",summ[[sid]]$iou),cex.main=0.95)
}
dev.off(); cat(sprintf("[98] DONE -> if_to_bf_register.pdf (%d sections)\n", length(pages)))
