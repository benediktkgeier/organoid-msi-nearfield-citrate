#!/usr/bin/env Rscript
# ============================================================================
# 12_dataset_pairs_hq.R - per-MSI-dataset, FULL-QUALITY side-by-side report.
#   For each MSI dataset, crop its brightfield region EXTENDED by 30% and show:
#     LEFT  = IF 3-channel composite (ZO-1=red ch2, b-catenin=green ch3,
#             DAPI=cyan ch4; b-catenin + F-actin/Cy5 EXCLUDED), read at NATIVE 30x resolution,
#             reprojected into the brightfield frame via the landmark transform,
#             at FULL opacity (no brightfield underneath).
#     RIGHT = native brightfield (full resolution) with the MSI-measured
#             rectangle drawn as a thin DASHED line.
#   Both panels = same region/orientation, rendered as large raster (full quality).
#
# Usage : Rscript R/06_if_registration/12_dataset_pairs_hq.R [test|all|sid_msi ...]
# Output: figures/if_registration/dataset_pairs_hq.pdf
#         figures/if_registration/crops/if_optical_<sid_msi>.png  (the LEFT panel's
#           native-res IF composite, per dataset, standalone PNG; no brightfield)
# ============================================================================

if (Sys.getenv("JAVA_HOME") == "") Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jre1.8.0_491")
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/if_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register_if.R"))
suppressPackageStartupMessages({ library(RBioFormats); library(Cardinal); library(viridisLite); library(EBImage) })

EXT<-0.30; OUT_CAP<-2600L; ION_MZ<-191.0217; CIT_ALPHA<-1.0; Dg<-2L
rect<-read.csv(file.path(IF_RES,"bf_rectangles.csv"), stringsAsFactors=FALSE)
mse<-readRDS(TISSUE_MSE); pd<-as.data.frame(pixelData(mse)); mzs<-mz(mse)
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_citrate.R"))
val_ion<-citrate_onto_pd(pd); ramp<-viridisLite::viridis(64)   # anchored citrate (raw imzML, TIC-norm); was grid feat ION_MZ
clip01<-function(m,p=0.995){ hi<-quantile(m[is.finite(m)&m>0],p,na.rm=TRUE); if(!is.finite(hi)||hi<=0)hi<-1; pmin(m/hi,1) }
slide_b<-function(A) if(A=="sl6A")"sl6b" else "sl4b"; grp_of<-function(A) if(A=="sl6A")"0h" else "20h"
msi_slide_of<-function(sid) if(grepl("sl6A",sid))"sl6A" else "sl4A"
inv_aff<-function(B,q) sweep(q,2,B[3,]) %*% solve(B[1:2,])
lut_apply<-function(ch,L) pmin(pmax((ch-L["lo"])/(L["hi"]-L["lo"]),0),1)^(1/L["g"])
bfdim<-list(); for(A in unique(rect$slide)){ cm<-coreMetadata(read.metadata(IF_SLIDES$msi_bf[IF_SLIDES$msi_slide==A]),series=1); bfdim[[A]]<-c(SX=cm$sizeX,SY=cm$sizeY) }
if_sid_for<-function(A,cx,cy){ cm<-bfdim[[A]]; rs<-rect[rect$slide==A,]
  for(i in seq_len(nrow(rs))) if(cx>=rs$fx0[i]*cm["SX"]&&cx<=rs$fx1[i]*cm["SX"]&&cy>=rs$fy0[i]*cm["SY"]&&cy<=rs$fy1[i]*cm["SY"]) return(sprintf("AO_%s_%s_sec%d",grp_of(A),slide_b(A),rs$section[i])); NULL }
hp_of<-function(sid_if){ slb<-if(grepl("sl6b",sid_if))"sl6b" else "sl4b"; n<-as.integer(sub(".*sec(\\d+)$","\\1",sid_if)); file.path(IF_DIR, paste0(IF_SLIDES$hr_prefix[IF_SLIDES$slide==slb], n, ".nd2")) }

args<-commandArgs(trailingOnly=TRUE)
SIDS<-sort(unique(as.character(pd$sample_id)))
if(length(args) && args[1]=="test") SIDS<-"AO_0h_sl6A_sec5b" else if(length(args) && args[1]!="all") SIDS<-args

pdf(file.path(IF_FIG, Sys.getenv("HQ_OUT","dataset_pairs_hq.pdf")), width=20, height=7.2)
for(sid_msi in SIDS){
  A<-msi_slide_of(sid_msi); xf<-file.path(REG_CACHE,sprintf("nd2final_%s.rds",sid_msi)); if(!file.exists(xf)) next
  Bm<-readRDS(xf)$B_msi_nd2; sub<-pd[as.character(pd$sample_id)==sid_msi,c("x","y")]; W<-max(sub$x); H<-max(sub$y)
  SX<-bfdim[[A]]["SX"]; SY<-bfdim[[A]]["SY"]
  nd<-apply_affine(Bm,cbind(sub$x,sub$y))                       # measured pixels -> BF native
  rx<-range(nd[,1]); ry<-range(nd[,2]); ex<-diff(rx)*EXT/2; ey<-diff(ry)*EXT/2
  cx0<-max(1,floor(rx[1]-ex)); cy0<-max(1,floor(ry[1]-ey)); cx1<-min(SX,ceiling(rx[2]+ex)); cy1<-min(SY,ceiling(ry[2]+ey))
  img<-read.image(IF_SLIDES$msi_bf[IF_SLIDES$msi_slide==A],series=1,subset=list(x=cx0:cx1,y=cy0:cy1),normalize=TRUE)
  a<-as.array(img); om<-if(length(dim(a))==2)a else a[,,1]; bfd<-t(om)/max(om); cw<-ncol(bfd); ch<-nrow(bfd)
  # MSI-measured rectangle corners (raster bbox) -> BF native -> crop-local
  mc<-apply_affine(Bm, cbind(c(0.5,W+0.5,W+0.5,0.5), c(0.5,0.5,H+0.5,H+0.5)))
  mcx<-mc[,1]-cx0+1; mcy<-mc[,2]-cy0+1

  # ---- IF composite (native, reprojected into BF frame) ----
  ctr<-colMeans(nd); sid_if<-if_sid_for(A,ctr[1],ctr[2]); IFrgb<-NULL
  if(!is.null(sid_if)){ tf<-file.path(IF_CACHE,sprintf("iftobf_lm_%s.rds",sid_if))
    if(file.exists(tf)){ Bif<-readRDS(tf)$B_if_bf; slb<-slide_b(A); LUT<-IF_LUT[[slb]]; hp<-hp_of(sid_if)
      hm<-coreMetadata(read.metadata(hp),series=1); ifNX<-hm$sizeX; ifNY<-hm$sizeY
      # output grid in BF-crop px at ~IF native resolution (capped)
      scl<-BF_UMPX/HR_UMPX; outW<-min(OUT_CAP, round(cw*scl)); outH<-round(outW*ch/cw)
      ogx<-seq(1,cw,length.out=outW); ogy<-seq(1,ch,length.out=outH)
      NXo<-cx0+rep(ogx,times=outH)-1; NYo<-cy0+rep(ogy,each=outW)-1
      ip<-inv_aff(Bif,cbind(NXo,NYo))                            # -> IF native px
      ib0<-floor(min(ip[,1]));ib1<-ceiling(max(ip[,1]));jb0<-floor(min(ip[,2]));jb1<-ceiling(max(ip[,2]))
      ib0<-max(1,ib0);jb0<-max(1,jb0);ib1<-min(ifNX,ib1);jb1<-min(ifNY,jb1)
      crp<-as.array(read.image(hp,series=1,subset=list(x=ib0:ib1,y=jb0:jb1),normalize=FALSE))  # [x,y,4]
      xi<-round(ip[,1])-ib0+1; yi<-round(ip[,2])-jb0+1; inb<-xi>=1&xi<=dim(crp)[1]&yi<=dim(crp)[2]&yi>=1
      smp<-function(c){ v<-rep(NA_real_,length(NXo)); if(any(inb)) v[inb]<-crp[cbind(xi[inb],yi[inb],c)]; v }
      Rr<-lut_apply(smp(IF_CH["ZO1"]),LUT$ZO1); Bb<-lut_apply(smp(IF_CH["DAPI"]),LUT$DAPI)   # ZO-1 + DAPI only (b-catenin dropped)
      f<-function(v){ v[is.na(v)]<-0; matrix(v,outH,outW,byrow=TRUE) }
      IFrgb<-array(c(f(Rr),f(Bb),f(Bb)),c(outH,outW,3)); rm(crp); gc(verbose=FALSE)   # R=ZO-1, G=B=DAPI -> DAPI cyan
      # standalone high-res IF composite crop (same raster as P1: ZO-1 red / DAPI cyan, no brightfield)
      png::writePNG(IFrgb, file.path(IF_FIG, "crops", sprintf("if_optical_%s.png", sid_msi))) } }

  # ---- citrate ion projected into BF crop frame ----
  gx<-seq(1,cw,by=Dg); gy<-seq(1,ch,by=Dg); NXc<-cx0+rep(gx,times=length(gy))-1; NYc<-cy0+rep(gy,each=length(gx))-1
  secval<-val_ion[as.character(pd$sample_id)==sid_msi]; lutm<-matrix(NA_real_,W,H); lutm[cbind(sub$x,sub$y)]<-secval
  msc<-inv_aff(Bm,cbind(NXc,NYc)); inbc<-msc[,1]>=0.5&msc[,1]<=W+0.5&msc[,2]>=0.5&msc[,2]<=H+0.5
  vv<-rep(NA_real_,length(NXc)); if(any(inbc)) vv[inbc]<-lutm[cbind(pmin(pmax(round(msc[inbc,1]),1),W),pmin(pmax(round(msc[inbc,2]),1),H))]
  vm<-matrix(vv,length(gy),length(gx),byrow=TRUE); vs<-clip01(vm)

  smn<-sqrt(abs(Bm[1,1]*Bm[2,2]-Bm[1,2]*Bm[2,1])); bp<-200/(MSI_PIXEL_UM/smn)
  layout(matrix(1:3,nrow=1)); par(mar=c(1,1,3,1))
  # P1: IF composite (full opacity, no brightfield)
  plot.new(); plot.window(c(0,1),c(1,0),asp=ch/cw)
  if(!is.null(IFrgb)){ rasterImage(IFrgb,0,1,1,0,interpolate=TRUE); title(sprintf("%s  -  IF: ZO-1(red)/DAPI(cyan)",disp_id(sid_msi)),cex.main=0.9)
  } else { rasterImage(matrix(0,2,2),0,1,1,0); title(sprintf("%s  -  IF (no transform)",disp_id(sid_msi)),cex.main=0.9) }
  # P2: brightfield + dashed MSI rectangle
  plot.new(); plot.window(c(1,cw),c(ch,1),asp=1); rasterImage(bfd,1,ch,cw,1)
  polygon(mcx,mcy,border="red",lty=2,lwd=1); title("brightfield (native) + MSI area (dashed)",cex.main=0.9)
  segments(cw-bp-10,ch-12,cw-10,ch-12,lwd=3,col="black"); text(cw-bp/2-10,ch-26,"200 um",cex=0.7)
  # P3: citrate ion on brightfield
  plot.new(); plot.window(c(1,cw),c(ch,1),asp=1); rasterImage(bfd,1,ch,cw,1)
  finc<-is.finite(vs); if(any(finc)){ cidx<-pmin(pmax(round(vs*63)+1,1),64); cm3<-matrix("#00000000",length(gy),length(gx)); cm3[finc]<-adjustcolor(ramp[cidx[finc]],alpha.f=CIT_ALPHA); rasterImage(as.raster(cm3),1,ch,cw,1,interpolate=FALSE) }
  # white dotted organoid segmentation outline (is_tissue boundary = signed-dist 0)
  tis<-pd[as.character(pd$sample_id)==sid_msi & pd$is_tissue, c("x","y")]
  if(nrow(tis)>5){ Mt<-matrix(0,W,H); Mt[cbind(tis$x,tis$y)]<-1; Mt<-fillHull(closing(Mt>0,makeBrush(3,"disc")))
    for(co in ocontour(bwlabel(Mt))){ if(nrow(co)<6) next; p<-apply_affine(Bm,cbind(co[,1]+0.5,co[,2]+0.5)); lines(c(p[,1]-cx0+1,p[1,1]-cx0+1),c(p[,2]-cy0+1,p[1,2]-cy0+1),col="white",lty=3,lwd=1.2) } }
  polygon(mcx,mcy,border="red",lty=2,lwd=1); title(sprintf("citrate [M-H]- %.4f ion (%.0f%%) + organoid seg (white dotted)",ION_MZ,CIT_ALPHA*100),cex.main=0.9)
  segments(cw-bp-10,ch-12,cw-10,ch-12,lwd=3,col="black"); text(cw-bp/2-10,ch-26,"200 um",cex=0.7)
  cat(sprintf("[101] %s: BF %dx%d  IF=%s\n",sid_msi,cw,ch,ifelse(is.null(IFrgb),"none",sid_if)))
}
dev.off(); cat("[101] DONE -> dataset_pairs_hq.pdf\n")
