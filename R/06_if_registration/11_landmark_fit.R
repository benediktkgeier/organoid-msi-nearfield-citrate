#!/usr/bin/env Rscript
# ============================================================================
# 11_landmark_fit.R - fit per-section IF -> brightfield transforms from the
#   user's organoid-center landmark marks on landmark_sheets.pdf.
#
#   Workflow: user places matching NUMBERED point marks at organoid centers on
#   BOTH panels of each page (same number = same organoid). A PyMuPDF step writes
#   results/if_registration/landmark_points_raw.csv (page, number, ptx, pty in PDF pt,
#   top-down). This script maps each point to image pixels via the captured panel
#   boxes (landmark_panels.rds), pairs IF<->BF by number per section, and fits:
#       >=3 pairs  -> affine (handles rotation/scale/reflection/shear)
#        2 pairs   -> similarity (translation + rotation + uniform scale)
#        1 pair    -> translation only at the known optical scale (0.197)
#   Output: cache/register_if/iftobf_lm_<sid>.rds (B_if_bf) + a validation overlay.
#
# Run AFTER annotation:
#   python ... -> landmark_points_raw.csv   (extractor; see run note)
#   Rscript R/06_if_registration/11_landmark_fit.R
# ============================================================================

if (Sys.getenv("JAVA_HOME") == "") Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jre1.8.0_491")
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/if_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register_if.R"))
suppressPackageStartupMessages({ library(RBioFormats); library(EBImage); library(Cardinal) })

LP <- readRDS(file.path(IF_CACHE,"landmark_panels.rds")); panels <- LP$panels
pts <- read.csv(file.path(IF_RES,"landmark_points_raw.csv"), stringsAsFactors=FALSE)  # page,number,ptx,pty (top-down)
.mse<-readRDS(TISSUE_MSE); .pd<-as.data.frame(pixelData(.mse))
msi_pts<-function(A){ sids<-grep(sprintf("_%s_",A),unique(as.character(.pd$sample_id)),value=TRUE)
  do.call(rbind,lapply(sids,function(s){B<-readRDS(file.path(REG_CACHE,sprintf("nd2final_%s.rds",s)))$B_msi_nd2; sub<-.pd[as.character(.pd$sample_id)==s&.pd$is_tissue,c("x","y")]; apply_affine(B,cbind(sub$x,sub$y))})) }
lut_apply<-function(ch,L) pmin(pmax((ch-L["lo"])/(L["hi"]-L["lo"]),0),1)^(1/L["g"])
inv_aff<-function(B,q) sweep(q,2,B[3,]) %*% solve(B[1:2,])
page2sid <- setNames(sapply(panels,function(p)p$sid), sapply(panels,function(p)p$page))

# map a page point -> (panel, native px) using the captured panel boxes
map_point <- function(P, ptx, pty){
  inbox<-function(b) ptx>=min(b["px0"],b["px1"]) && ptx<=max(b["px0"],b["px1"]) && pty>=min(b["pyt"],b["pyb"]) && pty<=max(b["pyt"],b["pyb"])
  if (inbox(P$IF)){ fx<-(ptx-P$IF["px0"])/(P$IF["px1"]-P$IF["px0"]); fy<-(pty-P$IF["pyt"])/(P$IF["pyb"]-P$IF["pyt"])
    return(list(panel="IF", x=fx*P$IF["ifNX"], y=fy*P$IF["ifNY"])) }
  if (inbox(P$BF)){ fx<-(ptx-P$BF["px0"])/(P$BF["px1"]-P$BF["px0"]); fy<-(pty-P$BF["pyt"])/(P$BF["pyb"]-P$BF["pyt"])
    return(list(panel="BF", x=P$BF["cx0"]+fx*(P$BF["cx1"]-P$BF["cx0"]), y=P$BF["cy0"]+fy*(P$BF["cy1"]-P$BF["cy0"]))) }
  NULL
}
fit_lm <- function(src, dst){ n<-nrow(src); s0<-HR_UMPX/BF_UMPX
  if (n>=3) return(list(B=fit_affine(src,dst), mode="affine"))
  if (n==2){ ds<-src[2,]-src[1,]; dd<-dst[2,]-dst[1,]; s<-sqrt(sum(dd^2)/sum(ds^2))
    th<-atan2(dd[2],dd[1])-atan2(ds[2],ds[1]); R<-rbind(c(cos(th),sin(th)),c(-sin(th),cos(th))); Mr<-s*R
    return(list(B=rbind(Mr, colMeans(dst)-colMeans(src)%*%Mr), mode="similarity")) }
  Mr<-diag(s0,2); list(B=rbind(Mr, dst[1,]-src[1,]%*%Mr), mode="translation")
}

summ<-list(); pageinfo<-list()
for (pg in sort(unique(pts$page))){ sid<-page2sid[[as.character(pg)]]; if(is.null(sid)) next; P<-panels[[sid]]
  pp<-pts[pts$page==pg,]; IFm<-list(); BFm<-list()
  for (j in seq_len(nrow(pp))){ m<-map_point(P, pp$ptx[j], pp$pty[j]); if(is.null(m)) next
    if(m$panel=="IF") IFm[[as.character(pp$number[j])]]<-c(m$x,m$y) else BFm[[as.character(pp$number[j])]]<-c(m$x,m$y) }
  common<-intersect(names(IFm),names(BFm)); if(length(common)<1){ cat(sprintf("[99b] %s: no paired marks - skip\n",sid)); next }
  src<-do.call(rbind,IFm[common]); dst<-do.call(rbind,BFm[common])
  fit<-fit_lm(src,dst); B<-fit$B; res<-if(nrow(src)>=2) sqrt(mean(rowSums((apply_affine(B,src)-dst)^2)))*BF_UMPX else NA
  saveRDS(list(sid=sid,A=P$A,n=P$n,B_if_bf=B,mode=fit$mode,n_pairs=length(common),rmse_um=res,
               src=src,dst=dst,rect=P$BF[c("cx0","cy0","cx1","cy1")]), file.path(IF_CACHE,sprintf("iftobf_lm_%s.rds",sid)))
  summ[[sid]]<-data.frame(sid=sid,A=P$A,n=P$n,n_pairs=length(common),mode=fit$mode,rmse_um=round(res,1))
  pageinfo[[sid]]<-list(P=P,B=B,src=src,dst=dst)
  cat(sprintf("[99b] %s: %d pairs, %s, RMSE=%.1f um\n",sid,length(common),fit$mode,ifelse(is.na(res),0,res)))
}
if(length(summ)) write.csv(do.call(rbind,summ), file.path(IF_RES,"if_to_bf_landmark_summary.csv"), row.names=FALSE)

# ---- validation overlay: IF (DAPI+ZO-1) projected on brightfield + MSI + marks ----
pdf(file.path(IF_FIG,"if_to_bf_landmark.pdf"), width=13, height=6.5)
for (sid in names(pageinfo)){ G<-pageinfo[[sid]]; P<-G$P; A<-P$A; LUT<-IF_LUT[[P$slb]]; B<-G$B
  cx0<-P$BF["cx0"]; cy0<-P$BF["cy0"]; cx1<-P$BF["cx1"]; cy1<-P$BF["cy1"]
  img<-read.image(IF_SLIDES$msi_bf[IF_SLIDES$msi_slide==A],series=1,subset=list(x=cx0:cx1,y=cy0:cy1),normalize=TRUE)
  a<-as.array(img); om<-if(length(dim(a))==2)a else a[,,1]; om<-t(om); bfd<-om/max(om); cw<-ncol(bfd); ch<-nrow(bfd)
  ALL<-readRDS(file.path(IF_CACHE,sprintf("hr4raw_%s.rds",sid)))
  Dn<-3L; gx<-seq(1,cw,by=Dn); gy<-seq(1,ch,by=Dn); NX<-cx0+rep(gx,times=length(gy))-1; NY<-cy0+rep(gy,each=length(gx))-1
  ifp<-inv_aff(B,cbind(NX,NY)); fx<-ifp[,1]/F_HR+0.5; fy<-ifp[,2]/F_HR+0.5; W5<-ncol(ALL); H5<-nrow(ALL); inb<-fx>=1&fx<=W5&fy>=1&fy<=H5
  sd<-function(M){ vv<-rep(NA_real_,length(NX)); if(any(inb)) vv[inb]<-M[cbind(pmin(pmax(round(fy[inb]),1),H5),pmin(pmax(round(fx[inb]),1),W5))]; matrix(vv,length(gy),length(gx),byrow=TRUE) }
  DAPI<-lut_apply(sd(ALL[,,IF_CH["DAPI"]]),LUT$DAPI); ZO1<-lut_apply(sd(ALL[,,IF_CH["ZO1"]]),LUT$ZO1)
  MP<-msi_pts(A); ib<-MP[,1]>=cx0&MP[,1]<=cx1&MP[,2]>=cy0&MP[,2]<=cy1
  layout(matrix(1:2,nrow=1)); par(mar=c(1,1,3,1))
  plot.new(); plot.window(c(1,cw),c(ch,1),asp=1); rasterImage(bfd,1,ch,cw,1)
  if(any(ib)) points(MP[ib,1]-cx0+1,MP[ib,2]-cy0+1,pch=15,cex=0.3,col=adjustcolor("#39ff14",0.5))
  dn<-apply_affine(B,G$dst*0+G$src); points(G$dst[,1]-cx0+1,G$dst[,2]-cy0+1,col="red",pch=3,cex=2,lwd=2)
  title(sprintf("%s sec%d  brightfield + MSI(green) + your BF marks(red)",A,P$n),cex.main=0.9)
  plot.new(); plot.window(c(1,cw),c(ch,1),asp=1); rasterImage(bfd,1,ch,cw,1)
  fin<-is.finite(DAPI); cmat<-matrix("#00000000",length(gy),length(gx)); cmat[fin]<-rgb(pmin(ZO1[fin],1),pmin(DAPI[fin],1),pmin(DAPI[fin],1),alpha=0.55); rasterImage(as.raster(cmat),1,ch,cw,1,interpolate=FALSE)
  if(any(ib)) points(MP[ib,1]-cx0+1,MP[ib,2]-cy0+1,pch=15,cex=0.22,col=adjustcolor("#39ff14",0.35))
  ifm<-apply_affine(B,G$src); points(ifm[,1]-cx0+1,ifm[,2]-cy0+1,col="yellow",pch=1,cex=2,lwd=2)
  title(sprintf("+ IF registered (DAPI cyan+ZO-1 red)  [%s, %d pairs]",summ[[sid]]$mode,summ[[sid]]$n_pairs),cex.main=0.9)
}
dev.off(); cat(sprintf("[99b] DONE -> if_to_bf_landmark.pdf (%d sections)\n", length(pageinfo)))
