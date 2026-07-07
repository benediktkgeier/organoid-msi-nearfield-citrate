#!/usr/bin/env Rscript
# ============================================================================
# 10_landmark_sheets.R - per-section LANDMARK-MARKING sheets for manual IF ->
#   brightfield registration. One section per page, TWO EQUAL-SIZE panels (each
#   filling half the page for visibility - display scale differs, that's fine):
#     LEFT  = IF DAPI (ch4) in GREYSCALE (structures clearest for organoid centres)
#     RIGHT = brightfield rectangle (native) + registered MSI outline (green)
#   The user places matching NUMBERED point marks at organoid CENTRES on BOTH
#   panels (same number = same organoid). R/06_if_registration/11_landmark_fit.R maps them to image pixels via the
#   panel boxes captured here (grconvert) and fits a per-section affine IF->BF.
#   The fit is independent of display scale (maps panel fraction -> native px).
#
# Output: figures/if_registration/landmark_sheets.pdf, cache/register_if/landmark_panels.rds
# ============================================================================

Sys.setenv(JAVA_HOME="C:/Program Files/Java/jre1.8.0_491")
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/if_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register_if.R"))
suppressPackageStartupMessages({ library(RBioFormats); library(EBImage); library(Cardinal) })

rect <- read.csv(file.path(IF_RES,"bf_rectangles.csv"), stringsAsFactors=FALSE)
slide_b <- function(A) if (A=="sl6A") "sl6b" else "sl4b"; grp_of <- function(A) if (A=="sl6A") "0h" else "20h"
.mse<-readRDS(TISSUE_MSE); .pd<-as.data.frame(pixelData(.mse))
msi_pts<-function(A){ sids<-grep(sprintf("_%s_",A),unique(as.character(.pd$sample_id)),value=TRUE)
  do.call(rbind,lapply(sids,function(s){B<-readRDS(file.path(REG_CACHE,sprintf("nd2final_%s.rds",s)))$B_msi_nd2; sub<-.pd[as.character(.pd$sample_id)==s&.pd$is_tissue,c("x","y")]; apply_affine(B,cbind(sub$x,sub$y))})) }
nat<-function(ch,lo=0.005,hi=0.995,g=0.75){ q<-quantile(ch[is.finite(ch)],c(lo,hi),na.rm=TRUE); if(!is.finite(q[2])||q[2]<=q[1]) q<-range(ch[is.finite(ch)]); pmin(pmax((ch-q[1])/(q[2]-q[1]),0),1)^g }
ds_disp<-function(m,target=1500){ f<-max(1,round(max(dim(m))/target)); if(f>1) ds_mean(m,f) else m }

PAGE_W<-1500; PAGE_H<-780
panels<-list(); ndc2ptx<-function(xn) xn*PAGE_W; ndc2pty<-function(yn) (1-yn)*PAGE_H
panel_box <- function(W,H){ list(px0=ndc2ptx(grconvertX(1,"user","ndc")), px1=ndc2ptx(grconvertX(W,"user","ndc")),
                                 pyt=ndc2pty(grconvertY(1,"user","ndc")), pyb=ndc2pty(grconvertY(H,"user","ndc"))) }
pdf(file.path(IF_FIG, Sys.getenv("LM_OUT","landmark_sheets.pdf")), width=PAGE_W/72, height=PAGE_H/72)
for (i in seq_len(nrow(rect))){ A<-rect$slide[i]; slb<-slide_b(A); n<-rect$section[i]; sid<-sprintf("AO_%s_%s_sec%d",grp_of(A),slb,n)
  hf<-file.path(IF_CACHE,sprintf("hr4raw_%s.rds",sid)); if(!file.exists(hf)){ cat("skip",sid,"\n"); next }
  ALL<-readRDS(hf); ifNX<-ncol(ALL)*F_HR; ifNY<-nrow(ALL)*F_HR
  IFg<-ds_disp(nat(ALL[,,IF_CH["DAPI"]]))          # greyscale DAPI, brightened
  cm<-coreMetadata(read.metadata(IF_SLIDES$msi_bf[IF_SLIDES$msi_slide==A]),series=1); SX<-cm$sizeX; SY<-cm$sizeY
  cx0<-max(1,round(rect$fx0[i]*SX)); cx1<-min(SX,round(rect$fx1[i]*SX)); cy0<-max(1,round(rect$fy0[i]*SY)); cy1<-min(SY,round(rect$fy1[i]*SY))
  img<-read.image(IF_SLIDES$msi_bf[IF_SLIDES$msi_slide==A],series=1,subset=list(x=cx0:cx1,y=cy0:cy1),normalize=TRUE)
  a<-as.array(img); om<-if(length(dim(a))==2)a else a[,,1]; om<-t(om); bfd<-om/max(om)

  layout(matrix(1:2,nrow=1)); par(mar=c(1,1,3,1), oma=c(0,0,2,0))
  # LEFT: IF greyscale, fills cell (asp=1)
  Wi<-ncol(IFg); Hi<-nrow(IFg); plot.new(); plot.window(c(1,Wi),c(Hi,1),asp=1); rasterImage(IFg,1,Hi,Wi,1)
  segments(Wi*0.03,Hi*0.97,Wi*0.03+1000/(HR_UMPX*F_HR),Hi*0.97,lwd=4,col="white"); text(Wi*0.03+500/(HR_UMPX*F_HR),Hi*0.93,"1 mm",col="white",cex=0.8)
  title("IF DAPI (greyscale) - 30x field",cex.main=1)
  IFb<-panel_box(Wi,Hi)
  # RIGHT: brightfield + MSI, fills cell
  cw<-ncol(bfd); ch<-nrow(bfd); plot.new(); plot.window(c(1,cw),c(ch,1),asp=1); rasterImage(bfd,1,ch,cw,1)
  MP<-msi_pts(A); ib<-MP[,1]>=cx0&MP[,1]<=cx1&MP[,2]>=cy0&MP[,2]<=cy1
  if(any(ib)) points(MP[ib,1]-cx0+1,MP[ib,2]-cy0+1,pch=15,cex=0.3,col=adjustcolor("#39ff14",0.15))  # faint fill
  if(sum(ib)>5){ lm<-MP[ib,,drop=FALSE]; bin<-6L; Wn<-max(3L,ceiling(cw/bin)); Hn<-max(3L,ceiling(ch/bin)); M<-matrix(0,Wn,Hn)
    xi<-pmin(pmax(round((lm[,1]-cx0)/bin)+1,1),Wn); yi<-pmin(pmax(round((lm[,2]-cy0)/bin)+1,1),Hn); M[cbind(xi,yi)]<-1
    M<-fillHull(closing(M>0,makeBrush(5,"disc"))); for(co in ocontour(bwlabel(M))){ if(nrow(co)<4) next; lines(c(co[,1]*bin,co[1,1]*bin),c(co[,2]*bin,co[1,2]*bin),col=adjustcolor("#39ff14",0.7),lwd=1.3) } }
  segments(cw*0.03,ch*0.97,cw*0.03+1000/BF_UMPX,ch*0.97,lwd=4,col="black"); text(cw*0.03+500/BF_UMPX,ch*0.93,"1 mm",cex=0.8)
  title("brightfield + MSI organoids (green)",cex.main=1)
  BFb<-panel_box(cw,ch)
  mtext(sprintf("%s sec %d   -   place matching NUMBERED marks at organoid CENTRES on BOTH panels (same number = same organoid)", A, n), outer=TRUE, line=0.2, font=2, cex=1)

  panels[[sid]]<-list(sid=sid,A=A,slb=slb,n=n,
                      IF=c(unlist(IFb), ifNX=ifNX, ifNY=ifNY),
                      BF=c(unlist(BFb), cx0=cx0,cy0=cy0,cx1=cx1,cy1=cy1), page=length(panels))
  cat(sprintf("[99] %s page %d\n", sid, length(panels)))
}
dev.off()
saveRDS(list(panels=panels, PAGE_W=PAGE_W, PAGE_H=PAGE_H), file.path(IF_CACHE,"landmark_panels.rds"))
cat(sprintf("[99] DONE -> landmark_sheets.pdf (%d pages) + landmark_panels.rds\n", length(panels)))
