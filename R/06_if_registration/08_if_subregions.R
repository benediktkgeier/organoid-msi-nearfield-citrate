#!/usr/bin/env Rscript
# ============================================================================
# 08_if_subregions.R - zoom report of the user-drawn subregions on roi_pairs_hq.
#   The user drew a small rectangle on each DAPI panel and each brightfield panel
#   (results/if_registration/subregions_raw.csv). This script:
#     1. reproduces the roi_pairs_hq layout and captures each panel's page-pt box
#        (grconvert) so the annotations map back to exact image pixels;
#     2. reads each IF subregion at NATIVE resolution (DAPI ch1 + ZO-1 ch3),
#        applies the slide's NIS LUTs, and overlays DAPI=cyan + ZO-1=red;
#     3. reads each brightfield subregion at native res (registration_native look)
#        with the registered MSI is_tissue outline (green);
#     4. renders ONE section pair per page, same physical scale, as large as
#        possible side by side.
#
# Input : subregions_raw.csv, bf_rectangles.csv, IF + BF .nd2, hrdapihq cache, TISSUE_MSE
# Output: figures/if_registration/roi_pairs_subregions.pdf
# ============================================================================

Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jre1.8.0_491")
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/if_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register_if.R"))
suppressPackageStartupMessages({ library(RBioFormats); library(EBImage); library(Cardinal) })

ann   <- read.csv(file.path(IF_RES, "subregions_raw.csv"), stringsAsFactors = FALSE)  # page,ax0,ay0,ax1,ay1 (pt, top-down)
rect  <- read.csv(file.path(IF_RES, "bf_rectangles.csv"), stringsAsFactors = FALSE)
slide_b <- function(A) if (A=="sl6A") "sl6b" else "sl4b"
grp_of  <- function(A) if (A=="sl6A") "0h" else "20h"
SLIDES  <- c("sl6A","sl4A")                       # page 0, page 1
HQF <- 2L; HR_PX_UM <- HR_UMPX*HQF; BF_PX_UM <- BF_UMPX; GAP_UM <- 400
PAGE_W <- 792; PAGE_H <- 1008                     # 11x14 in pt (roi_pairs_hq)
lut_apply <- function(ch, L) pmin(pmax((ch - L["lo"])/(L["hi"]-L["lo"]), 0), 1)^(1/L["g"])

.mse <- readRDS(TISSUE_MSE); .pd <- as.data.frame(pixelData(.mse))
msi_bf_points <- function(A){ sids <- grep(sprintf("_%s_",A), unique(as.character(.pd$sample_id)), value=TRUE)
  do.call(rbind, lapply(sids, function(s){ xf<-file.path(REG_CACHE,sprintf("nd2final_%s.rds",s)); if(!file.exists(xf)) return(NULL)
    B<-readRDS(xf)$B_msi_nd2; sub<-.pd[as.character(.pd$sample_id)==s & .pd$is_tissue,c("x","y")]; if(!nrow(sub)) return(NULL); apply_affine(B,cbind(sub$x,sub$y)) })) }

# geometry identical to roi_pairs_hq (R/06_if_registration/07_if_roi_pairs.R)
bfdim <- list(); for (A in SLIDES){ cm<-coreMetadata(read.metadata(IF_SLIDES$msi_bf[IF_SLIDES$msi_slide==A]),series=1); bfdim[[A]]<-c(SX=cm$sizeX,SY=cm$sizeY) }
hqref <- readRDS(file.path(IF_CACHE, sprintf("hrdapihq_%s.rds", paste0(IF_SLIDES$sid_prefix[1],1))))$m
HRNC <- ncol(hqref); HRNR <- nrow(hqref)            # HQ DAPI dims (F=2)
hr_W <- HRNC*HR_PX_UM; hr_H <- HRNR*HR_PX_UM
rect$bfW <- (rect$fx1-rect$fx0)*sapply(rect$slide,function(A)bfdim[[A]]["SX"])*BF_PX_UM
rect$bfH <- (rect$fy1-rect$fy0)*sapply(rect$slide,function(A)bfdim[[A]]["SY"])*BF_PX_UM
XSPAN <- hr_W+GAP_UM+max(rect$bfW); YSPAN <- max(hr_H, max(rect$bfH))

# ---- Part A: capture each panel's page-pt bbox by reproducing the layout ----
panels <- list()
ndc2pt <- function(xn,yn) c(x=xn*PAGE_W, y=(1-yn)*PAGE_H)   # ndc -> pt, top-down
pdf(tempfile(fileext=".pdf"), width=11, height=14)
for (pageidx in seq_along(SLIDES)){ A<-SLIDES[pageidx]; rs<-rect[rect$slide==A,]; rs<-rs[order(rs$section),]; SX<-bfdim[[A]]["SX"]; SY<-bfdim[[A]]["SY"]
  layout(matrix(1:6,nrow=6,byrow=TRUE)); par(mar=c(1,1,2.2,1))
  for (i in seq_len(nrow(rs))){ n<-rs$section[i]
    cx0<-max(1,round(rs$fx0[i]*SX)); cx1<-min(SX,round(rs$fx1[i]*SX)); cy0<-max(1,round(rs$fy0[i]*SY)); cy1<-min(SY,round(rs$fy1[i]*SY))
    bw<-(cx1-cx0+1)*BF_PX_UM; bh<-(cy1-cy0+1)*BF_PX_UM
    plot.new(); plot.window(c(0,XSPAN),c(0,YSPAN),asp=1)
    d0<-ndc2pt(grconvertX(0,"user","ndc"),       grconvertY(YSPAN,"user","ndc"))
    d1<-ndc2pt(grconvertX(hr_W,"user","ndc"),    grconvertY(YSPAN-hr_H,"user","ndc"))
    b0<-ndc2pt(grconvertX(hr_W+GAP_UM,"user","ndc"),    grconvertY(YSPAN,"user","ndc"))
    b1<-ndc2pt(grconvertX(hr_W+GAP_UM+bw,"user","ndc"), grconvertY(YSPAN-bh,"user","ndc"))
    panels[[sprintf("%d.%d.DAPI",pageidx-1,n)]] <- list(slide=A,sec=n,panel="DAPI",px0=d0["x"],px1=d1["x"],pyt=d0["y"],pyb=d1["y"])
    panels[[sprintf("%d.%d.BF",pageidx-1,n)]]   <- list(slide=A,sec=n,panel="BF",px0=b0["x"],px1=b1["x"],pyt=b0["y"],pyb=b1["y"],
                                                        cx0=cx0,cy0=cy0,cx1=cx1,cy1=cy1,bw=bw,bh=bh)
  }
}
dev.off()

# ---- Part B: assign each annotation to a panel + map -> native-image subregion ----
sub <- list()   # key slide.sec -> list(DAPI=c(nx0,ny0,nx1,ny1) hi-res native, BF=c(...) bf native)
for (r in seq_len(nrow(ann))){ pg<-ann$page[r]; acx<-(ann$ax0[r]+ann$ax1[r])/2; acy<-(ann$ay0[r]+ann$ay1[r])/2
  hit<-NULL; for (k in names(panels)){ P<-panels[[k]]; if (as.integer(strsplit(k,"\\.")[[1]][1])!=pg) next
    if (acx>=min(P$px0,P$px1) && acx<=max(P$px0,P$px1) && acy>=min(P$pyt,P$pyb) && acy<=max(P$pyt,P$pyb)){ hit<-P; break } }
  if (is.null(hit)){ cat(sprintf("[97] WARN annotation %d (page %d) not in any panel\n", r, pg)); next }
  key<-sprintf("%s.%d",hit$slide,hit$sec); if (is.null(sub[[key]])) sub[[key]]<-list()
  # fractions of the annotation rect within the panel pt-box
  fx0<-(ann$ax0[r]-hit$px0)/(hit$px1-hit$px0); fx1<-(ann$ax1[r]-hit$px0)/(hit$px1-hit$px0)
  fyt<-(ann$ay0[r]-hit$pyt)/(hit$pyb-hit$pyt); fyb<-(ann$ay1[r]-hit$pyt)/(hit$pyb-hit$pyt)
  fx<-sort(c(fx0,fx1)); fy<-sort(c(fyt,fyb))
  if (hit$panel=="DAPI"){ nx0<-round(fx[1]*HRNC*HQF); nx1<-round(fx[2]*HRNC*HQF); ny0<-round(fy[1]*HRNR*HQF); ny1<-round(fy[2]*HRNR*HQF)
    sub[[key]]$DAPI<-c(max(1,nx0),max(1,ny0),nx1,ny1) }
  else { nx0<-hit$cx0+round(fx[1]*(hit$cx1-hit$cx0)); nx1<-hit$cx0+round(fx[2]*(hit$cx1-hit$cx0)); ny0<-hit$cy0+round(fy[1]*(hit$cy1-hit$cy0)); ny1<-hit$cy0+round(fy[2]*(hit$cy1-hit$cy0))
    sub[[key]]$BF<-c(nx0,ny0,nx1,ny1) }
}

# ---- Part C: render one section pair per page, same scale, as large as possible ----
pdf(file.path(IF_FIG, Sys.getenv("SUB_OUT","roi_pairs_subregions.pdf")), width=14, height=8)
GAP2_UM <- 250
for (A in SLIDES){ slb<-slide_b(A); LUT<-IF_LUT[[slb]]
  hp_of <- function(n) file.path(IF_DIR, paste0(IF_SLIDES$hr_prefix[IF_SLIDES$slide==slb], n, ".nd2"))
  for (n in 1:6){ key<-sprintf("%s.%d",A,n); S<-sub[[key]]; if (is.null(S)||is.null(S$DAPI)||is.null(S$BF)){ cat(sprintf("[97] skip %s (missing sub)\n",key)); next }
    # IF subregion: native DAPI(ch1)+ZO-1(ch3) -> LUT -> cyan+red
    d<-S$DAPI; img<-read.image(hp_of(n), series=1, subset=list(x=d[1]:d[3], y=d[2]:d[4]), normalize=FALSE); a<-as.array(img)
    DAPI<-lut_apply(t(a[,,IF_CH["DAPI"]]), LUT$DAPI); ZO1<-lut_apply(t(a[,,IF_CH["ZO1"]]), LUT$ZO1)  # ch4=DAPI, ch2=ZO-1
    IFrgb<-array(0,c(nrow(DAPI),ncol(DAPI),3)); IFrgb[,,1]<-ZO1; IFrgb[,,2]<-DAPI; IFrgb[,,3]<-DAPI   # red=ZO-1, cyan=DAPI
    ifW<-(d[3]-d[1]+1)*HR_UMPX; ifH<-(d[4]-d[2]+1)*HR_UMPX
    # BF subregion: native brightfield + MSI outline
    b<-S$BF; bimg<-read.image(IF_SLIDES$msi_bf[IF_SLIDES$msi_slide==A], series=1, subset=list(x=b[1]:b[3], y=b[2]:b[4]), normalize=TRUE)
    ba<-as.array(bimg); bom<-if(length(dim(ba))==2)ba else ba[,,1]; bom<-t(bom); bfd<-bom/max(bom)
    bW<-(b[3]-b[1]+1)*BF_PX_UM; bH<-(b[4]-b[2]+1)*BF_PX_UM
    # layout in physical um, same scale, side by side
    Xs<-ifW+GAP2_UM+bW; Ys<-max(ifH,bH)
    par(mar=c(1,1,2.4,1)); plot.new(); plot.window(c(0,Xs),c(0,Ys),asp=1)
    rasterImage(IFrgb, 0, Ys-ifH, ifW, Ys, interpolate=TRUE)
    rasterImage(bfd, ifW+GAP2_UM, Ys-bH, ifW+GAP2_UM+bW, Ys, interpolate=TRUE)
    rect(ifW+GAP2_UM, Ys-bH, ifW+GAP2_UM+bW, Ys, border="red", lwd=1)
    # MSI outline in BF subregion
    MP<-msi_bf_points(A); if(!is.null(MP)){ ib<-MP[,1]>=b[1]&MP[,1]<=b[3]&MP[,2]>=b[2]&MP[,2]<=b[4]
      if(any(ib)){ Wn<-max(3L,ceiling((b[3]-b[1]+1)/8)); Hn<-max(3L,ceiling((b[4]-b[2]+1)/8)); M<-matrix(0,Wn,Hn)
        xi<-pmin(pmax(round((MP[ib,1]-b[1])/8)+1,1),Wn); yi<-pmin(pmax(round((MP[ib,2]-b[2])/8)+1,1),Hn); M[cbind(xi,yi)]<-1
        M<-fillHull(closing(M>0, makeBrush(5,"disc"))); for(cc in ocontour(bwlabel(M))){ if(nrow(cc)<4) next
          lines(c(ifW+GAP2_UM+cc[,1]*8*BF_PX_UM, ifW+GAP2_UM+cc[1,1]*8*BF_PX_UM), c(Ys-cc[,2]*8*BF_PX_UM, Ys-cc[1,2]*8*BF_PX_UM), col="#39ff14", lwd=2) } } }
    sb<-100; segments(0,Ys*0.02,sb,Ys*0.02,lwd=4); text(sb/2,Ys*0.05,"100 um",cex=0.8)
    title(sprintf("%s sec %d subregion  |  IF: DAPI(cyan)+ZO-1(red), NIS LUT (30x native)   |   brightfield + MSI outline (green)", A, n), cex.main=1.0)
    cat(sprintf("[97] %s sec%d: IF %dx%d px, BF %dx%d px\n", A, n, d[3]-d[1]+1, d[4]-d[2]+1, b[3]-b[1]+1, b[4]-b[2]+1))
  }
}
dev.off(); cat("[97] DONE -> roi_pairs_subregions.pdf\n")
