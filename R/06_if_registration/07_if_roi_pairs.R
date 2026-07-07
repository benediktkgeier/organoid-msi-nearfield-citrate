#!/usr/bin/env Rscript
# ============================================================================
# 07_if_roi_pairs.R - side-by-side, SAME-SCALE report of, per organoid section:
#     LEFT  = high-res IF organoid section, 4-CHANNEL COMPOSITE (B-slide, 30x)
#             DAPI=cyan, b-catenin(GFP)=green, ZO-1(mCherry)=red, F-actin(Cy5)=gray
#     RIGHT = the brightfield ROI the user selected (rectangle on the A-slide),
#             cropped from the native brightfield .nd2 and shown like
#             registration_native (om/max(om), no high-pass / no extreme contrast).
#
#   Both panels are drawn at the SAME physical scale (um/px isotropic, shared
#   extent across every panel) so sizes are directly comparable even though the
#   ROIs vary in size. Rectangles from results/if_registration/bf_rectangles.csv.
#
# Input : bf_rectangles.csv, hi-res .nd2 (cached hr4_<sid>.rds), BF .nd2, nd2thumb
# Output: figures/if_registration/roi_pairs.pdf
# Usage : Rscript R/06_if_registration/07_if_roi_pairs.R          (builds hr4 cache on first run)
# ============================================================================

Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jre1.8.0_491")
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/if_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register_if.R"))
suppressPackageStartupMessages({ library(RBioFormats); library(EBImage); library(Cardinal) })

rect    <- read.csv(file.path(IF_RES, "bf_rectangles.csv"), stringsAsFactors = FALSE)

# registered MSI is_tissue footprints (measured organoids) projected to BF native px
.mse <- readRDS(TISSUE_MSE); .pd <- as.data.frame(pixelData(.mse)); .msi_cache <- list()
msi_bf_points <- function(A) {
  if (!is.null(.msi_cache[[A]])) return(.msi_cache[[A]])
  sids <- grep(sprintf("_%s_", A), unique(as.character(.pd$sample_id)), value = TRUE)
  P <- do.call(rbind, lapply(sids, function(s) {
    xf <- file.path(REG_CACHE, sprintf("nd2final_%s.rds", s)); if (!file.exists(xf)) return(NULL)
    B <- readRDS(xf)$B_msi_nd2; sub <- .pd[as.character(.pd$sample_id)==s & .pd$is_tissue, c("x","y")]
    if (!nrow(sub)) return(NULL); apply_affine(B, cbind(sub$x, sub$y))
  }))
  .msi_cache[[A]] <<- P; P
}
slide_b <- function(A) if (A == "sl6A") "sl6b" else "sl4b"
grp_of  <- function(A) if (A == "sl6A") "0h"   else "20h"
HQF      <- 2L                    # high-quality DAPI cache factor (0.72 um/px)
HR_PX_UM <- HR_UMPX * HQF         # hi-res DAPI display scale
BF_PX_UM <- BF_UMPX               # brightfield native ~1.833 um/px
GAP_UM   <- 400; SB_UM <- 1000    # gap between panels; scale bar length
# native per-channel display: linear scale between low/high percentiles of the
# RAW counts (no saturating clip -> not blown out), matching the acquisition look.
nat <- function(ch, lo = 0.002, hi = 0.999) {
  q <- quantile(ch[is.finite(ch)], c(lo, hi), na.rm = TRUE)
  if (!is.finite(q[2]) || q[2] <= q[1]) q <- range(ch[is.finite(ch)])
  pmin(pmax((ch - q[1]) / (q[2] - q[1]), 0), 1)
}
ds_disp  <- function(m, target = 2200) { f <- max(1, round(max(dim(m))/target)); if (f > 1) ds_mean(m, f) else m }

# hi-res DAPI display (single channel) from RAW counts, native linear scaling.
dapi_disp <- function(m) ds_disp(nat(m))

# draw the OUTLINE of the registered MSI footprint that falls within a crop,
# into the brightfield panel placed with top-left at (x0_um, ytop_um).
draw_msi_outline <- function(MP, cx0, cy0, cx1, cy1, x0_um, ytop_um, bin = 10L, col = "#39ff14") {
  if (is.null(MP)) return(invisible())
  ib <- MP[,1]>=cx0 & MP[,1]<=cx1 & MP[,2]>=cy0 & MP[,2]<=cy1; if (!any(ib)) return(invisible())
  cw <- cx1-cx0+1; ch <- cy1-cy0+1; Wn <- max(3L, ceiling(cw/bin)); Hn <- max(3L, ceiling(ch/bin))
  M <- matrix(0, Wn, Hn)   # EBImage orientation [x, y]
  xi <- pmin(pmax(round((MP[ib,1]-cx0)/bin)+1, 1), Wn); yi <- pmin(pmax(round((MP[ib,2]-cy0)/bin)+1, 1), Hn)
  M[cbind(xi, yi)] <- 1
  M <- EBImage::fillHull(EBImage::closing(M > 0, EBImage::makeBrush(5, "disc")))
  oc <- EBImage::ocontour(EBImage::bwlabel(M))
  for (cc in oc) { if (nrow(cc) < 4) next
    px <- x0_um + (cc[,1]*bin)*BF_PX_UM; py <- ytop_um - (cc[,2]*bin)*BF_PX_UM
    lines(c(px, px[1]), c(py, py[1]), col = col, lwd = 1.8) }
}

# apply a NIS-Elements LUT (lo,hi,g) to one raw channel -> [0,1]
lut_apply <- function(ch, L) pmin(pmax((ch - L["lo"]) / (L["hi"] - L["lo"]), 0), 1)^(1/L["g"])
# 4-channel composite from RAW counts using the slide's NIS LUTs + channel colours
lut_composite <- function(ALL, LUT) {
  H <- dim(ALL)[1]; W <- dim(ALL)[2]; RGB <- array(0, c(H, W, 3))
  for (k in seq_len(min(dim(ALL)[3], length(LUT)))) {
    d <- lut_apply(ALL[,,k], LUT[[k]]); col <- IF_LUT_COL[[k]]
    for (j in 1:3) if (col[j] > 0) RGB[,,j] <- RGB[,,j] + d * col[j]
  }
  RGB <- pmin(RGB, 1)
  array(c(ds_disp(RGB[,,1]), ds_disp(RGB[,,2]), ds_disp(RGB[,,3])), c(dim(ds_disp(RGB[,,1])), 3))
}

# BF native px dims per A-slide
bfdim <- list(); for (A in unique(rect$slide)) { cm <- coreMetadata(read.metadata(IF_SLIDES$msi_bf[IF_SLIDES$msi_slide==A]), series=1); bfdim[[A]] <- c(SX=cm$sizeX, SY=cm$sizeY) }

# ---- pass 1: build hr4 cache (slow) + gather BF crop sizes for global extent ----
hr_ref_dim <- NULL
for (i in seq_len(nrow(rect))) {
  A <- rect$slide[i]; sid <- sprintf("AO_%s_%s_sec%d", grp_of(A), slide_b(A), rect$section[i])
  hq <- file.path(IF_CACHE, sprintf("hrdapihq_%s.rds", sid))  # high-quality RAW DAPI (F=2)
  hp <- file.path(IF_DIR, paste0(IF_SLIDES$hr_prefix[IF_SLIDES$slide==slide_b(A)], rect$section[i], ".nd2"))
  if (!file.exists(hq)) { cat(sprintf("[96] building HQ DAPI %s ...\n", sid)); saveRDS(list(m=nd2_channel_block_mean(hp, ch=DAPI_CH_DEFAULT, F=HQF, normalize=FALSE)$m, F=HQF), hq) }
  if (is.null(hr_ref_dim)) hr_ref_dim <- dim(readRDS(hq)$m)
}
hr_W <- hr_ref_dim[2]*HR_PX_UM; hr_H <- hr_ref_dim[1]*HR_PX_UM
rect$bfW <- (rect$fx1-rect$fx0)*sapply(rect$slide,function(A)bfdim[[A]]["SX"])*BF_PX_UM
rect$bfH <- (rect$fy1-rect$fy0)*sapply(rect$slide,function(A)bfdim[[A]]["SY"])*BF_PX_UM
XSPAN <- hr_W + GAP_UM + max(rect$bfW); YSPAN <- max(hr_H, max(rect$bfH))

# ---- pass 2: draw (same scale everywhere) ----
pdf(file.path(IF_FIG, Sys.getenv("ROI_OUT", "roi_pairs.pdf")), width = 11, height = 14)
for (A in c("sl6A","sl4A")) {
  rs <- rect[rect$slide==A,]; rs <- rs[order(rs$section),]; SX <- bfdim[[A]]["SX"]; SY <- bfdim[[A]]["SY"]
  layout(matrix(1:6, nrow=6, byrow=TRUE)); par(mar=c(1,1,2.2,1))
  for (i in seq_len(nrow(rs))) {
    n <- rs$section[i]; slb <- slide_b(A); sid <- sprintf("AO_%s_%s_sec%d", grp_of(A), slb, n)
    comp <- dapi_disp(readRDS(file.path(IF_CACHE, sprintf("hrdapihq_%s.rds", sid)))$m); hr_lab <- "IF DAPI (HQ)"
    cx0<-max(1,round(rs$fx0[i]*SX)); cx1<-min(SX,round(rs$fx1[i]*SX)); cy0<-max(1,round(rs$fy0[i]*SY)); cy1<-min(SY,round(rs$fy1[i]*SY))
    img <- read.image(IF_SLIDES$msi_bf[IF_SLIDES$msi_slide==A], series=1, subset=list(x=cx0:cx1,y=cy0:cy1), normalize=TRUE)
    a<-as.array(img); om<-if(length(dim(a))==2)a else a[,,1]; om<-t(om); bfd<-om/max(om)
    bw <- (cx1-cx0+1)*BF_PX_UM; bh <- (cy1-cy0+1)*BF_PX_UM
    plot.new(); plot.window(c(0, XSPAN), c(0, YSPAN), asp=1)
    rasterImage(comp, 0, YSPAN-hr_H, hr_W, YSPAN, interpolate=TRUE)
    rasterImage(bfd,  hr_W+GAP_UM, YSPAN-bh, hr_W+GAP_UM+bw, YSPAN, interpolate=TRUE)
    rect(hr_W+GAP_UM, YSPAN-bh, hr_W+GAP_UM+bw, YSPAN, border="red", lwd=1.5)
    # outline the registered MSI is_tissue footprint (measured organoid) within this crop
    draw_msi_outline(msi_bf_points(A), cx0, cy0, cx1, cy1, hr_W+GAP_UM, YSPAN)
    segments(0, 60, SB_UM, 60, lwd=4); text(SB_UM/2, 150, "1 mm", cex=0.8)
    title(sprintf("sec %d   |   hi-res %s (%s, 30x)   vs   brightfield ROI (%s, native)", n, hr_lab, slb, A), cex.main=0.95, adj=0)
  }
  mtext(sprintf("%s organoid sections - same scale (1 mm) | green outline = registered MSI-measured organoid (is_tissue)", A), outer=TRUE, line=-1.4, font=2, cex=0.85)
}
dev.off()
cat("[96] DONE -> roi_pairs.pdf\n")
