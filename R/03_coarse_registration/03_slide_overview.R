#!/usr/bin/env Rscript
# ============================================================================
# 03_slide_overview.R - STEP 1: slide-level co-registration from Bruker teach
#   marks ONLY. For each slide, place every section's measured Area (.mis ROI,
#   which exactly equals the raster extent) onto the original slide brightfield
#   (_small.jpg) via the OriginalImageTeachPoint affine, and overlay the MSI
#   is_tissue mask inside each area so MSI organoids can be compared to the BF.
#   No .nd2, no refinement - just the teach-mark foundation, for validation.
#
# Output: figures/registration/slide_overview_<slide>.pdf
# Usage : Rscript R/03_coarse_registration/03_slide_overview.R
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
suppressPackageStartupMessages({ library(Cardinal); library(jpeg) })
MSI_DIR <- file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"MSI")
REG_FIG <- file.path(FIG_DIR, "registration"); dir.create(REG_FIG, showWarnings = FALSE, recursive = TRUE)

slide_of <- function(sid) if (grepl("sl6A", sid)) "sl6A" else "sl4A"
slide_folder <- function(sl) if (sl=="sl6A") file.path(MSI_DIR,"06102026_AO_0h_sl6A") else file.path(MSI_DIR,"06102026_AO_20h_sl4A")
slide_jpg <- function(sl) if (sl=="sl6A") file.path(MSI_DIR,"06102026_AO_0h_sl6A_small.jpg") else file.path(MSI_DIR,"06102026_AO_20h_sl4A_small.jpg")
sec_code <- function(sid) sub(".*_sec","",sid)
find_mis <- function(folder, code) { f <- list.files(folder, pattern=sprintf("(_sec%s|_%s)\\.mis$", code, code), full.names=TRUE); f<-f[!grepl("\\.bak$",f)]; if(!length(f)) NA_character_ else f[which.min(nchar(basename(f)))] }

mse <- readRDS(TISSUE_MSE); pd <- as.data.frame(pixelData(mse))
SIDS <- levels(pixelData(mse)$sample_id); if (is.null(SIDS)) SIDS <- sort(unique(as.character(pd$sample_id)))

# map MSI(x,y) linearly into the flex Area bbox for orientation o, then ->jpg
place_jpg <- function(sub, bx, by, W, H, B_fo, o) {
  fx <- if (o[1]>0) bx[1]+(sub$x-1)/(W-1)*(bx[2]-bx[1]) else bx[2]-(sub$x-1)/(W-1)*(bx[2]-bx[1])
  fy <- if (o[2]>0) by[1]+(sub$y-1)/(H-1)*(by[2]-by[1]) else by[2]-(sub$y-1)/(H-1)*(by[2]-by[1])
  apply_affine(B_fo, cbind(fx, fy))
}
CAND <- list(c(1,1),c(1,-1),c(-1,1),c(-1,-1))

for (sl in c("sl6A","sl4A")) {
  sids <- SIDS[vapply(SIDS, slide_of, "")==sl]
  j <- jpeg::readJPEG(slide_jpg(sl)); JH <- nrow(j); JW <- ncol(j)
  jg <- if (length(dim(j))==3) (j[,,1]+j[,,2]+j[,,3])/3 else j
  dark_at <- function(ox,oy){ ox<-pmin(pmax(round(ox),1),JW); oy<-pmin(pmax(round(oy),1),JH); 1-jg[cbind(oy,ox)] }

  # collect per-section area + transform; choose ONE orientation per slide
  info <- list(); cands <- matrix(0, length(sids), 4)
  for (i in seq_along(sids)) {
    sid <- sids[i]; mis <- find_mis(slide_folder(sl), sec_code(sid)); if (is.na(mis)) next
    M <- parse_mis(mis); if (nrow(M$orig)<3 || nrow(M$area)<2) next
    B_fo <- fit_affine(M$orig[,1:2], M$orig[,3:4]); bx<-range(M$area[,1]); by<-range(M$area[,2])
    sub <- pd[as.character(pd$sample_id)==sid, c("x","y","is_tissue")]; W<-max(sub$x); H<-max(sub$y)
    for (k in seq_along(CAND)) { oc<-place_jpg(sub,bx,by,W,H,B_fo,CAND[[k]]); d<-dark_at(oc[,1],oc[,2]); cands[i,k]<-mean(d[sub$is_tissue])-mean(d[!sub$is_tissue]) }
    info[[sid]] <- list(sid=sid, B_fo=B_fo, bx=bx, by=by, sub=sub, W=W, H=H, area=M$area)
  }
  o <- CAND[[ which.max(colSums(cands)) ]]
  cat(sprintf("[31] %s global orientation (%+d,%+d)\n", sl, o[1], o[2]))

  pdf(file.path(REG_FIG, sprintf("slide_overview_%s.pdf", sl)), width=16, height=9)
  # ---- Page 1: whole slide + all measured areas ----
  par(mar=c(1,1,3,1)); plot.new(); plot.window(c(1,JW), c(JH,1), asp=1)
  rasterImage(jg, 1, JH, JW, 1)
  pal <- grDevices::hcl.colors(length(info), "Dark 3"); k<-0
  for (sid in names(info)) { k<-k+1; I<-info[[sid]]
    corners <- place_jpg(data.frame(x=c(1,I$W,I$W,1), y=c(1,1,I$H,I$H)), I$bx,I$by,I$W,I$H,I$B_fo,o)
    polygon(corners[,1], corners[,2], border="yellow", lwd=2)
    tn <- place_jpg(I$sub[I$sub$is_tissue,], I$bx,I$by,I$W,I$H,I$B_fo,o)
    points(tn[,1], tn[,2], pch=15, cex=0.12, col=adjustcolor(pal[k],0.55))
    cc <- colMeans(corners); text(cc[1], cc[2], sec_code(sid), col="red", font=2, cex=1.1)
  }
  title(sprintf("%s - slide-level co-registration from Bruker teach marks (measured areas + MSI tissue mask on slide BF)", sl))

  # ---- Per-section zoom pages (jpg crop + area + tissue mask) ----
  par(mfrow=c(2,3), mar=c(1,1,2.5,1))
  for (sid in names(info)) { I<-info[[sid]]
    corners <- place_jpg(data.frame(x=c(1,I$W,I$W,1), y=c(1,1,I$H,I$H)), I$bx,I$by,I$W,I$H,I$B_fo,o)
    mg <- 0.6*max(diff(range(corners[,1])), diff(range(corners[,2])))
    bx2 <- range(corners[,1])+c(-mg,mg); by2 <- range(corners[,2])+c(-mg,mg)
    cols <- max(1,floor(bx2[1])):min(JW,ceiling(bx2[2])); rows <- max(1,floor(by2[1])):min(JH,ceiling(by2[2]))
    plot.new(); plot.window(range(cols), rev(range(rows)), asp=1)
    rasterImage(jg[rows,cols], min(cols), max(rows), max(cols), min(rows))
    polygon(corners[,1], corners[,2], border="yellow", lwd=2)
    tn <- place_jpg(I$sub[I$sub$is_tissue,], I$bx,I$by,I$W,I$H,I$B_fo,o)
    points(tn[,1], tn[,2], pch=15, cex=0.5, col=adjustcolor("red",0.5))
    title(sprintf("%s (%dx%d px)", sec_code(sid), I$W, I$H), cex.main=1)
  }
  dev.off(); cat(sprintf("[31] wrote slide_overview_%s.pdf\n", sl))
}
cat("[31] DONE\n")
