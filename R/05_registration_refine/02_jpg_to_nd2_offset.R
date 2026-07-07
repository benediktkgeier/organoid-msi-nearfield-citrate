#!/usr/bin/env Rscript
# ============================================================================
# 02_jpg_to_nd2_offset.R - STEP 3: transfer refined MSI->JPG to native .nd2.
#   The JPG is the whole-slide .nd2 downscaled (x scale = nd2_W/jpg_W) and
#   VERTICALLY flipped. Scale is exact from dims; the only free parameter is the
#   global (tx,ty) OFFSET. Solve it ONCE per slide by maximising the TOTAL
#   tissue-on-organoid darkness contrast of all 20 refined masks against the
#   .nd2 (block-mean thumbnail from R/03_coarse_registration/02_jpg_to_nd2.R). Then render the slide-level .nd2
#   overlay for validation.
#
# Input : cache/register/jpgxform_<sid>.rds (R/05_registration_refine/01_refine_jpg.R), cache/register/nd2thumb_<slide>.rds (R/03_coarse_registration/02_jpg_to_nd2.R)
# Output: cache/register/jpg2nd2off_<slide>.rds (scale, jpg_H, tx, ty)
#         figures/registration/slide_nd2_<slide>.pdf
# Usage : Rscript R/05_registration_refine/02_jpg_to_nd2_offset.R
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
suppressPackageStartupMessages({ library(Cardinal); library(jpeg) })
MSI_DIR <- file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"MSI")
REG_CACHE <- file.path(CACHE_DIR,"register"); REG_FIG <- file.path(FIG_DIR,"registration")
slide_of <- function(sid) if (grepl("sl6A",sid)) "sl6A" else "sl4A"
slide_jpg <- function(sl) if (sl=="sl6A") file.path(MSI_DIR,"06102026_AO_0h_sl6A_small.jpg") else file.path(MSI_DIR,"06102026_AO_20h_sl4A_small.jpg")
sec_code <- function(sid) sub(".*_sec","",sid)

mse <- readRDS(TISSUE_MSE); pd <- as.data.frame(pixelData(mse))
SIDS <- levels(pixelData(mse)$sample_id); if (is.null(SIDS)) SIDS <- sort(unique(as.character(pd$sample_id)))

for (sl in c("sl6A","sl4A")) {
  sids <- SIDS[vapply(SIDS, slide_of, "")==sl]
  nb <- readRDS(file.path(REG_CACHE, sprintf("nd2thumb_%s.rds", sl)))   # list(m[Hn x Wn], F, SX, SY)
  Ft <- nb$F; Dn <- 1 - nb$m/max(nb$m); Hn <- nrow(Dn); Wn <- ncol(Dn)
  j <- jpeg::readJPEG(slide_jpg(sl)); jpg_H <- nrow(j); jpg_W <- ncol(j)
  scale <- nb$SX / jpg_W                                                # nd2 px per jpg px (~2.895)

  # per-section refined tissue + bg in JPG coords
  T <- list()
  for (sid in sids) { R <- readRDS(file.path(REG_CACHE, sprintf("jpgxform_%s.rds", sid)))
    sub <- pd[as.character(pd$sample_id)==sid, c("x","y","is_tissue")]; tis<-sub[sub$is_tissue,]; bg<-sub[!sub$is_tissue,]
    if (nrow(bg)>300) bg<-bg[sample(nrow(bg),300),]
    T[[sid]] <- list(jt=apply_affine(R$B_msi_jpg, cbind(tis$x,tis$y)), jb=apply_affine(R$B_msi_jpg, cbind(bg$x,bg$y)))
  }
  # JPG -> nd2 thumb pixel given offset (tx,ty) in nd2 px:  V-flip y, scale, +offset, /F
  samp <- function(jp, tx, ty){ tcol <- (scale*jp[,1] + tx)/Ft; trow <- (scale*(jpg_H - jp[,2]) + ty)/Ft
    c<-pmin(pmax(round(tcol),1),Wn); r<-pmin(pmax(round(trow),1),Hn); Dn[cbind(r,c)] }
  total_contrast <- function(tx,ty){ s<-0; for (sid in sids){ s<-s+ mean(samp(T[[sid]]$jt,tx,ty))-mean(samp(T[[sid]]$jb,tx,ty)) }; s }

  # coarse then fine offset search (nd2 px)
  best <- list(sc=-Inf,tx=0,ty=0)
  for (tx in seq(-600,600,Ft)) for (ty in seq(-600,600,Ft)){ sc<-total_contrast(tx,ty); if(sc>best$sc) best<-list(sc=sc,tx=tx,ty=ty) }
  for (tx in best$tx+seq(-Ft,Ft,4)) for (ty in best$ty+seq(-Ft,Ft,4)){ sc<-total_contrast(tx,ty); if(sc>best$sc) best<-list(sc=sc,tx=tx,ty=ty) }
  cat(sprintf("[33] %s: scale=%.4f offset nd2=(%d,%d) total_contrast=%.3f (mean %.3f/section)\n",
              sl, scale, best$tx, best$ty, best$sc, best$sc/length(sids)))
  saveRDS(list(slide=sl, scale=scale, jpg_H=jpg_H, jpg_W=jpg_W, tx=best$tx, ty=best$ty, nd2_W=nb$SX, nd2_H=nb$SY),
          file.path(REG_CACHE, sprintf("jpg2nd2off_%s.rds", sl)))

  # ---- slide-level nd2 overlay (thumbnail) ----
  pdf(file.path(REG_FIG, sprintf("slide_nd2_%s.pdf", sl)), width=16, height=9)
  par(mar=c(1,1,3,1)); plot.new(); plot.window(c(1,Wn),c(Hn,1), asp=1)
  rasterImage(as.raster(1 - Dn), 1, Hn, Wn, 1)              # show BF (light bg, dark organoids)
  pal <- grDevices::hcl.colors(length(sids),"Dark 3"); k<-0
  for (sid in sids) { k<-k+1
    jt <- T[[sid]]$jt; tcol <- (scale*jt[,1]+best$tx)/Ft; trow <- (scale*(jpg_H-jt[,2])+best$ty)/Ft
    points(tcol, trow, pch=15, cex=0.25, col=adjustcolor(pal[k],0.6))
    text(mean(tcol), mean(trow), sec_code(sid), col="red", font=2, cex=0.9)
  }
  title(sprintf("%s - refined MSI tissue masks on NATIVE .nd2 (block-mean thumb), offset=(%d,%d)", sl, best$tx, best$ty))
  dev.off(); cat(sprintf("[33] wrote slide_nd2_%s.pdf\n", sl))
}
cat("[33] DONE\n")
