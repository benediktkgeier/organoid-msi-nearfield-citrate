#!/usr/bin/env Rscript
# ============================================================================
# 01_refine_jpg.R - STEP 2: within-box refinement on the slide JPG.
#   Slide-level teach placement is correct (user-confirmed: boxes on the right
#   organoids). Refine each section's is_tissue mask onto the JPG organoid with
#   a SMALL local search - mainly translation (+/-12 MSI px), minor rotation
#   (+/-4 deg) and scale (+/-3%) only if it helps - scored by tissue-on-organoid
#   darkness contrast. No flips, no large moves (foundation already right).
#
# Input : cache/register/teach_<sid>.rds (coarse MSI->jpg from R/03_coarse_registration/01_teach_msi_to_jpg.R)
# Output: cache/register/jpgxform_<sid>.rds (refined MSI->jpg + score)
#         figures/registration/refine_jpg_<slide>.pdf (coarse red vs refined green)
#         results/registration/refine_jpg_summary.csv
# Usage : Rscript R/05_registration_refine/01_refine_jpg.R
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
suppressPackageStartupMessages({ library(Cardinal); library(jpeg) })
MSI_DIR <- file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"MSI")
REG_CACHE <- file.path(CACHE_DIR,"register"); REG_RES <- file.path(RES_DIR,"registration"); REG_FIG <- file.path(FIG_DIR,"registration")
slide_of <- function(sid) if (grepl("sl6A",sid)) "sl6A" else "sl4A"
slide_jpg <- function(sl) if (sl=="sl6A") file.path(MSI_DIR,"06102026_AO_0h_sl6A_small.jpg") else file.path(MSI_DIR,"06102026_AO_20h_sl4A_small.jpg")
sec_code <- function(sid) sub(".*_sec","",sid)

TRANS_MSI <- 12; ROTS <- c(-4,-2,0,2,4)*pi/180; SCALES <- c(0.97,1.0,1.03)

mse <- readRDS(TISSUE_MSE); pd <- as.data.frame(pixelData(mse))
SIDS <- levels(pixelData(mse)$sample_id); if (is.null(SIDS)) SIDS <- sort(unique(as.character(pd$sample_id)))

jpgc <- list(); getj <- function(sl){ if(is.null(jpgc[[sl]])){ j<-jpeg::readJPEG(slide_jpg(sl)); g<-if(length(dim(j))==3)(j[,,1]+j[,,2]+j[,,3])/3 else j; jpgc[[sl]]<<-list(g=g,H=nrow(g),W=ncol(g))}; jpgc[[sl]] }

summ <- list(); byslide <- list()
for (sid in SIDS) {
  sl <- slide_of(sid); tf <- file.path(REG_CACHE, sprintf("teach_%s.rds", sid)); if (!file.exists(tf)) next
  T <- readRDS(tf); B0 <- T$B_msi_jpg; J <- getj(sl)
  dark <- function(p){ c<-pmin(pmax(round(p[,1]),1),J$W); r<-pmin(pmax(round(p[,2]),1),J$H); 1 - J$g[cbind(r,c)] }
  sub <- pd[as.character(pd$sample_id)==sid, c("x","y","is_tissue")]
  tis <- sub[sub$is_tissue,]; bg <- sub[!sub$is_tissue,]; if (nrow(bg)>600) bg<-bg[sample(nrow(bg),600),]
  jt0 <- apply_affine(B0, cbind(tis$x,tis$y)); jb0 <- apply_affine(B0, cbind(bg$x,bg$y))
  jpm <- sqrt(abs(B0[1,1]*B0[2,2]-B0[1,2]*B0[2,1]))       # jpg px per MSI px (~1.85)
  cX <- mean(jt0[,1]); cY <- mean(jt0[,2])
  score <- function(tx,ty,s,th){ pt<-rs_transform(jt0,cX,cY,s,th,tx,ty); pb<-rs_transform(jb0,cX,cY,s,th,tx,ty); mean(dark(pt))-mean(dark(pb)) }
  # hierarchical small search (translation first, then minor rot/scale, then fine translation)
  bT <- list(sc=score(0,0,1,0),tx=0,ty=0)
  for (tx in seq(-TRANS_MSI,TRANS_MSI,1)*jpm) for (ty in seq(-TRANS_MSI,TRANS_MSI,1)*jpm){ s<-score(tx,ty,1,0); if(s>bT$sc) bT<-list(sc=s,tx=tx,ty=ty) }
  bRS <- list(sc=-Inf,s=1,th=0); for (s in SCALES) for (th in ROTS){ sc<-score(bT$tx,bT$ty,s,th); if(sc>bRS$sc) bRS<-list(sc=sc,s=s,th=th) }
  bF <- list(sc=bRS$sc,tx=bT$tx,ty=bT$ty); for (tx in bT$tx+seq(-3,3,1)*jpm) for (ty in bT$ty+seq(-3,3,1)*jpm){ sc<-score(tx,ty,bRS$s,bRS$th); if(sc>bF$sc) bF<-list(sc=sc,tx=tx,ty=ty) }
  refine <- function(msi){ j<-apply_affine(B0,msi); rs_transform(j,cX,cY,bRS$s,bRS$th,bF$tx,bF$ty) }
  Bref <- fit_affine(cbind(sub$x,sub$y), refine(cbind(sub$x,sub$y)))
  c0 <- score(0,0,1,0)
  saveRDS(list(sid=sid, slide=sl, B_msi_jpg=Bref, B_coarse=B0, trans_msi=c(bF$tx,bF$ty)/jpm,
               rot_deg=bRS$th*180/pi, scale=bRS$s, contrast_coarse=c0, contrast_refine=bF$sc, jpm=jpm),
          file.path(REG_CACHE, sprintf("jpgxform_%s.rds", sid)))
  summ[[sid]] <- data.frame(sample_id=sid, slide=sl, trans_x_msi=round(bF$tx/jpm,1), trans_y_msi=round(bF$ty/jpm,1),
    rot_deg=round(bRS$th*180/pi,1), scale=bRS$s, contrast_coarse=round(c0,4), contrast_refine=round(bF$sc,4), stringsAsFactors=FALSE)
  byslide[[sl]] <- c(byslide[[sl]], sid)
  cat(sprintf("[32] %s: trans(%.1f,%.1f)MSIpx rot=%.0f s=%.2f contrast %.3f->%.3f\n",
    sid, bF$tx/jpm, bF$ty/jpm, bRS$th*180/pi, bRS$s, c0, bF$sc))
}
summ_df <- do.call(rbind, summ); write.csv(summ_df, file.path(REG_RES,"refine_jpg_summary.csv"), row.names=FALSE)

# ---- QC per slide: per-section zoom, coarse (red) vs refined (green) ----
for (sl in names(byslide)) {
  J <- getj(sl)
  pdf(file.path(REG_FIG, sprintf("refine_jpg_%s.pdf", sl)), width=15, height=9)
  par(mfrow=c(2,3), mar=c(1,1,2.5,1))
  for (sid in byslide[[sl]]) {
    R <- readRDS(file.path(REG_CACHE, sprintf("jpgxform_%s.rds", sid)))
    sub <- pd[as.character(pd$sample_id)==sid, c("x","y","is_tissue")]; tis<-sub[sub$is_tissue,]
    j0 <- apply_affine(R$B_coarse, cbind(tis$x,tis$y)); j1 <- apply_affine(R$B_msi_jpg, cbind(tis$x,tis$y))
    allp <- rbind(j0,j1); mg <- 0.5*max(diff(range(allp[,1])),diff(range(allp[,2])))
    bx<-range(allp[,1])+c(-mg,mg); by<-range(allp[,2])+c(-mg,mg)
    cols<-max(1,floor(bx[1])):min(J$W,ceiling(bx[2])); rows<-max(1,floor(by[1])):min(J$H,ceiling(by[2]))
    plot.new(); plot.window(range(cols), rev(range(rows)), asp=1); rasterImage(J$g[rows,cols], min(cols),max(rows),max(cols),min(rows))
    points(j0[,1],j0[,2],pch=15,cex=0.5,col=adjustcolor("red",0.4))
    points(j1[,1],j1[,2],pch=15,cex=0.5,col=adjustcolor("green3",0.5))
    title(sprintf("%s  c:%.3f->%.3f  d(%.0f,%.0f)px", sec_code(sid), R$contrast_coarse, R$contrast_refine, R$trans_msi[1], R$trans_msi[2]), cex.main=0.95)
  }
  plot.new(); legend("center", legend=c("coarse (teach)","refined"), pch=15, col=c("red","green3"), bty="n", cex=1.3)
  dev.off(); cat(sprintf("[32] wrote refine_jpg_%s.pdf\n", sl))
}
cat("\n[32] DONE\n"); print(summ_df[,c("sample_id","trans_x_msi","trans_y_msi","rot_deg","scale","contrast_coarse","contrast_refine")], row.names=FALSE)
