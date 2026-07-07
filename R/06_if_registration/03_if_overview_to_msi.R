#!/usr/bin/env Rscript
# ============================================================================
# 03_if_overview_to_msi.R - STEP 2: fit the overview -> MSI-brightfield (.nd2)
#   transform per slide from CORRESPONDENCES
#       {hi-res section center in overview (R/06_if_registration/02_if_overview_to_bf.R)}  <->  {MSI section-N footprint
#        centroid in BF (existing nd2final + is_tissue)}
#   then compose the final per-section hi-res -> BF affine B_hr_bf, with an
#   optional small translation polish against the MSI footprint. Reports the fit.
#
#   Serial B-section vs A-slide => the fit is an APPROXIMATE rigid+scale bridge;
#   residuals are reported. sec_n <-> MSI section n (user). Sections with no MSI
#   counterpart (e.g. 20h sec6) are carried but excluded from the fit.
#
# Input : cache/register_if/locate_<sid_if>.rds, cache/register/nd2final_*.rds, TISSUE_MSE
# Output: cache/register_if/ov_to_bf_<slide>.rds, hr_to_bf_<sid_if>.rds
#         figures/if_registration/step2_overview_to_msi.pdf, results/if_registration/fit_summary.csv
# Usage : Rscript R/06_if_registration/03_if_overview_to_msi.R [all | slide ...]   (default = all slides)
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/if_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register_if.R"))
suppressPackageStartupMessages({ library(Cardinal); library(EBImage) })

args   <- commandArgs(trailingOnly = TRUE)
SEC    <- if_sections()
slides <- if (length(args) == 0 || identical(args, "all")) IF_SLIDES$slide else args

mse <- readRDS(TISSUE_MSE); pd <- as.data.frame(pixelData(mse))
# MSI section footprint points (BF native) + centroid, union of secNa/b
msi_section_pts <- function(msi_slide, secn) {
  sids <- msi_sids_for_section(msi_slide, secn); if (!length(sids)) return(NULL)
  pts <- do.call(rbind, lapply(sids, function(sid) {
    B <- readRDS(file.path(REG_CACHE, sprintf("nd2final_%s.rds", sid)))$B_msi_nd2
    sub <- pd[as.character(pd$sample_id) == sid & pd$is_tissue, c("x","y")]
    apply_affine(B, cbind(sub$x, sub$y))
  }))
  list(pts = pts, centroid = colMeans(pts), sids = sids)
}

# hi-res DAPI organoid centroid in hi-res NATIVE px: peak of blurred DAPI density
# (the compact organoid core), NOT the loose threshold (which fills the whole FOV).
hr_dapi_centroid <- function(h, blur = 51L, qcore = 0.985) {
  v  <- norm01(h$m5)
  bl <- EBImage::filter2(v, EBImage::makeBrush(blur, "disc")/sum(EBImage::makeBrush(blur, "disc")), boundary = "replicate")
  core <- which(bl > quantile(bl, qcore), arr.ind = TRUE)              # densest core region
  lab  <- EBImage::bwlabel(bl > quantile(bl, qcore)); tb <- table(lab[lab > 0])
  big  <- as.integer(names(which.max(tb))); core <- which(lab == big, arr.ind = TRUE)
  list(centroid = c(mean((core[,2]-0.5)*h$F5), mean((core[,1]-0.5)*h$F5)),
       pts = cbind((core[,2]-0.5)*h$F5, (core[,1]-0.5)*h$F5))           # organoid-core footprint
}

fit_rows <- list(); pages <- list()
for (sl in slides) {
  S <- IF_SLIDES[IF_SLIDES$slide == sl, ]; secs <- SEC[SEC$slide == sl, ]
  rec <- list(); used <- c()
  for (j in seq_len(nrow(secs))) {
    sid_if <- secs$sid_if[j]; lf <- file.path(IF_CACHE, sprintf("locate_%s.rds", sid_if))
    if (!file.exists(lf)) next
    L <- readRDS(lf); M <- msi_section_pts(S$msi_slide, secs$secn[j]); if (is.null(M)) next
    h <- readRDS(file.path(IF_CACHE, sprintf("hrthumb_%s.rds", sid_if))); HC <- hr_dapi_centroid(h)
    ov_c <- apply_affine(L$B_hr_ov, matrix(HC$centroid, 1))[1, ]   # IF organoid centroid in overview px
    rec[[sid_if]] <- list(L = L, M = M, hc = HC$centroid, hr_pts = HC$pts, ov_c = ov_c, secn = secs$secn[j])
    used <- c(used, sid_if)
  }
  if (length(used) < 3) { cat(sprintf("[92] %s: only %d correspondences - skip\n", sl, length(used))); next }
  OV <- do.call(rbind, lapply(rec[used], `[[`, "ov_c"))
  BF <- do.call(rbind, lapply(rec[used], function(z) z$M$centroid))
  B_ov_bf <- fit_affine(OV, BF)                                    # gives rotation+scale bridge
  res <- sqrt(rowSums((apply_affine(B_ov_bf, OV) - BF)^2))
  saveRDS(list(slide = sl, msi_slide = S$msi_slide, B_ov_bf = B_ov_bf, used = used,
               rmse = sqrt(mean(res^2)), res = res), file.path(IF_CACHE, sprintf("ov_to_bf_%s.rds", sl)))
  cat(sprintf("[92] %s: fit on %d sections, RMSE=%.1f BF px (%.1f um)\n", sl, length(used),
              sqrt(mean(res^2)), sqrt(mean(res^2)) * BF_UMPX))

  # per-section: compose then SNAP translation so IF organoid centroid == MSI centroid
  page_outlines <- list()
  for (sid_if in used) {
    z <- rec[[sid_if]]; B_hr_bf <- compose_affine(z$L$B_hr_ov, B_ov_bf)
    cur  <- apply_affine(B_hr_bf, matrix(z$hc, 1))[1, ]
    snap <- z$M$centroid - cur; B_hr_bf[3, ] <- B_hr_bf[3, ] + snap
    saveRDS(list(sid_if = sid_if, slide = sl, msi_slide = S$msi_slide, B_hr_bf = B_hr_bf,
                 B_hr_ov = z$L$B_hr_ov, snap = snap, SX_hr = z$L$SX_hr, SY_hr = z$L$SY_hr),
            file.path(IF_CACHE, sprintf("hr_to_bf_%s.rds", sid_if)))
    corners <- cbind(c(1,z$L$SX_hr,z$L$SX_hr,1), c(1,1,z$L$SY_hr,z$L$SY_hr))
    page_outlines[[sid_if]] <- list(corner = apply_affine(B_hr_bf, corners),
                                    iffoot = apply_affine(B_hr_bf, z$hr_pts),
                                    foot = z$M$pts, secn = z$secn)
    fit_rows[[sid_if]] <- data.frame(sid_if = sid_if, slide = sl,
                                     snap_x_um = round(snap[1]*BF_UMPX), snap_y_um = round(snap[2]*BF_UMPX),
                                     fit_res_um = round(res[match(sid_if, used)]*BF_UMPX))
  }
  pages[[sl]] <- list(B_ov_bf=B_ov_bf, res=res, used=used, outlines=page_outlines, msi_slide=S$msi_slide)
}
if (length(fit_rows)) write.csv(do.call(rbind, fit_rows), file.path(IF_RES, "fit_summary.csv"), row.names = FALSE)

# ---- report: BF frame with MSI footprints + projected IF section frames ----
pdf(file.path(IF_FIG, "step2_overview_to_msi.pdf"), width = 13, height = 6)
for (sl in names(pages)) {
  G <- pages[[sl]]
  allpts <- do.call(rbind, lapply(G$outlines, function(o) rbind(o$foot, o$corner)))
  rx <- range(allpts[,1]); ry <- range(allpts[,2])
  par(mar = c(1,1,3,1))
  plot.new(); plot.window(rx, rev(ry), asp = 1)
  pal <- rainbow(length(G$outlines)); k <- 1
  for (sid_if in names(G$outlines)) {
    o <- G$outlines[[sid_if]]
    points(o$foot[,1], o$foot[,2], pch=15, cex=0.4, col=adjustcolor("cyan",0.6))          # MSI tissue
    points(o$iffoot[,1], o$iffoot[,2], pch=15, cex=0.25, col=adjustcolor(pal[k],0.35))     # IF DAPI tissue
    polygon(o$corner[,1], o$corner[,2], border=pal[k], lwd=1.5)
    text(mean(o$corner[,1]), mean(o$corner[,2]), sprintf("sec%d", o$secn), col=pal[k], font=2, cex=0.9)
    k <- k+1
  }
  title(sprintf("STEP2  %s->%s  MSI tissue (cyan) vs IF DAPI tissue (color, centroid-snapped)   fit RMSE=%.0f um",
                sl, G$msi_slide, sqrt(mean(G$res^2))*BF_UMPX), cex.main=0.9)
}
dev.off()
cat(sprintf("[92] DONE -> step2_overview_to_msi.pdf (%d slides)\n", length(pages)))
