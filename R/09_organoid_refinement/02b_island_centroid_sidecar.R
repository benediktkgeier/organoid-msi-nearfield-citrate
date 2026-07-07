#!/usr/bin/env Rscript
# ============================================================================
# 02b_island_centroid_sidecar.R - regenerate the island-cleanup sidecar with
#   each organoid's CENTROID position in PDF points (in addition to the numeric
#   label position). Replays the EXACT, deterministic render geometry of
#   R/09_organoid_refinement/02_island_cleanup_canvas.R to a
#   throwaway PDF so the centroid PDF coords match the already-annotated
#   figures/annotation/organoid_island_cleanup.pdf 1:1. The user's annotated PDF
#   is NOT touched.
#
#   Why: apical comments (in/out/mixed) are placed ON the organoid body, but the
#   original sidecar only stores the numeric label position (placed OUTSIDE the
#   organoid via a leader line). Matching comments to the label position misses
#   ~70% of them; matching to the centroid is the correct, robust target.
#
# Out: cache/organoid_island_label_positions_centroid.rds
# Usage: Rscript R/09_organoid_refinement/02b_island_centroid_sidecar.R
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
suppressPackageStartupMessages({ library(Cardinal); library(png); library(EBImage) })

REG_CACHE <- file.path(CACHE_DIR, "register")
CROP_DIR  <- file.path(FIG_DIR, "registration", "crops")
PAGE_W_IN <- 14; PAGE_H_IN <- 10
PAGE_W_PT <- PAGE_W_IN * 72; PAGE_H_PT <- PAGE_H_IN * 72

# identical helpers to R/09_organoid_refinement/02_island_cleanup_canvas.R (geometry-relevant ones) -------------------------
instance_outlines <- function(M, B, cx0, cy0) {
  M <- EBImage::fillHull(EBImage::closing(M > 0, EBImage::makeBrush(3, "disc")))
  out <- list()
  for (co in EBImage::ocontour(EBImage::bwlabel(M))) {
    if (nrow(co) < 6) next
    p <- apply_affine(B, cbind(co[, 1] + 0.5, co[, 2] + 0.5))
    out[[length(out) + 1L]] <- list(x = c(p[, 1] - cx0 + 1, p[1, 1] - cx0 + 1),
                                    y = c(p[, 2] - cy0 + 1, p[1, 2] - cy0 + 1))
  }
  out
}

mse <- readRDS(TISSUE_MSE)
pd  <- as.data.frame(pixelData(mse))
SIDS20  <- levels(pixelData(mse)$sample_id); if (is.null(SIDS20)) SIDS20 <- sort(unique(as.character(pd$sample_id)))
GROUP20 <- setNames(rep("incubated_5min", length(SIDS20)), SIDS20)
ord20   <- sort(SIDS20)
rm(mse); gc(verbose = FALSE)

RENDER_SIDS <- ord20
OUT_PDF     <- file.path(CACHE_DIR, "_tmp_island_geom.pdf")   # throwaway
OUT_SIDECAR <- file.path(CACHE_DIR, "organoid_island_label_positions_centroid.rds")

pdf(OUT_PDF, width = PAGE_W_IN, height = PAGE_H_IN)
sidecar <- list()
# page 1 = cover (matches R/09_organoid_refinement/02_island_cleanup_canvas.R has_cover=TRUE so section pages start at page 2)
par(mar = c(2, 2, 2, 2)); plot.new()

for (sid in RENDER_SIDS) {
  xf <- file.path(REG_CACHE, sprintf("nd2final_%s.rds", sid))
  pf <- file.path(CROP_DIR, sprintf("optical_%s.png", sid))
  inf_final <- file.path(CACHE_DIR, sprintf("instances_final_%s.rds", sid))
  inf_clean <- file.path(CACHE_DIR, sprintf("instances_clean_%s.rds", sid))
  inf_split <- file.path(CACHE_DIR, sprintf("instances_split_%s.rds", sid))
  inf <- if (file.exists(inf_final)) inf_final else if (file.exists(inf_clean)) inf_clean else
         if (file.exists(inf_split)) inf_split else file.path(CACHE_DIR, sprintf("instances_%s.rds", sid))
  grp <- GROUP20[sid]
  if (!file.exists(xf) || !file.exists(pf) || !file.exists(inf)) {
    par(mar = c(2, 2, 3, 2)); plot.new(); next     # keep page count aligned
  }
  Xr <- readRDS(xf); B <- Xr$B_msi_nd2
  cx0 <- Xr$crop[1]; cy0 <- Xr$crop[2]
  om <- png::readPNG(pf); if (length(dim(om)) == 3) om <- om[, , 1]
  cw <- ncol(om); ch <- nrow(om)

  inst <- readRDS(inf); inst <- inst[inst$instance > 0, ]
  ids <- sort(unique(as.integer(inst$instance)))

  # EXACT same device setup as R/09_organoid_refinement/02_island_cleanup_canvas.R section page (geometry only)
  par(mar = c(1, 1, 3, 1))
  plot.new(); plot.window(c(1, cw), c(ch, 1), asp = 1)

  for (k in ids) {
    cen <- apply_affine(B, matrix(c(mean(inst$x[inst$instance == k]),
                                    mean(inst$y[inst$instance == k])), nrow = 1))
    cx <- cen[1] - cx0 + 1; cy <- cen[2] - cy0 + 1
    ndc_x <- grconvertX(cx, "user", "ndc"); ndc_y <- grconvertY(cy, "user", "ndc")
    sidecar[[length(sidecar) + 1L]] <- data.frame(
      sid = sid, group = grp, instance = k,
      pdf_x = ndc_x * PAGE_W_PT, pdf_y_from_top = (1 - ndc_y) * PAGE_H_PT,
      n_px = sum(inst$instance == k), stringsAsFactors = FALSE)
  }
}
dev.off()
unlink(OUT_PDF)

sc <- do.call(rbind, sidecar)
page_of <- setNames(seq_along(RENDER_SIDS) + 1L, RENDER_SIDS)   # +1 for cover
sc$page <- as.integer(page_of[sc$sid])
sc <- sc[, c("page", "sid", "group", "instance", "pdf_x", "pdf_y_from_top", "n_px")]
attr(sc, "page_w_pt") <- PAGE_W_PT; attr(sc, "page_h_pt") <- PAGE_H_PT
attr(sc, "render_sids") <- RENDER_SIDS; attr(sc, "has_cover") <- TRUE
saveRDS(sc, OUT_SIDECAR)
cat(sprintf("[42b] centroid sidecar -> %s  (%d organoids on %d pages)\n",
            OUT_SIDECAR, nrow(sc), length(unique(sc$page))))
