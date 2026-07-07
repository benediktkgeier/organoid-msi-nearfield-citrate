#!/usr/bin/env Rscript
# ============================================================================
# 01_segment_organoids.R - Phase 08 step 1: per-section organoid instance segmentation.
# ----------------------------------------------------------------------------
# The on-tissue mask `is_tissue` (R/04_ssc_ontissue/04_ssc_tissue_mask.R floor80 rule) IS the organoid+edge
# footprint; no cluster_decisions needed. Per section, label connected
# components of is_tissue (EBImage::bwlabel, 4-conn), drop instances < 50 px,
# and mark each instance's SURFACE = on-tissue pixels with a non-tissue
# 4-neighbour (or grid edge). Runs on the FOUR gradient datasets only.
#
# Input : cache/peaks_tissue_combined.rds (pixelData has x,y,sample_id,is_tissue)
# Output: cache/instances_<sid>.rds  - data.frame(gidx,x,y,instance,is_surface)
#         results/gradient/instances_summary.csv
#         figures/gradient/instances_<sid>.pdf (label + surface overlay)
# Usage : Rscript R/08_organoid_gradient_survey/01_segment_organoids.R
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
suppressPackageStartupMessages({
  library(Cardinal)
  library(EBImage)
  library(viridisLite)
})

stopifnot(file.exists(TISSUE_MSE))
cat("[40] Loading tissue MSE...\n")
mse <- readRDS(TISSUE_MSE)
pd  <- as.data.frame(pixelData(mse))
pd$gidx <- seq_len(nrow(pd))
stopifnot(all(c("x", "y", "sample_id", "is_tissue") %in% names(pd)))

summ <- list()
for (sid in GRAD_SIDS) {
  sel <- as.character(pd$sample_id) == sid
  if (!any(sel)) stop("Section not found in MSE: ", sid)
  sub <- pd[sel, c("gidx", "x", "y", "is_tissue")]
  g   <- section_grid(sub$x, sub$y)

  # Binary tissue mask on the section grid
  M <- matrix(FALSE, g$W, g$H)
  M[cbind(g$ix, g$iy)] <- ifelse(is.na(sub$is_tissue), FALSE, sub$is_tissue)

  # Connected components (EBImage::bwlabel = 4-connectivity in 2D)
  lab <- EBImage::bwlabel(M)
  tab <- table(lab[lab > 0])
  if (length(tab)) {
    small <- as.integer(names(tab))[as.integer(tab) < MIN_INSTANCE_PX]
    lab[lab %in% small] <- 0L
  }
  # Relabel sequentially 1..K
  ids <- sort(unique(as.integer(lab[lab > 0])))
  remap <- setNames(seq_along(ids), as.character(ids))
  lab2 <- lab
  if (length(ids)) lab2[lab > 0] <- remap[as.character(lab[lab > 0])]

  # Surface = tissue pixel with <4 tissue 4-neighbours (grid edge counts as non-tissue)
  T <- lab2 > 0
  pad <- matrix(FALSE, g$W + 2L, g$H + 2L)
  pad[2:(g$W + 1L), 2:(g$H + 1L)] <- T
  n4 <- pad[1:g$W, 2:(g$H + 1L)] + pad[3:(g$W + 2L), 2:(g$H + 1L)] +
        pad[2:(g$W + 1L), 1:g$H]  + pad[2:(g$W + 1L), 3:(g$H + 2L)]
  surf_mat <- T & (n4 < 4)

  # Map matrix -> per-pixel
  sub$instance   <- as.integer(lab2[cbind(g$ix, g$iy)])
  sub$is_surface <- as.logical(surf_mat[cbind(g$ix, g$iy)])
  sub$is_surface[is.na(sub$is_surface)] <- FALSE

  n_inst <- length(ids)
  inst_px <- if (n_inst) as.integer(table(factor(sub$instance[sub$instance > 0],
                                                  levels = seq_len(n_inst)))) else integer(0)
  cat(sprintf("[40] %s: %d tissue px -> %d instance(s) [%s], %d surface px\n",
              sid, sum(T), n_inst,
              if (n_inst) paste(inst_px, collapse = ",") else "-",
              sum(sub$is_surface)))

  saveRDS(sub, file.path(GRAD_CACHE, sprintf("instances_%s.rds", sid)))

  for (k in seq_len(n_inst)) {
    summ[[length(summ) + 1]] <- data.frame(
      sample_id = sid, instance = k,
      n_px = inst_px[k],
      n_surface = sum(sub$is_surface & sub$instance == k),
      stringsAsFactors = FALSE)
  }

  # QC PDF: instance label map + surface overlay
  pdf(file.path(GRAD_FIG, sprintf("instances_%s.pdf", sid)), width = 11, height = 6)
  par(mfrow = c(1, 2), mar = c(2, 2, 3, 1))
  pal <- c("grey15", grDevices::hcl.colors(max(n_inst, 1), "Dark 3"))
  image(seq_len(g$W), seq_len(g$H), lab2 + 1L, col = pal, asp = 1,
        axes = FALSE, xlab = "", ylab = "", main = sprintf("%s - instances (n=%d)", sid, n_inst))
  base <- matrix(0L, g$W, g$H); base[T] <- 1L
  image(seq_len(g$W), seq_len(g$H), base, col = c("grey15", "grey55"), asp = 1,
        axes = FALSE, xlab = "", ylab = "", main = "surface (red)")
  sp <- which(surf_mat, arr.ind = TRUE)
  if (nrow(sp)) points(sp[, 1], sp[, 2], pch = 15, cex = 0.35, col = "red")
  dev.off()
}

summ_df <- do.call(rbind, summ)
write.csv(summ_df, file.path(GRAD_RES, "instances_summary.csv"), row.names = FALSE)
cat(sprintf("\n[40] DONE. %d instances across %d sections.\n",
            nrow(summ_df), length(GRAD_SIDS)))
print(summ_df)
