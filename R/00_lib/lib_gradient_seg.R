#!/usr/bin/env Rscript
# ============================================================================
# lib_gradient_seg.R - reusable segmentation + signed-distance ring builder.
# Extracted from R/08_organoid_gradient_survey/01_segment_organoids.R +
# 02_buffer_rings.R so the all-20-dataset per-organoid analysis
# (R/11_per_organoid_final/01_zones_curated.R) can build instances_<sid>.rds /
# zones_<sid>.rds for any section WITHOUT touching the base gradient-survey
# scripts. Same logic, same outputs.
# ============================================================================
suppressPackageStartupMessages({ library(EBImage); library(RANN) })

# Segment one section's is_tissue mask -> data.frame(gidx,x,y,instance,is_surface).
gseg_segment <- function(sub) {           # sub: gidx,x,y,is_tissue for one section
  g <- section_grid(sub$x, sub$y)
  M <- matrix(FALSE, g$W, g$H)
  M[cbind(g$ix, g$iy)] <- ifelse(is.na(sub$is_tissue), FALSE, sub$is_tissue)
  lab <- EBImage::bwlabel(M)
  tab <- table(lab[lab > 0])
  if (length(tab)) {
    small <- as.integer(names(tab))[as.integer(tab) < MIN_INSTANCE_PX]
    lab[lab %in% small] <- 0L
  }
  ids <- sort(unique(as.integer(lab[lab > 0])))
  lab2 <- lab
  if (length(ids)) {
    remap <- setNames(seq_along(ids), as.character(ids))
    lab2[lab > 0] <- remap[as.character(lab[lab > 0])]
  }
  T <- lab2 > 0
  pad <- matrix(FALSE, g$W + 2L, g$H + 2L); pad[2:(g$W + 1L), 2:(g$H + 1L)] <- T
  n4 <- pad[1:g$W, 2:(g$H + 1L)] + pad[3:(g$W + 2L), 2:(g$H + 1L)] +
        pad[2:(g$W + 1L), 1:g$H]  + pad[2:(g$W + 1L), 3:(g$H + 2L)]
  surf_mat <- T & (n4 < 4)
  sub$instance   <- as.integer(lab2[cbind(g$ix, g$iy)])
  sub$is_surface <- as.logical(surf_mat[cbind(g$ix, g$iy)])
  sub$is_surface[is.na(sub$is_surface)] <- FALSE
  sub
}

# Add signed-distance zones to a segmented section (RANN nn2, Voronoi catchment).
gseg_rings <- function(sec) {             # sec: ...,is_tissue,instance,is_surface
  is_tis <- !is.na(sec$is_tissue) & sec$is_tissue
  surf <- sec[sec$is_surface, ]
  if (!nrow(surf)) return(NULL)
  nn <- RANN::nn2(data = surf[, c("x", "y")], query = sec[, c("x", "y")], k = 1)
  dist_um <- as.numeric(nn$nn.dists[, 1]) * MSI_PIXEL_UM
  sec$instance_catch <- surf$instance[as.integer(nn$nn.idx[, 1])]
  sec$dist_um <- dist_um
  sec$signed_dist_um <- ifelse(is_tis, -dist_um, dist_um)
  max_out <- max(OUT_BREAKS_UM)
  zo <- rep(NA_integer_, nrow(sec)); op <- !is_tis & dist_um <= max_out
  zo[op] <- as.integer(cut(dist_um[op], breaks = OUT_BREAKS_UM, include.lowest = TRUE, labels = FALSE))
  sec$zone_out <- zo
  zi <- rep(NA_integer_, nrow(sec))
  zi[is_tis] <- as.integer(cut(dist_um[is_tis], breaks = IN_BREAKS_UM, include.lowest = TRUE, labels = FALSE))
  sec$zone_in <- zi
  sec
}

# Ensure cache/instances_<sid>.rds + zones_<sid>.rds exist; build if missing.
gseg_ensure <- function(sid, pd) {
  # read-check upstream reuse first; only build (and write Final-local) if missing
  if (file.exists(cache_in(sprintf("instances_%s.rds", sid))) &&
      file.exists(cache_in(sprintf("zones_%s.rds", sid)))) return(invisible(FALSE))
  fi <- file.path(CACHE_DIR, sprintf("instances_%s.rds", sid))
  fz <- file.path(CACHE_DIR, sprintf("zones_%s.rds", sid))
  sub <- pd[as.character(pd$sample_id) == sid, c("gidx", "x", "y", "is_tissue")]
  inst <- gseg_segment(sub)
  saveRDS(inst, fi)
  sec <- merge(pd[as.character(pd$sample_id) == sid, c("gidx", "x", "y", "is_tissue")],
               inst[, c("gidx", "instance", "is_surface")], by = "gidx", all.x = TRUE)
  sec$instance[is.na(sec$instance)] <- 0L
  sec$is_surface[is.na(sec$is_surface)] <- FALSE
  z <- gseg_rings(sec)
  if (!is.null(z)) saveRDS(z, fz)
  invisible(TRUE)
}
