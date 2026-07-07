#!/usr/bin/env Rscript
# ============================================================================
# 01_zones_curated.R - signed-distance outward/inward zones built from the
#   CURATED organoid instances (final > clean > split > base), for ALL 20
#   sections. Unlike R/08_organoid_gradient_survey/02_buffer_rings.R (base segmentation), the instance IDs
#   here match the curated IDs the apical annotations are keyed to, so the
#   gradient can be stratified by apical_class via merge(by=c("sid","instance")).
#
#   Surface is RECOMPUTED per curated instance (a pixel is surface if any
#   4-neighbour has a DIFFERENT instance, incl. 0) so touching split organoids
#   each get their own boundary. Zones then via lib_gradient_seg::gseg_rings.
#
# In : cache/peaks_tissue_combined.rds (pixelData: full grid incl. gel, is_tissue)
#      cache/instances_{final,clean,split,}_<sid>.rds
# Out: cache/zones_<sid>.rds (gidx,x,y,is_tissue,instance,instance_catch,
#        is_surface,dist_um,signed_dist_um,zone_out,zone_in)  -- one per section
#      results/gradient/zones_curated_summary.csv
# Usage: Rscript R/11_per_organoid_final/01_zones_curated.R
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_gradient_seg.R"))
suppressPackageStartupMessages({ library(Cardinal); library(RANN) })

inst_file <- function(sid) {
  for (suf in c("final", "clean", "split")) {
    f <- cache_in(sprintf("instances_%s_%s.rds", suf, sid))
    if (file.exists(f)) return(f)
  }
  cache_in(sprintf("instances_%s.rds", sid))
}

# surface of EACH curated instance: instance>0 pixel touching a different instance
recompute_surface <- function(sec) {
  g <- section_grid(sec$x, sec$y)
  L <- matrix(0L, g$W, g$H); L[cbind(g$ix, g$iy)] <- as.integer(sec$instance)
  pad <- matrix(-1L, g$W + 2L, g$H + 2L); pad[2:(g$W+1L), 2:(g$H+1L)] <- L
  up <- pad[2:(g$W+1L), 1:g$H]; dn <- pad[2:(g$W+1L), 3:(g$H+2L)]
  lf <- pad[1:g$W, 2:(g$H+1L)]; rt <- pad[3:(g$W+2L), 2:(g$H+1L)]
  diff_nb <- (up != L) | (dn != L) | (lf != L) | (rt != L)
  surf_mat <- (L > 0) & diff_nb
  sec$is_surface <- as.logical(surf_mat[cbind(g$ix, g$iy)])
  sec$is_surface[is.na(sec$is_surface)] <- FALSE
  sec
}

mse <- readRDS(TISSUE_MSE)
pd  <- as.data.frame(pixelData(mse)); pd$gidx <- seq_len(nrow(pd)); rm(mse); gc(verbose = FALSE)
SIDS20 <- sort(unique(as.character(pd$sample_id)))

summ <- list()
for (sid in SIDS20) {
  fi <- inst_file(sid)
  if (!file.exists(fi)) { cat(sprintf("[51] %s: no instance file, skip\n", sid)); next }
  inst <- readRDS(fi)
  sec  <- pd[as.character(pd$sample_id) == sid, c("gidx", "x", "y", "is_tissue")]
  sec  <- merge(sec, inst[, c("gidx", "instance")], by = "gidx", all.x = TRUE)
  sec$instance[is.na(sec$instance)] <- 0L
  sec  <- sec[order(sec$gidx), ]
  sec  <- recompute_surface(sec)
  if (!any(sec$is_surface)) { cat(sprintf("[51] %s: no surface px, skip\n", sid)); next }
  z <- gseg_rings(sec)                                   # adds instance_catch, dist, zones
  if (is.null(z)) { cat(sprintf("[51] %s: gseg_rings NULL, skip\n", sid)); next }
  saveRDS(z, file.path(CACHE_DIR, sprintf("zones_%s.rds", sid)))
  ot <- table(factor(z$zone_out, levels = seq_along(OUT_ZONE_UM)))
  n_org <- length(unique(z$instance[z$instance > 0]))
  cat(sprintf("[51] %s: %d organoids, outward[%s]\n", sid, n_org, paste(as.integer(ot), collapse = ",")))
  summ[[sid]] <- data.frame(sample_id = sid, n_organoid = n_org,
                            t(as.integer(ot)), stringsAsFactors = FALSE)
  names(summ[[sid]]) <- c("sample_id", "n_organoid", paste0("out_", OUT_ZONE_UM, "um"))
}
summ_df <- do.call(rbind, summ)
dir.create(GRAD_RES, showWarnings = FALSE, recursive = TRUE)
write.csv(summ_df, file.path(GRAD_RES, "zones_curated_summary.csv"), row.names = FALSE)
cat(sprintf("\n[51] DONE -> %d zones_<sid>.rds files\n", length(summ)))
print(summ_df, row.names = FALSE)
