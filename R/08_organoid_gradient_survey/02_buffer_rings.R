#!/usr/bin/env Rscript
# ============================================================================
# 02_buffer_rings.R - Phase 08 step 2: signed-distance zones around each organoid.
# ----------------------------------------------------------------------------
# For every acquired pixel in a section, distance (um) to the NEAREST organoid
# surface pixel (RANN::nn2, k=1, x10 um/px). The nearest surface pixel's
# instance is the pixel's catchment (Voronoi assignment -> pixels nearer another
# organoid belong to that organoid, never double-counted).
#   OUTWARD (off-tissue / CMC gel): zone by OUT_BREAKS_UM 10/20/50/80/160 um.
#   INWARD  (on-tissue / organoid): zone by IN_BREAKS_UM 10 um steps (capped).
# Pixels farther than the outermost break are dropped (NA zone).
#
# Input : cache/instances_<sid>.rds (gidx,x,y,instance,is_surface)
# Output: cache/zones_<sid>.rds - data.frame(gidx,x,y,instance_catch,dist_um,
#           signed_dist_um, zone_out, zone_in)
#         results/gradient/zones_summary.csv
#         figures/gradient/zones_<sid>.pdf
# Usage : Rscript R/08_organoid_gradient_survey/02_buffer_rings.R
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
suppressPackageStartupMessages({
  library(Cardinal)
  library(RANN)
  library(viridisLite)
})

# Need all acquired pixels per section (incl. off-tissue gel) -> read MSE pixelData.
mse <- readRDS(TISSUE_MSE)
pd  <- as.data.frame(pixelData(mse))
pd$gidx <- seq_len(nrow(pd))

MAX_OUT_UM <- max(OUT_BREAKS_UM)   # 160
summ <- list()

for (sid in GRAD_SIDS) {
  inst <- readRDS(file.path(GRAD_CACHE, sprintf("instances_%s.rds", sid)))
  sel  <- as.character(pd$sample_id) == sid
  sec  <- pd[sel, c("gidx", "x", "y", "is_tissue")]
  # join instance/surface onto all acquired pixels of the section
  sec  <- merge(sec, inst[, c("gidx", "instance", "is_surface")], by = "gidx", all.x = TRUE)
  sec$instance[is.na(sec$instance)]     <- 0L
  sec$is_surface[is.na(sec$is_surface)] <- FALSE
  is_tis <- !is.na(sec$is_tissue) & sec$is_tissue

  surf <- sec[sec$is_surface, ]
  if (!nrow(surf)) { cat(sprintf("[50] %s: no surface px, skipped\n", sid)); next }

  # Nearest surface pixel (in grid px) for every acquired pixel
  nn <- RANN::nn2(data = surf[, c("x", "y")], query = sec[, c("x", "y")], k = 1)
  dist_px <- as.numeric(nn$nn.dists[, 1])
  dist_um <- dist_px * MSI_PIXEL_UM
  sec$instance_catch <- surf$instance[as.integer(nn$nn.idx[, 1])]   # Voronoi catchment

  # Signed distance: negative inward (tissue), positive outward (gel)
  sec$signed_dist_um <- ifelse(is_tis, -dist_um, dist_um)
  sec$dist_um <- dist_um

  # Outward zones (off-tissue only, <= MAX_OUT_UM)
  zone_out <- rep(NA_integer_, nrow(sec))
  out_px <- !is_tis & dist_um <= MAX_OUT_UM
  zone_out[out_px] <- as.integer(cut(dist_um[out_px], breaks = OUT_BREAKS_UM,
                                     include.lowest = TRUE, labels = FALSE))
  sec$zone_out <- zone_out

  # Inward zones (tissue only)
  zone_in <- rep(NA_integer_, nrow(sec))
  in_px <- is_tis
  zone_in[in_px] <- as.integer(cut(dist_um[in_px], breaks = IN_BREAKS_UM,
                                   include.lowest = TRUE, labels = FALSE))
  sec$zone_in <- zone_in

  saveRDS(sec, file.path(GRAD_CACHE, sprintf("zones_%s.rds", sid)))

  ot <- table(factor(zone_out, levels = seq_along(OUT_ZONE_UM)))
  it <- table(factor(zone_in,  levels = seq_along(IN_ZONE_LAB)))
  cat(sprintf("[50] %s: outward[%s] inward[%s]\n", sid,
              paste(as.integer(ot), collapse = ","),
              paste(as.integer(it), collapse = ",")))
  summ[[sid]] <- data.frame(sample_id = sid,
                            t(as.integer(ot)), t(as.integer(it)),
                            stringsAsFactors = FALSE)
  names(summ[[sid]]) <- c("sample_id",
                          paste0("out_", OUT_ZONE_UM, "um"),
                          paste0("in_", IN_ZONE_LAB, "um"))

  # QC PDF: signed-distance map + outward/inward zone maps
  g <- section_grid(sec$x, sec$y)
  fill_mat <- function(vals) { M <- matrix(NA_real_, g$W, g$H); M[cbind(g$ix, g$iy)] <- vals; M }
  pdf(file.path(GRAD_FIG, sprintf("zones_%s.pdf", sid)), width = 13, height = 5)
  par(mfrow = c(1, 3), mar = c(2, 2, 3, 4))
  sdm <- fill_mat(sec$signed_dist_um)
  image(seq_len(g$W), seq_len(g$H), sdm, col = viridisLite::magma(256), asp = 1,
        axes = FALSE, xlab = "", ylab = "", main = sprintf("%s signed dist (um)", sid))
  zo <- fill_mat(sec$zone_out)
  image(seq_len(g$W), seq_len(g$H), zo, col = viridisLite::viridis(length(OUT_ZONE_UM)),
        asp = 1, axes = FALSE, xlab = "", ylab = "", main = "outward zones 10/20/50/80/160")
  zi <- fill_mat(sec$zone_in)
  image(seq_len(g$W), seq_len(g$H), zi, col = viridisLite::mako(length(IN_ZONE_LAB)),
        asp = 1, axes = FALSE, xlab = "", ylab = "", main = "inward zones (10um steps)")
  dev.off()
}

summ_df <- do.call(rbind, summ)
write.csv(summ_df, file.path(GRAD_RES, "zones_summary.csv"), row.names = FALSE)
cat("\n[50] DONE.\n"); print(summ_df)
