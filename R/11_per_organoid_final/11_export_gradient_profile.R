#!/usr/bin/env Rscript
# ============================================================================
# 11_export_gradient_profile.R
# Export the per-organoid / per-zone outward gradient profile (`grad_prof`)
# that backs PAGE 6 of apical_citrate_dha_report.pdf. The report computes this
# in memory but never writes it to disk, so Page 6 cannot otherwise be
# reproduced in Prism. This script mirrors the report's computation EXACTLY
# (03_apical_report.R lines 60-205) but skips PDF / overlay rendering.
# Out: results/annotation/apical_gradient_profile_long.csv
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
suppressPackageStartupMessages(library(Cardinal))

CIT_MZ <- 191.0217; DHA_MZ <- 327.2330
TIC_PPM   <- 10
ANNOT_RES <- file.path(RES_DIR, "annotation")
APICAL_CSV  <- file.path(ANNOT_RES, "apical_map_consensus.csv")
CURATED_CSV <- file.path(RES_DIR, "curated_feature_table.csv")
CLASSES <- c("basolateral_out", "apical_out", "mixed")
stopifnot(file.exists(APICAL_CSV), file.exists(CURATED_CSV))

# ---- apical map ------------------------------------------------------------
ap <- read.csv(APICAL_CSV, stringsAsFactors = FALSE)
ap_map <- setNames(ap$apical_class, paste(ap$sid, ap$instance))

# ---- MSE + ion vectors (identical to report lines 60-79) -------------------
mse <- readRDS(TISSUE_MSE)
pd  <- as.data.frame(pixelData(mse)); pd$gidx <- seq_len(nrow(pd))
mzs <- mz(mse); di <- which.min(abs(mzs - DHA_MZ))
SIDS20 <- levels(pixelData(mse)$sample_id)
if (is.null(SIDS20)) SIDS20 <- sort(unique(as.character(pd$sample_id)))
ord20  <- sort(SIDS20)
SP <- as.matrix(spectra(mse)); rm(mse); gc(verbose = FALSE)
val_dha <- SP[di, ]
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_citrate.R"))
val_cit <- citrate_onto_pd(pd)
cur <- read.csv(CURATED_CSV, stringsAsFactors = FALSE)
F_idx <- vapply(cur$mz_curated, function(q) {
  j <- which.min(abs(mzs - q)); if (abs(mzs[j] - q)/q*1e6 <= TIC_PPM) j else NA_integer_
}, integer(1))
F_idx <- sort(unique(F_idx[!is.na(F_idx)]))
metTIC <- colSums(SP[F_idx, , drop = FALSE])
cat(sprintf("[grad] %d pixels; %d curated features mapped; zones = %s um\n",
            nrow(pd), length(F_idx), paste(OUT_ZONE_UM, collapse = "/")))

# ---- outward gradient loop (identical to report lines 152-205) -------------
MIN_ZONE_PX <- 3L
NZ          <- length(OUT_ZONE_UM)
FAR_ZONES   <- which(OUT_ZONE_UM >= 80)
# annotated organoids only (== report's `org` restricted to CLASSES)
ap_class_of <- ap_map[ap_map %in% CLASSES]

prof_rows <- list()
for (sid in ord20) {
  fz <- cache_in(sprintf("zones_%s.rds", sid)); if (!file.exists(fz)) next
  z <- readRDS(fz)
  for (k in sort(unique(z$instance[z$instance > 0]))) {
    key <- paste(sid, k); cls <- ap_class_of[key]
    if (is.na(cls)) next
    surf_g <- z$gidx[z$instance == k & !is.na(z$zone_in) & z$zone_in == 1]
    int_g  <- z$gidx[z$instance == k]
    cit_surf <- mean(val_cit[surf_g]); dha_surf <- mean(val_dha[surf_g])
    dha_int  <- mean(val_dha[int_g])
    cit_z <- dha_z <- mtic_z <- rep(NA_real_, NZ)
    for (zz in seq_len(NZ)) {
      gg <- z$gidx[z$instance_catch == k & !is.na(z$zone_out) & z$zone_out == zz]
      if (length(gg) >= MIN_ZONE_PX) {
        cit_z[zz] <- mean(val_cit[gg]); dha_z[zz] <- mean(val_dha[gg])
        gm <- gg[metTIC[gg] > 0]; if (length(gm)) mtic_z[zz] <- mean(val_cit[gm] / metTIC[gm]) * 100
      }
    }
    surf_norm <- if (is.finite(cit_surf) && cit_surf > 0) cit_z / cit_surf else cit_z * NA
    citdha    <- if (is.finite(dha_int)  && dha_int  > 0) cit_z / dha_int  else cit_z * NA
    dha_sn    <- if (is.finite(dha_surf) && dha_surf > 0) dha_z / dha_surf else dha_z * NA
    for (zz in seq_len(NZ))
      prof_rows[[length(prof_rows)+1L]] <- data.frame(
        sid = sid, instance = k, key = key, apical_class = unname(cls),
        zone = zz, zone_um = OUT_ZONE_UM[zz],
        surf = surf_norm[zz], abs = cit_z[zz], mtic = mtic_z[zz], citdha = citdha[zz],
        dha_sn = dha_sn[zz], dha_abs = dha_z[zz], stringsAsFactors = FALSE)
  }
}
grad_prof <- do.call(rbind, prof_rows)
out <- file.path(ANNOT_RES, "apical_gradient_profile_long.csv")
write.csv(grad_prof, out, row.names = FALSE)
cat(sprintf("[grad] %d rows (%d organoids x %d zones) -> %s\n",
            nrow(grad_prof), length(unique(grad_prof$key)), NZ, out))
print(table(unique(grad_prof[, c("key","apical_class")])$apical_class))
