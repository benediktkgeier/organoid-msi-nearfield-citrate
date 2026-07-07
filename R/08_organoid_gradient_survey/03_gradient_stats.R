#!/usr/bin/env Rscript
# ============================================================================
# 03_gradient_stats.R - per-ion outward/inward gradient ranking (POOLED).
# ----------------------------------------------------------------------------
# Single-condition study (all sections incubated in CMC >=5 min, treated
# identically). This is a POOLED descriptive survey across ALL sections: for
# each ion, pool every section's pixels and take the MEAN intensity in each
# distance zone, then:
#   rho_out = Spearman(zone idx, mean intensity) over the OUTWARD zones
#             (10/20/50/80/160/250/500 um). More positive = shallower outward
#             decay (more signal carried into the gel).
#   rho_in  = same over the INWARD zones (surface -> core).
# Descriptive only (organoids within a section are pseudo-replicates); no group
# split, no Delta-rho, no p-values.
#
# Input : cache/peaks_tissue_combined.rds, cache/zones_<sid>.rds
# Output: results/gradient/gradient_stats.csv      (one row per ion, ranked)
#         results/gradient/zone_profiles_long.csv  (ion x zone means)
# Usage : Rscript R/08_organoid_gradient_survey/03_gradient_stats.R
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
suppressPackageStartupMessages(library(Cardinal))

mse <- readRDS(TISSUE_MSE)
pd  <- as.data.frame(pixelData(mse))
fd  <- as.data.frame(featureData(mse))
mzs <- mz(mse)
nf  <- nrow(mse)
X   <- as.matrix(spectra(mse))                       # nf x n_pixels (all sections)
cat(sprintf("[03] %d ions x %d px (all sections, pooled)\n", nf, ncol(X)))

# Per-pixel zone lookup (z$gidx is the global pixel index into pd row order)
zone_out <- rep(NA_integer_, ncol(X))
zone_in  <- rep(NA_integer_, ncol(X))
for (sid in GRAD_SIDS) {
  zf <- cache_in(sprintf("zones_%s.rds", sid)); if (!file.exists(zf)) next
  z <- readRDS(zf)
  zone_out[z$gidx] <- z$zone_out
  zone_in[z$gidx]  <- z$zone_in
}

n_out <- length(OUT_ZONE_UM)   # 7
n_in  <- length(IN_ZONE_LAB)   # 5

# mean intensity matrix: ion x zone, pooled over ALL sections' pixels
zone_means <- function(zone_vec, nz) {
  M <- matrix(NA_real_, nf, nz); ncnt <- integer(nz)
  for (z in seq_len(nz)) {
    cols <- which(zone_vec == z); ncnt[z] <- length(cols)
    if (length(cols)) M[, z] <- rowMeans(X[, cols, drop = FALSE])
  }
  list(M = M, n = ncnt)
}
spearman_row <- function(v) {
  if (all(is.na(v))) return(NA_real_)
  idx <- seq_along(v)[!is.na(v)]
  if (length(idx) < 3 || stats::sd(v[idx]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(idx, v[idx], method = "spearman"))
}

out <- zone_means(zone_out, n_out); inn <- zone_means(zone_in, n_in)
cat(sprintf("[03] outward zone n: [%s]\n", paste(out$n, collapse = ",")))
cat(sprintf("[03] inward  zone n: [%s]\n", paste(inn$n, collapse = ",")))

rho_out <- apply(out$M, 1, spearman_row)
rho_in  <- apply(inn$M, 1, spearman_row)

res <- data.frame(
  feat = seq_len(nf), mz = round(mzs, 4),
  is_known = fd$is_known, is_on_tissue = fd$is_on_tissue,
  metabolite_name = fd$metabolite_name, metabolite_adduct = fd$metabolite_adduct,
  metabolite_score = fd$metabolite_score,
  rho_out = round(rho_out, 4), rho_in = round(rho_in, 4),
  stringsAsFactors = FALSE)
res$label <- ifelse(!is.na(res$metabolite_name) & nzchar(res$metabolite_name),
                    sprintf("%s (%.4f)", res$metabolite_name, res$mz),
                    sprintf("m/z %.4f (unknown)", res$mz))
res <- res[order(-res$rho_out), ]
write.csv(res, file.path(GRAD_RES, "gradient_stats.csv"), row.names = FALSE)

# Long zone-profile table (for the report)
long <- list()
add_long <- function(dir, zmlist, zlabs) {
  M <- zmlist$M
  for (z in seq_along(zlabs)) {
    long[[length(long) + 1]] <<- data.frame(
      feat = seq_len(nf), mz = round(mzs, 4), direction = dir,
      zone = z, zone_label = zlabs[z], n_px = zmlist$n[z],
      mean_intensity = M[, z], stringsAsFactors = FALSE)
  }
}
add_long("outward", out, paste0(OUT_ZONE_UM, "um"))
add_long("inward",  inn, IN_ZONE_LAB)
long_df <- do.call(rbind, long)
write.csv(long_df, file.path(GRAD_RES, "zone_profiles_long.csv"), row.names = FALSE)

cat("\n[03] DONE. Top outward ions (largest rho_out):\n")
top <- head(res[!is.na(res$rho_out), c("label", "is_known", "rho_out", "rho_in")], 15)
print(top, row.names = FALSE)