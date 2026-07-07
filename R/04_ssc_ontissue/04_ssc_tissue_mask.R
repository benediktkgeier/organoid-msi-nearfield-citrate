#!/usr/bin/env Rscript
# ============================================================================
# 04_ssc_tissue_mask.R - Phase 04 step 4: PER-SECTION two-pass
#   spatialShrunkenCentroids + cross-section CHEMOTYPE harmonization.
# ============================================================================
# DESIGN (user decision 2026-06-16): cluster each of the 20 sections INDEPENDENTLY
# (both passes), then HARMONIZE the per-section substructure clusters into a
# controlled set of "chemotypes" by hierarchical meta-clustering of their mean
# spectra. This avoids large/intense sections dominating shared centroids and
# keeps each section's chemistry distinct, while the chemotype label restores
# cross-section comparability ("color clusters by chemical composition").
#
# Input (R/04_ssc_ontissue/03_build_curated_set.R, 20-section self-contained MSE):
#   cache/peaks_curated.rds  - the 348 curated knowns+unknowns ions; the ONLY
#       data used here, for BOTH clustering and the tissueness proxy.
#
# Per section (each pixelData$sample_id):
#   Pass 1: k=4, s={6,9,12}  -> tissue vs CMC vs edge; on-tissue = clusters with
#                               mean curated signal >=50% of the richest, UNION
#                               pixels >= section 80th-pct curated signal (floor80:
#                               recovers epithelium edges, keeps the lumen excluded)
#   Pass 2: k=10, s={6,9,12} -> within-tissue substructure (local cluster ids)
# Harmonization (s=9): pool every (section, Pass-2 cluster) mean spectrum ->
#   1-Pearson(log1p) distance -> hclust(ward.D2) -> silhouette-chosen K cut ->
#   chemotype id per local cluster -> per-pixel pixelData$chemotype.
#
# Parallelism: register_parallel() SnowParam; BPPARAM passed explicitly.
#   NEVER setCardinalBPPARAM() (v3 bug). SerialParam fallback on error.
#
# Outputs:
#   cache/peaks_tissue_combined.rds  - clean MSE + pixelData: ssc_k4_sec(+_s*),
#       tissueness_px, is_tissue, ssc_k10_sec(+_s*, local), chemotype
#   results/ssc_clusters.csv         - per (sample_id,pass,s,local_cluster) + role
#   results/chemotype_map.csv        - chemotype -> members, top m/z, role/color
#   results/ssc_log.csv
#   figures/ssc/ssc_mask_<sid>.pdf          (per section)
#   figures/ssc/chemotype_harmonization.pdf (global)
#
# Optional override: cache/chemotype_k_override.rds (single integer K).
# Usage: Rscript R/04_ssc_ontissue/04_ssc_tissue_mask.R
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_paths.R"))
source("D:/R/Projects/BIG_MSI/HP_Comparative_2026/R/lib_ion_image.R")
suppressPackageStartupMessages({
  library(Cardinal)
  library(viridisLite)
})

# Everything in this analysis uses ONLY the 348 curated ions (R/04_ssc_ontissue/03_build_curated_set.R): both the
# clustering AND the per-pixel tissueness proxy. The 8,772 annot MSE is not used.
CLEAN_MSE  <- file.path(CACHE_DIR, "peaks_curated.rds")
FIG_SSC    <- file.path(FIG_DIR, "ssc")
OUT_RDS    <- file.path(CACHE_DIR, "peaks_tissue_combined.rds")
OUT_CLUST  <- file.path(RES_DIR, "ssc_clusters.csv")
OUT_CHEMO  <- file.path(RES_DIR, "chemotype_map.csv")
OUT_LOG    <- file.path(RES_DIR, "phase3_ssc_log.csv")
K_OVERRIDE <- file.path(CACHE_DIR, "chemotype_k_override.rds")
dir.create(FIG_SSC, recursive = TRUE, showWarnings = FALSE)

SSC_R          <- 2
S_GRID         <- c(6, 9, 12)
S_DEFAULT      <- 9
K_PASS1        <- 4
K_PASS2        <- 10
OFF_LABELS     <- c("off_tissue", "matrix")
NOISE_LABEL    <- "noise"
TISSUENESS_THR <- 0.5         # cluster cut: keep clusters >= 50% of richest cluster's mean
FLOOR_PCT      <- 0.80        # per-pixel floor: also keep pixels >= section 80th-pct curated
                             # signal (the 'floor80' rule: catches epithelium edges the cluster
                             # cut lumps into gel clusters, without flooding the low-signal lumen)
K_CHEMO_RANGE  <- 4:10        # silhouette search window for chemotype count
set.seed(42)

t_start <- Sys.time()

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Run SSC trying SnowParam (bp), falling back to SerialParam on error.
run_ssc <- function(x, k, bp) {
  call_ssc <- function(PP) {
    spatialShrunkenCentroids(x, r = SSC_R, k = k, s = S_GRID,
                             weights = "adaptive", BPPARAM = PP)
  }
  res <- tryCatch(call_ssc(bp), error = function(e) e)
  if (inherits(res, "error")) {
    cat(sprintf("[20]   SnowParam SSC failed (%s); retry SerialParam...\n",
                conditionMessage(res)))
    res <- call_ssc(BiocParallel::SerialParam())
    attr(res, "bpused") <- "SerialParam"
  } else attr(res, "bpused") <- "SnowParam"
  res
}

# ResultsList (one per s) -> list keyed by S_GRID, each a length-n cluster vec.
extract_by_s <- function(ssc, n) {
  one <- function(obj) {
    cid <- NULL
    cl <- try(obj$class, silent = TRUE)
    if (!inherits(cl, "try-error") && !is.null(cl)) cid <- as.integer(cl)
    if (is.null(cid)) {
      m <- try(slot(obj, "model"), silent = TRUE)
      if (!inherits(m, "try-error") && !is.null(m$class)) cid <- as.integer(m$class)
    }
    if (is.null(cid)) {
      fm <- try(fitted(obj), silent = TRUE)
      if (!inherits(fm, "try-error") && !is.null(fm)) {
        fmt <- if (is.matrix(fm)) fm else matrix(as.numeric(fm), nrow = n)
        cid <- as.integer(apply(fmt, 1, which.max))
      }
    }
    if (is.null(cid) || length(cid) != n)
      stop(sprintf("Could not extract %d-pixel clusters from SSC result", n))
    cid
  }
  out <- if (is(ssc, "ResultsList")) {
    stopifnot(length(ssc) == length(S_GRID))
    lapply(seq_along(S_GRID), function(i) one(ssc[[i]]))
  } else list(one(ssc))
  names(out) <- paste0("s", S_GRID)
  out
}

# render_ion_image() img list for ONE sample from a value vector.
build_img <- function(x, y, vals) {
  xs <- sort(unique(x)); ys <- sort(unique(y))
  mat <- matrix(NA_real_, nrow = length(xs), ncol = length(ys),
                dimnames = list(xs, ys))
  mat[cbind(match(x, xs), match(y, ys))] <- as.numeric(vals)
  list(matrix = mat, x_coord = xs, y_coord = ys, mz_matched = NA, ppm_err = NA,
       n_nonzero = sum(vals > 0, na.rm = TRUE), n_total = length(vals))
}
cat_palette <- function(k) grDevices::hcl.colors(max(k, 2), palette = "Dark 3")

# Mean silhouette width for a labeling, given a full symmetric distance matrix.
sil_avg <- function(dm, lab) {
  n <- length(lab); s <- numeric(n); ug <- unique(lab)
  if (length(ug) < 2) return(-Inf)
  for (i in seq_len(n)) {
    same <- which(lab == lab[i]); same <- same[same != i]
    if (!length(same)) { s[i] <- 0; next }
    a <- mean(dm[i, same])
    other <- ug[ug != lab[i]]
    b <- min(vapply(other, function(g) mean(dm[i, which(lab == g)]), numeric(1)))
    s[i] <- if (max(a, b) > 0) (b - a) / max(a, b) else 0
  }
  mean(s)
}

# ===========================================================================
# 1. Load inputs + guards
# ===========================================================================
if (!file.exists(CLEAN_MSE))
  stop("Missing input: ", CLEAN_MSE, "\n  Run R/21_build_curated_set.R first.")
cat("[20] Loading curated (348-ion) SSC input MSE...\n")
clean <- readRDS(CLEAN_MSE)
npix  <- ncol(clean)
cat(sprintf("[20] curated: %d feat x %d px\n", nrow(clean), npix))

stopifnot(nlevels(run(clean)) == 20)
pc <- as.data.frame(pixelData(clean))
sample_ids <- levels(pixelData(clean)$sample_id)
if (is.null(sample_ids)) sample_ids <- sort(unique(as.character(pc$sample_id)))

# Self-containment sanity: a mid-matrix spectra read must succeed (was the crash).
sp_clean <- spectra(clean)
cat(sprintf("[20] clean spectra class: %s\n", class(sp_clean)[1]))
.probe <- as.matrix(sp_clean[seq_len(min(3, nrow(clean))),
                             round(npix * 0.6) + 0:4, drop = FALSE])
cat(sprintf("[20] mid-matrix read OK (%s) -> self-contained\n",
            paste(dim(.probe), collapse = "x")))

# ===========================================================================
# 2. Per-pixel "tissueness" from the 348 CURATED ions only
# ===========================================================================
# All curated ions are real/on-tissue signal, so a pixel's TOTAL curated-ion
# intensity is a direct tissue-richness proxy: organoid pixels are rich in these
# metabolites (high), empty CMC/matrix pixels carry little (low). Normalize to
# [0,1] by the global p99.5 so the per-cluster tissue-cut is on a stable scale.
# NOTE: this analysis uses ONLY the 348 curated ions -- the 8,772 annot set is
# not loaded or used anywhere.
cur_sig <- as.numeric(colSums(spectra(clean)))
hi_sig  <- as.numeric(quantile(cur_sig[cur_sig > 0], IMG_CLIP_HI, na.rm = TRUE))
tissueness_px <- pmin(cur_sig / hi_sig, 1)
cat(sprintf("[20] curated-signal tissueness: median %.3f (p99.5 ref %.3g)\n",
            median(tissueness_px), hi_sig))

bp <- register_parallel()

# Global pixel-aligned outputs (filled per section)
k4_g  <- setNames(lapply(S_GRID, function(s) rep(NA_integer_, npix)), paste0("s", S_GRID))
k10_g <- setNames(lapply(S_GRID, function(s) rep(NA_integer_, npix)), paste0("s", S_GRID))
is_tissue <- rep(FALSE, npix)

centroids   <- list()  # one per (section, pass2 s9 local cluster)
cent_meta   <- list()
clust_rows  <- list()
cut_log     <- list()
bpused      <- NA_character_

# ===========================================================================
# 3. Per-section two-pass SSC
# ===========================================================================
for (sid in sample_ids) {
  gidx <- which(as.character(pc$sample_id) == sid)   # global pixel indices
  cat(sprintf("\n[20] === Section %s (%d px) ===\n", sid, length(gidx)))
  sec <- clean[, gidx]

  # ---- Pass 1: k=4 ----
  ssc1 <- run_ssc(sec, K_PASS1, bp)
  if (is.na(bpused)) bpused <- attr(ssc1, "bpused")
  k4_by_s <- extract_by_s(ssc1, length(gidx))
  for (s in S_GRID) k4_g[[paste0("s", s)]][gidx] <- k4_by_s[[paste0("s", s)]]
  k4_def <- k4_by_s[[paste0("s", S_DEFAULT)]]

  # ---- Per-cluster tissueness (k4 s9) -> auto tissue-cut ----
  tn_sec <- tissueness_px[gidx]
  agg <- aggregate(tn ~ cl, data = data.frame(cl = k4_def, tn = tn_sec), FUN = mean)
  agg$n <- as.integer(table(k4_def)[as.character(agg$cl)])
  agg <- agg[order(agg$tn), ]
  keep <- agg$cl[agg$tn >= TISSUENESS_THR * max(agg$tn)]   # relative: >=50% of richest cluster
  if (length(keep) == 0) keep <- agg$cl[nrow(agg)]          # guard: always keep the richest
  cut <- setdiff(agg$cl, keep)
  in_cluster <- k4_def %in% keep
  # Per-pixel floor (the 'floor80' rule): also keep high-curated-signal pixels the
  # cluster cut missed (epithelium edges lumped into gel-dominated clusters), WITHOUT
  # flooding the low-signal lumen the way a morphological buffer would.
  by_floor <- tn_sec >= as.numeric(quantile(tn_sec, FLOOR_PCT))
  sec_is_tissue <- in_cluster | by_floor
  is_tissue[gidx] <- sec_is_tissue
  cat(sprintf("[20]   keep k4 {%s} / cut {%s}; cluster=%d +floor(p%.0f)=%d -> %d/%d tissue px\n",
              paste(sort(keep), collapse = ","), paste(sort(cut), collapse = ","),
              sum(in_cluster), FLOOR_PCT * 100, sum(by_floor & !in_cluster),
              sum(sec_is_tissue), length(gidx)))
  cut_log[[sid]] <- data.frame(sample_id = sid, cut_clusters = paste(sort(cut), collapse = ","),
                               n_cluster = sum(in_cluster), n_floor = sum(by_floor & !in_cluster),
                               n_tissue = sum(sec_is_tissue), n_px = length(gidx))

  # ---- Pass 2: k=10 within tissue ----
  gidx_t <- gidx[sec_is_tissue]
  if (length(gidx_t) >= K_PASS2) {
    sec_t <- clean[, gidx_t]
    ssc2  <- run_ssc(sec_t, K_PASS2, bp)
    k10_by_s <- extract_by_s(ssc2, length(gidx_t))
    for (s in S_GRID) k10_g[[paste0("s", s)]][gidx_t] <- k10_by_s[[paste0("s", s)]]
    k10_def <- k10_by_s[[paste0("s", S_DEFAULT)]]
  } else {
    cat(sprintf("[20]   WARNING: only %d tissue px (<%d) - Pass 2 skipped\n",
                length(gidx_t), K_PASS2))
    k10_def <- integer(0)
  }

  # ---- Centroids per Pass-2 (s9) local cluster (mean spectrum, clean feats) ----
  if (length(gidx_t) >= K_PASS2) {
    for (cl in sort(unique(k10_def))) {
      px <- gidx_t[k10_def == cl]
      cen <- as.numeric(rowMeans(sp_clean[, px, drop = FALSE]))
      key <- sprintf("%s|c%d", sid, cl)
      centroids[[key]] <- cen
      cent_meta[[key]] <- data.frame(sample_id = sid, local_cluster = cl,
                                     n_pixels = length(px),
                                     mean_tissueness = round(mean(tissueness_px[px]), 4),
                                     stringsAsFactors = FALSE)
    }
  }

  # ---- Cluster-table rows (both passes, all s) ----
  add_rows <- function(pass, by_s, idx_space, is_k4) {
    for (s in S_GRID) {
      cid <- by_s[[paste0("s", s)]]
      if (!length(cid)) next
      for (cl in sort(unique(cid))) {
        sel <- which(cid == cl)
        clust_rows[[length(clust_rows) + 1]] <<- data.frame(
          sample_id = sid, pass = pass, s = s, local_cluster = cl,
          n_pixels = length(sel),
          mean_tissueness = round(mean(tissueness_px[idx_space[sel]]), 4),
          is_tissue = if (is_k4 && s == S_DEFAULT) (cl %in% keep) else NA,
          chemotype = NA_integer_, role = "", stringsAsFactors = FALSE)
      }
    }
  }
  add_rows("k4", k4_by_s, gidx, TRUE)
  if (length(gidx_t) >= K_PASS2) add_rows("k10", k10_by_s, gidx_t, FALSE)
}

# ===========================================================================
# 4. Chemotype harmonization (hierarchical meta-clustering on s9 centroids)
# ===========================================================================
cent_keys <- names(centroids)
ncent <- length(cent_keys)
cat(sprintf("\n[20] Harmonizing %d per-section Pass-2 centroids into chemotypes...\n", ncent))
chemo_of_key <- setNames(rep(NA_integer_, ncent), cent_keys)
K_chemo <- NA_integer_; hc <- NULL; cm <- NULL
if (ncent >= 3) {
  cmat <- do.call(rbind, centroids[cent_keys])     # ncent x nfeat
  lc   <- log1p(cmat)
  cm   <- cor(t(lc)); cm[is.na(cm)] <- 0           # ncent x ncent
  dm   <- 1 - cm
  hc   <- hclust(as.dist(dm), method = "ward.D2")
  if (file.exists(K_OVERRIDE)) {
    K_chemo <- as.integer(readRDS(K_OVERRIDE))
    cat(sprintf("[20]   K override = %d\n", K_chemo))
  } else {
    krange <- K_CHEMO_RANGE[K_CHEMO_RANGE < ncent]
    sils <- vapply(krange, function(k) sil_avg(dm, cutree(hc, k = k)), numeric(1))
    K_chemo <- krange[which.max(sils)]
    cat(sprintf("[20]   silhouette-chosen K = %d (avg sil %.3f over K in [%d,%d])\n",
                K_chemo, max(sils), min(krange), max(krange)))
  }
  chemo_of_key[cent_keys] <- cutree(hc, k = K_chemo)
} else {
  cat("[20]   <3 centroids - assigning each its own chemotype\n")
  K_chemo <- ncent
  chemo_of_key[cent_keys] <- seq_len(ncent)
}

# Per-pixel chemotype: (section, k10 s9 local cluster) -> chemotype
chemotype_px <- rep(NA_integer_, npix)
k10_s9 <- k10_g[[paste0("s", S_DEFAULT)]]
for (key in cent_keys) {
  parts <- strsplit(key, "\\|c")[[1]]; sid <- parts[1]; cl <- as.integer(parts[2])
  gidx <- which(as.character(pc$sample_id) == sid & k10_s9 == cl)
  chemotype_px[gidx] <- chemo_of_key[[key]]
}
# Back-fill chemotype into cluster table (k10 s9 rows)
for (i in seq_along(clust_rows)) {
  r <- clust_rows[[i]]
  if (r$pass == "k10" && r$s == S_DEFAULT) {
    key <- sprintf("%s|c%d", r$sample_id, r$local_cluster)
    if (key %in% names(chemo_of_key)) clust_rows[[i]]$chemotype <- chemo_of_key[[key]]
  }
}

# ===========================================================================
# 5. Attach pixelData + save
# ===========================================================================
out <- clean
pixelData(out)$tissueness_px <- tissueness_px
pixelData(out)$is_tissue     <- is_tissue
pixelData(out)$ssc_k4_sec    <- k4_g[[paste0("s", S_DEFAULT)]]
pixelData(out)$ssc_k10_sec   <- k10_s9
pixelData(out)$chemotype     <- chemotype_px
for (s in S_GRID) {
  pixelData(out)[[sprintf("ssc_k4_sec_s%d", s)]]  <- k4_g[[paste0("s", s)]]
  pixelData(out)[[sprintf("ssc_k10_sec_s%d", s)]] <- k10_g[[paste0("s", s)]]
}
saveRDS(out, OUT_RDS)
cat(sprintf("[20] Saved: %s\n", basename(OUT_RDS)))

# ===========================================================================
# 6. results/ssc_clusters.csv + chemotype_map.csv
# ===========================================================================
clust_df <- do.call(rbind, clust_rows)
write.csv(clust_df, OUT_CLUST, row.names = FALSE)
cat(sprintf("[20] Wrote %s (%d rows)\n", basename(OUT_CLUST), nrow(clust_df)))

mz_vec <- mz(clean)
chemo_rows <- list()
for (ch in sort(unique(na.omit(chemotype_px)))) {
  keys <- cent_keys[chemo_of_key[cent_keys] == ch]
  members <- paste(keys, collapse = ";")
  px_tot  <- sum(chemotype_px == ch, na.rm = TRUE)
  persamp <- sapply(sample_ids, function(s)
    sum(chemotype_px == ch & as.character(pc$sample_id) == s, na.rm = TRUE))
  # top discriminating m/z: mean centroid of this chemotype vs grand mean
  cen_ch <- rowMeans(do.call(cbind, centroids[keys]), na.rm = TRUE)
  cen_all <- rowMeans(do.call(cbind, centroids), na.rm = TRUE)
  top_mz <- mz_vec[order(cen_ch - cen_all, decreasing = TRUE)[1:5]]
  r <- data.frame(chemotype = ch, n_members = length(keys), n_pixels = px_tot,
                  mean_tissueness = round(mean(tissueness_px[chemotype_px == ch], na.rm = TRUE), 4),
                  top_mz = paste(sprintf("%.4f", top_mz), collapse = ";"),
                  members = members, role = "", color = cat_palette(K_chemo)[ch],
                  stringsAsFactors = FALSE)
  for (s in sample_ids) r[[paste0("npx_", s)]] <- persamp[[s]]
  chemo_rows[[length(chemo_rows) + 1]] <- r
}
chemo_df <- do.call(rbind, chemo_rows)
write.csv(chemo_df, OUT_CHEMO, row.names = FALSE)
cat(sprintf("[20] Wrote %s (%d chemotypes)\n", basename(OUT_CHEMO), nrow(chemo_df)))

# ===========================================================================
# 7. Per-section QC PDFs
# ===========================================================================
chemo_pal <- cat_palette(K_chemo)
render_cluster_panel <- function(x, y, cid, main, pal, k) {
  img <- build_img(x, y, cid)
  render_ion_image(img, clip_hi = 1.0, palette = if (missing(pal)) cat_palette(k) else pal,
                   main = main, add_colorbar = FALSE)
}
for (sid in sample_ids) {
  sel <- as.character(pc$sample_id) == sid
  x <- pc$x[sel]; y <- pc$y[sel]
  pdf(file.path(FIG_SSC, sprintf("ssc_mask_%s.pdf", sid)), width = 14, height = 8.5)

  par(mfrow = c(1, 3), oma = c(1, 1, 3, 1))
  for (s in S_GRID)
    render_cluster_panel(x, y, k4_g[[paste0("s", s)]][sel], sprintf("k4 s=%d", s), k = K_PASS1)
  mtext(sprintf("%s - Pass 1 SSC k=4 (per-section) - cluster maps by s", sid),
        outer = TRUE, font = 2)

  par(mfrow = c(1, 2), oma = c(1, 1, 3, 1))
  render_ion_image(build_img(x, y, tissueness_px[sel]), clip_hi = IMG_CLIP_HI,
                   gamma = IMG_GAMMA, palette = viridis(256),
                   main = "curated-signal (tissue richness)", add_colorbar = TRUE)
  render_ion_image(build_img(x, y, as.integer(is_tissue[sel]) + 1L), clip_hi = 1.0,
                   palette = c("grey30", "white"),
                   main = "on-tissue pixels (curated-signal cut)",
                   add_colorbar = FALSE)
  mtext(sprintf("%s - on-tissue delineation (348 curated ions)", sid), outer = TRUE, font = 2)

  par(mfrow = c(1, 3), oma = c(1, 1, 3, 1))
  for (s in S_GRID)
    render_cluster_panel(x, y, k10_g[[paste0("s", s)]][sel], sprintf("k10 s=%d (local)", s),
                         k = K_PASS2)
  mtext(sprintf("%s - Pass 2 SSC k=10 within tissue - LOCAL cluster maps by s", sid),
        outer = TRUE, font = 2)

  # Pass-2 colored by harmonized CHEMOTYPE (shared palette across sections)
  par(mfrow = c(1, 1), oma = c(1, 1, 3, 1))
  render_ion_image(build_img(x, y, chemotype_px[sel]), clip_hi = 1.0, palette = chemo_pal,
                   main = "Pass 2 colored by CHEMOTYPE (shared across sections)",
                   add_colorbar = FALSE)
  legend("topright", legend = sprintf("chemo %d", seq_len(K_chemo)),
         fill = chemo_pal, bty = "n", cex = 0.8, xpd = NA)
  mtext(sprintf("%s - chemotype map", sid), outer = TRUE, font = 2)

  dev.off()
  cat(sprintf("[20] Wrote PDF: ssc_mask_%s.pdf\n", sid))
}

# ===========================================================================
# 8. Global chemotype harmonization PDF
# ===========================================================================
if (!is.null(hc)) {
  pdf(file.path(FIG_SSC, "chemotype_harmonization.pdf"), width = 12, height = 9)
  # dendrogram with chemotype cut
  par(mar = c(8, 4, 3, 1))
  plot(hc, labels = cent_keys, main = sprintf("Per-section centroid dendrogram (K=%d chemotypes)", K_chemo),
       xlab = "", sub = "", cex = 0.7)
  rect.hclust(hc, k = K_chemo, border = chemo_pal)
  # correlation heatmap (ordered by dendrogram)
  par(mar = c(8, 8, 3, 4))
  ord <- hc$order
  image(seq_len(ncent), seq_len(ncent), cm[ord, ord], col = viridis(256),
        axes = FALSE, xlab = "", ylab = "", main = "Centroid Pearson correlation (log1p)")
  axis(1, at = seq_len(ncent), labels = cent_keys[ord], las = 2, cex.axis = 0.6)
  axis(2, at = seq_len(ncent), labels = cent_keys[ord], las = 2, cex.axis = 0.6)
  # chemotype composition by section
  par(mar = c(7, 4, 3, 6))
  comp <- sapply(sample_ids, function(s)
    sapply(seq_len(K_chemo), function(ch)
      sum(chemotype_px == ch & as.character(pc$sample_id) == s, na.rm = TRUE)))
  barplot(comp, col = chemo_pal, las = 2, ylab = "tissue pixels",
          main = "Chemotype composition by section", legend.text = sprintf("chemo %d", seq_len(K_chemo)),
          args.legend = list(x = "topright", bty = "n", cex = 0.7, inset = c(-0.12, 0), xpd = NA))
  dev.off()
  cat("[20] Wrote PDF: chemotype_harmonization.pdf\n")
}

# ===========================================================================
# 9. Log
# ===========================================================================
elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
log_df <- data.frame(
  n_pixels = npix, n_features = nrow(clean), n_sections = length(sample_ids),
  n_centroids = ncent, K_chemotypes = K_chemo,
  n_tissue_px = sum(is_tissue), bpparam = bpused,
  elapsed_min = round(elapsed, 2), stringsAsFactors = FALSE)
write.csv(log_df, OUT_LOG, row.names = FALSE)
print(log_df)
cat(sprintf("[20] DONE in %.1f min.\n", elapsed))
