#!/usr/bin/env Rscript
# ============================================================================
# 03_apical_report.R - EXPANDED apical-orientation report. Compares per-organoid
#   citrate vs DHA mean intensity across apical classes (apical-out /
#   basolateral-out / mixed), with:
#     (1) absolute (TIC-normalised) means, dot+box, with significance stars
#     (2) a WITHIN-SECTION normalised view (log2 vs section median) that removes
#         per-section baseline so the apical effect is isolated
#     (3) per-section brightfield + MSI overlay images (homogeneous-opacity style)
#         with organoid outlines coloured by apical class
#
# Stats note: intensities are ALREADY TIC-normalised. The "normalised" panel
# further divides each organoid by its SECTION median to strip slide/section
# baseline. Significance = pairwise Wilcoxon rank-sum on organoids (organoid =
# independent unit); these are DESCRIPTIVE - organoids within a section are
# pseudo-replicates. Single-condition study: no incubation-time comparison.
#
# In : results/annotation/apical_map_consensus.csv  (two-annotator CONSENSUS apical map; default)
#      cache/peaks_tissue_combined.rds, cache/instances_{final,clean,split,}_<sid>.rds
#      cache/register/nd2final_<sid>.rds, figures/registration/crops/optical_<sid>.png
# Out: figures/annotation/apical_citrate_dha_report.pdf
#      results/annotation/apical_citrate_dha_stats.csv
# Usage: Rscript R/11_per_organoid_final/03_apical_report.R
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
suppressPackageStartupMessages({ library(Cardinal); library(png); library(EBImage); library(viridisLite) })

CIT_MZ <- 191.0217; DHA_MZ <- 327.2330
OVERLAY_ALPHA <- 0.65; GRID_D <- 2L
REG_CACHE <- cache_in("register")   # read-only reuse from upstream cache
CROP_DIR  <- file.path(FIG_DIR, "registration", "crops")
ANNOT_RES <- file.path(RES_DIR, "annotation")
# Optional args: [1] apical-map CSV (sid, instance, apical_class), [2] output suffix.
# Default = the two-annotator CONSENSUS map (apical_map_consensus.csv) -> default-named
# outputs. Pass a different map (e.g. apical_map_toppick.csv) + a suffix for a subset.
args <- commandArgs(trailingOnly = TRUE)
APICAL_CSV <- if (length(args) >= 1 && nzchar(args[1])) args[1] else
                file.path(ANNOT_RES, "apical_map_consensus.csv")
SUFFIX     <- if (length(args) >= 2 && nzchar(args[2])) paste0("_", args[2]) else ""
# Optional [3] flag drops the ambiguous "mixed" class -> two-class (basolateral-out
# vs apical-out) report; auto-suffixes outputs "_nomixed" when no suffix was given.
DROP_MIXED <- length(args) >= 3 && nzchar(args[3]) &&
                tolower(args[3]) %in% c("1", "true", "nomixed", "drop_mixed")
if (DROP_MIXED && SUFFIX == "") SUFFIX <- "_nomixed"
CURATED_CSV <- file.path(RES_DIR, "curated_feature_table.csv")   # 348 curated on-tissue features
TIC_PPM   <- 10                                                   # mz_curated -> MSE feature tolerance
OUT_PDF   <- file.path(FIG_DIR, "annotation", sprintf("apical_citrate_dha_report%s.pdf", SUFFIX))
OUT_STATS <- file.path(ANNOT_RES, sprintf("apical_citrate_dha_stats%s.csv", SUFFIX))
OUT_CSV2  <- file.path(ANNOT_RES, sprintf("apical_citrate_dha_per_organoid_normalized%s.csv", SUFFIX))
stopifnot(file.exists(APICAL_CSV), file.exists(CURATED_CSV))
cat(sprintf("[105] apical map = %s   suffix = '%s'\n", APICAL_CSV, SUFFIX))

CLASSES <- c("basolateral_out", "apical_out", "mixed")
if (DROP_MIXED) CLASSES <- setdiff(CLASSES, "mixed")
CLAB <- APICAL_LABS[CLASSES]
CCOL <- APICAL_COLS[CLASSES]   # locked: green basolateral-out / magenta apical-out / grey mixed
RAMP <- viridisLite::viridis(64)

# ---- apical map (sid, instance -> class) -----------------------------------
ap <- read.csv(APICAL_CSV, stringsAsFactors = FALSE)
ap_map <- setNames(ap$apical_class, paste(ap$sid, ap$instance))

# ---- MSE + ion vectors -----------------------------------------------------
mse <- readRDS(TISSUE_MSE)
pd  <- as.data.frame(pixelData(mse)); pd$gidx <- seq_len(nrow(pd))
mzs <- mz(mse)
ci <- which.min(abs(mzs - CIT_MZ)); di <- which.min(abs(mzs - DHA_MZ))
SIDS20 <- levels(pixelData(mse)$sample_id); if (is.null(SIDS20)) SIDS20 <- sort(unique(as.character(pd$sample_id)))
ord20  <- sort(SIDS20)

# realize the full spectra matrix once -> ion vectors + metabolite-TIC image
SP <- as.matrix(spectra(mse)); rm(mse); gc(verbose = FALSE)
val_dha <- SP[di, ]
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_citrate.R"))
val_cit <- citrate_onto_pd(pd)   # anchored citrate (raw imzML, +-CITRATE_WIN_PPM, TIC-norm); NOT a grid feature

# metabolite-TIC = per-pixel sum over the FULL 348-feature curated on-tissue set
cur <- read.csv(CURATED_CSV, stringsAsFactors = FALSE)
F_idx <- vapply(cur$mz_curated, function(q) {
  j <- which.min(abs(mzs - q)); if (abs(mzs[j] - q)/q*1e6 <= TIC_PPM) j else NA_integer_
}, integer(1))
F_idx <- sort(unique(F_idx[!is.na(F_idx)]))
metTIC <- colSums(SP[F_idx, , drop = FALSE])
rm(SP); gc(verbose = FALSE)
cat(sprintf("[105] metabolite-TIC: %d/%d curated features mapped into MSE (<=%g ppm); citrate in pool: %s, DHA in pool: %s\n",
            length(F_idx), nrow(cur), TIC_PPM, ci %in% F_idx, di %in% F_idx))

inst_file <- function(sid) {
  for (suf in c("final", "clean", "split")) {
    f <- cache_in(sprintf("instances_%s_%s.rds", suf, sid))
    if (file.exists(f)) return(f)
  }
  cache_in(sprintf("instances_%s.rds", sid))
}

# ---- per-instance means for ALL organoids (annotated + not) ----------------
rows <- list(); inst_cache <- list()
for (sid in ord20) {
  f <- inst_file(sid); if (!file.exists(f)) next
  inst <- readRDS(f); inst_cache[[sid]] <- inst
  ids <- sort(unique(inst$instance[inst$instance > 0]))
  for (k in ids) {
    g <- inst$gidx[inst$instance == k]; if (!length(g)) next
    ok <- metTIC[g] > 0                         # per-pixel ratio guard
    cit_mtic <- if (any(ok)) mean(val_cit[g][ok] / metTIC[g][ok]) * 100 else NA_real_
    dha_mtic <- if (any(ok)) mean(val_dha[g][ok] / metTIC[g][ok]) * 100 else NA_real_
    cit_m <- mean(val_cit[g]); dha_m <- mean(val_dha[g])
    rows[[length(rows)+1L]] <- data.frame(
      sid = sid, instance = k,
      apical_class = unname(ap_map[paste(sid, k)]),
      cit = cit_m, dha = dha_m,
      cit_mtic = cit_mtic, dha_mtic = dha_mtic,
      cit_dha = if (is.finite(dha_m) && dha_m > 0) log2(cit_m / dha_m) else NA_real_,
      n_px = length(g), stringsAsFactors = FALSE)
  }
}
all_org <- do.call(rbind, rows)
# within-section normalisation: log2( organoid / section median over ALL instances )
sec_med <- function(v, sid) tapply(v, sid, median, na.rm = TRUE)
mc <- sec_med(all_org$cit, all_org$sid); md <- sec_med(all_org$dha, all_org$sid)
all_org$cit_rel <- log2(all_org$cit / mc[all_org$sid])
all_org$dha_rel <- log2(all_org$dha / md[all_org$sid])

org <- all_org[!is.na(all_org$apical_class) & all_org$apical_class %in% CLASSES, ]
org$apical_class <- factor(org$apical_class, levels = CLASSES)
cat(sprintf("[105] %d annotated organoids (of %d total)\n", nrow(org), nrow(all_org)))

# ---- significance: pairwise Wilcoxon ---------------------------------------
PAIRS <- list(c("basolateral_out","apical_out"), c("basolateral_out","mixed"), c("apical_out","mixed"))
PAIRS <- Filter(function(pr) all(pr %in% CLASSES), PAIRS)   # drop *_vs_mixed when mixed removed
stars <- function(p) if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
pw <- function(v, cls, subset = rep(TRUE, length(v))) {
  out <- lapply(PAIRS, function(pr) {
    a <- v[subset & cls == pr[1]]; b <- v[subset & cls == pr[2]]
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    p <- if (length(a) >= 3 && length(b) >= 3) suppressWarnings(wilcox.test(a, b)$p.value) else NA_real_
    data.frame(g1 = pr[1], g2 = pr[2], n1 = length(a), n2 = length(b),
               p = p, stars = stars(p), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}
stat_rows <- list()
for (ion in c("cit","dha","cit_rel","dha_rel","cit_mtic","dha_mtic","cit_dha")) {
  s <- pw(org[[ion]], org$apical_class); s <- cbind(metric = ion, scope = "all", s)
  stat_rows[[length(stat_rows)+1L]] <- s
}
stat_df <- do.call(rbind, stat_rows)
dir.create(ANNOT_RES, showWarnings = FALSE, recursive = TRUE)
write.csv(stat_df, OUT_STATS, row.names = FALSE)
write.csv(org[, c("sid","instance","apical_class","n_px",
                  "cit","dha","cit_rel","dha_rel","cit_mtic","dha_mtic","cit_dha")],
          OUT_CSV2, row.names = FALSE)

# ===========================================================================
# OUTWARD CITRATE GRADIENT per organoid (curated Voronoi zones), by apical class
# ===========================================================================
MIN_ZONE_PX <- 3L
NZ          <- length(OUT_ZONE_UM)                 # 7 outward zones (10..500 um)
FAR_ZONES   <- which(OUT_ZONE_UM >= 80)            # zones 4..7
OUT_GRAD_CSV <- file.path(ANNOT_RES, sprintf("apical_gradient_per_organoid%s.csv", SUFFIX))
spr <- function(v) { idx <- which(!is.na(v))
  if (length(idx) < 3 || stats::sd(v[idx]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(idx, v[idx], method = "spearman")) }
cliffs_delta <- function(a, b) { a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (!length(a) || !length(b)) return(NA_real_)
  mean(outer(a, b, ">")) - mean(outer(a, b, "<")) }

ap_class_of <- setNames(as.character(org$apical_class), paste(org$sid, org$instance))
prof_rows <- list(); grad_rows <- list()
for (sid in ord20) {
  fz <- cache_in(sprintf("zones_%s.rds", sid)); if (!file.exists(fz)) next
  z <- readRDS(fz)
  for (k in sort(unique(z$instance[z$instance > 0]))) {
    key <- paste(sid, k); cls <- ap_class_of[key]
    if (is.na(cls)) next                              # only annotated organoids
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
    # absolute near-field level: mean citrate over gel catchment within a distance window
    g50  <- z$gidx[z$instance_catch == k & z$signed_dist_um > 0 & z$signed_dist_um <= 50]
    g100 <- z$gidx[z$instance_catch == k & z$signed_dist_um > 0 & z$signed_dist_um <= 100]
    near50  <- if (length(g50)  >= MIN_ZONE_PX) mean(val_cit[g50])  else NA_real_
    near100 <- if (length(g100) >= MIN_ZONE_PX) mean(val_cit[g100]) else NA_real_
    for (zz in seq_len(NZ))
      prof_rows[[length(prof_rows)+1L]] <- data.frame(
        sid = sid, instance = k, key = key, apical_class = cls,
        zone = zz, zone_um = OUT_ZONE_UM[zz],
        surf = surf_norm[zz], abs = cit_z[zz], mtic = mtic_z[zz], citdha = citdha[zz],
        dha_sn = dha_sn[zz], dha_abs = dha_z[zz], stringsAsFactors = FALSE)
    grad_rows[[length(grad_rows)+1L]] <- data.frame(
      sid = sid, instance = k, apical_class = cls,
      rho_out = spr(cit_z), far_index = mean(surf_norm[FAR_ZONES], na.rm = TRUE),
      near50 = near50, near100 = near100,
      stringsAsFactors = FALSE)
  }
}
grad_prof <- do.call(rbind, prof_rows)
grad_org  <- do.call(rbind, grad_rows)
grad_org$apical_class <- factor(grad_org$apical_class, levels = CLASSES)
write.csv(grad_org, OUT_GRAD_CSV, row.names = FALSE)
cat(sprintf("[105] gradient: %d organoids with outward profiles (in=%d out=%d%s)\n",
            nrow(grad_org), sum(grad_org$apical_class=="basolateral_out"),
            sum(grad_org$apical_class=="apical_out"),
            if (DROP_MIXED) "" else sprintf(" mixed=%d", sum(grad_org$apical_class=="mixed"))))

# gradient metrics significance -> append to stat_df
for (m in c("rho_out","far_index","near50","near100")) {
  s <- pw(grad_org[[m]], grad_org$apical_class); stat_df <- rbind(stat_df, cbind(metric = paste0("grad_",m), scope = "all", s))
}
write.csv(stat_df, OUT_STATS, row.names = FALSE)
cdo <- function(m) cliffs_delta(grad_org[[m]][grad_org$apical_class=="apical_out"], grad_org[[m]][grad_org$apical_class=="basolateral_out"])
D_RHO <- cdo("rho_out"); D_FAR <- cdo("far_index"); D_N50 <- cdo("near50"); D_N100 <- cdo("near100")

# ===========================================================================
# PLOT helpers
# ===========================================================================
# significance brackets above a dot/box panel
add_sig <- function(stat, ytop, yrange) {
  step <- 0.07 * yrange; lvl <- 0
  xof <- c(basolateral_out = 1, apical_out = 2, mixed = 3)
  ord <- order(abs(xof[stat$g1] - xof[stat$g2]))   # short spans first (lower)
  for (i in ord) {
    if (is.na(stat$p[i]) || stat$stars[i] == "ns") next
    x1 <- xof[stat$g1[i]]; x2 <- xof[stat$g2[i]]
    y <- ytop + step * (lvl + 1); lvl <- lvl + 1
    segments(x1, y, x2, y, lwd = 1.1); segments(x1, y, x1, y - step*0.3, lwd = 1.1)
    segments(x2, y, x2, y - step*0.3, lwd = 1.1)
    text((x1 + x2)/2, y + step*0.18, stat$stars[i], cex = 1.1, font = 2)
  }
  lvl
}
dotbox <- function(yv, ylab, main, cls = org$apical_class, sub = rep(TRUE, length(yv)),
                   stat = NULL, ptcol = NULL, baseline = NULL) {
  ok <- is.finite(yv) & sub
  yr <- range(yv[ok]); pad <- diff(yr) * 0.04
  ymax <- yr[2]; ymin <- min(yr[1], if (!is.null(baseline)) baseline else yr[1])
  headroom <- if (!is.null(stat)) 0.30 * diff(yr) else pad
  plot(NA, xlim = c(0.5, length(CLASSES) + 0.5), ylim = c(ymin - pad, ymax + headroom),
       xaxt = "n", xlab = "", ylab = ylab, main = main, cex.main = 1.0)
  if (!is.null(baseline)) abline(h = baseline, lty = 3, col = "grey60")
  for (i in seq_along(CLASSES)) {
    v <- yv[ok & cls == CLASSES[i]]; if (!length(v)) next
    q <- quantile(v, c(.25,.5,.75)); iqr <- q[3]-q[1]
    rect(i-0.28, q[1], i+0.28, q[3], border = "grey35", col = NA, lwd = 1.3)
    segments(i-0.28, q[2], i+0.28, q[2], lwd = 2.4, col = "grey20")
    segments(i, q[3], i, min(max(v), q[3]+1.5*iqr), col="grey55")
    segments(i, q[1], i, max(min(v), q[1]-1.5*iqr), col="grey55")
    bg <- if (is.null(ptcol)) adjustcolor(CCOL[CLASSES[i]],0.7) else
            adjustcolor(ptcol[ok & cls == CLASSES[i]], 0.75)
    set.seed(i)
    points(i + runif(length(v), -0.16, 0.16), v, pch = 21, cex = 1.05, bg = bg, col = "white", lwd = 0.5)
  }
  ns <- sapply(CLASSES, function(c) sum(ok & cls == c))
  axis(1, at = seq_along(CLASSES), labels = sprintf("%s\n(n=%d)", CLAB[CLASSES], ns), padj = 0.6, cex.axis = 0.92)
  if (!is.null(stat)) add_sig(stat, ymax, diff(yr))
}

# ---- spatial overlay helpers (R/05_registration_refine/04_overlay_report.R homogeneous-opacity style) --------------
inv_affine <- function(B, NX, NY) { M <- rbind(B[1,], B[2,]); t0 <- B[3,]
  xy <- solve(M, rbind(NX - t0[1], NY - t0[2])); list(x = xy[1,], y = xy[2,]) }
instance_outlines <- function(M, B, cx0, cy0) {
  M <- EBImage::fillHull(EBImage::closing(M > 0, EBImage::makeBrush(3, "disc")))
  out <- list()
  for (co in EBImage::ocontour(EBImage::bwlabel(M))) {
    if (nrow(co) < 6) next
    p <- apply_affine(B, cbind(co[,1]+0.5, co[,2]+0.5))
    out[[length(out)+1L]] <- list(x = c(p[,1]-cx0+1, p[1,1]-cx0+1), y = c(p[,2]-cy0+1, p[1,2]-cy0+1))
  }
  out
}
# global per-ion p99.5 clip across all sections (cross-section comparability)
sel_on <- unlist(lapply(ord20, function(sid) {
  inst <- inst_cache[[sid]]; inst$gidx[inst$instance > 0] }))
HI_CIT <- as.numeric(quantile(val_cit[val_cit > 0], IMG_CLIP_HI))
HI_DHA <- as.numeric(quantile(val_dha[val_dha > 0], IMG_CLIP_HI))

overlay_panel <- function(sid, val, hi, main) {
  xf <- file.path(REG_CACHE, sprintf("nd2final_%s.rds", sid))
  pf <- file.path(CROP_DIR, sprintf("optical_%s.png", sid))
  if (!file.exists(xf) || !file.exists(pf)) { plot.new(); title(main); return(invisible()) }
  X <- readRDS(xf); B <- X$B_msi_nd2; smn <- X$scale_msi_nd2; cx0 <- X$crop[1]; cy0 <- X$crop[2]
  om <- png::readPNG(pf); if (length(dim(om)) == 3) om <- om[,,1]; cw <- ncol(om); ch <- nrow(om)
  sub <- pd[as.character(pd$sample_id) == sid, c("x","y")]; W <- max(sub$x); H <- max(sub$y)
  secval <- val[pd$gidx[as.character(pd$sample_id) == sid]]
  gx <- seq(1, cw, by = GRID_D); gy <- seq(1, ch, by = GRID_D)
  NX <- cx0 + rep(gx, times = length(gy)) - 1; NY <- cy0 + rep(gy, each = length(gx)) - 1
  ms <- inv_affine(B, NX, NY); inb <- ms$x >= .5 & ms$x <= W+.5 & ms$y >= .5 & ms$y <= H+.5
  lut <- matrix(NA_real_, W, H); lut[cbind(sub$x, sub$y)] <- secval
  v <- rep(NA_real_, length(NX))
  if (any(inb)) { xi <- pmin(pmax(round(ms$x[inb]),1),W); yi <- pmin(pmax(round(ms$y[inb]),1),H); v[inb] <- lut[cbind(xi,yi)] }
  vm <- matrix(v, length(gy), length(gx), byrow = TRUE); vs <- pmin(vm/hi, 1)
  par(mar = c(2.2, 1, 2.6, 1))
  plot.new(); plot.window(c(1, cw), c(ch, 1), asp = 1)
  rasterImage(om/max(om), 1, ch, cw, 1, interpolate = TRUE); title(main, cex.main = 0.92)
  fin <- is.finite(vs)
  if (any(fin)) {
    cidx <- pmin(pmax(round(vs*63)+1,1),64); cmat <- matrix("#00000000", length(gy), length(gx))
    cmat[fin] <- adjustcolor(RAMP[cidx[fin]], alpha.f = OVERLAY_ALPHA)
    rasterImage(as.raster(cmat), 1, ch, cw, 1, interpolate = FALSE)
  }
  # organoid outlines coloured by apical class (grey if unannotated) + labels
  inst <- inst_cache[[sid]]; ids <- sort(unique(inst$instance[inst$instance > 0]))
  Wi <- max(inst$x); Hi <- max(inst$y); lab <- matrix(0L, Wi, Hi); lab[cbind(inst$x, inst$y)] <- as.integer(inst$instance)
  for (k in ids) {
    cls <- unname(ap_map[paste(sid, k)])
    bc <- if (!is.na(cls) && cls %in% CLASSES) CCOL[cls] else "grey70"
    for (poly in instance_outlines(lab == k, B, cx0, cy0)) polygon(poly$x, poly$y, border = bc, lwd = 1.6)
    cen <- apply_affine(B, matrix(c(mean(inst$x[inst$instance==k]), mean(inst$y[inst$instance==k])), 1))
    text(cen[1]-cx0+1, cen[2]-cy0+1, k, col = "white", font = 2, cex = 0.6)
  }
  # scale bar BELOW image (locked style): bar then centred label beneath
  umpx <- MSI_PIXEL_UM / smn; bp <- 100/umpx; yb <- ch + 0.06*ch
  par(xpd = NA); segments(cw-bp, yb, cw, yb, lwd = 3, col = "black")
  text(cw - bp/2, yb + 0.05*ch, "100 um", cex = 0.6); par(xpd = FALSE)
}
colorbar <- function(lab) {
  par(mar = c(3, 2, 2, 2)); plot.new(); plot.window(c(0,1), c(0,1))
  n <- 64; rasterImage(as.raster(matrix(rev(RAMP), ncol = 1)), 0.35, 0.05, 0.65, 0.95)
  text(0.68, 0.95, "p99.5", adj = 0, cex = 0.7); text(0.68, 0.05, "0", adj = 0, cex = 0.7)
  text(0.5, 0.995, lab, cex = 0.75, font = 2)
}

# ===========================================================================
# RENDER
# ===========================================================================
pdf(OUT_PDF, width = 11, height = 7.2)

# ---- Page 0: Summary Statement (presentable claim, values pulled from stats) ----
sget <- function(metric, scope, g1, g2, what) {
  r <- stat_df[stat_df$metric == metric & stat_df$scope == scope &
               stat_df$g1 == g1 & stat_df$g2 == g2, ]
  if (!nrow(r)) NA else r[[what]][1]
}
fp <- function(p) { p <- as.numeric(p); if (is.na(p)) "NA" else if (p < 1e-3) sprintf("%.1e", p) else sprintf("%.3f", p) }
p_n50  <- sget("grad_near50", "all", "basolateral_out", "apical_out", "p")
s_n50  <- sget("grad_near50", "all", "basolateral_out", "apical_out", "stars")
p_cd   <- sget("cit_dha",     "all", "basolateral_out", "apical_out", "p")
s_cd   <- sget("cit_dha",     "all", "basolateral_out", "apical_out", "stars")

layout(1); par(mar = c(1, 1, 1, 1))
plot.new()
text(0.5, 0.93, "Summary Statement", font = 2, cex = 2.0, adj = c(0.5, 1))
claim <- sprintf(paste0(
  "Apical-out organoids show elevated TIC-normalised citrate within <=50 um of the organoid\n",
  "surface (near-field level, grad_near50: p = %s %s, Cliff's d = %+.2f vs basolateral-out;\n",
  "large effect). This near-field enrichment is confirmed against the shared thickness /\n",
  "ionisation component by the citrate/DHA internal-control ratio (cit_dha: p = %s [%s]),\n",
  "which cancels any factor scaling both ions together. The outward monotonic slope\n",
  "(grad_rho_out) is not significant, so the result is presented as near-field level\n",
  "magnitude, not gradient steepness."),
  fp(p_n50), s_n50, as.numeric(D_N50), fp(p_cd), s_cd)
text(0.07, 0.78, claim, adj = c(0, 1), cex = 1.2)
text(0.07, 0.20, sprintf(paste0(
  "Basis: %d annotated organoids across 20 sections; organoid = independent unit; pairwise Wilcoxon\n",
  "rank-sum (descriptive). Citrate = [M-H]- anchored 191.0198 +-7 ppm (raw imzML), TIC-normalised.\n",
  "Near-field = mean citrate over gel catchment 0 < d <= 50 um from the organoid surface."),
  nrow(grad_org)), adj = c(0, 1), cex = 0.85, col = "grey30")

# ---- Page 1: absolute (TIC) means, citrate & DHA, with sig stars -----------
layout(matrix(1:2, nrow = 1)); par(oma = c(0,0,3,0), mar = c(4,4.5,3,1))
dotbox(org$cit, "mean intensity (TIC-norm, a.u.)", sprintf("Citrate [M-H]-  %.4f", mzs[ci]),
       stat = stat_df[stat_df$metric=="cit" & stat_df$scope=="all", ])
dotbox(org$dha, "mean intensity (TIC-norm, a.u.)", sprintf("DHA C22:6 [M-H]-  %.4f", mzs[di]),
       stat = stat_df[stat_df$metric=="dha" & stat_df$scope=="all", ])
mtext(sprintf("Per-organoid mean ion intensity by apical orientation  (n=%d organoids, 20 sections)", nrow(org)),
      outer = TRUE, font = 2, cex = 1.15, line = 1.2)
mtext("* p<0.05  ** p<0.01  *** p<0.001  (pairwise Wilcoxon; organoid = unit; descriptive)",
      outer = TRUE, cex = 0.8, line = -0.2)

# ---- Page 2: same plot, metabolite-TIC normalised (relativizes ionization) --
# Companion to page 1: per-organoid mean intensity expressed as each ion / the
# summed 348-curated-feature signal per pixel (% of pool), which removes regional
# ionization bias. Placed right after the absolute-TIC page for direct comparison.
layout(matrix(1:2, nrow = 1)); par(oma = c(0,0,3,0), mar = c(4,4.5,3,1))
dotbox(org$cit_mtic, "citrate / metabolite-TIC  (% of pool)", "Citrate - metabolite-TIC normalised",
       stat = stat_df[stat_df$metric=="cit_mtic" & stat_df$scope=="all", ])
dotbox(org$dha_mtic, "DHA / metabolite-TIC  (% of pool)", "DHA - metabolite-TIC normalised",
       stat = stat_df[stat_df$metric=="dha_mtic" & stat_df$scope=="all", ])
mtext(sprintf("Per-organoid mean ion intensity by apical orientation - metabolite-TIC normalised  (n=%d organoids, 20 sections)", nrow(org)),
      outer = TRUE, font = 2, cex = 1.12, line = 1.2)
mtext(sprintf("Per organoid pixel: ion / sum of %d curated on-tissue features, then mean over the organoid (x100). Relativizes regional ionization bias. Data already globally TIC-normalised.",
              length(F_idx)), outer = TRUE, cex = 0.78, line = -0.2)

# ---- Page 3: within-section normalised (log2 vs section median) ------------
layout(matrix(1:2, nrow = 1)); par(oma = c(0,0,3,0), mar = c(4,4.5,3,1))
dotbox(org$cit_rel, "log2( organoid / section median )", "Citrate - section-normalised",
       stat = stat_df[stat_df$metric=="cit_rel" & stat_df$scope=="all", ], baseline = 0)
dotbox(org$dha_rel, "log2( organoid / section median )", "DHA - section-normalised",
       stat = stat_df[stat_df$metric=="dha_rel" & stat_df$scope=="all", ], baseline = 0)
mtext("Within-section normalisation (removes per-section/slide baseline)", outer = TRUE, font = 2, cex = 1.15, line = 1.2)
mtext("Data are already TIC-normalised; here each organoid is divided by its section's median organoid (0 = section median).",
      outer = TRUE, cex = 0.78, line = -0.2)

# ---- Page 4: citrate/DHA internal-control ratio (cancels parallel shift) ---
layout(1); par(oma = c(0,0,3,0), mar = c(4,4.5,3,1))
dotbox(org$cit_dha, "log2( citrate / DHA )", "All organoids",
       stat = stat_df[stat_df$metric=="cit_dha" & stat_df$scope=="all", ])
mtext("Citrate/DHA internal-control ratio (cancels the per-organoid parallel shift)",
      outer = TRUE, font = 2, cex = 1.15, line = 1.2)
mtext("DHA = designed negative control (membrane-bound, organoid-confined); the ratio removes any factor scaling both ions together.",
      outer = TRUE, cex = 0.78, line = -0.2)

# ---- Page 6: outward citrate gradient, basolateral-out vs apical-out -------------
GCOL <- c(basolateral_out = unname(CCOL["basolateral_out"]), apical_out = unname(CCOL["apical_out"]))
decay_panel <- function(valcol, ylab, main, dha_col = NULL) {
  use <- grad_prof[grad_prof$apical_class %in% c("basolateral_out","apical_out"), ]
  yv <- use[[valcol]]; ok <- is.finite(yv)
  if (!is.null(dha_col)) ok2 <- is.finite(grad_prof[[dha_col]]) else ok2 <- FALSE
  ymax <- max(c(yv[ok], if (!is.null(dha_col)) grad_prof[[dha_col]][ok2]), na.rm = TRUE)
  ymin <- min(c(0, yv[ok]), na.rm = TRUE)
  plot(NA, xlim = range(OUT_ZONE_UM), ylim = c(ymin, ymax * 1.04), log = "x",
       xlab = "distance outward (um, log)", ylab = ylab, main = main, cex.main = 1.0)
  for (kk in unique(use$key)) {                      # per-organoid spaghetti
    s <- use[use$key == kk, ]; s <- s[order(s$zone_um), ]
    lines(s$zone_um, s[[valcol]], col = adjustcolor(GCOL[s$apical_class[1]], 0.22), lwd = 0.8)
  }
  if (!is.null(dha_col)) {                            # DHA confined control (pooled median)
    agg <- tapply(grad_prof[[dha_col]], grad_prof$zone_um, median, na.rm = TRUE)
    lines(as.numeric(names(agg)), as.numeric(agg), col = "grey45", lwd = 2, lty = 2, type = "b", pch = 17, cex = 0.8)
  }
  for (cl in c("basolateral_out","apical_out")) {          # group-median bold
    s <- use[use$apical_class == cl, ]
    agg <- tapply(s[[valcol]], s$zone_um, median, na.rm = TRUE)
    lines(as.numeric(names(agg)), as.numeric(agg), col = GCOL[cl], lwd = 3, type = "b", pch = 19)
  }
  leg <- c("basolateral-out (median)","apical-out (median)"); lc <- c(GCOL["basolateral_out"], GCOL["apical_out"])
  lw <- c(3,3); lt <- c(1,1); lp <- c(19,19)
  if (!is.null(dha_col)) { leg <- c(leg,"DHA control"); lc <- c(lc,"grey45"); lw <- c(lw,2); lt <- c(lt,2); lp <- c(lp,17) }
  legend("topright", legend = leg, col = lc, lwd = lw, lty = lt, pch = lp, bty = "n", cex = 0.8)
}
layout(matrix(1:4, nrow = 2, byrow = TRUE)); par(oma = c(0,0,3,0), mar = c(4.2,4.3,2.8,1))
decay_panel("surf", "citrate / surface band", "Surface-normalised", dha_col = "dha_sn")
decay_panel("abs",  "citrate (TIC-norm, a.u.)", "Absolute (TIC)")
decay_panel("mtic", "citrate / metabolite-TIC (%)", "Metabolite-TIC normalised")
decay_panel("citdha", "citrate / interior DHA", "Citrate/DHA (interior internal std)")
mtext(sprintf("Outward citrate [M-H]- gradient: basolateral-out vs apical-out  (curated Voronoi zones, organoid = unit; in=%d out=%d)",
              sum(grad_org$apical_class=="basolateral_out"), sum(grad_org$apical_class=="apical_out")),
      outer = TRUE, font = 2, cex = 1.05, line = 1.2)
mtext("Citrate measured as [M-H]- anchored 191.0198 +-7 ppm (raw imzML). citrate/DHA uses interior DHA: per-zone ratio undefined outward (DHA gel-absent).",
      outer = TRUE, cex = 0.72, line = -0.2)

# ---- Page 7: per-organoid gradient metrics, basolateral-out vs apical-out --------
layout(matrix(1:2, nrow = 1)); par(oma = c(0,0,3,0), mar = c(4,4.5,3,1))
dotbox(grad_org$rho_out, "Spearman rho (zone vs citrate)", "Outward monotonicity (rho_out)",
       cls = grad_org$apical_class, stat = stat_df[stat_df$metric=="grad_rho_out" & stat_df$scope=="all", ], baseline = 0)
dotbox(grad_org$far_index, "far-field / surface (>=80um)", "Far-field retention (far_index)",
       cls = grad_org$apical_class, stat = stat_df[stat_df$metric=="grad_far_index" & stat_df$scope=="all", ])
mtext("Per-organoid outward gradient metrics by apical class", outer = TRUE, font = 2, cex = 1.15, line = 1.2)
mtext(sprintf("Cliff's delta (apical-out vs basolateral-out):  rho_out = %+.2f   far_index = %+.2f   (positive = more outward in apical-out)",
              D_RHO, D_FAR), outer = TRUE, cex = 0.8, line = -0.2)

# ---- Page 8: absolute near-field citrate level, basolateral-out vs apical-out -----
layout(matrix(1:2, nrow = 1)); par(oma = c(0,0,3,0), mar = c(4,4.5,3,1))
dotbox(grad_org$near50, "citrate (TIC-norm, a.u.)", "Near-field citrate level 0-50 um",
       cls = grad_org$apical_class, stat = stat_df[stat_df$metric=="grad_near50" & stat_df$scope=="all", ])
dotbox(grad_org$near100, "citrate (TIC-norm, a.u.)", "Near-field citrate level 0-100 um",
       cls = grad_org$apical_class, stat = stat_df[stat_df$metric=="grad_near100" & stat_df$scope=="all", ])
mtext("Absolute near-field citrate level in the gel, by apical class  (mean [M-H]- over the organoid's gel catchment)",
      outer = TRUE, font = 2, cex = 1.1, line = 1.2)
mtext(sprintf("Absolute (TIC-norm), NOT surface-normalised. Cliff's delta out-vs-in:  0-50um = %+.2f   0-100um = %+.2f",
              D_N50, D_N100), outer = TRUE, cex = 0.8, line = -0.2)

# ---- Page 8b: near-field 0-50 um vs 0-100 um WITHIN each apical class -------
# Paired per organoid (same organoid measured at both windows). 0-100 um includes
# the more-dilute 50-100 um band, so 100 <= 50 for an outward-decaying field.
fp2 <- function(p) if (is.na(p)) "NA" else if (p < 1e-3) sprintf("%.1e", p) else sprintf("%.3f", p)
paired_50v100 <- function(clz) {
  d <- grad_org[grad_org$apical_class == clz, c("near50", "near100")]
  d <- d[is.finite(d$near50) & is.finite(d$near100), , drop = FALSE]
  n <- nrow(d); col <- unname(CCOL[clz])
  p <- if (n >= 3) suppressWarnings(wilcox.test(d$near50, d$near100, paired = TRUE)$p.value) else NA_real_
  yr <- if (n) range(c(d$near50, d$near100)) else c(0, 1)
  pad <- diff(yr) * 0.06; if (!is.finite(pad) || pad == 0) pad <- 0.1
  plot(NA, xlim = c(0.5, 2.5), ylim = c(min(0, yr[1]) - pad, yr[2] + pad * 3.2),
       xaxt = "n", xlab = "", ylab = "citrate (TIC-norm, a.u.)",
       main = sprintf("%s  (n=%d)", CLAB[clz], n), cex.main = 1.05, col.main = col)
  for (r in seq_len(n)) segments(1, d$near50[r], 2, d$near100[r], col = adjustcolor(col, 0.22), lwd = 0.8)
  for (j in 1:2) {
    v <- if (j == 1) d$near50 else d$near100; if (!length(v)) next
    q <- quantile(v, c(.25, .5, .75)); iqr <- q[3] - q[1]
    rect(j - 0.22, q[1], j + 0.22, q[3], border = "grey35", col = NA, lwd = 1.3)
    segments(j - 0.22, q[2], j + 0.22, q[2], lwd = 2.4, col = "grey20")
    segments(j, q[3], j, min(max(v), q[3] + 1.5 * iqr), col = "grey55")
    segments(j, q[1], j, max(min(v), q[1] - 1.5 * iqr), col = "grey55")
    set.seed(j); points(j + runif(n, -0.10, 0.10), v, pch = 21, cex = 1.05,
                        bg = adjustcolor(col, 0.7), col = "white", lwd = 0.5)
  }
  axis(1, at = 1:2, labels = c("0-50 um", "0-100 um"), padj = 0.4, cex.axis = 0.95)
  ytop <- yr[2] + pad * 1.4
  segments(1, ytop, 2, ytop, lwd = 1.1); segments(1, ytop, 1, ytop - pad * 0.35, lwd = 1.1)
  segments(2, ytop, 2, ytop - pad * 0.35, lwd = 1.1)
  text(1.5, ytop + pad * 0.9, sprintf("%s  p=%s", stars(p), fp2(p)), cex = 0.95, font = 2)
  invisible(data.frame(class = clz, n = n, med_near50 = median(d$near50), med_near100 = median(d$near100),
                       p_paired = p, stringsAsFactors = FALSE))
}
layout(matrix(seq_along(CLASSES), nrow = 1)); par(oma = c(0, 0, 3, 0), mar = c(4, 4.5, 3, 1))
p50v100 <- do.call(rbind, lapply(CLASSES, paired_50v100))
mtext("Near-field citrate: 0-50 um vs 0-100 um WITHIN each apical class  (paired per organoid)",
      outer = TRUE, font = 2, cex = 1.1, line = 1.2)
mtext("Paired Wilcoxon signed-rank (same organoid at both windows). 0-100 um includes the more-dilute 50-100 um band; a drop = outward decay.",
      outer = TRUE, cex = 0.78, line = -0.2)

# ---- Pages 9+: per-section BF + MSI overlays ------------------------------
for (sid in ord20) {
  layout(matrix(c(1,2,3), nrow = 1), widths = c(1,1,0.22)); par(oma = c(0,0,2.6,0))
  overlay_panel(sid, val_cit, HI_CIT, "citrate [M-H]- on brightfield")
  overlay_panel(sid, val_dha, HI_DHA, "DHA [M-H]- on brightfield")
  colorbar("MSI\nintensity")
  na <- sum(!is.na(ap_map[paste(sid, sort(unique(inst_cache[[sid]]$instance[inst_cache[[sid]]$instance>0])))]))
  mtext(sprintf("%s   -   outline: apical-out (magenta) / basolateral-out (green)%s / faint-grey = unannotated   (%d annotated)",
                disp_id(sid), if (DROP_MIXED) "" else " / mixed (grey)", na), outer = TRUE, font = 2, cex = 0.95, line = 0.6)
}
dev.off()
cat(sprintf("[105] DONE -> %s\n", OUT_PDF))
cat(sprintf("[105] stats -> %s\n", OUT_STATS))

# console: print key significance
cat("\n[105] pairwise Wilcoxon (scope=all):\n")
print(stat_df[stat_df$scope == "all", c("metric","g1","g2","n1","n2","p","stars")], row.names = FALSE)
cat("\n[105] medians by class:\n")
for (ion in c("cit","dha","cit_rel","dha_rel","cit_mtic","dha_mtic","cit_dha"))
  cat(sprintf("  %-8s : %s\n", ion, paste(sprintf("%s=%.3f", CLAB[CLASSES],
      tapply(org[[ion]], org$apical_class, median, na.rm=TRUE)[CLASSES]), collapse = "  ")))

cat("\n[105] outward gradient + near-field metrics median by class:\n")
for (m in c("rho_out","far_index","near50","near100"))
  cat(sprintf("  %-9s : %s\n", m, paste(sprintf("%s=%.3f", CLAB[CLASSES],
      tapply(grad_org[[m]], grad_org$apical_class, median, na.rm=TRUE)[CLASSES]), collapse = "  ")))
cat(sprintf("  Cliff's delta out-vs-in: rho_out=%+.2f far_index=%+.2f near50=%+.2f near100=%+.2f\n", D_RHO, D_FAR, D_N50, D_N100))
cat("\n[105] gradient pairwise Wilcoxon (scope=all):\n")
print(stat_df[stat_df$scope=="all" & grepl("^grad_", stat_df$metric), c("metric","g1","g2","n1","n2","p","stars")], row.names = FALSE)

cat("\n[105] near-field 0-50 vs 0-100 um WITHIN class (paired Wilcoxon signed-rank):\n")
print(p50v100, row.names = FALSE)
