#!/usr/bin/env Rscript
# ============================================================================
# 05_ssc_report.R - VISUAL report of the SSC on-tissue clustering + organoid
#   segmentation ("how the datasets were delineated"). First SSC figure in the
#   Analysis_R_Final fork (figures/ssc/ was empty).
#
#   All inputs are ALREADY CACHED in peaks_tissue_combined.rds pixelData
#   (tissueness_px, is_tissue, ssc_k4_sec, ssc_k10_sec, chemotype) + zones_<sid>.rds
#   (instance, is_surface, signed_dist_um) - NO SSC recompute. Reads through
#   cache_in() (upstream Analysis_R/cache), writes only to the fork's figures/ssc.
#
#   Page 1  = methods + composition: SSC procedure writeup, chemotype-composition
#             stacked bar across the 20 sections, on-tissue pixel counts.
#   Pages 2+ = per dataset, compact 4-panel (shared-frame style of 10_..._final.R):
#     1 tissueness (curated-ion signal, viridis) + white-dotted on-tissue outline
#     2 SSC pass-1 k=4 clusters (r=2, s=9, adaptive)
#     3 on-tissue mask (white = tissue; k4>=50% richest  U  floor80)
#     4 organoid instance segmentation (per-instance outline + surface + id labels)
#
# In : cache/peaks_tissue_combined.rds, cache/zones_<sid>.rds, cache/register/nd2final_<sid>.rds
# Out: figures/ssc/ssc_clustering_segmentation_report[_TEST_<sid>].pdf
# Usage: Rscript R/04_ssc_ontissue/05_ssc_report.R [all|test|<sid>]   (env FINAL_OUT=<name>.pdf overrides)
#        all (default) = front page + 20 dataset pages;  test = front page + one
#        representative dataset (-> ..._TEST_<sid>.pdf);  <sid> = that one dataset.
# ============================================================================

ROOT <- file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final")
source(file.path(ROOT, "R/00_lib/gradient_config.R"))
source(file.path(ROOT, "R/00_lib/lib_register.R"))
source(file.path(ROOT, "R/00_lib/lib_report_frame.R"))   # PMAR / frame_box / scalebar_bottom / colorbar_img / SCALE_UM
suppressPackageStartupMessages({ library(Cardinal); library(viridisLite) })

log_msg  <- function(...) message(sprintf("[ssc-report] %s", sprintf(...)))
SSC_FIG  <- file.path(FIG_DIR, "ssc"); dir.create(SSC_FIG, showWarnings = FALSE, recursive = TRUE)
REG_CACHE<- cache_in("register")
PX_UM    <- MSI_PIXEL_UM
pal      <- viridisLite::viridis(256)
cat_palette <- function(k) grDevices::hcl.colors(max(k, 2), palette = "Dark 3")

# ============================================================================
# SETUP
# ============================================================================
mse <- readRDS(TISSUE_MSE); pd <- as.data.frame(pixelData(mse))
need <- c("tissueness_px","is_tissue","ssc_k4_sec","chemotype")
stopifnot(all(need %in% names(pd)))
Kc <- max(pd$chemotype, na.rm = TRUE)              # harmonized chemotype count (=4)
ord20 <- sort(levels(pd$sample_id)); if (is.null(ord20)) ord20 <- sort(unique(as.character(pd$sample_id)))

# Args:  (none)|all -> full 20-dataset report;  test -> quick TEST report (front page +
#        one representative dataset);  <sid> -> that one dataset (also a TEST report).
a <- commandArgs(trailingOnly = TRUE)
mode <- if (length(a) == 0) "all" else tolower(a[1])
TEST_SID <- if ("AO_0h_sl6A_sec2a" %in% ord20) "AO_0h_sl6A_sec2a" else ord20[1]
if (mode == "all")        { RENDER_SIDS <- ord20;    test_mode <- FALSE
} else if (mode == "test"){ RENDER_SIDS <- TEST_SID; test_mode <- TRUE
} else                    { RENDER_SIDS <- a[1];     test_mode <- TRUE }   # explicit <sid>
RENDER_SIDS <- RENDER_SIDS[RENDER_SIDS %in% ord20]; stopifnot(length(RENDER_SIDS) >= 1)
OUT_PDF <- if (nzchar(Sys.getenv("FINAL_OUT"))) file.path(SSC_FIG, Sys.getenv("FINAL_OUT")) else
  if (test_mode) file.path(SSC_FIG, sprintf("ssc_clustering_segmentation_report_TEST_%s.pdf", RENDER_SIDS[1])) else
  file.path(SSC_FIG, "ssc_clustering_segmentation_report.pdf")

# section grid + y-orientation (match the gradient report: B[2,2]<0 -> range(ys))
grid_of  <- function(sid) { m <- which(as.character(pd$sample_id) == sid); x <- pd$x[m]; y <- pd$y[m]
  xs <- sort(unique(x)); ys <- sort(unique(y)); list(m = m, x = x, y = y, xs = xs, ys = ys, ix = match(x,xs), iy = match(y,ys)) }
fill_g   <- function(g, v) { mm <- matrix(NA_real_, length(g$xs), length(g$ys)); mm[cbind(g$ix, g$iy)] <- v; mm }
ylim_of  <- function(sid, ys) { xf <- file.path(REG_CACHE, sprintf("nd2final_%s.rds", sid))
  if (file.exists(xf)) { B <- readRDS(xf)$B_msi_nd2; if (B[2,2] < 0) range(ys) else rev(range(ys)) } else range(ys) }

# ============================================================================
# PANELS
# ============================================================================
tissueness_panel <- function(g, yl, main, cap) {
  tn <- fill_g(g, pd$tissueness_px[g$m])                      # already 0..1 (p99.5-normalised)
  par(mar = PMAR)
  image(g$xs, g$ys, pmin(tn, 1), col = pal, asp = 1, useRaster = TRUE, zlim = c(0,1), ylim = yl,
        axes = FALSE, xlab = "", ylab = "", main = main, cex.main = 0.95)
  frame_box(); colorbar_img(pal, 1, min(g$ys), max(g$ys))
  mt <- fill_g(g, as.integer(pd$is_tissue[g$m])); mt[is.na(mt)] <- 0   # on-tissue boundary
  contour(g$xs, g$ys, mt, levels = 0.5, add = TRUE, drawlabels = FALSE, col = "white", lty = 3, lwd = 0.8)
  scalebar_bottom(1/PX_UM, cap = cap)
}
cat_panel <- function(g, vals, k, yl, main, cap, leg = NULL) {
  m <- fill_g(g, vals); cp <- cat_palette(k)
  par(mar = PMAR)
  image(g$xs, g$ys, m, col = cp, breaks = seq(0.5, k+0.5, 1), asp = 1, useRaster = TRUE, ylim = yl,
        axes = FALSE, xlab = "", ylab = "", main = main, cex.main = 0.95)
  frame_box()
  if (!is.null(leg)) legend("topright", legend = leg, fill = cp, cex = 0.52, border = NA,
                            bg = adjustcolor("white", 0.6), box.col = NA, inset = c(0.012, 0.012))
  scalebar_bottom(1/PX_UM, cap = cap)
}
ssc_mask_panel <- function(g, yl, main, cap) {
  m <- fill_g(g, as.integer(pd$is_tissue[g$m]) + 1L)          # 1 off -> grey30, 2 tissue -> white
  par(mar = PMAR)
  image(g$xs, g$ys, m, col = c("grey30","white"), breaks = c(0.5,1.5,2.5), asp = 1, useRaster = TRUE, ylim = yl,
        axes = FALSE, xlab = "", ylab = "", main = main, cex.main = 0.95)
  frame_box(); scalebar_bottom(1/PX_UM, cap = cap)
}
instances_panel <- function(sid, yl_dir_rev, main) {
  fz <- cache_in(sprintf("zones_%s.rds", sid))
  par(mar = PMAR)
  if (!file.exists(fz)) { plot.new(); title(main, cex.main = 0.95); text(0.5,0.5,"(no zones)",col="grey50"); frame_box(); return(invisible(0)) }
  z <- readRDS(fz); xs <- sort(unique(z$x)); ys <- sort(unique(z$y)); ix <- match(z$x,xs); iy <- match(z$y,ys)
  fz2 <- function(v) { mm <- matrix(NA_real_, length(xs), length(ys)); mm[cbind(ix,iy)] <- v; mm }
  ids <- sort(unique(z$instance[z$instance > 0])); k <- length(ids); cp <- cat_palette(max(k,2))
  yl <- if (yl_dir_rev) rev(range(ys)) else range(ys)
  tis <- fz2(as.integer(z$is_tissue))                        # grey tissue base
  image(xs, ys, tis, col = c(NA, "grey80"), breaks = c(-0.5,0.5,1.5), asp = 1, useRaster = TRUE, ylim = yl,
        axes = FALSE, xlab = "", ylab = "", main = main, cex.main = 0.95)
  sdm <- fz2(z$signed_dist_um); instm <- fz2(as.numeric(z$instance))
  for (j in seq_along(ids)) { sd_k <- sdm; sd_k[is.na(instm) | instm != ids[j]] <- NA
    if (any(is.finite(sd_k))) contour(xs, ys, sd_k, levels = 0, add = TRUE, drawlabels = FALSE, col = cp[((j-1) %% length(cp))+1], lwd = 1.7) }
  sf <- z[z$is_surface, ]; if (nrow(sf)) points(sf$x, sf$y, pch = 15, cex = 0.18, col = adjustcolor("red", 0.55))
  for (j in seq_along(ids)) { cc <- z[z$instance == ids[j], ]; text(mean(cc$x), mean(cc$y), ids[j], col = "black", cex = 0.62, font = 2) }
  frame_box(); scalebar_bottom(1/PX_UM, cap = sprintf("organoid instances - connected comp >=%d px - N=%d", MIN_INSTANCE_PX, k))
  invisible(k)
}

# ============================================================================
# PAGE 1 - METHODS + COMPOSITION
# ============================================================================
methods_lines <- c(
  "SSC ON-TISSUE DELINEATION + ORGANOID SEGMENTATION  (Phase 04 + 08/09)",
  "",
  "1. Feature set: 348 curated on-tissue ions (peaks_curated.rds).",
  "2. Tissueness proxy: per-pixel curated-signal sum, clipped at its p99.5 (0..1).",
  "3. SSC pass 1 (per section): spatialShrunkenCentroids r=2, k=4,",
  "   s in {6,9,12}, weights=adaptive (s=9 default).",
  "4. On-tissue cut: keep k4 clusters with mean tissueness >= 50% of the",
  "   richest cluster,  UNION  floor80 (pixels >= section 80th-pct tissueness)",
  "   -> is_tissue. floor80 recovers epithelium edges without flooding lumen.",
  "5. SSC pass 2 (within tissue): k=10 local substructure; per-section",
  "   centroids pooled and harmonized by ward.D2 on log1p correlation,",
  "   K chosen by max mean silhouette -> chemotype (shared across sections).",
  "6. Organoid instances: connected components (4-conn) of is_tissue,",
  sprintf("   drop < %d px, refined (one connected ROI = one id) -> instances_final.", MIN_INSTANCE_PX),
  "   Surface = tissue pixel touching a different instance / non-tissue;",
  "   signed distance (RANN nn2, 10 um/px): - inward, + outward.",
  "",
  "Per-dataset pages: 1 tissueness+outline | 2 SSC k=4 | 3 on-tissue mask |",
  "4 organoid instances (coloured outline + red surface + id).",
  "All maps read directly from the cached SSC columns - no recompute."
)
front_page <- function() {
  layout(matrix(c(1,1,2,3), nrow = 2)); par(oma = c(0.5, 0.5, 3.4, 0.5))
  # left: methods text
  par(mar = c(1,1,1,1)); plot.new(); plot.window(c(0,1), c(0,1))
  text(0.0, 1.0, paste(methods_lines, collapse = "\n"), adj = c(0,1), cex = 0.82, family = "mono")
  # per-section tissue-pixel row indices (one factor pass, reused by both bars)
  idx <- split(seq_len(nrow(pd))[pd$is_tissue], pd$sample_id[pd$is_tissue])[ord20]
  # right-top: chemotype composition stacked bar (fraction of tissue px per chemotype)
  comp <- vapply(ord20, function(s) tabulate(pd$chemotype[idx[[s]]], nbins = Kc), numeric(Kc))
  prop <- sweep(comp, 2, pmax(colSums(comp),1), "/")
  par(mar = c(6.5, 4, 3, 1))
  barplot(prop, col = cat_palette(Kc), names.arg = disp_id(ord20), las = 2, cex.names = 0.5, border = NA,
          ylab = "chemotype fraction", main = "chemotype composition per section (within tissue)", cex.main = 0.9)
  legend("topright", legend = paste("chemo", seq_len(Kc)), fill = cat_palette(Kc), bty = "n", cex = 0.6, border = NA, horiz = TRUE, inset = c(0,-0.08), xpd = NA)
  # right-bottom: on-tissue pixel counts
  ntis <- vapply(ord20, function(s) length(idx[[s]]), numeric(1))
  par(mar = c(6.5, 4, 3, 1))
  barplot(ntis, names.arg = disp_id(ord20), las = 2, cex.names = 0.5, col = "grey45", border = NA,
          ylab = "on-tissue pixels", main = "on-tissue pixel count per section", cex.main = 0.9)
  mtext("SSC on-tissue clustering + organoid segmentation - methods & composition (20 sections, single condition)",
        outer = TRUE, line = 1.2, cex = 1.25, font = 2)
}

# ============================================================================
# PER-DATASET PAGE
# ============================================================================
dataset_page <- function(sid) {
  g <- grid_of(sid); yl <- ylim_of(sid, g$ys); rev_y <- identical(yl, rev(range(g$ys)))
  ntis <- sum(pd$is_tissue[g$m])
  layout(rbind(c(1,2), c(3,4))); par(oma = c(1.5, 1, 4.2, 1))
  tissueness_panel(g, yl, "1) tissueness (curated-ion signal)", "curated signal - viridis - p99.5 + on-tissue outline")
  cat_panel(g, pd$ssc_k4_sec[g$m], 4, yl, "2) SSC pass-1 k=4 clusters", "SSC r=2 s=9 adaptive - k=4",
            leg = paste("cluster", 1:4))
  ssc_mask_panel(g, yl, "3) on-tissue mask", "k4 >=50% richest  U  floor80 - white=tissue")
  k <- instances_panel(sid, rev_y, "4) organoid instance segmentation")
  mtext(sprintf("%s   -   %d organoid%s, %s on-tissue pixels", disp_id(sid), k, if (k==1) "" else "s", format(ntis, big.mark=",")),
        outer = TRUE, line = 2.2, cex = 1.5, font = 2)
  mtext("SSC on-tissue delineation -> organoid segmentation   (all maps from cached peaks_tissue_combined.rds / zones_<sid>.rds)",
        outer = TRUE, line = 0.6, cex = 0.82, col = "grey25")
}

# ============================================================================
# RENDER
# ============================================================================
pdf(OUT_PDF, width = 16, height = 10)
log_msg("front page (methods + composition)"); front_page()
for (sid in RENDER_SIDS) { log_msg("dataset page: %s", sid); dataset_page(sid) }
dev.off()
log_msg("DONE -> %s (front + %d dataset page%s)", OUT_PDF, length(RENDER_SIDS), if (length(RENDER_SIDS)==1) "" else "s")