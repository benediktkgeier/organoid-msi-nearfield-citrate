#!/usr/bin/env Rscript
# ============================================================================
# 14_citrate_gradient_report_final_msigrid.R - VARIANT of the FINAL per-dataset
#   citrate-gradient report (10_citrate_gradient_report_final.R) that ALSO draws
#   the organoid segmentation outline in NATIVE MSI-GRID space, directly on the
#   MSI citrate ion panel (panel 1), via the native = FALSE branch of
#   class_outlines() (lib_nearfield_viz_whole.R:54-57).
#
#   Motivation: panels 5 & 6 of the final report draw the SAME zones_<sid>.rds
#   segmentation, but as instance_outlines_native() (binary-mask -> closing +
#   fillHull -> ocontour -> MSI->BF affine). That path flips (B[2,2]<0), scales
#   (~5.4x) and slightly rotates the outline, and rounds/holefills the shape, so
#   the boundary looks differently shaped than the raw MSI-grid contour used in
#   09_citrate_gradient_perdataset_v3.R. This variant adds the un-warped MSI-grid
#   outline on panel 1 so the true segmentation shape is visible in the SAME
#   report, directly comparable to the BF-warped outlines in panels 5 & 6.
#
#   ONLY DIFFERENCE vs 10_citrate_gradient_report_final.R:
#     - panel 1 (MSI citrate ion) now overlays class_outlines(native = FALSE):
#       per-instance signed-distance = 0 contour on the MSI grid, coloured by
#       consensus apical class (green basolateral-out / magenta apical-out /
#       grey mixed|unannotated) - same colour convention as panels 5 & 6.
#     - output file -> citrate_gradient_report_final_msigrid.pdf
#   Panels 2-6, the overview page, titles and legend are otherwise identical.
#
# In : cache/peaks_tissue_combined.rds, cache/zones_<sid>.rds,
#      cache/register/nd2final_<sid>.rds, figures/registration/crops/optical_<sid>.png,
#      figures/if_registration/crops/if_optical_<sid>.png,
#      results/annotation/apical_map_consensus.csv
# Out: figures/gradient/citrate_gradient_report_final_msigrid[_TEST_<sid>].pdf
# Usage: Rscript R/11_per_organoid_final/14_citrate_gradient_report_final_msigrid.R [all|test|<sid>]
#        all (default) = 20 dataset pages + overview;  test = front/overview + one
#        representative dataset (-> ..._TEST_<sid>.pdf);  <sid> = that one dataset.
# ============================================================================

ROOT <- file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final")
source(file.path(ROOT, "R/00_lib/gradient_config.R"))
source(file.path(ROOT, "R/00_lib/lib_register.R"))
source(file.path(ROOT, "R/00_lib/lib_nearfield_viz.R"))         # LOCKED 4-view toolkit (also sources lib_citrate)
source(file.path(ROOT, "R/00_lib/lib_nearfield_viz_whole.R"))   # LOCKED whole-region variants (class_outlines)
source(file.path(ROOT, "R/00_lib/if_config.R"))                 # BF_UMPX / HR_UMPX for the IF crop scale bar
source(file.path(ROOT, "R/00_lib/lib_report_frame.R"))          # PMAR / frame_box / scalebar_bottom / colorbar_img / SCALE_UM
source(file.path(ROOT, "R/00_lib/lib_gradient_report.R"))       # shared per-dataset panels (sec_data + 6 panels)
suppressPackageStartupMessages({ library(Cardinal); library(png); library(viridisLite); library(EBImage) })

log_msg <- function(...) message(sprintf("[final-msigrid] %s", sprintf(...)))

# ---- display constants -----------------------------------------------------
PX_UM      <- MSI_PIXEL_UM
MIN_ZONE_PX<- 3L
OV_MAX_UM  <- 100                        # overview mini-curves: cap outward x-axis at 100 um
EXT_IF     <- 0.30                       # IF crop = MSI bbox extended 30% (matches 12_dataset_pairs_hq.R)
IF_CROPS <- file.path(IF_FIG, "crops")
BF_CROPS <- file.path(FIG_DIR, "registration", "crops")
OUTLINE_LWD<- 1.6                        # MSI-grid class-outline line width on panel 1

# ============================================================================
# SHARED SETUP (mirrors 08 + 09 + 10)
# ============================================================================
nfviz_load_ion(CIT_MZ)                    # sets globals: pd (with gidx, is_tissue), val_cit, HI_CIT (citrate anchored)
HI <- HI_CIT                              # global p99.5 clip for the citrate ion image
stopifnot("is_tissue" %in% names(pd))

CONSENS_CSV <- file.path(RES_DIR, "annotation", "apical_map_consensus.csv")
stopifnot(file.exists(CONSENS_CSV))
ap     <- read.csv(CONSENS_CSV, stringsAsFactors = FALSE)
ap_map <- setNames(ap$apical_class, paste(ap$sid, ap$instance))

ord20 <- levels(pd$sample_id); if (is.null(ord20)) ord20 <- unique(as.character(pd$sample_id))
ord20 <- sort(ord20)

# Args:  (none)|all -> full 20-dataset report;  test -> quick TEST report (front/overview +
#        one representative dataset page);  <sid> -> that one dataset (also a TEST report).
a <- commandArgs(trailingOnly = TRUE)
mode <- if (length(a) == 0) "all" else tolower(a[1])
TEST_SID <- if ("AO_0h_sl6A_sec2a" %in% ord20) "AO_0h_sl6A_sec2a" else ord20[1]
if (mode == "all")       { RENDER_SIDS <- ord20;    test_mode <- FALSE
} else if (mode == "test"){ RENDER_SIDS <- TEST_SID; test_mode <- TRUE
} else                   { RENDER_SIDS <- a[1];     test_mode <- TRUE }   # explicit <sid>
RENDER_SIDS <- RENDER_SIDS[RENDER_SIDS %in% ord20]
stopifnot(length(RENDER_SIDS) >= 1)
OUT_PDF <- if (nzchar(Sys.getenv("FINAL_OUT"))) file.path(GRAD_FIG, Sys.getenv("FINAL_OUT")) else
  if (test_mode)
  file.path(GRAD_FIG, sprintf("citrate_gradient_report_final_msigrid_TEST_%s.pdf", RENDER_SIDS[1])) else
  file.path(GRAD_FIG, "citrate_gradient_report_final_msigrid.pdf")

# GLOBAL weather-heatmap clip (view 6): pooled p99.9 across ALL 20 (comparable page-to-page).
HEAT_HI <- heatmap_global_clip(ord20)
log_msg("HI(p99.5)=%.3g | HEAT_HI(global p99.9)=%.3g", HI, HEAT_HI)

# Per-dataset panels (sec_data + the 6 panel renderers) come from
# R/00_lib/lib_gradient_report.R (sourced above), shared with 10_..._final.R.

# ============================================================================
# OVERVIEW mini plot (A): per-organoid outward curves + class-mean trend lines
# ============================================================================
class_trend_mini <- function(sid) {
  par(mar = c(2.6, 2.3, 1.6, 0.6))
  fz <- cache_in(sprintf("zones_%s.rds", sid))
  if (!file.exists(fz)) { plot.new(); title(disp_id(sid), cex.main = 0.7); text(0.5,0.5,"(no zones)",col="grey50",cex=0.6); return(invisible()) }
  z <- readRDS(fz); ids <- sort(unique(z$instance[z$instance > 0]))
  zsel <- which(OUT_ZONE_UM <= OV_MAX_UM); nz <- length(zsel)   # cap the overview x-axis at OV_MAX_UM
  prof_k <- function(k) {
    surf <- z$gidx[z$instance == k & !is.na(z$zone_in) & z$zone_in == 1]
    rc <- mean(val_cit[surf], na.rm = TRUE); inside <- mean(val_cit[z$gidx[z$instance == k]], na.rm = TRUE)
    out <- vapply(zsel, function(zz){ g <- z$gidx[z$instance_catch == k & !is.na(z$zone_out) & z$zone_out == zz]
      if (length(g) >= MIN_ZONE_PX) mean(val_cit[g], na.rm = TRUE) else NA_real_ }, numeric(1))
    v <- c(inside, out); if (is.finite(rc) && rc > 0) v/rc else v*NA
  }
  if (!length(ids)) { plot.new(); title(disp_id(sid), cex.main = 0.7); return(invisible()) }
  P <- sapply(ids, prof_k); if (is.null(dim(P))) P <- matrix(P, ncol = length(ids))
  cls <- vapply(ids, function(k){ c <- unname(ap_map[paste(sid,k)]); if (is.na(c)) "unannotated" else c }, character(1))
  xp <- 0:nz; ymax <- max(1.2, quantile(P, 0.97, na.rm = TRUE), na.rm = TRUE); if (!is.finite(ymax)) ymax <- 2
  plot(NA, xlim = c(0,nz), ylim = c(0,ymax), xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = disp_id(sid), cex.main = 0.7)
  axis(1, at = xp, labels = c("in", OUT_ZONE_UM[zsel]), cex.axis = 0.42, tcl = -0.2, mgp = c(0,0.02,0), las = 2)
  axis(2, at = c(0,1), labels = c("0","1"), cex.axis = 0.5, las = 1, tcl = -0.2, mgp = c(0,0.3,0))
  abline(h = 1, lty = 3, col = "grey75"); abline(v = 0.5, lty = 3, col = "grey70")
  for (j in seq_along(ids)) lines(xp, P[,j], col = adjustcolor(.class_col(sid, ids[j], ap_map), 0.30), lwd = 0.7)
  for (cc in c("basolateral_out","apical_out")) {
    selc <- which(cls == cc); if (!length(selc)) next
    mp <- rowMeans(P[, selc, drop = FALSE], na.rm = TRUE)
    lines(xp, mp, col = CCOL_WHOLE[[cc]], lwd = 2.4, type = "b", pch = 19, cex = 0.5)
  }
}

# ============================================================================
# PAGES
# ============================================================================
dataset_page <- function(sid) {
  sd  <- sec_data(sid); sec <- prep_sec(sid); z <- sec$z   # prep_sec already loaded zones_<sid>.rds
  inst_ids <- if (is.null(z)) integer(0) else sort(unique(z$instance[z$instance > 0]))
  n_annot  <- sum(!is.na(ap_map[paste(sid, inst_ids)]))
  layout(rbind(c(1,2,3), c(4,5,6))); par(oma = c(2.0, 1, 5, 1))
  ion_panel(sd$cit, sd$xs, sd$ys, HI, "1) MSI citrate ion + segmentation outline", "#444444", TRUE,
            "citrate m/z 191.02 - viridis - p99.5 - MSI-grid class outline", ylim = ylim_grid(sec))
  # MSI-GRID outline (native = FALSE): per-instance signed-dist=0 contour on the MSI
  # grid, class-coloured; drawn onto panel 1 (still the active region after ion_panel).
  class_outlines(sec, sid, ap_map, native = FALSE, lwd = OUTLINE_LWD)
  native_bf_panel(sid, "2) native brightfield", "brightfield - .nd2 native crop")
  ssc_mask_panel(sid, sec, "3) SSC on-tissue mask", "on-tissue - SSC+floor80 - white=tissue")
  if_panel(sid, "4) matched IF", "IF - ZO-1 red / DAPI cyan - 30x +30%")
  native_overlay_norings(sec, sid, ap_map, "5) citrate on brightfield (overlay)",
                         "citrate on BF - a0.60 - class outlines")
  weather_panel(sec, sid, ap_map, "6) citrate gradient map (weather, GLOBAL)",
                "weather - gblur 12um - global p99.9")
  mtext(sprintf("%s   -   %d organoid%s (%d annotated)", disp_id(sid), length(inst_ids),
                if (length(inst_ids)==1) "" else "s", n_annot), outer = TRUE, line = 3.1, cex = 1.6, font = 2)
  mtext("1 MSI citrate + MSI-grid outline | 2 native BF | 3 SSC on-tissue | 4 matched IF (ZO-1/DAPI) | 5 MSI<->BF overlay (no rings) | 6 weather gradient map (global)",
        outer = TRUE, line = 1.3, cex = 0.8, col = "grey25")
  legend_cols <- c(CCOL_WHOLE[["basolateral_out"]], CCOL_WHOLE[["apical_out"]], CCOL_WHOLE[["mixed"]], UNANNOT_COL)
  par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE); plot.new()
  legend("bottom", horiz = TRUE, bty = "n", cex = 0.8, lwd = 3, inset = c(0, 0.003),
         legend = c("basolateral-out","apical-out","mixed","unannotated"), col = legend_cols, text.col = "grey15")
}

overview_page <- function() {
  layout(matrix(1:40, nrow = 5, byrow = TRUE)); par(oma = c(0.5, 1, 4.5, 1))
  for (sid in ord20) {
    sec <- prep_sec(sid)
    draw_gradmap_whole(sec, sid, ap_map, disp_id(sid))
    class_trend_mini(sid)
  }
  mtext("OVERVIEW - all 20 datasets: citrate ion image + 50/100 um rings, paired with per-organoid outward curves",
        outer = TRUE, line = 2.7, cex = 1.3, font = 2)
  mtext("curves normalized to organoid surface (=1.0); faint = organoids (coloured by apical class); BOLD = class mean trend; MAGENTA=apical-out, GREEN=basolateral-out; x = 'in' + outward zones (um)",
        outer = TRUE, line = 1.0, cex = 0.82, col = "grey25")
}

# ============================================================================
# RENDER
# ============================================================================
dir.create(GRAD_FIG, showWarnings = FALSE, recursive = TRUE)
pdf(OUT_PDF, width = 16, height = 10)
for (sid in RENDER_SIDS) { log_msg("dataset page: %s", sid); dataset_page(sid) }
log_msg("overview page (all 20)"); overview_page()
dev.off()
log_msg("DONE -> %s (%d dataset page%s + overview)", OUT_PDF, length(RENDER_SIDS), if (length(RENDER_SIDS)==1) "" else "s")
