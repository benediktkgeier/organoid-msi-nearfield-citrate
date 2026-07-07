#!/usr/bin/env Rscript
# ============================================================================
# 08b_view3_heatmap_hires.R - HIGH-RES export of VIEW 3 ONLY, one PNG per dataset.
#   Companion to 08_nearfield_wholeregion.R. That driver renders a 20-page 2x2 PDF
#   (views 1-4) per section; here we extract ONLY view 3 - the interpolated
#   weather-rainbow citrate heatmap on brightfield - for ALL 20 datasets, each as a
#   standalone high-resolution PNG (600 dpi) in its own folder.
#
#   Thin DRIVER over the LOCKED toolkit: uses draw_heatmap_whole() unchanged, on the
#   SAME global p99 heatmap clip (heatmap_global_clip over all 20 datasets) as the
#   apical_nearfield_emission_figure_all_globalheat.pdf report - so every PNG is on
#   one comparable colour scale (+ colorbar), page-to-page identical to that report.
#
# In : same inputs as 08 (cache/peaks_tissue_combined.rds, cache/register/*,
#      figures/registration/crops/*, results/annotation/apical_map_consensus.csv)
# Out: figures/annotation/nearfield_view3_heatmap_hires/view3_heatmap_<sid>.png  (x20)
# Usage: Rscript R/11_per_organoid_final/08b_view3_heatmap_hires.R [all|<sid>]
# ============================================================================

ROOT <- file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final")
source(file.path(ROOT, "R/00_lib/gradient_config.R"))
source(file.path(ROOT, "R/00_lib/lib_register.R"))
source(file.path(ROOT, "R/00_lib/lib_nearfield_viz.R"))         # LOCKED 4-view toolkit
source(file.path(ROOT, "R/00_lib/lib_nearfield_viz_whole.R"))   # whole-region variants

OUT_DIR <- file.path(FIG_DIR, "annotation", "nearfield_view3_heatmap_hires")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
CONSENS_CSV <- file.path(RES_DIR, "annotation", "apical_map_consensus.csv")
stopifnot(file.exists(CONSENS_CSV))

DPI <- 600L                                     # high-resolution raster export

nfviz_load_ion(CIT_MZ)                          # sets val_cit, HI_CIT, pd (citrate [M-H]-)

ap <- read.csv(CONSENS_CSV, stringsAsFactors = FALSE)
ap_map <- setNames(ap$apical_class, paste(ap$sid, ap$instance))

# ---- dataset list (all sections, alphabetical) -----------------------------
SIDS20 <- levels(pd$sample_id); if (is.null(SIDS20)) SIDS20 <- sort(unique(as.character(pd$sample_id)))
ord20  <- sort(SIDS20)

a <- commandArgs(trailingOnly = TRUE)
a_sid <- a[!(tolower(a) %in% c("all"))]
# optional numeric arg = HEAT_CLIP_Q override (calibration); optional "suffix=<str>" tags the filename
Q_OVERRIDE <- suppressWarnings(as.numeric(a_sid)); Q_OVERRIDE <- Q_OVERRIDE[!is.na(Q_OVERRIDE)][1]
FNAME_SUFFIX <- sub("^suffix=", "", a_sid[grepl("^suffix=", a_sid)])[1]; if (is.na(FNAME_SUFFIX)) FNAME_SUFFIX <- ""
a_sid <- a_sid[is.na(suppressWarnings(as.numeric(a_sid))) & !grepl("^suffix=", a_sid)]
RENDER_SIDS <- if (length(a_sid) == 0) ord20[ord20 %in% SIDS20] else {
  s <- a_sid[1][a_sid[1] %in% SIDS20]; stopifnot(length(s) == 1); s
}

# ---- GLOBAL heatmap clip (view 3): one clip pooled over ALL 20 datasets ------
# ALWAYS global (this export mirrors the *_globalheat report), even for a single sid.
# HEAT_CLIP_Q controls display intensity: the locked report uses p99 (0.99), which
# saturates the top 1% of pixels red (over-bright). To DIM the map (per user), raise
# the percentile toward the pooled global MAX. q = 1.0 eliminates the percentile
# entirely (scale to global max); ~0.999 gives a slighter reduction. Data unchanged.
HEAT_CLIP_Q <- if (!is.null(Q_OVERRIDE) && !is.na(Q_OVERRIDE)) Q_OVERRIDE else 0.999
HEAT_HI  <- heatmap_global_clip(ord20[ord20 %in% SIDS20], q = HEAT_CLIP_Q)
HEAT_P99 <- heatmap_global_clip(ord20[ord20 %in% SIDS20], q = 0.99)   # for reference/log only
cat(sprintf("[08b] heatmap clip: q=%.4g -> HI=%.3g   (locked-report p99=%.3g; ratio %.2fx dimmer)\n",
            HEAT_CLIP_Q, HEAT_HI, HEAT_P99, HEAT_HI / HEAT_P99))

# ---- one high-res view-3 PNG for a section ---------------------------------
render_view3 <- function(sid) {
  sec <- prep_sec(sid)
  ids <- sort(unique(sec$z$instance[sec$z$instance > 0]))
  n_annot <- sum(!is.na(ap_map[paste(sid, ids)]))
  cat(sprintf("[08b] %s: %d organoids (%d consensus-annotated)\n", sid, length(ids), n_annot))

  # size the canvas to the native-crop aspect (asp=1 in the panel => no distortion;
  # this just minimises letterboxing). Long dimension ~9in + margins for title/bar.
  asp_img <- if (!is.null(sec$cw) && !is.null(sec$ch)) sec$cw / sec$ch else 1
  if (asp_img >= 1) { w_in <- 9.0 + 1.4; h_in <- 9.0 / asp_img + 1.6 }
  else              { h_in <- 9.0 + 1.6; w_in <- 9.0 * asp_img + 1.4 }

  out_png <- file.path(OUT_DIR, sprintf("view3_heatmap_%s.png", sid))
  png(out_png, width = w_in, height = h_in, units = "in", res = DPI, type = "cairo")
  par(oma = c(2.2, 0.3, 3.0, 0.3))
  draw_heatmap_whole(sec, sid, ap_map,
                     "3) citrate heatmap on BF (weather rainbow, GLOBAL scale)", hi = HEAT_HI)
  mtext(sprintf("%s - whole-region citrate emission heatmap (view 3)   (%d organoids, %d annotated)",
                disp_id(sid), length(ids), n_annot),
        outer = TRUE, font = 2, cex = 1.15, line = 1.1)
  mtext(sprintf("weather-rainbow interpolated citrate [M-H]- on native brightfield; GLOBAL heatmap scale 0..%.2g (%s across all 20 datasets); 50/100 um rings; outlines by consensus apical class",
                HEAT_HI, if (HEAT_CLIP_Q >= 1) "global max" else sprintf("p%.4g", HEAT_CLIP_Q*100)),
        outer = TRUE, cex = 0.72, line = 0.0, col = "grey25")
  legend_cols <- c(CCOL_WHOLE[["basolateral_out"]], CCOL_WHOLE[["apical_out"]], CCOL_WHOLE[["mixed"]], UNANNOT_COL)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot.new()
  legend("bottom", horiz = TRUE, bty = "n", cex = 0.85, lwd = 3, inset = c(0, 0.004),
         legend = c("basolateral-out", "apical-out", "mixed", "unannotated"), col = legend_cols,
         text.col = "grey15")
  dev.off()
  cat(sprintf("[08b]   -> %s (%.1f x %.1f in @ %d dpi)\n", out_png, w_in, h_in, DPI))
}

for (sid in RENDER_SIDS) render_view3(sid)
cat(sprintf("[08b] DONE -> %s (%d PNG%s @ %d dpi, GLOBAL heat scale %.2f)\n",
            OUT_DIR, length(RENDER_SIDS), if (length(RENDER_SIDS) == 1) "" else "s", DPI, HEAT_HI))