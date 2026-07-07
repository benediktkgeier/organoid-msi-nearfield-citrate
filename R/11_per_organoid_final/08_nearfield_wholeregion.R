#!/usr/bin/env Rscript
# ============================================================================
# 08_nearfield_wholeregion.R - WHOLE-REGION near-field citrate emission figure.
#   *** LOCKED *** - frozen visual spec; see docs/nearfield_wholeregion.md.
#   Companion to the per-organoid 04_nearfield_figure.R: same four views, but each
#   is projected across the ENTIRE measurement region of a section instead of being
#   cropped per organoid. EVERY organoid outline is coloured by its consensus apical
#   class (apical_map_consensus.csv): basolateral-out=green / apical-out=magenta /
#   mixed=grey / unannotated=grey70 (locked APICAL_COLS). Thin DRIVER over the
#   LOCKED toolkit (lib_nearfield_viz.R) +
#   whole-region variants (lib_nearfield_viz_whole.R, also LOCKED).
#
#   Four views (2x2 page, identical order to 04):
#     1) citrate overlay on native brightfield (viridis, constant alpha 0.60)
#     2) citrate ion image on MSI grid + 50/100 um signed-distance rings
#     3) interpolated weather-rainbow heatmap on brightfield
#     4) outward emission vectors (length=absolute, thickness=GLOBAL relative)
#   ORIENTATION: MSI-grid panels (2,4) are y-oriented to MATCH the native-BF panels
#   (1,3) via ylim_grid() in the lib - the MSI->BF affine is vertically flipped
#   (B[2,2] < 0). DO NOT revert to rev(range(ys)) (would flip views 2 & 4).
#
#   DEFAULT (no args / "all"): render ALL 20 datasets, one whole-region page each,
#   into a single multi-page PDF -> apical_nearfield_emission_figure_all_globalheat.pdf.
#   Arrow THICKNESS uses ONE GLOBAL p10..p90 emission range across all sections so
#   thickness is comparable page-to-page (matches the already-global HI_CIT clip).
#   A single <sid> arg renders just that section -> apical_nearfield_wholeregion_<sid>.pdf.
#
#   HEATMAP scaling (view 3) - *** LOCKED, UNIFORM ***:
#     ALWAYS ONE pooled GLOBAL p99.9 clip over all 20 datasets (dimmer than the old p99;
#     hotspots stay warm, surround cools), so the weather heatmap is comparable across
#     pages (+ colorbar). The old per-section p99 auto-scale mode has been REMOVED; the
#     'globalheat' arg is still accepted but redundant.
#     Out -> apical_nearfield_emission_figure_all_globalheat.pdf.
#
# In : cache/peaks_tissue_combined.rds, cache/zones_<sid>.rds,
#      cache/register/nd2final_<sid>.rds, figures/registration/crops/optical_<sid>.png,
#      results/annotation/apical_map_consensus.csv
# Out: figures/annotation/apical_nearfield_emission_figure_all_globalheat.pdf  (all 20)
#      figures/annotation/apical_nearfield_wholeregion_<sid>_globalheat.pdf     (single sid)
# Usage: Rscript R/11_per_organoid_final/08_nearfield_wholeregion.R [all|<sid>] [globalheat]
# ============================================================================

ROOT <- file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final")
source(file.path(ROOT, "R/00_lib/gradient_config.R"))
source(file.path(ROOT, "R/00_lib/lib_register.R"))
source(file.path(ROOT, "R/00_lib/lib_nearfield_viz.R"))         # LOCKED 4-view toolkit
source(file.path(ROOT, "R/00_lib/lib_nearfield_viz_whole.R"))   # whole-region variants

ANNOT_FIG <- file.path(FIG_DIR, "annotation")
dir.create(ANNOT_FIG, showWarnings = FALSE, recursive = TRUE)
CONSENS_CSV <- file.path(RES_DIR, "annotation", "apical_map_consensus.csv")
stopifnot(file.exists(CONSENS_CSV))

nfviz_load_ion(CIT_MZ)                          # sets val_cit, HI_CIT, pd (citrate [M-H]-)

# consensus class per (sid, instance)
ap <- read.csv(CONSENS_CSV, stringsAsFactors = FALSE)
ap_map <- setNames(ap$apical_class, paste(ap$sid, ap$instance))

# ---- dataset list (all sections, alphabetical) -----------------------------
SIDS20  <- levels(pd$sample_id); if (is.null(SIDS20)) SIDS20 <- sort(unique(as.character(pd$sample_id)))
ord20   <- sort(SIDS20)

a <- commandArgs(trailingOnly = TRUE)
GLOBAL_HEAT <- TRUE   # LOCKED: view-3 heatmap is ALWAYS pooled global p99.9 (old per-section auto-scale removed); the 'globalheat' arg is still accepted but redundant
a_sid <- a[!(tolower(a) %in% c("globalheat", "all"))]
if (length(a_sid) == 0) {
  RENDER_SIDS <- ord20[ord20 %in% SIDS20]
  OUT_PDF <- file.path(ANNOT_FIG, sprintf("apical_nearfield_emission_figure_all%s.pdf",
                                          if (GLOBAL_HEAT) "_globalheat" else ""))
} else {
  RENDER_SIDS <- a_sid[1][a_sid[1] %in% SIDS20]; stopifnot(length(RENDER_SIDS) == 1)
  OUT_PDF <- file.path(ANNOT_FIG, sprintf("apical_nearfield_wholeregion_%s%s.pdf",
                                          RENDER_SIDS, if (GLOBAL_HEAT) "_globalheat" else ""))
}

# ---- GLOBAL arrow-thickness range: p10..p90 over ALL organoids in RENDER_SIDS
allvex <- do.call(rbind, lapply(RENDER_SIDS, function(s) {
  sec <- prep_sec(s); ids <- sort(unique(sec$z$instance[sec$z$instance > 0]))
  if (length(ids)) data.frame(sid = s, instance = ids, stringsAsFactors = FALSE) else NULL
}))
nfviz_arrow_range(allvex)                       # sets global GLO, GHI (used by all pages)

# ---- LOCKED GLOBAL heatmap clip (view 3): one pooled p99.9 over all datasets (always)
HEAT_HI <- heatmap_global_clip(RENDER_SIDS)   # q = 0.999 (locked default in heatmap_global_clip)

# ---- one whole-region 2x2 page for a section -------------------------------
render_page <- function(sid) {
  sec <- prep_sec(sid)
  ids <- sort(unique(sec$z$instance[sec$z$instance > 0]))
  n_annot <- sum(!is.na(ap_map[paste(sid, ids)]))
  cat(sprintf("[08] %s: %d organoids (%d consensus-annotated)\n", sid, length(ids), n_annot))
  layout(matrix(1:4, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 4.6, 0))
  heat3 <- "3) citrate heatmap on BF (weather rainbow, GLOBAL scale)"
  draw_overlay_whole(sec, sid, ap_map, "1) citrate on brightfield (overlay)")
  draw_gradmap_whole(sec, sid, ap_map, "2) citrate ion image + 50/100 um rings")
  draw_heatmap_whole(sec, sid, ap_map, heat3, hi = HEAT_HI)
  draw_vectors_whole(sec, sid, ap_map, "4) emission vectors (length=absolute, thickness=relative)")
  mtext(sprintf("%s - whole-region near-field citrate emission   (%d organoids, %d annotated)",
                disp_id(sid), length(ids), n_annot),
        outer = TRUE, font = 2, cex = 1.12, line = 2.6)
  heat_note <- sprintf("view 3 heatmap on a GLOBAL scale (0..%.2g, p99.9 across all 20 datasets) - comparable page-to-page", HEAT_HI)
  mtext(sprintf("projected across the entire measurement region (not per organoid); outlines coloured by consensus apical class; arrows on a GLOBAL thickness scale; %s", heat_note),
        outer = TRUE, cex = 0.74, line = 1.4, col = "grey25")
  legend_cols <- c(CCOL_WHOLE[["basolateral_out"]], CCOL_WHOLE[["apical_out"]], CCOL_WHOLE[["mixed"]], UNANNOT_COL)
  par(fig = c(0, 1, 0, 1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE); plot.new()
  legend("bottom", horiz = TRUE, bty = "n", cex = 0.8, lwd = 3, inset = c(0, 0.003),
         legend = c("basolateral-out", "apical-out", "mixed", "unannotated"), col = legend_cols,
         text.col = "grey15")
}

# ===========================================================================
# RENDER multi-page PDF (one whole-region page per dataset)
# ===========================================================================
pdf(OUT_PDF, width = 11, height = 8.5)
for (sid in RENDER_SIDS) render_page(sid)
dev.off()
cat(sprintf("[08] DONE -> %s (%d page%s)\n", OUT_PDF, length(RENDER_SIDS),
            if (length(RENDER_SIDS) == 1) "" else "s"))
