#!/usr/bin/env Rscript
# ============================================================================
# 01_apical_annotate.R - per-section ANNOTATION PDF for scoring each
#   organoid's apical orientation (basolateral-out / apical-out / mixed) in Adobe.
#
#   One page per section: native-resolution brightfield crop (optical_<sid>.png,
#   from R/05_registration_refine/03_native_crops.R) with a clean coloured +
#   numbered outline per segmented organoid (instances_<sid>.rds, from
#   R/08_organoid_gradient_survey/01_segment_organoids.R), projected MSI->native
#   via B_msi_nd2.
#
#   Also writes a SIDECAR of every organoid's label position in PDF points so
#   R/10_apical_annotation/02_apical_parse.R can map each Adobe comment back to a
#   (sid, instance). Coordinate
#   convention matches D:/R/PeakMe/PeakMe_GCPL/phaseZC_pdf_to_pairs.R:
#     pdf_y_from_top = PAGE_H_PT - y_pdfcoord    (PAGE_H_PT = height_in * 72)
#
#   Fixed page size (14 x 10 in) on EVERY page so PAGE_H_PT is constant; each
#   crop is letterboxed via asp=1 and exact label coords are captured with
#   grconvertX/Y, so letterboxing needs no special handling.
#
# Usage:
#   Rscript R/10_apical_annotation/01_apical_annotate.R              # all 20 sections
#   Rscript R/10_apical_annotation/01_apical_annotate.R <sample_id>  # one section (TEST)
#
# In : cache/register/nd2final_<sid>.rds, figures/registration/crops/optical_<sid>.png,
#      cache/instances_<sid>.rds
# Out: figures/annotation/organoid_apical_annotation[_TEST_<sid>].pdf
#      cache/organoid_apical_label_positions[_TEST_<sid>].rds
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
suppressPackageStartupMessages({ library(Cardinal); library(png); library(EBImage) })

REG_CACHE <- file.path(CACHE_DIR, "register")
CROP_DIR  <- file.path(FIG_DIR, "registration", "crops")
# dedicated pre-step folder (run right after brightfield alignment, before the
# gradient analysis); kept separate from figures/gradient on purpose.
ANNOT_FIG <- file.path(FIG_DIR, "annotation")
SCALE_UM  <- 100; PX_UM <- MSI_PIXEL_UM
PAGE_W_IN <- 14; PAGE_H_IN <- 10
PAGE_W_PT <- PAGE_W_IN * 72; PAGE_H_PT <- PAGE_H_IN * 72
# per-organoid palette so IDs/colours cross-reference with the
# citrate_gradient_perdataset report
ORG_PAL <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628","#f781bf",
             "#1b9e77","#d95f02","#7570b3","#66a61e","#e7298a")
org_colors <- function(ids) setNames(ORG_PAL[((seq_along(ids)-1) %% length(ORG_PAL))+1], as.character(ids))
log_msg <- function(...) message(sprintf("[102] %s", sprintf(...)))

# ---- section ordering (by slide block) -------------------
mse <- readRDS(TISSUE_MSE)
pd  <- as.data.frame(pixelData(mse))
SIDS20  <- levels(pixelData(mse)$sample_id); if (is.null(SIDS20)) SIDS20 <- sort(unique(as.character(pd$sample_id)))
GROUP20 <- setNames(rep("incubated_5min", length(SIDS20)), SIDS20)
FC      <- setNames(rep("#444444", length(SIDS20)), SIDS20)
ord20   <- sort(SIDS20)
rm(mse); gc(verbose = FALSE)

a <- commandArgs(trailingOnly = TRUE)
RENDER_SIDS <- if (length(a) == 0) ord20 else a[a %in% SIDS20]
stopifnot(length(RENDER_SIDS) >= 1)
TEST <- length(RENDER_SIDS) == 1 && length(a) >= 1

OUT_PDF <- if (TEST) file.path(ANNOT_FIG, sprintf("organoid_apical_annotation_TEST_%s.pdf", RENDER_SIDS)) else
                     file.path(ANNOT_FIG, "organoid_apical_annotation.pdf")
OUT_SIDECAR <- if (TEST) file.path(CACHE_DIR, sprintf("organoid_apical_label_positions_TEST_%s.rds", RENDER_SIDS)) else
                         file.path(CACHE_DIR, "organoid_apical_label_positions.rds")

# ============================================================================
# helpers
# ============================================================================
# clean per-instance binary mask -> list of native-space closed outline polygons
instance_outlines <- function(M, B, cx0, cy0) {
  M <- EBImage::fillHull(EBImage::closing(M > 0, EBImage::makeBrush(3, "disc")))
  out <- list()
  for (co in EBImage::ocontour(EBImage::bwlabel(M))) {
    if (nrow(co) < 6) next
    p <- apply_affine(B, cbind(co[, 1] + 0.5, co[, 2] + 0.5))   # MSI grid -> native nd2
    out[[length(out) + 1L]] <- list(x = c(p[, 1] - cx0 + 1, p[1, 1] - cx0 + 1),
                                    y = c(p[, 2] - cy0 + 1, p[1, 2] - cy0 + 1))
  }
  out
}

# draw a high-contrast organoid ID label (white halo behind coloured number)
label_halo <- function(x, y, txt, col, cex = 1.5) {
  off <- 0.6 * (cex / 1.5) * strwidth("0", cex = cex) / 4 + 0.4
  for (dx in c(-1, 1)) for (dy in c(-1, 1))
    text(x + dx * off, y + dy * off, txt, col = "white", cex = cex, font = 2, xpd = NA)
  text(x, y, txt, col = col, cex = cex, font = 2, xpd = NA)
}

# choose a label anchor JUST OUTSIDE the organoid outline, in the direction with
# the most free space, so the number never covers identification-critical
# structures (lumen / apical face). Scans compass directions; for each, takes the
# farthest boundary point along it and offsets the number a small gap beyond,
# penalising off-image positions and proximity to OTHER organoids' outlines.
# Coords are crop-local user (native px). Returns list(lx, ly, bx, by).
place_label <- function(cx, cy, self_pts, other_pts, cw, ch, labW, labH) {
  gap <- 0.7 * labH; rad <- 1.3 * labH
  ang <- seq(0, 2 * pi, length.out = 17)[-17]
  best <- NULL; best_pen <- Inf
  for (i in seq_along(ang)) {
    d    <- c(cos(ang[i]), sin(ang[i]))
    proj <- (self_pts[, 1] - cx) * d[1] + (self_pts[, 2] - cy) * d[2]
    b    <- self_pts[which.max(proj), ]
    half <- abs(d[1]) * labW / 2 + abs(d[2]) * labH / 2
    L    <- c(b[1] + d[1] * (gap + half), b[2] + d[2] * (gap + half))
    pen  <- i * 1e-6 +
      60 * ((L[1] - labW/2 < 2) + (L[1] + labW/2 > cw - 1) +
            (L[2] - labH/2 < 2) + (L[2] + labH/2 > ch - 1))
    if (nrow(other_pts))
      pen <- pen + sum((other_pts[, 1] - L[1])^2 + (other_pts[, 2] - L[2])^2 < rad^2)
    if (pen < best_pen) { best_pen <- pen; best <- list(lx = L[1], ly = L[2], bx = b[1], by = b[2]) }
  }
  best
}

# ============================================================================
# cover page
# ============================================================================
cover_page <- function() {
  par(mar = c(2, 2, 2, 2)); plot.new()
  lines <- c(
    "Organoid apical-orientation annotation",
    "",
    "One page per section. Each segmented organoid (from the MSI data) is drawn",
    "as a coloured outline with a number placed JUST OUTSIDE it (a thin leader line",
    "points back to the organoid), so the number never hides the lumen / apical face.",
    "",
    "HOW TO ANNOTATE (Adobe Acrobat / Reader):",
    "  1. For each numbered organoid, add a text/sticky-note comment ON or NEXT TO it.",
    "  2. The comment text must be one of:   in    out    mixed",
    "       in    = basolateral-out        out  = apical-out        mixed = mixed",
    "     (synonyms accepted: 'apical in'/'ai', 'apical out'/'ao', 'mix'.)",
    "  3. For safety you may prefix the organoid number, e.g.  '3 out'  or  '3: mixed'.",
    "     A number-prefixed comment maps by ID; an un-prefixed one maps to the",
    "     nearest organoid on the page (place it close to the number).",
    "  4. Save the PDF, then tell Claude \"apical annotation done\" to run R/103.",
    "",
    "Outline colours/IDs match the citrate_gradient_perdataset report for cross-reference.",
    "Scale bar = 100 um."
  )
  text(0, 1, paste(lines, collapse = "\n"), adj = c(0, 1), cex = 1.05, family = "mono", xpd = NA)
}

# ============================================================================
# render
# ============================================================================
dir.create(ANNOT_FIG, showWarnings = FALSE, recursive = TRUE)
pdf(OUT_PDF, width = PAGE_W_IN, height = PAGE_H_IN)
sidecar <- list()

if (!TEST) cover_page()

for (sid in RENDER_SIDS) {
  xf <- file.path(REG_CACHE, sprintf("nd2final_%s.rds", sid))
  pf <- file.path(CROP_DIR, sprintf("optical_%s.png", sid))
  # prefer final (09/04) > cleaned (09/03) > split (09/01) > raw (08/01) segmentation
  inf_final <- file.path(CACHE_DIR, sprintf("instances_final_%s.rds", sid))
  inf_clean <- file.path(CACHE_DIR, sprintf("instances_clean_%s.rds", sid))
  inf_split <- file.path(CACHE_DIR, sprintf("instances_split_%s.rds", sid))
  inf <- if (file.exists(inf_final)) inf_final else if (file.exists(inf_clean)) inf_clean else
         if (file.exists(inf_split)) inf_split else file.path(CACHE_DIR, sprintf("instances_%s.rds", sid))
  grp <- GROUP20[sid]; fc <- FC[sid]
  if (!file.exists(xf) || !file.exists(pf) || !file.exists(inf)) {
    par(mar = c(2, 2, 3, 2)); plot.new()
    title(sprintf("%s  (%s)  -  missing inputs", sub("AO_", "", sid), grp), col.main = fc)
    text(0.5, 0.5, "(no native crop / registration / instances cache)", col = "grey50")
    log_msg("%s: SKIP (missing inputs)", sid); next
  }
  Xr <- readRDS(xf); B <- Xr$B_msi_nd2; smn <- Xr$scale_msi_nd2
  cx0 <- Xr$crop[1]; cy0 <- Xr$crop[2]
  om <- png::readPNG(pf); if (length(dim(om)) == 3) om <- om[, , 1]
  cw <- ncol(om); ch <- nrow(om)

  inst <- readRDS(inf)                                   # x, y, instance, is_surface
  inst <- inst[inst$instance > 0, ]
  W <- max(inst$x); H <- max(inst$y)
  lab <- matrix(0L, W, H); lab[cbind(inst$x, inst$y)] <- as.integer(inst$instance)
  ids <- sort(unique(as.integer(inst$instance)))
  org_col <- org_colors(ids)

  par(mar = c(1, 1, 3, 1))
  plot.new(); plot.window(c(1, cw), c(ch, 1), asp = 1)
  rasterImage(om / max(om), 1, ch, cw, 1, interpolate = TRUE)
  title(sprintf("%s   (%s)   -   %d organoid%s", sub("AO_", "", sid), grp,
                length(ids), if (length(ids) == 1) "" else "s"),
        cex.main = 1.2, col.main = fc, font.main = 2)

  # precompute each organoid's outline polygons + centroid (crop-local native px)
  polys_of <- list(); cen_of <- list(); pts_of <- list()
  for (k in ids) {
    polys <- instance_outlines(lab == k, B, cx0, cy0)
    cen   <- apply_affine(B, matrix(c(mean(inst$x[inst$instance == k]),
                                      mean(inst$y[inst$instance == k])), nrow = 1))
    pts   <- if (length(polys)) do.call(rbind, lapply(polys, function(p) cbind(p$x, p$y)))
             else matrix(c(cen[1] - cx0 + 1, cen[2] - cy0 + 1), 1, 2)
    polys_of[[as.character(k)]] <- polys; cen_of[[as.character(k)]] <- cen; pts_of[[as.character(k)]] <- pts
  }
  # outlines first, so leader lines/labels sit on top
  for (k in ids) { kc <- org_col[as.character(k)]
    for (poly in polys_of[[as.character(k)]]) polygon(poly$x, poly$y, border = kc, lwd = 1.8) }
  # numbers placed just OUTSIDE each organoid in the freest direction + leader
  labH <- strheight("0", cex = 1.6)
  for (k in ids) {
    kc  <- org_col[as.character(k)]; ck <- as.character(k)
    cen <- cen_of[[ck]]; cx <- cen[1] - cx0 + 1; cy <- cen[2] - cy0 + 1
    others <- if (length(ids) > 1) do.call(rbind, pts_of[setdiff(as.character(ids), ck)]) else matrix(numeric(0), 0, 2)
    pl <- place_label(cx, cy, pts_of[[ck]], others, cw, ch, strwidth(ck, cex = 1.6), labH)
    segments(pl$bx, pl$by, pl$lx, pl$ly, col = kc, lwd = 1.1, xpd = NA)   # leader to organoid
    label_halo(pl$lx, pl$ly, ck, kc, cex = 1.6)
    ndc_x <- grconvertX(pl$lx, "user", "ndc"); ndc_y <- grconvertY(pl$ly, "user", "ndc")
    sidecar[[length(sidecar) + 1L]] <- data.frame(
      sid = sid, group = grp, instance = k, original_text = ck,
      pdf_x = ndc_x * PAGE_W_PT, pdf_y_from_top = (1 - ndc_y) * PAGE_H_PT,
      cx_native = cen[1], cy_native = cen[2], stringsAsFactors = FALSE)
  }
  # scale bar (native px)
  umpx <- PX_UM / smn; bp <- SCALE_UM / umpx
  segments(cw - bp - 10, ch - 14, cw - 10, ch - 14, lwd = 3, col = "black")
  text(cw - bp / 2 - 10, ch - 28, sprintf("%d um", SCALE_UM), cex = 0.7)
  log_msg("%s: BF %dx%d, %d organoid(s)", sid, cw, ch, length(ids))
}
dev.off()

# ---- sidecar: page index = render order (+1 for cover unless TEST) ----------
sc <- if (length(sidecar)) do.call(rbind, sidecar) else
        data.frame(sid = character(), group = character(), instance = integer(),
                   original_text = character(), pdf_x = numeric(),
                   pdf_y_from_top = numeric(), cx_native = numeric(),
                   cy_native = numeric(), stringsAsFactors = FALSE)
page_of <- setNames(seq_along(RENDER_SIDS) + (if (TEST) 0L else 1L), RENDER_SIDS)
sc$page <- as.integer(page_of[sc$sid])
sc <- sc[, c("page", "sid", "group", "instance", "original_text",
             "pdf_x", "pdf_y_from_top", "cx_native", "cy_native")]
attr(sc, "page_w_pt") <- PAGE_W_PT; attr(sc, "page_h_pt") <- PAGE_H_PT
attr(sc, "render_sids") <- RENDER_SIDS; attr(sc, "has_cover") <- !TEST
saveRDS(sc, OUT_SIDECAR)

log_msg("DONE -> %s (%d section page%s, %d organoid labels)",
        OUT_PDF, length(RENDER_SIDS), if (length(RENDER_SIDS) == 1) "" else "s", nrow(sc))
log_msg("sidecar -> %s", OUT_SIDECAR)
