#!/usr/bin/env Rscript
# ============================================================================
# 02_island_cleanup_canvas.R - Phase 09 canvas: per-section CLEAN numbered PDF
#   for the island CLEANUP step (mark segmented islands to DELETE or MERGE).
# ----------------------------------------------------------------------------
# Identical render to R/10_apical_annotation/01_apical_annotate.R (native BF +
# coloured numbered organoid outlines, the
# split segmentation when present), but with NO split/cut lines and a cleanup
# cover. The user marks islands in Adobe; a sibling parser (E2, built later) will
# read the comments via the sidecar this script writes.
#
# Proposed marking convention (refine with Claude before parsing):
#   DELETE an island : comment  del  (or x / delete / remove) on it.
#   MERGE islands     : comment  merge <tag>  on EACH island to combine; same tag
#                       = one merged organoid, e.g.  merge A  on both #4 and #7.
#   (Prefix the organoid number for safety, e.g.  4 del  or  7 merge A.)
#
# In : cache/instances_split_<sid>.rds (else instances_<sid>.rds; from
#      R/08_organoid_gradient_survey/01_segment_organoids.R and
#      R/09_organoid_refinement/01_organoid_split_apply.R),
#      cache/register/nd2final_<sid>.rds, figures/registration/crops/optical_<sid>.png
# Out: figures/annotation/organoid_island_cleanup.pdf
#      cache/organoid_island_label_positions.rds  (sidecar for the E2 parser)
# Usage: Rscript R/09_organoid_refinement/02_island_cleanup_canvas.R [sample_id]   (1 arg = TEST page)
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
suppressPackageStartupMessages({ library(Cardinal); library(png); library(EBImage) })

REG_CACHE <- file.path(CACHE_DIR, "register")
CROP_DIR  <- file.path(FIG_DIR, "registration", "crops")
ANNOT_FIG <- file.path(FIG_DIR, "annotation")
SCALE_UM  <- 100; PX_UM <- MSI_PIXEL_UM
PAGE_W_IN <- 14; PAGE_H_IN <- 10
PAGE_W_PT <- PAGE_W_IN * 72; PAGE_H_PT <- PAGE_H_IN * 72
ORG_PAL <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628","#f781bf",
             "#1b9e77","#d95f02","#7570b3","#66a61e","#e7298a")
org_colors <- function(ids) setNames(ORG_PAL[((seq_along(ids)-1) %% length(ORG_PAL))+1], as.character(ids))
log_msg <- function(...) message(sprintf("[42] %s", sprintf(...)))

# ---- shared render helpers (kept identical to R/10_apical_annotation/01_apical_annotate.R so the sidecar matches) --
instance_outlines <- function(M, B, cx0, cy0) {
  M <- EBImage::fillHull(EBImage::closing(M > 0, EBImage::makeBrush(3, "disc")))
  out <- list()
  for (co in EBImage::ocontour(EBImage::bwlabel(M))) {
    if (nrow(co) < 6) next
    p <- apply_affine(B, cbind(co[, 1] + 0.5, co[, 2] + 0.5))
    out[[length(out) + 1L]] <- list(x = c(p[, 1] - cx0 + 1, p[1, 1] - cx0 + 1),
                                    y = c(p[, 2] - cy0 + 1, p[1, 2] - cy0 + 1))
  }
  out
}
label_halo <- function(x, y, txt, col, cex = 1.5) {
  off <- 0.6 * (cex / 1.5) * strwidth("0", cex = cex) / 4 + 0.4
  for (dx in c(-1, 1)) for (dy in c(-1, 1))
    text(x + dx * off, y + dy * off, txt, col = "white", cex = cex, font = 2, xpd = NA)
  text(x, y, txt, col = col, cex = cex, font = 2, xpd = NA)
}
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

# ---- section ordering (sorted sids), matching R/10_apical_annotation/01_apical_annotate.R -------------------
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

OUT_PDF <- if (TEST) file.path(ANNOT_FIG, sprintf("organoid_island_cleanup_TEST_%s.pdf", RENDER_SIDS)) else
                     file.path(ANNOT_FIG, "organoid_island_cleanup.pdf")
OUT_SIDECAR <- if (TEST) file.path(CACHE_DIR, sprintf("organoid_island_label_positions_TEST_%s.rds", RENDER_SIDS)) else
                         file.path(CACHE_DIR, "organoid_island_label_positions.rds")

cover_page <- function() {
  par(mar = c(2, 2, 2, 2)); plot.new()
  lines <- c(
    "Organoid island cleanup  (delete / merge)",
    "",
    "One page per section. Each segmented organoid (after the morphological split)",
    "is a coloured outline with a number placed just outside it (leader line points",
    "to it). Use this to flag spurious islands to DELETE and pieces to MERGE.",
    "",
    "PROPOSED MARKING (Adobe Acrobat / Reader) - confirm with Claude before parsing:",
    "  DELETE an island : comment  del   (also: x / delete / remove)  on/next to it.",
    "  MERGE islands    : comment  merge <tag>  on EACH island to combine; the SAME",
    "                     tag means one merged organoid, e.g.  merge A  on #4 and #7.",
    "  For safety prefix the organoid number, e.g.  4 del   or   7 merge A.",
    "",
    "Numbers/colours match the split segmentation and the apical-annotation report.",
    "Save the PDF, then tell Claude \"island cleanup done\".",
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
  inf_final <- file.path(CACHE_DIR, sprintf("instances_final_%s.rds", sid))
  inf_clean <- file.path(CACHE_DIR, sprintf("instances_clean_%s.rds", sid))
  inf_split <- file.path(CACHE_DIR, sprintf("instances_split_%s.rds", sid))
  inf <- if (file.exists(inf_final)) inf_final else if (file.exists(inf_clean)) inf_clean else
         if (file.exists(inf_split)) inf_split else file.path(CACHE_DIR, sprintf("instances_%s.rds", sid))
  grp <- GROUP20[sid]; fc <- FC[sid]
  if (!file.exists(xf) || !file.exists(pf) || !file.exists(inf)) {
    par(mar = c(2, 2, 3, 2)); plot.new()
    title(sprintf("%s  (%s)  -  missing inputs", sub("AO_", "", sid), grp), col.main = fc)
    log_msg("%s: SKIP (missing inputs)", sid); next
  }
  Xr <- readRDS(xf); B <- Xr$B_msi_nd2; smn <- Xr$scale_msi_nd2
  cx0 <- Xr$crop[1]; cy0 <- Xr$crop[2]
  om <- png::readPNG(pf); if (length(dim(om)) == 3) om <- om[, , 1]
  cw <- ncol(om); ch <- nrow(om)

  inst <- readRDS(inf); inst <- inst[inst$instance > 0, ]
  W <- max(inst$x); H <- max(inst$y)
  lab <- matrix(0L, W, H); lab[cbind(inst$x, inst$y)] <- as.integer(inst$instance)
  ids <- sort(unique(as.integer(inst$instance))); org_col <- org_colors(ids)

  par(mar = c(1, 1, 3, 1))
  plot.new(); plot.window(c(1, cw), c(ch, 1), asp = 1)
  rasterImage(om / max(om), 1, ch, cw, 1, interpolate = TRUE)
  title(sprintf("%s   (%s)   -   %d organoid%s", sub("AO_", "", sid), grp,
                length(ids), if (length(ids) == 1) "" else "s"),
        cex.main = 1.2, col.main = fc, font.main = 2)

  polys_of <- list(); cen_of <- list(); pts_of <- list()
  for (k in ids) {
    polys <- instance_outlines(lab == k, B, cx0, cy0)
    cen   <- apply_affine(B, matrix(c(mean(inst$x[inst$instance == k]),
                                      mean(inst$y[inst$instance == k])), nrow = 1))
    pts   <- if (length(polys)) do.call(rbind, lapply(polys, function(p) cbind(p$x, p$y)))
             else matrix(c(cen[1] - cx0 + 1, cen[2] - cy0 + 1), 1, 2)
    polys_of[[as.character(k)]] <- polys; cen_of[[as.character(k)]] <- cen; pts_of[[as.character(k)]] <- pts
  }
  for (k in ids) { kc <- org_col[as.character(k)]
    for (poly in polys_of[[as.character(k)]]) polygon(poly$x, poly$y, border = kc, lwd = 1.8) }
  labH <- strheight("0", cex = 1.6)
  for (k in ids) {
    kc  <- org_col[as.character(k)]; ck <- as.character(k)
    cen <- cen_of[[ck]]; cx <- cen[1] - cx0 + 1; cy <- cen[2] - cy0 + 1
    others <- if (length(ids) > 1) do.call(rbind, pts_of[setdiff(as.character(ids), ck)]) else matrix(numeric(0), 0, 2)
    pl <- place_label(cx, cy, pts_of[[ck]], others, cw, ch, strwidth(ck, cex = 1.6), labH)
    segments(pl$bx, pl$by, pl$lx, pl$ly, col = kc, lwd = 1.1, xpd = NA)
    label_halo(pl$lx, pl$ly, ck, kc, cex = 1.6)
    ndc_x <- grconvertX(pl$lx, "user", "ndc"); ndc_y <- grconvertY(pl$ly, "user", "ndc")
    sidecar[[length(sidecar) + 1L]] <- data.frame(
      sid = sid, group = grp, instance = k, original_text = ck,
      pdf_x = ndc_x * PAGE_W_PT, pdf_y_from_top = (1 - ndc_y) * PAGE_H_PT,
      cx_native = cen[1], cy_native = cen[2], stringsAsFactors = FALSE)
  }
  umpx <- PX_UM / smn; bp <- SCALE_UM / umpx
  segments(cw - bp - 10, ch - 14, cw - 10, ch - 14, lwd = 3, col = "black")
  text(cw - bp / 2 - 10, ch - 28, sprintf("%d um", SCALE_UM), cex = 0.7)
  log_msg("%s: BF %dx%d, %d organoid(s)", sid, cw, ch, length(ids))
}
dev.off()

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
