#!/usr/bin/env Rscript
# ============================================================================
# lib_report_frame.R - shared base-graphics FRAME primitives for the non-locked
#   cache-only report drivers:
#     R/04_ssc_ontissue/05_ssc_report.R
#     R/11_per_organoid_final/10_citrate_gradient_report_final.R
#   Extracted so both drivers share ONE definition of the panel frame (the copies
#   had already begun to drift - colorbar_img gained a `lo` arg in one only).
#   Pure plotting, no SSC/citrate/registration coupling. Output-identical to the
#   former inline copies.
#
#   PMAR             unified symmetric panel margins (L==R -> asp=1 images centre;
#                    identical across panels -> equal boxes + equal row top gap)
#   frame_box()      thin uniform box around the plot region
#   scalebar_bottom  black bar anchored at a fixed npc y BELOW the plot region
#                    (so bars align across a row) + optional caption on the same
#                    line; px_per_um = image pixels per micron along x
#   colorbar_img     vertical colour bar in the right gutter, sized to the image's
#                    y-extent (never taller than the image), labelled lo..hi
# Depends on: nothing (SCALE_UM defined here).
# ============================================================================

SCALE_UM <- 100                                   # scale-bar length (um)
PMAR     <- c(3.4, 3.4, 2.4, 3.4)                 # bottom, left, top, right

frame_box <- function() box("plot", col = "#444444", lwd = 0.8)

scalebar_bottom <- function(px_per_um, cap = NULL) {
  usr <- par("usr"); bar_u <- SCALE_UM * px_per_um
  x1 <- usr[1] + 0.985*(usr[2]-usr[1]); x0 <- x1 - bar_u
  yb <- grconvertY(-0.055, "npc", "user")
  segments(x0, yb, x1, yb, col = "black", lwd = 3, xpd = NA)
  text((x0+x1)/2, grconvertY(-0.125, "npc", "user"), sprintf("%d um", SCALE_UM), col = "black", cex = 0.6, adj = c(0.5, 1), xpd = NA)
  if (!is.null(cap))
    text(usr[1] + 0.015*(usr[2]-usr[1]), yb, cap, col = "grey30", cex = 0.62, adj = c(0, 0.5), xpd = NA)
}

colorbar_img <- function(pal, hi, ylo, yhi, lo = 0) {
  yl <- grconvertY(min(ylo,yhi), "user","npc"); yh <- grconvertY(max(ylo,yhi), "user","npc")
  pad <- 0.05*(yh - yl)                                     # inset so the labelled bar <= image height
  cy0 <- grconvertY(yl + pad, "npc","user"); cy1 <- grconvertY(yh - pad, "npc","user")
  cx0 <- grconvertX(1.02, "npc","user"); cx1 <- grconvertX(1.05, "npc","user")
  yb  <- seq(cy0, cy1, length.out = length(pal)+1); rect(cx0, head(yb,-1), cx1, tail(yb,-1), col = pal, border = NA, xpd = NA)
  text(cx1, c(cy0, cy1), sprintf("%.2g", c(lo, hi)), pos = 4, cex = 0.55, xpd = NA, offset = 0.12)
}
