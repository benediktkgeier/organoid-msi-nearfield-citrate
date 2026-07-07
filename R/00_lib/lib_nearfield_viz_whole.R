#!/usr/bin/env Rscript
# ============================================================================
# lib_nearfield_viz_whole.R - *** LOCKED *** WHOLE-REGION variants of the locked
#   4-view near-field organoid emission figure (lib_nearfield_viz.R).
#   Frozen visual spec; see docs/nearfield_wholeregion.md. Driver: 08_nearfield_wholeregion.R.
# ----------------------------------------------------------------------------
# The locked toolkit renders each view PER ORGANOID (cropped to the organoid
# bbox + MARGIN_UM, masked to that organoid + its Voronoi catchment). This helper
# adds whole-section variants that project all four views across the ENTIRE
# measurement region of one section, drawing EVERY organoid outline coloured by
# its consensus apical class (basolateral-out / apical-out / mixed / unannotated).
#
# It SOURCES the locked lib_nearfield_viz.R and reuses, unchanged:
#   prep_sec(), instance_outlines_native(), arrow_emissions(), inv_affine(),
#   scalebar_native(), scalebar_grid(), nfviz_load_ion(), nfviz_arrow_range(),
#   and constants PAL / WEATHER / CIN / COUT / MARGIN_UM / RINGS_UM.
# The locked lib is NOT modified.
#
# REQUIRED GLOBALS (set by the caller, exactly as for the locked views):
#   val_cit, HI_CIT  via nfviz_load_ion();  GLO, GHI via nfviz_arrow_range().
# Depends on gradient_config.R + lib_register.R + lib_nearfield_viz.R (caller-sourced).
# ============================================================================
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_nearfield_viz.R"))

# consensus apical-class palette (reuses locked CIN/COUT); unannotated -> grey70
CCOL_WHOLE   <- APICAL_COLS   # locked: green basolateral-out / magenta apical-out / grey mixed
UNANNOT_COL  <- "grey70"
WGRID_D      <- 2L                              # native-crop sampling stride (speed)

# MSI-grid y-axis orientation that MATCHES the native-brightfield panels (1 & 3).
# Native panels draw row 1 at top; the MSI->BF affine B[2,2] sign tells how MSI y maps
# onto native Y. B[2,2] < 0 (vertical flip, the registration convention here) -> large
# MSI y sits at the TOP in the native panels, so the grid panel must put max(ys) at top
# (ylim = range(ys)); B[2,2] > 0 -> rev(range(ys)). No registration -> default to range(ys).
ylim_grid <- function(sec) {
  if (is.null(sec$B) || sec$B[2,2] < 0) range(sec$ys) else rev(range(sec$ys))
}

# class colour for one (sid,instance) via an ap_map keyed by paste(sid, instance)
.class_col <- function(sid, k, ap_map) {
  cls <- unname(ap_map[paste(sid, k)])
  if (!is.na(cls) && cls %in% names(CCOL_WHOLE)) CCOL_WHOLE[[cls]] else UNANNOT_COL
}

# ---- draw every organoid outline, coloured by consensus class --------------
# native = TRUE : NATIVE brightfield coords (instance_outlines_native);
# native = FALSE: MSI grid coords (per-instance signed-distance-0 contour).
class_outlines <- function(sec, sid, ap_map, native = TRUE, lwd = 1.8) {
  ids <- sort(unique(sec$z$instance[sec$z$instance > 0]))
  for (k in ids) {
    col <- .class_col(sid, k, ap_map)
    if (native) {
      for (poly in instance_outlines_native(sec, k)) lines(poly$x, poly$y, col = col, lwd = lwd)
    } else {
      sdm_k <- sec$sdm; sdm_k[is.na(sec$instm) | sec$instm != k] <- NA
      if (any(is.finite(sdm_k)))
        contour(sec$xs, sec$ys, sdm_k, levels = 0, drawlabels = FALSE, add = TRUE, col = col, lwd = lwd)
    }
  }
}

# ============================================================================
# WHOLE-REGION NATIVE BF (overlay + heatmap) - generalises locked native_overlay
# ============================================================================
# Shared sampler: project the citrate ion onto the whole native-crop grid and
# (optionally) gblur-smooth it. Returns the value matrix + the sampling context
# needed for ring overlays. Used by BOTH native_overlay_whole (render) and
# heatmap_global_clip (cross-dataset pre-pass) so the two see identical values.
whole_region_vm <- function(sec, smooth_um) {
  if (is.null(sec$om)) return(NULL)
  B <- sec$B; cx0 <- sec$cx0; cy0 <- sec$cy0; smn <- sec$smn; om <- sec$om; cw <- sec$cw; ch <- sec$ch
  W <- max(sec$z$x); H <- max(sec$z$y)
  gx <- seq(1, cw, by = WGRID_D); gy <- seq(1, ch, by = WGRID_D)
  NX <- cx0 + rep(gx, times = length(gy)) - 1; NY <- cy0 + rep(gy, each = length(gx)) - 1
  ms <- inv_affine(B, NX, NY)
  inb <- ms$x >= .5 & ms$x <= W + .5 & ms$y >= .5 & ms$y <= H + .5
  lut <- matrix(NA_real_, W, H); lut[cbind(sec$z$x, sec$z$y)] <- val_cit[sec$z$gidx]
  v <- rep(NA_real_, length(NX))
  if (any(inb)) { xi <- pmin(pmax(round(ms$x[inb]),1),W); yi <- pmin(pmax(round(ms$y[inb]),1),H); v[inb] <- lut[cbind(xi,yi)] }
  vm <- matrix(v, length(gy), length(gx), byrow = TRUE); mm <- is.finite(vm)
  if (smooth_um > 0 && any(mm)) {
    sig <- max(1, smooth_um * smn / (MSI_PIXEL_UM * WGRID_D)); v0 <- vm; v0[!mm] <- 0
    vm <- as.matrix(EBImage::gblur(v0, sigma = sig)); vm[!mm] <- NA
  }
  list(vm = vm, gx = gx, gy = gy, om = om, cw = cw, ch = ch, smn = smn,
       NX = NX, ms = ms, inb = inb, W = W, H = H, z = sec$z)
}

# GLOBAL heatmap clip: pool the smoothed whole-region values across `sids` and
# take a single percentile so the view-3 heatmap is comparable page-to-page.
# LOCKED clip = pooled p99.9 (q = 0.999): dimmer than the old p99 (which over-saturated
# the top 1% to red); keeps genuine hotspots warm while cooling the surround.
heatmap_global_clip <- function(sids, smooth_um = 12, q = 0.999) {
  vals <- unlist(lapply(sids, function(s) {
    r <- whole_region_vm(prep_sec(s), smooth_um); if (is.null(r)) return(numeric(0))
    r$vm[is.finite(r$vm)]
  }))
  hi <- if (length(vals)) as.numeric(quantile(vals, q)) else NA_real_
  cat(sprintf("[nfviz] GLOBAL heatmap clip (p%.1f over %d datasets, smoothed): %.3f\n",
              q*100, length(sids), hi))
  hi
}

# pal_cbar: small vertical colour bar in the right gutter (0..hi), for the global heatmap
.pal_cbar <- function(pal, hi, lab = "citrate") {
  cx0<-grconvertX(1.015,"npc","user"); cx1<-grconvertX(1.05,"npc","user")
  cy0<-grconvertY(0.06,"npc","user");  cy1<-grconvertY(0.94,"npc","user")
  yb<-seq(cy0,cy1,length.out=length(pal)+1)
  rect(cx0,head(yb,-1),cx1,tail(yb,-1),col=pal,border=NA,xpd=NA)
  text(cx1,grconvertY(c(0.06,0.94),"npc","user"),sprintf("%.2g",c(0,hi)),pos=4,cex=0.5,xpd=NA,offset=0.1)
  text(grconvertX(1.03,"npc","user"),grconvertY(0.99,"npc","user"),lab,cex=0.5,xpd=NA,adj=c(0.5,0))
}

native_overlay_whole <- function(sec, sid, ap_map, pal, alpha, smooth_um, main = "", hi = HI_CIT, cbar = FALSE) {
  if (is.null(sec$om)) { plot.new(); title(main); return(invisible()) }
  r <- whole_region_vm(sec, smooth_um)
  vm <- r$vm; gx <- r$gx; gy <- r$gy; om <- r$om; cw <- r$cw; ch <- r$ch; smn <- r$smn
  NX <- r$NX; ms <- r$ms; inb <- r$inb; W <- r$W; H <- r$H
  hi_use <- if (is.na(hi)) { fv <- vm[is.finite(vm)]; if (length(fv)) as.numeric(quantile(fv, 0.999)) else 1 } else hi   # LOCKED p99.9 (was 0.99)
  vs <- pmin(vm / hi_use, 1)
  par(mar = c(2.0, 1, 2.4, if (cbar) 3.4 else 1)); plot.new(); plot.window(c(1, cw), c(ch, 1), asp = 1)
  rasterImage(om / max(om), 1, ch, cw, 1, interpolate = TRUE); title(main, cex.main = 0.95)
  fin <- is.finite(vs)
  if (any(fin)) { cidx <- pmin(pmax(round(vs*255)+1,1),256); cmat <- matrix("#00000000", length(gy), length(gx))
    cmat[fin] <- adjustcolor(pal[cidx[fin]], alpha.f = alpha)
    rasterImage(as.raster(cmat), 1, ch, cw, 1, interpolate = (smooth_um > 0)) }
  # whole-region 50/100 um signed-distance rings (0 um = organoid outlines, drawn separately)
  lutd <- matrix(NA_real_, W, H); lutd[cbind(r$z$x, r$z$y)] <- r$z$signed_dist_um
  dn <- rep(NA_real_, length(NX)); if (any(inb)) { xi <- pmin(pmax(round(ms$x[inb]),1),W); yi <- pmin(pmax(round(ms$y[inb]),1),H); dn[inb] <- lutd[cbind(xi,yi)] }
  dmat <- matrix(dn, length(gy), length(gx), byrow = TRUE)
  if (any(is.finite(dmat)))
    contour(gx, gy, t(dmat), levels = RINGS_UM[RINGS_UM > 0], add = TRUE, drawlabels = FALSE,
            col = adjustcolor("white", 0.7), lwd = 0.7, lty = 3)
  class_outlines(sec, sid, ap_map, native = TRUE, lwd = 1.8)
  if (cbar && is.finite(hi_use)) .pal_cbar(pal, hi_use)
  scalebar_native(1, cw, 1, ch, smn)
}
draw_overlay_whole <- function(sec, sid, ap_map, main = "")
  native_overlay_whole(sec, sid, ap_map, PAL, 0.60, 0, main, hi = HI_CIT)
# hi = NA  -> per-panel p99 auto-scale (version 1, each section to its own max)
# hi = num -> fixed clip shared across pages (version 2, global cross-dataset scale) + colorbar
draw_heatmap_whole <- function(sec, sid, ap_map, main = "", hi = NA)
  native_overlay_whole(sec, sid, ap_map, WEATHER, 0.72, 12, main, hi = hi, cbar = is.finite(hi))

# ============================================================================
# WHOLE-REGION GRADIENT MAP (MSI grid)
# ============================================================================
draw_gradmap_whole <- function(sec, sid, ap_map, main = "") {
  sc <- pmin(sec$citm / HI_CIT, 1)
  par(mar = c(2.0, 1, 2.4, 1)); plot.new(); plot.window(range(sec$xs), ylim_grid(sec), asp = 1)
  rect(min(sec$xs)-.5, min(sec$ys)-.5, max(sec$xs)+.5, max(sec$ys)+.5, col = "#101010", border = NA)
  image(sec$xs, sec$ys, sc, col = PAL, zlim = c(0,1), add = TRUE, useRaster = TRUE)
  title(main, cex.main = 0.95)
  if (any(is.finite(sec$sdm)))
    contour(sec$xs, sec$ys, sec$sdm, levels = RINGS_UM[RINGS_UM > 0], drawlabels = FALSE, add = TRUE,
            col = adjustcolor("white", 0.8), lwd = 1.0, lty = 3)
  class_outlines(sec, sid, ap_map, native = FALSE, lwd = 1.6)
  scalebar_grid(sec$xs, sec$ys)
}

# ============================================================================
# WHOLE-REGION EMISSION VECTORS (MSI grid) - arrows from EVERY organoid surface
# ============================================================================
draw_vectors_whole <- function(sec, sid, ap_map, main = "") {
  sc <- pmin(sec$citm / HI_CIT, 1)
  par(mar = c(2.0, 1, 2.4, 1)); plot.new(); plot.window(range(sec$xs), ylim_grid(sec), asp = 1)
  rect(min(sec$xs)-.5, min(sec$ys)-.5, max(sec$xs)+.5, max(sec$ys)+.5, col = "#101010", border = NA)
  image(sec$xs, sec$ys, sc, col = PAL, zlim = c(0,1), add = TRUE, useRaster = TRUE)
  title(main, cex.main = 0.95)
  class_outlines(sec, sid, ap_map, native = FALSE, lwd = 1.4)
  ids <- sort(unique(sec$z$instance[sec$z$instance > 0]))
  ae <- do.call(rbind, lapply(ids, function(k) arrow_emissions(sec, k)))
  if (!is.null(ae) && nrow(ae)) {
    es <- pmin(pmax(ae$emis/HI_CIT, 0), 1)                                   # absolute -> LENGTH
    L  <- (0.35 + 0.65*es) * (MARGIN_UM*0.75/MSI_PIXEL_UM)
    srel <- if (is.finite(GHI) && GHI > GLO) pmin(pmax((ae$emis - GLO)/(GHI - GLO), 0), 1) else rep(0.5, nrow(ae))
    f  <- 0.25 + 1.75*srel                                                   # global relative score -> THICKNESS (0.25x..2.0x)
    for (j in seq_len(nrow(ae))) {
      x1 <- ae$x0[j] + ae$nx[j]*L[j]; y1 <- ae$y0[j] + ae$ny[j]*L[j]
      arrows(ae$x0[j], ae$y0[j], x1, y1, length = 0.05, lwd = 4.2*f[j], col = "black")
      arrows(ae$x0[j], ae$y0[j], x1, y1, length = 0.04, lwd = 2.2*f[j], col = "white")
    }
  }
  scalebar_grid(sec$xs, sec$ys)
}
