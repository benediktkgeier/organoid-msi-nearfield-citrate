#!/usr/bin/env Rscript
# ============================================================================
# lib_gradient_report.R - shared per-dataset PANEL renderers for the final
#   citrate-gradient report family:
#     R/11_per_organoid_final/10_citrate_gradient_report_final.R   (mixed as grey context)
#     R/11_per_organoid_final/13_citrate_gradient_report_3class.R  (mixed own trend / nomixed)
#   Extracted so the two drivers share ONE definition of the 6 per-dataset panels
#   (+ sec_data). The drivers differ only in the overview mini-plot and the page
#   title/legend, which stay in each driver. Output-identical to the former inline
#   copies (which mirrored the LOCKED 09 panels).
#
# NOT self-contained: the CALLER must, before rendering, have sourced
#   lib_report_frame.R (PMAR / frame_box / scalebar_bottom / colorbar_img / SCALE_UM),
#   lib_nearfield_viz*.R (PAL / WEATHER / RINGS_UM / class_outlines / whole_region_vm /
#   ylim_grid), lib_register.R (apply_affine), if_config.R (BF_UMPX);
#   and set the globals: pd, val_cit, HI_CIT, HEAT_HI, PX_UM, REG_CACHE, BF_CROPS,
#   IF_CROPS, EXT_IF. (All are resolved at call time, so source order is flexible.)
# ============================================================================

# citrate image grid for one section
sec_data <- function(sid) {
  m <- which(as.character(pd$sample_id) == sid); x <- pd$x[m]; y <- pd$y[m]
  xs <- sort(unique(x)); ys <- sort(unique(y)); ix <- match(x, xs); iy <- match(y, ys)
  mm <- matrix(0, length(xs), length(ys)); mm[cbind(ix, iy)] <- val_cit[m]
  list(xs = xs, ys = ys, cit = mm)
}

# Panel 1: MSI citrate ion image (viridis, colorbar sized to image, aligned scale bar)
ion_panel <- function(mat, xs, ys, hi, main, fc, colorbar, cap, ylim = range(ys)) {
  sc <- pmin(mat / hi, 1); par(mar = PMAR)
  image(xs, ys, sc, col = PAL, asp = 1, useRaster = TRUE, zlim = c(0,1), ylim = ylim,
        axes = FALSE, xlab = "", ylab = "", main = main, cex.main = 0.95, col.main = fc)
  frame_box()
  if (colorbar) colorbar_img(PAL, hi, min(ys), max(ys))
  scalebar_bottom(1/PX_UM, cap = cap)
}

# Panel 2: native brightfield crop (optical_<sid>.png)
native_bf_panel <- function(sid, main, cap) {
  xf <- file.path(REG_CACHE, sprintf("nd2final_%s.rds", sid)); pf <- file.path(BF_CROPS, sprintf("optical_%s.png", sid))
  par(mar = PMAR)
  if (!file.exists(pf)) { plot.new(); title(main, cex.main = 0.95); text(0.5,0.5,"(no native crop)",col="grey50"); frame_box(); return(invisible()) }
  om <- png::readPNG(pf); if (length(dim(om))==3) om <- om[,,1]; cw <- ncol(om); ch <- nrow(om)
  plot.new(); plot.window(c(1,cw), c(ch,1), asp = 1); rasterImage(om/max(om), 1, ch, cw, 1); title(main, cex.main = 0.95)
  frame_box()
  smn <- if (file.exists(xf)) readRDS(xf)$scale_msi_nd2 else 1
  scalebar_bottom(smn/PX_UM, cap = cap)          # native px per um = smn/PX_UM
}

# Panel 3: SSC on-tissue mask (white = tissue) - is_tissue column of TISSUE_MSE
ssc_mask_panel <- function(sid, sec, main, cap) {
  sel <- which(as.character(pd$sample_id) == sid); x <- pd$x[sel]; y <- pd$y[sel]
  xs <- sort(unique(x)); ys <- sort(unique(y))
  m <- matrix(NA_real_, length(xs), length(ys)); m[cbind(match(x,xs), match(y,ys))] <- as.integer(pd$is_tissue[sel]) + 1L
  par(mar = PMAR)
  image(xs, ys, m, col = c("grey30","white"), zlim = c(0.5,2.5), asp = 1, useRaster = TRUE, ylim = ylim_grid(sec),
        axes = FALSE, xlab = "", ylab = "", main = main, cex.main = 0.95)
  frame_box(); scalebar_bottom(1/PX_UM, cap = cap)
}

# Panel 4: matched, clipped IF composite (pre-rendered if_optical_<sid>.png), +30% brightness
if_panel <- function(sid, main, cap) {
  par(mar = PMAR)
  pf <- file.path(IF_CROPS, sprintf("if_optical_%s.png", sid))
  if (!file.exists(pf)) { plot.new(); title(main, cex.main = 0.95); text(0.5,0.5,"(no IF)",col="grey50"); frame_box(); return(invisible()) }
  im <- png::readPNG(pf); im <- pmin(im * 1.30, 1)             # +30% brightness gain (background stays ~black)
  outH <- dim(im)[1]; outW <- dim(im)[2]
  plot.new(); plot.window(c(1,outW), c(outH,1), asp = 1); rasterImage(im, 1, outH, outW, 1, interpolate = TRUE); title(main, cex.main = 0.95)
  frame_box()
  # scale bar: reconstruct the BF-crop geometry (Bm + MSI bbox + 30% extension), um/px = crop_um_width/outW
  xf <- file.path(REG_CACHE, sprintf("nd2final_%s.rds", sid)); px_per_um <- NA_real_
  if (file.exists(xf)) {
    Bm <- readRDS(xf)$B_msi_nd2; sub <- pd[as.character(pd$sample_id)==sid, c("x","y")]
    nd <- apply_affine(Bm, cbind(sub$x, sub$y)); rx <- range(nd[,1]); ex <- diff(rx)*EXT_IF/2
    cx0 <- max(1, floor(rx[1]-ex)); cx1 <- ceiling(rx[2]+ex); cw_bf <- cx1 - cx0 + 1
    umpx <- cw_bf*BF_UMPX/outW; px_per_um <- 1/umpx
  }
  if (is.finite(px_per_um)) scalebar_bottom(px_per_um, cap = cap) else
    text(outW/2, grconvertY(-0.055,"npc","user"), cap, col = "grey30", cex = 0.62, adj = c(0.5,0.5), xpd = NA)
}

# Panel 5: MSI<->BF overlay WITHOUT the 50/100 um rings (copy of native_overlay_whole minus the ring block)
native_overlay_norings <- function(sec, sid, apm, main, cap) {
  if (is.null(sec$om)) { par(mar = PMAR); plot.new(); title(main, cex.main = 0.95); frame_box(); return(invisible()) }
  r <- whole_region_vm(sec, 0)                                  # no smoothing (overlay)
  vm <- r$vm; om <- r$om; cw <- r$cw; ch <- r$ch; smn <- r$smn
  vs <- pmin(vm / HI_CIT, 1)
  par(mar = PMAR); plot.new(); plot.window(c(1, cw), c(ch, 1), asp = 1)
  rasterImage(om / max(om), 1, ch, cw, 1, interpolate = TRUE); title(main, cex.main = 0.95)
  fin <- is.finite(vs)
  if (any(fin)) { cidx <- pmin(pmax(round(vs*255)+1,1),256); cmat <- matrix("#00000000", length(r$gy), length(r$gx))
    cmat[fin] <- adjustcolor(PAL[cidx[fin]], alpha.f = 0.60); rasterImage(as.raster(cmat), 1, ch, cw, 1, interpolate = FALSE) }
  class_outlines(sec, sid, apm, native = TRUE, lwd = 1.8)       # NO rings drawn
  frame_box(); scalebar_bottom(smn/PX_UM, cap = cap)
}

# Panel 6: weather-rainbow gradient heatmap on BF (gblur 12um, GLOBAL clip) + 50/100 um rings.
weather_panel <- function(sec, sid, apm, main, cap) {
  if (is.null(sec$om)) { par(mar = PMAR); plot.new(); title(main, cex.main = 0.95); frame_box(); return(invisible()) }
  r <- whole_region_vm(sec, 12)                                 # 12 um gblur (locked heatmap smoothing)
  vm <- r$vm; om <- r$om; cw <- r$cw; ch <- r$ch; smn <- r$smn
  vs <- pmin(vm / HEAT_HI, 1)
  par(mar = PMAR); plot.new(); plot.window(c(1, cw), c(ch, 1), asp = 1)
  rasterImage(om / max(om), 1, ch, cw, 1, interpolate = TRUE); title(main, cex.main = 0.95)
  fin <- is.finite(vs)
  if (any(fin)) { cidx <- pmin(pmax(round(vs*255)+1,1),256); cmat <- matrix("#00000000", length(r$gy), length(r$gx))
    cmat[fin] <- adjustcolor(WEATHER[cidx[fin]], alpha.f = 0.72); rasterImage(as.raster(cmat), 1, ch, cw, 1, interpolate = TRUE) }
  # 50/100 um signed-distance rings (0 um = class-coloured organoid outlines, drawn after)
  lutd <- matrix(NA_real_, r$W, r$H); lutd[cbind(r$z$x, r$z$y)] <- r$z$signed_dist_um
  dn <- rep(NA_real_, length(r$NX)); if (any(r$inb)) { xi <- pmin(pmax(round(r$ms$x[r$inb]),1),r$W); yi <- pmin(pmax(round(r$ms$y[r$inb]),1),r$H); dn[r$inb] <- lutd[cbind(xi,yi)] }
  dmat <- matrix(dn, length(r$gy), length(r$gx), byrow = TRUE)
  if (any(is.finite(dmat)))
    contour(r$gx, r$gy, t(dmat), levels = RINGS_UM[RINGS_UM > 0], add = TRUE, drawlabels = FALSE, col = adjustcolor("white", 0.7), lwd = 0.7, lty = 3)
  class_outlines(sec, sid, apm, native = TRUE, lwd = 1.8)
  frame_box(); colorbar_img(WEATHER, HEAT_HI, 1, ch); scalebar_bottom(smn/PX_UM, cap = cap)
}
