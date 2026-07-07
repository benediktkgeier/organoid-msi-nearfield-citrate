#!/usr/bin/env Rscript
# ============================================================================
# lib_nearfield_viz.R - LOCKED near-field organoid emission visualization toolkit
# ----------------------------------------------------------------------------
# Reusable, frozen implementation of the four-view near-field figure used to show
# per-organoid emission of an ion into the surrounding gel (developed for citrate
# [M-H]- 191.0217, apical-out vs basolateral-out; see docs/apical_nearfield.md).
# Source this from any report so the visual style stays identical across figures.
#
# THE FOUR LOCKED VIEWS (each cropped to the organoid bbox + MARGIN_UM, and MASKED
# to the organoid interior + its Voronoi catchment = the region near-field metrics
# are measured over; neighbours are excluded):
#   draw_overlay()  VIEW 1  ion MSI on native brightfield, viridis, constant alpha
#                            0.60; white organoid outline + dotted 50/100 um rings.
#   draw_gradmap()  VIEW 2  cropped ion image (viridis, global p99.5 clip = HI_CIT)
#                            + solid 0 / dotted 50,100 um signed-distance rings.
#   draw_heatmap()  VIEW 3  interpolated "weather-map" heatmap on brightfield:
#                            EBImage::gblur-smoothed ion, WEATHER rainbow
#                            (LOW=blue/cold -> HIGH=red/hot). LOCKED clip = p99.9
#                            (per-organoid: its own p99.9; whole-region: pooled global
#                            p99.9) - dimmer than the old p99, hotspots stay warm.
#   draw_vectors()  VIEW 4  outward arrows from surface points (normal = grad of
#                            signed distance). LENGTH ~ absolute near-field ion
#                            (sampled 20-50 um out); THICKNESS ~ relative score
#                            scaled GLOBALLY across exemplars (p10..p90 of emission)
#                            spanning 0.25x..2.0x (= -75%..+100%) of base lwd.
# Conventions: viridis for ion images; scale bar BELOW image, centred label
# beneath; basolateral-out = CIN (green), apical-out = COUT (magenta).
#
# REQUIRED GLOBALS the caller must set before drawing (helpers provided):
#   val_cit  - numeric ion-intensity vector indexed by MSE pixel gidx  (nfviz_load_ion)
#   HI_CIT   - global p99.5 clip for the ion images                    (nfviz_load_ion)
#   GLO,GHI  - global emission range for arrow thickness               (nfviz_arrow_range)
# Depends on gradient_config.R (CACHE_DIR, FIG_DIR, MSI_PIXEL_UM, IMG_CLIP_HI,
# TISSUE_MSE) and lib_register.R (apply_affine), both sourced by the caller.
# ============================================================================
suppressPackageStartupMessages({ library(Cardinal); library(png); library(EBImage); library(viridisLite) })

# ---- LOCKED visual constants ----------------------------------------------
CIT_MZ    <- 191.0217; DHA_MZ <- 327.2330   # CIT_MZ = legacy label; citrate is extracted ANCHORED (raw imzML, +-CITRATE_WIN_PPM around CITRATE_ANCHOR_MZ) via lib_citrate, NOT this grid feature
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_citrate.R"))
MARGIN_UM <- 120; RINGS_UM <- c(0, 50, 100)
PAL       <- viridisLite::viridis(256)                                   # ion images
WEATHER   <- colorRampPalette(c("#00008b","blue","cyan","green3","yellow","orange","red","#7f0000"))(256)  # heatmap: low=blue/cold -> high=red/hot
CIN       <- unname(APICAL_COLS["basolateral_out"]); COUT <- unname(APICAL_COLS["apical_out"])  # basolateral-out / apical-out accents (locked)
REG_CACHE <- cache_in("register"); CROP_DIR <- file.path(FIG_DIR, "registration", "crops")

inv_affine <- function(B, NX, NY) { M <- rbind(B[1,], B[2,]); t0 <- B[3,]
  xy <- solve(M, rbind(NX - t0[1], NY - t0[2])); list(x = xy[1,], y = xy[2,]) }

# ---- initialisers (set the required globals) -------------------------------
# Load an ion from the working MSE into the globals the views need.
nfviz_load_ion <- function(mz_target = CIT_MZ, mse_path = TISSUE_MSE) {
  mse <- readRDS(mse_path); pd <<- as.data.frame(pixelData(mse)); pd$gidx <<- seq_len(nrow(pd))
  mzs <- mz(mse); idx <- which.min(abs(mzs - mz_target))
  if (abs((mz_target - CITRATE_ANCHOR_MZ)/CITRATE_ANCHOR_MZ*1e6) < 50) {   # CITRATE -> anchored raw imzML window
    rm(mse); gc(verbose = FALSE)
    val_cit <<- citrate_onto_pd(pd)            # +-CITRATE_WIN_PPM around the standard anchor, TIC-norm
    cat(sprintf("[nfviz] citrate ANCHORED %.5f +-%d ppm (raw imzML, TIC-norm); ", CITRATE_ANCHOR_MZ, CITRATE_WIN_PPM))
  } else {                                                                 # any other ion -> grid feature
    val_cit <<- as.numeric(spectra(mse)[idx, ]); rm(mse); gc(verbose = FALSE)
    cat(sprintf("[nfviz] grid ion idx=%d mz=%.4f (%.1f ppm); ", idx, mzs[idx], (mzs[idx]-mz_target)/mz_target*1e6))
  }
  HI_CIT <<- as.numeric(quantile(val_cit[val_cit > 0], IMG_CLIP_HI, na.rm = TRUE))
  cat(sprintf("HI_CIT(p%.1f)=%.3f\n", IMG_CLIP_HI*100, HI_CIT))
  invisible(idx)
}
# Set the GLOBAL arrow-thickness emission range from a set of exemplar organoids
# (data.frame with sid, instance). Robust p10..p90 so one hot pixel can't compress it.
nfviz_arrow_range <- function(vex) {
  ae_all <- unlist(lapply(seq_len(nrow(vex)), function(i) {
    ae <- arrow_emissions(prep_sec(vex$sid[i]), vex$instance[i]); if (nrow(ae)) ae$emis else numeric(0) }))
  GLO <<- if (length(ae_all)) as.numeric(quantile(ae_all, 0.10)) else 0
  GHI <<- if (length(ae_all)) as.numeric(quantile(ae_all, 0.90)) else 1
  cat(sprintf("[nfviz] global arrow emission range (thickness scale, p10..p90): %.3f .. %.3f\n", GLO, GHI))
  invisible(c(GLO, GHI))
}

# ---- per-section loader (cached) + organoid crop ---------------------------
SEC <- new.env()
prep_sec <- function(sid) {
  if (!is.null(SEC[[sid]])) return(SEC[[sid]])
  z <- readRDS(cache_in(sprintf("zones_%s.rds", sid)))
  xs <- sort(unique(z$x)); ys <- sort(unique(z$y)); ix <- match(z$x, xs); iy <- match(z$y, ys)
  fill <- function(v) { m <- matrix(NA_real_, length(xs), length(ys)); m[cbind(ix, iy)] <- v; m }
  obj <- list(z = z, xs = xs, ys = ys,
              citm = fill(val_cit[z$gidx]), sdm = fill(z$signed_dist_um),
              instm = fill(as.numeric(z$instance)), catchm = fill(as.numeric(z$instance_catch)),
              surfm = fill(as.numeric(z$is_surface)))
  xf <- file.path(REG_CACHE, sprintf("nd2final_%s.rds", sid)); pf <- file.path(CROP_DIR, sprintf("optical_%s.png", sid))
  if (file.exists(xf) && file.exists(pf)) {
    X <- readRDS(xf); om <- png::readPNG(pf); if (length(dim(om)) == 3) om <- om[,,1]
    obj$B <- X$B_msi_nd2; obj$smn <- X$scale_msi_nd2; obj$cx0 <- X$crop[1]; obj$cy0 <- X$crop[2]
    obj$om <- om; obj$cw <- ncol(om); obj$ch <- nrow(om)
  }
  SEC[[sid]] <- obj; obj
}
prep_org <- function(sec, k) {
  mp <- MARGIN_UM / MSI_PIXEL_UM
  oc <- which(sec$z$instance == k); ox <- sec$z$x[oc]; oy <- sec$z$y[oc]
  xi <- which(sec$xs >= min(ox)-mp & sec$xs <= max(ox)+mp); yi <- which(sec$ys >= min(oy)-mp & sec$ys <= max(oy)+mp)
  instc <- sec$instm[xi, yi, drop = FALSE]; catchc <- sec$catchm[xi, yi, drop = FALSE]
  maskc <- (!is.na(instc) & instc == k) | (!is.na(catchc) & catchc == k)
  list(k = k, xi = xi, yi = yi, xsc = sec$xs[xi], ysc = sec$ys[yi],
       citc = sec$citm[xi, yi, drop = FALSE], sdc = sec$sdm[xi, yi, drop = FALSE],
       instc = instc, catchc = catchc, surfc = sec$surfm[xi, yi, drop = FALSE], maskc = maskc)
}
scalebar_grid <- function(xsc, ysc, um = 50) {
  bp <- um / MSI_PIXEL_UM; x1 <- max(xsc) - bp; yb <- min(ysc) - 0.06*diff(range(ysc)) - 0.4
  par(xpd = NA); segments(x1, yb, max(xsc), yb, lwd = 3, col = "black")
  text((x1 + max(xsc))/2, yb - 0.07*diff(range(ysc)) - 0.4, sprintf("%d um", um), cex = 0.6); par(xpd = FALSE)
}
scalebar_native <- function(nx0, nx1, ny0, ny1, smn, um = 50) {
  bp <- um * smn / MSI_PIXEL_UM; x1 <- nx1 - bp; yb <- ny1 + 0.06*(ny1-ny0)
  par(xpd = NA); segments(x1, yb, nx1, yb, lwd = 3, col = "black")
  text((x1+nx1)/2, yb + 0.06*(ny1-ny0), sprintf("%d um", um), cex = 0.6); par(xpd = FALSE)
}
instance_outlines_native <- function(sec, k) {
  M <- (sec$instm == k); M[is.na(M)] <- FALSE
  M <- EBImage::fillHull(EBImage::closing(M, EBImage::makeBrush(3, "disc")))
  out <- list()
  for (co in EBImage::ocontour(EBImage::bwlabel(M))) { if (nrow(co) < 6) next
    gx <- sec$xs[pmin(pmax(co[,1],1),length(sec$xs))]; gy <- sec$ys[pmin(pmax(co[,2],1),length(sec$ys))]
    p <- apply_affine(sec$B, cbind(gx, gy))
    out[[length(out)+1L]] <- list(x = c(p[,1]-sec$cx0+1, p[1,1]-sec$cx0+1), y = c(p[,2]-sec$cy0+1, p[1,2]-sec$cy0+1)) }
  out
}
draw_rings_native <- function(sec, k, levels, cx0, cy0, B) {
  hm <- (!is.na(sec$instm) & sec$instm == k) | (!is.na(sec$catchm) & sec$catchm == k)
  sdm_h <- sec$sdm; sdm_h[!hm] <- NA
  cl <- grDevices::contourLines(sec$xs, sec$ys, sdm_h, levels = levels[levels > 0])
  for (L in cl) { p <- apply_affine(B, cbind(L$x, L$y))
    lines(p[,1]-cx0+1, p[,2]-cy0+1, col = adjustcolor("white", 0.8), lwd = 1.0, lty = 3) }
}

# ============================ VIEW 1 + VIEW 3 (native BF) ====================
native_overlay <- function(sec, k, pal, alpha, smooth_um, main = "", hi = HI_CIT) {
  if (is.null(sec$om)) { plot.new(); title(main); return(invisible()) }
  B <- sec$B; cx0 <- sec$cx0; cy0 <- sec$cy0; smn <- sec$smn; om <- sec$om; cw <- sec$cw; ch <- sec$ch
  oc <- which(sec$z$instance == k)
  np <- apply_affine(B, cbind(sec$z$x[oc], sec$z$y[oc]))
  mum <- MARGIN_UM * smn / MSI_PIXEL_UM
  nx0 <- max(1, floor(min(np[,1]) - cx0 + 1 - mum)); nx1 <- min(cw, ceiling(max(np[,1]) - cx0 + 1 + mum))
  ny0 <- max(1, floor(min(np[,2]) - cy0 + 1 - mum)); ny1 <- min(ch, ceiling(max(np[,2]) - cy0 + 1 + mum))
  omc <- om[ny0:ny1, nx0:nx1, drop = FALSE]
  gx <- seq(nx0, nx1, by = 1); gy <- seq(ny0, ny1, by = 1)
  NX <- cx0 + rep(gx, times = length(gy)) - 1; NY <- cy0 + rep(gy, each = length(gx)) - 1
  ms <- inv_affine(B, NX, NY); W <- max(sec$z$x); H <- max(sec$z$y)
  inb <- ms$x >= .5 & ms$x <= W+.5 & ms$y >= .5 & ms$y <= H+.5
  lut <- matrix(NA_real_, W, H); lut[cbind(sec$z$x, sec$z$y)] <- val_cit[sec$z$gidx]
  mlut <- matrix(FALSE, W, H); mlut[cbind(sec$z$x, sec$z$y)] <- (sec$z$instance == k) | (sec$z$instance_catch == k)
  v <- rep(NA_real_, length(NX))
  if (any(inb)) { xi <- pmin(pmax(round(ms$x[inb]),1),W); yi <- pmin(pmax(round(ms$y[inb]),1),H)
    vv <- lut[cbind(xi,yi)]; vv[!mlut[cbind(xi,yi)]] <- NA; v[inb] <- vv }
  vm <- matrix(v, length(gy), length(gx), byrow = TRUE); mm <- is.finite(vm)
  if (smooth_um > 0 && any(mm)) {
    sig <- max(1, smooth_um * smn / MSI_PIXEL_UM); v0 <- vm; v0[!mm] <- 0
    vm <- as.matrix(EBImage::gblur(v0, sigma = sig)); vm[!mm] <- NA
  }
  hi_use <- if (is.na(hi)) { fv <- vm[is.finite(vm)]; if (length(fv)) as.numeric(quantile(fv, 0.999)) else 1 } else hi   # LOCKED p99.9 own-scale (was 0.99)
  vs <- pmin(vm/hi_use, 1)
  par(mar = c(2.0, 1, 2.4, 1)); plot.new(); plot.window(c(nx0, nx1), c(ny1, ny0), asp = 1)
  rasterImage(omc/max(omc), nx0, ny1, nx1, ny0, interpolate = TRUE); title(main, cex.main = 0.95)
  fin <- is.finite(vs)
  if (any(fin)) { cidx <- pmin(pmax(round(vs*255)+1,1),256); cmat <- matrix("#00000000", length(gy), length(gx))
    cmat[fin] <- adjustcolor(pal[cidx[fin]], alpha.f = alpha)
    rasterImage(as.raster(cmat), nx0, ny1, nx1, ny0, interpolate = (smooth_um > 0)) }
  for (poly in instance_outlines_native(sec, k)) lines(poly$x, poly$y, col = "white", lwd = 1.8)
  draw_rings_native(sec, k, levels = RINGS_UM, cx0, cy0, B)
  scalebar_native(nx0, nx1, ny0, ny1, smn)
}
draw_overlay <- function(sec, k, main = "") native_overlay(sec, k, PAL, 0.60, 0, main, hi = HI_CIT)
draw_heatmap <- function(sec, k, main = "") native_overlay(sec, k, WEATHER, 0.72, 12, main, hi = NA)

# ============================ VIEW 2: GRADIENT MAP (MSI grid) ===============
draw_gradmap <- function(sec, k, main = "") {
  o <- prep_org(sec, k); sc <- pmin(o$citc/HI_CIT, 1); sc[!o$maskc] <- NA
  par(mar = c(2.0, 1, 2.4, 1)); plot.new(); plot.window(range(o$xsc), rev(range(o$ysc)), asp = 1)
  rect(min(o$xsc)-.5, min(o$ysc)-.5, max(o$xsc)+.5, max(o$ysc)+.5, col = "#101010", border = NA)
  image(o$xsc, o$ysc, sc, col = PAL, zlim = c(0,1), add = TRUE, useRaster = TRUE)
  title(main, cex.main = 0.95)
  sdm_c <- o$sdc; sdm_c[!o$maskc] <- NA
  contour(o$xsc, o$ysc, sdm_c, levels = RINGS_UM, drawlabels = FALSE, add = TRUE,
          col = adjustcolor("white", 0.85), lwd = c(1.6,1.0,1.0), lty = c(1,3,3))
  scalebar_grid(o$xsc, o$ysc)
}

# ============================ VIEW 4: VECTORS ==============================
arrow_emissions <- function(sec, k) {
  sd <- sec$sdm; gx <- matrix(NA_real_, nrow(sd), ncol(sd)); gy <- gx
  gx[2:(nrow(sd)-1), ] <- (sd[3:nrow(sd), ] - sd[1:(nrow(sd)-2), ]) / 2
  gy[, 2:(ncol(sd)-1)] <- (sd[, 3:ncol(sd)] - sd[, 1:(ncol(sd)-2)]) / 2
  cit <- sec$citm
  sp <- which(sec$surfm == 1 & sec$instm == k, arr.ind = TRUE)
  if (nrow(sp) <= 2) return(data.frame())
  cx <- mean(sp[,1]); cy <- mean(sp[,2]); ang <- atan2(sp[,2]-cy, sp[,1]-cx)
  sp <- sp[order(ang), , drop = FALSE]
  sel <- unique(round(seq(1, nrow(sp), length.out = min(28, nrow(sp))))); sp <- sp[sel, , drop = FALSE]
  arr_um <- c(2, 3, 4, 5)                            # sample 20-50 um outward
  rows <- list()
  for (i in seq_len(nrow(sp))) {
    a <- sp[i,1]; b <- sp[i,2]; n <- c(gx[a,b], gy[a,b]); nn <- sqrt(sum(n^2)); if (!is.finite(nn) || nn == 0) next
    n <- n/nn
    pts <- sapply(arr_um, function(d) { ia <- round(a + n[1]*d); ib <- round(b + n[2]*d)
      if (ia>=1 && ia<=nrow(cit) && ib>=1 && ib<=ncol(cit)) cit[ia,ib] else NA })
    e <- mean(pts, na.rm = TRUE)
    if (is.finite(e) && e > 0) rows[[length(rows)+1L]] <- data.frame(x0 = sec$xs[a], y0 = sec$ys[b], nx = n[1], ny = n[2], emis = e)
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}
draw_vectors <- function(sec, k, main = "") {
  o <- prep_org(sec, k); sc <- pmin(o$citc/HI_CIT, 1); sc[!o$maskc] <- NA
  par(mar = c(2.0, 1, 2.4, 1)); plot.new(); plot.window(range(o$xsc), rev(range(o$ysc)), asp = 1)
  rect(min(o$xsc)-.5, min(o$ysc)-.5, max(o$xsc)+.5, max(o$ysc)+.5, col = "#101010", border = NA)
  image(o$xsc, o$ysc, sc, col = PAL, zlim = c(0,1), add = TRUE, useRaster = TRUE)
  sdm_c <- o$sdc; sdm_c[!o$maskc] <- NA
  contour(o$xsc, o$ysc, sdm_c, levels = 0, drawlabels = FALSE, add = TRUE, col = "white", lwd = 1.4)
  title(main, cex.main = 0.95)
  ae <- arrow_emissions(sec, k)
  if (nrow(ae)) {
    es <- pmin(pmax(ae$emis/HI_CIT, 0), 1)                                  # absolute -> LENGTH
    L  <- (0.35 + 0.65*es) * (MARGIN_UM*0.75/MSI_PIXEL_UM)
    srel <- if (is.finite(GHI) && GHI > GLO) pmin(pmax((ae$emis - GLO)/(GHI - GLO), 0), 1) else rep(0.5, nrow(ae))
    f  <- 0.25 + 1.75*srel                                                  # global relative score -> THICKNESS (0.25x..2.0x)
    for (j in seq_len(nrow(ae))) {
      x1 <- ae$x0[j] + ae$nx[j]*L[j]; y1 <- ae$y0[j] + ae$ny[j]*L[j]
      arrows(ae$x0[j], ae$y0[j], x1, y1, length = 0.06, lwd = 4.2*f[j], col = "black")
      arrows(ae$x0[j], ae$y0[j], x1, y1, length = 0.05, lwd = 2.2*f[j], col = "white")
    }
  }
  scalebar_grid(o$xsc, o$ysc)
}

# ---- radial profile (ion vs signed distance) for a set of organoids --------
radial_df <- function(rows) {
  br <- seq(-40, 160, by = 10); mids <- (br[-1] + br[-length(br)])/2
  out <- list()
  for (r in seq_len(nrow(rows))) {
    sec <- prep_sec(rows$sid[r]); k <- rows$instance[r]; z <- sec$z
    g <- z$instance == k | z$instance_catch == k
    bin <- cut(z$signed_dist_um[g], br, labels = FALSE); m <- tapply(val_cit[z$gidx[g]], bin, mean, na.rm = TRUE)
    prof <- rep(NA_real_, length(mids)); prof[as.integer(names(m))] <- as.numeric(m)
    out[[r]] <- data.frame(sid = rows$sid[r], instance = k, apical_class = rows$apical_class[r],
                           dist = mids, cit = prof, stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}
