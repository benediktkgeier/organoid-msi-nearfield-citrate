#!/usr/bin/env Rscript
# ============================================================================
# lib_register_if.R - helpers for registering IF (B-section) .nd2 images onto
#   the MSI brightfield / MSI grid. Builds on R/lib_register.R (must be sourced
#   first: fit_affine/apply_affine/compose_affine/phase_corr/ds_mean/rs_transform).
#
# Chain (see docs/if_registration.md / plan):
#   hi-res .nd2 (0.362 um/px, 4ch) --B_hr_ov--> overview .nd2 (1.833 um/px, DAPI)
#   overview .nd2 --B_ov_bf--> MSI BF .nd2 (1.833 um/px) --inv(B_msi_nd2)--> MSI grid
# ============================================================================

suppressPackageStartupMessages({ library(RBioFormats) })

# ---- low-memory block-MEAN of ONE channel of a (possibly huge) .nd2 --------
# Reads horizontal row-chunks; selects channel `ch` via read.image subset.
# Returns list(m=[Hn x Wn] matrix, F, SX, SY, ch, nC).
nd2_channel_block_mean <- function(path, ch = 1L, F = 8L, chunk_blocks = 60L, normalize = TRUE) {
  cm <- coreMetadata(read.metadata(path), series = 1)
  SX <- cm$sizeX; SY <- cm$sizeY; nC <- cm$sizeC
  ch <- max(1L, min(as.integer(ch), nC))
  F  <- max(1L, as.integer(F))
  Wn <- SX %/% F; Hn <- SY %/% F
  out <- matrix(0, Hn, Wn); yb <- 1L
  while (yb <= Hn) {
    ye <- min(yb + chunk_blocks - 1L, Hn)
    y0 <- (yb - 1L) * F + 1L; y1 <- ye * F
    sub <- if (nC > 1L) list(x = 1:(Wn*F), y = y0:y1, c = ch) else list(x = 1:(Wn*F), y = y0:y1)
    img <- read.image(path, series = 1, subset = sub, normalize = normalize)
    a <- as.array(img)
    m <- if (length(dim(a)) == 2) a else a[, , 1]    # [x, y-chunk]
    m <- t(m)                                        # [y, x]
    nb <- ye - yb + 1L
    rid <- rep(1:nb, each = F); cid <- rep(1:Wn, each = F)
    bm <- rowsum(t(rowsum(m, rid)), cid)             # sum over F x F blocks
    out[yb:ye, ] <- t(bm) / (F * F)
    yb <- ye + 1L
  }
  list(m = out, F = F, SX = SX, SY = SY, ch = ch, nC = nC)
}

# block-MIN downsample of one .nd2 channel (preserves thin DARK tissue edges that
# block-mean averages away - good for faint post-MALDI brightfield organoids).
nd2_block_min <- function(path, F = 6L, chunk_blocks = 60L, ch = 1L, normalize = TRUE) {
  cm <- coreMetadata(read.metadata(path), series = 1)
  SX <- cm$sizeX; SY <- cm$sizeY; nC <- cm$sizeC; F <- max(1L, as.integer(F))
  Wn <- SX %/% F; Hn <- SY %/% F; out <- matrix(0, Hn, Wn); yb <- 1L
  while (yb <= Hn) {
    ye <- min(yb + chunk_blocks - 1L, Hn); y0 <- (yb-1L)*F + 1L; y1 <- ye*F
    sub <- if (nC > 1L) list(x = 1:(Wn*F), y = y0:y1, c = ch) else list(x = 1:(Wn*F), y = y0:y1)
    img <- read.image(path, series = 1, subset = sub, normalize = normalize)
    a <- as.array(img); m <- if (length(dim(a)) == 2) a else a[, , 1]; m <- t(m)   # [y, x]
    nb <- ye - yb + 1L
    ar <- array(m, c(F, nb, ncol(m))); rmin <- ar[1, , ]; for (k in 2:F) rmin <- pmin(rmin, ar[k, , ])
    tr <- t(rmin); ac <- array(tr, c(F, Wn, nb)); cmin <- ac[1, , ]; for (k in 2:F) cmin <- pmin(cmin, ac[k, , ])
    out[yb:ye, ] <- t(cmin); yb <- ye + 1L
  }
  list(m = out, F = F, SX = SX, SY = SY)
}

# block-MEAN ALL channels of an .nd2 in a single row-chunked pass -> [Hn,Wn,nC]
nd2_block_mean_allch <- function(path, F = 5L, chunk_blocks = 60L, normalize = TRUE) {
  cm <- coreMetadata(read.metadata(path), series = 1)
  SX <- cm$sizeX; SY <- cm$sizeY; nC <- cm$sizeC; F <- max(1L, as.integer(F))
  Wn <- SX %/% F; Hn <- SY %/% F; out <- array(0, c(Hn, Wn, nC)); yb <- 1L
  while (yb <= Hn) {
    ye <- min(yb + chunk_blocks - 1L, Hn); y0 <- (yb-1L)*F + 1L; y1 <- ye*F
    img <- read.image(path, series = 1, subset = list(x = 1:(Wn*F), y = y0:y1), normalize = normalize)
    a <- as.array(img); if (length(dim(a)) == 2) a <- array(a, c(dim(a), 1))
    nb <- ye - yb + 1L; rid <- rep(1:nb, each = F); cid <- rep(1:Wn, each = F)
    for (k in 1:nC) { m <- t(a[, , k]); bm <- rowsum(t(rowsum(m, rid)), cid); out[yb:ye, , k] <- t(bm)/(F*F) }
    yb <- ye + 1L
  }
  out
}

# normalize a matrix to [0,1] robustly (p1..p99), NA-safe
norm01 <- function(m, lo = 0.01, hi = 0.995) {
  q <- quantile(m[is.finite(m)], c(lo, hi), na.rm = TRUE)
  if (!is.finite(q[1]) || !is.finite(q[2]) || q[2] <= q[1]) { r <- range(m[is.finite(m)]); q <- if (diff(r) > 0) r else c(0, 1) }
  pmin(pmax((m - q[1]) / (q[2] - q[1]), 0), 1)
}

# ---- tissue masks ----------------------------------------------------------
# DAPI thumb: nuclei/tissue = signal above background. Otsu-ish quantile thr.
dapi_mask <- function(m, q = 0.80) {
  v <- norm01(m); thr <- quantile(v[is.finite(v)], q, na.rm = TRUE)
  (v >= thr) & is.finite(v)
}
# BF thumb (brightfield, tissue = DARK): invert then threshold.
bf_mask <- function(m, q = 0.55) {
  v <- norm01(m); d <- 1 - v; thr <- quantile(d[is.finite(d)], q, na.rm = TRUE)
  (d >= thr) & is.finite(d)
}

# ---- mask overlap (Intersection-over-Union) on common HxW -----------------
mask_iou <- function(A, B) {
  H <- min(nrow(A), nrow(B)); W <- min(ncol(A), ncol(B))
  A <- A[1:H, 1:W]; B <- B[1:H, 1:W]
  i <- sum(A & B); u <- sum(A | B); if (u == 0) 0 else i / u
}

# flip a logical/numeric matrix: "none","H" (left-right),"V" (up-down),"HV"
flip_mat <- function(m, f) {
  if (f == "H")  return(m[, ncol(m):1, drop = FALSE])
  if (f == "V")  return(m[nrow(m):1, , drop = FALSE])
  if (f == "HV") return(m[nrow(m):1, ncol(m):1, drop = FALSE])
  m
}

# rotate a matrix by k*90 deg (k=0..3), exact (no interpolation)
rot90_mat <- function(m, k) {
  k <- k %% 4L
  if (k == 0L) return(m)
  if (k == 1L) return(t(m[nrow(m):1, , drop = FALSE]))            # 90 CCW
  if (k == 2L) return(m[nrow(m):1, ncol(m):1, drop = FALSE])      # 180
  t(m[, ncol(m):1, drop = FALSE])                                 # 270 CCW
}

# ---- coarse flip x rot90 x translation search (mask IoU) -------------------
# Move maskB (moving) onto maskA (fixed). Same pixel scale assumed.
# Returns best flip, rot90 k, integer (dy,dx) shift of B, and IoU.
search_flip_rot_trans <- function(maskA, maskB, flips = c("none","H","V","HV"),
                                  ks = 0:3, refine = 6L) {
  best <- list(flip = "none", k = 0L, dy = 0, dx = 0, iou = -1)
  Af <- maskA * 1.0
  for (fl in flips) for (k in ks) {
    Bf <- rot90_mat(flip_mat(maskB, fl), k) * 1.0
    pc <- phase_corr(Af, Bf)                                      # coarse shift via phase corr
    # build shifted B and score IoU; small local refine around the pc peak
    for (ddy in seq(-refine, refine, by = 2L)) for (ddx in seq(-refine, refine, by = 2L)) {
      dy <- pc$dy + ddy; dx <- pc$dx + ddx
      Bs <- shift_mask(Bf > 0.5, dy, dx, dim(Af))
      io <- mask_iou(maskA, Bs)
      if (io > best$iou) best <- list(flip = fl, k = k, dy = dy, dx = dx, iou = io)
    }
  }
  best
}

# shift a logical mask by (dy,dx) into a canvas of dims `dd` (rows,cols)
shift_mask <- function(M, dy, dx, dd) {
  H <- dd[1]; W <- dd[2]
  out <- matrix(FALSE, H, W)
  idx <- which(M, arr.ind = TRUE)
  r <- idx[, 1] + dy; c <- idx[, 2] + dx
  ok <- r >= 1 & r <= H & c >= 1 & c <= W
  if (any(ok)) out[cbind(r[ok], c[ok])] <- TRUE
  out
}

# ---- pick the DAPI channel of a hi-res .nd2 by max phase-corr to overview ---
# Given a coarse overview-DAPI thumb window (already located), test each channel
# downsampled to the same scale; return the channel with the strongest peak.
# Fallback: default `default_ch` if all peaks are weak.
auto_dapi_channel <- function(hr_path, ov_window, F_hr, default_ch = 1L, nC = 4L) {
  best <- list(ch = default_ch, peak = -1)
  for (ch in 1:nC) {
    th <- nd2_channel_block_mean(hr_path, ch = ch, F = F_hr)$m
    pc <- phase_corr(norm01(ov_window), norm01(th))
    if (pc$peak > best$peak) best <- list(ch = ch, peak = pc$peak)
  }
  best
}

# draw a binary mask outline as points (for QC overlays)
mask_outline_pts <- function(M, step = 1L) {
  idx <- which(M, arr.ind = TRUE)
  if (step > 1L) idx <- idx[seq(1, nrow(idx), by = step), , drop = FALSE]
  idx  # columns: row(=y), col(=x)
}

# ---- similarity transform on continuous (x=col,y=row) points ---------------
# par = list(flip, theta(rad), s, dx, dy); src/target thumbs share um/px scale.
# Flip is taken w.r.t. the source-thumb dims `dd`=c(H,W); rotation about its center.
apply_sim <- function(xy, par, dd) {
  x <- xy[, 1]; y <- xy[, 2]; H <- dd[1]; W <- dd[2]
  if (!is.null(par$flip)) {
    if (par$flip %in% c("H", "HV")) x <- W + 1 - x
    if (par$flip %in% c("V", "HV")) y <- H + 1 - y
  }
  cx <- (W + 1) / 2; cy <- (H + 1) / 2
  ct <- cos(par$theta); st <- sin(par$theta); s <- par$s
  X <- x - cx; Y <- y - cy
  cbind(cx + s * (ct * X - st * Y) + par$dx, cy + s * (st * X + ct * Y) + par$dy)
}

# rasterize continuous (x,y) foreground points onto a logical canvas of dims dd
rasterize_pts <- function(xy, dd) {
  H <- dd[1]; W <- dd[2]; out <- matrix(FALSE, H, W)
  c <- round(xy[, 1]); r <- round(xy[, 2])
  ok <- r >= 1 & r <= H & c >= 1 & c <= W
  if (any(ok)) out[cbind(r[ok], c[ok])] <- TRUE
  out
}

# normalized cross-correlation of two equal-size matrices (mean-subtracted)
ncc_score <- function(a, b) {
  a <- a - mean(a); b <- b - mean(b)
  d <- sqrt(sum(a * a) * sum(b * b)); if (d <= 0) return(0)
  sum(a * b) / d
}

# windowed template match: slide `tmpl` over `img` for row/col offsets in
# `rseq`/`cseq` (offset of tmpl's top-left in img). Returns best offset + ncc.
locate_template <- function(tmpl, img, rseq, cseq) {
  th <- nrow(tmpl); tw <- ncol(tmpl); H <- nrow(img); W <- ncol(img)
  best <- list(r0 = NA, c0 = NA, ncc = -2)
  for (r0 in rseq) for (c0 in cseq) {
    if (r0 < 1 || c0 < 1 || r0 + th - 1 > H || c0 + tw - 1 > W) next
    sub <- img[r0:(r0 + th - 1), c0:(c0 + tw - 1)]
    sc <- ncc_score(tmpl, sub)
    if (sc > best$ncc) best <- list(r0 = r0, c0 = c0, ncc = sc)
  }
  best
}

# rotate a matrix by theta (rad) about center, bilinear (needs EBImage)
rotate_mat <- function(m, theta) {
  if (abs(theta) < 1e-6) return(m)
  EBImage::rotate(m, angle = theta * 180 / pi, output.dim = dim(m), bg.col = mean(m))
}

# local-texture (std) map - emphasizes punctate nuclei, suppresses smooth bubbles
texture_map <- function(m, k = 2L) {
  ker <- matrix(1, 2*k+1, 2*k+1) / ((2*k+1)^2)
  v <- norm01(m)
  mu  <- EBImage::filter2(v, ker, boundary = "replicate")
  mu2 <- EBImage::filter2(v*v, ker, boundary = "replicate")
  sqrt(pmax(mu2 - mu*mu, 0))
}

# Locate a hi-res DAPI section within its overview by FFT-NCC on texture maps,
# at the coarse (F16/F20) scale, with a small rotation search. Returns the
# similarity affine B_hr_ov (hi-res native px -> overview native px), the matched
# center in overview-native px, the NCC peak, and the rotation angle (deg).
#   ov : ovthumb list (m16,F16,SX,SY);  h : hrthumb list (m20,F20,SX,SY)
locate_hires_in_overview <- function(ov, h, angs = seq(-8, 8, 2)) {
  imupx <- OV_UMPX * ov$F16                       # overview coarse um/px (~29.3)
  sc    <- (HR_UMPX * h$F20) / imupx              # resize hi-res(7.24)->29.3
  Traw  <- EBImage::resize(texture_map(h$m20), round(nrow(h$m20)*sc), round(ncol(h$m20)*sc))
  Iimg  <- as.matrix(texture_map(ov$m16))
  best  <- list(ncc = -2)
  for (ang in angs) {
    Tt  <- if (ang == 0) Traw else EBImage::rotate(Traw, ang, output.dim = dim(Traw), bg.col = 0)
    res <- fast_ncc(Iimg, as.matrix(Tt))
    if (res$ncc > best$ncc) best <- c(res, list(ang = ang, th = nrow(Tt), tw = ncol(Tt)))
  }
  # matched center in overview-native px
  cc_coarse <- best$c0 + best$tw/2 - 0.5; rr_coarse <- best$r0 + best$th/2 - 0.5
  center_ov <- c((cc_coarse - 0.5) * ov$F16 + 0.5, (rr_coarse - 0.5) * ov$F16 + 0.5)
  # build B_hr_ov (similarity) from control points
  s   <- HR_UMPX / OV_UMPX                        # overview-native px per hi-res-native px (~0.1975)
  th  <- best$ang * pi / 180
  cxh <- (h$SX + 1)/2; cyh <- (h$SY + 1)/2
  gx  <- rep(seq(1, h$SX, length.out = 6), times = 6)
  gy  <- rep(seq(1, h$SY, length.out = 6), each  = 6)
  X <- gx - cxh; Y <- gy - cyh
  ovx <- center_ov[1] + s*( cos(th)*X - sin(th)*Y)
  ovy <- center_ov[2] + s*( sin(th)*X + cos(th)*Y)
  B_hr_ov <- fit_affine(cbind(gx, gy), cbind(ovx, ovy))
  list(B_hr_ov = B_hr_ov, center_ov = center_ov, ncc = best$ncc, ang = best$ang,
       box = c(r0 = best$r0, c0 = best$c0, th = best$th, tw = best$tw), F16 = ov$F16)
}
integral_img <- function(m) {
  S <- apply(m, 2, cumsum); S <- t(apply(S, 1, cumsum))
  rbind(0, cbind(0, S))
}
# window sums over h x w boxes from an integral image; result [H-h+1, W-w+1]
window_sum <- function(II, h, w) {
  H <- nrow(II) - 1L; W <- ncol(II) - 1L
  r2 <- (h):(H); c2 <- (w):(W); r1 <- r2 - h; c1 <- c2 - w
  II[r2 + 1L, c2 + 1L, drop = FALSE] - II[r1 + 1L, c2 + 1L, drop = FALSE] -
    II[r2 + 1L, c1 + 1L, drop = FALSE] + II[r1 + 1L, c1 + 1L, drop = FALSE]
}

# Fast normalized cross-correlation of template T over image I (FFT, Lewis 1995).
# Returns NCC map [H-h+1, W-w+1] (value at top-left window origin) + best location.
fast_ncc <- function(I, T) {
  H <- nrow(I); W <- ncol(I); h <- nrow(T); w <- ncol(T)
  Tz <- T - mean(T); tnorm <- sqrt(sum(Tz * Tz)); if (tnorm <= 0) tnorm <- 1
  Tp <- matrix(0, H, W); Tp[1:h, 1:w] <- Tz
  num <- Re(fft(fft(I) * Conj(fft(Tp)), inverse = TRUE)) / (H * W)   # cross-corr; [r,c]=window origin
  II  <- integral_img(I); II2 <- integral_img(I * I)
  ws  <- window_sum(II, h, w); ws2 <- window_sum(II2, h, w)
  n   <- h * w
  denom <- sqrt(pmax(ws2 - ws * ws / n, 0)) * tnorm
  nr <- H - h + 1L; nc <- W - w + 1L
  map <- num[1:nr, 1:nc] / ifelse(denom > 0, denom, NA)
  pk <- which.max(replace(map, is.na(map), -Inf))
  r0 <- ((pk - 1) %% nr) + 1L; c0 <- ((pk - 1) %/% nr) + 1L
  list(map = map, r0 = r0, c0 = c0, ncc = map[r0, c0])
}
