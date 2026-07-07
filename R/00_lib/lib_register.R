#!/usr/bin/env Rscript
# ============================================================================
# lib_register.R - helpers for MSI <-> optical (brightfield) registration.
# Chain: MSI pixel -> stage (poslog) -> flexImage -> original _small.jpg
#        (teach points), then jpg -> .nd2 (scale+flip) handled in R/03_coarse_registration/02_jpg_to_nd2.R.
# ============================================================================

# ---- parse a Bruker flexImaging .mis (imaging sequence) --------------------
parse_mis <- function(path) {
  L <- readLines(path, warn = FALSE)
  grab_pairs <- function(tag) {
    ln <- grep(sprintf("<%s>", tag), L, value = TRUE, fixed = TRUE)
    out <- lapply(ln, function(s) {
      s <- sub(sprintf(".*<%s>", tag), "", s); s <- sub(sprintf("</%s>.*", tag), "", s)
      ab <- strsplit(s, ";")[[1]]
      a <- as.numeric(strsplit(ab[1], ",")[[1]]); b <- as.numeric(strsplit(ab[2], ",")[[1]])
      c(a, b)
    })
    if (!length(out)) return(matrix(numeric(0), 0, 4))
    do.call(rbind, out)
  }
  teach <- grab_pairs("TeachPoint")              # flexImgX,flexImgY, stageX,stageY
  orig  <- grab_pairs("OriginalImageTeachPoint") # flexImgX,flexImgY, origX,origY
  # Area polygon points (flexImg coords)
  pl <- grep("<Point>", L, value = TRUE, fixed = TRUE)
  area <- if (length(pl)) do.call(rbind, lapply(pl, function(s) {
    s <- sub(".*<Point>", "", s); s <- sub("</Point>.*", "", s); as.numeric(strsplit(s, ",")[[1]]) })) else matrix(numeric(0),0,2)
  oimg <- sub(".*<OriginalImage>", "", grep("<OriginalImage>", L, value = TRUE, fixed = TRUE)[1]); oimg <- sub("</OriginalImage>.*", "", oimg)
  list(teach = teach, orig = orig, area = area, orig_image = oimg)
}

# ---- parse poslog -> raster index + stage coords ---------------------------
parse_poslog <- function(path) {
  L <- readLines(path, warn = FALSE)
  L <- L[!grepl("^#", L)]
  m <- regmatches(L, regexpr("R\\d+X\\d+Y\\d+", L))
  keep <- grepl("R\\d+X\\d+Y\\d+", L)
  L <- L[keep]
  rx <- as.integer(sub(".*X(\\d+)Y\\d+.*", "\\1", regmatches(L, regexpr("X\\d+Y\\d+", L))))
  ry <- as.integer(sub(".*Y(\\d+).*", "\\1", regmatches(L, regexpr("Y\\d+", L))))
  # tokens after the Pos id: first two numbers = commanded stage X,Y (regular 10um grid)
  rest <- sub(".*R\\d+X\\d+Y\\d+\\s+", "", L)
  nums <- lapply(strsplit(trimws(rest), "\\s+"), as.numeric)
  sx <- vapply(nums, `[`, numeric(1), 1); sy <- vapply(nums, `[`, numeric(1), 2)
  data.frame(rasterX = rx, rasterY = ry, stageX = sx, stageY = sy)
}

# ---- least-squares affine: dst(Nx2) ~ src(Nx2); returns 3x2 matrix B -------
fit_affine <- function(src, dst) {
  Xa <- cbind(src, 1)
  B <- solve(t(Xa) %*% Xa, t(Xa) %*% dst)   # 3x2
  B
}
apply_affine <- function(B, pts) { cbind(pts, 1) %*% B }
affine_rmse <- function(B, src, dst) sqrt(mean(rowSums((apply_affine(B, src) - dst)^2)))

# compose two affines: first apply B1 (a->b), then B2 (b->c) -> a->c (3x2)
compose_affine <- function(B1, B2) {
  M1 <- rbind(t(B1), c(0,0,1)); M2 <- rbind(t(B2), c(0,0,1))
  M <- M2 %*% M1; t(M[1:2, ])
}

# ---- image helpers (registration) -----------------------------------------
# block-mean downsample a matrix by integer factor f
ds_mean <- function(m, f) {
  f <- max(1L, as.integer(f)); if (f == 1L) return(m)
  H <- (nrow(m) %/% f) * f; W <- (ncol(m) %/% f) * f; m <- m[1:H, 1:W, drop = FALSE]
  rid <- rep(1:(H/f), each = f); cid <- rep(1:(W/f), each = f)
  t(rowsum(t(rowsum(m, rid)), cid)) / (f*f)
}

# phase correlation: integer shift (dy,dx) to move b onto a, + normalized peak
phase_corr <- function(a, b) {
  H <- min(nrow(a), nrow(b)); W <- min(ncol(a), ncol(b))
  a <- a[1:H, 1:W]; b <- b[1:H, 1:W]; a <- a - mean(a); b <- b - mean(b)
  if (sd(a) == 0 || sd(b) == 0) return(list(dy = 0, dx = 0, peak = 0))
  Fa <- fft(a); Fb <- fft(b); R <- Fa * Conj(Fb); R <- R / (Mod(R) + 1e-9)
  r <- Re(fft(R, inverse = TRUE)) / length(R)
  pk <- which.max(r); peak <- max(r)
  dy <- (pk - 1) %% H; dx <- (pk - 1) %/% H
  if (dy > H/2) dy <- dy - H; if (dx > W/2) dx <- dx - W
  list(dy = dy, dx = dx, peak = peak)
}

# rigid+scale transform of points about a center: scale s, rotation theta(rad)
rs_transform <- function(pts, cx, cy, s, theta, tx = 0, ty = 0) {
  ct <- cos(theta); st <- sin(theta)
  X <- pts[,1] - cx; Y <- pts[,2] - cy
  cbind(cx + s*(ct*X - st*Y) + tx, cy + s*(st*X + ct*Y) + ty)
}
