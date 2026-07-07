#!/usr/bin/env Rscript
# ============================================================================
# RETAINED ROLE: this script's lasting output is the .nd2 block-MEAN thumbnail
#   cache/register/nd2thumb_<slide>.rds (consumed by R/05_registration_refine/02_jpg_to_nd2_offset.R). Its flip/offset
#   estimate (jpg2nd2_<slide>.rds) is SUPERSEDED by R/05_registration_refine/02_jpg_to_nd2_offset.R's tissue-mask offset.
#   See docs/registration.md for the working pipeline order.
# ============================================================================
# 02_jpg_to_nd2.R - Phase R2: slide-JPG -> native .nd2 transform (per slide).
#   The _small.jpg is the whole-slide field downscaled ~2.895x and FLIPPED vs the
#   Nikon .nd2 (user). Scale is known from dims; FLIP + offset are found by phase-
#   correlating a block-MEAN .nd2 thumbnail (preserves thin organoids) against the
#   JPG darkness image under each of 4 flips. Refinement (R/05_registration_refine/01_refine_jpg.R) absorbs residual.
#
# Output: cache/register/nd2thumb_<slide>.rds  (block-mean thumb + factor)
#         cache/register/jpg2nd2_<slide>.rds   (scale, flip, offset, peak)
#         figures/registration/jpg2nd2_<slide>.pdf
# Usage : Rscript R/03_coarse_registration/02_jpg_to_nd2.R
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
suppressPackageStartupMessages({ library(RBioFormats); library(jpeg); library(png) })
REG_CACHE <- file.path(CACHE_DIR, "register"); REG_FIG <- file.path(FIG_DIR, "registration")

NV <- list(sl6A=list(nd2=file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"06102026_AO_0h_sl6A.nd2"),
                     jpg=file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"MSI/06102026_AO_0h_sl6A_small.jpg")),
           sl4A=list(nd2=file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"06102026_AO_20h_sl4A.nd2"),
                     jpg=file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"MSI/06102026_AO_20h_sl4A_small.jpg")))
F_THUMB <- 16L   # nd2 block-mean factor

# block-mean downsample of a huge .nd2 by reading horizontal row-chunks
nd2_block_mean <- function(path, F, chunk_blocks = 40L) {
  cm <- coreMetadata(read.metadata(path), series = 1); SX <- cm$sizeX; SY <- cm$sizeY
  Wn <- SX %/% F; Hn <- SY %/% F
  out <- matrix(0, Hn, Wn); rows_per <- F * chunk_blocks
  yb <- 1L
  while (yb <= Hn) {
    ye <- min(yb + chunk_blocks - 1L, Hn)
    y0 <- (yb-1L)*F + 1L; y1 <- ye*F
    img <- read.image(path, series = 1, subset = list(x = 1:(Wn*F), y = y0:y1), normalize = TRUE)
    a <- as.array(img); m <- if (length(dim(a))==2) a else a[,,1]   # [x, y] (W x H-chunk)
    m <- t(m)                                                       # [y, x]
    nb <- ye - yb + 1L
    # block-mean nb*F x Wn*F  ->  nb x Wn
    rid <- rep(1:nb, each = F); cid <- rep(1:Wn, each = F)
    bm <- rowsum(t(rowsum(m, rid)), cid)                            # sum over blocks
    out[yb:ye, ] <- t(bm) / (F*F)
    yb <- ye + 1L
  }
  list(m = out, F = F, SX = SX, SY = SY)
}

# phase correlation: returns best integer shift (dy,dx) of b to match a + peak
phase_corr <- function(a, b) {
  H <- min(nrow(a), nrow(b)); W <- min(ncol(a), ncol(b))
  a <- a[1:H, 1:W]; b <- b[1:H, 1:W]
  a <- a - mean(a); b <- b - mean(b)
  Fa <- fft(a); Fb <- fft(b); R <- Fa * Conj(Fb); R <- R / (Mod(R) + 1e-9)
  r <- Re(fft(R, inverse = TRUE)) / length(R)
  pk <- which.max(r); peak <- max(r)
  dy <- (pk - 1) %% H; dx <- (pk - 1) %/% H
  if (dy > H/2) dy <- dy - H; if (dx > W/2) dx <- dx - W
  list(dy = dy, dx = dx, peak = peak)
}

ds_mean <- function(m, f) {   # block-mean a matrix by integer factor f
  H <- (nrow(m) %/% f) * f; W <- (ncol(m) %/% f) * f; m <- m[1:H, 1:W]
  rid <- rep(1:(H/f), each = f); cid <- rep(1:(W/f), each = f)
  t(rowsum(t(rowsum(m, rid)), cid)) / (f*f)
}

for (sl in names(NV)) {
  cat(sprintf("[30c] %s: block-mean nd2 (factor %d)...\n", sl, F_THUMB))
  nb <- nd2_block_mean(NV[[sl]]$nd2, F_THUMB)
  saveRDS(nb, file.path(REG_CACHE, sprintf("nd2thumb_%s.rds", sl)))
  Dn <- 1 - nb$m / max(nb$m)                       # nd2 darkness thumb (Hn x Wn)
  j <- jpeg::readJPEG(NV[[sl]]$jpg); jg <- if (length(dim(j))==3) (j[,,1]+j[,,2]+j[,,3])/3 else j
  scale_jn <- nb$SX / ncol(jg)                     # ~2.895 (nd2 px per jpg px)
  fj <- F_THUMB / scale_jn                         # jpg downsample to nd2-thumb scale
  Dj <- 1 - ds_mean(jg, max(1, round(fj)))
  Dj <- Dj - min(Dj); Dn <- Dn - min(Dn)
  flips <- list(none=function(x) x, V=function(x) x[nrow(x):1,], H=function(x) x[,ncol(x):1],
                VH=function(x) x[nrow(x):1, ncol(x):1])
  best <- NULL
  for (fn in names(flips)) {
    Djf <- flips[[fn]](Dj); pc <- phase_corr(Dn, Djf)
    cat(sprintf("[30c]   flip %-4s peak=%.4f shift(dy=%d,dx=%d)\n", fn, pc$peak, pc$dy, pc$dx))
    if (is.null(best) || pc$peak > best$peak) best <- c(pc, flip = fn)
  }
  # Build jpg->nd2 mapping at thumb scale, then to full nd2:
  # thumb: nd2thumb_pt = flip(jpg_thumb_pt) + (dx,dy); full nd2 = thumb * F_THUMB
  Hn <- nrow(Dn); Wn <- ncol(Dn); Hj <- nrow(Dj); Wj <- ncol(Dj)
  saveRDS(list(slide = sl, scale_jpg_nd2 = scale_jn, F_thumb = F_THUMB, flip = best$flip,
               dy_thumb = best$dy, dx_thumb = best$dx, peak = best$peak,
               nd2_W = nb$SX, nd2_H = nb$SY, jpg_W = ncol(jg), jpg_H = nrow(jg),
               thumb_Wn = Wn, thumb_Hn = Hn, jpgthumb_W = Wj, jpgthumb_H = Hj, jpg_ds = max(1, round(fj))),
          file.path(REG_CACHE, sprintf("jpg2nd2_%s.rds", sl)))
  cat(sprintf("[30c] %s -> flip=%s peak=%.4f scale=%.4f\n", sl, best$flip, best$peak, scale_jn))
  # QC: nd2 thumb with flipped jpg darkness contours overlaid at found shift
  pdf(file.path(REG_FIG, sprintf("jpg2nd2_%s.pdf", sl)), width = 12, height = 7)
  par(mfrow = c(1,2), mar = c(2,2,3,1))
  image(t(Dn)[, nrow(Dn):1], col = grey.colors(64), main = sprintf("%s nd2 block-mean (darkness)", sl), axes = FALSE)
  Djf <- flips[[best$flip]](Dj)
  image(t(Djf)[, nrow(Djf):1], col = grey.colors(64), main = sprintf("jpg darkness (flip=%s)", best$flip), axes = FALSE)
  dev.off()
}
cat("[30c] DONE\n")
