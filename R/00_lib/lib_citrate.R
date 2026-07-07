#!/usr/bin/env Rscript
# lib_citrate.R  -- shared citrate extraction (gate phase + downstream)
#
# THE citrate definition for the whole pipeline: integrate a +-CITRATE_WIN_PPM
# window around the standard-anchored mass CITRATE_ANCHOR_MZ (lib_paths.R), read
# PER PIXEL from the RAW centroid imzML/.ibd. This deliberately bypasses the
# processed 25-ppm peak grid, whose single ~191 feature merges citrate with a
# +16 ppm co-isobar shoulder and cannot represent a +-7 ppm citrate window.
#
# Houses the low-level imzML reader (parse_imzml/read_spots/mz_win, moved here
# from 02_CitrateStandard/00_config.R so the gate and downstream share one copy)
# plus the per-sample helpers used by 02_CitrateStandard/05_build_citrate_cache.R
# and every downstream citrate consumer.
#
# CAVEAT (centroided data): each pixel has ONE merged 191 centroid; a narrow
# window SELECTS citrate-dominant pixels, it does not integrate citrate area.
# Documented limit of this instrument (R~2400); see docs / project memory.

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_paths.R"))

WIN_PPM <- 10        # default half-window for the gate's QC reads (uM spectra etc.)
HIST_BY <- 1.0e-4    # high-res spectrum bin (Da) for the gate's histogram reads

mz_win <- function(mz, ppm = WIN_PPM) c(mz*(1 - ppm/1e6), mz*(1 + ppm/1e6))
ppm_of <- function(a, b) (a - b)/b * 1e6

# ---- low-level imzML/.ibd reader -------------------------------------------
parse_imzml <- function(imz) {
  x <- readChar(imz, file.info(imz)$size, useBytes=TRUE)
  pat <- paste0('ref="(mzArray|intensities)"/>\\s*<cvParam[^>]*external array length" value="([0-9]+)"',
                '[^>]*/>\\s*<cvParam[^>]*external encoded length" value="[0-9]+"[^>]*/>\\s*',
                '<cvParam[^>]*external offset" value="([0-9]+)"')
  m <- gregexpr(pat,x,perl=TRUE)[[1]]; st<-attr(m,"capture.start"); ln<-attr(m,"capture.length")
  kind<-substring(x,st[,1],st[,1]+ln[,1]-1); alen<-as.numeric(substring(x,st[,2],st[,2]+ln[,2]-1))
  off<-as.numeric(substring(x,st[,3],st[,3]+ln[,3]-1)); mzi<-which(kind=="mzArray"); ini<-which(kind=="intensities")
  gx<-regmatches(x,gregexpr('name="position x" value="[0-9]+"',x,perl=TRUE))[[1]]
  gy<-regmatches(x,gregexpr('name="position y" value="[0-9]+"',x,perl=TRUE))[[1]]
  list(ibd=sub("\\.imzML$",".ibd",imz),mz_off=off[mzi],mz_len=alen[mzi],in_off=off[ini],
       x=as.integer(sub('.*"([0-9]+)"$',"\\1",gx)), y=as.integer(sub('.*"([0-9]+)"$',"\\1",gy)))
}

# one-pass reader. Returns per pixel: window sums M [n x length(wins)], TIC [n];
# per window: intensity-weighted centroid accumulators cnum/cden; optional
# high-res histogram over hrange by HIST_BY. wins: list of c(lo,hi) m/z windows.
read_spots <- function(p, wins, hrange = NULL) {
  con <- file(p$ibd, "rb"); on.exit(close(con)); n <- length(p$mz_off)
  M <- matrix(0, n, length(wins)); tic <- numeric(n)
  cnum <- numeric(length(wins)); cden <- numeric(length(wins))
  use_h <- !is.null(hrange)
  if (use_h) { hbrk <- seq(hrange[1], hrange[2], by = HIST_BY)
               hctr <- (head(hbrk,-1)+tail(hbrk,-1))/2; hwl <- numeric(length(hctr)) }
  glo <- min(sapply(wins,`[`,1), if(use_h) hrange[1] else Inf)
  ghi <- max(sapply(wins,`[`,2), if(use_h) hrange[2] else -Inf)
  for (i in seq_len(n)) {
    if (p$mz_len[i] == 0) next
    seek(con,p$mz_off[i]); mz <- readBin(con,"double",p$mz_len[i],8,endian="little")
    seek(con,p$in_off[i]); ic <- readBin(con,"numeric",p$mz_len[i],4,endian="little")
    tic[i] <- sum(ic)
    if (!any(mz>=glo & mz<=ghi)) next
    if (use_h) { sel<-which(mz>=hrange[1]&mz<=hrange[2])
      if (length(sel)) { b<-findInterval(mz[sel],hbrk,rightmost.closed=TRUE); ok<-b>=1&b<=length(hwl)
        hwl[b[ok]] <- hwl[b[ok]] + ic[sel][ok] } }
    for (w in seq_along(wins)) {
      s <- mz>=wins[[w]][1] & mz<=wins[[w]][2]
      if (any(s)) { M[i,w]<-sum(ic[s]); cnum[w]<-cnum[w]+sum(ic[s]*mz[s]); cden[w]<-cden[w]+sum(ic[s]) }
    }
  }
  out <- list(x=p$x, y=p$y, M=M, tic=tic, cnum=cnum, cden=cden, n=n)
  if (use_h) { out$hctr<-hctr; out$hw<-hwl }
  out
}

read_file <- function(dir, file, wins, hrange = NULL)
  read_spots(parse_imzml(file.path(dir, file)), wins, hrange)

# ============================================================================
# Anchored citrate per-pixel extraction (the canonical pipeline citrate value)
# ============================================================================

# raw +-CITRATE_WIN_PPM citrate sum + TIC per pixel for one sample.
# sid = inventory sample_id; imz overrides the inventory imzml_path if given.
citrate_raw_pixels <- function(sid, imz = NULL) {
  if (is.null(imz)) imz <- inventory_row(sid)$imzml_path
  p   <- parse_imzml(imz)
  out <- read_spots(p, wins = list(mz_win(CITRATE_ANCHOR_MZ, CITRATE_WIN_PPM)))
  data.frame(x = out$x, y = out$y, cit_raw = out$M[,1], tic = out$tic)
}

CITRATE_CACHE <- function(sid) cache_in(sprintf("citrate_anchored_%s.rds", sid))

# load the precomputed per-sample citrate (built by 02_CitrateStandard/05);
# falls back to a live imzML read if the cache is absent.
load_citrate_cached <- function(sid) {
  f <- CITRATE_CACHE(sid)
  if (file.exists(f)) readRDS(f) else citrate_raw_pixels(sid)
}

# Build a full-length per-pixel citrate vector aligned to an MSE's pixel order.
# pd = as.data.frame(pixelData(mse)) with columns sample_id, x, y (gidx = row).
# normalize=TRUE -> cit_raw/tic (TIC-normalised, matching the old grid feature's
# scale, since the processed grid was TIC-normalised); FALSE -> raw sum.
# Left-joined on (sample_id, x, y); NA where a raw pixel is missing.
citrate_onto_pd <- function(pd, sids = NULL, normalize = TRUE) {
  if (is.null(sids)) sids <- unique(as.character(pd$sample_id))
  v <- rep(NA_real_, nrow(pd))
  for (sid in sids) {
    df  <- load_citrate_cached(sid)
    val <- if (normalize) df$cit_raw / ifelse(df$tic > 0, df$tic, NA_real_) else df$cit_raw
    sel <- which(as.character(pd$sample_id) == sid)
    key <- match(paste(pd$x[sel], pd$y[sel]), paste(df$x, df$y))
    v[sel] <- val[key]
  }
  v
}
