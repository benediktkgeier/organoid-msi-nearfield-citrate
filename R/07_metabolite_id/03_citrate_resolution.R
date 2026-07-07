#!/usr/bin/env Rscript
# 03_citrate_resolution.R   (revamped)
# Shows that the timsTOF fleX cannot cleanly resolve citrate [M-H]- from the
# co-isobaric ion(s) at nominal m/z 191, and that the SCiLS feature the metabolite
# report labels "Citrate" is centred ~+8 ppm above the true citrate mass.
#
# Both pages are built from FOUR sections chosen for genuine 191 signal:
#   AO_0h_sl6A_sec2a, AO_0h_sl6A_sec5a, AO_20h_sl4A_sec3a, AO_20h_sl4A_sec3b
#
# Page 1  Publication figure: intensity-weighted reconstructed spectrum around
#         m/z 191 + per-pixel centroid distribution (signal-bearing pixels only),
#         vs theoretical citrate (191.0197) and the report grid ion (191.0217);
#         annotates the blended-peak envelope and the unresolved 10-ppm gap.
# Page 2  Ion images: the full "citrate" feature, then its m/z sub-bands and the
#         dominant +138 ppm isobar, across the 4 sections, with the spatial
#         correlation of each band vs the citrate-mass band.
#
# Reads SCiLS centroids straight from the .ibd (processed/centroid imzML) -> no
# rebinning. Out: figures/metabolites/citrate_resolution_report.pdf (standalone).
# Usage: Rscript R/07_metabolite_id/03_citrate_resolution.R

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_paths.R"))
suppressPackageStartupMessages(library(viridisLite))

CIT  <- 191.019726          # citrate [M-H]- C6H7O7  (theoretical)
GRID <- 191.021704          # feature m/z the report matched & labelled citrate
OUT_PDF <- file.path(FIG_DIR, "metabolites", "citrate_resolution_report.pdf")
IMZDIR  <- file.path(CACHE_DIR, "imzml")
ppm <- function(a, b = CIT) (a - b) / b * 1e6
log_msg <- function(...) message(sprintf("[82] %s", sprintf(...)))

# sections with proper citrate signal (display order: 0h, 0h, 20h, 20h)
SECT <- list(
  list(sid="AO_0h_sl6A_sec2a",  grp="0h",  imz="06102026_ao_0h_sl6a_sec2a.imzML"),
  list(sid="AO_0h_sl6A_sec5a",  grp="0h",  imz="06102026_ao_0h_sl6a_sec5a.imzML"),
  list(sid="AO_20h_sl4A_sec3a", grp="20h", imz="06102026_ao_20h_sl4a_3a.imzML"),
  list(sid="AO_20h_sl4A_sec3b", grp="20h", imz="06102026_ao_20h_sl4a_3b.imzML"))

# m/z bands. Band 1 = the whole "citrate" feature the report renders; bands 2-4
# are disjoint thirds of that blend; band 5 = the dominant non-citrate isobar.
BANDS <- list(
  list(name="report \"citrate\" feature", c=191.0217, lo=191.0188, hi=191.0246, ref=FALSE, full=TRUE),
  list(name="citrate region",            c=191.0197, lo=191.0188, hi=191.0206, ref=TRUE,  full=FALSE),
  list(name="report apex",               c=191.0217, lo=191.0208, hi=191.0226, ref=FALSE, full=FALSE),
  list(name="high shoulder",             c=191.0237, lo=191.0228, hi=191.0246, ref=FALSE, full=FALSE),
  list(name="dominant isobar",           c=191.0461, lo=191.0440, hi=191.0482, ref=FALSE, full=FALSE))
wins   <- lapply(BANDS, function(b) c(b$lo, b$hi))
REFB   <- which(vapply(BANDS, `[[`, logical(1), "ref"))

WLO <- 190.95; WHI <- 191.10                    # context window for spectrum
brk <- seq(WLO, WHI, by = 0.00010); ctr <- (head(brk,-1)+tail(brk,-1))/2
CL_LO <- 191.013; CL_HI <- 191.027              # the citrate "blend" cluster

# ---------- imzML (processed/centroid) low-level reader ----------------------
parse_imzml <- function(imz) {
  x <- readChar(imz, file.info(imz)$size, useBytes = TRUE)
  pat <- paste0('ref="(mzArray|intensities)"/>\\s*<cvParam[^>]*external array length" ',
                'value="([0-9]+)"[^>]*/>\\s*<cvParam[^>]*external encoded length" ',
                'value="[0-9]+"[^>]*/>\\s*<cvParam[^>]*external offset" value="([0-9]+)"')
  m  <- gregexpr(pat, x, perl = TRUE)[[1]]
  st <- attr(m, "capture.start"); ln <- attr(m, "capture.length")
  kind <- substring(x, st[,1], st[,1]+ln[,1]-1)
  alen <- as.numeric(substring(x, st[,2], st[,2]+ln[,2]-1))
  off  <- as.numeric(substring(x, st[,3], st[,3]+ln[,3]-1))
  mzi <- which(kind == "mzArray"); ini <- which(kind == "intensities")
  gx <- regmatches(x, gregexpr('name="position x" value="[0-9]+"', x, perl=TRUE))[[1]]
  gy <- regmatches(x, gregexpr('name="position y" value="[0-9]+"', x, perl=TRUE))[[1]]
  list(ibd = sub("\\.imzML$", ".ibd", imz),
       mz_off = off[mzi], mz_len = alen[mzi], in_off = off[ini],
       x = as.integer(sub('.*"([0-9]+)"$', "\\1", gx)),
       y = as.integer(sub('.*"([0-9]+)"$', "\\1", gy)))
}

# one pass: band sums (pixel x band) + spectrum histogram + per-pixel strongest
# centroid in the cluster (m/z + intensity).
read_dataset <- function(p) {
  con <- file(p$ibd, "rb"); on.exit(close(con))
  n <- length(p$mz_off)
  M <- matrix(0, n, length(wins)); hwl <- numeric(length(ctr))
  pkm <- rep(NA_real_, n); pki <- numeric(n)
  for (i in seq_len(n)) {
    seek(con, p$mz_off[i]); mz <- readBin(con, "double",  p$mz_len[i], 8, endian="little")
    sel <- which(mz >= WLO & mz <= WHI); if (!length(sel)) next
    seek(con, p$in_off[i]); ic <- readBin(con, "numeric", p$mz_len[i], 4, endian="little")
    mw <- mz[sel]; iw <- ic[sel]
    b <- findInterval(mw, brk, rightmost.closed=TRUE); ok <- b>=1 & b<=length(hwl)
    hwl[b[ok]] <- hwl[b[ok]] + iw[ok]
    for (w in seq_along(wins)) { s <- mw>=wins[[w]][1] & mw<=wins[[w]][2]; if (any(s)) M[i,w] <- sum(iw[s]) }
    inc <- which(mw >= CL_LO & mw <= CL_HI)
    if (length(inc)) { j <- inc[which.max(iw[inc])]; pkm[i] <- mw[j]; pki[i] <- iw[j] }
  }
  list(x=p$x, y=p$y, M=M, hw=hwl, pk_mz=pkm, pk_in=pki)
}

log_msg("Reading 4 signal-bearing datasets ...")
SEC <- lapply(SECT, function(s) {
  d <- read_dataset(parse_imzml(file.path(IMZDIR, s$imz)))
  d$cor <- sapply(seq_along(BANDS), function(w)
             suppressWarnings(cor(d$M[,REFB], d$M[,w], use="complete.obs")))
  c(s, d)
})
names(SEC) <- vapply(SECT, `[[`, character(1), "sid")

# ---------- pooled distribution metrics --------------------------------------
hw <- Reduce(`+`, lapply(SEC, `[[`, "hw"))
pk_mz <- unlist(lapply(SEC, `[[`, "pk_mz")); pk_in <- unlist(lapply(SEC, `[[`, "pk_in"))
keep  <- is.finite(pk_mz) & pk_in > 0
thr   <- as.numeric(quantile(pk_in[keep], 0.50))     # signal-bearing pixels
sig   <- keep & pk_in >= thr
pk    <- pk_mz[sig]; pk_med <- median(pk)

cl <- ctr >= CL_LO & ctr <= CL_HI
iwmean <- sum(ctr[cl]*hw[cl]) / sum(hw[cl])
sm <- as.numeric(stats::filter(hw, rep(1/5,5))); sm[!is.finite(sm)] <- 0
smc <- sm[cl]; ctc <- ctr[cl]
apex <- ctc[which.max(smc)]; hm <- max(smc)/2
above <- which(smc >= hm); lo_i <- min(above); hi_i <- max(above)
fwhm_mz <- ctc[hi_i]-ctc[lo_i]; fwhm_ppm <- fwhm_mz/apex*1e6; Rpow <- apex/fwhm_mz
log_msg("apex %.5f (%+.2f ppm) | iw-mean %.5f (%+.2f ppm) | FWHM %.1f ppm | R %.0f | per-px(signal) median %.5f (%+.2f ppm) | n_sig %d",
        apex, ppm(apex), iwmean, ppm(iwmean), fwhm_ppm, Rpow, pk_med, ppm(pk_med), length(pk))

# ================= RENDER ====================================================
pal <- viridis(256)
dir.create(dirname(OUT_PDF), showWarnings = FALSE, recursive = TRUE)
pdf(OUT_PDF, width = 11, height = 8.5)

## ---- PAGE 1 -----------------------------------------------------------------
XL <- c(191.010, 191.031); inXL <- ctr>=XL[1] & ctr<=XL[2]; yt <- max(hw[inXL])
par(fig=c(0,1,0.42,1), mar=c(3.6,4.8,5.0,1.6))
plot(ctr, hw, type="n", xlim=XL, xaxs="i", yaxs="i", ylim=c(0, yt*1.30),
     xlab="", ylab="summed ion intensity  (a.u.)", axes=FALSE)
segments(ctr, 0, ctr, hw, col="#bcbcbc", lwd=1.2)
lines(ctr, sm, col="#222222", lwd=2.0)
axis(1, at=seq(191.010,191.030,0.005)); axis(2, las=1)
pp <- seq(-40,60,20); axis(3, at=CIT*(1+pp/1e6), labels=sprintf("%+d", pp), col.axis="grey35", col="grey55")
mtext("ppm relative to citrate", side=3, line=1.9, cex=0.82, col="grey35")
mtext("m/z", side=1, line=2.3, cex=0.95)
abline(v=CIT, col="#c0392b", lwd=2.4); abline(v=iwmean, col="#1f5fa8", lwd=2.2, lty=2)
abline(v=GRID, col="#1f5fa8", lwd=1.1, lty=3)
text(CIT,  yt*1.27, "citrate [M-H]-\n191.0197 (theo.)", col="#c0392b", cex=0.74, font=2, adj=c(1.08,1))
text(iwmean, yt*1.27, sprintf("intensity-weighted centroid %.4f\n(%+.1f ppm; pulled up by co-isobars)", iwmean, ppm(iwmean)),
     col="#1f5fa8", cex=0.74, font=2, adj=c(-0.05,1))
yb <- yt*0.5
segments(ctc[lo_i], yb, ctc[hi_i], yb, col="grey20", lwd=1.4)
segments(ctc[lo_i], yb*0.92, ctc[lo_i], yb*1.08, col="grey20")
segments(ctc[hi_i], yb*0.92, ctc[hi_i], yb*1.08, col="grey20")
text((ctc[lo_i]+ctc[hi_i])/2, yb, sprintf("blended-peak envelope FWHM ~ %.0f ppm  (R ~ %s)",
     fwhm_ppm, format(round(Rpow,-3), big.mark=",")), pos=3, cex=0.74, col="grey20")
ya <- yt*0.20
arrows(CIT, ya, GRID, ya, length=0.05, angle=90, code=3, col="#7d3c98", lwd=1.8)
text((CIT+GRID)/2, ya, sprintf("citrate <-> report ion = %.1f ppm  <<  peak width  =>  NOT resolved", ppm(GRID)),
     pos=3, cex=0.76, col="#7d3c98", font=2)
title("How well can this timsTOF separate citrate at m/z 191?", cex.main=1.3, line=3.4, adj=0)

par(fig=c(0.085,0.40,0.69,0.92), mar=c(1.8,0.6,1.4,0.6), new=TRUE)
plot(ctr, log1p(hw), type="h", col="#9a9a9a", xlim=c(190.965,191.09), axes=FALSE, xlab="", ylab="")
abline(v=CIT, col="#c0392b", lwd=1.4)
axis(1, at=c(190.98,191.02,191.06), cex.axis=0.62, mgp=c(1,0.25,0), tck=-0.04); box(col="grey70")
text(191.046, log1p(max(hw)), "dominant\nisobar (+138 ppm)", col="#6f6f6f", cex=0.6, adj=c(0.5,1.05))
text(190.967, log1p(max(hw))*0.96, "full window (log)", cex=0.62, col="grey40", adj=0)

par(fig=c(0,1,0,0.42), mar=c(5.6,4.8,1.4,1.6), new=TRUE)
h2 <- hist(pk, breaks=seq(CL_LO,CL_HI,0.0002), plot=FALSE)
plot(h2, xlim=XL, col="#1f5fa855", border="#1f5fa8", main="", xaxs="i", xlab="", ylab="number of pixels", las=1)
abline(v=CIT, col="#c0392b", lwd=2.4); abline(v=pk_med, col="#1f5fa8", lwd=2.2, lty=2)
mtext("per-pixel strongest-centroid m/z in the 191 cluster", side=1, line=2.4, cex=0.95)
legend("topright", bty="n", cex=0.80,
       legend=c("citrate (theoretical) 191.0197",
                sprintf("per-pixel median %.4f (%+.1f ppm)", pk_med, ppm(pk_med)),
                sprintf("signal-bearing pixels (top 50%% by intensity), n = %s", format(length(pk), big.mark=","))),
       col=c("#c0392b","#1f5fa8",NA), lwd=c(2.4,2.2,NA), lty=c(1,2,NA))
mtext(sprintf("Even in strong-signal pixels, per-pixel centroids scatter across one ~%.0f-ppm-wide peak and sit ~%+.0f ppm above citrate on average -> citrate is not cleanly separated from heavier co-isobars.",
              fwhm_ppm, ppm(pk_med)), side=1, line=4.2, cex=0.74, col="grey25")

## ---- PAGE 2: sub-band ion images --------------------------------------------
SCALE_UM <- 500; PX_UM <- MSI_PIXEL_UM
render <- function(s, w, hi) {
  v <- s$M[, w]; xs <- sort(unique(s$x)); ys <- sort(unique(s$y))
  mat <- matrix(0, length(xs), length(ys)); mat[cbind(match(s$x,xs), match(s$y,ys))] <- v
  par(mar=c(0.5,0.5,1.7,0.5))
  image(xs, ys, pmin(mat/hi,1), col=pal, asp=1, useRaster=TRUE, zlim=c(0,1), axes=FALSE, xlab="", ylab="", main="")
  box(col="#444444", lwd=1.8)
  x0 <- xs[1]+0.05*(max(xs)-min(xs)); y0 <- ys[1]+0.06*(max(ys)-min(ys))
  segments(x0,y0,x0+SCALE_UM/PX_UM,y0, col="white", lwd=2.5, xpd=NA)
}
layout(matrix(1:25, nrow=5, byrow=TRUE), widths=c(0.66,1,1,1,1))
par(oma=c(3.4,0.5,4.8,0.5))
for (bi in seq_along(BANDS)) {
  b <- BANDS[[bi]]
  pos <- unlist(lapply(SEC, function(s){ vv<-s$M[,bi]; vv[vv>0] }))
  hi  <- if (length(pos)) as.numeric(quantile(pos, IMG_CLIP_HI, na.rm=TRUE)) else 1
  if (!is.finite(hi) || hi<=0) hi <- 1
  par(mar=c(0.5,0.4,1.7,0.2)); plot.new()
  text(0.02, 0.66, b$name, adj=c(0,0.5), font=2, cex=0.92,
       col=if(b$full) "#222222" else if(b$ref) "#c0392b" else "#1f5fa8")
  text(0.02, 0.44, sprintf("m/z %.4f  (%+.0f ppm)", b$c, ppm(b$c)), adj=c(0,0.5), cex=0.74, col="grey30")
  text(0.02, 0.26, sprintf("band %.4f-%.4f", b$lo, b$hi), adj=c(0,0.5), cex=0.64, col="grey45")
  if (b$full) text(0.02, 0.08, "= what the report sums as citrate", adj=c(0,0.5), cex=0.62, col="grey45", font=3)
  if (b$ref)  text(0.02, 0.08, "spatial reference (r=1)", adj=c(0,0.5), cex=0.62, col="grey45", font=3)
  for (si in seq_along(SEC)) {
    s <- SEC[[si]]; render(s, bi, hi)
    if (bi == 1) mtext(sprintf("%s (%s)", sub("^AO_","",s$sid), s$grp), side=3, line=0.25, cex=0.62,
                       col="#444444")
    if (!b$full && !b$ref) mtext(sprintf("r=%.2f vs citrate", s$cor[bi]), side=3, line=0.15, cex=0.6, col="grey35")
  }
}
mtext("Ion images of the \"citrate\" feature and its m/z sub-bands (4 sections with strong 191 signal)",
      outer=TRUE, line=2.8, cex=1.18, font=2)
mtext("Locked render: viridis, linear, gamma 1.0, per-band global p99.5 clip, 500 um scale bar.   r = spatial Pearson correlation of each band vs the citrate-mass band.",
      outer=TRUE, line=1.3, cex=0.72, col="grey30")
mtext("One molecule would give spatially identical sub-bands (r~1). Divergent maps / low r => the report's 191 \"citrate\" feature blends distinct ions the instrument cannot separate.",
      outer=TRUE, side=1, line=1.5, cex=0.74, col="grey25")

## ---- PAGE 3: with vs without (citrate band shifts, noise stays flat) --------
# Logic: if the 191.0197 signal is a real compound (citrate) its intensity should
# track sections that show it, while an empty m/z window (pure noise) does not.
WW <- list(
  list(sid="20h_sl4A_sec3b", grp="20h", g="with",    imz="06102026_ao_20h_sl4a_3b.imzML"),
  list(sid="0h_sl6A_sec2a",  grp="0h",  g="with",    imz="06102026_ao_0h_sl6a_sec2a.imzML"),
  list(sid="20h_sl4A_sec3a", grp="20h", g="with",    imz="06102026_ao_20h_sl4a_3a.imzML"),
  list(sid="0h_sl6A_sec5a",  grp="0h",  g="with",    imz="06102026_ao_0h_sl6a_sec5a.imzML"),
  list(sid="20h_sl4A_sec2a", grp="20h", g="with",    imz="06102026_ao_20h_sl4a_2a.imzML"),
  list(sid="0h_sl6A_sec1b",  grp="0h",  g="without", imz="06102026_ao_0h_sl6a_sec1b.imzML"),
  list(sid="0h_sl6A_sec4b",  grp="0h",  g="without", imz="06102026_ao_0h_sl6a_sec4b.imzML"),
  list(sid="20h_sl4A_sec4a", grp="20h", g="without", imz="06102026_ao_20h_sl4a_4a.imzML"),
  list(sid="20h_sl4A_sec4b", grp="20h", g="without", imz="06102026_ao_20h_sl4a_4b.imzML"),
  list(sid="20h_sl4A_sec5a", grp="20h", g="without", imz="06102026_ao_20h_sl4a_5a.imzML"))
CITW <- c(191.0188,191.0206)   # citrate-specific band (the "0 ppm" image)
NOIW <- c(190.9450,190.9468)   # empty m/z window, equal width = noise floor
ISOW <- c(191.0440,191.0482)   # dominant isobar (different compound, specificity ctrl)
FHL <- 191.005; FHH <- 191.030; fbrk <- seq(FHL,FHH,0.0002); fctr <- (head(fbrk,-1)+tail(fbrk,-1))/2

read_ww <- function(p) {
  con <- file(p$ibd, "rb"); on.exit(close(con)); n <- length(p$mz_off)
  cit <- numeric(n); noi <- numeric(n); iso <- numeric(n); tic <- 0; fh <- numeric(length(fctr))
  for (i in seq_len(n)) {
    seek(con, p$mz_off[i]); mz <- readBin(con,"double",p$mz_len[i],8,endian="little")
    seek(con, p$in_off[i]); ic <- readBin(con,"numeric",p$mz_len[i],4,endian="little")
    tic <- tic + sum(ic)
    s<-mz>=CITW[1]&mz<=CITW[2]; if(any(s)) cit[i]<-sum(ic[s])
    s<-mz>=NOIW[1]&mz<=NOIW[2]; if(any(s)) noi[i]<-sum(ic[s])
    s<-mz>=ISOW[1]&mz<=ISOW[2]; if(any(s)) iso[i]<-sum(ic[s])
    s<-which(mz>=FHL&mz<=FHH); if(length(s)){b<-findInterval(mz[s],fbrk,rightmost.closed=TRUE);ok<-b>=1&b<=length(fh);fh[b[ok]]<-fh[b[ok]]+ic[s][ok]}
  }
  list(cit=cit, noi=noi, iso=iso, tic=tic, fh=fh)
}
log_msg("Reading 10 with/without datasets ...")
WD <- lapply(WW, function(w){ p<-parse_imzml(file.path(IMZDIR,w$imz)); d<-read_ww(p); c(w, p["ibd"], list(x=p$x,y=p$y), d) })
gW <- vapply(WW,`[[`,character(1),"g")
citf <- sapply(WD,function(d) sum(d$cit)/d$tic)      # TIC-fraction per dataset
noif <- sapply(WD,function(d) sum(d$noi)/d$tic)
isof <- sapply(WD,function(d) sum(d$iso)/d$tic)
fold <- function(v) mean(v[gW=="with"])/mean(v[gW=="without"])
pval <- function(v) suppressWarnings(wilcox.test(v[gW=="with"], v[gW=="without"])$p.value)
log_msg("citrate fold %.2f p=%.4f | noise fold %.2f p=%.3f | isobar fold %.2f p=%.3f",
        fold(citf),pval(citf),fold(noif),pval(noif),fold(isof),pval(isof))
# shared image clip for the 10 citrate-band thumbnails (raw, comparable scale)
allcit <- unlist(lapply(WD,function(d) d$cit[d$cit>0]))
HI <- as.numeric(quantile(allcit, IMG_CLIP_HI, na.rm=TRUE))

frame()                                    # new page
thumb <- function(d, xr, yr) {
  xs<-sort(unique(d$x)); ys<-sort(unique(d$y)); m<-matrix(0,length(xs),length(ys))
  m[cbind(match(d$x,xs),match(d$y,ys))]<-d$cit
  par(fig=c(xr,yr), mar=c(0.3,0.3,1.45,0.3), new=TRUE)
  image(xs,ys,pmin(m/HI,1),col=pal,asp=1,useRaster=TRUE,zlim=c(0,1),axes=FALSE,xlab="",ylab="",main="")
  box(col="#444444", lwd=1.6)
  segments(xs[1]+0.05*diff(range(xs)), ys[1]+0.06*diff(range(ys)),
           xs[1]+0.05*diff(range(xs))+SCALE_UM/PX_UM, ys[1]+0.06*diff(range(ys)), col="white", lwd=2.2, xpd=NA)
}
xL<-0.075; xR<-0.625; cw<-(xR-xL)/5
yW<-c(0.50,0.76); yWO<-c(0.155,0.415)
for (i in 1:5) { d<-WD[[i]]; xr<-c(xL+(i-1)*cw, xL+(i-1)*cw+cw*0.94)
  thumb(d,xr,yW); par(fig=c(xr,yW),new=TRUE)
  mtext(sub("^.*?_","",d$sid), side=3, line=0.45, cex=0.58, col="#444444")
  mtext(sprintf("%.1f e-5 TIC", citf[i]*1e5), side=3, line=-0.15, cex=0.52, col="grey30") }
for (i in 1:5) { d<-WD[[i+5]]; xr<-c(xL+(i-1)*cw, xL+(i-1)*cw+cw*0.94)
  thumb(d,xr,yWO); par(fig=c(xr,yWO),new=TRUE)
  mtext(sub("^.*?_","",d$sid), side=3, line=0.45, cex=0.58, col="#444444")
  mtext(sprintf("%.1f e-5 TIC", citf[i+5]*1e5), side=3, line=-0.15, cex=0.52, col="grey30") }
# row labels + shared colourbar
par(fig=c(0,1,0,1), mar=c(0,0,0,0), new=TRUE); plot.new()
text(0.038, mean(yW),  "WITH\ncitrate",       srt=90, font=2, cex=0.9, col="#1a7a3a")
text(0.038, mean(yWO), "WITHOUT\n/ low",      srt=90, font=2, cex=0.9, col="#8a8a8a")
cbx<-seq(xL,xR,length.out=257); cby<-c(0.105,0.122)
rect(head(cbx,-1),cby[1],tail(cbx,-1),cby[2],col=pal,border=NA,xpd=NA)
text(xL,cby[1]-0.012,"0",cex=0.6,xpd=NA); text(xR,cby[1]-0.012,sprintf("%.2g (p99.5)",HI),cex=0.6,pos=2,xpd=NA)
text((xL+xR)/2, cby[2]+0.012, "citrate band 191.0188-191.0206  (shared raw intensity scale)", cex=0.62, col="grey30", xpd=NA)

# right-top: group-mean spectra overlay
par(fig=c(0.69,0.985,0.55,0.80), mar=c(3.0,3.4,2.0,0.8), new=TRUE)
ff <- sapply(WD,function(d) d$fh/d$tic)
mw<-rowMeans(ff[,gW=="with"]); mo<-rowMeans(ff[,gW=="without"])
plot(fctr,mw,type="n",xlim=c(FHL,FHH),ylim=c(0,max(mw)*1.1),axes=FALSE,xlab="",ylab="")
polygon(c(FHL,fctr,FHH),c(0,mw,0),col="#1a7a3a33",border=NA); lines(fctr,mw,col="#1a7a3a",lwd=2)
lines(fctr,mo,col="#8a8a8a",lwd=2)
abline(v=CIT,col="#c0392b",lwd=1.6,lty=2)
axis(1,at=c(191.01,191.02,191.03),cex.axis=0.7,mgp=c(2,0.4,0)); axis(2,las=1,cex.axis=0.7,mgp=c(2,0.5,0))
mtext("m/z",1,line=1.7,cex=0.72); mtext("mean intensity / TIC",2,line=2.3,cex=0.72)
title("Group-mean spectra @ m/z 191",cex.main=0.92,line=0.6)
legend("topright",bty="n",cex=0.7,legend=c("with citrate","without / low","citrate 191.0197"),
       col=c("#1a7a3a","#8a8a8a","#c0392b"),lwd=c(2,2,1.6),lty=c(1,1,2))

# right-bottom: quantification dotplot (log y): citrate vs noise, with vs without
par(fig=c(0.69,0.985,0.09,0.46), mar=c(3.4,3.8,2.0,0.8), new=TRUE)
xs4 <- c(1,2,3.2,4.2); vals<-list(citf[gW=="with"],citf[gW=="without"],noif[gW=="with"],noif[gW=="without"])
cols4<-c("#1a7a3a","#8a8a8a","#1a7a3a","#8a8a8a")
yl<-range(c(citf,noif[noif>0],min(citf)/3)); yl[1]<-max(yl[1],min(c(citf,noif[noif>0]))/2)
plot(NA,xlim=c(0.5,4.7),ylim=yl,log="y",axes=FALSE,xlab="",ylab="")
for(k in 1:4){ set.seed(k); jx<-xs4[k]+runif(length(vals[[k]]),-0.12,0.12)
  points(jx,pmax(vals[[k]],yl[1]),pch=21,bg=cols4[k],col="grey20",cex=1.2)
  segments(xs4[k]-0.18,mean(vals[[k]]),xs4[k]+0.18,mean(vals[[k]]),lwd=2,col=cols4[k]) }
axis(2,las=1,cex.axis=0.7,mgp=c(2,0.5,0)); box(col="grey70")
axis(1,at=c(1.5,3.7),labels=c("citrate band","noise band"),tick=FALSE,line=0.6,cex.axis=0.85,font=2)
axis(1,at=xs4,labels=c("with","w/o","with","w/o"),tick=FALSE,line=-0.4,cex.axis=0.65)
abline(v=2.6,col="grey85")
mtext("intensity / TIC  (log)",2,line=2.6,cex=0.72)
text(1.5, yl[2], sprintf("%.1fx, p=%.3f", fold(citf), pval(citf)), col="#1a7a3a", font=2, cex=0.78, pos=1)
text(3.7, yl[2], sprintf("%.1fx, p=%.2f", fold(noif), pval(noif)), col="grey40", font=2, cex=0.78, pos=1)
title("Relative abundance shift",cex.main=0.92,line=0.6)

par(fig=c(0,1,0,1),mar=c(0,0,0,0),new=TRUE); plot.new()
text(0.5,0.975,"Validation: the citrate-specific band tracks sections that show it, while the noise floor does not",
     font=2,cex=1.15)
text(0.5,0.945,sprintf("citrate band rises %.1fx in 'with' sections (p=%.3f); equal-width empty noise window is unchanged (%.1fx, p=%.2f); the unrelated +138 ppm isobar does not track citrate grouping (%.1fx, p=%.2f).",
     fold(citf),pval(citf),fold(noif),pval(noif),fold(isof),pval(isof)), cex=0.74, col="grey25")
text(0.5,0.045,"Constant noise + proportional citrate-band change = the 191.0197 signal is a real, abundance-variable compound consistent with citrate (not noise / not a flat artefact).",
     cex=0.78, col="grey25", font=3)

dev.off()
log_msg("DONE -> %s", OUT_PDF)
cat(sprintf("\ncitrate theo %.4f | apex %.4f (%+.2f ppm) | per-px(signal) median %.4f (%+.2f ppm) | FWHM %.0f ppm | R~%.0f\n",
            CIT, apex, ppm(apex), pk_med, ppm(pk_med), fwhm_ppm, Rpow))
for (s in SEC) cat(sprintf("  %-20s r(apex,high,isobar vs citrate) = %.2f, %.2f, %.2f\n",
                           s$sid, s$cor[3], s$cor[4], s$cor[5]))
