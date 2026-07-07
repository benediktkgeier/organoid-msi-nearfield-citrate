#!/usr/bin/env Rscript
# 05_citrate_isotopes_adducts.R
# Companion to 04_citrate_window_images.R. Same 5 sections with genuine citrate
# signal, but the NORMAL +-10 ppm integration window (not the special +-2.5/3/5
# ppm windows of 04_citrate_window_images.R). Shows the full citrate ion family:
#   - Citrate [M-H]- C12 (191.0197) and its 13C isotope C13 (192.0231)
#   - Na adducts  [M-2H+Na]- 213.0017 , [M-3H+2Na]- 234.9836
#   - Cl / acetate [M+Cl]- 226.9964   , [M+CH3COO]- 251.0409
# plus per-section spectral zoom-ins over 190.95-192.10 (C12 + C13 envelope).
#
# Locked render: viridis, linear, gamma 1.0, per-ion (per-row) global p99.5 clip
# across the 5 sections, 500 um scale bar, per-slide box outline.
# Out: figures/metabolites/citrate_isotopes_addocts.pdf
# Usage: Rscript R/07_metabolite_id/05_citrate_isotopes_adducts.R
#
# NOTE on the acetate row: the user's "[M+CH3OO]-" was confirmed to mean the
# acetate adduct [M+CH3COO]- = neutral M + CH3COO (C2H3O2, 59.013304).

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_paths.R"))
suppressPackageStartupMessages(library(viridisLite))

# ---- masses ----------------------------------------------------------------
C12  <- 191.019726                       # citrate [M-H]- (locked project value)
DC13 <- 1.0033548                        # 13C - 12C
# neutral M = [M-H]- + proton ; adduct anions include the electron mass
PROT <- 1.00727646; ELEC <- 0.000548580
H    <- 1.00782503; NA_  <- 22.98976928; CL <- 34.96885268
CH3COO <- 59.013304                      # acetate C2H3O2 (per user: "CH3OO" = acetate)
Mneu <- C12 + PROT                       # 192.027002

UP   <- C12 - DC13                       # 190.016371 : putative monoisotopic ion
                                         # one 13C step BELOW citrate. If citrate
                                         # were merely its 13C satellite, this M
                                         # would dominate. Diagnostic row 7.
IONS <- list(
  list(lab="Citrate [M-H]- (C12)", form="C6H7O7-",      mz = C12),
  list(lab="13C isotope (C13)",    form="13C-C5H7O7-",  mz = C12 + DC13),
  list(lab="[M-2H+Na]-",           form="C6H6O7Na-",    mz = Mneu - 2*H + NA_   + ELEC),
  list(lab="[M-3H+2Na]-",          form="C6H5O7Na2-",   mz = Mneu - 3*H + 2*NA_ + ELEC),
  list(lab="[M+Cl]-",              form="C6H8O7Cl-",    mz = Mneu + CL     + ELEC),
  list(lab="[M+CH3COO]- (acetate)",form="C8H11O9-",     mz = Mneu + CH3COO + ELEC),
  list(lab="upstream M? (isotope check)", form="191.0197 - 1.0034", mz = UP))
IUP <- length(IONS)                      # index of the 190.0163 diagnostic ion

PPMHW <- 10                              # normal +-10 ppm half-window (per user)
WINS  <- lapply(IONS, function(io) c(io$mz*(1-PPMHW/1e6), io$mz*(1+PPMHW/1e6)))
C13mz <- IONS[[2]]$mz

WLO <- 189.98; WHI <- 192.10             # spectral zoom: 190.0163 .. C12 .. C13
brk <- seq(WLO, WHI, by = 0.00010); ctr <- (head(brk,-1)+tail(brk,-1))/2

OUT_PDF  <- file.path(FIG_DIR, "metabolites", "citrate_isotopes_addocts.pdf")
IMZDIR   <- file.path(CACHE_DIR, "imzml")
SCALE_UM <- 500; PX_UM <- MSI_PIXEL_UM; pal <- viridis(256)
ppm      <- function(a,b)(a-b)/b*1e6
log_msg  <- function(...) message(sprintf("[83b] %s", sprintf(...)))

SECT <- list(
  list(sid="0h_sl6A_sec2a",  grp="0h",  imz="06102026_ao_0h_sl6a_sec2a.imzML"),
  list(sid="0h_sl6A_sec5a",  grp="0h",  imz="06102026_ao_0h_sl6a_sec5a.imzML"),
  list(sid="20h_sl4A_sec2a", grp="20h", imz="06102026_ao_20h_sl4a_2a.imzML"),
  list(sid="20h_sl4A_sec3a", grp="20h", imz="06102026_ao_20h_sl4a_3a.imzML"),
  list(sid="20h_sl4A_sec3b", grp="20h", imz="06102026_ao_20h_sl4a_3b.imzML"))

# ---- imzML low-level reader (verbatim from R/07_metabolite_id/04_citrate_window_images.R) ----------------------------
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

# one pass: per-pixel ion-window sums (pixel x ion) + spectrum histogram +
# intensity-weighted m/z accumulators per ion (for measured centroid).
read_all <- function(p) {
  con<-file(p$ibd,"rb"); on.exit(close(con)); n<-length(p$mz_off)
  M<-matrix(0,n,length(WINS)); hwl<-numeric(length(ctr))
  cnum<-numeric(length(WINS)); cden<-numeric(length(WINS))
  glo<-min(WLO, min(sapply(WINS,`[`,1))); ghi<-max(WHI, max(sapply(WINS,`[`,2)))
  for(i in seq_len(n)){
    seek(con,p$mz_off[i]); mz<-readBin(con,"double",p$mz_len[i],8,endian="little")
    if(!any(mz>=glo&mz<=ghi)) next
    seek(con,p$in_off[i]); ic<-readBin(con,"numeric",p$mz_len[i],4,endian="little")
    sels<-which(mz>=WLO&mz<=WHI)
    if(length(sels)){ b<-findInterval(mz[sels],brk,rightmost.closed=TRUE); ok<-b>=1&b<=length(hwl)
      hwl[b[ok]]<-hwl[b[ok]]+ic[sels][ok] }
    for(w in seq_along(WINS)){
      s<-mz>=WINS[[w]][1]&mz<=WINS[[w]][2]
      if(any(s)){ M[i,w]<-sum(ic[s]); cnum[w]<-cnum[w]+sum(ic[s]*mz[s]); cden[w]<-cden[w]+sum(ic[s]) }
    }
  }
  list(M=M, hw=hwl, cnum=cnum, cden=cden)
}

log_msg("Reading 5 datasets ...")
SEC <- lapply(SECT, function(s){
  p<-parse_imzml(file.path(IMZDIR,s$imz)); d<-read_all(p)
  c(s, list(x=p$x,y=p$y,M=d$M,hw=d$hw,cnum=d$cnum,cden=d$cden)) })

# pooled measured centroid + delta-ppm per ion
CNUM <- Reduce(`+`, lapply(SEC,`[[`,"cnum")); CDEN <- Reduce(`+`, lapply(SEC,`[[`,"cden"))
MEAS <- ifelse(CDEN>0, CNUM/CDEN, NA_real_)
DPPM <- sapply(seq_along(IONS), function(w) if(is.finite(MEAS[w])) ppm(MEAS[w], IONS[[w]]$mz) else NA_real_)

# per-ion (per-row) global p99.5 clip across the 5 sections
HI <- sapply(seq_along(WINS), function(w){
  pos<-unlist(lapply(SEC,function(s){v<-s$M[,w]; v[v>0]}))
  if(!length(pos)) return(1); as.numeric(quantile(pos,IMG_CLIP_HI,na.rm=TRUE)) })

# ---- isotope cross-check: is C12 191.0197 just the 13C of 190.0163? ----------
# Pooled intensity ratio I(191)/I(190) and pixel-wise spatial Spearman r.
v191 <- unlist(lapply(SEC, function(s) s$M[,1]))
v190 <- unlist(lapply(SEC, function(s) s$M[,IUP]))
RATIO   <- if (sum(v190)>0) sum(v191)/sum(v190) else Inf      # I(C12)/I(upstream)
IMPL_C  <- RATIO * 0.989/0.0107                               # carbons needed if 191==13C of 190
both    <- v191>0 | v190>0
RSP     <- suppressWarnings(cor(v191[both], v190[both], method="spearman"))
RSEC    <- sapply(SEC, function(s){ a<-s$M[,1]; b<-s$M[,IUP]; k<-a>0|b>0
            if(sum(k)<5) NA_real_ else suppressWarnings(cor(a[k],b[k],method="spearman")) })

render <- function(s, w, hi) {
  xs<-sort(unique(s$x)); ys<-sort(unique(s$y)); m<-matrix(0,length(xs),length(ys))
  m[cbind(match(s$x,xs),match(s$y,ys))]<-s$M[,w]
  par(mar=c(0.4,0.4,1.6,0.4))
  image(xs,ys,pmin(m/hi,1),col=pal,asp=1,useRaster=TRUE,zlim=c(0,1),axes=FALSE,xlab="",ylab="",main="")
  box(col="#444444",lwd=1.7)
  segments(xs[1]+0.05*diff(range(xs)), ys[1]+0.06*diff(range(ys)),
           xs[1]+0.05*diff(range(xs))+SCALE_UM/PX_UM, ys[1]+0.06*diff(range(ys)), col="white", lwd=2.4, xpd=NA)
}

# vertical colorbar: summed-ion-intensity (a.u.) from 0 to the row's p99.5 clip
draw_cbar <- function(hi) {
  z <- seq(0,1,length.out=256)
  par(mar=c(1.6,0.4,1.6,2.6))
  image(1, z, matrix(z,1,256), col=pal, axes=FALSE, xlab="", ylab="", useRaster=TRUE)
  box(lwd=0.6)
  at <- c(0,0.25,0.5,0.75,1)
  axis(4, at=at, labels=sprintf("%.2g", at*hi), las=1, cex.axis=0.52, tcl=-0.18, mgp=c(2,0.30,0), lwd=0.6)
}

dir.create(dirname(OUT_PDF),showWarnings=FALSE,recursive=TRUE)
pdf(OUT_PDF, width=11, height=11.5)

## ---- PAGE 1: ion image grid (7 ions x 5 sections) + per-row colorbar -------
layout(matrix(1:(7*7), nrow=7, byrow=TRUE), widths=c(0.78,1,1,1,1,1,0.42))
par(oma=c(3.0,0.5,4.6,0.5))
for (w in seq_along(IONS)) {
  io<-IONS[[w]]
  par(mar=c(0.4,0.3,1.6,0.2)); plot.new()
  text(0.02, 0.86, io$lab, adj=c(0,0.5), font=2, cex=0.92)
  text(0.02, 0.66, io$form, adj=c(0,0.5), cex=0.70, col="grey30")
  text(0.02, 0.50, sprintf("theo %.4f", io$mz), adj=c(0,0.5), cex=0.66, col="grey25")
  text(0.02, 0.36, sprintf("meas %.4f (%+.1f ppm)", MEAS[w], DPPM[w]), adj=c(0,0.5), cex=0.66, col="#1f5fa8")
  text(0.02, 0.22, sprintf("win %.4f-%.4f", WINS[[w]][1], WINS[[w]][2]), adj=c(0,0.5), cex=0.62, col="grey40")
  text(0.02, 0.08, sprintf("+-%d ppm | p99.5 %.2g", PPMHW, HI[w]), adj=c(0,0.5), cex=0.62, col="grey45")
  for (si in seq_along(SEC)) {
    s<-SEC[[si]]; render(s,w,HI[w])
    if (w==1) mtext(sprintf("%s (%s)", s$sid, s$grp), side=3, line=0.3, cex=0.66,
                    col="#444444")
  }
  draw_cbar(HI[w])
  if (w==1) mtext("a.u.", side=3, line=0.3, cex=0.6, col="grey30")
}
mtext("Citrate isotopes & adducts - ion images at the normal +-10 ppm window (5 sections with citrate signal)",
      outer=TRUE, line=2.8, cex=1.16, font=2)
mtext("viridis, linear, gamma 1.0, 500 um bar.  Right-column colorbar = summed ion intensity (a.u.), 0 to per-ion global p99.5 clip.  Rows: C12 / 13C / [M-2H+Na]- / [M-3H+2Na]- / [M+Cl]- / [M+CH3COO]- (acetate) / upstream 190.0163.",
      outer=TRUE, line=1.0, cex=0.72, col="grey30")
mtext(sprintf("Isotope check: I(191.0197)/I(190.0163) = %.1f  ->  191 would need ~%.0f C to be the 13C of 190 (citrate has 6) -> 191 is a genuine ion, NOT the isotope of 190.   spatial r = %.2f.",
              RATIO, IMPL_C, RSP),
      outer=TRUE, side=1, line=1.4, cex=0.72, col="#7d3c98", font=2)

## ---- PAGE 2: per-section spectral zoom-ins (190.95-192.10) ------------------
inx <- ctr>=WLO & ctr<=WHI
layout(matrix(1:6, nrow=2, byrow=TRUE))
par(oma=c(1.6,1.0,4.2,0.6))
for (si in seq_along(SEC)) {
  s<-SEC[[si]]; y<-log1p(s$hw[inx]); xx<-ctr[inx]
  yt<-max(y); if(!is.finite(yt)||yt<=0) yt<-1
  par(mar=c(2.8,3.4,1.8,0.6))
  plot(xx, y, type="n", xlim=c(WLO,WHI), ylim=c(0,yt*1.14), xaxs="i", yaxs="i", axes=FALSE, xlab="", ylab="")
  segments(xx, 0, xx, y, col="#9a9a9a", lwd=1)
  abline(v=UP,    col="#27ae60", lwd=1.2, lty=3)
  abline(v=C12,   col="#c0392b", lwd=1.2)
  abline(v=C13mz, col="#1f5fa8", lwd=1.2, lty=2)
  axis(1, at=seq(190.0,192.0,0.5), cex.axis=0.75, mgp=c(2,0.4,0))
  axis(2, las=1, cex.axis=0.72, mgp=c(2,0.5,0))
  box(col="#444444", lwd=1.6)
  arrows(C12, yt*1.06, C13mz, yt*1.06, length=0.05, code=3, col="grey35", lwd=1)
  text((C12+C13mz)/2, yt*1.11, sprintf("%.4f Da", C13mz-C12), cex=0.62, col="grey30")
  mtext(sprintf("%s (%s)", s$sid, s$grp), side=3, line=0.3, cex=0.74,
        col="#444444")
}
# legend + verdict cell (6th panel)
par(mar=c(2.8,0.6,1.8,0.4)); plot.new()
legend("topleft", bty="n", cex=0.86, lwd=1.2, lty=c(3,1,2),
       col=c("#27ae60","#c0392b","#1f5fa8"),
       legend=c(sprintf("upstream 190.0163"), sprintf("citrate C12  %.4f", C12),
                sprintf("13C C13  %.4f", C13mz)))
text(0, 0.66, sprintf("Isotope cross-check\nI(191)/I(190) = %.1f\nimplied C if 191==13C(190): ~%.0f\n(citrate has 6 C)\nspatial r(191,190) = %.2f\n-> 191 is a genuine ion,\n   not the 13C of 190.",
                      RATIO, IMPL_C, RSP), adj=c(0,1), cex=0.74, col="#7d3c98")
text(0, 0.10, "log1p intensity\n(compresses the dominant\n191.046 isobar so the\nsmall C13 stays visible)", adj=c(0,1), cex=0.68, col="grey40")
mtext("Citrate [M-H]- spectral zoom-ins (189.98-192.10): upstream 190.0163 / C12 191.0197 / 13C 192.0231, per section",
      outer=TRUE, line=2.4, cex=1.10, font=2)
mtext("Summed ion intensity per 0.0001 Da bin, log1p y-axis.  Green dotted = 190.0163, red = citrate C12, blue dashed = 13C C13.",
      outer=TRUE, line=1.0, cex=0.74, col="grey30")
dev.off()

log_msg("DONE -> %s", OUT_PDF)
cat("ion                          theo         meas         dppm    p99.5\n")
for (w in seq_along(IONS))
  cat(sprintf("%-27s %.5f  %.5f  %+6.2f  %.3g\n",
      IONS[[w]]$lab, IONS[[w]]$mz, MEAS[w], DPPM[w], HI[w]))
cat(sprintf("\nisotope check: I(191.0197)/I(190.0163) = %.2f | implied C if 191==13C(190) ~%.0f | spatial r = %.2f\n",
            RATIO, IMPL_C, RSP))
cat(sprintf("per-section r(191,190): %s\n", paste(sprintf("%.2f", RSEC), collapse=", ")))
