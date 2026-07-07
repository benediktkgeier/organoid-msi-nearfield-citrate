#!/usr/bin/env Rscript
# 04_citrate_window_images.R
# Small report: citrate [M-H]- (191.0197) ion images for the 5 sections with
# genuine citrate signal, integrated over three m/z windows: +-2.5, +-3, +-5 ppm.
# Lets the user judge how the citrate image looks as the integration window widens
# toward the isobar at 191.0211 (midpoint 191.0204 = +3.6 ppm).
#
# Locked render: viridis, linear, gamma 1.0, per-window global p99.5 clip across
# the 5 sections, 500 um scale bar.
# Out: figures/metabolites/citrate_window_images.pdf (+ png)
# Usage: Rscript R/07_metabolite_id/04_citrate_window_images.R

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_paths.R"))
suppressPackageStartupMessages(library(viridisLite))

CIT <- 191.019726; ISO <- 191.0211
OUT_PDF <- file.path(FIG_DIR, "metabolites", "citrate_window_images.pdf")
IMZDIR  <- file.path(CACHE_DIR, "imzml")
SCALE_UM <- 500; PX_UM <- MSI_PIXEL_UM; pal <- viridis(256)
ppm <- function(a,b=CIT)(a-b)/b*1e6
log_msg <- function(...) message(sprintf("[83] %s", sprintf(...)))

SECT <- list(
  list(sid="0h_sl6A_sec2a",  grp="0h",  imz="06102026_ao_0h_sl6a_sec2a.imzML"),
  list(sid="0h_sl6A_sec5a",  grp="0h",  imz="06102026_ao_0h_sl6a_sec5a.imzML"),
  list(sid="20h_sl4A_sec2a", grp="20h", imz="06102026_ao_20h_sl4a_2a.imzML"),
  list(sid="20h_sl4A_sec3a", grp="20h", imz="06102026_ao_20h_sl4a_3a.imzML"),
  list(sid="20h_sl4A_sec3b", grp="20h", imz="06102026_ao_20h_sl4a_3b.imzML"))
PPMW <- c(2.5, 3, 5)
WINS <- lapply(PPMW, function(p) c(CIT*(1-p/1e6), CIT*(1+p/1e6)))

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
read_wins <- function(p) {
  con<-file(p$ibd,"rb"); on.exit(close(con)); n<-length(p$mz_off); M<-matrix(0,n,length(WINS))
  glo<-min(sapply(WINS,`[`,1)); ghi<-max(sapply(WINS,`[`,2))
  for(i in seq_len(n)){
    seek(con,p$mz_off[i]); mz<-readBin(con,"double",p$mz_len[i],8,endian="little")
    if(!any(mz>=glo&mz<=ghi)) next
    seek(con,p$in_off[i]); ic<-readBin(con,"numeric",p$mz_len[i],4,endian="little")
    for(w in seq_along(WINS)){s<-mz>=WINS[[w]][1]&mz<=WINS[[w]][2]; if(any(s)) M[i,w]<-sum(ic[s])}
  }
  M
}
log_msg("Reading 5 datasets ...")
SEC <- lapply(SECT, function(s){ p<-parse_imzml(file.path(IMZDIR,s$imz)); c(s, list(x=p$x,y=p$y,M=read_wins(p))) })

# per-window global p99.5 clip across the 5 sections
HI <- sapply(seq_along(WINS), function(w){
  pos<-unlist(lapply(SEC,function(s){v<-s$M[,w]; v[v>0]})); as.numeric(quantile(pos,IMG_CLIP_HI,na.rm=TRUE)) })

render <- function(s, w, hi) {
  xs<-sort(unique(s$x)); ys<-sort(unique(s$y)); m<-matrix(0,length(xs),length(ys))
  m[cbind(match(s$x,xs),match(s$y,ys))]<-s$M[,w]
  par(mar=c(0.4,0.4,1.6,0.4))
  image(xs,ys,pmin(m/hi,1),col=pal,asp=1,useRaster=TRUE,zlim=c(0,1),axes=FALSE,xlab="",ylab="",main="")
  box(col="#444444",lwd=1.7)
  segments(xs[1]+0.05*diff(range(xs)), ys[1]+0.06*diff(range(ys)),
           xs[1]+0.05*diff(range(xs))+SCALE_UM/PX_UM, ys[1]+0.06*diff(range(ys)), col="white", lwd=2.4, xpd=NA)
}

dir.create(dirname(OUT_PDF),showWarnings=FALSE,recursive=TRUE)
pdf(OUT_PDF, width=11, height=7)
layout(matrix(1:18, nrow=3, byrow=TRUE), widths=c(0.55,1,1,1,1,1))
par(oma=c(2.6,0.5,4.4,0.5))
for (w in seq_along(WINS)) {
  par(mar=c(0.4,0.3,1.6,0.2)); plot.new()
  text(0.04, 0.66, sprintf("+-%.1f ppm", PPMW[w]), adj=c(0,0.5), font=2, cex=1.15)
  text(0.04, 0.44, sprintf("%.4f-%.4f", WINS[[w]][1], WINS[[w]][2]), adj=c(0,0.5), cex=0.72, col="grey30")
  text(0.04, 0.26, sprintf("width %.4f Da", diff(WINS[[w]])), adj=c(0,0.5), cex=0.66, col="grey45")
  text(0.04, 0.08, sprintf("p99.5 clip %.2g", HI[w]), adj=c(0,0.5), cex=0.64, col="grey45")
  # viridis intensity colorbar for this window-row (shared global p99.5 clip)
  np<-128; xb<-seq(0.04,0.60,length.out=np+1)
  rect(xb[-(np+1)],0.80,xb[-1],0.90,col=viridis(np),border=NA,xpd=NA)
  rect(0.04,0.80,0.60,0.90,border="black",lwd=0.5,xpd=NA)
  text(0.04,0.955,"0",adj=c(0,0.5),cex=0.56,xpd=NA)
  text(0.60,0.955,sprintf("%.2g",HI[w]),adj=c(1,0.5),cex=0.56,xpd=NA)
  text(0.04,0.985,"intensity (viridis)",adj=c(0,0.5),cex=0.5,col="grey40",xpd=NA)
  for (si in seq_along(SEC)) {
    s<-SEC[[si]]; render(s,w,HI[w])
    if (w==1) mtext(sprintf("%s (%s)", s$sid, s$grp), side=3, line=0.3, cex=0.7,
                    col="#444444")
  }
}
mtext("Citrate [M-H]- 191.0197 ion images at three integration windows (5 sections with citrate signal)",
      outer=TRUE, line=2.6, cex=1.18, font=2)
mtext("viridis, linear, gamma 1.0, per-window global p99.5 clip, 500 um bar.   Isobar at 191.0211 (+7.2 ppm); midpoint 191.0204 (+3.6 ppm) -> +-5 ppm window starts to include it.",
      outer=TRUE, line=1.1, cex=0.74, col="grey30")
mtext("Wider window = more citrate signal but more isobar bleed-through (upper edge of +-5 ppm = 191.0207, past the midpoint).",
      outer=TRUE, side=1, line=1.2, cex=0.74, col="grey25")
dev.off()
log_msg("DONE -> %s", OUT_PDF)
cat(sprintf("windows: +-2.5 [%.4f-%.4f] | +-3 [%.4f-%.4f] | +-5 [%.4f-%.4f]\n",
    WINS[[1]][1],WINS[[1]][2],WINS[[2]][1],WINS[[2]][2],WINS[[3]][1],WINS[[3]][2]))
cat(sprintf("p99.5 clips: %.3g, %.3g, %.3g\n", HI[1],HI[2],HI[3]))
