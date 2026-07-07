#!/usr/bin/env Rscript
# 02_metabolite_report.R
# One metabolite per page x 4 display datasets (2 from each slide block). Ion
# images for every reference/HMDB compound matched within 20 ppm (from
# R/07_metabolite_id/01_metabolite_match.R).
# v3 locked rendering: viridis, linear, GLOBAL p99.5 clip across the 4 panels
# per ion, gamma 1.0, NA->0, scale bar (500 um). Page title carries the accurate
# theoretical m/z, matched data m/z, delta (mDa + ppm), score, HMDB, source.
#
# Out: figures/metabolites/metabolite_report.pdf
# Usage: Rscript R/07_metabolite_id/02_metabolite_report.R

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_paths.R"))
suppressPackageStartupMessages({ library(Cardinal); library(viridisLite) })

TARGET <- c("AO_0h_sl6A_sec1a"="0h", "AO_0h_sl6A_sec4b"="0h",
            "AO_20h_sl4A_sec5a"="20h", "AO_20h_sl4A_sec2b"="20h")
SIDS <- names(TARGET)
MATCH_CSV <- file.path(RES_DIR, "metabolites", "metabolite_match_table.csv")
OUT_PDF   <- file.path(FIG_DIR, "metabolites", "metabolite_report.pdf")
dir.create(dirname(OUT_PDF), showWarnings = FALSE, recursive = TRUE)
SCALE_UM <- 500; PX_UM <- MSI_PIXEL_UM   # 500 um bar, 10 um/px
log_msg <- function(...) message(sprintf("[81] %s", sprintf(...)))

mt <- read.csv(MATCH_CSV, stringsAsFactors = FALSE)
mt <- mt[mt$matched == TRUE, ]
mt <- mt[order(-mt$score, mt$theo_mz), ]
log_msg("Matched compounds to render: %d", nrow(mt))

log_msg("Loading peaks_after_freq.rds + slicing 4 datasets...")
mse <- readRDS(file.path(CACHE_DIR, "peaks_after_freq.rds"))
# peaks_after_freq.rds predates the sample_id annotation in 01_preprocess.R, so
# derive sample_id from the run() factor + inventory run_name mapping.
inv <- load_inventory()
inv$run_name <- sub("\\.imzML$", "", basename(inv$imzml_path), ignore.case = TRUE)
run_to_sid <- setNames(inv$sample_id, inv$run_name)
px_sid <- unname(run_to_sid[as.character(run(mse))])
co <- as.data.frame(coord(mse)); stopifnot(all(c("x","y") %in% names(co)))
cols <- which(px_sid %in% SIDS)
feat <- sort(unique(mt$feat_idx))
sp  <- as.matrix(spectra(mse)[feat, cols, drop = FALSE])   # features x pixels(4 sec)
rownames(sp) <- as.character(feat)
sid4 <- px_sid[cols]; x4 <- co$x[cols]; y4 <- co$y[cols]   # aligned to sp columns
ps <- lapply(SIDS, function(sid) {
  m <- which(sid4 == sid)
  x <- x4[m]; y <- y4[m]; xs <- sort(unique(x)); ys <- sort(unique(y))
  list(sid=sid, grp=unname(TARGET[sid]), cols=m, xs=xs, ys=ys,
       ix=match(x,xs), iy=match(y,ys), nx=length(xs), ny=length(ys))
})
names(ps) <- SIDS
log_msg("Section pixel counts: %s",
        paste(sprintf("%s=%d", SIDS, sapply(ps, function(p) length(p$cols))), collapse=", "))

pal <- viridis(256)
panel <- function(p, vals_all, hi, title) {
  v <- vals_all[p$cols]
  mat <- matrix(0, p$nx, p$ny); mat[cbind(p$ix, p$iy)] <- v
  if (!is.finite(hi) || hi <= 0) hi <- 1
  sc <- pmin(mat/hi, 1)
  par(mar = c(1.5, 1, 2.5, 3.5))
  image(p$xs, p$ys, sc, col = pal, asp = 1, useRaster = TRUE, zlim = c(0,1),
        axes = FALSE, xlab = "", ylab = "",
        main = sprintf("%s  (%s)", p$sid, p$grp), cex.main = 0.85,
        col.main = "#444444")
  box(col = "#444444", lwd = 2)
  # colorbar (shared clip)
  cx0<-grconvertX(1.03,"npc","user"); cx1<-grconvertX(1.08,"npc","user")
  cy0<-grconvertY(0.05,"npc","user"); cy1<-grconvertY(0.95,"npc","user")
  yb<-seq(cy0,cy1,length.out=257)
  rect(cx0,head(yb,-1),cx1,tail(yb,-1),col=pal,border=NA,xpd=NA)
  text(cx1,grconvertY(c(0.05,0.95),"npc","user"),sprintf("%.2g",c(0,hi)),
       pos=4,cex=0.55,xpd=NA,offset=0.15)
  # scale bar: 500 um, bottom-left, white, label centered below
  blen <- SCALE_UM/PX_UM
  x0 <- p$xs[1] + 0.04*(max(p$xs)-min(p$xs)); y0 <- p$ys[1] + 0.06*(max(p$ys)-min(p$ys))
  segments(x0, y0, x0+blen, y0, col="white", lwd=3, xpd=NA)
  text(x0+blen/2, y0 - 0.045*(max(p$ys)-min(p$ys)),
       sprintf("%d um", SCALE_UM), col="white", cex=0.55, xpd=NA)
}

log_msg("Rendering PDF...")
pdf(OUT_PDF, width = 11, height = 8.5)

# ---- index pages -----------------------------------------------------------
idx <- mt[, c("name","adduct","theo_mz","data_mz","delta_ppm","score","in_analysis_set","source")]
per_pg <- 40
np <- ceiling(nrow(idx)/per_pg)
for (pg in seq_len(np)) {
  par(mfrow=c(1,1), mar=c(1,1,2,1)); plot.new()
  title(sprintf("Metabolite report - index (%d matched within 20 ppm)  [page %d/%d]",
                nrow(mt), pg, np), cex.main=0.95)
  rows <- idx[((pg-1)*per_pg+1):min(pg*per_pg, nrow(idx)), ]
  hdr <- sprintf("%-34s %-11s %9s %9s %7s %5s %4s %-6s",
                 "name","adduct","theo m/z","data m/z","d_ppm","score","ana","src")
  lines_txt <- apply(rows, 1, function(r) sprintf("%-34.34s %-11s %9.4f %9.4f %7.1f %5.0f %4s %-6s",
       r["name"], r["adduct"], as.numeric(r["theo_mz"]), as.numeric(r["data_mz"]),
       as.numeric(r["delta_ppm"]), as.numeric(r["score"]),
       ifelse(r["in_analysis_set"]=="TRUE","Y","-"), sub("_.*","",r["source"])))
  text(0, 1, paste(c(hdr, "", lines_txt), collapse="\n"), adj=c(0,1), family="mono", cex=0.62)
}

# ---- one metabolite per page -----------------------------------------------
for (i in seq_len(nrow(mt))) {
  r <- mt[i, ]
  vals <- sp[as.character(r$feat_idx), ]
  pos <- vals[vals > 0]
  hi <- if (length(pos)) as.numeric(quantile(pos, IMG_CLIP_HI, na.rm=TRUE)) else 1
  par(mfrow=c(2,2), oma=c(0.5,0.5,5.2,0.5))
  for (sid in SIDS) panel(ps[[sid]], vals, hi, sid)
  iso <- if (r$isobaric_collision) "  [isobaric collision]" else ""
  ana <- if (r$in_analysis_set) "in analysis set" else "below freq-floor (recall grid only)"
  mtext(sprintf("%s   %s", r$name, r$adduct), outer=TRUE, line=3.4, cex=1.15, font=2)
  mtext(sprintf("HMDB %s  |  %s  |  %s", r$hmdb, r$formula, r$hmdb_name),
        outer=TRUE, line=2.3, cex=0.75)
  mtext(sprintf("theoretical m/z %.4f   |   measured %.4f   |   delta  %+.1f mDa  (%+.1f ppm)   |   score %.0f/100",
                r$theo_mz, r$data_mz, r$delta_mDa, r$delta_ppm, r$score),
        outer=TRUE, line=1.3, cex=0.85)
  mtext(sprintf("source: %s   |   %s%s   |   global p99.5 clip = %.3g",
                r$source, ana, iso, hi),
        outer=TRUE, line=0.4, cex=0.7, col="grey30")
  if (i %% 25 == 0) log_msg("  %d / %d pages", i, nrow(mt))
}
dev.off()
log_msg("DONE -> %s (%d metabolite pages + %d index)", OUT_PDF, nrow(mt), np)
