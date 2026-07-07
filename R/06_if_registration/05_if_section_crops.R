#!/usr/bin/env Rscript
# ============================================================================
# 05_if_section_crops.R - fresh start (the NCC/overview->BF alignment was wrong;
#   only the step0 DAPI foundation is kept). Detect & CROP the 6 organoid
#   sections from each whole-slide OVERVIEW DAPI and report them, so we can
#   rebuild registration on clean per-section crops.
#
#   Detection: DAPI texture (nuclei) -> blur -> top-6 density blobs -> boxes.
#   The boxes are PROVISIONAL: a clean annotatable overview PDF is also written
#   per slide so the user can draw the 6 rectangles by hand if any box is off.
#
# Input : cache/register_if/ovthumb_<slide>.rds, hrthumb_<sid_if>.rds
# Output: figures/if_registration/section_crops_<slide>.pdf  (boxes + 6 crops + hi-res)
#         figures/if_registration/overview_<slide>_annotate.pdf  (clean, for drawing)
#         results/if_registration/section_boxes_<slide>.csv  (editable box coords, native px)
# Usage : Rscript R/06_if_registration/05_if_section_crops.R [all | slide ...]   (default = all slides)
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/if_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register_if.R"))
suppressPackageStartupMessages({ library(EBImage) })

args   <- commandArgs(trailingOnly = TRUE)
slides <- if (length(args) == 0 || identical(args, "all")) IF_SLIDES$slide else args
SEC    <- if_sections()
ANNOT_BF_F <- 4L   # brightfield annotate downsample (block-min, ~7.3 um/px)

# detect the n strongest organoid (DAPI-texture) blobs -> bounding boxes (m16 px)
detect_organoids <- function(m16, n = 6, blur = 12L, qd = 0.55, min_area = 25L) {
  tex <- texture_map(m16, 2L)
  bk  <- makeBrush(blur*2+1, "disc"); bl <- filter2(tex, bk/sum(bk), boundary = "replicate")
  bm  <- closing(bl > quantile(bl, qd), makeBrush(9, "disc")); bm <- fillHull(bm)
  lab <- bwlabel(bm); tb <- table(lab[lab > 0]); comps <- as.integer(names(tb))
  area <- as.integer(tb[as.character(comps)]); comps <- comps[area >= min_area]
  score <- sapply(comps, function(l) sum(tex[lab == l]))
  keep  <- comps[order(score, decreasing = TRUE)][seq_len(min(n, length(comps)))]
  bx <- lapply(keep, function(l) { idx <- which(lab == l, arr.ind = TRUE)
    c(r0 = min(idx[,1]), r1 = max(idx[,1]), c0 = min(idx[,2]), c1 = max(idx[,2])) })
  do.call(rbind, bx)
}
# order boxes column-major (down each column, columns left->right) -> sec1..6
order_colmajor <- function(B, ncol_grid = 3L) {
  cx <- (B[,"c0"]+B[,"c1"])/2; cy <- (B[,"r0"]+B[,"r1"])/2
  col <- cut(rank(cx, ties.method = "first"), breaks = ncol_grid, labels = FALSE)
  order(col, cy)
}

# box source: "ncc" = step-1 NCC field boxes (validated on organoids), "detect" = fresh detection
BOX_SOURCE <- Sys.getenv("IF_BOX_SOURCE", "ncc")

for (sl in slides) {
  ov <- readRDS(file.path(IF_CACHE, sprintf("ovthumb_%s.rds", sl)))
  m16 <- ov$m16; F16 <- ov$F16
  secs <- SEC[SEC$slide == sl, ]
  if (BOX_SOURCE == "ncc" && all(file.exists(file.path(IF_CACHE, sprintf("locate_%s.rds", secs$sid_if))))) {
    # use the validated NCC field boxes, in sec1..6 order
    B <- do.call(rbind, lapply(secs$sid_if, function(s){ b <- readRDS(file.path(IF_CACHE, sprintf("locate_%s.rds", s)))$box
      c(r0 = b["r0"], r1 = b["r0"]+b["th"]-1, c0 = b["c0"], c1 = b["c0"]+b["tw"]-1) }))
    colnames(B) <- c("r0","r1","c0","c1")
  } else {
    B <- detect_organoids(m16); B <- B[order_colmajor(B), , drop = FALSE]
  }
  pad <- function(b, f = 0.05) { dr <- (b["r1"]-b["r0"])*f; dc <- (b["c1"]-b["c0"])*f
    c(r0 = max(1, floor(b["r0"]-dr)), r1 = min(nrow(m16), ceiling(b["r1"]+dr)),
      c0 = max(1, floor(b["c0"]-dc)), c1 = min(ncol(m16), ceiling(b["c1"]+dc))) }
  Bp <- t(apply(B, 1, pad))

  # editable CSV in overview-native px
  csv <- data.frame(section = seq_len(nrow(Bp)),
                    x0 = round((Bp[,"c0"]-1)*F16+1), y0 = round((Bp[,"r0"]-1)*F16+1),
                    x1 = round(Bp[,"c1"]*F16),       y1 = round(Bp[,"r1"]*F16))
  write.csv(csv, file.path(IF_RES, sprintf("section_boxes_%s.csv", sl)), row.names = FALSE)

  # ---- report: overview + boxes, then the 6 crops (+ paired hi-res DAPI) ----
  pdf(file.path(IF_FIG, sprintf("section_crops_%s.pdf", sl)), width = 12, height = 7)
  # page 1: whole overview + numbered boxes
  par(mar = c(1,1,3,1))
  plot.new(); plot.window(c(1, ncol(m16)), c(nrow(m16), 1), asp = 1)
  rasterImage(as.raster(norm01(m16)), 1, nrow(m16), ncol(m16), 1)
  for (j in seq_len(nrow(Bp))) { b <- Bp[j, ]; rect(b["c0"], b["r0"], b["c1"], b["r1"], border = "red", lwd = 2.5)
    text((b["c0"]+b["c1"])/2, (b["r0"]+b["r1"])/2, j, col = "yellow", cex = 1.6, font = 2) }
  title(sprintf("%s overview DAPI - 6 organoid sections (boxes = hi-res DAPI cross-corr fields; editable CSV)", sl), cex.main = 1.0)
  # page 2: the 6 overview crops
  par(mfrow = c(2, 3), mar = c(1,1,2.5,1))
  for (j in seq_len(nrow(Bp))) { b <- Bp[j, ]; cr <- m16[b["r0"]:b["r1"], b["c0"]:b["c1"]]
    plot.new(); plot.window(c(1, ncol(cr)), c(nrow(cr), 1), asp = 1); rasterImage(as.raster(norm01(cr)), 1, nrow(cr), ncol(cr), 1)
    title(sprintf("section %d  (overview DAPI crop)", j), cex.main = 1) }
  # page 3: paired hi-res DAPI sections (step0 foundation), sec1..6 order
  secs <- SEC[SEC$slide == sl, ]
  par(mfrow = c(2, 3), mar = c(1,1,2.5,1))
  for (j in seq_len(nrow(secs))) { hf <- file.path(IF_CACHE, sprintf("hrthumb_%s.rds", secs$sid_if[j]))
    if (!file.exists(hf)) { plot.new(); next }; h <- readRDS(hf)
    plot.new(); plot.window(c(1, ncol(h$m5)), c(nrow(h$m5), 1), asp = 1); rasterImage(as.raster(norm01(h$m5)), 1, nrow(h$m5), ncol(h$m5), 1)
    title(sprintf("hi-res file sec%d DAPI", secs$secn[j]), cex.main = 1) }
  dev.off()

  # ---- clean annotatable BRIGHTFIELD A-slide (no boxes) for hand-drawn rects --
  # user draws the 6 organoid rectangles on the MSI A-slide brightfield (the side
  # the MSI data sits in). Contrast-stretched thumbnail for visibility.
  msi_slide <- IF_SLIDES$msi_slide[IF_SLIDES$slide == sl]
  msi_bf <- IF_SLIDES$msi_bf[IF_SLIDES$slide == sl]
  bmf <- file.path(IF_CACHE, sprintf("bfmin_%s_F%d.rds", msi_slide, ANNOT_BF_F))
  bm  <- if (file.exists(bmf)) readRDS(bmf) else { x <- nd2_block_min(msi_bf, F = ANNOT_BF_F)$m; saveRDS(x, bmf); x }
  # block-min (keeps faint tissue edges) + high-pass (removes the MALDI grid/illumination)
  v   <- norm01(bm); kb <- EBImage::makeBrush(51, "disc")
  bfd <- norm01(v - EBImage::filter2(v, kb/sum(kb), boundary = "replicate"), 0.01, 0.99)
  pdf(file.path(IF_FIG, sprintf("annotate_brightfield_%s.pdf", msi_slide)),
      width = 16, height = 16*nrow(bfd)/ncol(bfd))
  par(mar = c(0,0,0,0)); plot.new(); plot.window(c(1, ncol(bfd)), c(nrow(bfd), 1), asp = 1)
  rasterImage(as.raster(bfd), 1, nrow(bfd), ncol(bfd), 1)
  dev.off()
  cat(sprintf("[94] %s: %d sections -> section_crops_%s.pdf (+ brightfield %s annotate PDF, boxes CSV)\n",
              sl, nrow(Bp), sl, msi_slide))
}
cat("[94] DONE\n")
