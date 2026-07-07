#!/usr/bin/env Rscript
# ============================================================================
# 01_teach_msi_to_jpg.R - Phase R1: deterministic MSI -> slide-JPG transform
#   from Bruker teach points (.mis) + stage log (poslog.txt).
#   raster -> stage -> flexImage -> original _small.jpg  (3 affines composed).
#   Orientation of the MSI grid vs raster is auto-resolved by which flip makes
#   the is_tissue mask land on DARK (organoid) jpg pixels (max contrast).
#
# Output: cache/register/teach_<sid>.rds  (B_msi_jpg 3x2, W,H, bbox, score...)
#         results/registration/teach_summary.csv
#         figures/registration/teach_msi2jpg_<slide>.pdf
# Usage : Rscript R/03_coarse_registration/01_teach_msi_to_jpg.R
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
suppressPackageStartupMessages({ library(Cardinal); library(jpeg) })

MSI_DIR <- file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"MSI")
REG_CACHE <- file.path(CACHE_DIR, "register"); REG_RES <- file.path(RES_DIR, "registration"); REG_FIG <- file.path(FIG_DIR, "registration")
for (d in c(REG_CACHE, REG_RES, REG_FIG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

slide_of <- function(sid) if (grepl("sl6A", sid)) "sl6A" else "sl4A"
slide_folder <- function(sl) if (sl == "sl6A") file.path(MSI_DIR, "06102026_AO_0h_sl6A") else file.path(MSI_DIR, "06102026_AO_20h_sl4A")
slide_jpg <- function(sl) if (sl == "sl6A") file.path(MSI_DIR, "06102026_AO_0h_sl6A_small.jpg") else file.path(MSI_DIR, "06102026_AO_20h_sl4A_small.jpg")
sec_code <- function(sid) sub(".*_sec", "", sid)   # AO_0h_sl6A_sec1a -> 1a

# find .mis / poslog for a section (sl6A uses _sec1a, sl4A uses _1a)
find_file <- function(folder, code, suffix) {
  f <- list.files(folder, pattern = sprintf("(_sec%s|_%s)%s$", code, code, suffix), full.names = TRUE)
  f <- f[!grepl("\\.bak$", f)]
  if (!length(f)) return(NA_character_); f[which.min(nchar(basename(f)))]
}

mse <- readRDS(TISSUE_MSE); pd <- as.data.frame(pixelData(mse))
SIDS <- levels(pixelData(mse)$sample_id); if (is.null(SIDS)) SIDS <- sort(unique(as.character(pd$sample_id)))

jpg_cache <- list()
get_jpg <- function(sl) {
  if (!is.null(jpg_cache[[sl]])) return(jpg_cache[[sl]])
  j <- jpeg::readJPEG(slide_jpg(sl)); g <- if (length(dim(j)) == 3) (j[,,1]+j[,,2]+j[,,3])/3 else j
  jpg_cache[[sl]] <<- list(gray = g, H = nrow(g), W = ncol(g)); jpg_cache[[sl]]
}
# darkness (organoid=high) at jpg pixel coords (ox=col, oy=row)
dark_at <- function(J, ox, oy) {
  ox <- pmin(pmax(round(ox), 1), J$W); oy <- pmin(pmax(round(oy), 1), J$H)
  1 - J$gray[cbind(oy, ox)]
}

CAND <- list(c(1,1), c(1,-1), c(-1,1), c(-1,-1))
place <- function(sub, bx, by, W, H, o) {   # MSI(x,y) -> flex coords for orientation o
  fx <- if (o[1] > 0) bx[1] + (sub$x - 1)/(W - 1)*(bx[2]-bx[1]) else bx[2] - (sub$x - 1)/(W - 1)*(bx[2]-bx[1])
  fy <- if (o[2] > 0) by[1] + (sub$y - 1)/(H - 1)*(by[2]-by[1]) else by[2] - (sub$y - 1)/(H - 1)*(by[2]-by[1])
  cbind(fx, fy)
}

# ---- Pass 1: collect per-section geometry + per-orientation contrast ----
secs <- list(); per_slide <- list()
for (sid in SIDS) {
  sl <- slide_of(sid); code <- sec_code(sid); folder <- slide_folder(sl)
  mis <- find_file(folder, code, "\\.mis")
  if (is.na(mis)) { cat(sprintf("[30b] %s: .mis not found (code %s)\n", sid, code)); next }
  M <- parse_mis(mis)
  if (nrow(M$orig) < 3 || nrow(M$area) < 2) { cat(sprintf("[30b] %s: insufficient teach/area\n", sid)); next }
  B_fo <- fit_affine(M$orig[, 1:2], M$orig[, 3:4]); rmse_fo <- affine_rmse(B_fo, M$orig[,1:2], M$orig[,3:4])
  bx <- range(M$area[, 1]); by <- range(M$area[, 2])
  sub <- pd[as.character(pd$sample_id) == sid, c("x","y","is_tissue")]; W <- max(sub$x); H <- max(sub$y)
  J <- get_jpg(sl)
  contrasts <- vapply(CAND, function(o) {
    oc <- apply_affine(B_fo, place(sub, bx, by, W, H, o)); dk <- dark_at(J, oc[,1], oc[,2])
    mean(dk[sub$is_tissue], na.rm = TRUE) - mean(dk[!sub$is_tissue], na.rm = TRUE)
  }, numeric(1))
  secs[[sid]] <- list(sid=sid, slide=sl, B_fo=B_fo, rmse_fo=rmse_fo, bx=bx, by=by, sub=sub, W=W, H=H, contrasts=contrasts)
  per_slide[[sl]] <- c(per_slide[[sl]], sid)
}

# ---- Choose ONE orientation per slide (max total contrast across its sections) ----
global_o <- list()
for (sl in names(per_slide)) {
  tot <- rowSums(vapply(per_slide[[sl]], function(s) secs[[s]]$contrasts, numeric(length(CAND))))
  global_o[[sl]] <- CAND[[which.max(tot)]]
  cat(sprintf("[30b] slide %s global orientation (%+d,%+d) | per-cand totals: %s\n",
              sl, global_o[[sl]][1], global_o[[sl]][2], paste(sprintf("%.3f", tot), collapse=" ")))
}

# ---- Pass 2: save transforms with the slide-global orientation ----
summ <- list()
for (sid in names(secs)) {
  S <- secs[[sid]]; o <- global_o[[S$slide]]
  oc <- apply_affine(S$B_fo, place(S$sub, S$bx, S$by, S$W, S$H, o))
  dk <- dark_at(get_jpg(S$slide), oc[,1], oc[,2])
  contrast <- mean(dk[S$sub$is_tissue], na.rm=TRUE) - mean(dk[!S$sub$is_tissue], na.rm=TRUE)
  B_msi_jpg <- fit_affine(cbind(S$sub$x, S$sub$y), oc); bbox <- apply(oc, 2, range)
  saveRDS(list(sid = sid, slide = S$slide, B_msi_jpg = B_msi_jpg, W = S$W, H = S$H,
               orientation = o, contrast = contrast, rmse_flex_orig = S$rmse_fo,
               jpg_bbox = bbox, jpg = slide_jpg(S$slide)),
          file.path(REG_CACHE, sprintf("teach_%s.rds", sid)))
  summ[[sid]] <- data.frame(sample_id = sid, slide = S$slide, W = S$W, H = S$H,
                            orient_x = o[1], orient_y = o[2], contrast = round(contrast, 4),
                            rmse_flex_orig = round(S$rmse_fo, 3), stringsAsFactors = FALSE)
  cat(sprintf("[30b] %s: orient(%+d,%+d) contrast=%.3f bbox x[%.0f-%.0f] y[%.0f-%.0f]\n",
              sid, o[1], o[2], contrast, bbox[1,1], bbox[2,1], bbox[1,2], bbox[2,2]))
}
summ_df <- do.call(rbind, summ); write.csv(summ_df, file.path(REG_RES, "teach_summary.csv"), row.names = FALSE)

# ---- QC PDF per slide: tissue mask overlaid on jpg ----
for (sl in names(per_slide)) {
  J <- get_jpg(sl)
  pdf(file.path(REG_FIG, sprintf("teach_msi2jpg_%s.pdf", sl)), width = 12, height = 7)
  # whole slide with all section footprints
  par(mar = c(1,1,3,1)); plot.new(); plot.window(c(1, J$W), c(J$H, 1), asp = 1)
  rasterImage(J$gray, 1, J$H, J$W, 1)
  pal <- grDevices::hcl.colors(length(per_slide[[sl]]), "Dark 3")
  for (i in seq_along(per_slide[[sl]])) {
    sid <- per_slide[[sl]][i]; r <- readRDS(file.path(REG_CACHE, sprintf("teach_%s.rds", sid)))
    sub <- pd[as.character(pd$sample_id) == sid, c("x","y","is_tissue")]
    oc <- apply_affine(r$B_msi_jpg, cbind(sub$x, sub$y))
    points(oc[sub$is_tissue,1], oc[sub$is_tissue,2], pch = 15, cex = 0.12, col = adjustcolor(pal[i], 0.6))
    cen <- colMeans(oc); text(cen[1], cen[2], sec_code(sid), col = "red", cex = 0.8, font = 2)
  }
  title(sprintf("%s - all MSI section tissue masks on slide JPG (teach-point registration)", sl))
  # per-section zoomed panels
  par(mfrow = c(2,3), mar = c(1,1,2,1))
  for (sid in per_slide[[sl]]) {
    r <- readRDS(file.path(REG_CACHE, sprintf("teach_%s.rds", sid)))
    sub <- pd[as.character(pd$sample_id) == sid, c("x","y","is_tissue")]
    oc <- apply_affine(r$B_msi_jpg, cbind(sub$x, sub$y))
    bx <- range(oc[,1]) + c(-30,30); by <- range(oc[,2]) + c(-30,30)
    bx <- pmin(pmax(bx,1),J$W); by <- pmin(pmax(by,1),J$H)
    cols <- max(1,floor(bx[1])):min(J$W,ceiling(bx[2])); rows <- max(1,floor(by[1])):min(J$H,ceiling(by[2]))
    plot.new(); plot.window(range(cols), rev(range(rows)), asp = 1)
    rasterImage(J$gray[rows, cols], min(cols), max(rows), max(cols), min(rows))
    points(oc[sub$is_tissue,1], oc[sub$is_tissue,2], pch = 15, cex = 0.3, col = adjustcolor("red", 0.45))
    title(sprintf("%s  contrast=%.2f orient(%+d,%+d)", sec_code(sid), r$contrast, r$orientation[1], r$orientation[2]), cex.main = 0.9)
  }
  dev.off(); cat(sprintf("[30b] wrote QC: teach_msi2jpg_%s.pdf\n", sl))
}
cat("\n[30b] DONE\n"); print(summ_df, row.names = FALSE)
