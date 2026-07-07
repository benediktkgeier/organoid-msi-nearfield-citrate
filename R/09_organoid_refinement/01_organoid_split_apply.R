#!/usr/bin/env Rscript
# ============================================================================
# 01_organoid_split_apply.R - Phase 09 step 1: MORPHOLOGICAL SPLIT of clustered
#   organoids, applied RIGHT AFTER segmentation
#   (R/08_organoid_gradient_survey/01_segment_organoids.R).
# ----------------------------------------------------------------------------
# Plain connected-component labelling in R/08_organoid_gradient_survey/
# 01_segment_organoids.R merges touching organoids into one instance. The user
# marks the cuts by drawing GREEN freehand (/Ink) lines through the clustered
# organoids on the R/10_apical_annotation/01_apical_annotate.R annotation PDF
# (any copy). This
# script reads those green lines, lays each as a 1px-dilated barrier over the
# instance footprint, removes it, re-runs connected components, and splits the
# merged instance -- mirroring the proven gastric-gland workflow in
# D:/R/PeakMe/PeakMe_GCPL/phaseR_apply_splits.R.
#
# ACCUMULATION (multi-round): cuts are drawn over several rounds, each on the
#   freshly re-rendered canvas (which shows the prior splits). So the green
#   strokes live in DIFFERENT PDFs and a re-render wipes a PDF's strokes. To make
#   this durable, every stroke is mapped to MSI coords and accumulated into a
#   persistent store cache/organoid_split_strokes.rds (deduped by geometry). The
#   store -- not any PDF -- is the source of truth; re-rendering never loses cuts.
#   Each run reads the input PDF(s), adds any NEW strokes to the store, then
#   re-applies ALL accumulated strokes to the untouched originals.
#
# ID rule (preserve + append): an instance untouched by a green line keeps its
#   id; a split instance keeps its LARGEST piece as the original id and new
#   pieces get appended ids. Barrier pixels stay with the kept piece (no px lost).
#
# Coordinate map (PDF point -> MSI grid): green-line points share the
#   R/10_apical_annotation/01_apical_annotate.R sidecar's PDF-point frame. We
#   invert that render exactly by REPLAYING its
#   plot setup (14x10 in, plot.window(c(1,cw),c(ch,1),asp=1)) and grconvertX/Y
#   ("ndc","user") to recover crop-local native px, then native -> MSI via the
#   inverse of B_msi_nd2. Needs NO organoid correspondences (works on 1-organoid
#   pages), and is render-independent (same physical cut -> same MSI coords).
#
# Storage: split result -> SEPARATE file cache/instances_split_<sid>.rds (same
#   columns as instances_<sid>.rds); originals untouched.
#   R/10_apical_annotation/01_apical_annotate.R prefers the split file when
#   present. Idempotent: always re-derived from the untouched originals.
#
# Usage:
#   Rscript R/09_organoid_refinement/01_organoid_split_apply.R [pdf1 pdf2 ...]
#     default inputs: organoid_apical_annotation_BG.pdf + organoid_apical_annotation.pdf
#   Rscript R/09_organoid_refinement/01_organoid_split_apply.R --reset [pdf ...]   # clear the store first
#
# Out: cache/instances_split_<sid>.rds (sids with >=1 split),
#      cache/organoid_split_strokes.rds (accumulated stroke store),
#      figures/annotation/organoid_split_curation.pdf,
#      results/annotation/organoid_split_report.csv
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
suppressPackageStartupMessages({ library(png); library(EBImage) })

REG_CACHE <- file.path(CACHE_DIR, "register")
CROP_DIR  <- file.path(FIG_DIR, "registration", "crops")
ANNOT_FIG <- file.path(FIG_DIR, "annotation")
ANNOT_RES <- file.path(RES_DIR, "annotation")
dir.create(ANNOT_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(ANNOT_RES, showWarnings = FALSE, recursive = TRUE)

ANNOT_PY     <- file.path(PROJECT_ROOT, "py", "extract_pdf_lines.py")
ANNOT_PY_BIN <- "C:/Users/bened/.virtualenvs/r-reticulate/Scripts/python.exe"
STORE_RDS    <- file.path(CACHE_DIR, "organoid_split_strokes.rds")
DROPPED_RDS  <- file.path(CACHE_DIR, "organoid_split_dropped.rds")  # keys to NEVER re-add
SIDECAR_RDS  <- file.path(CACHE_DIR, "organoid_apical_label_positions.rds")
SCALE_UM <- 100; PX_UM <- MSI_PIXEL_UM
PAGE_W_IN <- 14; PAGE_H_IN <- 10
PAGE_W_PT <- PAGE_W_IN * 72; PAGE_H_PT <- PAGE_H_IN * 72
ORG_PAL <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628","#f781bf",
             "#1b9e77","#d95f02","#7570b3","#66a61e","#e7298a")
org_colors <- function(ids) setNames(ORG_PAL[((seq_along(ids)-1) %% length(ORG_PAL))+1], as.character(ids))
group_of <- function(sid) "incubated_5min"
log_msg <- function(...) message(sprintf("[41] %s", sprintf(...)))

args <- commandArgs(trailingOnly = TRUE)
RESET <- "--reset" %in% args; args <- args[args != "--reset"]
# --apply-only: re-derive instances_split from the EXISTING stroke store without
# reading any PDF (used after a targeted store edit, e.g. dropping a section's
# superseded strokes so it can be re-split differently).
APPLY_ONLY <- "--apply-only" %in% args; args <- args[args != "--apply-only"]
# --any-color: accept ink of ANY colour as a split (default: GREEN only). Use when
# a split line was drawn in a non-green colour; pass ONLY the relevant PDF so other
# coloured marks (e.g. red pointers) aren't swept in.
ANY_COLOR <- "--any-color" %in% args; args <- args[args != "--any-color"]
default_pdfs <- c(file.path(ANNOT_FIG, "organoid_apical_annotation_BG.pdf"),
                  file.path(ANNOT_FIG, "organoid_apical_annotation.pdf"))
INPUT_PDFS <- if (APPLY_ONLY) character(0) else if (length(args)) args else default_pdfs[file.exists(default_pdfs)]
INPUT_PDFS <- INPUT_PDFS[file.exists(INPUT_PDFS)]
stopifnot(file.exists(ANNOT_PY), file.exists(SIDECAR_RDS))

# ---- helpers ---------------------------------------------------------------
is_green <- function(cs) {
  v <- suppressWarnings(as.numeric(strsplit(cs, ",")[[1]]))
  length(v) == 3 && !any(is.na(v)) && v[2] > 0.4 && v[1] < 0.3 && v[3] < 0.3
}
parse_path <- function(s) {
  pts <- strsplit(s, ";")[[1]]
  do.call(rbind, lapply(pts, function(p) as.numeric(strsplit(p, ",")[[1]])))
}
instance_outlines <- function(M, B, cx0, cy0) {       # same idiom as R/10_apical_annotation/01_apical_annotate.R
  M <- EBImage::fillHull(EBImage::closing(M > 0, EBImage::makeBrush(3, "disc")))
  out <- list()
  for (co in EBImage::ocontour(EBImage::bwlabel(M))) {
    if (nrow(co) < 6) next
    p <- apply_affine(B, cbind(co[, 1] + 0.5, co[, 2] + 0.5))
    out[[length(out) + 1L]] <- list(x = c(p[, 1] - cx0 + 1, p[1, 1] - cx0 + 1),
                                    y = c(p[, 2] - cy0 + 1, p[1, 2] - cy0 + 1))
  }
  out
}
label_halo <- function(x, y, txt, col, cex = 1.4) {
  off <- 0.6 * (cex / 1.5) * strwidth("0", cex = cex) / 4 + 0.4
  for (dx in c(-1, 1)) for (dy in c(-1, 1))
    text(x + dx * off, y + dy * off, txt, col = "white", cex = cex, font = 2, xpd = NA)
  text(x, y, txt, col = col, cex = cex, font = 2, xpd = NA)
}
raster_polyline <- function(v, W, H) {
  if (nrow(v) < 2) return(matrix(integer(0), 0, 2))
  out <- list()
  for (k in seq_len(nrow(v) - 1)) {
    n  <- max(ceiling(max(abs(v[k+1,1]-v[k,1]), abs(v[k+1,2]-v[k,2]))), 1L)
    xx <- round(seq(v[k,1], v[k+1,1], length.out = n + 1))
    yy <- round(seq(v[k,2], v[k+1,2], length.out = n + 1))
    out[[k]] <- cbind(xx, yy)
  }
  p <- do.call(rbind, out)
  p[p[,1] >= 1 & p[,1] <= W & p[,2] >= 1 & p[,2] <= H, , drop = FALSE]
}
native_to_msi <- function(B, N) {                     # invert apply_affine
  M2 <- matrix(c(B[1,1], B[1,2], B[2,1], B[2,2]), nrow = 2)
  sweep(N, 2, B[3, ], "-") %*% t(solve(M2))
}
surface_of <- function(T) {                           # 4-neighbour surface (matches R/08_organoid_gradient_survey/01_segment_organoids.R)
  W <- nrow(T); H <- ncol(T)
  pad <- matrix(FALSE, W + 2L, H + 2L); pad[2:(W + 1L), 2:(H + 1L)] <- T
  n4 <- pad[1:W, 2:(H + 1L)] + pad[3:(W + 2L), 2:(H + 1L)] +
        pad[2:(W + 1L), 1:H]  + pad[2:(W + 1L), 3:(H + 2L)]
  T & (n4 < 4)
}

# ---- page -> sid via sidecar -----------------------------------------------
sc <- readRDS(SIDECAR_RDS)
PG_H <- attr(sc, "page_h_pt"); if (is.null(PG_H)) PG_H <- PAGE_H_PT
PG_W <- attr(sc, "page_w_pt"); if (is.null(PG_W)) PG_W <- PAGE_W_PT
sid_of_page <- function(pg) { u <- unique(sc$sid[sc$page == pg]); if (length(u)) u[1] else NA_character_ }

# cache per-sid registration + crop dims (read once)
reg_cache <- new.env()
get_reg <- function(sid) {
  if (!is.null(reg_cache[[sid]])) return(reg_cache[[sid]])
  xf <- file.path(REG_CACHE, sprintf("nd2final_%s.rds", sid))
  pf <- file.path(CROP_DIR, sprintf("optical_%s.png", sid))
  if (!file.exists(xf) || !file.exists(pf)) return(NULL)
  Xr <- readRDS(xf); om <- png::readPNG(pf); if (length(dim(om)) == 3) om <- om[, , 1]
  r <- list(B = Xr$B_msi_nd2, smn = Xr$scale_msi_nd2, cx0 = Xr$crop[1], cy0 = Xr$crop[2],
            cw = ncol(om), ch = nrow(om), om = om)
  reg_cache[[sid]] <- r; r
}

# ---- extract candidate strokes from the input PDFs -------------------------
candidates <- list()
for (pdf in INPUT_PDFS) {
  tsv <- tempfile(fileext = ".tsv")
  out <- system2(ANNOT_PY_BIN, c(shQuote(ANNOT_PY), shQuote(pdf), shQuote(tsv)),
                 stdout = TRUE, stderr = TRUE)
  if (!file.exists(tsv)) { log_msg("WARN extract failed for %s", basename(pdf)); next }
  ann <- read.table(tsv, sep = "\t", header = TRUE, quote = "", comment.char = "",
                    stringsAsFactors = FALSE, fill = TRUE); unlink(tsv)
  keep_col <- if (ANY_COLOR) rep(TRUE, nrow(ann)) else vapply(ann$color, is_green, logical(1))
  ann <- ann[ann$subtype %in% c("/Ink","/Line","/PolyLine") & nzchar(ann$path) & keep_col, , drop = FALSE]
  if (!nrow(ann)) next
  for (i in seq_len(nrow(ann))) {
    sid <- sid_of_page(ann$page[i]); if (is.na(sid)) next
    pts <- parse_path(ann$path[i])
    key <- paste0(sid, "|", ann$page[i], "|", paste(sprintf("%.1f", c(t(pts))), collapse = ","))
    candidates[[length(candidates) + 1L]] <- list(sid = sid, page = ann$page[i],
                                                  src = basename(pdf), pts = pts, key = key)
  }
  log_msg("%s: %d green stroke(s)", basename(pdf), sum(!is.na(vapply(seq_len(nrow(ann)),
          function(i) sid_of_page(ann$page[i]), character(1)))))
}

# ---- map candidates -> MSI (replay R/10_apical_annotation/01_apical_annotate.R render), on a throwaway device -----
if (length(candidates)) {
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp, width = PAGE_W_IN, height = PAGE_H_IN)
  by_sid <- split(seq_along(candidates), vapply(candidates, `[[`, character(1), "sid"))
  for (sid in names(by_sid)) {
    rg <- get_reg(sid); if (is.null(rg)) { log_msg("%s: SKIP (missing registration)", sid); next }
    par(mar = c(1, 1, 3, 1)); plot.new(); plot.window(c(1, rg$cw), c(rg$ch, 1), asp = 1)
    for (idx in by_sid[[sid]]) {
      m <- candidates[[idx]]$pts
      ndc_x <- m[, 1] / PG_W; ndc_y <- 1 - m[, 2] / PG_H
      ux <- grconvertX(ndc_x, "ndc", "user"); uy <- grconvertY(ndc_y, "ndc", "user")
      candidates[[idx]]$msi <- native_to_msi(rg$B, cbind(ux + rg$cx0 - 1, uy + rg$cy0 - 1))
    }
  }
  dev.off(); unlink(tmp)
}

# ---- accumulate into the persistent store (dedup by geometry key) ----------
store <- if (!RESET && file.exists(STORE_RDS)) readRDS(STORE_RDS) else list()
# keys explicitly dropped by store surgery -- never re-added even if a source PDF
# still carries the stroke (e.g. a section re-split differently after a merge)
DROPPED <- if (!RESET && file.exists(DROPPED_RDS)) readRDS(DROPPED_RDS) else character(0)
if (RESET) { if (file.exists(DROPPED_RDS)) file.remove(DROPPED_RDS) }
have  <- vapply(store, `[[`, character(1), "key")
n_new <- 0L; n_skip_dropped <- 0L
for (cd in candidates) {
  if (is.null(cd$msi)) next
  if (cd$key %in% DROPPED) { n_skip_dropped <- n_skip_dropped + 1L; next }
  if (cd$key %in% have) next
  store[[length(store) + 1L]] <- list(key = cd$key, sid = cd$sid, page = cd$page,
                                      src = cd$src, msi = cd$msi)
  have <- c(have, cd$key); n_new <- n_new + 1L
}
saveRDS(store, STORE_RDS)
log_msg("%d stroke(s) in store (%d new this run, %d skipped as dropped, reset=%s)",
        length(store), n_new, n_skip_dropped, RESET)
if (!length(store)) stop("no split strokes accumulated; nothing to do")

# ============================================================================
# apply ALL stored strokes to the untouched originals, per affected section
# ============================================================================
CURATION_PDF <- file.path(ANNOT_FIG, "organoid_split_curation.pdf")
pdf(CURATION_PDF, width = PAGE_W_IN, height = PAGE_H_IN)
report <- list(); duds <- list(); n_written <- 0L
store_by_sid <- split(store, vapply(store, `[[`, character(1), "sid"))

for (sid in names(store_by_sid)) {
  rg  <- get_reg(sid); inf <- file.path(CACHE_DIR, sprintf("instances_%s.rds", sid))
  if (is.null(rg) || !file.exists(inf)) { log_msg("%s: SKIP (missing inputs)", sid); next }
  B <- rg$B; cx0 <- rg$cx0; cy0 <- rg$cy0; cw <- rg$cw; ch <- rg$ch; smn <- rg$smn; om <- rg$om
  strokes_msi <- lapply(store_by_sid[[sid]], `[[`, "msi")

  inst <- readRDS(inf); pos <- inst$instance > 0
  W <- max(inst$x); H <- max(inst$y)
  lab <- matrix(0L, W, H); lab[cbind(inst$x[pos], inst$y[pos])] <- as.integer(inst$instance[pos])
  ids0 <- sort(unique(as.integer(inst$instance[pos])))

  # dilated barrier from all stored strokes of this section
  bar <- matrix(FALSE, W, H)
  for (v in strokes_msi) { pc <- raster_polyline(v, W, H); if (nrow(pc)) bar[pc] <- TRUE }
  bar <- as.array(EBImage::dilate(EBImage::as.Image(bar * 1), EBImage::makeBrush(3, "box"))) > 0.5

  labN <- lab; next_id <- if (length(ids0)) max(ids0) + 1L else 1L
  split_detail <- integer(0)
  for (k in ids0) {
    foot <- lab == k
    comp <- as.array(EBImage::bwlabel(EBImage::as.Image(((foot & !bar)) * 1)))
    npc  <- max(comp); if (npc <= 1) next
    keep <- which.max(tabulate(comp[comp > 0]))
    for (cc in setdiff(seq_len(npc), keep)) { labN[foot & comp == cc] <- next_id; next_id <- next_id + 1L }
    split_detail <- c(split_detail, setNames(npc, as.character(k)))
  }
  ids1 <- sort(unique(as.integer(labN[labN > 0]))); did_split <- length(ids1) > length(ids0)

  # per-stroke diagnostic: does each stroke ALONE disconnect an instance? a stroke
  # that doesn't (clipped an edge / didn't cross) is flagged so it can be redrawn
  for (si in seq_along(store_by_sid[[sid]])) {
    bs <- matrix(FALSE, W, H); pc <- raster_polyline(strokes_msi[[si]], W, H)
    if (nrow(pc)) bs[pc] <- TRUE
    bs <- as.array(EBImage::dilate(EBImage::as.Image(bs * 1), EBImage::makeBrush(3, "box"))) > 0.5
    eff <- any(vapply(ids0, function(k) { f <- lab == k; if (!any(f & bs)) return(FALSE)
      max(as.array(EBImage::bwlabel(EBImage::as.Image(((f & !bs)) * 1)))) > 1 }, logical(1)))
    if (!eff) duds[[length(duds) + 1L]] <- data.frame(sid = sid,
      page = store_by_sid[[sid]][[si]]$page, src = store_by_sid[[sid]][[si]]$src,
      stringsAsFactors = FALSE)
  }

  report[[length(report) + 1L]] <- data.frame(
    sid = sid, group = group_of(sid), n_strokes = length(strokes_msi),
    n_before = length(ids0), n_after = length(ids1),
    split_ids = paste(names(split_detail), collapse = ";"),
    pieces = paste(unname(split_detail), collapse = ";"), stringsAsFactors = FALSE)

  if (did_split) {
    T <- labN > 0; surf <- surface_of(T)
    out <- inst
    out$instance <- as.integer(labN[cbind(out$x, out$y)]); out$instance[is.na(out$instance)] <- 0L
    out$is_surface <- as.logical(surf[cbind(out$x, out$y)]); out$is_surface[is.na(out$is_surface)] <- FALSE
    saveRDS(out, file.path(CACHE_DIR, sprintf("instances_split_%s.rds", sid)))
    n_written <- n_written + 1L
    log_msg("%s: %d -> %d organoid(s); wrote instances_split_%s.rds", sid, length(ids0), length(ids1), sid)
  } else log_msg("%s: %d stroke(s) but NO split", sid, length(strokes_msi))

  # curation page: BF + original (grey) + new (coloured + id) + cut lines (green)
  par(mar = c(1, 1, 3, 1)); plot.new(); plot.window(c(1, cw), c(ch, 1), asp = 1)
  rasterImage(om / max(om), 1, ch, cw, 1, interpolate = TRUE)
  for (k in ids0) for (poly in instance_outlines(lab == k, B, cx0, cy0))
    polygon(poly$x, poly$y, border = "grey70", lwd = 1.0)
  oc <- org_colors(ids1)
  for (k in ids1) {
    kc <- oc[as.character(k)]
    for (poly in instance_outlines(labN == k, B, cx0, cy0)) polygon(poly$x, poly$y, border = kc, lwd = 1.8)
    cen <- colMeans(which(labN == k, arr.ind = TRUE)); pc <- apply_affine(B, matrix(cen, nrow = 1))
    label_halo(pc[1] - cx0 + 1, pc[2] - cy0 + 1, as.character(k), kc, cex = 1.4)
  }
  for (v in strokes_msi) { p <- apply_affine(B, v); lines(p[, 1] - cx0 + 1, p[, 2] - cy0 + 1, col = "#16a020", lwd = 2.2) }
  title(sprintf("%s  (%s)   split %d -> %d organoid%s%s", sub("AO_", "", sid), group_of(sid),
                length(ids0), length(ids1), if (length(ids1) == 1) "" else "s",
                if (!did_split) "  [no split]" else ""),
        cex.main = 1.15, col.main = if (did_split) "#16a020" else "grey40", font.main = 2)
  umpx <- PX_UM / smn; bp <- SCALE_UM / umpx
  segments(cw - bp - 10, ch - 14, cw - 10, ch - 14, lwd = 3, col = "black")
  text(cw - bp / 2 - 10, ch - 28, sprintf("%d um", SCALE_UM), cex = 0.7)
}
dev.off()

rep_df <- if (length(report)) do.call(rbind, report) else
  data.frame(sid=character(), group=character(), n_strokes=integer(), n_before=integer(),
             n_after=integer(), split_ids=character(), pieces=character())
write.csv(rep_df, file.path(ANNOT_RES, "organoid_split_report.csv"), row.names = FALSE)

cat("\n[41] ===== SPLIT REPORT =====\n"); print(rep_df, row.names = FALSE)
if (length(duds)) {
  dud_df <- do.call(rbind, duds)
  write.csv(dud_df, file.path(ANNOT_RES, "organoid_split_nonsplitting_strokes.csv"), row.names = FALSE)
  cat(sprintf("\n[41] %d NON-SPLITTING stroke(s) (drawn but did not disconnect any organoid -- redundant\n      redraw of an existing cut, or did not cross the organoid; redraw across it if intended):\n", nrow(dud_df)))
  print(dud_df, row.names = FALSE)
}
log_msg("wrote %d instances_split_*.rds; curation -> %s", n_written, CURATION_PDF)
log_msg("NEXT: re-run R/102 to regenerate the annotation PDF on the split organoids.")
