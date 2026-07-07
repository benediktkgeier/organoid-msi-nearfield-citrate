#!/usr/bin/env Rscript
# ============================================================================
# 03_apical_consensus.R - consensus of TWO independent apical-orientation
#   annotations (in / out / mixed) made on the SAME island-cleanup canvas.
# ----------------------------------------------------------------------------
# Two scientists each placed Adobe in/out/mixed comments on the identical
# 21-page island-cleanup PDF geometry, so both share the centroid sidecar from
# R/09_organoid_refinement/02b_island_centroid_sidecar.R. We parse each comment
# (python extractor -> classify keyword -> nearest organoid centroid <= 60pt),
# join per (sid, instance), and render a per-section
# report so the two can SEE where they agree vs disagree and converge.
#
#   Annotator A ("orig"): figures/annotation/organoid_island_cleanup.pdf
#                         (the file the apical_citrate_dha report was built from)
#   Annotator B ("JM")  : figures/annotation/organoid_island_cleanup_INOUT_JM_0623.pdf
#
# Outline colour by status (user-chosen):
#   cyan    = both annotated, SAME class        (agree)
#   magenta = both annotated, DIFFERENT class   (disagree)
#   yellow  = only ONE annotator scored it      (call still missing)
#   grey    = neither annotated                 (thin)
#
# In : the two PDFs above; cache/organoid_island_label_positions_centroid.rds
#      (R/09_organoid_refinement/02b_island_centroid_sidecar.R);
#      cache/register/nd2final_<sid>.rds; cache/instances_{final,clean,
#      split,}_<sid>.rds; figures/registration/crops/optical_<sid>.png
# Out: results/annotation/apical_consensus_per_organoid.csv
#      figures/annotation/apical_consensus_report.pdf
# Usage: Rscript R/10_apical_annotation/03_apical_consensus.R
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
suppressPackageStartupMessages({ library(png); library(EBImage) })

# ---- inputs ----------------------------------------------------------------
REG_CACHE <- file.path(CACHE_DIR, "register")
CROP_DIR  <- file.path(FIG_DIR, "registration", "crops")
ANNOT_FIG <- file.path(FIG_DIR, "annotation")
ANNOT_RES <- file.path(RES_DIR, "annotation")
SCALE_UM  <- 100; PX_UM <- MSI_PIXEL_UM
PAGE_W_IN <- 14; PAGE_H_IN <- 10
PAGE_W_PT <- PAGE_W_IN * 72; PAGE_H_PT <- PAGE_H_IN * 72

PDF_ORIG <- file.path(ANNOT_FIG, "organoid_island_cleanup.pdf")
PDF_JM   <- file.path(ANNOT_FIG, "organoid_island_cleanup_INOUT_JM_0623.pdf")
SIDECAR  <- file.path(CACHE_DIR, "organoid_island_label_positions_centroid.rds")
ANNOT_PY     <- file.path(PROJECT_ROOT, "py", "extract_pdf_annots.py")
ANNOT_PY_BIN <- "C:/Users/bened/.virtualenvs/r-reticulate/Scripts/python.exe"
ACCEPT_TYPES <- c("/FreeText", "/Text", "/Square", "/Circle", "/Highlight",
                  "/StrikeOut", "/Underline", "/Squiggly", "/Popup")
MATCH_RADIUS_PT <- 60          # generous: comment-on-body -> centroid (R/42b)

OUT_CSV <- file.path(ANNOT_RES, "apical_consensus_per_organoid.csv")
OUT_PDF <- file.path(ANNOT_FIG, "apical_consensus_report.pdf")
stopifnot(file.exists(PDF_ORIG), file.exists(PDF_JM), file.exists(SIDECAR),
          file.exists(ANNOT_PY))
dir.create(ANNOT_RES, showWarnings = FALSE, recursive = TRUE)

CLS    <- c("basolateral_out", "apical_out", "mixed")           # canonical class order
ABBR   <- c(basolateral_out = "in", apical_out = "out", mixed = "mix")
abbr   <- function(x) ifelse(is.na(x), "-", ABBR[x])
log_msg <- function(...) message(sprintf("[108] %s", sprintf(...)))

# ---- classify a comment string into an apical class ------------------------
# (legacy 'in'/'apical in'/'ai' input tokens resolve to the basolateral-out class)
classify <- function(s) {
  t <- tolower(trimws(s)); t <- gsub("^[0-9]+[ \t:.)-]*", "", t); t <- trimws(t)
  if (grepl("\\b(apical[ _-]?out|out|ao)\\b", t)) return("apical_out")
  if (grepl("\\b(apical[ _-]?in|in|ai)\\b", t))   return("basolateral_out")
  if (grepl("\\b(mixed|mix)\\b", t))              return("mixed")
  NA_character_
}

# ---- shared render helpers (kept identical to R/09_organoid_refinement/02_island_cleanup_canvas.R
#      / R/10_apical_annotation/01_apical_annotate.R) ---------------------------
instance_outlines <- function(M, B, cx0, cy0) {
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
label_halo <- function(x, y, txt, col, cex = 1.0) {
  off <- 0.6 * (cex / 1.5) * strwidth("0", cex = cex) / 4 + 0.4
  for (dx in c(-1, 1)) for (dy in c(-1, 1))
    text(x + dx * off, y + dy * off, txt, col = "white", cex = cex, font = 2, xpd = NA)
  text(x, y, txt, col = col, cex = cex, font = 2, xpd = NA)
}
place_label <- function(cx, cy, self_pts, other_pts, cw, ch, labW, labH) {
  gap <- 0.7 * labH; rad <- 1.3 * labH
  ang <- seq(0, 2 * pi, length.out = 17)[-17]
  best <- NULL; best_pen <- Inf
  for (i in seq_along(ang)) {
    d    <- c(cos(ang[i]), sin(ang[i]))
    proj <- (self_pts[, 1] - cx) * d[1] + (self_pts[, 2] - cy) * d[2]
    b    <- self_pts[which.max(proj), ]
    half <- abs(d[1]) * labW / 2 + abs(d[2]) * labH / 2
    L    <- c(b[1] + d[1] * (gap + half), b[2] + d[2] * (gap + half))
    pen  <- i * 1e-6 +
      60 * ((L[1] - labW/2 < 2) + (L[1] + labW/2 > cw - 1) +
            (L[2] - labH/2 < 2) + (L[2] + labH/2 > ch - 1))
    if (nrow(other_pts))
      pen <- pen + sum((other_pts[, 1] - L[1])^2 + (other_pts[, 2] - L[2])^2 < rad^2)
    if (pen < best_pen) { best_pen <- pen; best <- list(lx = L[1], ly = L[2], bx = b[1], by = b[2]) }
  }
  best
}

# ============================================================================
# 1. parse each annotator -> data.frame(sid, instance, class)
# ============================================================================
sc <- readRDS(SIDECAR)
PAGE_H_SC <- attr(sc, "page_h_pt"); if (is.null(PAGE_H_SC)) PAGE_H_SC <- PAGE_H_PT
RENDER_SIDS <- attr(sc, "render_sids")
if (is.null(RENDER_SIDS)) RENDER_SIDS <- unique(sc$sid[order(sc$page)])
GROUP_OF <- setNames(sc$group[!duplicated(sc$sid)], sc$sid[!duplicated(sc$sid)])

parse_pdf <- function(pdf, tag) {
  tsv <- file.path(ANNOT_RES, sprintf("consensus_%s_raw.tsv", tag))
  st <- system2(ANNOT_PY_BIN, c(shQuote(ANNOT_PY), shQuote(pdf), shQuote(tsv)),
                stdout = TRUE, stderr = TRUE)
  if (!file.exists(tsv)) stop("annotation extraction failed: ", paste(st, collapse = " | "))
  ann <- read.table(tsv, sep = "\t", header = TRUE, quote = "", comment.char = "",
                    stringsAsFactors = FALSE, fill = TRUE)
  ann <- ann[ann$type %in% ACCEPT_TYPES & nzchar(trimws(ann$contents)) & !is.na(ann$x_left), ]
  ann$apical_class <- vapply(ann$contents, classify, character(1))
  ann <- ann[!is.na(ann$apical_class), ]
  ann$cx <- (ann$x_left + ann$x_right) / 2
  ann$cy <- PAGE_H_SC - (ann$y_bot_pdfcoord + ann$y_top_pdfcoord) / 2

  hits <- list(); n_far <- 0L
  for (i in seq_len(nrow(ann))) {
    sp <- sc[sc$page == ann$page[i], , drop = FALSE]
    if (!nrow(sp)) { n_far <- n_far + 1L; next }
    d <- sqrt((sp$pdf_x - ann$cx[i])^2 + (sp$pdf_y_from_top - ann$cy[i])^2)
    j <- which.min(d)
    if (d[j] > MATCH_RADIUS_PT) { n_far <- n_far + 1L; next }
    hits[[length(hits) + 1L]] <- data.frame(
      sid = sp$sid[j], instance = sp$instance[j],
      apical_class = ann$apical_class[i], stringsAsFactors = FALSE)
  }
  hit <- if (length(hits)) do.call(rbind, hits) else
    data.frame(sid = character(), instance = integer(), apical_class = character())
  log_msg("%s: %d comment(s) classified, matched %d/%d to organoids (%d beyond %dpt)",
          tag, nrow(ann), nrow(hit), nrow(ann), n_far, MATCH_RADIUS_PT)

  # resolve intra-annotator duplicates: same class -> keep; conflicting -> mixed
  if (!nrow(hit)) return(setNames(data.frame(sid=character(), instance=integer(),
                                             class=character()), c("sid","instance",tag)))
  key <- paste(hit$sid, hit$instance)
  agg <- lapply(split(hit, key), function(g) {
    cls <- unique(g$apical_class)
    data.frame(sid = g$sid[1], instance = g$instance[1],
               class = if (length(cls) == 1) cls else "mixed", stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, agg)
  names(out)[names(out) == "class"] <- tag
  out
}

orig <- parse_pdf(PDF_ORIG, "orig")
jm   <- parse_pdf(PDF_JM,   "jm")

# ============================================================================
# 2. consensus table over the FULL organoid roster (centroid sidecar)
# ============================================================================
roster <- unique(sc[, c("sid", "group", "instance")])
con <- merge(roster, orig, by = c("sid", "instance"), all.x = TRUE)
con <- merge(con,   jm,   by = c("sid", "instance"), all.x = TRUE)
names(con)[names(con) == "orig"] <- "class_orig"
names(con)[names(con) == "jm"]   <- "class_jm"

has_o <- !is.na(con$class_orig); has_j <- !is.na(con$class_jm)
con$status <- ifelse(has_o & has_j,
                     ifelse(con$class_orig == con$class_jm, "agree", "disagree"),
                ifelse(has_o & !has_j, "only_orig",
                ifelse(!has_o & has_j, "only_jm", "none")))
con$consensus_class <- ifelse(con$status == "agree", con$class_orig, NA_character_)
con <- con[order(match(con$sid, RENDER_SIDS), con$instance),
           c("sid", "group", "instance", "class_orig", "class_jm",
             "status", "consensus_class")]
write.csv(con, OUT_CSV, row.names = FALSE)
log_msg("wrote %s (%d organoids)", OUT_CSV, nrow(con))

# ---- overall agreement stats + 3x3 confusion matrix ------------------------
both <- con[con$status %in% c("agree", "disagree"), ]
cm   <- table(factor(both$class_orig, levels = CLS), factor(both$class_jm, levels = CLS))
n_both <- nrow(both); n_agree <- sum(con$status == "agree")
po <- if (n_both) sum(diag(cm)) / n_both else NA_real_
pe <- if (n_both) sum(rowSums(cm) * colSums(cm)) / n_both^2 else NA_real_
kappa <- if (n_both) (po - pe) / (1 - pe) else NA_real_
log_msg("both-annotated=%d  agreement=%.1f%%  kappa=%.3f", n_both,
        100 * (n_agree / max(n_both, 1)), kappa)
log_msg("status counts: %s",
        paste(sprintf("%s=%d", names(table(con$status)), table(con$status)), collapse = "  "))

# ============================================================================
# 3. render report
# ============================================================================
STAT_COL <- c(agree = "#00BFC4", disagree = "#FF00FF",
              only_orig = "#E6B800", only_jm = "#E6B800", none = "grey70")
STAT_LWD <- c(agree = 2.0, disagree = 2.4, only_orig = 2.0, only_jm = 2.0, none = 1.0)

dir.create(ANNOT_FIG, showWarnings = FALSE, recursive = TRUE)
pdf(OUT_PDF, width = PAGE_W_IN, height = PAGE_H_IN)

# ---- cover / summary page --------------------------------------------------
disagreements <- con[con$status == "disagree", ]
cover_page <- function() {
  par(mar = c(2, 2, 2, 2)); plot.new()
  cmtxt <- capture.output(print(cm))
  st <- table(factor(con$status, levels = c("agree","disagree","only_orig","only_jm","none")))
  head <- c(
    "Apical-orientation annotation CONSENSUS",
    "  orig = organoid_island_cleanup.pdf      JM = organoid_island_cleanup_INOUT_JM_0623.pdf",
    "",
    "Outline colours:  CYAN = agree   MAGENTA = disagree   YELLOW = one annotator only   grey = neither",
    "Each labelled organoid tag:  O:<orig> J:<JM>   (in / out / mix / - = no call)",
    "",
    sprintf("Organoids total        : %d", nrow(con)),
    sprintf("Both annotated         : %d", n_both),
    sprintf("  agree                : %d (%.1f%%)", st["agree"], 100*st["agree"]/max(n_both,1)),
    sprintf("  disagree             : %d (%.1f%%)", st["disagree"], 100*st["disagree"]/max(n_both,1)),
    sprintf("Only orig / only JM    : %d / %d", st["only_orig"], st["only_jm"]),
    sprintf("Neither                : %d", st["none"]),
    sprintf("Cohen's kappa (3-class): %.3f", kappa),
    "",
    "Confusion matrix  (rows = orig, cols = JM):",
    paste0("  ", cmtxt),
    "",
    sprintf("Disagreements to discuss: %d  (listed on the next page)", nrow(disagreements))
  )
  text(0, 1, paste(head, collapse = "\n"), adj = c(0, 1), cex = 0.9,
       family = "mono", xpd = NA)
}
# dedicated worklist page(s): every disagreement laid out in 3 columns
disagree_page <- function() {
  if (!nrow(disagreements)) return(invisible())
  dz <- sprintf("%-18s #%-3d  O:%-4s J:%-4s",
                sub("AO_", "", disagreements$sid), disagreements$instance,
                abbr(disagreements$class_orig), abbr(disagreements$class_jm))
  NCOL <- 3L; PER_COL <- 40L; PER_PAGE <- NCOL * PER_COL
  npage <- ceiling(length(dz) / PER_PAGE)
  for (pg in seq_len(npage)) {
    par(mar = c(2, 2, 3, 2)); plot.new(); plot.window(c(0, 1), c(0, 1))
    title(sprintf("Disagreements to discuss  (%d total%s)", length(dz),
                  if (npage > 1) sprintf(", page %d/%d", pg, npage) else ""),
          cex.main = 1.2, font.main = 2)
    chunk <- dz[((pg-1)*PER_PAGE + 1):min(pg*PER_PAGE, length(dz))]
    for (cc in seq_len(NCOL)) {
      col_lines <- chunk[((cc-1)*PER_COL + 1):min(cc*PER_COL, length(chunk))]
      col_lines <- col_lines[!is.na(col_lines)]
      if (!length(col_lines)) next
      text((cc - 1) / NCOL + 0.01, 0.97, paste(col_lines, collapse = "\n"),
           adj = c(0, 1), cex = 0.78, family = "mono", xpd = NA)
    }
  }
}
cover_page()
disagree_page()

# ---- one page per section --------------------------------------------------
inst_file <- function(sid) {
  for (suf in c("final", "clean", "split")) {
    f <- file.path(CACHE_DIR, sprintf("instances_%s_%s.rds", suf, sid))
    if (file.exists(f)) return(f)
  }
  file.path(CACHE_DIR, sprintf("instances_%s.rds", sid))
}
# (sid,instance) -> status / tag lookup
con$key <- paste(con$sid, con$instance)
STATUS_OF <- setNames(con$status, con$key)
TAG_OF    <- setNames(sprintf("O:%s J:%s", abbr(con$class_orig), abbr(con$class_jm)), con$key)

for (sid in RENDER_SIDS) {
  xf <- file.path(REG_CACHE, sprintf("nd2final_%s.rds", sid))
  pf <- file.path(CROP_DIR, sprintf("optical_%s.png", sid))
  inf <- inst_file(sid)
  grp <- GROUP_OF[sid]; fc <- "#444444"
  if (!file.exists(xf) || !file.exists(pf) || !file.exists(inf)) {
    par(mar = c(2, 2, 3, 2)); plot.new()
    title(sprintf("%s  (%s)  -  missing inputs", sub("AO_", "", sid), grp), col.main = fc)
    log_msg("%s: SKIP (missing inputs)", sid); next
  }
  Xr <- readRDS(xf); B <- Xr$B_msi_nd2; smn <- Xr$scale_msi_nd2
  cx0 <- Xr$crop[1]; cy0 <- Xr$crop[2]
  om <- png::readPNG(pf); if (length(dim(om)) == 3) om <- om[, , 1]
  cw <- ncol(om); ch <- nrow(om)

  inst <- readRDS(inf); inst <- inst[inst$instance > 0, ]
  W <- max(inst$x); H <- max(inst$y)
  lab <- matrix(0L, W, H); lab[cbind(inst$x, inst$y)] <- as.integer(inst$instance)
  ids <- sort(unique(as.integer(inst$instance)))

  par(mar = c(1, 1, 3, 1))
  plot.new(); plot.window(c(1, cw), c(ch, 1), asp = 1)
  rasterImage(om / max(om), 1, ch, cw, 1, interpolate = TRUE)

  # status counts for this section title
  sst <- STATUS_OF[paste(sid, ids)]
  ag <- sum(sst == "agree", na.rm = TRUE); ds <- sum(sst == "disagree", na.rm = TRUE)
  sg <- sum(sst %in% c("only_orig","only_jm"), na.rm = TRUE)
  title(sprintf("%s   (%s)   -   agree %d  ·  disagree %d  ·  one-only %d",
                sub("AO_", "", sid), grp, ag, ds, sg),
        cex.main = 1.15, col.main = fc, font.main = 2)

  # geometry per organoid
  polys_of <- list(); cen_of <- list(); pts_of <- list()
  for (k in ids) {
    polys <- instance_outlines(lab == k, B, cx0, cy0)
    cen   <- apply_affine(B, matrix(c(mean(inst$x[inst$instance == k]),
                                      mean(inst$y[inst$instance == k])), nrow = 1))
    pts   <- if (length(polys)) do.call(rbind, lapply(polys, function(p) cbind(p$x, p$y)))
             else matrix(c(cen[1] - cx0 + 1, cen[2] - cy0 + 1), 1, 2)
    polys_of[[as.character(k)]] <- polys; cen_of[[as.character(k)]] <- cen; pts_of[[as.character(k)]] <- pts
  }
  # outlines coloured by status
  for (k in ids) {
    sct <- STATUS_OF[paste(sid, k)]; if (is.na(sct)) sct <- "none"
    kc <- STAT_COL[sct]; klw <- STAT_LWD[sct]
    for (poly in polys_of[[as.character(k)]]) polygon(poly$x, poly$y, border = kc, lwd = klw)
  }
  # labels (two-call tag) on every annotated organoid (skip 'none')
  labH <- strheight("0", cex = 1.0)
  for (k in ids) {
    sct <- STATUS_OF[paste(sid, k)]; if (is.na(sct) || sct == "none") next
    ck  <- as.character(k); kc <- STAT_COL[sct]
    txt <- TAG_OF[paste(sid, k)]
    cen <- cen_of[[ck]]; cx <- cen[1] - cx0 + 1; cy <- cen[2] - cy0 + 1
    others <- if (length(ids) > 1) do.call(rbind, pts_of[setdiff(as.character(ids), ck)]) else matrix(numeric(0), 0, 2)
    pl <- place_label(cx, cy, pts_of[[ck]], others, cw, ch, strwidth(txt, cex = 0.85), labH)
    segments(pl$bx, pl$by, pl$lx, pl$ly, col = kc, lwd = 1.0, xpd = NA)
    label_halo(pl$lx, pl$ly, txt, kc, cex = 0.85)
  }
  # scale bar
  umpx <- PX_UM / smn; bp <- SCALE_UM / umpx
  segments(cw - bp - 10, ch - 14, cw - 10, ch - 14, lwd = 3, col = "black")
  text(cw - bp / 2 - 10, ch - 28, sprintf("%d um", SCALE_UM), cex = 0.7)
  log_msg("%s: BF %dx%d, %d organoid(s)  [agree %d disagree %d one-only %d]",
          sid, cw, ch, length(ids), ag, ds, sg)
}
dev.off()
log_msg("DONE -> %s", OUT_PDF)
log_msg("       -> %s", OUT_CSV)
