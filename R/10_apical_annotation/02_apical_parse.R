#!/usr/bin/env Rscript
# ============================================================================
# 02_apical_parse.R - parse the EDITED annotation PDF
#   (R/10_apical_annotation/01_apical_annotate.R output)
#   back into a machine-readable per-organoid apical-orientation table.
#
#   Mirrors D:/R/PeakMe/PeakMe_GCPL/phaseZC_pdf_to_pairs.R: calls the python
#   pypdf extractor, then maps each Adobe comment to a (sid, instance) using the
#   sidecar of label positions, and parses the class keyword from the comment.
#
#   Mapping per annotation:
#     (a) if the comment text starts with an integer matching an organoid id on
#         that page -> map by ID (unambiguous);
#     (b) else -> nearest sidecar label on the page within MATCH_RADIUS_PT.
#   Class keywords (lowercased): apical_out <- out|apical out|ao ;
#     basolateral_out <- in|apical in|ai (legacy annotator synonyms that resolve
#     to the basolateral-out class) ; mixed <- mix|mixed.
#
# Usage:
#   Rscript R/10_apical_annotation/02_apical_parse.R <edited.pdf> [sidecar.rds]
#   (sidecar defaults to cache/organoid_apical_label_positions.rds)
#
# Out: results/annotation/organoid_apical_annotations.csv  (sid, instance, group,
#        apical_class, note, page, match_mode, dist_pt)
#      results/annotation/organoid_apical_annotations.rds
#      results/annotation/organoid_apical_parse_log.csv     (audit: every annot)
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))

MATCH_RADIUS_PT  <- 25
ANNOT_PY         <- file.path(PROJECT_ROOT, "py", "extract_pdf_annots.py")
ANNOT_PY_BIN     <- "C:/Users/bened/.virtualenvs/r-reticulate/Scripts/python.exe"
ACCEPT_TYPES     <- c("/FreeText", "/Text", "/Highlight", "/StrikeOut",
                      "/Underline", "/Squiggly", "/Square", "/Circle", "/Popup")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript R/103_organoid_apical_parse.R <edited.pdf> [sidecar.rds]")
EDITED_PDF  <- args[1]
SIDECAR_RDS <- if (length(args) >= 2) args[2] else file.path(CACHE_DIR, "organoid_apical_label_positions.rds")
stopifnot(file.exists(EDITED_PDF), file.exists(SIDECAR_RDS), file.exists(ANNOT_PY))

# dedicated annotation pre-step folder (parallel to figures/annotation)
ANNOT_RES <- file.path(RES_DIR, "annotation")
dir.create(ANNOT_RES, showWarnings = FALSE, recursive = TRUE)
OUT_CSV <- file.path(ANNOT_RES, "organoid_apical_annotations.csv")
OUT_RDS <- file.path(ANNOT_RES, "organoid_apical_annotations.rds")
OUT_LOG <- file.path(ANNOT_RES, "organoid_apical_parse_log.csv")
OUT_TSV <- file.path(ANNOT_RES, "organoid_apical_annotations_raw.tsv")

# ---- classify a comment string into an apical class ------------------------
classify <- function(s) {
  t <- tolower(trimws(s))
  t <- gsub("^[0-9]+[ \t:.)-]*", "", t)          # drop a leading "<id>" prefix
  t <- trimws(t)
  if (grepl("\\b(apical[ _-]?out|out|ao)\\b", t))  return("apical_out")
  if (grepl("\\b(apical[ _-]?in|in|ai)\\b", t))    return("basolateral_out")
  if (grepl("\\b(mixed|mix)\\b", t))               return("mixed")
  NA_character_
}
# leading explicit organoid id, if any
lead_id <- function(s) {
  m <- regmatches(s, regexpr("^\\s*[0-9]+", s))
  if (length(m)) as.integer(trimws(m)) else NA_integer_
}

# ---- extract annotations via python helper ---------------------------------
status <- system2(ANNOT_PY_BIN, c(shQuote(ANNOT_PY), shQuote(EDITED_PDF), shQuote(OUT_TSV)),
                  stdout = TRUE, stderr = TRUE)
if (!file.exists(OUT_TSV)) stop("annotation extraction failed: ", paste(status, collapse = " | "))
ann <- read.table(OUT_TSV, sep = "\t", header = TRUE, quote = "", comment.char = "",
                  stringsAsFactors = FALSE, fill = TRUE)
ann <- ann[ann$type %in% ACCEPT_TYPES & nzchar(trimws(ann$contents)), , drop = FALSE]
ann <- ann[!is.na(ann$x_left), , drop = FALSE]
cat(sprintf("[103] %d candidate annotation(s) with text\n", nrow(ann)))

sc <- readRDS(SIDECAR_RDS)
PAGE_H_PT <- attr(sc, "page_h_pt"); if (is.null(PAGE_H_PT)) PAGE_H_PT <- 10 * 72
GROUP_OF  <- setNames(sc$group[!duplicated(sc$sid)], sc$sid[!duplicated(sc$sid)])

if (nrow(ann)) {
  ann$cx_pt          <- (ann$x_left + ann$x_right) / 2
  ann$cy_from_top_pt <- PAGE_H_PT - (ann$y_bot_pdfcoord + ann$y_top_pdfcoord) / 2
}

# ---- map each annotation -> (sid, instance) + class ------------------------
log_rows <- list(); hits <- list()
for (i in seq_len(nrow(ann))) {
  pg <- ann$page[i]; contents <- ann$contents[i]
  spg <- sc[sc$page == pg, , drop = FALSE]
  rec <- function(action, sid = NA, inst = NA, cls = NA, mode = NA, dist = NA)
    data.frame(page = pg, sid = sid, instance = inst, contents = contents,
               apical_class = cls, match_mode = mode, dist_pt = dist,
               action = action, stringsAsFactors = FALSE)
  if (!nrow(spg)) { log_rows[[length(log_rows)+1L]] <- rec("no-organoids-on-page"); next }
  sid <- spg$sid[1]
  cls <- classify(contents)
  lid <- lead_id(contents)
  mode <- NA_character_; dist <- NA_real_; inst <- NA_integer_
  if (!is.na(lid) && lid %in% spg$instance) {
    inst <- lid; mode <- "by-id"
    j <- which(spg$instance == lid)
    dist <- sqrt((spg$pdf_x[j] - ann$cx_pt[i])^2 + (spg$pdf_y_from_top[j] - ann$cy_from_top_pt[i])^2)
  } else {
    d <- sqrt((spg$pdf_x - ann$cx_pt[i])^2 + (spg$pdf_y_from_top - ann$cy_from_top_pt[i])^2)
    j <- which.min(d)
    if (length(j) && is.finite(d[j]) && d[j] <= MATCH_RADIUS_PT) {
      inst <- spg$instance[j]; mode <- "by-nearest"; dist <- d[j]
    }
  }
  if (is.na(inst)) { log_rows[[length(log_rows)+1L]] <- rec("no-match-within-radius", sid = sid); next }
  if (is.na(cls))  { log_rows[[length(log_rows)+1L]] <- rec("unparseable-class", sid = sid, inst = inst, mode = mode, dist = dist); next }
  hits[[length(hits)+1L]] <- data.frame(
    page = pg, sid = sid, group = unname(GROUP_OF[sid]), instance = inst,
    apical_class = cls, note = contents, match_mode = mode, dist_pt = round(dist, 1),
    stringsAsFactors = FALSE)
  log_rows[[length(log_rows)+1L]] <- rec("accept", sid = sid, inst = inst, cls = cls, mode = mode, dist = round(dist, 1))
}

hit_df <- if (length(hits)) do.call(rbind, hits) else
  data.frame(page = integer(), sid = character(), group = character(), instance = integer(),
             apical_class = character(), note = character(), match_mode = character(),
             dist_pt = numeric(), stringsAsFactors = FALSE)

# ---- conflicts: same (sid, instance) annotated more than once --------------
if (nrow(hit_df)) {
  key <- paste(hit_df$sid, hit_df$instance)
  dup <- key %in% key[duplicated(key)]
  if (any(dup)) {
    confl <- unique(key[dup]);
    cat(sprintf("[103] WARNING: %d organoid(s) have conflicting/duplicate annotations:\n", length(confl)))
    for (k in confl) {
      sub <- hit_df[key == k, ]
      cat(sprintf("   %s : %s\n", k, paste(sprintf("%s('%s')", sub$apical_class, sub$note), collapse = " ; ")))
    }
    cat("   -> kept all rows; resolve in PDF and re-run, or dedup downstream.\n")
  }
}

# ---- organoids with NO annotation (so nothing is silently missed) ----------
if (nrow(sc)) {
  keyed <- if (nrow(hit_df)) paste(hit_df$sid, hit_df$instance) else character(0)
  miss  <- sc[!(paste(sc$sid, sc$instance) %in% keyed), c("sid", "group", "instance")]
  if (nrow(miss)) cat(sprintf("[103] NOTE: %d/%d organoid(s) not yet annotated.\n", nrow(miss), nrow(sc)))
}

hit_df <- hit_df[order(hit_df$sid, hit_df$instance), ]
write.csv(hit_df[, c("sid", "instance", "group", "apical_class", "note", "page", "match_mode", "dist_pt")],
          OUT_CSV, row.names = FALSE)
saveRDS(hit_df, OUT_RDS)
if (length(log_rows)) write.csv(do.call(rbind, log_rows), OUT_LOG, row.names = FALSE)

cat(sprintf("[103] %d organoid annotation(s) accepted (by-id: %d, by-nearest: %d)\n",
            nrow(hit_df), sum(hit_df$match_mode == "by-id"), sum(hit_df$match_mode == "by-nearest")))
if (nrow(hit_df)) print(table(hit_df$group, hit_df$apical_class))
cat(sprintf("[103] DONE -> %s\n", OUT_CSV))
