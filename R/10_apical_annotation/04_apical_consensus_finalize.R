#!/usr/bin/env Rscript
# ============================================================================
# 04_apical_consensus_finalize.R - turn the discussed consensus sheet
#   (apical_consensus_report_2annotate.pdf) into analysis-ready apical maps.
# ----------------------------------------------------------------------------
# On the R/10_apical_annotation/03_apical_consensus.R report the two scientists, sitting together:
#   (a) wrote a RESOLVED in/out/mix call on each DISAGREEMENT organoid, and
#   (b) flagged "top pick" on organoids whose basolateral-out/out call is especially
#       clean (gold-standard examples).
#
# The consensus report PDF has the SAME per-section render geometry as the
# island-cleanup canvas (so the R/09_organoid_refinement/02b_island_centroid_sidecar.R sidecar applies), but its pages
# are shifted by +1 (page1 = cover, page2 = disagreement worklist, sections
# start at page 3 vs page 2 on the canvas) -> SECTION_PAGE_OFFSET.
#
# Final consensus class per organoid:
#   - resolved        : a call was written on the consensus sheet            -> use it
#   - agreed          : both annotators already agreed (03_apical_consensus.R status=agree)  -> consensus_class
#   - excluded        : unresolved disagreement / single-annotator-only / neither
#
# In : figures/annotation/apical_consensus_report_2annotate.pdf
#      cache/organoid_island_label_positions_centroid.rds   (R/09_organoid_refinement/02b_island_centroid_sidecar.R)
#      results/annotation/apical_consensus_per_organoid.csv (R/10_apical_annotation/03_apical_consensus.R)
# Out: results/annotation/apical_consensus_final.csv        (master, all organoids)
#      results/annotation/apical_map_consensus.csv          (analysis map: consensus set - the headline input)
#      results/annotation/apical_map_toppick.csv            (analysis map: top picks only)
# Usage: Rscript R/10_apical_annotation/04_apical_consensus_finalize.R
#   then: Rscript R/11_per_organoid_final/03_apical_report.R results/annotation/apical_map_consensus.csv consensus
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))

ANNOT_FIG <- file.path(FIG_DIR, "annotation")
ANNOT_RES <- file.path(RES_DIR, "annotation")
PDF      <- file.path(ANNOT_FIG, "apical_consensus_report_2annotate.pdf")
SIDECAR  <- file.path(CACHE_DIR, "organoid_island_label_positions_centroid.rds")
CONSENSUS_CSV <- file.path(ANNOT_RES, "apical_consensus_per_organoid.csv")   # from 03_apical_consensus.R
ANNOT_PY     <- file.path(PROJECT_ROOT, "py", "extract_pdf_annots.py")
ANNOT_PY_BIN <- Sys.getenv("MSI_PYTHON", "python")
ACCEPT_TYPES <- c("/FreeText", "/Text", "/Square", "/Circle", "/Highlight",
                  "/StrikeOut", "/Underline", "/Squiggly")
MATCH_RADIUS_PT    <- 60
SECTION_PAGE_OFFSET <- 1L   # consensus report has 1 extra page (worklist) before sections
stopifnot(file.exists(PDF), file.exists(SIDECAR), file.exists(CONSENSUS_CSV), file.exists(ANNOT_PY))

OUT_MASTER   <- file.path(ANNOT_RES, "apical_consensus_final.csv")
OUT_CONSENS  <- file.path(ANNOT_RES, "apical_map_consensus.csv")
OUT_TOPPICK  <- file.path(ANNOT_RES, "apical_map_toppick.csv")
OUT_TSV      <- file.path(ANNOT_RES, "consensus_2annotate_raw.tsv")
log_msg <- function(...) message(sprintf("[109] %s", sprintf(...)))

classify <- function(s) {
  t <- tolower(trimws(s)); t <- gsub("^[0-9]+[ \t:.)-]*", "", t); t <- trimws(t)
  if (grepl("\\b(apical[ _-]?out|out|ao)\\b", t)) return("apical_out")
  if (grepl("\\b(apical[ _-]?in|in|ai)\\b", t))   return("basolateral_out")
  if (grepl("\\b(mixed|mix)\\b", t))              return("mixed")
  NA_character_
}

# ---- 1. extract + classify comments on the consensus sheet -----------------
st <- system2(ANNOT_PY_BIN, c(shQuote(ANNOT_PY), shQuote(PDF), shQuote(OUT_TSV)),
              stdout = TRUE, stderr = TRUE)
if (!file.exists(OUT_TSV)) stop("annotation extraction failed: ", paste(st, collapse = " | "))
ann <- read.table(OUT_TSV, sep = "\t", header = TRUE, quote = "", comment.char = "",
                  stringsAsFactors = FALSE, fill = TRUE)
ann <- ann[ann$type %in% ACCEPT_TYPES & nzchar(trimws(ann$contents)) & !is.na(ann$x_left), ]
ann$is_top <- grepl("top[ _-]?pick", tolower(ann$contents))
ann$cls    <- ifelse(ann$is_top, NA_character_, vapply(ann$contents, classify, character(1)))
ann <- ann[ann$is_top | !is.na(ann$cls), ]                 # keep top-picks + class calls

# ---- 2. map each comment -> nearest organoid centroid (page offset) --------
sc <- readRDS(SIDECAR)
PAGE_H <- attr(sc, "page_h_pt"); if (is.null(PAGE_H)) PAGE_H <- 720
ann$cx <- (ann$x_left + ann$x_right) / 2
ann$cy <- PAGE_H - (ann$y_bot_pdfcoord + ann$y_top_pdfcoord) / 2
ann$sc_page <- ann$page - SECTION_PAGE_OFFSET
ann$sid <- NA_character_; ann$instance <- NA_integer_; ann$dist <- NA_real_
for (i in seq_len(nrow(ann))) {
  sp <- sc[sc$page == ann$sc_page[i], , drop = FALSE]; if (!nrow(sp)) next
  d <- sqrt((sp$pdf_x - ann$cx[i])^2 + (sp$pdf_y_from_top - ann$cy[i])^2)
  j <- which.min(d)
  if (d[j] <= MATCH_RADIUS_PT) { ann$sid[i] <- sp$sid[j]; ann$instance[i] <- sp$instance[j]; ann$dist[i] <- d[j] }
}
# self-validation: the +1 offset must give tight matches
mok <- !is.na(ann$dist)
log_msg("page offset %d -> matched %d/%d comments (median dist %.1f pt, max %.1f)",
        SECTION_PAGE_OFFSET, sum(mok), nrow(ann), median(ann$dist[mok]), max(ann$dist[mok]))
if (median(ann$dist[mok]) > 30) stop("matches are loose - SECTION_PAGE_OFFSET likely wrong")
ann <- ann[mok, ]
ann$key <- paste(ann$sid, ann$instance)

# resolved class per organoid (conflicting calls on one organoid -> mixed)
cc <- ann[!ann$is_top, ]
new_cls <- tapply(cc$cls, cc$key, function(v) { u <- unique(v); if (length(u) == 1) u else "mixed" })
n_conf  <- sum(tapply(cc$cls, cc$key, function(v) length(unique(v)) > 1))
top_keys <- unique(ann$key[ann$is_top])
log_msg("resolved class on %d organoid(s) (%d with conflicting calls -> mixed); %d top-pick organoid(s)",
        length(new_cls), n_conf, length(top_keys))

# ---- 3. merge with 03_apical_consensus.R output -> final class -------------
con <- read.csv(CONSENSUS_CSV, stringsAsFactors = FALSE)
con$key <- paste(con$sid, con$instance)
con$new_class <- unname(new_cls[con$key])
con$top_pick  <- con$key %in% top_keys

reviewed <- !is.na(con$new_class)
con$final_class <- ifelse(reviewed, con$new_class,
                   ifelse(con$status == "agree", con$consensus_class, NA_character_))
con$class_source <- ifelse(reviewed & con$status == "disagree", "resolved",
                    ifelse(reviewed,                          "resolved_single",
                    ifelse(con$status == "agree",             "agreed", "excluded")))

master <- con[, c("sid", "group", "instance", "class_orig", "class_jm", "status",
                  "new_class", "final_class", "class_source", "top_pick")]
master <- master[order(master$sid, master$instance), ]
write.csv(master, OUT_MASTER, row.names = FALSE)
log_msg("wrote %s (%d organoids)", OUT_MASTER, nrow(master))

# ---- 4. analysis maps (sid, group, instance, apical_class) -----------------
cons <- master[!is.na(master$final_class), c("sid", "group", "instance", "final_class")]
names(cons)[4] <- "apical_class"
write.csv(cons, OUT_CONSENS, row.names = FALSE)

tp <- master[master$top_pick & !is.na(master$final_class), c("sid", "group", "instance", "final_class")]
names(tp)[4] <- "apical_class"
write.csv(tp, OUT_TOPPICK, row.names = FALSE)

log_msg("consensus map -> %s (%d organoids)", OUT_CONSENS, nrow(cons))
log_msg("top-pick  map -> %s (%d organoids)", OUT_TOPPICK, nrow(tp))
cat("\n[109] class_source breakdown (all organoids):\n"); print(table(master$class_source))
cat("\n[109] consensus set apical_class:\n"); print(table(cons$apical_class))
cat("\n[109] top-pick set apical_class:\n");   print(table(tp$apical_class))
cat("\n[109] NEXT: the consensus map is the DEFAULT for the report; for the top-pick subset:\n")
cat("  Rscript R/11_per_organoid_final/03_apical_report.R                                              # consensus (default)\n")
cat("  Rscript R/11_per_organoid_final/03_apical_report.R results/annotation/apical_map_toppick.csv toppick\n")
