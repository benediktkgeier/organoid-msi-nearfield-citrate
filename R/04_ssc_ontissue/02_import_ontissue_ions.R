#!/usr/bin/env Rscript
# ============================================================================
# 02_import_ontissue_ions.R - Phase 04 step 2: attach on-tissue selection annotations to the
#                           20-section SCiLS combined MSE; build clean MSE.
# ============================================================================
# Single-pipeline (the DEV/BACKGROUND split is retired, 2026-06-16). Operates on
# the self-contained 20-section MSE cache/peaks_combined.rds. Labels come from
# the 3 annotated sections; the m/z grid is shared across all 20 sections, so a
# feature's label applies wherever that m/z appears (off_tissue/matrix are
# feature-level, treated as section-independent; unannotated -> keep).
#
# Reads the per-sample on-tissue selection CSVs from results/peakme_annotations/
#   (columns: mz_value, label_name[, starred, confidence, annotator, annotated_at]).
# Collapses to ONE label per m/z by KEEP-OVERRIDE precedence (NOT majority vote):
#   an ion is KEPT (on_tissue) if unannotated, or any sample said on_tissue/unclear;
#   a keep signal OVERRIDES off_tissue/matrix/noise on other samples. An ion is
#   DROPPED only when EVERY annotating sample gave a drop label and none kept it.
#
# Outputs:
#   cache/peaks_combined_annot.rds   - full 6845-feature MSE + annotation fData
#                                      (kept intact: tissueness denominator in 04_ssc_tissue_mask.R)
#   cache/peaks_combined_clean.rds   - off_tissue/matrix/noise features dropped
#                                      (SSC input); unannotated features KEPT
#   results/annotation_import_summary.csv
#
# Adapts match/coverage logic of
#   D:/R/Projects/BIG_MSI/HP_Comparative_2026/peakme/peakme_export.R (lines 200-334).
#
# Usage: Rscript R/04_ssc_ontissue/02_import_ontissue_ions.R
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_paths.R"))
suppressPackageStartupMessages({
  library(Cardinal)
})

MSE_IN     <- file.path(CACHE_DIR, "peaks_combined.rds")
ANNOT_DIR  <- file.path(RES_DIR, "peakme_annotations")
ANNOT_GLOB <- "peakme_*_annotations.csv"
# Same-grid matching: the annotations are exported FROM this exact MSE, so each
# annotation m/z equals an MSE feature m/z. We still match to the nearest feature
# within half the 25 ppm bin (12.5 ppm) for robustness; offsets are ~0 ppm.
NN_PPM_TOL <- 12.5

# ---- Cross-sample label resolution (KEEP-OVERRIDE precedence) --------------
# An ion is KEPT (resolved to "on_tissue") if ANY of these hold:
#   - it is unannotated in all CSVs, OR
#   - any sample labeled it KEEP_SIGNALS ("on_tissue" or "unclear").
# A KEEP signal OVERRIDES any drop label (incl. off_tissue/matrix/noise) on
# another sample. An ion is DROPPED only when every annotating sample gave a
# drop label AND none said on_tissue/unclear.
KEEP_SIGNALS <- c("on_tissue", "unclear")
DROP_LABELS  <- c("off_tissue", "matrix", "noise")  # affirmative drop categories
OFF_LABELS   <- c("off_tissue", "matrix")           # tissueness denominator (step 20)

OUT_ANNOT  <- file.path(CACHE_DIR, "peaks_combined_annot.rds")
OUT_CLEAN  <- file.path(CACHE_DIR, "peaks_combined_clean.rds")
OUT_SUMM   <- file.path(RES_DIR, "annotation_import_summary.csv")

# ---- Label normalization: lowercase, spaces/hyphens -> underscore ----------
norm_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[[:space:]\\-]+", "_", x)   # "off tissue"/"off-tissue" -> "off_tissue"
  x
}

# ===========================================================================
# 1. Load the 20-section MSE
# ===========================================================================
if (!file.exists(MSE_IN)) {
  stop("MSE not found: ", MSE_IN)
}
cat(sprintf("[04] Loading MSE: %s\n", basename(MSE_IN)))
mse <- readRDS(MSE_IN)
n_feat <- nrow(mse)
mz_vec <- mz(mse)
cat(sprintf("[04] MSE: %d features x %d pixels\n", n_feat, ncol(mse)))

# Guard: per-sample run factor must be intact (20 sections) for downstream SSC.
n_runs <- nlevels(run(mse))
cat(sprintf("[04] run() levels (samples): %d\n", n_runs))
stopifnot(n_runs == 20)

# ===========================================================================
# 2. Read all per-sample annotation CSVs
# ===========================================================================
csv_files <- list.files(ANNOT_DIR, pattern = utils::glob2rx(ANNOT_GLOB),
                         full.names = TRUE)
if (length(csv_files) == 0) {
  stop("No annotation CSVs found in ", ANNOT_DIR,
       "\n  Expected files matching '", ANNOT_GLOB, "'.",
       "\n  Drop your on-tissue selection exports there and re-run.")
}
cat(sprintf("[04] Found %d annotation CSV(s):\n", length(csv_files)))
for (f in csv_files) cat("       - ", basename(f), "\n", sep = "")

ann_list <- list()
for (f in csv_files) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  miss <- setdiff(c("mz_value", "label_name"), colnames(df))
  if (length(miss) > 0) {
    stop("CSV '", basename(f), "' is missing required column(s): ",
         paste(miss, collapse = ", "))
  }
  df$source_file <- basename(f)
  if (!"starred" %in% colnames(df))    df$starred    <- NA
  if (!"confidence" %in% colnames(df)) df$confidence <- NA
  if (!"annotator" %in% colnames(df))  df$annotator  <- NA_character_
  ann_list[[f]] <- df[, c("mz_value", "label_name", "starred",
                          "confidence", "annotator", "source_file")]
}
ann <- do.call(rbind, ann_list)
ann$label_norm <- norm_label(ann$label_name)
cat(sprintf("[04] Total annotation rows across CSVs: %d\n", nrow(ann)))

# ===========================================================================
# 3. Match each annotation row's mz_value to the NEAREST MSE feature within
#    NN_PPM_TOL (ppm). Summary-counted, not per-row warnings.
# ===========================================================================
feat_idx <- integer(nrow(ann))
off_ppm  <- numeric(nrow(ann))
for (i in seq_len(nrow(ann))) {
  mzi   <- ann$mz_value[i]
  d_ppm <- abs(mz_vec - mzi) / mzi * 1e6
  nn    <- which.min(d_ppm)
  off_ppm[i] <- d_ppm[nn]
  feat_idx[i] <- if (d_ppm[nn] <= NN_PPM_TOL) nn else NA_integer_
}
n_unmatched <- sum(is.na(feat_idx))
cat(sprintf("[04] Matching to rebuilt grid within %.1f ppm: matched %d / %d rows (%.1f%%); %d unmatched.\n",
            NN_PPM_TOL, sum(!is.na(feat_idx)), nrow(ann),
            100 * mean(!is.na(feat_idx)), n_unmatched))
cat(sprintf("[04] Matched-row offset ppm: median %.2f, 95th %.2f, max %.2f\n",
            median(off_ppm[!is.na(feat_idx)]),
            quantile(off_ppm[!is.na(feat_idx)], 0.95),
            max(off_ppm[!is.na(feat_idx)])))
if (n_unmatched > 0) {
  ann <- ann[!is.na(feat_idx), ]
  feat_idx <- feat_idx[!is.na(feat_idx)]
}
ann$feat_idx <- feat_idx
cat(sprintf("[04] Matched %d annotation rows to MSE features.\n", nrow(ann)))

# ===========================================================================
# 4. Resolve one label per feature via KEEP-OVERRIDE precedence
#    (see KEEP_SIGNALS / DROP_LABELS comment above)
# ===========================================================================
# Default EVERY feature to on_tissue -> covers the "unannotated -> keep" rule.
fd <- as.data.frame(featureData(mse))
fd$peakme_label          <- "on_tissue"   # default keep; overridden below if all-drop
fd$peakme_votes          <- NA_character_  # raw cross-sample tally (diagnostic)
fd$peakme_n_annotated    <- 0L
fd$peakme_label_disagree <- FALSE          # >1 distinct raw label across samples
fd$peakme_keep_override  <- FALSE          # TRUE when a KEEP signal beat a drop label
fd$peakme_starred        <- NA
fd$peakme_confidence     <- NA_real_
fd$peakme_annotator      <- NA_character_

# Drop-label representative: prefer background (off_tissue/matrix) over noise so
# the ion still feeds the tissueness denominator in step 20; most-frequent wins.
pick_drop <- function(tab) {
  bg <- tab[names(tab) %in% OFF_LABELS]
  if (length(bg) > 0) return(names(bg)[which.max(bg)])  # off_tissue/matrix
  names(tab)[which.max(tab)]                            # else noise (or other)
}

for (fi in sort(unique(ann$feat_idx))) {
  rows <- ann[ann$feat_idx == fi, , drop = FALSE]
  labs <- rows$label_norm
  tab  <- sort(table(labs), decreasing = TRUE)

  has_keep  <- any(labs %in% KEEP_SIGNALS)
  drop_labs <- labs[labs %in% DROP_LABELS]
  all_drop  <- length(drop_labs) == length(labs) && length(labs) > 0

  if (has_keep || !all_drop) {
    # any on_tissue/unclear, or any non-drop/unknown label -> KEEP
    resolved <- "on_tissue"
    fd$peakme_keep_override[fi] <- has_keep && length(drop_labs) > 0
  } else {
    resolved <- pick_drop(table(drop_labs))   # every annotating sample said drop
  }

  fd$peakme_label[fi]          <- resolved
  fd$peakme_votes[fi]          <- paste(sprintf("%s:%d", names(tab), as.integer(tab)),
                                        collapse = ";")
  fd$peakme_n_annotated[fi]    <- nrow(rows)
  fd$peakme_label_disagree[fi] <- length(unique(labs)) > 1
  st <- suppressWarnings(as.logical(rows$starred))
  fd$peakme_starred[fi]        <- if (any(st, na.rm = TRUE)) TRUE else FALSE
  cf <- suppressWarnings(as.numeric(rows$confidence))
  fd$peakme_confidence[fi]     <- if (all(is.na(cf))) NA_real_ else median(cf, na.rm = TRUE)
  an <- unique(rows$annotator[!is.na(rows$annotator) & nzchar(rows$annotator)])
  fd$peakme_annotator[fi]      <- if (length(an)) paste(an, collapse = ",") else NA_character_
}

featureData(mse)$peakme_label          <- fd$peakme_label
featureData(mse)$peakme_votes          <- fd$peakme_votes
featureData(mse)$peakme_n_annotated    <- fd$peakme_n_annotated
featureData(mse)$peakme_label_disagree <- fd$peakme_label_disagree
featureData(mse)$peakme_keep_override  <- fd$peakme_keep_override
featureData(mse)$peakme_starred        <- fd$peakme_starred
featureData(mse)$peakme_confidence     <- fd$peakme_confidence
featureData(mse)$peakme_annotator      <- fd$peakme_annotator

# ===========================================================================
# 5. Coverage + label breakdown
# ===========================================================================
n_annotated <- sum(fd$peakme_n_annotated > 0)
n_disagree  <- sum(fd$peakme_label_disagree, na.rm = TRUE)
n_override  <- sum(fd$peakme_keep_override, na.rm = TRUE)
cat("\n[04] ---- Coverage --------------------------------------------------\n")
cat(sprintf("       Total features         : %d\n", n_feat))
cat(sprintf("       Annotated (>=1 sample) : %d (%.1f%%)\n",
            n_annotated, 100 * n_annotated / n_feat))
cat(sprintf("       Unannotated -> on_tissue: %d\n", n_feat - n_annotated))
cat(sprintf("       Mixed raw labels        : %d\n", n_disagree))
cat(sprintf("       keep-override fired     : %d (on_tissue/unclear beat a drop label)\n",
            n_override))
cat("[04] ---- Resolved-label breakdown ---------------------------------\n")
lab_tab <- sort(table(fd$peakme_label), decreasing = TRUE)
for (nm in names(lab_tab)) {
  cat(sprintf("       %-20s %d\n", nm, lab_tab[[nm]]))
}

# ===========================================================================
# 6. Save annotated MSE + build clean MSE
# ===========================================================================
saveRDS(mse, OUT_ANNOT)
cat(sprintf("\n[04] Saved annotated MSE: %s\n", basename(OUT_ANNOT)))

remove_mask <- fd$peakme_label %in% DROP_LABELS   # only affirmative all-drop ions
n_remove    <- sum(remove_mask)
n_keep      <- n_feat - n_remove
if (n_keep == 0) stop("All features would be dropped - check labels_to_remove.")

mse_clean <- mse[!remove_mask, ]
saveRDS(mse_clean, OUT_CLEAN)
cat(sprintf("[04] Clean MSE: dropped %d (%s), kept %d features -> %s\n",
            n_remove, paste(DROP_LABELS, collapse = "/"), n_keep,
            basename(OUT_CLEAN)))

# ===========================================================================
# 7. Summary CSV
# ===========================================================================
summ <- data.frame(
  mz                    = mz_vec,
  peakme_label          = fd$peakme_label,
  peakme_votes          = fd$peakme_votes,
  peakme_n_annotated    = fd$peakme_n_annotated,
  peakme_label_disagree = fd$peakme_label_disagree,
  peakme_keep_override  = fd$peakme_keep_override,
  peakme_starred        = fd$peakme_starred,
  peakme_confidence     = fd$peakme_confidence,
  kept_in_clean         = !remove_mask,
  stringsAsFactors      = FALSE
)
write.csv(summ, OUT_SUMM, row.names = FALSE)
cat(sprintf("[04] Wrote summary: %s\n", basename(OUT_SUMM)))
cat("[04] DONE.\n")
