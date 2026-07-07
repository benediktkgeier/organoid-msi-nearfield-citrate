#!/usr/bin/env Rscript
# 03_build_curated_set.R
# Build the curated "knowns + unknowns" feature MSE for clustering + gradients:
#   knowns   = all reference-list-matched features (metabolite_match_table.csv, matched==TRUE)
#   unknowns = the user's "on tissue" on-tissue selection ions (results/peakme_annotations/)
# Union of features, subset from the 61,124 recall grid (peaks_after_freq.rds, so
# the 23 sub-floor reference knowns are included), all 90,760 pixels, then realized
# in-memory (self-contained). featureData carries known/unknown provenance tags.
# featureData metabolite cols: metabolite_name/adduct/score/hmdb.
#
# Out: cache/peaks_curated.rds, results/curated_feature_table.csv
# Usage: Rscript R/04_ssc_ontissue/03_build_curated_set.R

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_paths.R"))
suppressPackageStartupMessages({ library(Cardinal); library(matter) })

MATCH_CSV <- file.path(RES_DIR, "metabolites", "metabolite_match_table.csv")
ANNOT_DIR <- file.path(RES_DIR, "peakme_annotations")
OUT_RDS   <- file.path(CACHE_DIR, "peaks_curated.rds")
OUT_TBL   <- file.path(RES_DIR, "curated_feature_table.csv")
PPM_MATCH <- 5
log_msg <- function(...) message(sprintf("[21] %s", sprintf(...)))

# ---- 1. reference-list knowns (dedupe isobaric to one feature) -------------
mt <- read.csv(MATCH_CSV, stringsAsFactors = FALSE)
mt <- mt[mt$matched == TRUE, ]
mt <- mt[order(mt$feat_idx, -mt$score), ]
mt_u <- mt[!duplicated(mt$feat_idx), ]            # one row per data feature
sam_feat <- mt_u$feat_idx
log_msg("Reference matched rows: %d -> unique features: %d", nrow(mt), length(sam_feat))

# ---- 2. on-tissue + off-tissue m/z from on-tissue selection annotations -----------------
csvs <- list.files(ANNOT_DIR, pattern = "_annotations\\.csv$", full.names = TRUE)
stopifnot(length(csvs) >= 1)
ann <- do.call(rbind, lapply(csvs, read.csv, stringsAsFactors = FALSE))
lab <- tolower(trimws(ann$label_name))
ot_mz  <- unique(ann$mz_value[lab %in% c("on tissue", "on_tissue")])
off_mz <- unique(ann$mz_value[lab %in% c("off tissue","off_tissue","matrix","noise")])
log_msg("on-tissue selection (%d CSV): on-tissue mz=%d, off-tissue/matrix/noise mz=%d",
        length(csvs), length(ot_mz), length(off_mz))

# ---- 3. recall grid + map on-tissue mz -> feature index --------------------
mse <- readRDS(file.path(CACHE_DIR, "peaks_after_freq.rds"))
mzv <- mz(mse)
nn_feat <- function(q) { j <- which.min(abs(mzv - q)); if (abs(mzv[j]-q)/q*1e6 <= PPM_MATCH) j else NA_integer_ }
ot_feat <- na.omit(vapply(ot_mz, nn_feat, integer(1)))
log_msg("on-tissue mz mapped to grid within %g ppm: %d / %d", PPM_MATCH, length(ot_feat), length(ot_mz))

# ---- 4. union -> subset -----------------------------------------------------
keep <- sort(unique(c(sam_feat, ot_feat)))
log_msg("Curated union: %d features (knowns %d + on-tissue %d - overlap %d)",
        length(keep), length(sam_feat), length(ot_feat),
        length(sam_feat) + length(ot_feat) - length(keep))
cur <- mse[keep, ]
cmz <- mz(cur)

# ---- 5. pixelData sample_id (peaks_after_freq predates annotation) ----------
inv <- load_inventory()
inv$run_name <- sub("\\.imzML$", "", basename(inv$imzml_path), ignore.case = TRUE)
run_to_sid <- setNames(inv$sample_id, inv$run_name)
pixelData(cur)$sample_id <- factor(unname(run_to_sid[as.character(run(cur))]), levels = inv$sample_id)
stopifnot(!any(is.na(pixelData(cur)$sample_id)), nlevels(droplevels(pixelData(cur)$sample_id)) == 20)

# ---- 6. featureData provenance tags ----------------------------------------
mz8772 <- mz(readRDS(file.path(CACHE_DIR, "peaks_combined.rds")))
within5 <- function(q, grid) any(abs(grid - q)/q*1e6 <= PPM_MATCH)
mrow <- mt_u[match(keep, mt_u$feat_idx), ]
fd <- featureData(cur)
fd$mz_curated        <- cmz
fd$is_known          <- keep %in% sam_feat
fd$is_on_tissue      <- keep %in% ot_feat
fd$metabolite_name   <- mrow$name
fd$metabolite_adduct <- mrow$adduct
fd$metabolite_score  <- mrow$score
fd$metabolite_hmdb   <- mrow$hmdb
fd$off_tissue_flagged <- fd$is_known & vapply(cmz, within5, logical(1), grid = off_mz)
fd$in_8772_set       <- vapply(cmz, function(q) any(abs(mz8772 - q) < 1e-6), logical(1))
featureData(cur) <- fd
log_msg("Tags: known=%d, on-tissue=%d, both=%d, off-flagged knowns=%d, in-8772=%d, sub-floor knowns=%d",
        sum(fd$is_known), sum(fd$is_on_tissue), sum(fd$is_known & fd$is_on_tissue),
        sum(fd$off_tissue_flagged), sum(fd$in_8772_set), sum(fd$is_known & !fd$in_8772_set))

# ---- 7. realize self-contained ---------------------------------------------
log_msg("Realizing spectra in-memory (sparse)...")
dense <- as.matrix(spectra(cur))
sp <- matter::sparse_mat(dense, rowMaj = FALSE); rm(dense); gc()
stopifnot(is.null(attr(sp, "path")))
spectra(cur) <- sp
probe <- as.matrix(spectra(cur)[seq_len(min(3, nrow(cur))), round(ncol(cur)*0.6) + 0:4, drop = FALSE])
log_msg("In-memory mid-matrix read OK (%s)", paste(dim(probe), collapse = "x"))

# ---- 8. save ----------------------------------------------------------------
saveRDS(cur, OUT_RDS)
write.csv(as.data.frame(featureData(cur)), OUT_TBL, row.names = FALSE)
log_msg("DONE -> %s (%d feat x %d px) + %s", basename(OUT_RDS), nrow(cur), ncol(cur), basename(OUT_TBL))
