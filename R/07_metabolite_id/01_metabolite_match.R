#!/usr/bin/env Rscript
# 01_metabolite_match.R
# Match a published reference metabolite m/z list (negative-mode spatial
# metabolomics; supplement sheets ST1 liver + ST3 small intestine) to our SCiLS
# feature grid, using ACCURATE theoretical masses computed from HMDB molecular
# formulas (NOT the supplement's 2-decimal m/z, which is too coarse for 20 ppm).
# Adds extra negative-adduct forms of citrate + lactate (user request).
# (Reference source: keep the proper DOI/PMID citation in the manuscript.)
#
# Theoretical m/z = nM * EXACT_MASS(HMDB) + adduct_offset (deterministic).
# Match each compound to the NEAREST data feature within 20 ppm; report
# delta m/z (mDa + ppm) and a score = 100*(1-|ppm|/20).
#
# Grid: peaks_after_freq.rds (61,124 feat, richest, best recall); flag which
# matched features survive into the 8,772 analysis set (peaks_combined.rds).
#
# Out: results/metabolites/metabolite_match_table.csv
# Usage: Rscript R/07_metabolite_id/01_metabolite_match.R

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_paths.R"))
suppressPackageStartupMessages({ library(Cardinal); library(readxl) })

XLSX <- file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Supplement_metabolites.xlsx")
HMDB <- "D:/R/PeakMe/PeakMe_GCPL/phaseG_db_cache/HMDB.tsv"
OUT_DIR <- file.path(RES_DIR, "metabolites"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
PPM_TOL <- 20
MZ_LO <- 100; MZ_HI <- 900
log_msg <- function(...) message(sprintf("[80] %s", sprintf(...)))

# ---- adduct offsets (monoisotopic, negative mode) --------------------------
mH <- 1.007825031898; me <- 0.000548579909; mO <- 15.994914619
H2O <- 2*mH + mO; Na <- 22.989769282; K <- 38.963706487; Cl <- 34.968852682
ADD <- list(
  "[M-H]-"      = list(nM=1, off = -mH + me),
  "[M-H2O-H]-"  = list(nM=1, off = -H2O - mH + me),
  "[M+Cl]-"     = list(nM=1, off = +Cl + me),
  "[M+Na-2H]-"  = list(nM=1, off = +Na - 2*mH + me),
  "[M+K-2H]-"   = list(nM=1, off = +K  - 2*mH + me),
  "[2M-H]-"     = list(nM=2, off = -mH + me))
norm_add <- function(a) gsub("\\s", "", a)

# ---- HMDB: accession -> EXACT_MASS, FORMULA --------------------------------
log_msg("Loading HMDB.tsv ...")
hd <- read.delim(HMDB, stringsAsFactors = FALSE, quote = "")
em <- setNames(as.numeric(hd$EXACT_MASS), hd$HMDB_ID)
fm <- setNames(hd$FORMULA, hd$HMDB_ID)
nm <- setNames(hd$GENERIC_NAME, hd$HMDB_ID)
log_msg("HMDB rows: %d", nrow(hd))

# ---- reference list ST1 + ST3 -----------------------------------------------
read_st <- function(sheet, src) {
  d <- as.data.frame(read_excel(XLSX, sheet = sheet))
  data.frame(source = src, name = d$Name, adduct = norm_add(d$Ion),
             published_mz = suppressWarnings(as.numeric(d$`m/z`)),
             hmdb = trimws(d$`HMDB Accession*`), stringsAsFactors = FALSE)
}
lst <- rbind(read_st("ST1","ST1_liver"), read_st("ST3","ST3_intestine"))
lst <- lst[!is.na(lst$hmdb) & nzchar(lst$hmdb), ]

# NOTE (2026-06-22): the previously-injected extra citrate/lactate adduct forms
# ([M+Na-2H]-/[M+K-2H]-/[M+Cl]-/[2M-H]-) were REMOVED. The "citrate adduct-switching
# / total-citrate" idea was an artifact (unprovable here; adduct ion images are noise).
# Citrate is matched as [M-H]- only via the reference supplement sheets.

# Normalize HMDB accessions: the supplement mixes old 5-digit IDs (HMDB02985)
# and IDs with internal spaces (HMDB 0011494) with the modern 7-digit form.
# Strip spaces, extract HMDB+digits, zero-pad to 7 -> recovers entries that an
# exact-string lookup misses. (METPA*/non-HMDB IDs return NA and stay unmatched.)
norm_hmdb <- function(x) {
  x <- gsub("[[:space:]]", "", x)
  vapply(x, function(s) {
    m <- regmatches(s, regexpr("HMDB[0-9]+", s))
    if (length(m) == 0) return(NA_character_)
    sprintf("HMDB%07d", as.integer(sub("HMDB", "", m)))
  }, character(1), USE.NAMES = FALSE)
}
nh <- norm_hmdb(lst$hmdb)
lst$hmdb <- ifelse(is.na(nh), lst$hmdb, nh)

# dedupe by (hmdb, adduct), keep first source
lst$key <- paste(lst$hmdb, lst$adduct)
lst <- lst[!duplicated(lst$key), ]

# ---- accurate theoretical m/z ----------------------------------------------
lst$exact_mass <- em[lst$hmdb]
lst$formula    <- fm[lst$hmdb]
lst$hmdb_name  <- nm[lst$hmdb]
ok_add <- lst$adduct %in% names(ADD)
lst$theo_mz <- NA_real_
for (i in which(ok_add & !is.na(lst$exact_mass))) {
  a <- ADD[[lst$adduct[i]]]; lst$theo_mz[i] <- a$nM * lst$exact_mass[i] + a$off
}
n_noemass <- sum(is.na(lst$exact_mass))
n_noadd   <- sum(!ok_add)
log_msg("Compounds: %d | no HMDB exact-mass: %d | unhandled adduct: %d",
        nrow(lst), n_noemass, n_noadd)
lst <- lst[!is.na(lst$theo_mz), ]
# keep only acquisition range
out_range <- lst$theo_mz < MZ_LO | lst$theo_mz > MZ_HI
log_msg("Dropped %d compounds outside m/z %d-%d (e.g. lactate [M-H]- 89).",
        sum(out_range), MZ_LO, MZ_HI)
lst <- lst[!out_range, ]

# ---- data grids -------------------------------------------------------------
log_msg("Loading feature grids...")
mse_full <- readRDS(file.path(CACHE_DIR, "peaks_after_freq.rds"))  # 61,124
mz_full  <- mz(mse_full)
mz_ana   <- mz(readRDS(file.path(CACHE_DIR, "peaks_combined.rds"))) # 8,772 analysis set

# ---- match within PPM_TOL ---------------------------------------------------
nearest <- function(q, grid) { j <- which.min(abs(grid - q)); c(idx=j, mz=grid[j]) }
res <- lst
res$data_mz <- NA_real_; res$delta_mDa <- NA_real_; res$delta_ppm <- NA_real_
res$score <- NA_real_; res$feat_idx <- NA_integer_; res$in_analysis_set <- FALSE
for (i in seq_len(nrow(res))) {
  nn <- nearest(res$theo_mz[i], mz_full)
  ppm <- (nn["mz"] - res$theo_mz[i]) / res$theo_mz[i] * 1e6
  if (abs(ppm) <= PPM_TOL) {
    res$feat_idx[i]   <- as.integer(nn["idx"])
    res$data_mz[i]    <- nn["mz"]
    res$delta_mDa[i]  <- (nn["mz"] - res$theo_mz[i]) * 1e3
    res$delta_ppm[i]  <- ppm
    res$score[i]      <- round(100 * (1 - abs(ppm)/PPM_TOL), 1)
    res$in_analysis_set[i] <- any(abs(mz_ana - nn["mz"]) < 1e-6)
  }
}
res$matched <- !is.na(res$feat_idx)

# ---- collisions: multiple compounds -> same data feature -------------------
res$isobaric_group <- NA_integer_
mres <- res[res$matched, ]
dup_feats <- mres$feat_idx[duplicated(mres$feat_idx) | duplicated(mres$feat_idx, fromLast=TRUE)]
res$isobaric_collision <- res$matched & res$feat_idx %in% unique(dup_feats)

# ---- write + summary --------------------------------------------------------
ord <- order(!res$matched, -res$score, res$theo_mz)
res <- res[ord, c("source","name","hmdb","hmdb_name","formula","adduct",
                  "exact_mass","theo_mz","published_mz","matched","data_mz",
                  "delta_mDa","delta_ppm","score","in_analysis_set",
                  "isobaric_collision","feat_idx")]
write.csv(res, file.path(OUT_DIR, "metabolite_match_table.csv"), row.names = FALSE)

cat(sprintf("\n==== SUMMARY ====\n"))
cat(sprintf("Compounds considered (in range, accurate mass): %d\n", nrow(res)))
cat(sprintf("MATCHED within %d ppm:                          %d\n", PPM_TOL, sum(res$matched)))
cat(sprintf("  ... also in 8,772 analysis set:               %d\n", sum(res$in_analysis_set)))
cat(sprintf("  ... isobaric collisions (shared feature):     %d\n", sum(res$isobaric_collision)))
cat(sprintf("Unmatched:                                      %d\n", sum(!res$matched)))
cat(sprintf("\nScore (ppm proximity) distribution of matches:\n"))
mm <- res[res$matched,]
cat(sprintf("  |ppm|<2 (score>90): %d | <5 (>75): %d | <10 (>50): %d | <20: %d\n",
    sum(abs(mm$delta_ppm)<2), sum(abs(mm$delta_ppm)<5), sum(abs(mm$delta_ppm)<10), nrow(mm)))
cat(sprintf("\nWrote: %s\n", file.path(OUT_DIR, "metabolite_match_table.csv")))
