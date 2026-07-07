# ============================================================================
# 12_joy_tables.R
# Focused hand-off subset: a "Joy_Tables" folder with the Excel workbook + a
# single Prism project for three panels of apical_citrate_dha_report.pdf:
#   Page 1   Per-organoid mean ion intensity by apical orientation   (cit, dha)
#   Page 8   Absolute near-field citrate level in the gel            (near50, near100)
#   Page 8b  Near-field citrate 0-50 vs 0-100 um WITHIN each class   (paired)
# All are Column / paired tables, so one .pzfx suffices (no XY).
# READ-ONLY w.r.t. analysis.
# ============================================================================

stopifnot(requireNamespace("writexl", quietly = TRUE), requireNamespace("pzfx", quietly = TRUE))

ROOT      <- file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final")
ANNOT_RES <- file.path(ROOT, "results", "annotation")
OUT_DIR   <- file.path(ANNOT_RES, "prism", "Joy_Tables")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

org  <- read.csv(file.path(ANNOT_RES, "apical_citrate_dha_per_organoid_normalized.csv"), stringsAsFactors = FALSE)
grad <- read.csv(file.path(ANNOT_RES, "apical_gradient_per_organoid.csv"), stringsAsFactors = FALSE)
stat <- read.csv(file.path(ANNOT_RES, "apical_citrate_dha_stats.csv"), stringsAsFactors = FALSE)

CLASSES <- c("basolateral_out", "apical_out", "mixed")
CLAB    <- c(basolateral_out = "Basolateral-out", apical_out = "Apical-out", mixed = "Mixed")

wide_by_class <- function(df, valcol) {
  cols <- lapply(CLASSES, function(cl) { v <- df[[valcol]][df$apical_class == cl]; v[is.finite(v)] })
  n <- max(vapply(cols, length, 0L))
  out <- lapply(cols, function(v) { length(v) <- n; v })
  setNames(as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE), CLAB[CLASSES])
}
paired_class <- function(df, cl, keep_id = TRUE) {
  d <- df[df$apical_class == cl, c("sid", "instance", "near50", "near100")]
  d <- d[is.finite(d$near50) & is.finite(d$near100), , drop = FALSE]
  base <- data.frame(`0-50 um` = d$near50, `0-100 um` = d$near100,
                     check.names = FALSE, stringsAsFactors = FALSE)
  if (keep_id) cbind(id = paste(d$sid, d$instance, sep = " | "), base) else base
}

# ---- stats reference (only the three figures) ------------------------------
tab_of <- c(cit = "P1_Citrate_abs", dha = "P1_DHA_abs",
            grad_near50 = "P8_near50", grad_near100 = "P8_near100")
sref <- stat[stat$metric %in% names(tab_of) & stat$scope == "all", ]
sref <- sref[order(match(sref$metric, names(tab_of))), ]
stats_ref <- data.frame(
  tab = unname(tab_of[sref$metric]), metric = sref$metric,
  comparison = paste(CLAB[sref$g1], "vs", CLAB[sref$g2]),
  n1 = sref$n1, n2 = sref$n2, p_value = sref$p, stars = sref$stars,
  prism_test = "Column -> unpaired Mann-Whitney (per pair); or Kruskal-Wallis + Dunn",
  check.names = FALSE, stringsAsFactors = FALSE)
pcount <- function(cl) sum(is.finite(grad$near50[grad$apical_class==cl]) &
                           is.finite(grad$near100[grad$apical_class==cl]))
stats_ref <- rbind(stats_ref, data.frame(
  tab = c("P8b_basolateral_out","P8b_apical_out","P8b_mixed"),
  metric = "near50 vs near100 (paired)", comparison = "0-50 um vs 0-100 um within class",
  n1 = c(pcount("basolateral_out"), pcount("apical_out"), pcount("mixed")),
  n2 = NA_integer_, p_value = NA_real_, stars = "",
  prism_test = "Column -> paired Wilcoxon matched-pairs",
  check.names = FALSE, stringsAsFactors = FALSE))

# ---- README tab ------------------------------------------------------------
readme <- data.frame(field=character(), value=character(), check.names=FALSE, stringsAsFactors=FALSE)
add <- function(k,v) readme[nrow(readme)+1L, ] <<- list(k,v)
add("Source report", "figures/annotation/apical_citrate_dha_report.pdf (pages 1, 8, 8b)")
add("Citrate m/z", "191.01976 [M-H]-  (anchored, +/-7 ppm)")
add("DHA m/z", "327.2330 [M-H]- (C22:6)")
add("Units", "TIC-normalized mean intensity, a.u.")
add("Replication unit", "one organoid")
add("n per class", sprintf("Basolateral-out=%d, Apical-out=%d, Mixed=%d",
    sum(org$apical_class=="basolateral_out"), sum(org$apical_class=="apical_out"), sum(org$apical_class=="mixed")))
add("", "")
add("P1_Citrate_abs / P1_DHA_abs", "PAGE 1: per-organoid mean citrate / DHA by apical class (Column)")
add("P8_near50 / P8_near100", "PAGE 8: absolute near-field citrate, gel 0-50 / 0-100 um (Column)")
add("P8b_*", "PAGE 8b: paired near50 vs near100 within each class (Column, before-after)")
add("Stats_reference", "report p-values + which Prism test to use")
add("", "")
add("Prism (P1, P8)", "Column table, scatter+box; Mann-Whitney per pair (or Kruskal-Wallis+Dunn)")
add("Prism (P8b)", "Column table, before-after; Wilcoxon matched-pairs (paired by row; id = row title)")
add("Blank cells", "unequal n per class (23/30/31); blanks are padding, not zeros")

# ---- Excel workbook --------------------------------------------------------
sheets <- list(
  README             = readme,
  P1_Citrate_abs     = wide_by_class(org,  "cit"),
  P1_DHA_abs         = wide_by_class(org,  "dha"),
  P8_near50          = wide_by_class(grad, "near50"),
  P8_near100         = wide_by_class(grad, "near100"),
  P8b_basolateral_out = paired_class(grad, "basolateral_out"),
  P8b_apical_out      = paired_class(grad, "apical_out"),
  P8b_mixed           = paired_class(grad, "mixed"),
  Stats_reference    = stats_ref
)
xlsx <- file.path(OUT_DIR, "Joy_Tables.xlsx")
writexl::write_xlsx(sheets, xlsx)

# ---- Prism project (all Column/paired -> one file; no id text col) ---------
pz_tables <- list(
  P1_Citrate_abs     = wide_by_class(org,  "cit"),
  P1_DHA_abs         = wide_by_class(org,  "dha"),
  P8_near50          = wide_by_class(grad, "near50"),
  P8_near100         = wide_by_class(grad, "near100"),
  P8b_basolateral_out = paired_class(grad, "basolateral_out", keep_id = FALSE),
  P8b_apical_out      = paired_class(grad, "apical_out",      keep_id = FALSE),
  P8b_mixed           = paired_class(grad, "mixed",           keep_id = FALSE)
)
pzf <- file.path(OUT_DIR, "Joy_Tables.pzfx")
pzfx::write_pzfx(pz_tables, path = pzf, row_names = FALSE)

cat(sprintf("[joy] xlsx: %d tabs -> %s\n", length(sheets), xlsx))
cat(sprintf("[joy] pzfx: %d tables -> %s\n", length(pz_tables), pzf))