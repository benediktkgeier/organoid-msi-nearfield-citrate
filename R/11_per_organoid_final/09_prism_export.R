# ============================================================================
# 09_prism_export.R
# Reshape ALL reproducible quantitative panels of apical_citrate_dha_report.pdf
# into a GraphPad Prism-ready Excel workbook (one tab per Prism data table).
# READ-ONLY w.r.t. analysis: only writes apical_prism_tables.xlsx.
#
# Report page -> tab(s):
#   pg1  absolute TIC means          cit, dha                 (Column, by class)
#   pg2  metabolite-TIC normalised   cit_mtic, dha_mtic       (Column, by class)
#   pg3  within-section normalised   cit_rel, dha_rel         (Column, by class)
#   pg4  citrate/DHA ratio           cit_dha                  (Column, by class)
#   pg6  outward decay curves        surf/abs/mtic/citdha     (XY, distance vs value)
#   pg7  gradient metrics            rho_out, far_index       (Column, by class)
#   pg8  near-field level            near50, near100          (Column, by class)
#   pg8b near50 vs near100 paired    per class                (Column, paired)
#   pg0  summary text  -> Stats_reference tab
#   pg9+ spatial BF+MSI overlays     -> NOT reproducible in Prism (raster images)
# Prereq: run 11_export_gradient_profile.R first (writes the pg6 long CSV).
# ============================================================================

ROOT      <- file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final")
ANNOT_RES <- file.path(ROOT, "results", "annotation")
OUT_DIR   <- file.path(ANNOT_RES, "prism")
FIG_DIR   <- file.path(ROOT, "figures", "annotation")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

F_ORG   <- file.path(ANNOT_RES, "apical_citrate_dha_per_organoid_normalized.csv")
F_GRAD  <- file.path(ANNOT_RES, "apical_gradient_per_organoid.csv")
F_STATS <- file.path(ANNOT_RES, "apical_citrate_dha_stats.csv")
F_PROF  <- file.path(ANNOT_RES, "apical_gradient_profile_long.csv")   # from 11_export_gradient_profile.R
stopifnot(file.exists(F_ORG), file.exists(F_GRAD), file.exists(F_STATS))
have_prof <- file.exists(F_PROF)

org  <- read.csv(F_ORG,   stringsAsFactors = FALSE)
grad <- read.csv(F_GRAD,  stringsAsFactors = FALSE)
stat <- read.csv(F_STATS, stringsAsFactors = FALSE)
prof <- if (have_prof) read.csv(F_PROF, stringsAsFactors = FALSE) else NULL

CLASSES <- c("basolateral_out", "apical_out", "mixed")
CLAB    <- c(basolateral_out = "Basolateral-out", apical_out = "Apical-out", mixed = "Mixed")

# ---- helpers ---------------------------------------------------------------
wide_by_class <- function(df, valcol) {          # Column table: one col per class
  cols <- lapply(CLASSES, function(cl) { v <- df[[valcol]][df$apical_class == cl]; v[is.finite(v)] })
  n <- max(vapply(cols, length, 0L))
  out <- lapply(cols, function(v) { length(v) <- n; v })
  setNames(as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE), CLAB[CLASSES])
}
paired_class <- function(df, cl) {               # paired near50/near100 within a class
  d <- df[df$apical_class == cl, c("sid", "instance", "near50", "near100")]
  d <- d[is.finite(d$near50) & is.finite(d$near100), , drop = FALSE]
  data.frame(id = paste(d$sid, d$instance, sep = " | "),
             `0-50 um` = d$near50, `0-100 um` = d$near100,
             check.names = FALSE, stringsAsFactors = FALSE)
}
# pg6 XY table: X = zone_um, group-median curves + per-organoid spaghetti
decay_xy <- function(gp, metric, with_dha = FALSE) {
  zones <- sort(unique(gp$zone_um))
  med <- function(cl) vapply(zones, function(z)
    stats::median(gp[[metric]][gp$zone_um == z & gp$apical_class == cl], na.rm = TRUE), 0.0)
  out <- data.frame(`distance_um` = zones,
                    `med_Basolateral-out` = med("basolateral_out"),
                    `med_Apical-out`      = med("apical_out"),
                    check.names = FALSE)
  if (with_dha)                                   # DHA control = median over ALL organoids
    out[["med_DHA-control"]] <- vapply(zones, function(z)
      stats::median(gp$dha_sn[gp$zone_um == z], na.rm = TRUE), 0.0)
  # per-organoid spaghetti (basolateral then apical, as plotted)
  sub <- gp[gp$apical_class %in% c("basolateral_out", "apical_out"), ]
  keys <- unique(sub[order(match(sub$apical_class, CLASSES), sub$key), "key"])
  for (kk in keys) {
    s <- sub[sub$key == kk, ]; cl <- s$apical_class[1]
    v <- s[[metric]][match(zones, s$zone_um)]
    out[[sprintf("%s | %s", CLAB[cl], kk)]] <- v
  }
  out
}

# ---- Column-table specs (report page -> tab) -------------------------------
col_specs <- list(
  list(tab="P1_Citrate_abs",    df=org,  col="cit"),
  list(tab="P1_DHA_abs",        df=org,  col="dha"),
  list(tab="P2_Citrate_metaTIC",df=org,  col="cit_mtic"),
  list(tab="P2_DHA_metaTIC",    df=org,  col="dha_mtic"),
  list(tab="P3_Citrate_secnorm",df=org,  col="cit_rel"),
  list(tab="P3_DHA_secnorm",    df=org,  col="dha_rel"),
  list(tab="P4_CitDHA_ratio",   df=org,  col="cit_dha"),
  list(tab="P7_rho_out",        df=grad, col="rho_out"),
  list(tab="P7_far_index",      df=grad, col="far_index"),
  list(tab="P8_near50",         df=grad, col="near50"),
  list(tab="P8_near100",        df=grad, col="near100")
)

# ---- Stats reference -------------------------------------------------------
tab_of <- c(cit="P1_Citrate_abs", dha="P1_DHA_abs", cit_mtic="P2_Citrate_metaTIC",
            dha_mtic="P2_DHA_metaTIC", cit_rel="P3_Citrate_secnorm", dha_rel="P3_DHA_secnorm",
            cit_dha="P4_CitDHA_ratio", grad_rho_out="P7_rho_out", grad_far_index="P7_far_index",
            grad_near50="P8_near50", grad_near100="P8_near100")
sref <- stat[stat$metric %in% names(tab_of) & stat$scope == "all", ]
sref <- sref[order(match(sref$metric, names(tab_of))), ]
stats_ref <- data.frame(
  tab = unname(tab_of[sref$metric]), metric = sref$metric,
  comparison = paste(CLAB[sref$g1], "vs", CLAB[sref$g2]),
  n1 = sref$n1, n2 = sref$n2, p_value = sref$p, stars = sref$stars,
  prism_test = "Column -> unpaired Mann-Whitney (per pair); or Kruskal-Wallis + Dunn across 3",
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
add("Source report", "figures/annotation/apical_citrate_dha_report.pdf")
add("Generated by", "R/11_per_organoid_final/09_prism_export.R (+ 11_export_gradient_profile.R for pg6)")
add("Citrate m/z", "191.01976 [M-H]- (anchored, +/-7 ppm, raw imzML)")
add("DHA m/z", "327.2330 [M-H]- (C22:6)")
add("Units", "TIC-normalized mean intensity, a.u.; *_rel = within-section log2; *_mtic = % of 341-feature pool")
add("Replication unit", "one organoid (per-pixel means aggregated per organoid)")
add("n per class", sprintf("Basolateral-out=%d, Apical-out=%d, Mixed=%d",
    sum(org$apical_class=="basolateral_out"), sum(org$apical_class=="apical_out"), sum(org$apical_class=="mixed")))
add("Class keys", "basolateral_out=Basolateral-out; apical_out=Apical-out; mixed=Mixed")
add("", "")
add("TAB -> report page / meaning", "")
add("P1_Citrate_abs / P1_DHA_abs", "pg1: per-organoid mean citrate / DHA (absolute TIC)")
add("P2_Citrate_metaTIC / P2_DHA_metaTIC", "pg2: ion / metabolite-TIC pool (% of pool)")
add("P3_Citrate_secnorm / P3_DHA_secnorm", "pg3: within-section log2( organoid / section median )")
add("P4_CitDHA_ratio", "pg4: log2( citrate / DHA ) internal-control ratio")
add("P6_surf / P6_abs / P6_mtic / P6_citdha", "pg6: outward decay curves (XY: distance vs value); medians + spaghetti")
add("P7_rho_out / P7_far_index", "pg7: outward monotonicity (Spearman rho) / far-field retention")
add("P8_near50 / P8_near100", "pg8: absolute near-field citrate 0-50 / 0-100 um")
add("P8b_*", "pg8b: paired near50 vs near100 within each apical class")
add("Stats_reference", "pg0 summary: report p-values + which Prism test to use")
add("", "")
add("Prism (Column tabs)", "Column table, scatter+box; Mann-Whitney per pair (or Kruskal-Wallis+Dunn)")
add("Prism (P8b tabs)", "Column table, before-after; Wilcoxon matched-pairs (paired by row)")
add("Prism (P6 tabs)", "XY table; col1=distance_um (X, log axis); med_* = bold group curves; rest = per-organoid spaghetti")
add("pg6 note", "pg6 uses basolateral-out + apical-out only (mixed omitted); DHA-control line only on P6_surf")
add("Not included", "pg9+ per-section brightfield+MSI overlay IMAGES are spatial rasters - not a Prism plot type")
add("Blank cells", "columns have unequal n; blanks are padding, not zeros")

# ---- assemble workbook -----------------------------------------------------
sheets <- list(README = readme)
for (s in col_specs) sheets[[s$tab]] <- wide_by_class(s$df, s$col)
sheets[["P8b_basolateral_out"]] <- paired_class(grad, "basolateral_out")
sheets[["P8b_apical_out"]]      <- paired_class(grad, "apical_out")
sheets[["P8b_mixed"]]           <- paired_class(grad, "mixed")
if (have_prof) {
  sheets[["P6_surf"]]   <- decay_xy(prof, "surf", with_dha = TRUE)
  sheets[["P6_abs"]]    <- decay_xy(prof, "abs")
  sheets[["P6_mtic"]]   <- decay_xy(prof, "mtic")
  sheets[["P6_citdha"]] <- decay_xy(prof, "citdha")
} else {
  cat("[prism] WARNING: pg6 long CSV missing; run 11_export_gradient_profile.R. Skipping P6_* tabs.\n")
}
sheets[["Stats_reference"]] <- stats_ref

out_xlsx <- file.path(OUT_DIR, "apical_prism_tables.xlsx")
fig_xlsx <- file.path(FIG_DIR, "apical_prism_tables.xlsx")
writer <- NULL
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(sheets, out_xlsx); writer <- "writexl"
} else if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(sheets, out_xlsx); writer <- "openxlsx"
} else {
  for (nm in names(sheets)) write.csv(sheets[[nm]], file.path(OUT_DIR, paste0(nm, ".csv")), row.names = FALSE)
  writer <- "csv-fallback"
}
if (writer != "csv-fallback") invisible(file.copy(out_xlsx, fig_xlsx, overwrite = TRUE))

cat(sprintf("[prism] writer=%s  tabs=%d\n", writer, length(sheets)))
cat(sprintf("[prism] %s\n", paste(names(sheets), collapse = ", ")))
cat(sprintf("[prism] out = %s\n", if (writer=="csv-fallback") OUT_DIR else out_xlsx))
