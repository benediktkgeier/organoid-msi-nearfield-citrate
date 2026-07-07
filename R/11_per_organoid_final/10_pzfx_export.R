# ============================================================================
# 10_pzfx_export.R
# Native GraphPad Prism projects (.pzfx) for the apical citrate/DHA report.
# Two files (a single .pzfx cannot mix Column and XY table types via write_pzfx):
#   apical_prism_tables.pzfx  -> Column + paired tables (pg1-4, 7, 8, 8b)
#   apical_gradient_decay.pzfx-> XY tables, outward decay curves (pg6)
# .pzfx tables are numeric-only: README / Stats_reference live in the .xlsx.
# Prereq: run 11_export_gradient_profile.R first (writes the pg6 long CSV).
# READ-ONLY w.r.t. analysis.
# ============================================================================

stopifnot(requireNamespace("pzfx", quietly = TRUE))

ROOT      <- file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final")
ANNOT_RES <- file.path(ROOT, "results", "annotation")
OUT_DIR   <- file.path(ANNOT_RES, "prism")
FIG_DIR   <- file.path(ROOT, "figures", "annotation")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

org  <- read.csv(file.path(ANNOT_RES, "apical_citrate_dha_per_organoid_normalized.csv"), stringsAsFactors = FALSE)
grad <- read.csv(file.path(ANNOT_RES, "apical_gradient_per_organoid.csv"), stringsAsFactors = FALSE)
F_PROF   <- file.path(ANNOT_RES, "apical_gradient_profile_long.csv")
have_prof <- file.exists(F_PROF)
prof <- if (have_prof) read.csv(F_PROF, stringsAsFactors = FALSE) else NULL

CLASSES <- c("basolateral_out", "apical_out", "mixed")
CLAB    <- c(basolateral_out = "Basolateral-out", apical_out = "Apical-out", mixed = "Mixed")

wide_by_class <- function(df, valcol) {
  cols <- lapply(CLASSES, function(cl) { v <- df[[valcol]][df$apical_class == cl]; v[is.finite(v)] })
  n <- max(vapply(cols, length, 0L))
  out <- lapply(cols, function(v) { length(v) <- n; v })
  setNames(as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE), CLAB[CLASSES])
}
paired_class <- function(df, cl) {
  d <- df[df$apical_class == cl, c("near50", "near100")]
  d <- d[is.finite(d$near50) & is.finite(d$near100), , drop = FALSE]
  setNames(data.frame(d$near50, d$near100, check.names = FALSE), c("0-50 um", "0-100 um"))
}
decay_xy_med <- function(gp, metric, with_dha = FALSE) {   # medians only (bold curves)
  zones <- sort(unique(gp$zone_um))
  med <- function(cl) vapply(zones, function(z)
    stats::median(gp[[metric]][gp$zone_um == z & gp$apical_class == cl], na.rm = TRUE), 0.0)
  out <- data.frame(distance_um = zones,
                    `Basolateral-out` = med("basolateral_out"),
                    `Apical-out`      = med("apical_out"), check.names = FALSE)
  if (with_dha) out[["DHA-control"]] <- vapply(zones, function(z)
    stats::median(gp$dha_sn[gp$zone_um == z], na.rm = TRUE), 0.0)
  out
}

# ---- file 1: Column + paired tables ----------------------------------------
col_tables <- list(
  P1_Citrate_abs     = wide_by_class(org,  "cit"),
  P1_DHA_abs         = wide_by_class(org,  "dha"),
  P2_Citrate_metaTIC = wide_by_class(org,  "cit_mtic"),
  P2_DHA_metaTIC     = wide_by_class(org,  "dha_mtic"),
  P3_Citrate_secnorm = wide_by_class(org,  "cit_rel"),
  P3_DHA_secnorm     = wide_by_class(org,  "dha_rel"),
  P4_CitDHA_ratio    = wide_by_class(org,  "cit_dha"),
  P7_rho_out         = wide_by_class(grad, "rho_out"),
  P7_far_index       = wide_by_class(grad, "far_index"),
  P8_near50          = wide_by_class(grad, "near50"),
  P8_near100         = wide_by_class(grad, "near100"),
  P8b_basolateral_out = paired_class(grad, "basolateral_out"),
  P8b_apical_out      = paired_class(grad, "apical_out"),
  P8b_mixed           = paired_class(grad, "mixed")
)
f1 <- file.path(OUT_DIR, "apical_prism_tables.pzfx")
pzfx::write_pzfx(col_tables, path = f1, row_names = FALSE)
invisible(file.copy(f1, file.path(FIG_DIR, "apical_prism_tables.pzfx"), overwrite = TRUE))
cat(sprintf("[pzfx] file1: %d Column/paired tables -> %s\n", length(col_tables), f1))

# ---- file 2: XY decay tables (pg6) -----------------------------------------
if (have_prof) {
  xy_tables <- list(
    P6_surf   = decay_xy_med(prof, "surf", with_dha = TRUE),
    P6_abs    = decay_xy_med(prof, "abs"),
    P6_mtic   = decay_xy_med(prof, "mtic"),
    P6_citdha = decay_xy_med(prof, "citdha")
  )
  f2 <- file.path(OUT_DIR, "apical_gradient_decay.pzfx")
  pzfx::write_pzfx(xy_tables, path = f2, row_names = FALSE, x_col = 1L)   # col1 = X (distance)
  invisible(file.copy(f2, file.path(FIG_DIR, "apical_gradient_decay.pzfx"), overwrite = TRUE))
  cat(sprintf("[pzfx] file2: %d XY decay tables (medians; spaghetti in .xlsx) -> %s\n", length(xy_tables), f2))
} else {
  cat("[pzfx] WARNING: pg6 long CSV missing; run 11_export_gradient_profile.R. Skipping decay .pzfx.\n")
}
