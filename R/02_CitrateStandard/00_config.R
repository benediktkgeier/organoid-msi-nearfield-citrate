#!/usr/bin/env Rscript
# 02_CitrateStandard / 00_config.R
# Shared config for the citrate-standard calibration + ID-confidence step.
#
# An authentic Sodium Citrate Tribasic Dihydrate dilution series was spotted in
# the SAME CMC matrix (4% in 150 uM ammonium bicarbonate) on the SAME slide as
# the organoid sections and sprayed with the SAME matrix. Because the spots are
# pure citrate in CMC (NO tissue), there is NO tissue co-isobar -> the standard
# pins the true on-instrument citrate mass / peak shape that the tissue feature
# (191.0217, +8 ppm, blended; see R/07_metabolite_id/03_citrate_resolution.R)
# cannot reveal on its own.
#
# Sourced by 01_/02_/03_. Defines: paths, spot<->concentration table, citrate +
# adduct masses (+ the new [2M-H]- dimer), the verbatim low-level imzML/.ibd
# reader, and the tissue / same-run organoid hook for standard-vs-tissue compare.

# lib_citrate.R provides the low-level imzML reader (parse_imzml/read_spots/
# read_file/mz_win/ppm_of, WIN_PPM, HIST_BY) and sources lib_paths.R.
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_citrate.R"))

# ---- standard data location ------------------------------------------------
STD_DIR <- file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"MSI/06102026_AO_0h_sl7A/imzml")

# spot <-> concentration table, in ascending concentration.  The acquisition
# filenames had been mislabeled (10<->100 uM and 10<->100 mM swapped); the raw
# .imzML/.ibd files were PHYSICALLY RENAMED to their TRUE concentration on
# 2026-06-25 (see MSI/06102026_AO_0h_sl7A/imzml/RENAME_LOG.txt).  So filename ==
# true concentration now and 'file' maps directly.  0um = citrate-free CMC blank.
SPOTS <- data.frame(
  label    = c("blank", "10 uM", "100 uM", "1 mM", "10 mM", "100 mM"),
  file     = c("06102026_ao_0h_sl7a_0um.imzML",  "06102026_ao_0h_sl7a_10um.imzML",
               "06102026_ao_0h_sl7a_100um.imzML", "06102026_ao_0h_sl7a_1mm.imzML",
               "06102026_ao_0h_sl7a_10mm.imzML",  "06102026_ao_0h_sl7a_100mm.imzML"),
  conc_M   = c(0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
  file_conc= c("0um", "10um", "100um", "1mm", "10mm", "100mm"),
  is_blank = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

# provenance footnote (acquisition mislabel was corrected on disk)
LABEL_NOTE <- "Acquisition filenames had 10<->100 uM and 10<->100 mM swapped; raw files physically renamed to TRUE concentration on 2026-06-25 (see RENAME_LOG.txt)."

# QC: with the corrected mapping the dose-response must increase with conc.
# returns list(monotonic = TRUE/FALSE, order = labels in ascending response).
verify_monotonic <- function(conc, resp) {
  k <- conc > 0 & is.finite(resp)
  o <- order(conc[k])
  list(monotonic = all(diff(resp[k][o]) >= 0), resp_sorted = resp[k][o])
}

# ---- masses (verbatim from R/07_metabolite_id/05_citrate_isotopes_adducts.R) -
C12   <- 191.019726                 # citrate [M-H]-  C6H7O7-   (locked project value)
DC13  <- 1.0033548                  # 13C - 12C
PROT  <- 1.00727646; ELEC <- 0.000548580
H     <- 1.00782503; NA_  <- 22.98976928; CL <- 34.96885268
CH3COO<- 59.013304                  # acetate C2H3O2
Mneu  <- C12 + PROT                 # neutral citrate = 192.027002

# adduct / isotope ion table (same construction as script 05) + [2M-H]- dimer.
# rel_to is the denominator ion for the abundance-ratio fingerprint (03_).
IONS <- list(
  list(key="MH",     lab="Citrate [M-H]- (C12)",  form="C6H7O7-",     mz = C12),
  list(key="C13",    lab="13C isotope (C13)",      form="13C-C5H7O7-", mz = C12 + DC13),
  list(key="Na",     lab="[M-2H+Na]-",             form="C6H6O7Na-",   mz = Mneu - 2*H + NA_   + ELEC),
  list(key="Na2",    lab="[M-3H+2Na]-",            form="C6H5O7Na2-",  mz = Mneu - 3*H + 2*NA_ + ELEC),
  list(key="Cl",     lab="[M+Cl]-",                form="C6H8O7Cl-",   mz = Mneu + CL     + ELEC),
  list(key="OAc",    lab="[M+CH3COO]- (acetate)",  form="C8H11O9-",    mz = Mneu + CH3COO + ELEC),
  list(key="dimer",  lab="[2M-H]- dimer",          form="C12H15O14-",  mz = 2*Mneu - PROT)
)
ION_MZ  <- sapply(IONS, `[[`, "mz"); names(ION_MZ) <- sapply(IONS, `[[`, "key")
MH_MZ   <- C12
DIMER_MZ<- ION_MZ[["dimer"]]
C13_MZ  <- ION_MZ[["C13"]]

# theoretical 13C M+1/M for citrate: 6 carbons dominate.  13C natural abundance
# 1.07% -> carbon term 6 * 0.0107/0.9893 = 6.49%; +17O / 2H terms ~0.3% -> ~6.8%.
ISO_THEO_C  <- 6 * 0.0107 / 0.9893      # carbon-only ~0.0649
ISO_THEO    <- 0.068                    # full envelope ~6.8% (label as "~6.5-7%")

# empty m/z window (equal width, far from any citrate ion) = per-pixel noise floor
NOISE_WIN_C <- 190.9459                 # centre of the empty window used in script 03
# NOTE: mz_win/ppm_of/WIN_PPM/HIST_BY + the imzML reader now live in lib_citrate.R

# pretty-print a molar concentration (used in figure annotations)
fmtM <- function(M) {
  if (!is.finite(M)) return("NA")
  unit <- function(v, u) paste0(formatC(v, format="g", digits=3), " ", u)
  if (M >= 1e-3) unit(M*1e3, "mM")
  else if (M >= 1e-6) unit(M*1e6, "uM")
  else unit(M*1e9, "nM")
}

# display label for a spot from its (swap-corrected) concentration
spot_label <- function(conc) if (!is.finite(conc) || conc == 0) "blank" else fmtM(conc)

# ---- render constants (locked, same as Phase 06) ---------------------------
SCALE_UM <- 500; PX_UM <- MSI_PIXEL_UM
OUT_FIG  <- file.path(FIG_DIR, "citrate_standard")
OUT_RES  <- file.path(RES_DIR, "citrate_standard")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_RES, showWarnings = FALSE, recursive = TRUE)

# ---- tissue comparison (standard is the matrix-matched reference) ----------
# EXPERIMENTAL DESIGN (user, 2026-06-25): on-sample / same-run co-acquisition of
# organoids with the standard is NOT possible. The spotted standard stands as the
# reference because it was spotted in the same CMC matrix, on the SAME GLASS SLIDE
# as organoid sections (sl7A), and sprayed identically -> matrix/preparation-
# matched to the organoid measurements.
#
# TISSUE_SECTIONS = the 5 genuine-citrate organoid sections available as imzML
# (sl6A/sl4A): same study + same matrix/spray prep, but a SEPARATE acquisition run
# and different slide, so absolute mass-cal / detector sensitivity can differ
# slightly (this is the only residual gap; the matrix is matched).
#
# LITERAL SAME-SLIDE option: sl7A also carries organoid sections, but their MSI
# data lives only in the SCiLS project (..._sections.sbd) -- not as imzML/.d. To
# compare on the SAME SLIDE as the standard, export those sl7A sections from SCiLS
# Lab to imzML and list them in SAMERUN_SECTIONS below (dir = STD_DIR); scripts
# 01 & 03 then prefer them and relabel "sep-run" -> "same-slide". No other change.
TISSUE_DIR <- file.path(CACHE_DIR, "imzml")
TISSUE_SECTIONS <- list(
  list(sid="0h_sl6A_sec2a",  grp="0h",  imz="06102026_ao_0h_sl6a_sec2a.imzML"),
  list(sid="0h_sl6A_sec5a",  grp="0h",  imz="06102026_ao_0h_sl6a_sec5a.imzML"),
  list(sid="20h_sl4A_sec2a", grp="20h", imz="06102026_ao_20h_sl4a_2a.imzML"),
  list(sid="20h_sl4A_sec3a", grp="20h", imz="06102026_ao_20h_sl4a_3a.imzML"),
  list(sid="20h_sl4A_sec3b", grp="20h", imz="06102026_ao_20h_sl4a_3b.imzML")
)
SAMERUN_SECTIONS <- list()        # <- optional: sl7A organoid sections (SCiLS->imzML)

# returns list(dir=, sections=, mode="same-slide"|"sep-run") -----------------
active_tissue <- function() {
  if (length(SAMERUN_SECTIONS)) {
    list(dir = STD_DIR, sections = SAMERUN_SECTIONS, mode = "same-slide")
  } else {
    list(dir = TISSUE_DIR, sections = TISSUE_SECTIONS, mode = "sep-run")
  }
}
TISSUE_CAVEAT <- paste(
  "Tissue side is %s.",
  "The standard is matrix-matched to the organoids (same CMC, same slide, same",
  "spray); on-sample co-acquisition was not possible. sep-run = organoid sections",
  "from a separate acquisition run/slide, so absolute mass-cal/sensitivity can",
  "differ slightly (matrix is matched). Export sl7A SCiLS sections to imzML for a",
  "literal same-slide comparison."
)

viridis_pal <- function() viridisLite::viridis(256)
log_cfg <- function(tag, ...) message(sprintf("[%s] %s", tag, sprintf(...)))
