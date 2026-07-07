#!/usr/bin/env Rscript
# ============================================================================
# if_config.R - registry + constants for IF (B-section) -> BF/MSI registration.
#   Sourced by R/90..93. Builds on gradient_config.R (paths) + lib_register*.R.
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))

IF_DIR    <- file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"LM_Bsections")
IF_CACHE  <- cache_in("register_if"); dir.create(IF_CACHE, showWarnings = FALSE, recursive = TRUE)
IF_FIG    <- file.path(FIG_DIR,  "if_registration"); dir.create(file.path(IF_FIG, "crops"), showWarnings = FALSE, recursive = TRUE)
IF_RES    <- file.path(RES_DIR,  "if_registration"); dir.create(IF_RES, showWarnings = FALSE, recursive = TRUE)
REG_CACHE <- cache_in("register")   # existing MSI<->BF transforms live here (read-only reuse)

# ---- pixel scales (um/px) --------------------------------------------------
OV_UMPX <- 1.833333   # overview .nd2 (6x) == MSI brightfield .nd2 scale
HR_UMPX <- 0.361835   # hi-res .nd2 (30x)
BF_UMPX <- 1.833333   # MSI brightfield .nd2 (6x)
HR_OV_RATIO <- OV_UMPX / HR_UMPX           # ~5.066 (hi-res px per overview px... downscale factor)

# ---- thumbnail block-mean factors ------------------------------------------
F_OV  <- 16L   # overview thumb: 1.833*16 = 29.3 um/px  (matches BF nd2thumb F=16)
F_BF  <- 16L   # BF nd2thumb_<slide>.rds is F=16
F_HR  <- 5L    # hi-res DAPI thumb: 0.362*5 = 1.81 um/px ~ overview native
# common coarse scale for hi-res<->overview windowed match (~7.3 um/px):
F_HR_C <- 20L  # 0.362*20 = 7.24 um/px
F_OV_C <- 4L   # 1.833*4  = 7.33 um/px

# ACTUAL hi-res .nd2 channel order (from NIS-Elements metadata, emission wl):
#   ch1 = Cy5 / F-actin (670), ch2 = mCherry / ZO-1 (642),
#   ch3 = FITC / b-catenin (524), ch4 = DAPI / DNA (405).
# (This is REVERSED from the earlier assumption; corrected 2026-06-18.)
# NOTE: the whole-slide OVERVIEW .nd2 is single-channel DAPI, so it uses ch=1.
IF_CH <- c(Factin = 1L, ZO1 = 2L, bcat = 3L, DAPI = 4L)   # hi-res channel index per marker
OV_DAPI_CH      <- 1L    # overview is single-channel DAPI
DAPI_CH_DEFAULT <- 1L    # (legacy) overview/registration thumbnails; hi-res DAPI = IF_CH["DAPI"]
IF_MARK_COL <- c(DAPI = "cyan", ZO1 = "red", bcat = "green", Factin = "white")  # DAPI cyan: better contrast on black

# ---- native NIS-Elements display LUTs per MARKER, per B-slide ----
# From the user's NIS LUT sliders (lo, hi raw counts; gamma 1.0).
# Display: clamp((raw - lo)/(hi - lo), 0, 1) ^ (1/gamma), then * marker colour.
IF_LUT <- list(
  sl6b = list(DAPI=c(lo=148,hi=4095,g=1), ZO1=c(lo=413,hi=2663,g=1),
              bcat=c(lo=420,hi=3324,g=1), Factin=c(lo=132,hi=3184,g=1)),
  sl4b = list(DAPI=c(lo=210,hi=2546,g=1), ZO1=c(lo=459,hi=2904,g=1),
              bcat=c(lo=335,hi=2771,g=1), Factin=c(lo=140,hi=3083,g=1))
)

# ---- slide-level registry (overviews) --------------------------------------
# msi_slide / msi_bf = the A-slide this B-slide is serial to (the MSI anchor).
IF_SLIDES <- data.frame(
  slide     = c("sl6b", "sl4b"),
  group     = c("0h",   "20h"),
  over      = file.path(IF_DIR, c("06172026_AO_0h_sl6b_over.nd2", "06172026_AO_20h_sl4b_over.nd2")),
  hr_prefix = c("06172026_AO_0h_sl6b_sec", "06172026_AO_20h_sl4_sec"),  # note: 20h hi-res files are 'sl4_'
  sid_prefix= c("AO_0h_sl6b_sec", "AO_20h_sl4b_sec"),
  msi_slide = c("sl6A", "sl4A"),
  msi_bf    = c(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"06102026_AO_0h_sl6A.nd2"),
                file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"06102026_AO_20h_sl4A.nd2")),
  stringsAsFactors = FALSE
)

# ---- section-level registry (hi-res, sec1..6 per slide) --------------------
if_sections <- function() {
  out <- do.call(rbind, lapply(seq_len(nrow(IF_SLIDES)), function(i) {
    s <- IF_SLIDES[i, ]
    data.frame(
      sid_if   = paste0(s$sid_prefix, 1:6),
      slide    = s$slide, group = s$group, secn = 1:6,
      hr_path  = file.path(IF_DIR, paste0(s$hr_prefix, 1:6, ".nd2")),
      over     = s$over, msi_slide = s$msi_slide, msi_bf = s$msi_bf,
      stringsAsFactors = FALSE)
  }))
  out
}

# MSI dataset id(s) for a given slide + section number (e.g. sl6A, 1 -> sec1a/1b)
msi_sids_for_section <- function(msi_slide, secn) {
  grp <- if (grepl("6A", msi_slide)) "0h" else "20h"
  base <- sprintf("AO_%s_%s_sec%d", grp, msi_slide, secn)
  cand <- paste0(base, c("a", "b"))
  cand[file.exists(file.path(REG_CACHE, sprintf("nd2final_%s.rds", cand)))]
}

# ---- pilot ----------------------------------------------------------------
PILOT_SIDS <- c("AO_0h_sl6b_sec1", "AO_20h_sl4b_sec5")
