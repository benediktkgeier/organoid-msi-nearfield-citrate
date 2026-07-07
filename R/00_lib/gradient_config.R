#!/usr/bin/env Rscript
# ============================================================================
# gradient_config.R - shared config for the GRADIENT pipeline.
# ----------------------------------------------------------------------------
# Single-condition study: all sections were incubated in CMC for at least 5 min
# before freezing and are treated identically (no incubation-time comparison).
# The gradient survey is a POOLED descriptive across all sections -- no group
# split and no Delta-rho.
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_paths.R"))

# ---- Datasets (all sections, single condition) ----------------------------
GRAD_SIDS <- load_inventory()$sample_id

# ---- Segmentation / ring constants ----------------------------------------
MIN_INSTANCE_PX <- 50L                       # drop connected components < 50 px
# Outward ring zones (um from organoid surface into CMC gel). 7 zones: the
# 250/500 um zones extend coverage to the full gel field of the ion images
# (user 2026-06-16) beyond the original 160 um diffusion window.
OUT_ZONE_UM     <- c(BUF_LEVELS_UM, 250, 500) # 10,20,50,80,160,250,500
OUT_BREAKS_UM   <- c(0, OUT_ZONE_UM)          # 0,10,20,50,80,160,250,500 -> zones 1..7
# Inward ring breaks (um from surface into the organoid). 10 um steps, capped.
IN_BREAKS_UM    <- c(0, 10, 20, 30, 40, Inf)  # zones 1..5 (last = core)
IN_ZONE_LAB     <- c("0-10", "10-20", "20-30", "30-40", "40+")

GRAD_CACHE <- CACHE_DIR
GRAD_RES   <- file.path(RES_DIR, "gradient")
GRAD_FIG   <- file.path(FIG_DIR, "gradient")
dir.create(GRAD_RES, showWarnings = FALSE, recursive = TRUE)
dir.create(GRAD_FIG, showWarnings = FALSE, recursive = TRUE)

TISSUE_MSE <- cache_in("peaks_tissue_combined.rds")   # read-only reuse from upstream cache

# ---- Helpers ---------------------------------------------------------------
# Build a section's (x,y)->matrix index map. Returns list with W,H,xmin,ymin and
# functions to go between global pixel order and matrix [ix,iy].
section_grid <- function(x, y) {
  xmin <- min(x); ymin <- min(y)
  W <- max(x) - xmin + 1L
  H <- max(y) - ymin + 1L
  ix <- x - xmin + 1L
  iy <- y - ymin + 1L
  list(W = W, H = H, xmin = xmin, ymin = ymin, ix = ix, iy = iy)
}
