if (Sys.getenv("JAVA_HOME") == "" && dir.exists("C:/Program Files/Java/jre1.8.0_491")) {
  Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jre1.8.0_491")
}

# ---- Parallelism (Windows: SnowParam, SOCK cluster) -------------------------
# Each worker is a separate R session; pick conservatively for memory.
# Dropped 8 -> 4 for the 20-section run (~90,760 px, 5.8x the 5-section pixels):
# 8 workers x 20 sections spiked RAM too high. Overnight runtime is acceptable.
N_WORKERS <- 4L
register_parallel <- function(n = N_WORKERS) {
  if (!requireNamespace("BiocParallel", quietly = TRUE)) return(invisible(NULL))
  bp <- BiocParallel::SnowParam(workers = n, type = "SOCK", progressbar = TRUE)
  BiocParallel::register(bp, default = TRUE)
  invisible(bp)
}

# ============================================================================
# Analysis_R_Final -- publication fork.
# Heavy intermediate caches are REUSED, read-only, from the original Analysis_R
# (single-condition study, so they carry no group split; produced by manual annotation gates).
# This fork only RE-RUNS the analysis/figure layer and writes its own outputs.
#   - Reads of upstream caches go through `cache_in()` (Final-local first, else
#     the read-only upstream cache).
#   - Writes always target CACHE_DIR (Final-local) so the original is untouched.
# ============================================================================
# ---- SINGLE PATH SWITCH -----------------------------------------------------
# Every absolute path in the pipeline is anchored to ONE data root. To run on a
# different machine, set the MSI_ROOT environment variable to your data root
# (the parent folder holding Analysis_R_Final/, Analysis_R/cache/, MSI/,
# LM_Bsections/, and the *.nd2 files) — e.g. in ~/.Renviron:  MSI_ROOT=/data/JoyMetabolGrad
# Leave it unset to use the reference-machine default. Every script's bootstrap
# uses Sys.getenv("MSI_ROOT", "<default>"), so this one variable relocates the
# whole pipeline (incl. the inventory.csv paths, rewritten in load_inventory()).
MSI_DEFAULT_ROOT <- "D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"
DATA_ROOT    <- Sys.getenv("MSI_ROOT", MSI_DEFAULT_ROOT)
PROJECT_ROOT <- file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final")
CACHE_SRC    <- file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R/cache")  # read-only upstream
INVENTORY    <- file.path(PROJECT_ROOT, "inventory.csv")
CACHE_DIR    <- file.path(PROJECT_ROOT, "cache")
RES_DIR      <- file.path(PROJECT_ROOT, "results")
FIG_DIR      <- file.path(PROJECT_ROOT, "figures")

# Resolve a cache file: prefer a Final-local copy (freshly recomputed), else
# fall back to the read-only upstream cache. Use for ALL upstream-cache READS.
cache_in <- function(name) {
  local <- file.path(CACHE_DIR, name)
  if (file.exists(local)) local else file.path(CACHE_SRC, name)
}

# Display label for a sample_id: strip the leading "AO_" and the incubation
# token so nothing surfaces 0h/20h in any figure/title (single-condition study).
disp_id <- function(sid) sub("^AO_", "", sub("_(0h|20h)_", "_", sid))

# ---- LOCKED apical/basolateral class convention (single source of truth) -----
# Used everywhere (outlines, dot/box plots, near-field overlays, legends,
# whole-region figure) so colors are consistent. Polarity classes:
#   apical-out      = apical membrane faces the gel
#   basolateral-out = basolateral faces the gel (apical faces the lumen)
APICAL_CLASSES <- c("apical_out", "basolateral_out", "mixed")
APICAL_COLS <- c(
  apical_out      = "#C2399A",  # magenta
  basolateral_out = "#2C9E4B",  # green
  mixed           = "#888888"   # grey
)
APICAL_LABS <- c(apical_out = "apical-out", basolateral_out = "basolateral-out", mixed = "mixed")

MSI_PIXEL_UM <- 10

ALIGN_ANCHORS <- c(
  taurine_M_H         = 124.0068,
  palmitate_FA16_0_MH = 255.2330,
  AMP_M_H             = 346.0558,
  ADP_M_H             = 426.0222,
  taurocholate_M_H    = 514.2844
)
ALIGN_TOL_PPM        <- 25   # widened to bin width + slack for post-binning search
LOCKMASS_PASS_PPM    <- 2

BUF_LEVELS_UM <- c(10, 20, 50, 80, 160)

PPM_TOL         <- 25   # user's empirical bin-width optimum for this instrument
SNR_THRESH      <- 3
REF_SAMPLE_SIZE <- 1000

IMG_CLIP_HI <- 0.995
IMG_GAMMA   <- 1.0

# ---- Citrate (locked from the authentic standard, 02_CitrateStandard) -------
# True on-instrument citrate [M-H]- centroid measured from the Sodium Citrate
# Tribasic Dihydrate standard (+0.16 ppm vs theoretical 191.019726). All citrate
# extraction (gate + downstream) integrates a +-CITRATE_WIN_PPM window here, read
# from the RAW centroid imzML -- NOT the 25-ppm-binned grid feature (which merges
# citrate with a +16 ppm co-isobar). 7 ppm chosen from a 5/7/10 sweep: max citrate
# capture at 0% shoulder contamination. See R/00_lib/lib_citrate.R.
CITRATE_ANCHOR_MZ <- 191.01976
CITRATE_WIN_PPM   <- 7

load_inventory <- function() {
  inv <- read.csv(INVENTORY, stringsAsFactors = FALSE)
  # relocate the absolute paths in inventory.csv via the same MSI_ROOT switch
  if (DATA_ROOT != MSI_DEFAULT_ROOT)
    for (col in intersect(c("imzml_path","bf_jpg_path","nd2_path"), names(inv)))
      inv[[col]] <- sub(MSI_DEFAULT_ROOT, DATA_ROOT, inv[[col]], fixed = TRUE)
  inv
}

inventory_row <- function(sample_id) {
  inv <- load_inventory()
  row <- inv[inv$sample_id == sample_id, , drop = FALSE]
  if (nrow(row) != 1) {
    stop(sprintf("sample_id '%s' not found in inventory (%s)", sample_id, INVENTORY))
  }
  row
}