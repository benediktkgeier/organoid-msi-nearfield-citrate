# 01b_build_hr4raw.R - build the 4-channel hi-res IF raster cache (STEP 0 helper).
#   Block-mean all channels of each hi-res section .nd2 (F=F_HR) to a compact
#   RAW-count array and cache it, so downstream IF steps read the cache not the .nd2.
# In : hi-res section .nd2 files (IF_DIR), IF_SLIDES table (if_config.R)
# Out: cache/register_if/hr4raw_<sid>.rds  (one per section of the slide)
# Usage: Rscript R/06_if_registration/01b_build_hr4raw.R [slide]   (default sl6b)
if (Sys.getenv("JAVA_HOME") == "") Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jre1.8.0_491")
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/if_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register_if.R"))
suppressPackageStartupMessages(library(RBioFormats))
args <- commandArgs(trailingOnly=TRUE)
slb  <- if (length(args)) args[1] else "sl6b"
S <- IF_SLIDES[IF_SLIDES$slide==slb,]
for (n in 1:6) {
  sid <- paste0(S$sid_prefix, n); out <- file.path(IF_CACHE, sprintf("hr4raw_%s.rds", sid))
  if (file.exists(out)) { cat("cached", sid, "\n"); next }
  hp <- file.path(IF_DIR, paste0(S$hr_prefix, n, ".nd2"))
  cat("building", sid, "...\n"); saveRDS(nd2_block_mean_allch(hp, F=F_HR, normalize=FALSE), out)
}
cat("DONE\n")
