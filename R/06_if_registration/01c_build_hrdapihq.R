# 01c_build_hrdapihq.R - build the high-quality DAPI raster cache (STEP 0 helper).
#   Block-mean the DAPI channel of every hi-res section .nd2 at F=HQF (finer than
#   the hr4raw F=F_HR) to a compact RAW-count matrix, for crisp DAPI displays.
# In : hi-res section .nd2 files (IF_DIR), IF_SLIDES table (if_config.R)
# Out: cache/register_if/hrdapihq_<sid>.rds  (all sections, both slides)
# Usage: Rscript R/06_if_registration/01c_build_hrdapihq.R
if (Sys.getenv("JAVA_HOME") == "") Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jre1.8.0_491")
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/if_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register_if.R"))
suppressPackageStartupMessages(library(RBioFormats))
HQF <- 2L
for (i in 1:nrow(IF_SLIDES)) {
  S <- IF_SLIDES[i,]
  for (n in 1:6) {
    sid <- paste0(S$sid_prefix, n); out <- file.path(IF_CACHE, sprintf("hrdapihq_%s.rds", sid))
    if (file.exists(out)) { cat("cached", sid, "\n"); next }
    hp <- file.path(IF_DIR, paste0(S$hr_prefix, n, ".nd2"))
    cat("building HQ DAPI", sid, "...\n")
    saveRDS(list(m=nd2_channel_block_mean(hp, ch=DAPI_CH_DEFAULT, F=HQF, normalize=FALSE)$m, F=HQF), out)
  }
}
cat("DONE\n")
