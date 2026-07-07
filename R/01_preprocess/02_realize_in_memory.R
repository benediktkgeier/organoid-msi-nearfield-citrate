#!/usr/bin/env Rscript
# 02_realize_in_memory.R
# Make the 20-section peaks_combined.rds SELF-CONTAINED (zero .ibd dependency).
#
# WHY: a Cardinal v3 MSE saved via saveRDS keeps featureData/pixelData in memory
# but its spectra is a matter `sparse_mat` whose values/index are lazily read
# from cache/imzml/*.ibd. The .rds stores only byte-offset pointers. So any
# spectra read (ion-image render, SSC) opens the .ibd -- and matter opens it
# READ-WRITE, so the file must be present, writable, and byte-identical. A
# parallel re-export already corrupted the frozen 5-section MSE this way
# (memory: project_organoid_ibd_dependency).
#
# FIX (verified on subsets): materialise spectra to a dense matrix once, convert
# to an IN-MEMORY matter sparse_mat (rowMaj=FALSE, path=NULL -> Cardinal-native
# class, no file backing, ~3 GB sparse vs ~8 GB dense), assign back, re-save.
# Result: spectra has no .ibd path; the object is fully portable.
#
# Writes a one-time backup of the lazy object, then overwrites peaks_combined.rds
# atomically (temp file + rename) after an on-disk reload read passes.
# Usage: Rscript R/01_preprocess/02_realize_in_memory.R

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_paths.R"))
suppressPackageStartupMessages({ library(Cardinal); library(matter) })

PEAKS   <- file.path(CACHE_DIR, "peaks_combined.rds")
BACKUP  <- file.path(CACHE_DIR, "peaks_combined_lazy_prerelease.rds")
TMP     <- file.path(CACHE_DIR, "peaks_combined.rds.tmp")
log_msg <- function(...) message(sprintf("[06 %s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

t0 <- Sys.time()
log_msg("Loading %s ...", basename(PEAKS))
mse <- readRDS(PEAKS)
nr <- nrow(mse); nc <- ncol(mse)
cls0 <- class(spectra(mse))[1]
log_msg("MSE: %d features x %d pixels; spectra class = %s", nr, nc, cls0)

if (cls0 == "matter_vec" || cls0 == "sparse_mat") {
  has_path <- !is.null(attr(spectra(mse), "path"))
  log_msg("spectra is matter-backed; file-backed (.ibd) = %s", has_path)
}

# One-time backup of the lazy object (cheap insurance; depends on .ibd, kept only
# until the self-contained version is verified).
if (!file.exists(BACKUP)) {
  file.copy(PEAKS, BACKUP)
  log_msg("Backed up lazy object -> %s", basename(BACKUP))
} else {
  log_msg("Backup already exists: %s (not overwritten)", basename(BACKUP))
}

# ---- Realize: dense (transient) -> in-memory matter sparse_mat --------------
log_msg("Materialising spectra to dense (transient ~%.1f GB)...", nr * nc * 8 / 1024^3)
td <- Sys.time()
dense <- as.matrix(spectra(mse))
log_msg("  dense %.2f GB in %.0fs; nonzero frac %.3f",
        object.size(dense) / 1024^3,
        as.numeric(difftime(Sys.time(), td, units = "secs")),
        mean(dense != 0))

log_msg("Converting to in-memory matter sparse_mat...")
sp_mem <- matter::sparse_mat(dense, rowMaj = FALSE)
rm(dense); gc()
stopifnot(is.null(attr(sp_mem, "path")))         # must be in-memory
spectra(mse) <- sp_mem
log_msg("  spectra now: class=%s, file-backed=%s, ~%.1f GB",
        class(spectra(mse))[1], !is.null(attr(spectra(mse), "path")),
        object.size(spectra(mse)) / 1024^3)

# ---- In-process self-containment proof: mid-matrix read --------------------
probe <- as.matrix(spectra(mse)[seq_len(min(3, nr)),
                                round(nc * 0.6) + 0:4, drop = FALSE])
log_msg("In-memory mid-matrix read OK (%s)", paste(dim(probe), collapse = "x"))

# ---- Save (temp) then verify on-disk reload, then atomic replace -----------
log_msg("Saving temp self-contained MSE...")
saveRDS(mse, TMP)
log_msg("Reloading temp from disk to verify...")
mse2 <- readRDS(TMP)
stopifnot(is.null(attr(spectra(mse2), "path")))  # reloaded object has no .ibd path
r2 <- as.matrix(spectra(mse2)[seq_len(min(3, nrow(mse2))),
                              round(ncol(mse2) * 0.6) + 0:4, drop = FALSE])
stopifnot(isTRUE(all.equal(r2, probe)))          # reloaded values match
log_msg("Reload read OK + values match in-process realize.")

file.rename(TMP, PEAKS)
log_msg("Replaced %s (self-contained).", basename(PEAKS))

# ---- Report ---------------------------------------------------------------
sz <- file.info(PEAKS)$size / 1024^2
bz <- file.info(BACKUP)$size / 1024^2
log_msg("DONE in %.1f min. peaks_combined.rds = %.0f MB (was lazy %.0f MB).",
        as.numeric(difftime(Sys.time(), t0, units = "mins")), sz, bz)
cat("\nspectra class is now in-memory matter sparse_mat: zero cache/imzml/*.ibd dependency.\n")
cat("The .ibd files can now be backed up / locked / moved without affecting this RDS.\n")
