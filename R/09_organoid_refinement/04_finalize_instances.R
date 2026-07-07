#!/usr/bin/env Rscript
# ============================================================================
# 04_finalize_instances.R - Phase 09 step 4: normalize so EACH connected ROI has its
#   own unique id. Runs last (after splits R/09_organoid_refinement/
#   01_organoid_split_apply.R, merges/deletes R/09_organoid_refinement/
#   03_island_cleanup_apply.R).
# ----------------------------------------------------------------------------
# Merges of non-adjacent pieces (R/09_organoid_refinement/03_island_cleanup_apply.R)
# and split-barrier remnants (R/09_organoid_refinement/01_organoid_split_apply.R) can
# leave one instance id with >1 connected component (two ROIs sharing a
# number/colour but not touching). Per the invariant "one connected ROI = one
# id", this step relabels every section:
#   - each instance's LARGEST connected component keeps the id;
#   - any OTHER component >= DISCONNECT_MIN_PX gets a NEW unique id (appended);
#   - components < DISCONNECT_MIN_PX are dropped as noise (split-barrier specks).
#
# Operates on the effective segmentation (instances_clean > instances_split >
#   instances) and writes a SEPARATE final layer cache/instances_final_<sid>.rds
#   for any section that changed. Renders prefer final > clean > split > orig.
#   Idempotent: re-derived each run; stale finals cleared first.
#
# Usage: Rscript R/09_organoid_refinement/04_finalize_instances.R [min_px]   (default DISCONNECT_MIN_PX)
# Out  : cache/instances_final_<sid>.rds,
#        results/annotation/organoid_finalize_report.csv
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
suppressPackageStartupMessages(library(EBImage))

ANNOT_RES <- file.path(RES_DIR, "annotation")
SIDECAR_RDS <- file.path(CACHE_DIR, "organoid_apical_label_positions.rds")
REMOVE_CSV  <- file.path(ANNOT_RES, "organoid_remove.csv")   # sid,id[,note]: ROIs to drop
DISCONNECT_MIN_PX <- 10L
log_msg <- function(...) message(sprintf("[44] %s", sprintf(...)))

a <- commandArgs(trailingOnly = TRUE)
if (length(a) >= 1) DISCONNECT_MIN_PX <- as.integer(a[1])

surface_of <- function(T) {
  W <- nrow(T); H <- ncol(T)
  pad <- matrix(FALSE, W + 2L, H + 2L); pad[2:(W + 1L), 2:(H + 1L)] <- T
  n4 <- pad[1:W, 2:(H + 1L)] + pad[3:(W + 2L), 2:(H + 1L)] +
        pad[2:(W + 1L), 1:H]  + pad[2:(W + 1L), 3:(H + 2L)]
  T & (n4 < 4)
}
effective_file <- function(sid) {
  for (lay in c("clean", "split")) {
    f <- file.path(CACHE_DIR, sprintf("instances_%s_%s.rds", lay, sid)); if (file.exists(f)) return(f)
  }
  file.path(CACHE_DIR, sprintf("instances_%s.rds", sid))
}

stopifnot(file.exists(SIDECAR_RDS))
sids <- attr(readRDS(SIDECAR_RDS), "render_sids")

# clear stale final layer so removed disconnections don't leave orphans
old <- list.files(CACHE_DIR, pattern = "^instances_final_.*\\.rds$", full.names = TRUE)
if (length(old)) { file.remove(old); log_msg("cleared %d stale instances_final_*.rds", length(old)) }

# explicit removals (spurious / off-tissue ROIs the user flagged), by final id
removes <- if (file.exists(REMOVE_CSV)) read.csv(REMOVE_CSV, stringsAsFactors = FALSE) else NULL

promos <- list(); drops <- list(); removed_rows <- list(); n_written <- 0L
for (sid in sids) {
  f <- effective_file(sid); inst <- readRDS(f); pos <- inst$instance > 0
  if (!any(pos)) next
  W <- max(inst$x); H <- max(inst$y)
  lab <- matrix(0L, W, H); lab[cbind(inst$x[pos], inst$y[pos])] <- as.integer(inst$instance[pos])
  ids <- sort(unique(as.integer(lab[lab > 0])))
  labN <- lab; next_id <- max(ids) + 1L; changed <- FALSE
  for (k in ids) {
    mask <- lab == k
    cc <- EBImage::bwlabel(mask * 1); n <- max(cc)
    if (n <= 1) next
    sizes <- tabulate(cc[cc > 0]); keep <- which.max(sizes)
    for (c in setdiff(seq_len(n), keep)) {
      sel <- mask & (cc == c)
      if (sizes[c] >= DISCONNECT_MIN_PX) {
        labN[sel] <- next_id
        promos[[length(promos) + 1L]] <- data.frame(sid = sid, old_id = k, new_id = next_id, px = sizes[c], stringsAsFactors = FALSE)
        next_id <- next_id + 1L
      } else {
        labN[sel] <- 0L
        drops[[length(drops) + 1L]] <- data.frame(sid = sid, id = k, px = sizes[c], stringsAsFactors = FALSE)
      }
      changed <- TRUE
    }
  }
  # explicit removals for this section (by id on the finalized canvas)
  rm_ids <- if (!is.null(removes)) as.integer(removes$id[removes$sid == sid]) else integer(0)
  rm_ids <- rm_ids[rm_ids %in% as.integer(labN[labN > 0])]
  if (length(rm_ids)) {
    labN[labN %in% rm_ids] <- 0L; changed <- TRUE
    removed_rows[[length(removed_rows) + 1L]] <- data.frame(sid = sid, removed_id = rm_ids, stringsAsFactors = FALSE)
  }
  if (changed) {
    T <- labN > 0; surf <- surface_of(T)
    out <- inst
    out$instance <- as.integer(labN[cbind(out$x, out$y)]); out$instance[is.na(out$instance)] <- 0L
    out$is_surface <- as.logical(surf[cbind(out$x, out$y)]); out$is_surface[is.na(out$is_surface)] <- FALSE
    saveRDS(out, file.path(CACHE_DIR, sprintf("instances_final_%s.rds", sid)))
    n_written <- n_written + 1L
    log_msg("%s: %d -> %d id(s)  (from %s)", sid, length(ids), length(unique(out$instance[out$instance > 0])), basename(f))
  }
}

promo_df <- if (length(promos)) do.call(rbind, promos) else data.frame(sid=character(), old_id=integer(), new_id=integer(), px=integer())
drop_df  <- if (length(drops))  do.call(rbind, drops)  else data.frame(sid=character(), id=integer(), px=integer())
write.csv(promo_df, file.path(ANNOT_RES, "organoid_finalize_report.csv"), row.names = FALSE)
cat(sprintf("\n[44] floor=%d px. Promoted %d disconnected component(s) to new ids:\n", DISCONNECT_MIN_PX, nrow(promo_df)))
print(promo_df, row.names = FALSE)
cat(sprintf("[44] Dropped %d sub-floor noise speck(s)", nrow(drop_df)))
if (nrow(drop_df)) cat(sprintf(" (px: %s)", paste(drop_df$px, collapse = ","))); cat("\n")
if (length(removed_rows)) { rmdf <- do.call(rbind, removed_rows)
  cat(sprintf("[44] Removed %d flagged ROI(s): %s\n", nrow(rmdf),
              paste(sprintf("%s#%d", sub("AO_","",rmdf$sid), rmdf$removed_id), collapse=", "))) }
log_msg("wrote %d instances_final_*.rds; re-run R/42 (and R/102) -- they prefer final > clean > split > orig.", n_written)
