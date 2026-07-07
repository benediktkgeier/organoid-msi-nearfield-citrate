#!/usr/bin/env Rscript
# ============================================================================
# 03_island_cleanup_apply.R - Phase 09 step 3: apply island DELETE / MERGE actions to
#   the split segmentation, producing the cleaned per-section instances.
# ----------------------------------------------------------------------------
# Reads an auditable action list (results/annotation/organoid_island_actions.csv)
# that Claude writes from the user's comments on organoid_island_cleanup.pdf
# (R/09_organoid_refinement/02_island_cleanup_canvas.R). Each row: sid, op (merge|delete), ids (";"-sep instance ids on the
# split canvas), note. Free-form Adobe notes are interpreted (with the user) into
# this explicit CSV rather than NLP-parsed, so every applied action is on record.
#
# ID rule: PRESERVE ids (so the action list stays valid across split rounds and
#   future cleanup rounds) -- a MERGE sets every member to the group's LOWEST id;
#   a DELETE removes the instance (its id simply disappears -> a gap). No
#   compaction. Targets are resolved on the ORIGINAL labels so merge/delete ids
#   refer to the same organoids regardless of action order.
#
# Operates on the SPLIT segmentation (cache/instances_split_<sid>.rds, else
#   instances_<sid>.rds). Output -> SEPARATE file cache/instances_clean_<sid>.rds
#   (split + originals untouched). R/09_organoid_refinement/02_island_cleanup_canvas.R
#   and R/10_apical_annotation/01_apical_annotate.R prefer clean > split > orig.
#   Idempotent: re-derived from instances_split + the action CSV every run.
#
# Usage: Rscript R/09_organoid_refinement/03_island_cleanup_apply.R [actions.csv]
# Out  : cache/instances_clean_<sid>.rds (sids in the action list),
#        figures/annotation/organoid_island_cleanup_curation.pdf,
#        results/annotation/organoid_island_cleanup_report.csv
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_register.R"))
suppressPackageStartupMessages({ library(png); library(EBImage) })

REG_CACHE <- file.path(CACHE_DIR, "register")
CROP_DIR  <- file.path(FIG_DIR, "registration", "crops")
ANNOT_FIG <- file.path(FIG_DIR, "annotation")
ANNOT_RES <- file.path(RES_DIR, "annotation")
SCALE_UM  <- 100; PX_UM <- MSI_PIXEL_UM
PAGE_W_IN <- 14; PAGE_H_IN <- 10
ORG_PAL <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628","#f781bf",
             "#1b9e77","#d95f02","#7570b3","#66a61e","#e7298a")
org_colors <- function(ids) setNames(ORG_PAL[((seq_along(ids)-1) %% length(ORG_PAL))+1], as.character(ids))
group_of <- function(sid) "incubated_5min"
log_msg <- function(...) message(sprintf("[43] %s", sprintf(...)))

args <- commandArgs(trailingOnly = TRUE)
ACTIONS_CSV <- if (length(args) >= 1) args[1] else file.path(ANNOT_RES, "organoid_island_actions.csv")
stopifnot(file.exists(ACTIONS_CSV))
act <- read.csv(ACTIONS_CSV, stringsAsFactors = FALSE)
act$op <- tolower(trimws(act$op))
parse_ids <- function(s) as.integer(strsplit(gsub("[^0-9;]", "", s), ";")[[1]])

# clear ALL stale cleaned files first, so a dropped action (e.g. an un-deleted or
# un-merged section) doesn't leave an orphan instances_clean overriding the split
old_clean <- list.files(CACHE_DIR, pattern = "^instances_clean_.*\\.rds$", full.names = TRUE)
if (length(old_clean)) { file.remove(old_clean); log_msg("cleared %d stale instances_clean_*.rds", length(old_clean)) }

instance_outlines <- function(M, B, cx0, cy0) {
  M <- EBImage::fillHull(EBImage::closing(M > 0, EBImage::makeBrush(3, "disc")))
  out <- list()
  for (co in EBImage::ocontour(EBImage::bwlabel(M))) {
    if (nrow(co) < 6) next
    p <- apply_affine(B, cbind(co[, 1] + 0.5, co[, 2] + 0.5))
    out[[length(out) + 1L]] <- list(x = c(p[, 1] - cx0 + 1, p[1, 1] - cx0 + 1),
                                    y = c(p[, 2] - cy0 + 1, p[1, 2] - cy0 + 1))
  }
  out
}
label_halo <- function(x, y, txt, col, cex = 1.4) {
  off <- 0.6 * (cex / 1.5) * strwidth("0", cex = cex) / 4 + 0.4
  for (dx in c(-1, 1)) for (dy in c(-1, 1))
    text(x + dx * off, y + dy * off, txt, col = "white", cex = cex, font = 2, xpd = NA)
  text(x, y, txt, col = col, cex = cex, font = 2, xpd = NA)
}
surface_of <- function(T) {
  W <- nrow(T); H <- ncol(T)
  pad <- matrix(FALSE, W + 2L, H + 2L); pad[2:(W + 1L), 2:(H + 1L)] <- T
  n4 <- pad[1:W, 2:(H + 1L)] + pad[3:(W + 2L), 2:(H + 1L)] +
        pad[2:(W + 1L), 1:H]  + pad[2:(W + 1L), 3:(H + 2L)]
  T & (n4 < 4)
}

CURATION_PDF <- file.path(ANNOT_FIG, "organoid_island_cleanup_curation.pdf")
pdf(CURATION_PDF, width = PAGE_W_IN, height = PAGE_H_IN)
report <- list()
sids <- unique(act$sid)

for (sid in sids) {
  xf  <- file.path(REG_CACHE, sprintf("nd2final_%s.rds", sid))
  pf  <- file.path(CROP_DIR, sprintf("optical_%s.png", sid))
  split_f <- file.path(CACHE_DIR, sprintf("instances_split_%s.rds", sid))
  base_f  <- if (file.exists(split_f)) split_f else file.path(CACHE_DIR, sprintf("instances_%s.rds", sid))
  if (!file.exists(xf) || !file.exists(pf) || !file.exists(base_f)) { log_msg("%s: SKIP (missing inputs)", sid); next }
  Xr <- readRDS(xf); B <- Xr$B_msi_nd2; smn <- Xr$scale_msi_nd2; cx0 <- Xr$crop[1]; cy0 <- Xr$crop[2]
  om <- png::readPNG(pf); if (length(dim(om)) == 3) om <- om[, , 1]; cw <- ncol(om); ch <- nrow(om)

  inst <- readRDS(base_f); pos <- inst$instance > 0
  W <- max(inst$x); H <- max(inst$y)
  lab0 <- matrix(0L, W, H); lab0[cbind(inst$x[pos], inst$y[pos])] <- as.integer(inst$instance[pos])
  ids0 <- sort(unique(as.integer(lab0[lab0 > 0])))

  rows <- act[act$sid == sid, , drop = FALSE]
  del_ids <- sort(unique(unlist(lapply(rows$ids[rows$op == "delete"], parse_ids))))
  merge_groups <- lapply(which(rows$op == "merge"), function(i) parse_ids(rows$ids[i]))
  del_ids <- del_ids[del_ids %in% ids0]
  merge_groups <- lapply(merge_groups, function(g) g[g %in% ids0])
  merge_groups <- merge_groups[vapply(merge_groups, length, 1L) >= 2]

  # apply on ORIGINAL-id reference: delete first, then merge (lowest id wins)
  labC <- lab0
  if (length(del_ids)) labC[lab0 %in% del_ids] <- 0L
  merge_map <- character(0)
  for (g in merge_groups) { tgt <- min(g); labC[lab0 %in% g] <- tgt
    merge_map <- c(merge_map, sprintf("{%s}->%d", paste(sort(g), collapse=","), tgt)) }
  ids1 <- sort(unique(as.integer(labC[labC > 0])))

  T <- labC > 0; surf <- surface_of(T)
  out <- inst
  out$instance   <- as.integer(labC[cbind(out$x, out$y)]); out$instance[is.na(out$instance)] <- 0L
  out$is_surface <- as.logical(surf[cbind(out$x, out$y)]); out$is_surface[is.na(out$is_surface)] <- FALSE
  saveRDS(out, file.path(CACHE_DIR, sprintf("instances_clean_%s.rds", sid)))

  report[[length(report) + 1L]] <- data.frame(
    sid = sid, group = group_of(sid), n_before = length(ids0), n_after = length(ids1),
    merged = paste(merge_map, collapse = " "), deleted = paste(del_ids, collapse = ";"),
    stringsAsFactors = FALSE)
  log_msg("%s: %d -> %d organoid(s) | merged: %s | deleted: %s", sid, length(ids0), length(ids1),
          if (length(merge_map)) paste(merge_map, collapse=" ") else "-",
          if (length(del_ids)) paste(del_ids, collapse=";") else "-")

  # curation page: BF + cleaned outlines (coloured + id) + deleted (red hatch) ----
  par(mar = c(1, 1, 3, 1)); plot.new(); plot.window(c(1, cw), c(ch, 1), asp = 1)
  rasterImage(om / max(om), 1, ch, cw, 1, interpolate = TRUE)
  oc <- org_colors(ids1)
  for (k in ids1) {
    kc <- oc[as.character(k)]
    for (poly in instance_outlines(labC == k, B, cx0, cy0)) polygon(poly$x, poly$y, border = kc, lwd = 1.8)
    cen <- colMeans(which(labC == k, arr.ind = TRUE)); pc <- apply_affine(B, matrix(cen, nrow = 1))
    label_halo(pc[1] - cx0 + 1, pc[2] - cy0 + 1, as.character(k), kc, cex = 1.4)
  }
  if (length(del_ids)) for (d in del_ids)        # deleted islands: red dashed outline
    for (poly in instance_outlines(lab0 == d, B, cx0, cy0)) {
      polygon(poly$x, poly$y, border = "red", lwd = 2.0, lty = 2)
      label_halo(mean(poly$x), mean(poly$y), "DEL", "red", cex = 1.1)
    }
  title(sprintf("%s  (%s)   cleanup  %d -> %d   merged: %s   deleted: %s",
                sub("AO_", "", sid), group_of(sid), length(ids0), length(ids1),
                if (length(merge_map)) paste(merge_map, collapse=" ") else "-",
                if (length(del_ids)) paste(del_ids, collapse=";") else "-"),
        cex.main = 1.05, col.main = "#6a3d9a", font.main = 2)
  umpx <- PX_UM / smn; bp <- SCALE_UM / umpx
  segments(cw - bp - 10, ch - 14, cw - 10, ch - 14, lwd = 3, col = "black")
  text(cw - bp / 2 - 10, ch - 28, sprintf("%d um", SCALE_UM), cex = 0.7)
}
dev.off()

rep_df <- if (length(report)) do.call(rbind, report) else
  data.frame(sid=character(), group=character(), n_before=integer(), n_after=integer(),
             merged=character(), deleted=character())
write.csv(rep_df, file.path(ANNOT_RES, "organoid_island_cleanup_report.csv"), row.names = FALSE)
cat("\n[43] ===== CLEANUP REPORT =====\n"); print(rep_df, row.names = FALSE)
log_msg("wrote %d instances_clean_*.rds; curation -> %s", length(sids), CURATION_PDF)
log_msg("NEXT: re-run R/42 (and R/102) -- they prefer instances_clean when present.")
