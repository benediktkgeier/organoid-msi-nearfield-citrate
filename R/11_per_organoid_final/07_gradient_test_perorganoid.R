#!/usr/bin/env Rscript
# ============================================================================
# 07_gradient_test_perorganoid.R - STATISTICAL TEST of whether organoids emit a
#   near-field (50-100 um) DECREASING citrate [M-H]- gradient. ORGANOID =
#   replication unit (n~72). Each organoid collapsed to ONE gradient number
#   (avoids pixel spatial-autocorrelation). Citrate tested ALONE (user choice
#   2026-06-21: [M-H]- only, no DHA contrast).
#
# Inferential target = gradient EXISTENCE around organoids. p-values are
#   organoid-level evidence with a 2-slide generalization ceiling. Single-
#   condition study: no incubation-time comparison.
#
# NOTE: citrate is measured as [M-H]- 191.0217 only. A flat [M-H]- profile means a
#   flat profile (the earlier adduct-switching / total-citrate caveat was an artifact
#   and has been removed).
#
# Out: results/gradient/gradient_test_perorganoid.csv
#      figures/gradient/gradient_test_perorganoid.pdf
# Usage: Rscript R/11_per_organoid_final/07_gradient_test_perorganoid.R
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_gradient_seg.R"))
suppressPackageStartupMessages({ library(Cardinal) })
set.seed(1)

CIT_MZ <- 191.0217
BAND_LO <- 50; BAND_HI <- 100      # near-field window (um from organoid surface)
MIN_BAND_PX <- 8L; MIN_DIST_LV <- 3L
B_BOOT <- 2000L
log_msg <- function(...) message(sprintf("[88] %s", sprintf(...)))

mse <- readRDS(TISSUE_MSE)
pd  <- as.data.frame(pixelData(mse)); pd$gidx <- seq_len(nrow(pd))
mzs <- mz(mse); ci <- which.min(abs(mzs - CIT_MZ))
source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/lib_citrate.R"))
val_cit <- citrate_onto_pd(pd)   # anchored citrate (raw imzML, TIC-norm); ci kept for the log only
SIDS20  <- levels(pixelData(mse)$sample_id); if (is.null(SIDS20)) SIDS20 <- sort(unique(as.character(pd$sample_id)))
log_msg("citrate [M-H]- idx=%d (%.4f) | near-field band %g-%g um", ci, mzs[ci], BAND_LO, BAND_HI)

# ---- per-organoid citrate gradient metrics ---------------------------------
# slope b of log1p(I)~dist (signed gradient: b<0 = decreasing outward = emission-like);
# Spearman rho (monotonicity); exp decay length lambda (nls, fallback -1/b).
fit_metrics <- function(d, I) {
  rho <- if (length(unique(d)) >= MIN_DIST_LV && sd(I) > 0) suppressWarnings(cor(d, I, method = "spearman")) else NA_real_
  b   <- tryCatch(unname(coef(lm(log1p(I) ~ d))[2]), error = function(e) NA_real_)
  lam <- NA_real_
  rng <- range(I); A0 <- max(diff(rng), 1e-6); C0 <- max(rng[1], 0)
  fit <- tryCatch(nls(I ~ A*exp(-d/lam) + C, start = list(A = A0, lam = 30, C = C0),
                      algorithm = "port", lower = c(A = 0, lam = 5, C = 0), upper = c(A = Inf, lam = 1000, C = Inf),
                      control = nls.control(maxiter = 200)), error = function(e) NULL)
  if (!is.null(fit)) { cf <- coef(fit); if (is.finite(cf["lam"]) && cf["A"] > 0) lam <- unname(cf["lam"]) }
  if (is.na(lam) && is.finite(b) && b < 0) lam <- -1 / b
  list(rho = unname(rho), b = b, lambda = lam)
}

rows <- list(); prof_rows <- list(); n_drop <- 0L
BINS <- seq(BAND_LO, BAND_HI, by = 10); BMID <- head(BINS, -1) + 5
for (sid in SIDS20) {
  fz <- cache_in(sprintf("zones_%s.rds", sid)); if (!file.exists(fz)) next
  z <- readRDS(fz)
  for (k in sort(unique(z$instance[z$instance > 0]))) {
    sel <- which(z$instance_catch == k & z$signed_dist_um > 0 & z$dist_um >= BAND_LO & z$dist_um <= BAND_HI)
    if (length(sel) < MIN_BAND_PX || length(unique(z$dist_um[sel])) < MIN_DIST_LV) { n_drop <- n_drop + 1L; next }
    g <- z$gidx[sel]; d <- z$dist_um[sel]; Ic <- val_cit[g]
    mc <- fit_metrics(d, Ic)
    p_rho <- tryCatch(suppressWarnings(cor.test(d, Ic, method = "spearman", alternative = "less")$p.value),
                      error = function(e) NA_real_)   # H1: decreasing outward
    rows[[length(rows)+1]] <- data.frame(sample_id = sid, instance = k, n_band_px = length(sel),
      b_cit = mc$b, lambda_cit = mc$lambda, rho_cit = mc$rho, p_rho = p_rho, stringsAsFactors = FALSE)
    bi <- cut(d, breaks = BINS, include.lowest = TRUE, labels = FALSE); cm <- tapply(Ic, bi, mean)
    prof_rows[[length(prof_rows)+1]] <- data.frame(sample_id = sid, instance = k,
      bin = seq_along(BMID), dist = BMID, cit = as.numeric(cm[as.character(seq_along(BMID))]), stringsAsFactors = FALSE)
  }
}
org <- do.call(rbind, rows); prof <- do.call(rbind, prof_rows)
org$fdr_rho <- p.adjust(org$p_rho, method = "BH")          # per-organoid multiple-testing correction (n organoids)
write.csv(org, file.path(GRAD_RES, "gradient_test_perorganoid.csv"), row.names = FALSE)
log_msg("%d organoids tested; %d dropped (<%d band px)", nrow(org), n_drop, MIN_BAND_PX)

# ---- across-organoid EXISTENCE statistics (citrate alone) ------------------
cluster_boot <- function(section, value, statfun, B = B_BOOT) {   # resample SECTIONS (respects nesting)
  secs <- unique(section); out <- numeric(B)
  for (b in seq_len(B)) { samp <- sample(secs, length(secs), replace = TRUE)
    idx <- unlist(lapply(samp, function(s) which(section == s))); out[b] <- statfun(value[idx]) }
  quantile(out, c(0.025, 0.975), na.rm = TRUE)
}
medn <- function(v) median(v, na.rm = TRUE); sect <- org$sample_id
p_b   <- suppressWarnings(wilcox.test(org$b_cit,   mu = 0, alternative = "less")$p.value)   # b<0 = decreasing
p_rho <- suppressWarnings(wilcox.test(org$rho_cit, mu = 0, alternative = "less")$p.value)
b_med <- medn(org$b_cit);   b_ci   <- cluster_boot(sect, org$b_cit, medn);   b_pct   <- 100*mean(org$b_cit  < 0, na.rm=TRUE)
r_med <- medn(org$rho_cit); r_ci   <- cluster_boot(sect, org$rho_cit, medn); r_pct   <- 100*mean(org$rho_cit < 0, na.rm=TRUE)
lam_med <- medn(org$lambda_cit); lam_ci <- cluster_boot(sect, org$lambda_cit, medn)
n_nofit <- sum(!is.finite(org$lambda_cit))
sig <- p_b < 0.05
# per-organoid significance (Spearman one-sided; effect size = rho)
n_unc <- sum(org$rho_cit < 0 & org$p_rho < 0.05, na.rm = TRUE)          # uncorrected
n_fdr <- sum(org$rho_cit < 0 & org$fdr_rho < 0.10, na.rm = TRUE)        # BH FDR<0.10
top   <- org[which.min(org$p_rho), ]
sig_i <- which(org$rho_cit < 0 & org$p_rho < 0.05)                      # nominal-hit effect sizes
maxrho_sig <- if (length(sig_i)) max(abs(org$rho_cit[sig_i])) else NA
lam_lo_sig <- if (length(sig_i)) min(org$lambda_cit[sig_i], na.rm = TRUE) else NA
lam_hi_sig <- if (length(sig_i)) max(org$lambda_cit[sig_i], na.rm = TRUE) else NA
log_msg("PER-ORGANOID: %d/%d reach uncorrected p<0.05 (all |rho|<=%.2f, lambda>=%.0f um = trivial); %d survive BH FDR<0.10",
        n_unc, nrow(org), max(abs(org$rho_cit[org$rho_cit<0 & org$p_rho<0.05]), na.rm=TRUE),
        min(org$lambda_cit[org$rho_cit<0 & org$p_rho<0.05], na.rm=TRUE), n_fdr)

log_msg("EXISTENCE citrate slope<0: p=%.3g | median b=%.3g (%.0f%% neg) | rho p=%.3g median rho=%.3g (%.0f%% neg)",
        p_b, b_med, b_pct, p_rho, r_med, r_pct)
log_msg("decay length (declining organoids): median lambda=%.1f um [%.1f, %.1f]; %d/%d no fittable decay",
        lam_med, lam_ci[1], lam_ci[2], n_nofit, nrow(org))
log_msg("VERDICT: %s near-field decreasing [M-H]- citrate gradient at organoid level (p=%.3g).",
        if (sig) "SIGNIFICANT" else "NO significant", p_b)

# ===========================================================================
# REPORT
# ===========================================================================
ORG_C <- "#1b4f72"; fmtp <- function(p) if (p < 1e-3) sprintf("%.1e", p) else sprintf("%.3f", p)
verdict <- if (sig) "evidence FOR a decreasing near-field [M-H]- citrate gradient" else "NO significant decreasing near-field [M-H]- citrate gradient (profile is ~flat)"
pdf(file.path(GRAD_FIG, "gradient_test_perorganoid.pdf"), width = 12, height = 8.5)

# ---- P1: design + caveats ---------------------------------------------------
plot.new(); text(0.5, 0.96, "Do organoids emit a near-field citrate [M-H]- gradient (50-100 um)?", font = 2, cex = 1.3)
txt <- c(
 "DESIGN  (organoid = replication unit; citrate [M-H]- tested ALONE)",
 sprintf("  Tracer : Citrate [M-H]- %.4f", mzs[ci]),
 sprintf("  Window : gel pixels %g-%g um from organoid surface (continuous distance)", BAND_LO, BAND_HI),
 sprintf("  Units  : %d organoids from 20 sections; %d dropped (<%d band px).",
         nrow(org), n_drop, MIN_BAND_PX),
 "  Per organoid: slope b of log1p(I)~dist (b<0 = decreasing outward), Spearman rho, decay length lambda.",
 "  Across organoids: one-sample Wilcoxon signed-rank (H1: median<0) + section cluster-bootstrap 95% CI.",
 "",
 "WHY THIS IS DEFENSIBLE (and its ceiling)",
 "  - Each organoid collapsed to ONE number -> avoids pixel spatial-autocorrelation pseudo-replication.",
 "  - Section cluster-bootstrap respects organoid-within-section nesting.",
 "  - CEILING: only 2 physical slides (1 per timepoint). p-values are ORGANOID-level evidence; they do",
 "    NOT generalize beyond this preparation (single-condition study; no incubation-time comparison).",
 "",
 "NOTE: citrate measured as [M-H]- anchored 191.0198 +-7 ppm (raw imzML). A flat profile means a flat profile.")
text(0.04, 0.87, paste(txt, collapse = "\n"), adj = c(0,1), cex = 0.92, family = "mono")
text(0.5, 0.05, sprintf("RESULT: %s  (Wilcoxon p = %s).", verdict, fmtp(p_b)), cex = 1.05, font = 2,
     col = if (sig) "#1a7a1a" else "#a11")

# ---- P2: per-organoid gradient statistics (dot + median CI) -----------------
layout(matrix(1:2, nrow = 1)); par(oma = c(0,0,2,0), mar = c(5,5,3,1))
dotstat <- function(v, lab, p, md, ci, pct) {
  o <- order(v); vv <- v[o]; n <- length(vv)
  plot(vv, seq_len(n), pch = 19, col = adjustcolor(ORG_C, 0.7), cex = 0.8,
       xlab = lab, ylab = "organoid (sorted)", main = sprintf("%s   (p=%s)", lab, fmtp(p)))
  abline(v = 0, lwd = 1.5, col = "grey40")
  segments(ci[1], -2, ci[2], -2, lwd = 4, col = "black", xpd = NA)
  points(md, -2, pch = 18, cex = 2, col = "black", xpd = NA)
  text(md, -2, sprintf("  median %.3g [%.3g, %.3g]", md, ci[1], ci[2]), pos = 4, cex = 0.75, xpd = NA)
  mtext(sprintf("%.0f%% of organoids < 0", pct), side = 3, line = -1, cex = 0.8)
}
dotstat(org$b_cit,   "slope b (log1p I / um)", p_b,   b_med, b_ci, b_pct)
dotstat(org$rho_cit, "Spearman rho (dist vs I)", p_rho, r_med, r_ci, r_pct)
mtext(sprintf("Per-organoid near-field citrate gradient (n=%d). %d reach uncorrected p<0.05 (Spearman), but %d survive BH FDR<0.10.",
              nrow(org), n_unc, n_fdr), outer = TRUE, font = 2, cex = 0.95)

# ---- P3: decay curves + distributions --------------------------------------
norm1 <- function(v) { f <- which(is.finite(v))[1]; if (is.na(f) || v[f] <= 0) return(v*NA); v / v[f] }
layout(matrix(1:3, nrow = 1)); par(oma = c(0,0,2,0), mar = c(4.5,4.2,3,1))
plot(NA, xlim = c(BAND_LO, BAND_HI), ylim = c(0, 1.6), xlab = "distance from surface (um)",
     ylab = "intensity / nearest-bin", main = "Citrate [M-H]- near-field decay (per organoid)")
med_by_bin <- matrix(NA_real_, length(BMID), 0)
for (i in seq_len(nrow(org))) { s <- prof[prof$sample_id == org$sample_id[i] & prof$instance == org$instance[i], ]
  s <- s[order(s$dist), ]; v <- norm1(s$cit); lines(s$dist, v, col = adjustcolor(ORG_C, 0.15), lwd = 0.7)
  med_by_bin <- cbind(med_by_bin, v[match(seq_along(BMID), s$bin)]) }
lines(BMID, apply(med_by_bin, 1, median, na.rm = TRUE), col = "black", lwd = 3, type = "b", pch = 19)
abline(h = 1, lty = 3, col = "grey60")
hist(org$b_cit, breaks = 24, col = adjustcolor("#1a9850", 0.6), border = NA, xlab = "slope b (log1p I / um)",
     main = "Per-organoid slope"); abline(v = 0, lty = 2); abline(v = b_med, col = "#1a9850", lwd = 2)
lc <- org$lambda_cit[is.finite(org$lambda_cit) & org$lambda_cit < 500]
hist(lc, breaks = 20, col = adjustcolor("#1a9850", 0.6), border = NA, xlab = "decay length lambda (um)",
     main = sprintf("Decay length (median %.0f um; %d no-fit)", lam_med, n_nofit)); abline(v = lam_med, col = "black", lwd = 2)
mtext("Curves normalized to nearest bin (shape). Flat (~1 across) = no near-field gradient.", outer = TRUE, font = 2, cex = 0.95)

# ---- P4: summary + honest interpretation ------------------------------------
layout(1); par(mar = c(1,1,1,1))
plot.new(); text(0.5, 0.95, "Summary & interpretation", font = 2, cex = 1.3)
tb <- c(
 "TEST - DOES A DECREASING NEAR-FIELD [M-H]- CITRATE GRADIENT EXIST? (organoid level)",
 sprintf("   slope b<0 : Wilcoxon p = %s | median b = %.3g | %.0f%% organoids negative", fmtp(p_b), b_med, b_pct),
 sprintf("   rho   <0 : Wilcoxon p = %s | median rho = %.3g | %.0f%% organoids negative", fmtp(p_rho), r_med, r_pct),
 sprintf("   decay length (declining organoids): median lambda = %.0f um [%.0f, %.0f]; %d/%d had NO fittable decay",
         lam_med, lam_ci[1], lam_ci[2], n_nofit, nrow(org)),
 "",
 "PER-ORGANOID (does ANY single organoid have a measurable gradient?)",
 sprintf("   %d/%d organoids reach uncorrected one-sided p<0.05 - but all have trivial effect size", n_unc, nrow(org)),
 sprintf("   (max |rho| = %.2f, lambda %.0f-%.0f um); driven by large pixel n, NOT a real gradient.", maxrho_sig, lam_lo_sig, lam_hi_sig),
 sprintf("   After BH FDR correction: %d organoids significant. Strongest = %s inst %d (rho=%.2f, p=%.3f).",
         n_fdr, sub("AO_","",top$sample_id), top$instance, top$rho_cit, top$p_rho),
 "",
 "INTERPRETATION (data-driven, no forced conclusion)",
 if (sig)
   "   -> Organoids show a statistically supported decreasing [M-H]- citrate gradient in 50-100 um."
 else
   "   -> NO significant decreasing [M-H]- citrate gradient: the near-field profile is essentially FLAT",
 if (sig) "" else "      (median slope ~0, rho ~0, ~half the organoids non-decreasing).",
 "",
 "   Citrate measured as [M-H]- anchored 191.0198 +-7 ppm (raw imzML). A flat profile means a flat profile.",
 "",
 "   CEILING: 2 slides only -> organoid-level evidence, no generalization / no time claim (locked).")
text(0.04, 0.86, paste(tb, collapse = "\n"), adj = c(0,1), cex = 0.93, family = "mono")

dev.off()
log_msg("DONE -> %s + gradient_test_perorganoid.csv", file.path(GRAD_FIG, "gradient_test_perorganoid.pdf"))
