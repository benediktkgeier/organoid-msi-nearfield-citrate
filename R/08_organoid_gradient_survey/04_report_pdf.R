#!/usr/bin/env Rscript
# ============================================================================
# 04_report_pdf.R - pooled gradient survey report (single condition).
# ----------------------------------------------------------------------------
# Supporting-context survey across ALL sections (no group split, no Delta-rho).
# Pages:
#   1. rho_out ranking barplot: most-outward and most-inward ions (known vs
#      unknown), from the pooled per-ion Spearman gradient.
#   2+. Per top-outward ion: pooled outward + inward zone-profile line plots.
#
# Input : results/gradient/gradient_stats.csv, zone_profiles_long.csv
# Output: figures/gradient/gradient_report.pdf
# Usage : Rscript R/08_organoid_gradient_survey/04_report_pdf.R
# ============================================================================

source(file.path(Sys.getenv("MSI_ROOT","D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad"),"Analysis_R_Final/R/00_lib/gradient_config.R"))

N_TOP <- 16L   # per-ion profile pages (top outward by rho_out)

stats <- read.csv(file.path(GRAD_RES, "gradient_stats.csv"), stringsAsFactors = FALSE)
long  <- read.csv(file.path(GRAD_RES, "zone_profiles_long.csv"), stringsAsFactors = FALSE)

# ---- pooled zone-profile line plot (single condition) ---------------------
profile_plot <- function(feat, direction, xlab) {
  sub <- long[long$feat == feat & long$direction == direction, ]
  sub <- sub[order(sub$zone), ]
  zlabs <- sub$zone_label; zk <- sub$zone
  par(mar = c(4, 4, 3, 1))
  ymax <- max(sub$mean_intensity, na.rm = TRUE); if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  plot(NA, xlim = range(zk), ylim = c(0, ymax * 1.05), xaxt = "n",
       xlab = xlab, ylab = "mean intensity (pooled)", main = sprintf("%s profile", direction))
  axis(1, at = zk, labels = zlabs, cex.axis = 0.8)
  lines(sub$zone, sub$mean_intensity, col = "#1b4f72", lwd = 2.4, type = "b", pch = 19)
}

pdf(file.path(GRAD_FIG, "gradient_report.pdf"), width = 12, height = 8.5)

# ===== Page 1: rho_out ranking =====
par(mfrow = c(1, 2), mar = c(4, 11, 3, 1))
rk   <- stats[!is.na(stats$rho_out), ]
topN <- head(rk[order(-rk$rho_out), ], 20)
botN <- head(rk[order(rk$rho_out), ], 20)
barplot(rev(topN$rho_out), horiz = TRUE, names.arg = rev(topN$label), las = 1,
        cex.names = 0.6, col = ifelse(rev(topN$is_known), "#d62728", "#7f7f7f"),
        main = "Most outward (rho_out > 0)", xlab = "rho_out (pooled)")
barplot(rev(botN$rho_out), horiz = TRUE, names.arg = rev(botN$label), las = 1,
        cex.names = 0.6, col = ifelse(rev(botN$is_known), "#1f77b4", "#7f7f7f"),
        main = "Most inward (rho_out < 0)", xlab = "rho_out (pooled)")
mtext("Pooled per-ion outward gradient across all sections (single condition; descriptive)",
      outer = TRUE, line = -1.5, font = 2, cex = 1.0)

# ===== Per-ion pages: pooled outward + inward profiles =====
top_ions <- head(stats[order(-stats$rho_out), ], N_TOP)
for (i in seq_len(nrow(top_ions))) {
  r <- top_ions[i, ]; feat <- r$feat
  layout(matrix(c(1, 2), nrow = 1)); par(oma = c(0, 0, 3, 0))
  profile_plot(feat, "outward", "distance from surface into CMC (um)")
  profile_plot(feat, "inward",  "distance from surface into organoid (um)")
  ktag <- if (isTRUE(r$is_known)) "KNOWN" else "unknown"
  mtext(sprintf("%s  |  %s  |  rho_out=%.2f   rho_in=%.2f", r$label, ktag, r$rho_out, r$rho_in),
        outer = TRUE, line = -1.0, font = 2, cex = 1.05)
}

dev.off()
cat(sprintf("[04] DONE -> %s (%d per-ion pages)\n",
            file.path(GRAD_FIG, "gradient_report.pdf"), nrow(top_ions)))