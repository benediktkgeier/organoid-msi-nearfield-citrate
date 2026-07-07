#!/usr/bin/env Rscript
# =============================================================================
# generate_methods_report.R
# Detailed methods / pipeline report for the organoid MSI near-field citrate
# study (apical-out vs basolateral-out). Renders a multi-page A4 PDF via base `pdf()`
# device + grid (no pandoc / LaTeX dependency).
#
#   Output: figures/methods_report/methods_pipeline_report.pdf
#
# Content is drawn from the verified on-disk pipeline (branch restructure-pipeline)
# and R/00_lib/lib_paths.R locked constants. Numbers without a verified source
# are intentionally omitted rather than guessed.
# =============================================================================

suppressPackageStartupMessages(library(grid))

OUT_DIR <- file.path("figures", "methods_report")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
OUT_PDF <- file.path(OUT_DIR, "methods_pipeline_report.pdf")

# ---- page geometry (A4 portrait, inches) -----------------------------------
PAGE_W <- 8.27; PAGE_H <- 11.69
ML <- 0.85; MR <- 0.85; MT <- 0.85; MB <- 0.90
TEXT_W <- PAGE_W - ML - MR

# ---- palette ---------------------------------------------------------------
COL_TEAL  <- "#0B6E6E"
COL_DARK  <- "#1A1A1A"
COL_GREY  <- "#555555"
COL_ACCENT<- "#B23A48"
COL_RULE  <- "#CCCCCC"
COL_HEADBG<- "#E4F0EF"
COL_ROWBG <- "#F4F8F8"

# ---- state -----------------------------------------------------------------
.cursor <- MT          # inches from top
.page   <- 0
DOC_TITLE <- "Organoid MSI Metabolic-Gradient Study - Methods & Pipeline"

lh_for <- function(fs) fs / 72 * 1.32   # line height (inches) for a font size

start_page <- function() {
  grid.newpage()
  .page <<- .page + 1
  # footer rule + text
  grid.lines(x = unit(c(ML, PAGE_W - MR), "inches"),
             y = unit(c(MB - 0.28, MB - 0.28), "inches"),
             gp = gpar(col = COL_RULE, lwd = 0.7))
  grid.text(DOC_TITLE, x = unit(ML, "inches"), y = unit(MB - 0.42, "inches"),
            just = c("left", "top"), gp = gpar(fontsize = 7, col = COL_GREY))
  grid.text(sprintf("p. %d", .page), x = unit(PAGE_W - MR, "inches"),
            y = unit(MB - 0.42, "inches"), just = c("right", "top"),
            gp = gpar(fontsize = 7, col = COL_GREY))
  .cursor <<- MT
}

room <- function(need) {            # ensure `need` inches available, else new page
  if (.cursor + need > PAGE_H - MB) start_page()
}

# wrap `text` to `width_in` at fontsize/face using real glyph metrics
wrap_text <- function(text, fs, width_in, face = "plain") {
  words <- strsplit(text, "\\s+")[[1]]
  if (length(words) == 0) return("")
  lines <- character(0); cur <- ""
  for (w in words) {
    test <- if (nchar(cur) == 0) w else paste(cur, w)
    wd <- convertWidth(grobWidth(textGrob(test, gp = gpar(fontsize = fs, fontface = face))),
                       "inches", valueOnly = TRUE)
    if (wd > width_in && nchar(cur) > 0) { lines <- c(lines, cur); cur <- w }
    else cur <- test
  }
  c(lines, cur)
}

draw_lines <- function(lines, fs, face = "plain", col = COL_DARK,
                       x_in = ML, lead = 1.0) {
  lh <- lh_for(fs) * lead
  for (ln in lines) {
    room(lh)
    grid.text(ln, x = unit(x_in, "inches"), y = unit(PAGE_H - .cursor, "inches"),
              just = c("left", "top"), gp = gpar(fontsize = fs, fontface = face, col = col))
    .cursor <<- .cursor + lh
  }
}

# ---- element renderers -----------------------------------------------------
add_space <- function(h = 0.10) .cursor <<- .cursor + h

add_title <- function(text, fs = 22) {
  lines <- wrap_text(text, fs, TEXT_W, "bold")
  draw_lines(lines, fs, "bold", COL_TEAL, lead = 1.05)
}

add_subtitle <- function(text, fs = 11) {
  draw_lines(wrap_text(text, fs, TEXT_W, "plain"), fs, "plain", COL_GREY, lead = 1.1)
}

add_h1 <- function(text, fs = 14) {
  add_space(0.16)
  room(lh_for(fs) * 2.2)
  draw_lines(wrap_text(text, fs, TEXT_W, "bold"), fs, "bold", COL_TEAL, lead = 1.05)
  grid.lines(x = unit(c(ML, PAGE_W - MR), "inches"),
             y = unit(PAGE_H - .cursor + 0.02, "inches"),
             gp = gpar(col = COL_TEAL, lwd = 1.1))
  add_space(0.10)
}

add_h2 <- function(text, fs = 11.5) {
  add_space(0.10)
  room(lh_for(fs) * 2.0)
  draw_lines(wrap_text(text, fs, TEXT_W, "bold"), fs, "bold", COL_ACCENT, lead = 1.05)
  add_space(0.02)
}

add_body <- function(text, fs = 9.5) {
  draw_lines(wrap_text(text, fs, TEXT_W, "plain"), fs, "plain", COL_DARK, lead = 1.16)
  add_space(0.06)
}

add_bullet <- function(text, fs = 9.5, indent = 0.26) {
  lh <- lh_for(fs) * 1.16
  lines <- wrap_text(text, fs, TEXT_W - indent, "plain")
  for (i in seq_along(lines)) {
    room(lh)
    if (i == 1)
      grid.circle(x = unit(ML + 0.10, "inches"),
                  y = unit(PAGE_H - .cursor - lh_for(fs) * 0.42, "inches"),
                  r = unit(0.022, "inches"),
                  gp = gpar(fill = COL_TEAL, col = NA))
    grid.text(lines[i], x = unit(ML + indent, "inches"),
              y = unit(PAGE_H - .cursor, "inches"), just = c("left", "top"),
              gp = gpar(fontsize = fs, col = COL_DARK))
    .cursor <<- .cursor + lh
  }
  add_space(0.02)
}

# three/n-column table with wrapping, shaded header, zebra rows, page-aware
add_table <- function(header, rows, widths, fs = 8.6, pad = 0.05) {
  stopifnot(length(header) == length(widths))
  widths <- widths / sum(widths)
  colw <- widths * TEXT_W
  xstarts <- ML + c(0, head(cumsum(colw), -1))
  lh <- lh_for(fs) * 1.12

  draw_row <- function(cells, face, bg = NA) {
    wrapped <- lapply(seq_along(cells), function(j)
      wrap_text(as.character(cells[j]), fs, colw[j] - 2 * pad, face))
    nlines <- max(vapply(wrapped, length, 1L))
    rh <- nlines * lh + 2 * pad
    room(rh)
    if (!is.na(bg))
      grid.rect(x = unit(ML, "inches"), y = unit(PAGE_H - .cursor, "inches"),
                width = unit(TEXT_W, "inches"), height = unit(rh, "inches"),
                just = c("left", "top"), gp = gpar(fill = bg, col = NA))
    for (j in seq_along(cells)) {
      yy <- .cursor + pad
      for (ln in wrapped[[j]]) {
        grid.text(ln, x = unit(xstarts[j] + pad, "inches"),
                  y = unit(PAGE_H - yy, "inches"), just = c("left", "top"),
                  gp = gpar(fontsize = fs, fontface = face,
                            col = if (identical(face, "bold")) COL_TEAL else COL_DARK))
        yy <- yy + lh
      }
    }
    grid.lines(x = unit(c(ML, ML + TEXT_W), "inches"),
               y = unit(PAGE_H - (.cursor + rh) + 0.005, "inches"),
               gp = gpar(col = COL_RULE, lwd = 0.5))
    .cursor <<- .cursor + rh
  }

  hdr <- function() draw_row(header, "bold", COL_HEADBG)
  room(lh * 3); hdr()
  for (i in seq_along(rows)) {
    if (.cursor + lh * 2 > PAGE_H - MB) { start_page(); hdr() }
    draw_row(rows[[i]], "plain", if (i %% 2 == 0) COL_ROWBG else NA)
  }
  add_space(0.10)
}

# ---- schematic pipeline diagram (own page) ---------------------------------
render_schematic <- function() {
  start_page()
  fsH <- 14
  draw_lines(wrap_text("Pipeline at a glance", fsH, TEXT_W, "bold"), fsH, "bold", COL_TEAL, lead = 1.05)
  grid.lines(x = unit(c(ML, PAGE_W - MR), "inches"),
             y = unit(PAGE_H - .cursor + 0.02, "inches"), gp = gpar(col = COL_TEAL, lwd = 1.1))
  add_space(0.06)
  draw_lines(wrap_text("From 20 raw MALDI-MSI datasets to the apical-out citrate result. Phase numbers refer to R/<NN> folders; Steps A-E match the manuscript methods (E - citrate identification - underpins the citrate measurement used in D).",
                       8.5, TEXT_W, "plain"), 8.5, "plain", COL_GREY, lead = 1.12)

  cx <- ML + TEXT_W / 2
  MW <- 4.9

  box <- function(top, h, fill, border, title, sub = NULL, tcol = "#FFFFFF", w = MW) {
    left <- cx - w / 2
    grid.rect(x = unit(left, "inches"), y = unit(PAGE_H - top, "inches"),
              width = unit(w, "inches"), height = unit(h, "inches"), just = c("left", "top"),
              gp = gpar(fill = fill, col = border, lwd = 1.3))
    grid.text(title, x = unit(cx, "inches"), y = unit(PAGE_H - (top + 0.12), "inches"),
              just = c("center", "top"), gp = gpar(fontsize = 10, fontface = "bold", col = tcol))
    if (!is.null(sub))
      grid.text(sub, x = unit(cx, "inches"), y = unit(PAGE_H - (top + 0.37), "inches"),
                just = c("center", "top"), gp = gpar(fontsize = 8, col = tcol))
  }
  arr <- function(y1, y2, x = cx) {
    grid.lines(x = unit(c(x, x), "inches"), y = unit(c(PAGE_H - y1, PAGE_H - y2), "inches"),
               arrow = arrow(angle = 22, length = unit(0.09, "inches"), type = "closed"),
               gp = gpar(col = COL_GREY, lwd = 1.4, fill = COL_GREY))
  }

  y <- .cursor + 0.22
  GAP <- 0.20
  step  <- function(h, fill, border, title, sub, tcol = "#FFFFFF") {
    box(y, h, fill, border, title, sub, tcol); y <<- y + h
  }
  gap   <- function() { arr(y, y + GAP); y <<- y + GAP }
  label <- function(txt, col = COL_ACCENT) {
    y <<- y + 0.10
    grid.text(txt, x = unit(cx, "inches"), y = unit(PAGE_H - y, "inches"),
              just = c("center", "top"), gp = gpar(fontsize = 9.3, fontface = "bold", col = col))
    y <<- y + 0.22
  }

  step(0.56, "#707070", "#555555", "Raw MALDI-MSI  -  20 datasets",
       "negative mode  |  m/z 100-900  |  10 um pixel"); gap()
  step(0.66, COL_TEAL, "#085656", "Phase 01  -  Spectral processing",
       "TIC norm | 25 ppm grid | dedup / deisotope  ->  8,772 x 90,760"); gap()
  step(0.58, "#0E8A8A", "#085656", "Phase 02  -  Citrate standard (gate)",
       "anchor m/z 191.01976  |  +/- 7 ppm"); gap()
  step(0.62, "#CDE7E6", COL_TEAL, "E  -  Citrate identification (Phases 02 & 07)",
       "accurate mass + authentic standard + LC-MS/MS  |  MSI Level 2  (underpins D)",
       tcol = COL_DARK); gap()

  label("Spatial & statistical analysis  -  Steps A-D")
  sc <- COL_ACCENT; sf <- "#FBEDEF"
  step(0.56, sf, sc, "A  -  Spatial clustering (Phase 04)",
       "SSC on-tissue mask | 348 ions | floor80", tcol = COL_DARK); gap()
  step(0.56, sf, sc, "B  -  Co-registration (Phases 03 / 05 / 06)",
       "MSI <-> brightfield <-> IF | affine + landmarks", tcol = COL_DARK); gap()
  step(0.56, sf, sc, "C  -  Apical annotation (Phases 08-10)",
       "segment + split | manual in / out / mixed", tcol = COL_DARK); gap()
  step(0.56, sf, sc, "D  -  Citrate statistics (Phase 11)",
       "near-field 50 um | Wilcoxon + Cliff's delta", tcol = COL_DARK); gap()
  step(0.60, COL_ACCENT, "#7E2833", "Result",
       "apical-out > basolateral-out near-field citrate  |  p = 6.1e-5, delta = +0.62")

  # footnote on citrate identity / orthogonal evidence
  y <- y + 0.26
  grid.lines(x = unit(c(ML, ML + 2.0), "inches"),
             y = unit(PAGE_H - y, "inches"), gp = gpar(col = COL_RULE, lwd = 0.7))
  .cursor <<- y + 0.07
  draw_lines(wrap_text(paste0(
    "Citrate is annotated at MSI Level 2 (m/z confirmed vs a matrix-matched authentic standard); ",
    "LC-MS/MS of conditioned media corroborates secretion but does not co-confirm the imaging feature, ",
    "and m/z 191.0217 does not resolve citrate from isocitrate. See Step E (Section 9)."),
    7.6, TEXT_W, "italic"), 7.6, "italic", COL_GREY, lead = 1.14)
}

# =============================================================================
# RENDER
# =============================================================================
pdf(OUT_PDF, width = PAGE_W, height = PAGE_H, onefile = TRUE)
start_page()

# ---------- TITLE BLOCK ----------
add_space(0.6)
add_title("MALDI-MSI Metabolic-Gradient Analysis")
add_space(0.05)
add_title("Methods & Pipeline Report", fs = 16)
add_space(0.18)
add_subtitle("Intestinal-organoid sections incubated in CMC for at least 5 min before freezing (all sections treated identically; no incubation-time comparison). MALDI-timsTOF flex imaging in negative mode; question: does organoid apical polarity (apical-out vs basolateral-out) relate to near-field citrate emission into the surrounding gel?")
add_space(0.20)
grid.lines(x = unit(c(ML, ML + 2.4), "inches"),
           y = unit(PAGE_H - .cursor, "inches"), gp = gpar(col = COL_ACCENT, lwd = 2))
add_space(0.18)
add_body("Dataset: 20 imaging datasets (10 sections x 2 slides, sl6A and sl4A; all treated identically). Final processed feature table: 8,772 features x 90,760 pixels.")
add_body("Software: R 4.4.2, Cardinal v3 (Bioconductor), matter, BiocParallel (SnowParam, 4 workers), RBioFormats/EBImage, RANN, viridisLite. Pipeline organized as 11 numbered phase folders under R/ plus shared libraries in R/00_lib/ (locked constants in lib_paths.R).")
add_body(paste0("Generated from branch 'restructure-pipeline'. This report documents spectral processing and the five analytical steps (A-E) described in the manuscript methods, mapped onto the 11 computational pipeline phases."))

# ---------- SCHEMATIC (page of its own) ----------
render_schematic()

# ---------- 1. OVERVIEW ----------
start_page()
add_h1("1.  Study design & data acquisition")
add_body("Intestinal organoids embedded in 4% carboxymethylcellulose (CMC) were sectioned and analyzed by matrix-assisted laser desorption/ionization mass spectrometry imaging (MALDI-MSI) on a Bruker timsTOF flex. Organoids were incubated in CMC for at least 5 minutes prior to freezing; all sections were treated identically (no incubation-time comparison). Two physical slides (sl6A, sl4A) carry 10 sections each, for 20 datasets in total.")
add_h2("Sample preparation & matrix")
add_bullet("Embedding medium: 4% carboxymethylcellulose (CMC).")
add_bullet("Matrix: 1,5-diaminonaphthalene (DAN), for negative-ion-mode small-metabolite detection.")
add_h2("MALDI imaging acquisition")
add_bullet("Polarity: negative ion mode.")
add_bullet("Recorded mass range: m/z 100-900.")
add_bullet("Spatial (raster) step size: 10 um pixel pitch.")
add_bullet("Laser: custom setting with 'Beam Scan' on, ~6 um scan range (spot size) and 10 um resulting field size, chosen to avoid oversampling.")
add_bullet("Ion optics / transfer settings as configured on the timsTOF flex method (MS1 only, no fragmentation; in-source CID energy 0 eV). See the instrument-parameter summary in Section 9.")

# ---------- 2. SPECTRAL PROCESSING ----------
add_h1("2.  Spectral processing  (Phase 01 - R/01_preprocess)")
add_body("All processing was performed in Cardinal v3. The 20 individual datasets were read as MS-imaging arrays, concatenated into a single combined experiment, and total-ion-current (TIC) normalized. A common reference m/z grid was then built, the data binned to that grid, and the feature list reduced by a sequence of filters to a final, self-contained peak table.")

add_h2("2.1  Combination & normalization")
add_bullet("Each .imzML dataset read with readMSIData() as an MSImagingArrays object.")
add_bullet("Datasets combined into one experiment (element-wise concatenation), then TIC-normalized so that per-pixel total signal is comparable across sections and slides.")

add_h2("2.2  Reference m/z grid")
add_bullet("Because Cardinal v3's automatic reference-peak estimation is unreliable on centroided data, the reference grid was built by manual single-linkage clustering of centroid m/z values.")
add_bullet("Density-invariant seeding: the grid is estimated from a fixed random subsample of 15,000 pixels (seeded for reproducibility) so the peak list does not collapse/over-chain as total pixel count grows.")
add_bullet("Clustering tolerance: 25 ppm single linkage defines peak boundaries.")

add_h2("2.3  Binning & alignment")
add_bullet("Spectra converted from arrays to a unified-m/z experiment at 25 ppm tolerance.")
add_bullet("Peaks aligned (peakAlign) at 25 ppm to merge the expanded grid into consistent feature centroids. Peak detection used a signal-to-noise threshold of 3.")

add_h2("2.4  Feature filtering")
add_bullet("Frequency filter (singletons): drop features detected in <= 0.1% of pixels (freq > 0.001 kept).")
add_bullet("Intensity cutoff: a data-driven elbow (Kneedle algorithm) on the ranked mean-intensity curve removes low-intensity noise features.")
add_bullet("Spatial-coverage floor: retain features present in >= 5% of pixels (freq >= 0.05), supporting per-zone gradient statistics.")

add_h2("2.5  Deduplication & deisotoping")
add_bullet("Feature deduplication: 5 ppm single-linkage with a 10 ppm diameter cap, applied before any per-ion ranking.")
add_bullet("13C deisotoping: an M+1 peak is removed when its intensity ratio to the monoisotopic peak is within [0.20, 0.80] and the two share a mean spatial correlation r > 0.85 (conservative).")

add_h2("2.6  Final processed dataset")
add_bullet("8,772 features x 90,760 pixels across all 20 sections.")
add_bullet("Stored as a self-contained in-memory sparse matrix (cache/peaks_combined.rds), with no dependency on the original .ibd cache.")

# ---------- 3. CITRATE STANDARD ----------
add_h1("3.  Citrate standard & mass anchor  (Phase 02 - R/02_CitrateStandard)")
add_body("Because citrate is the focal analyte, a matrix-matched citrate dilution series was run to validate its mass position, detection window, and quantitative behavior before any tissue analysis (an early gate in the pipeline).")
add_bullet("Target ion: citrate [M-H]-, theoretical m/z 191.0197.")
add_bullet("Measured/locked anchor: m/z 191.01976, used for all downstream citrate extraction.")
add_bullet("Extraction window: +/- 7 ppm (selected from a 5/7/10 ppm sweep to maximize citrate capture with no neighboring-ion contamination).")
add_bullet("Calibration: log-log slope ~1.00 over ~1-100 mM (R2 ~ 0.99); estimated limit of detection ~0.37 mM.")
add_bullet("Citrate is extracted per pixel from the raw centroid data (function citrate_onto_pd() in R/00_lib/lib_citrate.R) rather than from the 25-ppm-binned grid feature, to avoid blending citrate with a near-isobaric neighbor.")
add_body("Note: m/z 191.0217 (C6H7O7-) is shared by citrate and its structural isomer isocitrate, which are not separable by accurate mass. See Step E (Section 9) for the identity-confidence assessment and the orthogonal LC-MS/MS evidence.")

# ---------- 4. STEP A: SPATIAL CLUSTERING ----------
add_h1("4.  Step A - Spatial clustering to delineate organoid tissue  (Phase 04 - R/04_ssc_ontissue)")
add_body("Organoid tissue areas were separated from surrounding gel using spatial shrunken centroids (SSC) clustering, computed on a curated set of on-tissue ions rather than on the full feature table.")
add_h2("4.1  Curated ion set")
add_bullet("348 ions = 149 published 'known' metabolites union 237 manually selected on-tissue 'unknown' ions. Clustering uses only these 348 ions.")
add_h2("4.2  SSC parameters (Cardinal v3)")
add_bullet("Spatial smoothing radius r = 2 pixels (= 20 um at 10 um pitch), adaptive weights.")
add_bullet("Two passes: a coarse pass at k = 4 for the tissue/gel split, and a finer pass at k = 10 for substructure.")
add_bullet("Shrinkage parameter s tested over {6, 9, 12}; s = 9 used downstream.")
add_h2("4.3  On-tissue rule ('floor80')")
add_bullet("A pixel is on-tissue if its k=4 cluster mean curated signal is >= 50% of the richest cluster, OR its own curated signal exceeds the section's 80th percentile.")
add_bullet("Rationale: clusters alone miss low-signal epithelial edges; the 80th-percentile floor recovers them without flooding the lumen. A morphological-dilation buffer was tested and rejected (it floods the lumen).")
add_h2("4.4  Cross-section harmonization")
add_bullet("Per-section local clusters pooled and merged by hierarchical clustering (Ward.D2) on a 1 - Pearson distance of log1p mean spectra, with the cluster count chosen by silhouette over k = 4-10.")
add_bullet("Result: 4 tissue chemotypes plus an is_tissue mask, written to cache/peaks_tissue_combined.rds.")

# ---------- 5. STEP B: CO-REGISTRATION ----------
add_h1("5.  Step B - Co-registration of brightfield & MS images  (Phases 03 & 05)")
add_body("Each MS image was aligned to its high-resolution microscopy image through a chain of affine transforms, coarse first then refined, so that ion images and optical morphology share a common coordinate frame.")
add_h2("5.1  Coarse registration (Phase 03 - R/03_coarse_registration)")
add_bullet("The MSI raster is mapped to the slide-overview image using the stage Area bounding-box and teach points stored in the Bruker acquisition method (.mis XML).")
add_bullet("Global orientation (axis flips) is resolved automatically by maximizing overlap of the is_tissue mask with the optical tissue.")
add_h2("5.2  Refinement to native microscopy (Phase 05 - R/05_registration_refine)")
add_bullet("A within-box search refines translation/rotation/scale on the slide image, then the slide image is mapped to the native-resolution .nd2 brightfield (vertical flip, ~2.895x scale).")
add_bullet("A final small polish (<= 3 MSI pixels) aligns all 20 sections; per-section 3x2 affine transforms and native crops are cached.")
add_h2("5.3  Immunofluorescence registration (Phase 06 - R/06_if_registration)")
add_bullet("Serial immunofluorescence sections (channels: F-actin, ZO-1, beta-catenin, DAPI) are registered to the brightfield/MSI frame by manually placed organoid-center landmarks.")
add_bullet("10 of 12 sections registered (5-8 organoid landmark pairs each); organoid-level RMSE 28-222 um, the realistic ceiling for serial sections. Pixel sizes taken from .nd2 metadata (includes the 1.5x optical changer).")

# ---------- 6. METABOLITE ID ----------
add_h1("6.  Metabolite identification & confidence  (Phase 07 - R/07_metabolite_id)")
add_body("Processed features were matched against a published metabolite list (149 knowns) by accurate mass (m/z), with citrate-resolution QC confirming citrate is cleanly resolved within the +/- 7 ppm window. No on-tissue fragmentation (MS/MS) and no chromatographic separation were performed.")
add_h2("6.1  Identification confidence (Metabolomics Standards Initiative)")
add_body("Annotation confidence is reported on the four-level scheme of the Metabolomics Standards Initiative (Sumner et al., Metabolomics 2007; 3:211-221). Level 1 = identified compound (>= 2 orthogonal properties matched to an authentic standard analyzed under identical conditions); Level 2 = putatively annotated compound (accurate-mass or spectral match to a database, no authentic standard); Level 3 = putatively characterized compound class; Level 4 = unknown.")
add_bullet("Database-matched metabolites (the 149 knowns): MSI Level 2 - putative annotation by accurate mass alone (a single property; no MS/MS or chromatographic dimension).")
add_bullet("Citrate ([M-H]- m/z 191.0217): identified at MSI Level 2 (standard-confirmed). The full identification, the orthogonal LC-MS/MS evidence, and the citrate/isocitrate caveat are detailed in Step E (Section 9).")
add_bullet("The 237 manually selected on-tissue ions used for clustering: MSI Level 4 - unknown features, carried for spatial analysis only.")

# ---------- 7. STEP C: APICAL ANNOTATION ----------
add_h1("7.  Step C - Manual apical-orientation annotation  (Phases 08-10)")
add_body("Individual organoids were isolated, curated, and manually annotated for apical orientation (basolateral-out vs apical-out vs mixed) using the registered microscopy.")
add_h2("7.1  Segmentation & organoid curation (Phases 08-09)")
add_bullet("Organoids segmented as connected components on the tissue mask. Touching/merged organoids are split by user-drawn freehand lines on a review PDF; strokes accumulate in a persistent store so re-rendering never loses a cut. Curated ROIs are written to cache/instances_final_<sid>.rds.")
add_h2("7.2  Apical annotation (Phase 10 - R/10_apical_annotation)")
add_bullet("Organoid outlines are rendered on the native brightfield; the annotator adds a text comment (in / out / mixed) on or beside each numbered organoid.")
add_bullet("Comments are extracted from the PDF and matched to the nearest organoid centroid (<= 25 pt). Conflicting comments on one organoid are forced to 'mixed'; unmatched comments are logged, never silently dropped.")
add_bullet("Output: results/annotation/ with a per-organoid apical class.")

# ---------- 8. STEP D: STATISTICS ----------
add_h1("8.  Step D - Statistical analysis of citrate  (Phase 11 - R/11_per_organoid_final)")
add_body("Citrate abundance in the gel immediately surrounding each organoid was compared between apical-out and basolateral-out organoids, testing whether apical-out organoids release more citrate into the near-field.")
add_h2("8.1  Citrate measurement")
add_bullet("Citrate ([M-H]-, m/z 191.0217) extracted per pixel from raw centroids within +/- 7 ppm and TIC-normalized.")
add_bullet("Signed Euclidean distance of each gel pixel from the organoid surface computed with RANN::nn2.")
add_h2("8.2  Near-field comparison (headline)")
add_bullet("Per organoid, the mean near-field citrate within 50 um of the surface (grad_near50) is the headline metric. The replication unit is the individual organoid, collapsing pixels to one value per organoid to avoid spatial pseudo-replication.")
add_bullet("apical-out vs basolateral-out compared by Wilcoxon rank-sum test with Cliff's delta as effect size: p = 6.1e-5, delta = +0.62 (large).")
add_bullet("Corroborated by the log2 citrate/DHA ratio (an internal control that cancels shared thickness/ionization effects): p = 1.4e-8, significant at both 0 h and 20 h.")
add_h2("8.3  Per-organoid gradient (supporting)")
add_bullet("Outward gradients were also fit per organoid (slope of log1p intensity vs distance; Spearman rho; decay length) and tested across organoids with a one-sample Wilcoxon signed-rank test, with BH FDR < 0.10 for per-organoid calls.")
add_bullet("Scope (locked): the design supports organoid-level evidence within these slides. Single-condition study (all sections incubated in CMC >=5 min, treated identically); no incubation-time comparison.")

# ---------- 9. STEP E: CITRATE IDENTIFICATION ----------
add_h1("9.  Step E - Citrate identification  (Phases 02 & 07)")
add_body("Because citrate is the focal analyte, its identity was established before any quantitative comparison. Citrate was detected as the deprotonated ion [M-H]- at m/z 191.0217 in negative-ion mode. A matrix-matched authentic citrate standard dilution series was analyzed under identical MALDI conditions (Phase 02): this confirmed the exact mass (measured anchor m/z 191.01976, ~0.2 ppm from theory), fixed the +/- 7 ppm extraction window, and established a linear quantitative response (log-log slope ~1.0, R2 ~0.99, limit of detection ~0.37 mM). In parallel, quantitative LC-MS/MS of conditioned media from the same organoid line (with derivatization and fragmentation, a separate experiment) independently confirmed that this line secretes citrate.")
add_body("Within the Metabolomics Standards Initiative scheme (Sumner et al., Metabolomics 2007), the imaging annotation corresponds to Level 2 (putatively annotated, standard-confirmed): citrate is matched by accurate mass to an authentic standard run under identical conditions, but only one orthogonal property (m/z) is available on the imaging platform. It does not reach Level 1, which would require a second orthogonal identifier co-measured on the imaging feature itself - for example, on-tissue MS/MS fragmentation or ion-mobility (CCS) - matched to the standard. The media LC-MS/MS strengthens the biological case for citrate secretion but, being a different platform and sample, does not co-confirm the imaging feature.")
add_body("Finally, m/z 191.0217 (C6H7O7-) is shared by citrate and its structural isomer isocitrate and cannot distinguish them by mass alone (accurate mass and on-tissue MS/MS generally do not separate the isomers; chromatography does). The conditioned-media LC-MS/MS indicates citrate as the secreted species, so the ion is reported as citrate throughout, without excluding an isocitrate contribution to the imaging signal.")

# ---------- 10. PARAMETER SUMMARY TABLE ----------
add_h1("10.  Instrument & pipeline parameter summary")
add_h2("10.1  MALDI-MSI acquisition (timsTOF flex, MS1 only)")
add_table(
  header = c("Parameter", "Value"),
  widths = c(0.55, 0.45),
  rows = list(
    c("Polarity", "Negative ion mode"),
    c("Mass range (recorded)", "m/z 100-900"),
    c("In-source CID (isCID) energy", "0 eV (no fragmentation)"),
    c("MALDI plate offset", "50 V"),
    c("Deflection 1 delta", "-70 V"),
    c("Funnel 1 / Funnel 2 RF", "125 / 200 Vpp"),
    c("Multipole RF", "200 Vpp"),
    c("Quadrupole ion energy / low mass", "8.0 eV / 120 m/z"),
    c("Collision energy / collision RF", "5.0 eV / 500 Vpp"),
    c("Transfer time / pre-pulse storage", "50 us / 5 us"),
    c("High-Sensitivity Detection / Focus Mode", "Disabled"),
    c("Pixel (raster) size", "10 um"),
    c("Laser field / spot (Beam Scan)", "10 um field / ~6 um spot"),
    c("Matrix / embedding", "DAN / 4% CMC")
  ))

add_h2("10.2  Spectral processing & analysis")
add_table(
  header = c("Component", "Parameter", "Value"),
  widths = c(0.30, 0.40, 0.30),
  rows = list(
    c("Peak detection", "Signal-to-noise threshold", "3"),
    c("Reference grid", "Single-linkage tolerance", "25 ppm"),
    c("Reference grid", "Seed subsample", "15,000 pixels (fixed)"),
    c("Binning / align", "Tolerance", "25 ppm"),
    c("Filtering", "Singleton freq floor", "> 0.001"),
    c("Filtering", "Intensity cutoff", "Kneedle elbow (auto)"),
    c("Filtering", "Spatial coverage floor", ">= 0.05 (5% px)"),
    c("Dedup", "Linkage / diameter cap", "5 ppm / 10 ppm"),
    c("Deisotope", "13C ratio / spatial r", "[0.20, 0.80] / r > 0.85"),
    c("Final table", "Features x pixels", "8,772 x 90,760"),
    c("SSC (Step A)", "radius r / passes k", "2 px / 4 then 10"),
    c("SSC (Step A)", "shrinkage s", "{6, 9, 12}; 9 used"),
    c("SSC (Step A)", "on-tissue rule", ">=50% richest OR >=80th pct"),
    c("SSC (Step A)", "curated ions / chemotypes", "348 / 4"),
    c("Citrate", "anchor m/z / window", "191.01976 / +/- 7 ppm"),
    c("Citrate", "calibration / LOD", "slope ~1.0, R2 ~0.99 / ~0.37 mM"),
    c("Registration (B)", "MSI->BF transform", "per-section 3x2 affine"),
    c("Registration (B)", "refine polish / IF RMSE", "<=3 MSI px / 28-222 um"),
    c("Apical (C)", "scoring / matching", "manual in/out/mixed, <=25 pt"),
    c("Stats (D)", "near-field window", "50 um"),
    c("Stats (D)", "test / effect size", "Wilcoxon rank-sum / Cliff's delta"),
    c("Stats (D)", "headline result", "p = 6.1e-5, delta = +0.62"),
    c("Stats (D)", "replication unit", "individual organoid")
  ))

# ---------- 10. REPRODUCIBILITY ----------
add_h1("11.  Reproducibility")
add_bullet("Environment: R 4.4.2 with Cardinal v3; parallelism via SnowParam (4 workers). Locked constants centralized in R/00_lib/lib_paths.R.")
add_bullet("Runners: run_all.sh executes all computational phases in order; regen_reports.sh rebuilds every report PDF from cached .rds; run_apical_nearfield_pipeline.sh reproduces the apical near-field analysis and figure.")
add_bullet("Manual gates (on-tissue ion selection, organoid splitting, apical comments) are cached, so downstream steps run automatically; gates are marked '>>> MANUAL <<<'.")
add_bullet("This report: figures/methods_report/generate_methods_report.R -> methods_pipeline_report.pdf.")

invisible(dev.off())
cat("Wrote", OUT_PDF, "(", .page, "pages )\n")
