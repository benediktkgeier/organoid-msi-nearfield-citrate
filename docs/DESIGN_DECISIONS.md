# Locked design decisions

Canonical list of the project's locked methodological choices. These are settled — change only with
explicit reason.

## Acquisition

- MALDI-timsTOF flex, no ion mobility, negative mode, *m/z* 100–900, DAN matrix, 10 µm pixel,
  10 µm sections, 4% CMC embedding. Pre-calibration: red phosphorus 0.99 ppm.
- Online lock-mass calibration, 5 anchors: taurine 124.0068, palmitate FA16:0 [M−H] 255.2330,
  AMP 346.0558, ADP 426.0222, taurocholate 514.2844 (`ALIGN_ANCHORS` in `R/00_lib/lib_paths.R`).
- High-res microscopy: Nikon `.nd2` (transmitted light, single channel, ~23K × 11K px) for figure
  backgrounds at **native resolution** (never downsampled); low-res `_small.jpg` for landmark clicking.

## Feature pipeline (phase 01)

- Mass tolerance everywhere: **25 ppm** (`PPM_TOL`) — empirical bin-width optimum for this instrument.
- Manual single-linkage reference grid (Cardinal v3 `estimateReferencePeaks()` is profile-tuned and
  fails on centroid data) → `convertMSImagingArrays2Experiment` → `peakAlign(25 ppm)` →
  `summarizeFeatures`. **Never** `setCardinalBPPARAM()` (v3 bug); use `register_parallel()`.
- Feature dedup: 5 ppm single-linkage + 10 ppm diameter cap. 13C deisotope: M+1 drop when ratio ∈
  [0.20, 0.80] and mean spatial r > 0.85 (conservative). `freq ≥ 0.05` spatial floor on top of the
  Kneedle cut → **8,772 features** (`cache/peaks_combined.rds`, self-contained, no `.ibd` dep).

## On-tissue ion selection + delineation (phase 04)

- **Selected on-tissue ions:** a manually curated list of on-tissue ions (237) selected from a
  candidate ion sheet by a scientist. (Stored under `results/peakme_annotations/` — legacy dir name.)
- **Curated analysis set = 348 ions** = 149 published-list "known" metabolites (isobars collapse to
  145 features) ∪ 237 manually selected on-tissue "unknown" ions (34 overlap). Built from the recall
  grid to include 23 sub-floor knowns. → `cache/peaks_curated.rds` (self-contained).
- **On-tissue rule "floor80":** clustering uses ONLY the 348 curated ions. Per section: Pass-1 SSC
  k=4 (s=6/9/12, r=2). A pixel is on-tissue iff its k4 cluster mean curated-signal ≥ 50% of the
  richest cluster's (`TISSUENESS_THR=0.5`) **OR** its own curated-signal ≥ the section's 80th
  percentile (`FLOOR_PCT=0.80`). A morphological buffer was tried and **rejected** (floods lumen).
- Pass-2 SSC k=10 within on-tissue → per-section substructure → cross-section **chemotype**
  harmonization. Result: on-tissue 18,610 px (20.5%); **4 chemotypes**
  (silhouette), user-confirmed.

## Gradient (phases 08 / 11)

- Segmentation: connected components (`EBImage::bwlabel`, 4-conn) of `is_tissue`; drop < 50 px.
- Signed distance from surface via `RANN::nn2` × `MSI_PIXEL_UM` (10). Outward into CMC at
  `BUF_LEVELS_UM` = **10/20/50/80/160 µm** (Voronoi catchment: drop pixels nearer another instance);
  inward at 10 µm steps. Ring intensity = **mean per pixel = area-normalized**.
- ρ_out = Spearman(zone idx, mean intensity), POOLED across all sections (single condition; no group
  split, no Δρ). Near-field band 50–100 µm. Headline near-field metric = `grad_near50` (≤50 µm).
- Single-condition study: all sections were incubated in CMC for at least 5 min before freezing and
  are treated identically (no incubation-time comparison). Contrast of interest = apical-out vs
  basolateral-out polarity (colours locked in `lib_paths.R::APICAL_COLS`).

## Statistics

- organoid = unit → **descriptive only, NO p-values / NO LMM** for group contrasts; per-organoid /
  per-class work uses pairwise Wilcoxon + Cliff's δ effect size. The apical per-class analysis (Phase
  11) runs over the two-annotator consensus set, **n = 84** organoids (30 apical-out / 23
  basolateral-out / 31 mixed; `results/annotation/apical_map_consensus.csv`). (The earlier pooled
  20-dataset gradient survey covered 72 segmented organoids pre-refinement.)

## Citrate tracer (phases 02 / 07 / 08 / 11)

- Tracer = strict **[M−H]⁻**, **anchored** to the standard-measured mass
  `CITRATE_ANCHOR_MZ = 191.01976` with a **±7 ppm** window (`CITRATE_WIN_PPM`, both locked in
  `lib_paths.R` by the Phase 02 citrate-standard GATE). ALL downstream citrate extraction reads the
  anchored per-pixel citrate via `citrate_onto_pd()` in `R/00_lib/lib_citrate.R` (raw imzML, TIC-
  normalised), **not** the old processed-grid feature 191.0217. **[M−H]− only** — the earlier
  "adduct-switching / total-citrate" (Na/K/Cl gel adducts tracking outward diffusion) was found to be
  an **artifact and was removed**. DHA 327.2330 = confined-lipid negative control.
- Caveat (centroided data): a narrow window **selects** citrate-dominant pixels, it does not
  integrate citrate area.
- **191 QC:** the matched processed-grid feature is the centroid of an **unresolved blend** at +8–10 ppm; true
  citrate (191.0197) is a real, abundance-variable component (3.4× with-vs-without, p=0.008) but
  cannot be cleanly separated from the isobar at 191.0211. To isolate citrate, keep the integration
  upper edge ≤ 191.0204 (≈ 191.0197 ±3 ppm); ~25–30% residual contamination is overlap-limited.

## Rendering (locked)

- Ion images: **viridis, linear scale, global p99.5 clip per ion, gamma 1.0**, NA → 0 (darkest
  viridis). 100 µm scale bar below-right. **Every ion image carries a viridis intensity colorbar**
  labelled with its clip range.
- MSI-on-brightfield overlay: homogeneous **constant opacity (0.65)** across the whole MSI image (not
  intensity-thresholded); microscopy at native `.nd2` resolution.
