# SSC on-tissue delineation, chemotypes & organoid segmentation (Phase 04 + 08/09)

How the datasets are delineated ("what is tissue") and segmented into organoids, and the visual
report that documents it. All clustering is computed once by `R/04_ssc_ontissue/04_ssc_tissue_mask.R`
and persisted into `cache/peaks_tissue_combined.rds`; the visual report `05_ssc_report.R` only reads
those cached columns (no recompute).

## Method

**Feature set.** 348 curated on-tissue ions (`cache/peaks_curated.rds`; 149 published knowns + 237
on-tissue unknowns).

**Tissueness proxy.** Per pixel, sum of the 348 curated-ion intensities, clipped to its p99.5 → a 0..1
richness score (`tissueness_px`). Uses only curated ions so matrix/gel ions do not inflate it.

**SSC pass 1 (per section).** `Cardinal::spatialShrunkenCentroids(r = 2, k = 4, s ∈ {6,9,12},
weights = "adaptive")`, each of the 20 sections clustered independently (`set.seed(42)`). The s = 9
labelling drives the tissue cut. Column: `ssc_k4_sec` (and `ssc_k4_sec_s6/_s9/_s12`).

**On-tissue cut ("floor80").** Keep the k4 clusters whose mean `tissueness_px` is ≥ 50 % of the richest
cluster's mean (the richest is always kept), then **union** with *floor80* = pixels at or above the
section's 80th-percentile tissueness. The union recovers epithelium edges that get lumped into
gel-dominated clusters **without** flooding the low-signal lumen. Result → logical `is_tissue`.
(A morphological-buffer alternative was tested and **rejected** — it floods the lumen.)

**SSC pass 2 (within tissue).** k = 10 local substructure over on-tissue pixels only (`ssc_k10_sec`).
Per-section pass-2 centroids (mean spectrum over the 348 features) are pooled and **harmonized**
across sections by `hclust(ward.D2)` on `1 − cor(log1p(centroids))`; K is chosen by maximum mean
silhouette over 4..10 (override via `cache/chemotype_k_override.rds`) → `chemotype` (integer, NA
off-tissue; shared identity across sections).

**Organoid segmentation (Phase 08 → 09).** `is_tissue` is the organoid footprint: 4-connectivity
connected components (`EBImage::bwlabel`), drop components < 50 px, relabel → `cache/instances_<sid>.rds`.
Phase 09 refines ("one connected ROI = one id", manual split/cleanup) → `cache/instances_final_<sid>.rds`
(fields `gidx, x, y, is_tissue, instance, is_surface`). `11/01_zones_curated.R` then builds
`cache/zones_<sid>.rds` with the signed distance field (`signed_dist_um`, RANN nn2, 10 µm/px; negative
inward, positive outward), surface flag, and Voronoi catchment used by the near-field gradient work.

## Cached columns (in `peaks_tissue_combined.rds` pixelData)

`tissueness_px` (0..1), `is_tissue` (logical), `ssc_k4_sec` + `_s6/_s9/_s12`, `ssc_k10_sec` + `_s*`
(NA off-tissue), `chemotype` (NA off-tissue). Downstream analysis consumes only `is_tissue`; the
cluster/chemotype columns exist for QC and are read by the visual report below.

## Reports (`figures/ssc/`)

- **`ssc_mask_<sid>.pdf`** (per section, from `04_ssc_tissue_mask.R`): k4 maps (s-sweep), tissueness vs
  on-tissue mask, k10 maps, chemotype map. Plus `chemotype_harmonization.pdf` (dendrogram + centroid
  correlation heatmap + per-section composition).
- **`ssc_clustering_segmentation_report.pdf`** (from `05_ssc_report.R`): a single presentable report.
  Page 1 = methods writeup + chemotype-composition stacked bar across the 20 sections + on-tissue pixel
  counts. Pages 2–21 = per dataset, compact 4-panel (shared-frame style of the final gradient report):
  1. tissueness (viridis) with the white-dotted `is_tissue` outline,
  2. SSC pass-1 k = 4 clusters,
  3. on-tissue mask (white = tissue),
  4. organoid instance segmentation (per-instance coloured outline + red surface pixels + instance id).

Run (cache-only, R 4.4.2): `Rscript R/04_ssc_ontissue/05_ssc_report.R all` (or a single `<sid>` for one
test page; `FINAL_OUT=<name>.pdf` overrides the output filename).

See also [`DESIGN_DECISIONS.md`](DESIGN_DECISIONS.md) (floor80 + chemotypes), [`../PIPELINE.md`](../PIPELINE.md).
