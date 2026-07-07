# Phase 04 — SSC on-tissue delineation

Build the curated ion set and the on-tissue mask (with chemotype labels) used by
segmentation and all downstream analysis. On-tissue rule = "floor80" + 4 chemotypes.

## Scripts (run order)
| Script | In → Out |
|---|---|
| `01_export_ontissue_candidates.R` | `peaks_combined.rds` → candidate ion sheet (**MANUAL** on-tissue selection). |
| `02_import_ontissue_ions.R` | Manual selection (`results/peakme_annotations/`) → `cache/peaks_combined_annot.rds`. |
| `03_build_curated_set.R` | Annot MSE + published m/z list → **`cache/peaks_curated.rds`** (348 ions). |
| `04_ssc_tissue_mask.R` | `peaks_curated.rds` → **`cache/peaks_tissue_combined.rds`** (`is_tissue`, `ssc_k4_sec`, `ssc_k10_sec`, `chemotype`, `tissueness_px`) + per-section `figures/ssc/ssc_mask_<sid>.pdf`. |
| `05_ssc_report.R` | `peaks_tissue_combined.rds` + `zones_<sid>.rds` → **`figures/ssc/ssc_clustering_segmentation_report.pdf`** (visual SSC + segmentation report; reads cached columns, no recompute). |

## Inputs
- `cache/peaks_combined.rds`; manual on-tissue ion selection (`results/peakme_annotations/`); published m/z list.

## Outputs
- `cache/peaks_curated.rds` (348 ions), `cache/peaks_tissue_combined.rds` (headline mask + SSC columns),
  per-section `figures/ssc/ssc_mask_<sid>.pdf`, and `figures/ssc/ssc_clustering_segmentation_report.pdf`.

## Run
```bash
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/04_ssc_ontissue/03_build_curated_set.R
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/04_ssc_ontissue/04_ssc_tissue_mask.R
# visual SSC + segmentation report (cache-only):
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/04_ssc_ontissue/05_ssc_report.R all    # full 20-dataset report
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/04_ssc_ontissue/05_ssc_report.R test   # quick TEST report (front + 1 dataset) -> ..._TEST_<sid>.pdf
```

## Notes / gotchas
- Steps 01/02 are a **manual gate** (on-tissue ion selection); its marked output is reused/cached in this fork.
- `results/peakme_annotations/` is the legacy storage name for the manual on-tissue selection.
- `04` persists `tissueness_px / is_tissue / ssc_k4_sec / ssc_k10_sec / chemotype` into the mask cache; downstream only `is_tissue` is consumed, but `05_ssc_report.R` reads the cluster/chemotype columns directly to render the visual report (no recompute).
- Curated set = 149 published knowns + 237 on-tissue unknowns (348 total).

See also: [`../../docs/DESIGN_DECISIONS.md`](../../docs/DESIGN_DECISIONS.md) (floor80 + chemotypes), [`../../PIPELINE.md`](../../PIPELINE.md).
