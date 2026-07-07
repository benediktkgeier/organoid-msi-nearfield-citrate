# Phase 09 — Organoid separation refinement (ROI curation)

Turn the raw connected-component instances into curated, one-ROI-per-organoid
final instances using manual green-ink cuts and island delete/merge actions
(both cached). `00_run_roi_curation.R` chains the apply steps.

## Scripts (run order)
| Script | Role |
|---|---|
| `00_run_roi_curation.R` | Driver: chains `01` → `03` → `04`. |
| `01_organoid_split_apply.R` | Split merged organoids (**MANUAL** green-ink cuts, cached). |
| `02_island_cleanup_canvas.R` | Build the island-cleanup annotation canvas. |
| `02b_island_centroid_sidecar.R` | Centroid sidecar for parsing cleanup actions. |
| `03_island_cleanup_apply.R` | Apply delete/merge (**MANUAL** actions, cached). |
| `04_finalize_instances.R` | One connected ROI = one id → `cache/instances_final_<sid>.rds`. |

## Inputs
- `cache/instances_<sid>.rds` (Phase 08); manual split/cleanup markup (cached).

## Outputs
- `cache/instances_final_<sid>.rds` (curated ROIs); annotation canvas + centroid sidecar.

## Run
```bash
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/09_organoid_refinement/00_run_roi_curation.R
```

## Notes / gotchas
- Split cuts and island cleanup are **manual gates**; their marked artifacts are cached and reused in this fork.
- Apical classes are annotated later (Phase 10) on the `organoid_island_cleanup.pdf` canvas built here.

See also: [`../../docs/roi_curation.md`](../../docs/roi_curation.md), [`../../PIPELINE.md`](../../PIPELINE.md).
