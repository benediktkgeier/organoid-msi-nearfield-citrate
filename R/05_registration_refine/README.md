# Phase 05 — Registration refinement

Refine the coarse Phase 03 registration into publication-grade native-resolution
MSI↔`.nd2` affines and native brightfield crops.

## Scripts (run order)
| Script | Role |
|---|---|
| `01_refine_jpg.R` | Refine the MSI→small-JPG transform. |
| `02_jpg_to_nd2_offset.R` | JPG→`.nd2` refinement (V-flip, offset) → reusable MSI→`.nd2` affines `cache/register/nd2final_<sid>.rds`. |
| `03_native_crops.R` | Native-resolution brightfield crops per section. |
| `04_overlay_report.R` | Overlay QC → `figures/registration/registration_native.pdf`. |

## Inputs
- Coarse transforms `cache/register/` (Phase 03); `.nd2` optical images.

## Outputs
- Reusable affines `cache/register/nd2final_<sid>.rds`, native crops,
  `figures/registration/registration_native.pdf`.

## Run
```bash
for s in 01_refine_jpg 02_jpg_to_nd2_offset 03_native_crops 04_overlay_report; do
  "/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/05_registration_refine/$s.R
done
```

## Notes / gotchas
- All microscopy is rendered at NATIVE resolution (never downsampled to the MSI grid).
- MSI-on-BF overlays use homogeneous constant opacity (alpha ~0.5), not intensity-thresholded transparency.
- Register cache reused READ-ONLY in this fork via `cache_in()`.

See also: [`../../docs/registration.md`](../../docs/registration.md), [`../../PIPELINE.md`](../../PIPELINE.md).
