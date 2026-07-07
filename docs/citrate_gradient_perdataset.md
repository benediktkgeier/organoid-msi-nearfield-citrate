# Per-dataset citrate/DHA gradient report — v3, all 20 datasets (LOCKED)

Per-dataset combined report comparing **anchored citrate [M-H]- 191.0198 ±7 ppm** (signal;
`citrate_onto_pd()` from `lib_citrate.R`, not the old grid feature 191.0217) vs
**DHA C22:6 [M-H]- 327.2330** (confined-lipid control): overlays, single-ion images,
brightfield overlays with Voronoi zones, and per-organoid outward-gradient curves.

`v3` is the **all-20-datasets** report. (It was derived from an earlier single-dataset TEST
report — `06_citrate_gradient_perdataset` — which was **dropped in this fork**; `v3` is the
only per-dataset gradient script carried forward and defaults to all 20 sections.)

Project root: `D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad/Analysis_R_Final`
R version: **R 4.4.2** (`C:/Program Files/R/R-4.4.2/bin/Rscript.exe`).

## Run
```
RS="C:/Program Files/R/R-4.4.2/bin/Rscript.exe"
"$RS" R/11_per_organoid_final/09_citrate_gradient_perdataset_v3.R         # all 20 (default)
"$RS" R/11_per_organoid_final/09_citrate_gradient_perdataset_v3.R all     # all 20
"$RS" R/11_per_organoid_final/09_citrate_gradient_perdataset_v3.R <sid>   # one named section
```
Or run `./run_wholeregion_figures.sh` (rebuilds this + the whole-region figures).

## Outputs
| arg | output |
|---|---|
| (none) / `all` | `figures/gradient/citrate_gradient_perdataset_v3.pdf` (1 summary + 20 pages = **21**) |
| `<sid>` | `figures/gradient/citrate_gradient_perdataset_v3_TEST_<sid>.pdf` (summary + 1 page) |

## Inputs
`cache/peaks_tissue_combined.rds`, `cache/zones_<sid>.rds` (`R/11_per_organoid_final/01_zones_curated.R`),
`cache/register/nd2final_<sid>.rds`, `figures/registration/crops/optical_<sid>.png`.

## "Uses the consensus segmentation annotations"
No class recolour is applied. The consensus segmentation is **already baked into the
inputs**: `R/11_per_organoid_final/01_zones_curated.R` rebuilds every `cache/zones_<sid>.rds`
from the consensus-curated `instances_{final,clean,split}_<sid>.rds`, whose instance IDs
match the apical annotations. Re-rendering from the current zones therefore reflects the
consensus segmentation. Outlines stay white-dotted; per-organoid curves keep the v2 rainbow.

## Page layout (frozen)
- **Page 1 — SUMMARY:** all 20 datasets tiled, each as the citrate(yellow)+DHA(purple)
  overlay paired with its outward-decay curve (per-organoid faint + bold section median;
  x = 10–500 µm log).
- **Pages 2–21 — per dataset (2×3):** overlay (citrate=yellow, DHA=purple) · citrate ion
  image + white-dotted organoid outline · DHA ion image + outline · native brightfield ·
  brightfield + citrate overlay (constant opacity) + white zone contours · per-organoid
  outward gradient (citrate lines coloured per organoid, DHA grey-dashed control). Citrate
  uses the anchored [M−H]− definition (191.0198 ±7 ppm) from `lib_citrate.R`.

Global p99.5 clips for citrate and DHA are computed once across all 20 datasets.

## Scripts
| script | role |
|---|---|
| `R/11_per_organoid_final/09_citrate_gradient_perdataset_v3.R` | **LOCKED** v3 (all-20 default); the only per-dataset gradient script in the fork |

## Notes
- All sections render identically (single condition; no incubation-time comparison).
- Glob can't traverse `D:` here — use `ls`/Grep with explicit paths.
