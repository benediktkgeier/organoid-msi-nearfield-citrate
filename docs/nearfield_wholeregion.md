# Whole-region near-field citrate emission figure (LOCKED)

Whole-section companion to the per-organoid near-field figure
(`docs/apical_nearfield.md`). Instead of cropping each of the four views to one
organoid, every view is **projected across the entire measurement region** of a
section, and **every** organoid outline is coloured by its consensus apical class.
Produced for **all 20 datasets** as a single multi-page PDF.

Project root: `D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad/Analysis_R_Final`
R version: **R 4.4.2** (`C:/Program Files/R/R-4.4.2/bin/Rscript.exe`).

## Run
```
RS="C:/Program Files/R/R-4.4.2/bin/Rscript.exe"
"$RS" R/11_per_organoid_final/08_nearfield_wholeregion.R all              # 20 pages, per-section heatmap
"$RS" R/11_per_organoid_final/08_nearfield_wholeregion.R all globalheat   # 20 pages, GLOBAL heatmap scale
"$RS" R/11_per_organoid_final/08_nearfield_wholeregion.R AO_0h_sl6A_sec2a # single section (+ optional globalheat)
```
Or run `./run_wholeregion_figures.sh` (also rebuilds the perdataset v3 report).
No preprocessing/clustering/annotation is recomputed — reads finalized cache only.

## Outputs
| arg | output |
|---|---|
| `all` | `figures/annotation/apical_nearfield_emission_figure_all.pdf` (20 pages) |
| `all globalheat` | `figures/annotation/apical_nearfield_emission_figure_all_globalheat.pdf` (20 pages) |
| `<sid>` | `figures/annotation/apical_nearfield_wholeregion_<sid>[_globalheat].pdf` (1 page) |

The original per-organoid `apical_nearfield_emission_figure.pdf` (04) is **untouched**.

## Inputs
`cache/peaks_tissue_combined.rds`, `cache/zones_<sid>.rds` (`R/11_per_organoid_final/01_zones_curated.R`),
`cache/register/nd2final_<sid>.rds`, `figures/registration/crops/optical_<sid>.png`,
`results/annotation/apical_map_consensus.csv` (consensus class per `(sid,instance)`).

## Scripts
| script | role |
|---|---|
| `R/00_lib/lib_nearfield_viz.R` | **LOCKED** per-organoid 4-view toolkit (sourced, unchanged) |
| `R/00_lib/lib_nearfield_viz_whole.R` | **LOCKED** whole-region variants (`draw_*_whole`, `class_outlines`, `whole_region_vm`, `heatmap_global_clip`, `ylim_grid`) |
| `R/11_per_organoid_final/08_nearfield_wholeregion.R` | thin driver: dataset loop, global ranges, page layout, legend |

## LOCKED visual spec
2×2 landscape page per dataset, four views in this order:
1. **overlay** — citrate on native brightfield, **viridis**, constant alpha 0.60,
   global p99.5 clip `HI_CIT`; whole-region 50/100 µm signed-distance rings.
2. **gradmap** — citrate ion image on the MSI grid (viridis, `HI_CIT`) + 50/100 µm rings.
3. **heatmap** — `EBImage::gblur`-smoothed citrate on brightfield, **WEATHER rainbow**
   (LOW=blue/cold → HIGH=red/hot). Scale: **LOCKED pooled GLOBAL p99.9** (see below).
4. **vectors** — outward emission arrows from **every** organoid surface; **LENGTH ∝
   absolute** near-field citrate, **THICKNESS ∝ GLOBAL relative** score (p10..p90 of
   emission pooled across all rendered sections), 0.25×..2.0× base lwd, white + black halo.

**Outline colours (consensus apical class, from `lib_paths.R::APICAL_COLS`):**
apical-out magenta `#C2399A`, basolateral-out green `#2C9E4B`, mixed grey `#888888`,
unannotated `grey70`. Legend drawn along the bottom of each page; scale bars below each panel.

**Global comparability:** `HI_CIT` (views 1,2,4 ion images) and the arrow-thickness
range (view 4) are computed once across all rendered sections, so colour/length/thickness
mean the same thing on every page.

### Heatmap scaling (view 3) — *** LOCKED, UNIFORM ***
- **ALWAYS pooled GLOBAL p99.9.** `heatmap_global_clip()` (default `q = 0.999`) pools the
  smoothed whole-region values across all 20 datasets and takes one p99.9; every page uses that
  fixed clip and a small weather colorbar is drawn, so a colour = the same citrate level on every
  page. p99.9 is **dimmer** than the old p99 (which over-saturated the top 1% to red): genuine
  hotspots stay warm, the surround cools. The shared sampler `whole_region_vm()` guarantees the
  pre-pass and the render see identical values.
- The old **per-section auto-scale** default has been **removed** for uniformity; the
  `globalheat` CLI arg is still accepted but redundant. The same p99.9 percentile is used by the
  per-organoid hero heatmap (`04`, on its own scale) and the hi-res crops (`08b`).

## ORIENTATION (do not regress)
The native-BF panels (1,3) draw image row 1 at the top. The MSI→brightfield affine is
**vertically flipped** (`B[2,2] < 0`, ≈ −5.4 for every section), so the pure MSI-grid
panels (2,4) must put `max(ys)` at the top — `ylim_grid(sec)` returns `range(ys)` (not
`rev(range(ys))`) to match. Reverting to `rev(range(ys))` flips views 2 & 4 vertically
relative to 1 & 3. Arrows/contours/scale bars share the grid coords and flip consistently.

## Notes
- All sections render identically (single condition; no incubation-time comparison).
- A section with no consensus annotation (e.g. `AO_20h_sl4A_sec5a`, 0 annotated) still
  renders — all outlines grey.
- Glob can't traverse `D:` here — use `ls`/Grep with explicit paths.
