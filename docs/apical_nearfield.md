# Apical orientation → near-field citrate emission (LOCKED)

End-to-end, reproducible pipeline that asks: **do apical-out organoids emit more
citrate into the surrounding gel than basolateral-out organoids?** It scores each
organoid's apical orientation, normalizes per-organoid intensities several ways,
builds outward Voronoi zones, and renders a **locked four-view publication figure**.

Project root: `D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad/Analysis_R_Final`
R version: **R 4.4.2** (`C:/Program Files/R/R-4.4.2/bin/Rscript.exe`).

## Headline result (descriptive; organoid = unit)
With **absolute anchored citrate [M-H]- 191.0198 ±7 ppm** (no surface normalization),
apical-OUT organoids emit more citrate into the immediate gel than basolateral-out. The
locked headline number is the near-field ≤50 µm level (`grad_near50`, apical-out >
basolateral-out, **p ≈ 6.07e-5**), read live from `apical_citrate_dha_stats.csv` and shown
as the Summary Statement on page 1 of `apical_citrate_dha_report.pdf`. The consensus map has
n=84 organoids (30 apical-out / 23 basolateral-out / 31 mixed).

Caveats locked in: citrate is measured as **[M-H]- only** (the earlier
"adduct-switching / total-citrate" idea was an artifact and has been removed); the
*surface-normalized* outward gradient is a confounded readout (basolateral-out looks
"flat" only because its surface citrate ≈ gel background), so the **absolute
near-field** numbers are the headline. Stats are descriptive (organoids within a
section are pseudo-replicates; organoid = unit). This is a **single-condition** study —
all sections incubated in CMC for ≥5 min before freezing and treated identically (no
incubation-time comparison); the contrast is apical polarity only.

## Run order (each step depends on the previous data)
```
RS="C:/Program Files/R/R-4.4.2/bin/Rscript.exe"
# apical classes come from the CONSENSUS map (static committed input); zones/instances reused.
"$RS" R/11_per_organoid_final/03_apical_report.R      # headline report + apical_gradient_per_organoid.csv (consensus)
"$RS" R/11_per_organoid_final/04_nearfield_figure.R   # LOCKED 4-view publication figure + panels
```
Or just run `./run_apical_nearfield_pipeline.sh`.

### Two-class (`nomixed`) variant — drop the ambiguous `mixed` class
Both drivers take an optional flag that excludes the `mixed` annotator-conflict group,
leaving only the two clean polarity classes (basolateral-out n=23 / apical-out n=30). Outputs
are **suffixed `_nomixed`** and written alongside the 3-class originals (non-destructive).
```
"$RS" R/11_per_organoid_final/03_apical_report.R    "" "" nomixed   # args: [1]map [2]suffix [3]drop-mixed flag
"$RS" R/11_per_organoid_final/04_nearfield_figure.R nomixed         # arg:  [1]drop-mixed flag
```
Produces `apical_citrate_dha_report_nomixed.pdf`, `apical_nearfield_emission_figure_nomixed.pdf`,
`nearfield_panels_nomixed/`, and the matching `*_nomixed.csv` stats/gradient tables. Implementation:
a `DROP_MIXED` flag filters `CLASSES`/`PAIRS` and the dot-box group count derives from
`length(CLASSES)`; `04` reads the `_nomixed` gradient CSV that `03` writes. The headline
`grad_near50` (basolateral-out vs apical-out) is a two-class comparison, so **p ≈ 6.07e-5 is
unchanged** — dropping `mixed` only removes that group from the panels and the `*_vs_mixed`
pairwise Wilcoxon rows; in the per-section overlays mixed organoids become faint-grey unannotated.

## Where the apical annotations live (important)
The apical class per organoid is the **two-annotator CONSENSUS** map
**`results/annotation/apical_map_consensus.csv`** (`sid, instance, apical_class` ∈
{`apical_out`, `basolateral_out`, `mixed`}; n=84 = 30/23/31). This is a **static committed input**
in the fork — the consensus was finalized upstream (Phase 10) from two annotators' comments on the
island-cleanup canvas. The fork therefore does **not** re-parse any marked PDF (no python/reticulate
at runtime). `03_apical_report.R` reads this map by default; pass a different map + suffix as CLI
args to analyze a subset.

## Scripts
| script | role | key output |
|---|---|---|
| `R/11_per_organoid_final/01_zones_curated.R` | signed-distance Voronoi zones from **curated** instances, all 20 sections (surface recomputed per instance) | `cache/zones_<sid>.rds`, `results/gradient/zones_curated_summary.csv` |
| `R/11_per_organoid_final/03_apical_report.R` | headline report (default = consensus map): absolute-TIC / metabolite-TIC (companion, right after the absolute page) / within-section / citrate-DHA-ratio / outward-gradient / metrics / **near-field** / a **paired within-class page** comparing near-field citrate at 0–50 vs 0–100 µm inside each class / per-section overlays; page 1 = Summary Statement | `figures/annotation/apical_citrate_dha_report.pdf`, `apical_gradient_per_organoid.csv`, `apical_citrate_dha_stats.csv` |
| `R/00_lib/lib_nearfield_viz.R` | **LOCKED** four-view visualization toolkit (sourced by reports); colours from `lib_paths.R::APICAL_COLS` (apical-out = magenta `#C2399A`, basolateral-out = green `#2C9E4B`, mixed = grey `#888888`) | — |
| `R/11_per_organoid_final/04_nearfield_figure.R` | thin driver → publication figure + panels | `figures/annotation/apical_nearfield_emission_figure.pdf`, `figures/annotation/nearfield_panels/` |

## LOCKED visualization spec (`R/00_lib/lib_nearfield_viz.R`)
Every view is cropped to the organoid bbox + **120 µm** margin and **masked to the
organoid interior + its Voronoi catchment** (the exact region near-field metrics are
measured over; neighbours excluded). Distance rings/outlines are hero-only.
1. **`draw_overlay`** — ion MSI on native brightfield, **viridis**, constant alpha
   0.60; solid white organoid outline + **dotted 50/100 µm rings**.
2. **`draw_gradmap`** — cropped ion image (viridis, global p99.5 clip = `HI_CIT`)
   + solid 0 µm / dotted 50,100 µm signed-distance rings.
3. **`draw_heatmap`** — interpolated **weather-map heatmap** on brightfield:
   `EBImage::gblur`-smoothed ion, **WEATHER rainbow LOW=blue(cold) → HIGH=red(hot)**,
   **LOCKED p99.9 clip** (per-organoid: its own p99.9; whole-region: pooled global p99.9) —
   dimmer than the old p99, so genuine hotspots stay warm and the surround cools.
4. **`draw_vectors`** — outward arrows from surface points (normal = ∇ signed
   distance). **LENGTH ∝ absolute** near-field ion (sampled 20–50 µm out);
   **THICKNESS ∝ relative score scaled GLOBALLY** across exemplars (p10..p90 of
   emission) over **0.25×..2.0×** base lwd (= −75%..+100%). White arrow + black halo.

Conventions: viridis ion images; scale bar **below** the image with centred label; polarity
outline/legend colours come from `lib_paths.R::APICAL_COLS` (apical-out = magenta `#C2399A`,
basolateral-out = green `#2C9E4B`, mixed = grey `#888888`).
Exemplars = top-3 apical-out / bottom-3 basolateral-out by `near50`; hero = top apical-out.

## Reuse the locked views in a new report
```r
ROOT <- "D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad/Analysis_R_Final"
source(file.path(ROOT,"R/00_lib/gradient_config.R")); source(file.path(ROOT,"R/00_lib/lib_register.R"))
source(file.path(ROOT,"R/00_lib/lib_nearfield_viz.R"))
nfviz_load_ion(CIT_MZ)                       # sets val_cit, HI_CIT, pd (any m/z works)
nfviz_arrow_range(my_exemplars[,c("sid","instance")])   # sets GLO, GHI for arrow thickness
sec <- prep_sec("AO_20h_sl4A_sec3b")
draw_heatmap(sec, 3)                          # or draw_overlay / draw_gradmap / draw_vectors
```
Requires `cache/zones_<sid>.rds` (`R/11_per_organoid_final/01_zones_curated.R`) and the registration assets
(`cache/register/nd2final_<sid>.rds`, `figures/register/crops/optical_<sid>.png`).

## Reproducibility notes
- Tooling caveat: the Glob tool can't traverse `D:` here — use `ls`/Grep with explicit
  paths to check files.
- Cleaned data CSVs from the adduct purge keep `*.bak` backups
  (`curated_feature_table.csv`, `metabolite_match_table.csv`, `gradient_stats*.csv`).
- `peaks_tissue_combined.rds` is the curated **on-tissue** MSE (348 features incl. gel
  pixels with `is_tissue`); citrate adducts were removed from the metabolite-TIC pool.
