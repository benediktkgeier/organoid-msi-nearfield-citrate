# Phase 11 — Per-organoid final gradient + statistics

The headline analysis: per-organoid near-field citrate emission by apical class,
plus the locked publication figures. Apical classes come from the two-annotator
**consensus** map `results/annotation/apical_map_consensus.csv` (static committed
input) — no marked-PDF re-parse (and no python/reticulate) at runtime.

## Scripts (run order)
| Script | Role |
|---|---|
| `01_zones_curated.R` | Curated Voronoi zones on **refined** instances (all 20) → `cache/zones_<sid>.rds`. |
| `03_apical_report.R` | **Headline** report (default = consensus map) + `apical_gradient_per_organoid.csv`, `apical_citrate_dha_stats.csv`. |
| `04_nearfield_figure.R` | **LOCKED** four-view per-organoid near-field figure (`lib_nearfield_viz.R`). |
| `07_gradient_test_perorganoid.R` | Per-organoid Wilcoxon near-field test. |
| `08_nearfield_wholeregion.R` | **LOCKED** whole-region figure, all 20 datasets (`lib_nearfield_viz_whole.R`); `all globalheat` → global-scale variant. |
| `08b_view3_heatmap_hires.R` | Hi-res standalone render of the whole-region weather-heatmap view (view 3). |
| `09_citrate_gradient_perdataset_v3.R` | **LOCKED** per-dataset citrate gradient — all 20 datasets. |
| `09b_citrate_ion_hires.R` | Hi-res (600 dpi) standalone citrate ion-image PNGs per dataset. |
| `09_prism_export.R` | GraphPad Prism export (xlsx/pzfx) of the apical citrate/DHA stats. |
| `10_citrate_gradient_report_final.R` | **FINAL combined report**: per-dataset 6-panel page (MSI / native BF / SSC mask / matched IF / MSI↔BF overlay / weather gradient map) + capped ≤100 µm apical-class overview page → `figures/gradient/citrate_gradient_report_final.pdf`. Non-locked driver; sources the locked libs and copies their small helpers. |
| `10_pzfx_export.R` | Prism `.pzfx` export companion to `09_prism_export`. |
| `11_export_gradient_profile.R` | Tabular export of the per-organoid outward gradient profiles. |
| `12_joy_tables.R` | Curated "Joy_Tables" subset export for the collaborator. |
| `13_citrate_gradient_report_3class.R` | **Apical-class variant** of `10_` (non-destructive sibling). Overview minis treat "mixed" explicitly: default (`withmixed`) gives mixed its own bold grey trend line (3-group); `nomixed` drops mixed (2-group headline). → `figures/gradient/citrate_gradient_report_{withmixed,nomixed}.pdf`. |

## Inputs
- `cache/instances_final_<sid>.rds`; citrate cache (`lib_citrate.R`); **`apical_map_consensus.csv`**; native crops (Phase 05).

## Outputs
- `figures/annotation/apical_citrate_dha_report.pdf` (headline; page 1 = Summary Statement) +
  `results/annotation/apical_citrate_dha_stats.csv`; whole-region + per-dataset figures.

## Run
```bash
./run_apical_nearfield_pipeline.sh
./run_wholeregion_figures.sh
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/11_per_organoid_final/08_nearfield_wholeregion.R all globalheat
# final combined per-dataset report:
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/11_per_organoid_final/10_citrate_gradient_report_final.R all    # 20 dataset pages + overview
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/11_per_organoid_final/10_citrate_gradient_report_final.R test   # quick TEST report (overview + 1 dataset) -> ..._TEST_<sid>.pdf
```

## Notes / gotchas
- **Number order:** the two `09_` (`_v3` figure / `_prism_export`) and two `10_` (`_report_final` / `_pzfx_export`)
  scripts are independent — a figure/report and its data export share a number; run either directly.
- `10_citrate_gradient_report_final.R`, `13_citrate_gradient_report_3class.R`, and
  `04_ssc_ontissue/05_ssc_report.R` are non-locked cache-only reports; all accept `all` (default),
  `test` (a quick one-dataset TEST report → `..._TEST_<sid>.pdf`), or an explicit `<sid>`, and honor a
  `FINAL_OUT=<name>.pdf` env override. All are wired into `run_all.sh`.
- **The "mixed" apical class across the final-report family** (three explicit outputs, so nothing is
  ambiguous): `10_ → citrate_gradient_report_final.pdf` shows mixed as **grey context only** (no bold
  trend); `13_ (default) → citrate_gradient_report_withmixed.pdf` gives mixed **its own bold grey trend
  line** (3-group comparison); `13_ nomixed → citrate_gradient_report_nomixed.pdf` **drops mixed
  entirely** (2-group headline). `run_all.sh`/`regen_reports.sh` emit all three.
- Headline claim: `grad_near50` (≤50 µm) citrate apical-out > basolateral-out (p ≈ 6e-5, Cliff's δ),
  corroborated by citrate/DHA; outward slope `grad_rho_out` is n.s. — this is a near-field **level**, not steepness.
- Dropped in this fork (NOT present): `02_apical_citrate_dha` (single-annotator, superseded by consensus),
  `03b_apical_report_50um_nomixed`, `05_apical_toppick_report`, `06_citrate_gradient_perdataset`.

See also: [`../../docs/apical_nearfield.md`](../../docs/apical_nearfield.md), [`../../docs/nearfield_wholeregion.md`](../../docs/nearfield_wholeregion.md), [`../../docs/citrate_gradient_perdataset.md`](../../docs/citrate_gradient_perdataset.md), [`../../PIPELINE.md`](../../PIPELINE.md).
