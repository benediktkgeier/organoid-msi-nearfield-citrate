# Phase 08 — Organoid segmentation + pooled gradient survey

Segment organoids from the on-tissue mask, build signed-distance rings into the
surrounding CMC gel, and run a **pooled, single-condition descriptive** gradient
survey across ALL sections. Supporting context only — the headline result is the
Phase 11 apical near-field analysis.

## Scripts (run order)
| Script | Role |
|---|---|
| `01_segment_organoids.R` | Connected components of `is_tissue` → `cache/instances_<sid>.rds`. |
| `02_buffer_rings.R` | Signed-distance rings (outward `BUF_LEVELS_UM`, inward 10 µm) → `cache/zones_<sid>.rds`. |
| `03_gradient_stats.R` | Per-ion pooled ρ_out / ρ_in across all sections → `results/gradient/`. |
| `04_report_pdf.R` | `figures/gradient/gradient_report.pdf` (pooled ρ ranking + profiles). |

## Inputs
- `cache/peaks_tissue_combined.rds`; citrate cache (`lib_citrate.R`); `gradient_config.R`.

## Outputs
- `cache/instances_<sid>.rds`, `cache/zones_<sid>.rds`, `results/gradient/`, `figures/gradient/gradient_report.pdf`.

## Run
```bash
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/08_organoid_gradient_survey/03_gradient_stats.R
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/08_organoid_gradient_survey/04_report_pdf.R
```

## Notes / gotchas
- **Single-condition study: POOLED descriptive** — no group split, no Δρ.
- The former split-by-group scripts (`05_citrate_20datasets`, `06_citrate_dha_compare`,
  `07_citrate_dha_20datasets`, `08_gradient_20datasets`) were **dropped** in this fork — not present.
- Segmentation/rings caches reused READ-ONLY where available; refined ROIs come from Phase 09.

See also: [`../../PIPELINE.md`](../../PIPELINE.md).
