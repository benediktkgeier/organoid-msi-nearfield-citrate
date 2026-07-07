# Phase 07 — Metabolite identification + citrate QC

Match a published negative-mode reference m/z list to the feature grid using
accurate HMDB-formula theoretical masses, render per-metabolite ion images, and
run the citrate-resolution QC that documents the 191 co-isobar.

## Scripts (run order)
| Script | Role |
|---|---|
| `01_metabolite_match.R` | Published m/z list (ST1 liver + ST3 small intestine) → nearest feature within 20 ppm; **149 matched** + ppm/score/HMDB. |
| `02_metabolite_report.R` | One metabolite/page ion-image report (locked viridis/linear/p99.5/500 µm bar) → `figures/metabolites/metabolite_report.pdf`. |
| `02b_metabolite_report_single.R` | Single-section variant (default `AO_0h_sl6A_sec2a`); reuses the same match table. |
| `03_citrate_resolution.R` | Shows citrate [M−H]⁻ 191.0197 is not cleanly resolved from the 191 co-isobar (grid ion ~+8 ppm). |
| `04_citrate_window_images.R` | Citrate ion images over ±2.5/3/5 ppm windows (5 citrate-positive sections). |
| `05_citrate_isotopes_adducts.R` | Full citrate ion family (C12/C13, Na/Cl/acetate adducts) at ±10 ppm + spectral zoom-ins. |

## Inputs
- `cache/peaks_after_freq.rds` / `peaks_combined.rds`; published reference m/z list; RAW citrate windows via `lib_citrate.R`.

## Outputs
- `results/metabolites/`, `figures/metabolites/` (metabolite + citrate-resolution reports).

## Run (QC layer, off by default; enable with `./run_all.sh qc`)
```bash
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/07_metabolite_id/01_metabolite_match.R
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/07_metabolite_id/05_citrate_isotopes_adducts.R
```

## Notes / gotchas
- Use ACCURATE theoretical masses (HMDB formulas), not the supplement's 2-decimal m/z.
- Citrate is anchored/windowed via Phase 02 (`lib_citrate.R`), not the merged grid feature.
- Every tentative ID needs a PMID/DOI in the manuscript.

See also: [`../../docs/citrate_isotopes_adducts.md`](../../docs/citrate_isotopes_adducts.md), [`../../PIPELINE.md`](../../PIPELINE.md).
