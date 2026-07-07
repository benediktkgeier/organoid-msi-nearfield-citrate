# Phase 02 — Citrate standard GATE

Critical early gate. Validates citrate against an authentic Sodium Citrate
Tribasic Dihydrate dilution series (spotted in the same CMC matrix, same slide
sl7A) and **LOCKS** the citrate definition for the whole pipeline:
`CITRATE_ANCHOR_MZ = 191.01976` (+0.16 ppm vs theo), `CITRATE_WIN_PPM = 7`
(both in `lib_paths.R`). Only proceed downstream if the spectra make sense.

## Scripts (run order)
| Script | Role |
|---|---|
| `00_config.R` | Spot↔conc table, masses; sources `R/00_lib/lib_citrate.R` (reader + helpers). |
| `01_standard_spectra_mz.R` | Pure-standard centroid +0.16 ppm, FWHM ~7 ppm; tissue 191 = +14.6 ppm / 3.2× broader = co-isobar. |
| `02_calibration_curve.R` | log–log slope 1.00, R² 0.987 (1–100 mM); LOD ~0.37 mM; both decade pairs were mislabeled → mapped to true conc. |
| `03_id_fingerprint.R` | Standard ¹³C 7.3% ≈ 6.8% theo + clean adduct profile; tissue inflated = co-isobar. |
| `04_standard_anchored_citrate.R` | GATE self-check (measured anchor ≈ `CITRATE_ANCHOR_MZ`); 5/7/10 ppm sweep → **±7 ppm**. |
| `05_build_citrate_cache.R` | Writes `cache/citrate_anchored_<sid>.rds` (x,y,cit_raw,tic) for all 20 samples. |

## Inputs
- Standard imzML `MSI/06102026_AO_0h_sl7A/imzml` (0/blank, 10 µM … 100 mM); RAW centroid tissue imzML.

## Outputs
- Locked citrate anchor/window (via `lib_paths.R`), `cache/citrate_anchored_<sid>.rds`,
  reports in `results/citrate_standard/`, `figures/citrate_standard/`.

## Run (QC layer, off by default in `run_all.sh`)
```bash
./run_all.sh qc      # runs 01–04 (+ metabolite QC)
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/02_CitrateStandard/05_build_citrate_cache.R
```

## Notes / gotchas
- ALL downstream citrate (phases 05/07/08/11 + `lib_nearfield_viz.R`) calls
  `citrate_onto_pd()` from `lib_citrate.R`, NOT the merged grid feature.
- **Centroided data caveat**: a narrow window SELECTS citrate-dominant pixels; it does not integrate area.
- Comparison design: matrix-matched standard (same CMC/slide/spray); on-sample co-acquisition impossible.

See also: [`../../docs/citrate_isotopes_adducts.md`](../../docs/citrate_isotopes_adducts.md), [`../../PIPELINE.md`](../../PIPELINE.md).
