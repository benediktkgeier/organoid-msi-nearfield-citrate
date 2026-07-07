# Phase 01 — Spectral preprocessing

Build the self-contained, feature-aligned peak MSE that the rest of the pipeline
consumes. Cardinal v3 idiom: `readMSIData` → `c()` combine → `normalize(tic)` →
`process()` → manual single-linkage reference grid → `convertMSImagingArrays2Experiment`
→ `peakAlign` → `summarizeFeatures`.

## Scripts (run order)
| Script | Role |
|---|---|
| `01_preprocess.R` | imzML (`cache/imzml`) + `inventory.csv` → aligned MSE `cache/peaks_after_freq.rds`. |
| `02_realize_in_memory.R` | Realize peaks MSE in-memory → self-contained (breaks lazy `.ibd` dependency). |
| `03_refilter_freq_floor.R` | Frequency-floor refilter → **`cache/peaks_combined.rds`** (8,772 × 90,760). |
| `qc_render.R` | QC → `figures/preprocess/qc_combined.pdf` + feature-filter diagnostic. |
| `qc_top_ion_images.R` | Top-ion images → `figures/preprocess/top_ion_images.pdf`. |

## Inputs
- imzML/`.ibd` under `cache/imzml`, `inventory.csv` (sample registry).

## Outputs
- `cache/peaks_after_freq.rds`, `cache/peaks_combined.rds` (headline, self-contained MSE),
  QC PDFs under `figures/preprocess/`.

## Run
```bash
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/01_preprocess/01_preprocess.R
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/01_preprocess/02_realize_in_memory.R
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/01_preprocess/03_refilter_freq_floor.R
```

## Notes / gotchas
- R 4.4.2 / Cardinal v3 only; **never** `setCardinalBPPARAM()`.
- Cardinal RDS is NOT self-contained until step 02 realizes it in memory (lazy `.ibd` reads otherwise).
- In this fork `peaks_combined.rds` is reused READ-ONLY from upstream via `cache_in()`.

See also: [`../../PIPELINE.md`](../../PIPELINE.md).
