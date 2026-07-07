# Phase 03 — Coarse MSI→brightfield registration

Establish the coarse chain that maps MSI pixels onto the optical images:
MSI pixel → stage/flexImage → slide `_small.jpg` (teach points) → `.nd2`.

## Scripts (run order)
| Script | Role |
|---|---|
| `01_teach_msi_to_jpg.R` | MSI→slide-JPG teach transform (from Bruker `.mis` Area + teach points). |
| `02_jpg_to_nd2.R` | JPG→`.nd2` block-mean thumbs / coarse transform (scale + flip). |
| `03_slide_overview.R` | Slide-level QC overview. |

## Inputs
- Bruker `.mis`, slide `_small.jpg`, `.nd2` optical images; `inventory.csv`.

## Outputs
- Coarse MSI→`.nd2` transforms cached in `cache/register/`; slide overview QC figures.

## Run
```bash
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/03_coarse_registration/01_teach_msi_to_jpg.R
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/03_coarse_registration/02_jpg_to_nd2.R
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/03_coarse_registration/03_slide_overview.R
```

## Notes / gotchas
- Transforms are refined in Phase 05; keep the two phases' caches (`cache/register/`) coherent.
- In this fork the register cache is reused READ-ONLY via `cache_in("register")`.

See also: [`../../docs/registration.md`](../../docs/registration.md), [`../../PIPELINE.md`](../../PIPELINE.md).
