# Phase 06 — IF (B-section) → brightfield/MSI registration

Register serial immunofluorescence B-sections onto the brightfield / MSI grid via
manual landmarks (cached), and emit the paired IF↔MSI dataset report.
Chain: hi-res `.nd2` (0.362 µm/px, 4ch) → overview `.nd2` (1.833 µm/px, DAPI) →
MSI BF `.nd2` → inverse MSI affine → MSI grid.

## Scripts (run order)
| Script | Role |
|---|---|
| `01_if_thumbs.R` | IF thumbnails; `01b_build_hr4raw.R` / `01c_build_hrdapihq.R` = hi-res raster builds. |
| `02_if_overview_to_bf.R` … `09_if_to_bf_register.R` | Overview→BF / overview→MSI, native overlays, section crops, subregions, BF-guided register. |
| `10_landmark_sheets.R` | Emit manual-landmark sheets. |
| `11_landmark_fit.R` | Fit affine from cached manual landmarks. |
| `12_dataset_pairs_hq.R` | **Headline** paired IF↔MSI report → `figures/if_registration/dataset_pairs_hq.pdf`. |

## Inputs
- IF `.nd2` (`LM_Bsections`); MSI↔BF transforms `cache/register/` (read-only reuse); manual landmarks (cached).

## Outputs
- `cache/register_if/`, `figures/if_registration/dataset_pairs_hq.pdf`, `results/if_registration/`.

## Run
```bash
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/06_if_registration/12_dataset_pairs_hq.R
```

## Notes / gotchas
- Channel order: ch1=Cy5 / ch2=ZO-1 / ch3=b-catenin / ch4=DAPI.
- **DAPI = CYAN.** The pairs report (`12_`) shows **ZO-1 red + DAPI cyan only** (b-catenin dropped).
- Landmarks are a **manual gate**; cached and reused in this fork.

See also: [`../../docs/if_registration.md`](../../docs/if_registration.md), [`../../PIPELINE.md`](../../PIPELINE.md).
