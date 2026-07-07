# IF (immunofluorescence) ↔ brightfield / MSI registration

Registers immunofluorescence imaging of **consecutive (serial) B-sections** onto the
**MSI A-slide brightfield/MSI frame**, so IF markers can be related to the MSI metabolite
data per organoid. This is a **manual-landmark, organoid-level** registration — automatic
shape matching failed (serial sections + IF captures more tissue than MSI measures), so
the working approach uses user-marked organoid centers.

> **Status: complete.** Final deliverable: `figures/if_registration/dataset_pairs_hq.pdf`
> (produced by `R/06_if_registration/12_dataset_pairs_hq.R`).
>
> **Path note:** this doc predates the phase-folder restructure and uses the old flat script
> numbers (`R/90…R/101`) and the old figure dir `register_if/`. In the fork the scripts live
> under `R/06_if_registration/` (01…12 — see `PIPELINE.md` Phase 06 and the legacy→new table)
> and the output dir is `figures/if_registration/`. Configs/helpers are `R/00_lib/if_config.R`
> and `R/00_lib/lib_register_if.R`. The channel/pixel/LUT facts below are unchanged.

---

## 1. Data

| Item | Location | Notes |
|---|---|---|
| IF whole-slide overview | `LM_Bsections/06172026_AO_{0h_sl6b,20h_sl4b}_over.nd2` | **DAPI only**, 6× (4×·1.5), **1.833 µm/px**, ~29k×13k |
| IF high-res sections | `LM_Bsections/06172026_AO_0h_sl6b_sec{1..6}.nd2`, `..._20h_sl4_sec{1..6}.nd2` | **4 channels**, 30× (20×·1.5), **0.362 µm/px**, 13080×13071 |
| MSI A-slide brightfield | `06102026_AO_{0h_sl6A,20h_sl4A}.nd2` | 6× (4×·1.5), **1.833 µm/px** — the registration target frame |
| MSI data | `cache/peaks_tissue_combined.rds` (`TISSUE_MSE`) | `sample_id`, `is_tissue`, ion spectra |
| MSI↔BF transforms | `cache/register/nd2final_<sid>.rds` | `B_msi_nd2` (MSI grid → BF native); from the MSI-registration arc (`README_registration.md`) |

### Channel order (CRITICAL — verified from NIS metadata, emission wavelength)
The hi-res `.nd2` channel order is **reversed** from the naive assumption:

| ch | marker | colour | emission |
|---|---|---|---|
| **1** | Cy5 / **F-actin** | gray | 670 nm |
| **2** | mCherry / **ZO-1** | red | 642 nm |
| **3** | FITC / **β-catenin** | green | 524 nm |
| **4** | DAPI / **DNA** | cyan | 405 nm |

Codified in `R/if_config.R` as `IF_CH = c(Factin=1, ZO1=2, bcat=3, DAPI=4)` (marker-keyed,
so a channel/marker mix-up cannot recur). **DAPI is ch4, ZO-1 is ch2.**

### Pixel size & the 1.5× magnification
We use `dCalibration` from each `.nd2`, which **already includes** the 1.5× changer
(`dObjCalibration1to1 / dCalibration = 1.5000` exactly for every file). So:
`HR_UMPX = 0.36183` (30×, with 1.5×), `BF_UMPX = OV_UMPX = 1.83333` (6×). **Do not** apply
1.5× by hand and **do not** use `dObjCalibration1to1` (the objective-only value).

### Native NIS display LUTs (per marker, per slide) — `IF_LUT` in `R/if_config.R`
Linear, γ=1.0; display = `clamp((raw-lo)/(hi-lo),0,1)`.

| | DAPI | ZO-1 | β-catenin | F-actin |
|---|---|---|---|---|
| **sl6b** lo–hi | 148–4095 | 413–2663 | 420–3324 | 132–3184 |
| **sl4b** lo–hi | 210–2546 | 459–2904 | 335–2771 | 140–3083 |

---

## 2. Pipeline (run order) — all **R 4.4.2**, `JAVA_HOME=jre1.8.0_491`

Config/helpers: `R/if_config.R` (paths, `IF_CH`, `IF_LUT`, calib, `IF_SLIDES`),
`R/lib_register_if.R` (nd2 block-mean/min, FFT-NCC, masks, LUT helpers).

| Step | Script | Output |
|---|---|---|
| Caches: DAPI thumbs | `R/90_if_thumbs.R` | `cache/register_if/{ovthumb,hrthumb}_*.rds` |
| Caches: 4-ch raw (F=5) | `R/90b_build_hr4raw.R [sl6b\|sl4b]` | `hr4raw_<sid>.rds` (used by reports & fit) |
| Caches: HQ DAPI (F=2) | `R/90c_build_hrdapihq.R` | `hrdapihq_<sid>.rds` |
| **Brightfield rectangles** | `R/94_if_section_crops.R` (sheets) + `R/95_if_brightfield_guided.R` (MSI-guided) → user draws 6 boxes/slide → `extract_bf_rectangles.py` | `results/register_if/bf_rectangles.csv` |
| QC: ROI pairs (full field) | `R/96_if_roi_pairs.R` | `roi_pairs_hq.pdf` |
| QC: subregions | `R/97_if_subregions.R` (uses `roi_pairs_hq_annotate.pdf`) | `roi_pairs_subregions.pdf` |
| **Landmark sheets** | `R/99_landmark_sheets.R` | `landmark_sheets_marked.pdf` (user marks), `cache/register_if/landmark_panels.rds` |
| **Landmark fit** | `extract_landmarks.py` → `R/99b_landmark_fit.R` | `iftobf_lm_<sid>.rds` (IF→BF affine), `if_to_bf_landmark.pdf`, `if_to_bf_landmark_summary.csv` |
| **FINAL report** | `R/101_dataset_pairs_hq.R [all]` | **`dataset_pairs_hq.pdf`** |

### The manual-annotation loop (why there are PDFs to draw on)
1. **Brightfield rectangles** — the post-MALDI brightfield is faint, so `R/95` draws the
   MSI organoids as a guide. The user boxes the 6 organoids per slide in Adobe (the boxes
   were mis-numbered → `extract_bf_rectangles.py` applies a 1↔2/3↔4/5↔6 **pair-swap**, and
   corrects R's ~3.7% page-margin so box→pixel mapping is exact).
2. **Landmark marks** — on `landmark_sheets_marked.pdf` (grayscale IF | brightfield+MSI,
   equal-size panels) the user drops **numbered comments at organoid centers in BOTH panels**
   (same number = same organoid). `R/99` captured each panel's page-pt box via `grconvert`,
   so `extract_landmarks.py` + `R/99b` map marks → native pixels and fit a per-section affine.

### Final fit quality
`R/99b`: **10/12 sections** registered (affine, 5–8 organoid pairs each),
**RMSE 28–222 µm** — organoid-level (the realistic ceiling for serial sections).
`AO_20h_sl4b_sec1` and `sec5` were not marked → MSI datasets `sl4A_sec2a/2b` lack IF.

---

## 3. Final report — `dataset_pairs_hq.pdf` (`R/06_if_registration/12_dataset_pairs_hq.R`, legacy `R/101`)

One MSI dataset per page (20 pages), 3 panels, same region/scale (crop = MSI footprint
**+30%**), 200 µm bar:

1. **IF** — native-resolution 2-channel composite **ZO-1 (red) / DAPI (cyan)**
   (β-catenin and F-actin excluded), full opacity, reprojected into the brightfield frame
   via the landmark transform (not overlaid on brightfield).
2. **Brightfield** — native, with the **MSI-measured rectangle** as a thin red dashed line.
3. **Citrate** — MSI **anchored citrate [M−H]− 191.0198 ±7 ppm** ion image (viridis, 100% opacity) on the
   brightfield, with the **organoid segmentation outline** (`is_tissue` boundary = signed-
   distance-0 contour) as a **white dotted** line.

`R/100_dataset_pairs.R` is the earlier 2-panel version (superseded by R/101).

---

## 4. Re-running / extending

- **Fill the 2 missing IF datasets:** mark organoid centers on the sl4b sec1 & sec5 pages of
  `landmark_sheets_marked.pdf`, then:
  ```
  python extract_landmarks.py figures/register_if/landmark_sheets_marked.pdf
  Rscript R/99b_landmark_fit.R
  Rscript R/101_dataset_pairs_hq.R all
  ```
- **Different ion:** edit `ION_MZ` in `R/101`. **Different IF channels/colours:** the IF
  composite block in `R/101` (`IF_CH`, `IF_LUT`). **Crop size:** `EXT` (0.30). **Opacity:**
  `CIT_ALPHA`. Scripts honour an `HQ_OUT` env var for the output filename.

## 5. Superseded / abandoned (kept for history, NOT in the live path)
- `R/91_if_overview_to_bf.R` (FFT-NCC localize — only its boxes feed `R/94` section crops),
  `R/92_if_overview_to_msi.R`, `R/93_if_native_overlay.R` — the abandoned *automatic*
  overview→BF registration (failed: faint cross-modal BF + serial sections).
- `R/98_if_to_bf_register.R` — automatic organoid-shape fit (failed: IF tissue extent ≫ MSI
  subset; moment-scale collapsed). Replaced by the landmark approach above.
