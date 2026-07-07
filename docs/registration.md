# MSI ↔ brightfield registration pipeline

> **Path note:** this doc uses the old flat script numbers (`R/30b`…`R/35`). In the fork these
> live under **Phase 03** (`R/03_coarse_registration/`) and **Phase 05**
> (`R/05_registration_refine/`); shared helpers are `R/00_lib/lib_register.R`. See the
> legacy→new table in `PIPELINE.md` to map each number.

Registers each MSI section's data onto the native-resolution Nikon `.nd2` whole-slide brightfield,
using the Bruker flexImaging teach marks as the foundation and the MSI `is_tissue` mask for fine
refinement. Output = per-section native optical crops with MSI projected on top.

> **Status: working, validated.** Earlier one-shot attempt (`R/30d`, `R/30e`, now deleted) failed
> because a spurious JPG→`.nd2` offset and an over-large refinement search compounded. Rebuilt
> step-by-step with user validation; the key was refining in **small steps on validated ground**.

## Data

| Item | Location | Notes |
|---|---|---|
| Whole-slide BF (native) | `…/06102026_AO_0h_sl6A.nd2`, `…_20h_sl4A.nd2` | 23163×11784 uint16, ~1.83 µm/px. **One `.nd2` per slide**; its 10 sections share it. The 3 "series" are NOT a usable pyramid (series>1 read fails) → read `series=1` with `subset=list(x=,y=)`. |
| Slide BF (downscaled) | `MSI/06102026_AO_*_small.jpg` | 8000×4070; the `.nd2` field downscaled **×2.895** and **vertically flipped**. The Bruker teach marks reference this image. |
| Imaging sequence | `MSI/06102026_AO_*/*_sec*.mis` (sl6A) / `*_<code>.mis` (sl4A) | XML: `<TeachPoint>` flex↔stage, `<OriginalImageTeachPoint>` flex↔`_small.jpg`, `<Area>` = measured ROI in flex coords. |
| Stage log | `*_poslog.txt` | Per-spot stage positions. **NOT used** — its stage frame differs from the teach stage frame (different sign/offset). |

## Coordinate chain

```
MSI(x,y) ──teach (.mis Area + OriginalImageTeachPoint)──▶ slide _small.jpg
         ──V-flip + ×2.895 (≈0 offset)────────────────▶ native .nd2
```

Key facts (validated):
- The `.mis` `<Area>` bbox **exactly equals the raster extent** (sec1a: 224 flex-px = 60 raster ×
  10 µm ÷ 2.679 µm/flex-px). So MSI grid maps linearly into the Area bbox.
- MSI→raster orientation is **(+1,+1)** for both slides (global, by tissue-on-organoid contrast).
- JPG→`.nd2` is **vertical flip + scale 2.895, offset ≈ 0** (8/8 and 4/8 `.nd2` px). Scale = `nd2_W/jpg_W`.
- Within-box and `.nd2` refinement is **small translation only** (≤3 MSI px), minor rot/scale — the
  teach foundation is correct (user-confirmed: boxes land on the right organoids).

## Pipeline (run in order)

All with **R 4.4.2** (`C:/Program Files/R/R-4.4.2/bin/Rscript.exe`); needs `RBioFormats` +
`JAVA_HOME=C:/Program Files/Java/jre1.8.0_491`. Shared helpers in `R/lib_register.R`.

| # | Script | Input | Output | Purpose |
|---|---|---|---|---|
| 0a | `R/30b_teach_msi_to_jpg.R` | `peaks_tissue_combined.rds`, `.mis`, jpg | `cache/register/teach_<sid>.rds`, `results/registration/teach_summary.csv` | Coarse MSI→jpg from teach marks (Area bbox + OriginalImageTeachPoint), global orientation. |
| 0b | `R/30c_jpg_to_nd2.R` | `.nd2`, jpg | `cache/register/nd2thumb_<slide>.rds` | Block-MEAN `.nd2` thumbnail (factor 16; strided aliases thin organoids). *(Its flip/offset output is superseded by R/33.)* |
| 1 | `R/31_slide_overview.R` | `peaks_tissue_combined.rds`, `.mis`, jpg | `figures/registration/slide_overview_<slide>.pdf` | **Slide-level QC**: measured areas + tissue masks on slide BF (validation only). |
| 2 | `R/32_refine_jpg.R` | `teach_<sid>.rds`, jpg | `cache/register/jpgxform_<sid>.rds`, `results/registration/refine_jpg_summary.csv`, `figures/registration/refine_jpg_<slide>.pdf` | Within-box small-translation refine on JPG (±12 MSI px, rot ±4°, scale ±3%). |
| 3 | `R/33_jpg_to_nd2_offset.R` | `jpgxform_<sid>.rds`, `nd2thumb_<slide>.rds`, jpg | `cache/register/jpg2nd2off_<slide>.rds`, `figures/registration/slide_nd2_<slide>.pdf` | Global JPG→`.nd2` offset (V-flip + scale + tx,ty) by max total tissue contrast. |
| 4 | `R/34_native_crops.R` | `jpgxform`, `jpg2nd2off`, `.nd2` | `cache/register/nd2final_<sid>.rds`, `figures/registration/crops/optical_<sid>.png`, `figures/registration/registration_native.pdf`, `results/registration/native_summary.csv` | Compose MSI→`.nd2`, read native crop, tiny ±3 px polish, save crop + report. |
| 5 | `R/35_overlay_report.R` | `nd2final_<sid>.rds`, `optical_<sid>.png` | `figures/registration/registration_native.pdf` | **Report-only** (no `.nd2` re-read): native BF \| +`is_tissue` \| + MSI ion overlay. |

```bash
RS="C:/Program Files/R/R-4.4.2/bin/Rscript.exe"
$RS R/30b_teach_msi_to_jpg.R && $RS R/30c_jpg_to_nd2.R && $RS R/31_slide_overview.R \
 && $RS R/32_refine_jpg.R && $RS R/33_jpg_to_nd2_offset.R && $RS R/34_native_crops.R \
 && $RS R/35_overlay_report.R
```

## Tunables

- `R/35`: `OVERLAY_ALPHA` (constant MSI opacity, **0.65**); whole-MSI-image vs on-tissue (one line:
  `lut[cbind(sub$x,sub$y)]<-secval` for whole image, `sub$is_tissue` subset for on-tissue); `ION_MZ`
  (default citrate 191.0217); colormap (viridis). **Overlay style is locked: homogeneous opacity,
  whole MSI image** (see memory `feedback_msi_bf_overlay`).
- `R/32`: `TRANS_MSI` (±12), `ROTS`, `SCALES`. `R/34`: `POLISH_MSI` (±3).

## Per-section transform (downstream use)

`cache/register/nd2final_<sid>.rds` holds `B_msi_nd2` (3×2 affine MSI→`.nd2` px), `crop`
(nd2 bbox of the saved optical PNG), `scale_msi_nd2` (~5.4 nd2 px/MSI px). To project any MSI ion
onto native BF: invert `B_msi_nd2` per crop pixel → MSI (x,y) → sample ion (see `R/35`).

## Validation

All 20 sections registered; refinement polish ≤3 MSI px everywhere; tissue-on-organoid contrast
0.05–0.13. QC: `registration_native.pdf` (masks trace organoid rings at native res), plus the
step PDFs (`slide_overview_*`, `refine_jpg_*`, `slide_nd2_*`).
