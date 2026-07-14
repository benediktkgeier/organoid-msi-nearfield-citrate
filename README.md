# Organoid MSI — apical-out vs basolateral-out near-field citrate

MALDI-timsTOF flex MSI (negative mode, *m/z* 100–900, DAN matrix, 10 µm pixel, 4% CMC embedding)
of intestinal-organoid sections. All sections were **incubated in CMC for at least 5 min before
freezing and are treated identically** (no incubation-time comparison). Research question: does
organoid apical polarity (**apical-out** vs **basolateral-out**) relate to near-field citrate
emission into the surrounding CMC gel? Two physical slides carry 10 sections each → **20 datasets**.

This repository is the **publication fork** (`Analysis_R_Final`). It reuses the validated, timepoint-
agnostic intermediate caches from the working analysis (read-only) and re-runs only the analysis /
figure layer. It is organized as a numbered, reproducible pipeline under `R/`; see
**[`PIPELINE.md`](PIPELINE.md)** for the authoritative phase → script → input → output map, and
**[`docs/`](docs/)** for per-topic detail.

> **Reproducing this pipeline from the paper?** Start at **[`REPRODUCE.md`](REPRODUCE.md)** (data
> layout, run modes, validation) with the exact version pin in
> **[`docs/ENVIRONMENT.md`](docs/ENVIRONMENT.md)**.
>
> **Data availability:** the **MSI datasets (imzML)** are on METASPACE —
> **[McKenna-2026 project](https://metaspace2020.org/project/McKenna-2026?tab=datasets)**. This repo is
> code + docs + small inputs; the intermediate cache, `.nd2` microscopy, and rendered figures live in a
> separate data deposit — link to be added (see [`REPRODUCE.md`](REPRODUCE.md#data-availability)).

## Environment (locked)

- **R 4.4.2** — `C:/Program Files/R/R-4.4.2/bin/Rscript.exe` (NOT 4.5.3). Cardinal v3 (Bioconductor).
- Parallelism via `register_parallel()` (SnowParam SOCK, 4 workers) in `R/00_lib/lib_paths.R`.
  **Never** call `setCardinalBPPARAM()` (Cardinal v3 bug). `JAVA_HOME=C:/Program Files/Java/jre1.8.0_491`
  for RBioFormats.
- Packages: `Cardinal`, `matter`, `BiocParallel`, `RBioFormats`, `EBImage`, `RANN`, `viridisLite`,
  `jpeg`, `png`, `pdftools`.

## Data reuse (this fork)

`R/00_lib/lib_paths.R` defines `CACHE_SRC` (the original `Analysis_R/cache`, read-only) and a
`cache_in()` resolver: upstream caches (peaks, citrate, instances, zones, registration, IF) are read
in place; any recomputation writes to this fork's own `cache/`. The original `Analysis_R` is never
modified.

## Pipeline phases

| # | Phase | Folder | Headline output |
|---|---|---|---|
| 01 | Spectral preprocessing | `R/01_preprocess/` | `cache/peaks_combined.rds` (self-contained MSE) |
| 02 | Citrate standard GATE | `R/02_CitrateStandard/` | locked `CITRATE_ANCHOR_MZ` ±7 ppm; `cache/citrate_anchored_<sid>.rds` |
| 03 | Coarse MSI→brightfield registration | `R/03_coarse_registration/` | coarse MSI→`.nd2` transforms (`cache/register/`) |
| 04 | SSC on-tissue delineation | `R/04_ssc_ontissue/` | `cache/peaks_tissue_combined.rds` (`is_tissue`, `chemotype`) |
| 05 | Registration refinement | `R/05_registration_refine/` | native crops + reusable MSI→`.nd2` affines |
| 06 | IF↔brightfield/MSI registration | `R/06_if_registration/` | `figures/if_registration/dataset_pairs_hq.pdf` |
| 07 | Metabolite identification + citrate QC | `R/07_metabolite_id/` | `results/metabolites/`, citrate-resolution reports |
| 08 | Organoid segmentation + pooled gradient survey | `R/08_organoid_gradient_survey/` | `figures/gradient/gradient_report.pdf` (pooled, descriptive) |
| 09 | Organoid separation refinement | `R/09_organoid_refinement/` | `cache/instances_final_<sid>.rds` (curated ROIs) |
| 10 | Apical-orientation annotation | `R/10_apical_annotation/` | `results/annotation/apical_*` (apical-out / basolateral-out / mixed) |
| 11 | Per-organoid final + statistics | `R/11_per_organoid_final/` | `figures/annotation/apical_citrate_dha_report.pdf` (page 1 = Summary Statement) + `apical_nearfield_emission_figure.pdf` |

`R/00_lib/` holds `lib_paths.R` (locked constants + `cache_in()`/`CACHE_SRC` + `disp_id()` + `APICAL_COLS`),
`gradient_config.R`, `if_config.R`, `lib_register.R`, `lib_register_if.R`, `lib_nearfield_viz.R`,
`lib_nearfield_viz_whole.R`, `lib_gradient_seg.R`, `lib_citrate.R`.

Each phase folder (`R/NN_*/`) additionally carries its own `README.md`; **[`PIPELINE.md`](PIPELINE.md)**
is the cross-phase authority.

### Apical / colour convention (locked)

Polarity classes and their figure colours are defined once in `lib_paths.R` (`APICAL_COLS`) and reused
everywhere (outlines, dot/box plots, near-field overlays, legends): **apical-out = magenta `#C2399A`**,
**basolateral-out = green `#2C9E4B`**, **mixed = grey `#888888`**. ("apical-out" = apical membrane faces
the gel; "basolateral-out" = basolateral faces the gel, i.e. apical faces the lumen.)

## Running

- **`./run_all.sh`** — cache-reuse figure regenerator: runs the analysis/figure layer in phase order,
  reusing the upstream caches read-only. Manual annotation gates (on-tissue ion selection, organoid
  split/island markup, apical comments) are already cached/reused. Pass the **`qc`** arg
  (`./run_all.sh qc`) to also rebuild the citrate/metabolite QC reports.
- **`./run_apical_nearfield_pipeline.sh`** — reproduces the per-organoid apical → near-field analysis
  and the locked four-view publication figure.
- **`./run_wholeregion_figures.sh`** — whole-region near-field figures + the all-20 per-dataset citrate
  gradient report.

## Key result

Apical-out organoids show elevated TIC-normalised citrate within ≤50 µm of the organoid surface
(near-field level, `grad_near50`; apical-out > basolateral-out, p ≈ 6.07e-5, large Cliff's δ),
corroborated by the citrate/DHA internal-control ratio (which cancels thickness/ionisation). The outward monotonic slope
(`grad_rho_out`) is not significant, so the result is a near-field **level magnitude**, not gradient
steepness. Citrate = `[M−H]⁻` anchored 191.0198 ±7 ppm (raw imzML), TIC-normalised; DHA 327.2330 =
confined-lipid negative control. Detail in [`docs/apical_nearfield.md`](docs/apical_nearfield.md).

## Documentation

- **[`REPRODUCE.md`](REPRODUCE.md)** — reproduction protocol (data layout, run modes, validation).
- **[`docs/ENVIRONMENT.md`](docs/ENVIRONMENT.md)** — R/Bioconductor/package version pin + install.
- **[`PIPELINE.md`](PIPELINE.md)** — authoritative phase → script → I/O map.
- **[`docs/DESIGN_DECISIONS.md`](docs/DESIGN_DECISIONS.md)** — all locked design decisions.
- **[`docs/ssc_ontissue.md`](docs/ssc_ontissue.md)** — SSC on-tissue clustering, floor80 rule, chemotypes, organoid segmentation + the visual SSC report (phase 04).
- **[`docs/registration.md`](docs/registration.md)** — MSI↔brightfield registration (phases 03/05).
- **[`docs/if_registration.md`](docs/if_registration.md)** — IF↔brightfield/MSI registration (phase 06).
- **[`docs/roi_curation.md`](docs/roi_curation.md)** — organoid ROI curation (phase 09).
- **[`docs/apical_annotation.md`](docs/apical_annotation.md)** + **[`docs/apical_nearfield.md`](docs/apical_nearfield.md)** — phases 10/11.
- **[`docs/nearfield_wholeregion.md`](docs/nearfield_wholeregion.md)** — whole-region near-field figure (all 20).
- **[`docs/citrate_gradient_perdataset.md`](docs/citrate_gradient_perdataset.md)** — per-dataset citrate/DHA gradient (v3, all 20).
- **[`docs/citrate_isotopes_adducts.md`](docs/citrate_isotopes_adducts.md)** — citrate isotope/adduct analysis.
- **[`docs/HISTORY.md`](docs/HISTORY.md)** — project history / changelog.
- **`inventory.csv`** — sample registry (sample_id, imzml_path, bf_jpg_path, nd2_path).

## Methods (incubation)

> Organoids were incubated in carboxymethylcellulose (CMC) for at least 5 minutes prior to freezing;
> all sections were treated identically (no incubation-time comparison).

## Datasets (20)

20 organoid sections across two physical slides (sl6A, sl4A), 10 sections each, all treated identically.
Statistics are descriptive (organoid = independent unit; pairwise Wilcoxon, Cliff's δ); no group
contrasts.
