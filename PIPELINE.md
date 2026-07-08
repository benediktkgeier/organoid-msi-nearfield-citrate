# PIPELINE.md — authoritative phase → script → I/O map

R 4.4.2 only. All scripts source `R/00_lib/lib_paths.R` (locked constants). Run order is the folder
order (01→11); within a phase, the file-number order. `R/00_lib/` = shared libs/configs.

This is the **cache-reuse publication fork** (`Analysis_R_Final`): `lib_paths.R` defines `CACHE_SRC`
(the untouched upstream `Analysis_R/cache`, read-only) and a `cache_in()` resolver, so upstream
intermediate caches (peaks, citrate, instances, zones, registration, IF) are read in place and only
the analysis/figure layer is re-run. This fork's own `cache/` starts empty; `.gitignore` excludes
`cache/`, `figures/`, and large `results/`. Dropped/deprecated scripts (`R/_deprecated/`,
`99_scratch/`, the split-by-group Phase-08 `05–08`, and the superseded Phase-11 helpers) were **not**
carried into the fork — see the legacy→new table for what maps where.

**Per-phase READMEs:** each `R/NN_*/` folder has its own `README.md` describing that phase's scripts
and I/O; this file is the cross-phase authority.

> **Legacy script numbers in code comments:** many inline comments still cite the old flat numbers
> (e.g. "from R/40", "see R/87"). Use the **Legacy → new mapping** table at the bottom to decode them.

## Phase 01 — Spectral preprocessing
| Script | In | Out |
|---|---|---|
| `01_preprocess.R` | imzML (cache/imzml) + `inventory.csv` | aligned MSE → `cache/peaks_after_freq.rds` |
| `02_realize_in_memory.R` | peaks MSE | in-memory self-contained MSE |
| `03_refilter_freq_floor.R` | realized MSE | **`cache/peaks_combined.rds`** (8,772 × 90,760) |
| `qc_render.R` | `peaks_combined.rds` | `figures/preprocess/qc_combined.pdf`, feature-filter diagnostic |

Cardinal v3 idiom (locked): `readMSIData` → `c()` combine → `normalize(tic)` → `process()` →
manual single-linkage reference grid → `convertMSImagingArrays2Experiment` → `peakAlign` →
`summarizeFeatures`. **Never** `setCardinalBPPARAM()`.

## Phase 02 — Citrate standard GATE (lock anchored citrate; build per-sample cache)
**Critical early gate.** Validates citrate against an authentic Sodium Citrate Tribasic Dihydrate
dilution series and LOCKS the citrate definition for the whole pipeline: `CITRATE_ANCHOR_MZ =
191.01976` (+0.16 ppm vs theo), `CITRATE_WIN_PPM = 7` (lib_paths.R). Only proceed to downstream image
analysis if the spectra make sense.

| Script | Role |
|---|---|
| `00_config.R` | spot<->conc table, masses; sources `R/00_lib/lib_citrate.R` (reader + helpers) |
| `01_standard_spectra_mz.R` | pure-standard centroid +0.16 ppm, FWHM ~7 ppm; tissue 191 +14.6 ppm / 3.2x broader = co-isobar |
| `02_calibration_curve.R` | log-log slope 1.00, R2 0.987 (1-100 mM); LOD ~0.37 mM; both decade pairs were mislabeled (corrected on disk) |
| `03_id_fingerprint.R` | standard 13C 7.3% ~ 6.8% theo + clean adduct profile; tissue inflated = co-isobar |
| `04_standard_anchored_citrate.R` | GATE self-check (measured anchor ~ CITRATE_ANCHOR_MZ); 5/7/10 ppm sweep -> **+-7 ppm** (max citrate capture, 0% shoulder) |
| `05_build_citrate_cache.R` | writes `cache/citrate_anchored_<sid>.rds` (x,y,cit_raw,tic) for all 20 samples |

Why early + separate: the processed 25-ppm grid merges citrate (191.0198) with a +16 ppm co-isobar
into ONE feature (the old 191.0217), so it cannot represent a +-7 ppm citrate window. The gate reads
the RAW centroid imzML and caches the anchored citrate; **ALL downstream citrate extraction**
(phases 05/06/08/11 + `R/00_lib/lib_nearfield_viz.R`) now calls `citrate_onto_pd()` from
`R/00_lib/lib_citrate.R` instead of the grid feature (TIC-normalised; joined on sample_id+x,y;
verified 0% unmatched). CAVEAT (centroided data): a narrow window SELECTS citrate-dominant pixels, it
does not integrate citrate area — true integration needs profile spectra. Comparison design: standard
is matrix-matched (same CMC, same slide sl7A, same spray); on-sample co-acquisition not possible.
Outputs in `results/citrate_standard/`, `figures/citrate_standard/`.

## Phase 03 — Coarse MSI→brightfield registration
| Script | Out |
|---|---|
| `01_teach_msi_to_jpg.R` | MSI→slide-JPG teach transform |
| `02_jpg_to_nd2.R` | JPG→`.nd2` block-mean thumbs / coarse transform |
| `03_slide_overview.R` | slide-level QC overview |

Transforms cached in `cache/register/`. Detail: [`docs/registration.md`](docs/registration.md).

## Phase 04 — SSC on-tissue delineation (selected on-tissue ions)
| Script | In | Out |
|---|---|---|
| `01_export_ontissue_candidates.R` | `peaks_combined.rds` | candidate ion sheet (**MANUAL** on-tissue selection) |
| `02_import_ontissue_ions.R` | selection (`results/peakme_annotations/`) | `cache/peaks_combined_annot.rds` |
| `03_build_curated_set.R` | annot MSE + published m/z list | **`cache/peaks_curated.rds`** (348 ions) |
| `04_ssc_tissue_mask.R` | `peaks_curated.rds` | **`cache/peaks_tissue_combined.rds`** (`is_tissue`, `ssc_k4_sec`, `ssc_k10_sec`, `chemotype`, `tissueness_px`) + per-section `figures/ssc/ssc_mask_<sid>.pdf` |
| `05_ssc_report.R` | `peaks_tissue_combined.rds` + `zones_<sid>.rds` | **`figures/ssc/ssc_clustering_segmentation_report.pdf`** — visual SSC + organoid-segmentation report (methods/composition front page + 4-panel per dataset); reads cached columns, no recompute |

On-tissue rule "floor80" + 4 chemotypes — see [`docs/DESIGN_DECISIONS.md`](docs/DESIGN_DECISIONS.md).
`results/peakme_annotations/` is the legacy storage name for the manual on-tissue ion selection.

## Phase 05 — Registration refinement
`01_refine_jpg.R` → `02_jpg_to_nd2_offset.R` → `03_native_crops.R` → `04_overlay_report.R`
→ `figures/registration/registration_native.pdf` + native crops + reusable MSI→`.nd2` affines
(`cache/register/nd2final_<sid>.rds`). Detail: [`docs/registration.md`](docs/registration.md).

## Phase 06 — IF (B-section) → brightfield/MSI registration
`01_if_thumbs` (+`01b`,`01c` hi-res raster builds) → `02`–`09` overview/native/crops/subregions →
`10_landmark_sheets` → `11_landmark_fit` → **`12_dataset_pairs_hq.R`** →
`figures/if_registration/dataset_pairs_hq.pdf`. Manual landmarks (cached).
Detail: [`docs/if_registration.md`](docs/if_registration.md).

## Phase 07 — Metabolite identification + citrate QC
`01_metabolite_match.R` (published m/z list → 149 matched) → `02_metabolite_report.R`
(+ `02b_metabolite_report_single.R`, single-section variant) →
citrate QC: `03_citrate_resolution.R`, `04_citrate_window_images.R`, `05_citrate_isotopes_adducts.R`.
Outputs in `results/metabolites/`, `figures/metabolites/`.
Detail: [`docs/citrate_isotopes_adducts.md`](docs/citrate_isotopes_adducts.md).

### `R/02_CitrateStandard/` — citrate standard DETAIL (the GATE is **Phase 02** above; window LOCKED at +-7 ppm and DRIVES all downstream citrate via lib_citrate; the "+-5 ppm" / "standalone QC" notes below are superseded)
Uses the Sodium Citrate Tribasic Dihydrate dilution series spotted in the same CMC matrix
(`MSI/06102026_AO_0h_sl7A/imzml`: 0/blank, 10 µM, 100 µM, 1 mM, 10 mM, 100 mM).
`00_config.R` (spot↔conc table, masses, .ibd reader, tissue/same-run hook) →
`01_standard_spectra_mz.R` (pure-standard centroid **+0.2 ppm** vs theo, FWHM ~7 ppm; tissue overlay:
+14.6 ppm shift, 3.2× broader = co-isobar quantified) →
`02_calibration_curve.R` (log–log **slope 1.00, R² 0.987** over 1–100 mM; **both decade pairs were
mislabeled at acquisition — 10↔100 µM and 10↔100 mM — `00_config.R` maps each file to its TRUE
concentration**; LOD ~0.37 mM, LOQ ~1.1 mM, matrix-limited by a 191 background ion in citrate-free
CMC; dimer not detected) →
`03_id_fingerprint.R` (standard 13C **7.3% ≈ 6.8% theo** + clean adduct profile; tissue inflated =
co-isobar signature) →
`04_standard_anchored_citrate.R` (anchor tissue citrate on the standard-measured mass **191.01976
(+0.16 ppm)**, ±5 ppm window, vs the old 191.0217 feature). KEY: data are **centroided** — the
NEW/OLD windows select **disjoint pixels** (co-occurrence ≈ 0), so a narrow window **selects
citrate-dominant pixels, it does not integrate citrate area**; a true integration needs **profile**
spectra. Anchor + window saved to `results/citrate_standard/citrate_anchor.csv`.
Outputs in `results/citrate_standard/`, `figures/citrate_standard/`.
**Comparison design:** on-sample (same-run) organoid co-acquisition is NOT possible; the standard
stands as the reference because it is **matrix-matched** — same CMC, **same slide (sl7A)**, same spray
as the organoid sections. Default tissue side = the 5 citrate-positive organoid sections available as
imzML (sl6A/sl4A) — same prep, separate acquisition run ("sep-run"). For a literal **same-slide**
comparison, export sl7A's organoid sections from the SCiLS `_sections.sbd` to imzML and list them in
`SAMERUN_SECTIONS` (00_config.R) → comparisons relabel "sep-run" → "same-slide".
Standalone QC — does **not** modify any organoid/gradient script.

## Phase 08 — Organoid segmentation + gradient survey (pooled, all sections)
| Script | Role |
|---|---|
| `01_segment_organoids.R` | connected components of `is_tissue` → `cache/instances_<sid>.rds` |
| `02_buffer_rings.R` | signed-distance rings (outward `BUF_LEVELS_UM`, inward 10 µm) → `cache/zones_<sid>.rds` |
| `03_gradient_stats.R` | per-ion pooled ρ_out / ρ_in across ALL sections → `results/gradient/` |
| `04_report_pdf.R` | `figures/gradient/gradient_report.pdf` (pooled ρ ranking + profiles) |

Single-condition study: this is a POOLED descriptive survey (no group split, no Δρ). Supporting
context only — the headline result is the Phase 11 apical near-field analysis. (The former split-by-
group `05_citrate_20datasets` / `06_citrate_dha_compare` / `07_citrate_dha_20datasets` /
`08_gradient_20datasets` scripts were dropped in this fork.)

## Phase 09 — Organoid separation refinement
`00_run_roi_curation.R` chains: `01_organoid_split_apply` (split merged; **MANUAL** green-ink cuts,
cached) → `03_island_cleanup_apply` (delete/merge; **MANUAL** actions, cached) →
`04_finalize_instances` (one connected ROI = one id). `02_island_cleanup_canvas` + `02b_island_centroid_sidecar`
build the annotation canvas + centroid sidecar. → `cache/instances_final_<sid>.rds`.
Detail: [`docs/roi_curation.md`](docs/roi_curation.md).

## Phase 10 — Apical-orientation annotation
`01_apical_annotate.R` (emits markup PDF; **MANUAL** apical-out / basolateral-out / mixed) → `02_apical_parse.R` →
`03_apical_consensus.R` (two-annotator) → `04_apical_consensus_finalize.R` →
`results/annotation/apical_*` (per-organoid apical class keyed by sid+instance).
Detail: [`docs/apical_annotation.md`](docs/apical_annotation.md).

## Phase 11 — Per-organoid final gradient + statistics
Apical classes come from the **two-annotator CONSENSUS** map
`results/annotation/apical_map_consensus.csv` (a static committed input; the consensus was finalized
upstream). It is the SOLE annotation basis for every report/figure here — no marked-PDF re-parse
(and therefore no python/reticulate) is needed at runtime.

| Script | Role |
|---|---|
| `01_zones_curated.R` | curated Voronoi zones on **refined** instances (all 20) → `cache/zones_<sid>.rds` |
| `03_apical_report.R` | **headline** report (default = consensus map) + `apical_gradient_per_organoid.csv`, `apical_citrate_dha_stats.csv` |
| `04_nearfield_figure.R` | LOCKED four-view near-field publication figure (`lib_nearfield_viz.R`) |
| `07_gradient_test_perorganoid.R` | per-organoid Wilcoxon near-field test |
| `08_nearfield_wholeregion.R` | **LOCKED** whole-region near-field figure, all 20 datasets (`lib_nearfield_viz_whole.R`); `all globalheat` → the global-scale `..._globalheat.pdf` variant |
| `08b_view3_heatmap_hires.R` | hi-res standalone render of the whole-region weather-heatmap (view 3) |
| `09_citrate_gradient_perdataset_v3.R` | **LOCKED** per-dataset citrate gradient — **all 20 datasets** |
| `09b_citrate_ion_hires.R` | hi-res (600 dpi) citrate ion-image PNGs per dataset |
| `09_prism_export.R` | GraphPad Prism (xlsx/pzfx) export of the apical citrate/DHA stats |
| `10_citrate_gradient_report_final.R` | **FINAL combined report** → `figures/gradient/citrate_gradient_report_final.pdf`: 6-panel per dataset (MSI / native BF / SSC mask / matched IF / MSI↔BF overlay / weather gradient map) + ≤100 µm apical-class overview page. Non-locked driver; sources the locked libs + copies their small helpers |
| `10_pzfx_export.R` | Prism `.pzfx` export companion to `09_prism_export` |
| `11_export_gradient_profile.R` | tabular export of per-organoid outward gradient profiles |
| `12_joy_tables.R` | curated "Joy_Tables" subset export for the collaborator |
| `13_citrate_gradient_report_3class.R` | apical-CLASS variant of `10_` (untouched sibling): default → `citrate_gradient_report_withmixed.pdf` (mixed = its own grey trend line, 3-group); `nomixed` arg → `citrate_gradient_report_nomixed.pdf` (2-group, mixed dropped). Both publication options |
| `14_citrate_gradient_report_final_msigrid.R` | MSI-grid-outline variant of `10_`: panel 1 also overlays the per-instance segmentation outline in native MSI-grid space (`class_outlines(native=FALSE)`, class-coloured) so the un-warped shape is comparable to the BF-warped outlines in panels 5/6 → `citrate_gradient_report_final_msigrid.pdf` |

Headline: `figures/annotation/apical_citrate_dha_report.pdf` (page 1 = Summary Statement) +
`results/annotation/apical_citrate_dha_stats.csv` — the consensus source of the near-field claim
(grad_near50 apical-out > basolateral-out, p ≈ 6e-5).

(Dropped in this fork: `02_apical_citrate_dha` [single-annotator, superseded by the consensus map],
`03b_apical_report_50um_nomixed`, `05_apical_toppick_report`,
`06_citrate_gradient_perdataset` [TEST-scope, superseded by `09`].)

`10_citrate_gradient_report_final.R` is the non-locked **final combined** report (recombines the two
LOCKED scripts + the SSC/BF/IF modalities into one per-dataset deliverable); it and the `09b/10_pzfx/
11/12` export scripts are cache-only and run independently. Shared numbers (`09_`×2, `10_`×2) pair a
figure/report with its data export — run either directly.

Detail: [`docs/apical_nearfield.md`](docs/apical_nearfield.md) (per-organoid),
[`docs/nearfield_wholeregion.md`](docs/nearfield_wholeregion.md) (whole-region),
[`docs/citrate_gradient_perdataset.md`](docs/citrate_gradient_perdataset.md) (v3 all-20).
Runner: `./run_wholeregion_figures.sh`; the final combined report:
`Rscript R/11_per_organoid_final/10_citrate_gradient_report_final.R all`.

---

## Legacy → new mapping (old flat number → new path)

| Legacy | New path |
|---|---|
| R/01_preprocess | R/01_preprocess/01_preprocess.R |
| R/01b_qc_render | R/01_preprocess/qc_render.R |
| R/01c_top_ion_images_clean | *(dropped — top-ion QC not carried into fork)* |
| R/06_realize_in_memory_20section | R/01_preprocess/02_realize_in_memory.R |
| R/07_refilter_freq_floor | R/01_preprocess/03_refilter_freq_floor.R |
| R/30b_teach_msi_to_jpg | R/03_coarse_registration/01_teach_msi_to_jpg.R |
| R/30c_jpg_to_nd2 | R/03_coarse_registration/02_jpg_to_nd2.R |
| R/31_slide_overview | R/03_coarse_registration/03_slide_overview.R |
| R/05_peakme_export_3section_scils | R/04_ssc_ontissue/01_export_ontissue_candidates.R |
| R/04_import_annotations | R/04_ssc_ontissue/02_import_ontissue_ions.R |
| R/21_build_curated_set | R/04_ssc_ontissue/03_build_curated_set.R |
| R/20_ssc_tissue_mask | R/04_ssc_ontissue/04_ssc_tissue_mask.R |
| R/32_refine_jpg | R/05_registration_refine/01_refine_jpg.R |
| R/33_jpg_to_nd2_offset | R/05_registration_refine/02_jpg_to_nd2_offset.R |
| R/34_native_crops | R/05_registration_refine/03_native_crops.R |
| R/35_overlay_report | R/05_registration_refine/04_overlay_report.R |
| R/90_if_thumbs | R/06_if_registration/01_if_thumbs.R |
| R/90b_build_hr4raw / R/90c_build_hrdapihq | R/06_if_registration/01b_ / 01c_ |
| R/91…99,99b | R/06_if_registration/02_…11_ (overview→landmark_fit) |
| R/101_dataset_pairs_hq | R/06_if_registration/12_dataset_pairs_hq.R |
| R/80_metabolite_match | R/07_metabolite_id/01_metabolite_match.R |
| R/81_metabolite_report | R/07_metabolite_id/02_metabolite_report.R |
| R/82_citrate_resolution | R/07_metabolite_id/03_citrate_resolution.R |
| R/83_citrate_window_images | R/07_metabolite_id/04_citrate_window_images.R |
| R/83b_citrate_isotopes_adducts | R/07_metabolite_id/05_citrate_isotopes_adducts.R |
| R/40_segment_organoids | R/08_organoid_gradient_survey/01_segment_organoids.R |
| R/50_buffer_rings | R/08_organoid_gradient_survey/02_buffer_rings.R |
| R/60_gradient_stats | R/08_organoid_gradient_survey/03_gradient_stats.R |
| R/70_report_pdf | R/08_organoid_gradient_survey/04_report_pdf.R |
| R/82_citrate_20datasets | *(dropped — split-by-group, not carried into fork)* |
| R/84_citrate_dha_compare | *(dropped — split-by-group, not carried into fork)* |
| R/85_citrate_dha_20datasets | *(dropped — split-by-group, not carried into fork)* |
| R/86_gradient_20datasets | *(dropped — split-by-group, not carried into fork)* |
| R/45_run_roi_curation | R/09_organoid_refinement/00_run_roi_curation.R |
| R/41_organoid_split_apply | R/09_organoid_refinement/01_organoid_split_apply.R |
| R/42_island_cleanup_canvas | R/09_organoid_refinement/02_island_cleanup_canvas.R |
| R/42b_island_centroid_sidecar | R/09_organoid_refinement/02b_island_centroid_sidecar.R |
| R/43_island_cleanup_apply | R/09_organoid_refinement/03_island_cleanup_apply.R |
| R/44_finalize_instances | R/09_organoid_refinement/04_finalize_instances.R |
| R/102_organoid_apical_annotate | R/10_apical_annotation/01_apical_annotate.R |
| R/103_organoid_apical_parse | R/10_apical_annotation/02_apical_parse.R |
| R/108_apical_consensus | R/10_apical_annotation/03_apical_consensus.R |
| R/109_apical_consensus_finalize | R/10_apical_annotation/04_apical_consensus_finalize.R |
| R/51_zones_curated | R/11_per_organoid_final/01_zones_curated.R |
| R/104_apical_citrate_dha | *(dropped — single-annotator, superseded by the consensus map)* |
| R/105_apical_report | R/11_per_organoid_final/03_apical_report.R |
| R/106_nearfield_figure | R/11_per_organoid_final/04_nearfield_figure.R |
| R/110_apical_toppick_report | *(dropped — not carried into fork)* |
| R/87b_citrate_gradient_perdataset_v2 | *(dropped — TEST scope, superseded by 09_citrate_gradient_perdataset_v3.R)* |
| R/88_gradient_test_perorganoid | R/11_per_organoid_final/07_gradient_test_perorganoid.R |
| R/00_rebuild_5section, R/01c (non-clean), R/02*/03 export variants, R/05b, R/10_build_reference_peaks, R/22, R/87 (v1), R/89, R/100, probe_* | *(deprecated upstream; `_deprecated/` and `99_scratch/` were NOT carried into this fork)* |

Libraries (`R/00_lib/`): `lib_paths.R` (locked constants + `cache_in()`/`CACHE_SRC` + `disp_id()` +
`APICAL_COLS`), `lib_citrate.R` (`citrate_onto_pd()`), `gradient_config.R`, `if_config.R`,
`lib_register.R`, `lib_register_if.R`, `lib_nearfield_viz.R`, `lib_nearfield_viz_whole.R`,
`lib_gradient_seg.R`, `lib_report_frame.R` (shared report-frame primitives for `04/05` + `11/10`).

Non-phase generator (outside the `R/` tree, run by `run_all.sh`):
`figures/methods_report/generate_methods_report.R` → `figures/methods_report/methods_pipeline_report.pdf`
— a self-contained methods/pipeline PDF (base `pdf()` + grid, no pandoc/LaTeX).
