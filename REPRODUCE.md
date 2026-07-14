# REPRODUCE.md — reproducing the MSI / microscopy pipeline

Audience: a researcher who has read the paper (methods + results + data are provided separately)
and wants to reproduce the **imaging pipeline** — mass-spectrometry imaging (MSI) preprocessing →
on-tissue delineation → organoid segmentation → MSI↔brightfield↔immunofluorescence registration →
near-field citrate-gradient analysis and the publication figures.

This document is the entry point. The authoritative script-by-script map is
[`PIPELINE.md`](PIPELINE.md); exact package versions are in [`docs/ENVIRONMENT.md`](docs/ENVIRONMENT.md);
each phase folder `R/NN_*/` has its own `README.md`; locked methodological choices are in
[`docs/DESIGN_DECISIONS.md`](docs/DESIGN_DECISIONS.md).

---

## Data availability

This repository holds **code + documentation + small committed inputs** (the sample registry,
the manual-gate reference choices, and result tables). The heavy artifacts — the intermediate
**cache** (the Mode-A dependency: `peaks_*`, `citrate_anchored_*`, `instances_*`, `zones_*`,
`register*/`, `imzml/`), the **raw** imzML/`.nd2`, and the rendered **figures** — are archived in a
separate data deposit and are **`.gitignore`d** from the repo.

> **🔬 MSI datasets (imzML) → METASPACE:**
> **https://metaspace2020.org/project/McKenna-2026?tab=datasets** — the section datasets (plus the
> citrate-standard dilution series) with their METASPACE annotations. Download the imzML/`.ibd` here for
> the raw MSI layer (`imzml/` in the §2 layout).
>
> **🔬 Microscopy (`.nd2`) + Bruker MSI raw (`.d`) → figshare:**
> **https://doi.org/10.6084/m9.figshare.32979014** — the native Nikon `.nd2` brightfield/IF microscopy
> and the archival Bruker `.d`/TDF acquisitions. Place the `.nd2` (and `_over.nd2`) under `<MSI_ROOT>/`
> per the §2 layout. (The runnable MSI input is the METASPACE imzML above, **not** the Bruker `.d`.)
>
> **📦 Intermediate cache + rendered figures: _link to be added_** — the validated `cache/` tree (the
> Mode-A dependency) and the figure PDFs. Place the `cache/` tree at `<MSI_ROOT>/Analysis_R/cache`
> (see §2–§3). Not yet deposited; **Mode B** can rebuild the cache from the raw imzML above.

*(MSI imzML on METASPACE and microscopy/Bruker-raw on figshare are available now; the cache/figures
deposit placeholder will be filled with a DOI once published.)*

---

## 0. What this reproduces, and the two modes

The study is **single-condition** (all organoid sections incubated in CMC ≥5 min, treated
identically — there is *no* 0 h vs 20 h contrast; the tokens in file names are section labels only).
Headline result: near-field citrate emission is higher around **apical-out** than **basolateral-out**
organoids within ≤50 µm (`grad_near50`, p ≈ 6e-5, Cliff's δ).

There are two reproduction modes. **Read this before choosing.**

| Mode | What it does | Needs | Time |
|---|---|---|---|
| **A. Figure-layer (recommended)** | Regenerate every figure/result PDF + CSV from the validated intermediate **cache** | R env + Java + the `cache/` data package | ~15–40 min |
| **B. From-raw rebuild** | Re-run preprocessing → segmentation from the raw imzML/`.nd2` | A + raw imzML/`.nd2` + **the committed manual-gate artifacts** | hours; **partial** (see §6) |

Three steps in Mode B are **human-in-the-loop** (on-tissue ion selection, organoid ROI curation,
apical-orientation annotation). Each reproducer performs these **independently** — that is part of
reproducing the method. We provide **our exact choices** as committed reference artifacts (§6) so you
can either (a) reuse them to reproduce *our* figures bit-for-bit, or (b) make your own calls and
compare. Mode A (which reuses our choices via the cache) is how the paper's figures were produced and
is the faithful reproduction target.

---

## 1. Prerequisites

1. **R 4.4.2** and the **Bioconductor 3.20** package set — install exactly as in
   [`docs/ENVIRONMENT.md`](docs/ENVIRONMENT.md). (`Cardinal 3.8.3`, `matter 2.8.0`, `EBImage 4.48.0`,
   `RBioFormats 1.6.0`, `BiocParallel 1.40.2`, + CRAN `RANN/png/jpeg/viridisLite/writexl/readxl/pzfx`.)
2. **Java 8 JRE** (Bio-Formats backend for `RBioFormats`, required by every `.nd2` step). Set
   `JAVA_HOME` to it (§3.3).
3. **Disk**: ~2 GB for the cache data package (some intermediate `.rds` are ~800 MB each).
4. **The data package** — see §2.

Verify the toolchain before starting (R + every package the pipeline loads):
```bash
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" -e 'cat(R.version.string,"\n"); for(p in c("Cardinal","matter","EBImage","RBioFormats","BiocParallel","RANN","png","jpeg","viridisLite","writexl","readxl","pzfx")) cat(p, as.character(packageVersion(p)), "\n")'
```
Then verify Java 8 is visible (needed by every `.nd2` step — Phases 03/05/06/12; a missing JRE
otherwise only surfaces at the first `.nd2` read):
```bash
echo "$JAVA_HOME"        # should point at a Java 8 JRE
java -version            # should report 1.8.x
```

---

## 2. Data package & layout

The pipeline reads raw acquisitions and a set of intermediate caches. Provide them in this layout
(the reference root is `D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad/`; see §3 to relocate):

```
01062026_JoyMetabolGrad/
├─ Analysis_R_Final/            # this repository
│  ├─ inventory.csv             # 20-row sample registry (sample_id, group, imzml_path, bf_jpg_path, nd2_path, notes)
│  ├─ cache/                    # fork-local outputs (starts EMPTY; writes land here)
│  ├─ figures/  results/        # generated PDFs / CSVs
│  └─ R/ docs/ *.sh
├─ Analysis_R/cache/            # << THE DATA PACKAGE (see below) — read-only upstream cache
│  ├─ imzml/                    # MSI imzML + .ibd, all 20 sections (the effective MSI raw)
│  ├─ peaks_combined.rds, peaks_curated.rds, peaks_tissue_combined.rds, citrate_anchored_<sid>.rds
│  ├─ instances_final_<sid>.rds, zones_<sid>.rds
│  └─ register/  register_if/   # MSI↔BF and IF↔BF transforms (incl. nd2final_<sid>.rds)
├─ 06102026_AO_0h_sl6A.nd2, 06102026_AO_20h_sl4A.nd2   # brightfield .nd2 (project root)
├─ MSI/                         # slide QC JPGs + citrate-standard imzml (sl7A)
│  └─ 06102026_AO_0h_sl7A/imzml/   # standard dilution series imzML
└─ LM_Bsections/               # IF hi-res .nd2 (B-sections): 06172026_AO_0h_sl6b_sec{1..6}.nd2, 06172026_AO_20h_sl4_sec{1..6}.nd2, *_over.nd2
```

**Data provenance / cleanliness:**
- **MSI**: the pipeline consumes **already-converted centroid imzML** (SCiLS-Lab re-export of the
  Bruker `.d`/TDF acquisitions — TIMSCONVERT produced broken files, so SCiLS re-export is the clean
  source; see the `notes` column of `inventory.csv`). The Bruker `.d` raw is **not** used by any
  script; the runnable input is the imzML+ibd. These imzML datasets are published on **METASPACE**
  (<https://metaspace2020.org/project/McKenna-2026?tab=datasets>), and the archival Bruker `.d`/TDF raw
  is on **figshare** (<https://doi.org/10.6084/m9.figshare.32979014>).
- **Microscopy**: native **Nikon `.nd2`** (NIS-Elements) for both brightfield and IF, read directly
  via `RBioFormats`/Bio-Formats (`read.image`/`read.metadata`, series 1) — no conversion needed.
  IF channel order and per-slide display LUTs are locked in `R/00_lib/if_config.R`
  (`ch1=Cy5/F-actin, ch2=mCherry/ZO-1, ch3=FITC/β-catenin, ch4=DAPI`). The `.nd2` files are on
  **figshare**: <https://doi.org/10.6084/m9.figshare.32979014>.
- **Citrate standard**: a Sodium-Citrate dilution series imzML (matrix-matched, same slide sl7A).
  The decade concentrations were mislabeled at acquisition and physically renamed to the TRUE
  concentration (`MSI/.../imzml/RENAME_LOG.txt`; `R/02_CitrateStandard/00_config.R` maps each file).

**The cache is the real data dependency.** `run_all.sh` (Mode A) explicitly does *not* recompute
preprocessing/SSC/segmentation — it reads the caches above through `cache_in()`, which resolves to
`Analysis_R/cache` (read-only) when the fork-local `cache/` is empty. Ship/obtain that tree.

---

## 3. Configure paths — one environment variable

**There is a single switch: `MSI_ROOT`.** Every absolute path in the pipeline is anchored to one
**data root** — the parent folder that holds `Analysis_R_Final/`, `Analysis_R/cache/`, `MSI/`,
`LM_Bsections/`, and the `*.nd2` files (reference value:
`D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad`). Every script bootstraps via
`Sys.getenv("MSI_ROOT", "<reference default>")`, `R/00_lib/lib_paths.R` derives `PROJECT_ROOT` /
`CACHE_SRC` / `MSI_DIR` / `IF_DIR` / … from it, and `load_inventory()` rewrites the paths in
`inventory.csv` through the same switch. So to run on another machine you keep the layout (repo at
`<data root>/Analysis_R_Final`) and set **one** variable — **no source-code edits**:

```bash
# bash: for this session
export MSI_ROOT="/data/JoyMetabolGrad"
# ...or persist it for R in ~/.Renviron:
echo 'MSI_ROOT=/data/JoyMetabolGrad' >> ~/.Renviron
```
Leave `MSI_ROOT` unset to use the reference-machine default. (Verify: `Rscript -e
'source("R/00_lib/lib_paths.R"); cat(PROJECT_ROOT, CACHE_SRC)'` should print your root.)

A few paths are **outside** `MSI_ROOT`, each with its own env override:
- **Java** (`JAVA_HOME`) — see §3.3. Set it to your Java 8 JRE; the reference default is only
  applied when `JAVA_HOME` is unset.
- **Python** for the Phase-09/10 parse steps: `ANNOT_PY_BIN <- Sys.getenv("MSI_PYTHON", "python")`.
  Set `MSI_PYTHON` to your interpreter if `python` on `PATH` is not the right one. Mode-B only —
  not needed for Mode A or Phase 11.
- **`ION_LIB`** — an optional external QC helper (`file.exists`-guarded; falls back to a built-in
  renderer if absent). Edit the literal in `R/01_preprocess/01_preprocess.R` if you have it.

### 3.3 Java
`R/00_lib/lib_paths.R` sets `JAVA_HOME=C:/Program Files/Java/jre1.8.0_491` only if that folder exists
and `JAVA_HOME` is unset. On any other machine, set `JAVA_HOME` in your environment (or `.Renviron`)
to your Java 8 JRE **before** running the `.nd2` steps; ~10 Phase-06 IF scripts also set the
reference path explicitly — edit those if yours differs.

---

## 4. Mode A — reproduce the figure/result set (recommended)

With the environment installed (§1), paths configured (§3), and the cache present (§2):

```bash
cd Analysis_R_Final
./run_all.sh          # regenerates the publication figure/result set from cache
./run_all.sh qc       # also regenerates the upstream QC (citrate-standard gate, metabolite ID)
```
`run_all.sh` runs the analysis/figure layer in dependency order, logging each script to
`results/_regen_logs/<name>.log` (a failing script is logged and skipped; the run continues), and
prints a pass/fail summary to `results/_regen_logs/summary.txt`. The full set it runs (12 steps):
the SSC clustering+segmentation report, the pooled gradient survey (stats + report), the apical
near-field chain (headline stats + locked per-organoid figure + Wilcoxon test + whole-region figure),
the per-dataset v3 gradient, the **final combined per-dataset report**, a second whole-region
invocation for the GLOBAL-heatmap variant (`08_nearfield_wholeregion.R all globalheat`), the IF pairs
report, and the **methods/pipeline report** (`figures/methods_report/generate_methods_report.R` →
`figures/methods_report/methods_pipeline_report.pdf`; note this generator lives under `figures/`, not
in the `R/` phase tree). A clean run reports `TOTAL: 14 passed, 0 failed`.

**Outputs may already exist** (a prior copy of every figure/CSV ships in `figures/`/`results/`).
`run_all.sh` **overwrites** them in place — so a green run is confirmed by the `summary.txt` PASS
list and by fresh file mtimes, not by files merely being present. To be sure you're looking at your
own run, note the start time (or clear the targets first).

Report-only regeneration (no QC recompute) is also available via `./regen_reports.sh`.

Individual headline artifacts:
```bash
RS="/c/Program Files/R/R-4.4.2/bin/Rscript.exe"
$RS R/11_per_organoid_final/03_apical_report.R              # headline stats + figure
$RS R/11_per_organoid_final/08_nearfield_wholeregion.R all globalheat
$RS R/11_per_organoid_final/10_citrate_gradient_report_final.R all   # or `test` for a 1-dataset preview
$RS R/04_ssc_ontissue/05_ssc_report.R all                  # or `test`
```
(The two `…report_final.R` and `05_ssc_report.R` scripts accept `all` | `test` | `<sid>`; `test`
renders a single representative dataset to a `…_TEST_<sid>.pdf` — use it to sanity-check quickly.)

---

## 5. Mode B — rebuild from raw (phase order)

Run in folder order `01 → 11`; within a phase, file-number order. Full I/O per script is in
[`PIPELINE.md`](PIPELINE.md) and each phase `README.md`. Phases that hit a **manual gate** are
marked ⚠ — you cannot recompute the human markup; reuse the committed artifact (§6).

| Phase | Produces | Manual gate? |
|---|---|---|
| **01 preprocess** | `cache/peaks_combined.rds` (Cardinal v3 idiom; locked) | |
| **02 CitrateStandard** | validates + locks citrate; `cache/citrate_anchored_<sid>.rds` | |
| **03 coarse registration** | MSI→slide-JPG→`.nd2` transforms (`cache/register/`) | |
| **04 SSC on-tissue** | curated ions + `cache/peaks_tissue_combined.rds` (`is_tissue`, SSC clusters, chemotype) | ⚠ on-tissue ion selection |
| **05 registration refine** | native BF crops + `cache/register/nd2final_<sid>.rds` | |
| **06 IF registration** | IF↔BF/MSI landmark fits (`cache/register_if/`) → `dataset_pairs_hq.pdf` | ⚠ IF landmarks (manual) |
| **07 metabolite ID** | matched metabolites + citrate QC | |
| **08 organoid survey** | `cache/instances_<sid>.rds`, `zones_<sid>.rds`, pooled gradient | |
| **09 organoid refinement** | `cache/instances_final_<sid>.rds` | ⚠ ROI curation (green-ink cuts) |
| **10 apical annotation** | `results/annotation/apical_map_consensus.csv` | ⚠ two-annotator markup |
| **11 per-organoid final** | headline stats + all publication figures | |

Preprocessing (Phase 01) writes ~800 MB `.rds`; it uses `register_parallel()` (SnowParam SOCK, 4
workers). **Never** `setCardinalBPPARAM()`.

---

## 6. Manual gates — you perform these; our reference choices are provided

Three steps require human judgement and are **meant to be performed by each reproducer**. We ship
**our exact choices** as committed reference artifacts: reuse them to reproduce our figures exactly,
or substitute your own markup to test robustness. The reference artifacts are:

- **On-tissue ion selection** (Phase 04) → `results/peakme_annotations/peakme_*_annotations.csv`
  (curated on/off-tissue ion labels). Consumed by `02_import_ontissue_ions.R`.
- **Organoid ROI curation** (Phase 09) → `results/annotation/organoid_split_lines_raw.tsv`,
  `organoid_remove.csv`, `organoid_island_actions.csv` (from green-ink cut PDFs).
- **Apical-orientation consensus** (Phases 10/11) → **`results/annotation/apical_map_consensus.csv`**
  (two-annotator consensus; the SOLE annotation basis for every Phase-11 figure/stat — no marked-PDF
  re-parse and therefore no Python at Phase-11 runtime).

The parse steps that originally produced these (Phases 09/10) call a Python interpreter
(`ANNOT_PY_BIN <- Sys.getenv("MSI_PYTHON", "python")`, helpers in `py/`; set `MSI_PYTHON` to your
interpreter). You do **not** need Python for Mode A or for Phase 11 — only if you re-derive the
curation TSVs from fresh marked PDFs.

---

## 7. Validation — did it reproduce?

- **Headline stat**: `results/annotation/apical_citrate_dha_stats.csv` → `grad_near50`,
  `basolateral_out` vs `apical_out`, **p ≈ 6.07e-5** (Wilcoxon), apical-out higher (mean 1.70e-5 vs
  1.05e-5; **Cliff's δ ≈ +0.62**; n = 23 vs 30). Page 1 of
  `figures/annotation/apical_citrate_dha_report.pdf` is the Summary Statement.
- **Whole-region figure**: `figures/annotation/apical_nearfield_emission_figure_all_globalheat.pdf`
  (20 pages, outlines coloured by apical class; view-3 heatmap on a global p99.9 scale).
- **Final combined report**: `figures/gradient/citrate_gradient_report_final.pdf` (21 pages).
- **SSC + segmentation**: `figures/ssc/ssc_clustering_segmentation_report.pdf` (front page + 20
  per-dataset pages; tissueness → k4 clusters → on-tissue mask → organoid instances).
- **Run summary**: `results/_regen_logs/summary.txt` should read `TOTAL: 14 passed, 0 failed`.
- **Sanity (optional)**: every script should parse and the runner should be shell-clean:
  ```bash
  "/c/Program Files/R/R-4.4.2/bin/Rscript.exe" -e 'ok<-TRUE; for(f in list.files("R",pattern="[.]R$",recursive=TRUE,full.names=TRUE)) tryCatch(parse(f),error=function(e){cat("PARSE FAIL:",f,"\n");ok<<-FALSE}); cat(if(ok)"all parse OK\n" else "PARSE ERRORS\n")'
  bash -n run_all.sh
  ```
- **Expected benign warning**: logs show `package 'viridisLite' was built under R version 4.4.3` on
  most scripts — harmless (the pinned `viridisLite 0.4.3` binary was simply built on a later R; it is
  *not* a version-pin violation and does not affect output).

---

## 8. Reproducibility assessment (state of the code)

Honest status against the usual checklist, so you know what to expect:

- **Data cleanliness** — good. MSI = clean SCiLS re-exported centroid imzML (the broken TIMSCONVERT
  output was discarded); microscopy = native `.nd2`; the citrate-standard mislabeling was found and
  corrected on disk with a logged rename. The runnable inputs are imzML + `.nd2`, not Bruker `.d`.
- **Instructions** — the pipeline is a clear, linear, phase-numbered analysis (01→11) with a
  top-level `README.md`, an authoritative `PIPELINE.md`, per-phase `README.md`, per-topic `docs/`, and
  three runner scripts. This file is the added reproduction entry point.
- **Packages/versions** — now pinned in [`docs/ENVIRONMENT.md`](docs/ENVIRONMENT.md) (R 4.4.2,
  Bioconductor 3.20, exact package versions). There is **no `renv.lock`** — pin via the Bioconductor
  release and the version table.
- **Scripts** — all parse cleanly and the libraries source without error; the report drivers were
  re-rendered and verified. Locked figure scripts (`08`/`09`, `lib_nearfield_viz*`) are frozen visual
  specs — do not edit; the newer report drivers reuse them and a shared `lib_report_frame.R`.
- **Known frictions** (call-outs, not blockers): (1) paths relocate via the single **`MSI_ROOT`**
  environment variable (§3) — keep the layout, set one variable, no code edits; (2) the intermediate
  **cache is the real data dependency** for the faithful (Mode A) reproduction; (3) three **manual
  gates** are user-performed —
  our reference choices are provided as committed artifacts (reuse for exact figures, or substitute your
  own); (4) Java 8 path is Windows-specific; (5) the external
  `lib_ion_image.R` and Python parse helpers are optional (guarded / Mode-B-only).

---

## 9. Quick start (Mode A, reference machine layout)

```bash
cd /d/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad/Analysis_R_Final
# 1. verify env
"/c/Program Files/R/R-4.4.2/bin/Rscript.exe" -e 'cat(R.version.string)'
# 2. run the figure/result layer from cache
./run_all.sh
# 3. check the summary
cat results/_regen_logs/summary.txt
```