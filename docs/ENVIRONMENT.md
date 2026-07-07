# Environment & package versions (reproduction pin)

The pipeline has **no `renv.lock`/`DESCRIPTION`**; this file is the authoritative version
pin. It records the exact environment the results were produced in. Match it as closely as
possible — the Bioconductor **release** (3.20) pins the Bioc package versions, which is the
part most likely to affect numerical output (Cardinal/matter peak processing).

## Core

| Component | Version (used) | Notes |
|---|---|---|
| R | **4.4.2** (2024-10-31 "Pile of Leaves", ucrt) | NOT 4.5.x. `C:/Program Files/R/R-4.4.2/bin/Rscript.exe` on the reference machine |
| Platform | x86_64-w64-mingw32 (Windows 11) | Bio-Formats path config is Windows-specific (see Java) |
| Bioconductor | **3.20** | pins the Bioc packages below |
| Java | **JRE 1.8.0_491 (Java 8)** | required by `RBioFormats` (Bio-Formats) for all `.nd2` reads |

## Packages

**Bioconductor (release 3.20):**

| Package | Version | Role |
|---|---|---|
| Cardinal | **3.8.3** | MSI data model, preprocessing, peak alignment, SSC clustering |
| matter | **2.8.0** | out-of-memory backend for Cardinal |
| EBImage | **4.48.0** | morphology (segmentation outlines, gblur, connected components) |
| RBioFormats | **1.6.0** | reads Nikon `.nd2` brightfield + IF (wraps Java Bio-Formats) |
| BiocParallel | **1.40.2** | `SnowParam` SOCK parallelism (`register_parallel()`) |

**CRAN:**

| Package | Version | Role |
|---|---|---|
| RANN | **2.6.2** | nearest-neighbour signed-distance rings (`nn2`) |
| png | **0.1.8** | read/write crop PNGs |
| jpeg | **0.1.11** | read slide QC JPGs |
| viridisLite | **0.4.3** | ion-image / heatmap palettes |
| writexl | **1.5.4** | Prism/xlsx export (primary writer) |
| readxl | **1.4.5** | reads the published-metabolite supplement xlsx |
| pzfx | **0.3.1** | GraphPad Prism `.pzfx` export |
| openxlsx | *optional* (not installed on ref machine) | fallback xlsx writer only; `writexl` is used when present |
| pdftools | *optional* | listed in README; not required by any pipeline script (handy for rasterizing PDFs to inspect) |

Base R (bundled): `stats`, `grDevices`, `graphics`, `utils`, `parallel`.

> Most run logs print `package 'viridisLite' was built under R version 4.4.3`. This is **expected and
> harmless** — the pinned `viridisLite 0.4.3` binary happens to have been built on a later R; it is not
> a version mismatch and does not affect output. (Reproduced cleanly on R 4.4.2.)

## Install

```r
# 1. Bioconductor release 3.20 (pins Cardinal/matter/EBImage/RBioFormats/BiocParallel)
install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install(c("Cardinal", "matter", "EBImage", "RBioFormats", "BiocParallel"))

# 2. CRAN
install.packages(c("RANN", "png", "jpeg", "viridisLite", "writexl", "readxl", "pzfx"))
# optional: install.packages(c("openxlsx", "pdftools"))
```

Java 8 (for `RBioFormats`): install a Java 8 JRE and point `JAVA_HOME` at it. On the reference
machine that is `C:/Program Files/Java/jre1.8.0_491`; `R/00_lib/lib_paths.R` sets `JAVA_HOME`
to that path **only if it exists and is unset** — on any other machine set `JAVA_HOME` in your
shell/`.Renviron` before running the `.nd2` steps (Phases 03/05/06), or edit that path.

## Locked runtime conventions (do not change)

- **Never call `setCardinalBPPARAM()`** — triggers a Cardinal v3 bug. Parallelism is set once via
  `register_parallel()` (SnowParam SOCK, 4 workers) in `R/00_lib/lib_paths.R`.
- Cardinal v3 preprocessing idiom (locked): `readMSIData` → `c()` combine → `normalize(tic)` →
  `process()` → manual single-linkage reference grid → `convertMSImagingArrays2Experiment` →
  `peakAlign` → `summarizeFeatures`.
- Citrate is defined ONCE (locked): `citrate_onto_pd()` reads a ±7 ppm window around 191.01976 from
  the RAW centroid imzML — never the merged 25 ppm grid feature (`R/00_lib/lib_citrate.R`).

Captured on the reference machine with `Rscript` + `packageVersion()`; see the reproduction
protocol in [`../REPRODUCE.md`](../REPRODUCE.md).