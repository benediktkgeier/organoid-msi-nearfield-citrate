# Project history / changelog

Consolidated from the former `HANDOFF.md`, `SESSIONS.md`, `NEXT_SESSION_PROMPT.md`, and
`FINALIZE_RENAME.md` (all superseded by this file + `PIPELINE.md` + `docs/DESIGN_DECISIONS.md`).

> **Final-state note.** The study is **single-condition** (all sections incubated in CMC ≥5 min
> before freezing, treated identically; **no 0h-vs-20h / incubation-time comparison** — any
> "0h"/"20h" below or in filenames is only a raw-file / sample-ID token, never a study aim). The
> published fork `Analysis_R_Final` **reuses the upstream caches read-only** and re-runs only the
> analysis/figure layer. The restructure entry below lists "10 phase folders"; the fork adds an
> early **Phase 02 — Citrate standard GATE**, so the current structure is **11 phases**
> (`01_preprocess` … `11_per_organoid_final`) — see `PIPELINE.md`.

## 2026-06-24 — Pipeline restructure into 10 phase folders
Reorganized `R/` from a flat, organically-grown numbering into ten phase-named folders
(`01_preprocess` … `10_per_organoid_final`) plus `R/00_lib/` (shared libs/configs) and
`R/_deprecated/` (provenance). All `source()` paths repointed; figure/results dirs renamed off the
stale `phaseN_` scheme; empty placeholder folders and stray files removed; docs consolidated into
`README.md` + `PIPELINE.md` + `docs/`. "PeakMe" wording removed from the narrative — the on-tissue
ion list is presented as a manual **selected on-tissue ions** curation. Reports regenerated from
`cache/*.rds`. See `PIPELINE.md` for the legacy→new script mapping.

## 2026-06-23 — Near-field distances 40 → 50 µm
Near-field analysis reframed around 50 & 100 µm (was 40 & 100): emission pair `near40`→`near50`,
contour rings `0/40/100`→`0/50/100`, outward ladder `BUF_LEVELS_UM` …/40/… → …/50/…; per-organoid
test band 50–100 µm. Inward zones unchanged; metabolite reports untouched.

## 2026-06-22 — Phase F (IF↔brightfield/MSI registration) DONE
Serial IF B-sections (sl6b/sl4b) registered onto the MSI A-slide brightfield via manual-landmark,
organoid-level alignment (auto failed). Final deliverable
`figures/if_registration/dataset_pairs_hq.pdf` (20 datasets × 3 panels). IF `.nd2` channel order
ch1=Cy5, ch2=ZO-1, ch3=β-cat, ch4=DAPI; pixel size from `dCalibration` already includes the 1.5×
mag; 10/12 IF sections fit; RMSE 28–222 µm. Detail: `docs/if_registration.md`.

## 2026-06-22 — Metabolite naming finalized
The published reference m/z list was formerly labelled "samarah" (paper author); all
files/dirs/scripts/reports renamed to **metabolite(s)**. The featureData columns `samarah_* →
metabolite_*` were patched in place in `peaks_curated.rds` / `peaks_tissue_combined.rds`
(metadata-only — no SSC re-run); gradient CSVs/PDFs regenerated. Verified 0 "samarah" in results.

## 2026-06-21 — Citrate annotation QC
Established that the *m/z* 191 feature is an unresolved blend (+8–10 ppm above true citrate);
true citrate (191.0197) is real and abundance-variable (3.4×, p=0.008) but overlap-limited.
Window recommendation: upper edge ≤ 191.0204. See `docs/DESIGN_DECISIONS.md` + reports in
`figures/metabolites/`.

## 2026-06-17 — 20-dataset expansion + brightfield registration
20-dataset citrate/DHA + gradient survey (72 organoids). All 20 MSI sections registered to native
`.nd2` brightfield (`cache/register/nd2final_<sid>.rds`); report
`figures/registration/registration_native.pdf`. Detail: `docs/registration.md`.

## 2026-06-16 — Re-plan: single 20-section SCiLS pipeline
The earlier DEV/BACKGROUND split was retired (the TIMSCONVERT-based annotations were invalid; all
5-section annotation-derived work was discarded). One pipeline on the 20-section SCiLS MSE
(`cache/peaks_combined.rds`, 8,772 × 90,760, self-contained). Preprocess → curated 348-ion set →
SSC on-tissue mask → gradient pipeline.

## Notes / caveats
- Taurocholate (514.2844) sits +19.28 ppm in the binned grid — bin-edge effect, not real drift.
- `SnowParam(workers=4)` (`N_WORKERS` in `R/00_lib/lib_paths.R`) — conservative for the ~90,760-px run.
- All working MSEs are self-contained in-memory sparse; `cache/imzml/*.ibd` is a raw archive only.
- imzML conversion: 5/21 raw imzMLs were healthy at TIMSCONVERT v2.0.0; the rest hit a deterministic
  converter bug. Not blocking the 20-section SCiLS pipeline (which exported cleanly from SCiLS Lab).
