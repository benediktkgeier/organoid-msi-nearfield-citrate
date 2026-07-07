# Phase 00 — Shared libraries & locked config

Common libraries, locked constants, and helpers sourced by every phase. Not a
run-order phase; nothing here is executed standalone — the other phases
`source()` these files.

## Scripts / libraries
| File | Role |
|---|---|
| `lib_paths.R` | Locked constants + `cache_in()`/`CACHE_SRC` (read-only upstream cache reuse), `disp_id()` (strips 0h/20h tokens from titles), `APICAL_COLS`, `register_parallel()` (SnowParam SOCK, 4 workers), citrate anchor. Sourced by ALL scripts. |
| `lib_citrate.R` | THE citrate definition: `citrate_onto_pd()` reads a ±`CITRATE_WIN_PPM` (7 ppm) window around `CITRATE_ANCHOR_MZ` (191.01976) per-pixel from RAW centroid imzML; bypasses the merged 25-ppm grid feature. Houses the imzML reader. |
| `gradient_config.R` | Gradient-pipeline config: `GRAD_SIDS`, ring levels, `MIN_INSTANCE_PX`. Single-condition, pooled (no group split, no Δρ). |
| `if_config.R` | IF registration registry + pixel scales (overview/hi-res/BF µm/px), cache/figure/results dirs. Builds on `gradient_config.R`. |
| `lib_register.R` | MSI↔brightfield helpers (`.mis` parser, affine fit/apply/compose, phase-corr). |
| `lib_register_if.R` | IF `.nd2`→BF/MSI helpers (block-mean channel reader, affine chain). |
| `lib_gradient_seg.R` | Reusable segmentation + signed-distance ring builder (`gseg_segment`). |
| `lib_report_frame.R` | Shared base-graphics FRAME primitives (`PMAR`, `frame_box()`, `scalebar_bottom()`, `colorbar_img()`, `SCALE_UM`) for the non-locked cache-only report drivers (`04/05_ssc_report.R`, `11/10`, `11/13`). |
| `lib_gradient_report.R` | Shared per-dataset PANEL renderers (`sec_data` + the 6 panels: MSI ion / native BF / SSC mask / IF / overlay / weather heatmap) for the final-report family `11/10_citrate_gradient_report_final.R` + `11/13_citrate_gradient_report_3class.R`. Requires the caller to have set the render globals (see file header). |
| `lib_nearfield_viz.R` | **LOCKED** four-view per-organoid near-field emission figure toolkit. |
| `lib_nearfield_viz_whole.R` | **LOCKED** whole-region variants of the four-view figure (sources `lib_nearfield_viz.R` unchanged). |

## Notes / gotchas
- **R 4.4.2 only**, Cardinal v3; **never** `setCardinalBPPARAM()`.
- `Analysis_R_Final` is a cache-reuse publication fork: upstream `Analysis_R/cache`
  is read READ-ONLY via `cache_in()`; recomputation writes to this fork's `cache/`.
- Single-condition study — no 0h-vs-20h comparison; `disp_id()` enforces this in titles.
- Apical colours locked in `APICAL_COLS`: apical-out `#C2399A`, basolateral-out `#2C9E4B`, mixed `#888888`.

See also: [`../../PIPELINE.md`](../../PIPELINE.md), [`../../docs/DESIGN_DECISIONS.md`](../../docs/DESIGN_DECISIONS.md).
