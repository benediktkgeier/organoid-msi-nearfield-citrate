# Organoid ROI curation (Phase E) — run immediately after segmentation

> **Path note:** this doc uses the old flat script numbers (`R/40`…`R/45`, `R/102`). In the
> fork ROI curation is **Phase 09** (`R/09_organoid_refinement/`, runner `00_run_roi_curation.R`
> chaining `01`→`03`→`04` with the `02`/`02b` canvas + centroid sidecar); segmentation is
> `R/08_organoid_gradient_survey/01_segment_organoids.R` and the apical re-render is
> `R/10_apical_annotation/01_apical_annotate.R`. See the legacy→new table in `PIPELINE.md`.

**This step runs right after organoid segmentation (`R/40_segment_organoids.R`) and
before any downstream organoid analysis (apical scoring, gradient, etc.).** Plain
connected-component segmentation merges touching organoids, leaves over-split
fragments, spurious off-tissue blobs, and ROIs that share a number but aren't
connected. This phase curates the per-section instance maps into clean ROIs —
one connected ROI = one unique id — through an iterative Adobe-marking loop.

Everything is **layered and non-destructive**: the originals (`instances_<sid>.rds`)
are never modified; each step writes its own layer, and renders always read the
most-processed layer present.

```
orig  (R/40)            cache/instances_<sid>.rds
  └─ split   (R/41)     cache/instances_split_<sid>.rds     green-line cuts
       └─ clean (R/43)  cache/instances_clean_<sid>.rds     merges + deletes
            └─ final (R/44) cache/instances_final_<sid>.rds  1 connected ROI = 1 id; removals
render preference everywhere:  final > clean > split > orig   (R/42 canvas, R/102 apical)
```

## The loop

1. **Render the canvas:** `Rscript R/42_island_cleanup_canvas.R`
   → `figures/annotation/organoid_island_cleanup.pdf` (one page/section, native
   brightfield, each organoid a coloured outline with its number placed just
   outside it; deterministic colours: id1=red, id2=blue, id3=green, …).
2. **Mark in Adobe** (see conventions below), save the PDF.
3. **Apply:** `Rscript R/45_run_roi_curation.R figures/annotation/organoid_island_cleanup.pdf`
   (runs split → clean → finalize → re-render in order).
4. **Check** `organoid_island_cleanup.pdf` (+ the QC PDFs / reports) and **repeat**
   until the ROIs are correct. Then proceed to apical scoring (`R/102`, which
   auto-reads the `final` layer).

## Marking conventions

| Intent | How to mark | Applied by |
|---|---|---|
| **Split** a merged blob | draw a **GREEN** freehand line fully *across* the neck so it disconnects | `R/41` (accumulated in `cache/organoid_split_strokes.rds`) |
| **Delete** / **Merge** ROIs | sticky-note text near the ROI(s); Claude interprets into `results/annotation/organoid_island_actions.csv` (`sid,op,ids`) | `R/43` |
| **Disconnected ROI → own id** | a **RED** cross on the detached piece (a pointer); the general rule in `R/44` promotes *any* disconnected component ≥ `DISCONNECT_MIN_PX` to a new id automatically | `R/44` |
| **Remove a spurious / off-tissue ROI** | flag it; Claude adds `sid,id` to `results/annotation/organoid_remove.csv` | `R/44` |

Notes:
- A split line must **fully cross** the organoid or it won't disconnect (it is
  reported in `results/annotation/organoid_split_nonsplitting_strokes.csv`).
- Splits **accumulate** across rounds (the stroke store is the source of truth, so
  re-rendering never loses cuts). To drop a superseded cut, surgery on the store
  records its key in `cache/organoid_split_dropped.rds` so it is never re-added.
- `R/41 --any-color`: accept a split line drawn in a non-green colour (pass only
  the relevant PDF). `R/41 --apply-only`: re-derive from the store without reading PDFs.
- Free-form notes are interpreted **with the user** into the explicit action CSVs —
  not NLP-parsed — so every applied edit is on record and reproducible.

## Scripts & files

| script | role |
|---|---|
| `R/41_organoid_split_apply.R` | apply/accumulate green-line splits |
| `R/43_island_cleanup_apply.R` | apply merge/delete actions |
| `R/44_finalize_instances.R` | one connected ROI = one id; drop specks; apply removals |
| `R/42_island_cleanup_canvas.R` | render the clean numbered canvas (+ sidecar) |
| `R/45_run_roi_curation.R` | runner: 41 → 43 → 44 → 42 in order |
| `py/extract_pdf_lines.py` | extract ink geometry + colour (pypdf) |
| `py/extract_pdf_annots.py` | extract text comments (pypdf) |

Reports: `organoid_split_report.csv`, `organoid_split_nonsplitting_strokes.csv`,
`organoid_island_cleanup_report.csv`, `organoid_finalize_report.csv` (all in
`results/annotation/`). Action inputs: `organoid_island_actions.csv` (merge/delete),
`organoid_remove.csv` (removals).
