# Organoid apical-orientation annotation

> **Provenance doc (upstream annotation workflow).** This describes how the apical classes
> were originally produced (manual Adobe markup → parse → two-annotator consensus). In this
> publication **fork** the annotation is NOT re-run: the finalized two-annotator **consensus
> map** `results/annotation/apical_map_consensus.csv` is a static committed input, so the fork
> does **not** re-parse any marked PDF (no python/reticulate at runtime;
> `organoid_island_cleanup.pdf` is kept only as provenance). The three classes are
> **apical-out**, **basolateral-out**, and **mixed** — what was once typed "apical in"/"ai"
> during markup is recorded as **basolateral-out** (there is no "apical-in" class). Legacy flat
> script numbers below (`R/40`, `R/41`, `R/102`, `R/103`, …) decode via the Phase 09/10 tables
> in `PIPELINE.md`.

**Pipeline pre-step.** Run this **right after the MSI↔brightfield alignment**
(Phase 03/05) and **before** the gradient analysis (Phase 08+). It is kept in its own
`figures/annotation` / `results/annotation` folders (separate from
`figures/gradient`) because every run of the pipeline should do it once, early.

Score each segmented organoid as **basolateral-out**, **apical-out**, or **mixed** by
commenting on a per-section brightfield PDF in Adobe. The result is a
machine-readable table keyed by `(sid, instance)` that joins straight onto the
gradient analysis.

> **Do the morphological split FIRST.** Plain connected-component segmentation
> (R/40) merges touching organoids into one instance. Before apical scoring,
> separate them with the green-line split step below (`R/41`), then re-render this
> PDF so each organoid is its own numbered outline. See **[Step 0: morphological
> split](#step-0-morphological-split-of-clustered-organoids)**.

## Step 0: ROI curation (split / merge-delete / finalize)

**Full workflow + conventions: [`README_roi_curation.md`](../../README_roi_curation.md)**
(Phase E1–E3, run immediately after segmentation `R/40`; one-shot runner
`R/45_run_roi_curation.R` = `R/41` split → `R/43` merge/delete → `R/44` finalize →
`R/42` canvas). The summary below covers the green-line **split** part; merges,
deletes, and the "one connected ROI = one id" finalize are in the runbook.

This corrects the segmentation **right after R/40** (mirrors the gastric-gland
split in `D:/R/PeakMe/PeakMe_GCPL/phaseR_apply_splits.R`). It is **iterative** —
draw a round of cuts, apply, re-render, draw more on the updated canvas:

1. Render the annotation PDF (`R/102`, see below) — it is also the drawing canvas.
2. In Adobe, **draw a GREEN freehand line** *through* each pair/cluster of
   organoids that the segmentation wrongly merged (pencil/ink tool, colour green).
   The line must fully cross the blob so it disconnects into ≥2 pieces. Save (you
   can save over `organoid_apical_annotation.pdf` or as `*_BG.pdf` — R/41 reads
   both by default).
3. Apply the cuts:

   ```
   Rscript R/41_organoid_split_apply.R
   ```

   It reads the green `/Ink` lines, maps them onto the MSI grid (exact inverse of
   the R/102 render), lays each as a 1-px barrier, re-runs connected components,
   and splits the instance. **IDs: preserve + append** — untouched organoids keep
   their id; a split keeps its largest piece as the original id, new pieces get
   appended ids. Output → a **separate** file `cache/instances_split_<sid>.rds`
   (originals `cache/instances_<sid>.rds` untouched).

   **Accumulation (multi-round).** Every stroke is mapped to MSI coords and stored
   in `cache/organoid_split_strokes.rds` (deduped by geometry). The **store**, not
   any PDF, is the source of truth, so re-rendering the canvas never loses earlier
   cuts — each run adds only the *new* strokes and re-applies them all. To rebuild
   the cut set from scratch (e.g. to drop a misplaced cut), point R/41 at one clean
   PDF with `--reset`:  `Rscript R/41_organoid_split_apply.R --reset <clean.pdf>`.
4. **Check** `figures/annotation/organoid_split_curation.pdf`: each green cut sits
   on the drawn line, grey = original outline, coloured + numbered = the new
   pieces. `results/annotation/organoid_split_report.csv` lists per-section
   `n_before → n_after`; any stroke that **failed to disconnect** anything is
   listed in `organoid_split_nonsplitting_strokes.csv` (redraw it across the blob).
5. **Re-render** the annotation PDF so it shows the split organoids:

   ```
   Rscript R/102_organoid_apical_annotate.R
   ```

   R/102 automatically prefers `instances_split_<sid>.rds` when present. Draw the
   next round of cuts on this updated PDF, or proceed to apical scoring below.

## Files

| file | purpose | edit? |
|---|---|---|
| `figures/annotation/organoid_apical_annotation.pdf` | cover + 1 page/section: native BF with numbered, coloured organoid outlines | **YES — green split lines, then apical comments** |
| `figures/annotation/organoid_apical_annotation_BG.pdf` | optional copy carrying **green `/Ink` split lines** (Step 0); read by R/41 too | **YES — draw green cuts** |
| `R/41_organoid_split_apply.R` | accumulates + applies the green cuts → `cache/instances_split_<sid>.rds` | — |
| `py/extract_pdf_lines.py` | pulls Adobe ink/line geometry + colour → TSV (pypdf) | — |
| `cache/organoid_split_strokes.rds` | **accumulated cut store** (all rounds, deduped); source of truth | read-only |
| `figures/annotation/organoid_split_curation.pdf` | per-section split QC: green cut + grey (orig) + coloured (new) outlines | read-only |
| `results/annotation/organoid_split_report.csv` | per-section `n_before → n_after`, which ids split | read-only |
| `results/annotation/organoid_split_nonsplitting_strokes.csv` | strokes that failed to disconnect (redraw these); only when present | read-only |
| `cache/organoid_apical_label_positions.rds` | sidecar: each organoid's PDF-point position (for mapping comments back) | read-only |
| `R/102_organoid_apical_annotate.R` | renders the PDF + sidecar | — |
| `py/extract_pdf_annots.py` | pulls Adobe comments → TSV (pypdf) | — |
| `R/103_organoid_apical_parse.R` | edited PDF → keyed annotation table | — |
| `results/annotation/organoid_apical_annotations.csv` / `.rds` | **output**: `sid, instance, group, apical_class, note, …` | read-only |
| `results/annotation/organoid_apical_parse_log.csv` | audit of every comment (accepted / skipped / why) | read-only |

Outline colours and ID numbers match `citrate_gradient_perdataset.pdf`, so you
can cross-reference the same organoids between the two reports.

## How to annotate (Adobe Acrobat / Reader)

1. Open `organoid_apical_annotation.pdf`. Page 1 is instructions; pages 2–21 are
   the 20 sections.
2. For each numbered organoid, add a **text / sticky-note comment placed ON or
   NEXT TO it**. The comment text is one of:

   | type to record | meaning | also accepted |
   |---|---|---|
   | `in` | basolateral-out | `apical in`, `ai` |
   | `out` | apical-out | `apical out`, `ao` |
   | `mixed` | mixed | `mix` |

3. **Safest:** prefix the organoid number, e.g. `3 out` or `7: mixed`. A
   number-prefixed comment is matched by ID (unambiguous). An un-prefixed comment
   is matched to the **nearest** organoid on the page (≤ 25 pt), so place it close
   to the number.
4. **Save** the PDF.
5. Tell Claude **"apical annotation done"** (or run it yourself):

   ```
   Rscript R/103_organoid_apical_parse.R figures/annotation/organoid_apical_annotation.pdf
   ```

The parser reports counts per class, flags any organoid with conflicting/duplicate
comments, and notes any organoid not yet annotated. Nothing is silently dropped —
unmatched or unparseable comments land in `organoid_apical_parse_log.csv`. You can
re-edit and re-run as many times as you like (re-running overwrites the outputs).

## Re-rendering

```
Rscript R/102_organoid_apical_annotate.R                 # all 20 sections
Rscript R/102_organoid_apical_annotate.R AO_0h_sl6A_sec1a # single-page TEST
```

Re-rendering regenerates the sidecar; if organoid IDs change you must re-annotate.

## Feeding back into the analysis

`organoid_apical_annotations.csv` is keyed by `sid` + `instance`. The per-organoid
gradient outputs (`zones_<sid>.rds` `instance`; the per-organoid lines in `R/87`;
`R/88_gradient_test_perorganoid.R`) use the same `instance` IDs, so downstream
scripts can `merge(..., by = c("sid", "instance"))` (or join `sample_id`→`sid`) to
stratify citrate emission by `apical_class`.
