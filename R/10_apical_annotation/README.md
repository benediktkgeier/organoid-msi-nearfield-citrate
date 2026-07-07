# Phase 10 — Apical-orientation annotation

Assign each curated organoid a polarity class — **apical-out** vs
**basolateral-out** vs **mixed** — from manual markup, then reconcile the two
annotators into a consensus map keyed by sid+instance.

## Scripts (run order)
| Script | Role |
|---|---|
| `01_apical_annotate.R` | Emit markup PDF for **MANUAL** apical-out / basolateral-out / mixed scoring. |
| `02_apical_parse.R` | Parse the marked PDF (via centroid sidecar) into per-organoid classes. |
| `03_apical_consensus.R` | Reconcile the two annotators → consensus. |
| `04_apical_consensus_finalize.R` | Finalize → `results/annotation/apical_map_consensus.csv`. |

## Inputs
- `cache/instances_final_<sid>.rds` (Phase 09); island-cleanup canvas + centroid sidecar; two annotators' markup.

## Outputs
- `results/annotation/apical_*` including the committed **`apical_map_consensus.csv`** (Phase 11's sole annotation basis).

## Run
```bash
for s in 02_apical_parse 03_apical_consensus 04_apical_consensus_finalize; do
  "/c/Program Files/R/R-4.4.2/bin/Rscript.exe" R/10_apical_annotation/$s.R
done
```

## Notes / gotchas
- Polarity/colours locked (`APICAL_COLS`): apical-out `#C2399A`, basolateral-out `#2C9E4B`, mixed `#888888`.
  ("apical-in" was renamed to "basolateral-out".)
- Apical classes are annotated on `organoid_island_cleanup.pdf`; parsing uses the centroid sidecar.
- The consensus map is a static committed input; Phase 11 reads it directly (no runtime PDF re-parse / no python).

See also: [`../../docs/apical_annotation.md`](../../docs/apical_annotation.md), [`../../PIPELINE.md`](../../PIPELINE.md).
