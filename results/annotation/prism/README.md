# Apical citrate / DHA report — Prism-ready data

This folder reproduces **every quantitative panel** of
`figures/annotation/apical_citrate_dha_report.pdf` in GraphPad Prism.

| File | What it is |
|---|---|
| `apical_prism_tables.xlsx` | 20 tabs — all data tables + stats + this readme, one tab per Prism table |
| `apical_prism_tables.pzfx` | native Prism project: 14 **Column / paired** tables (pages 1–4, 7, 8, 8b) |
| `apical_gradient_decay.pzfx` | native Prism project: 4 **XY** decay-curve tables (page 6) |
| `apical_gradient_profile_long.csv` | tidy long form of the page-6 data (1 row per organoid × zone) |
| `README.md` | this file |

**No values were recomputed for the tables** — they are read verbatim from the CSVs the
report writes. The page-6 profile (`grad_prof`) is the one series the report computes but
never saved to disk, so `R/11_per_organoid_final/11_export_gradient_profile.R` re-derives
it with the *identical* code path (spot-checked: medians match to 4 dp).

Two `.pzfx` files exist because a single Prism project written by the `pzfx` package cannot
mix **Column** and **XY** table types.

---

## 1. The data

MALDI-MSI (10 µm pixel) of intestinal organoids in hydrogel (Stanford metabolic-gradient
study; 0 h + 20 h pooled). Each organoid is classified by **apical orientation**:

| `apical_class` key | Label | Meaning | n |
|---|---|---|---|
| `basolateral_out` | Basolateral-out | apical surface faces inward, basolateral faces the gel | 23 |
| `apical_out`      | Apical-out      | apical surface faces outward, toward the gel | 30 |
| `mixed`           | Mixed           | mixed / indeterminate | 31 |

**84 annotated organoids. Unit of replication = one organoid** (per-pixel means aggregated
per organoid; organoids, not pixels, are the plotted points).

**Ions:** citrate `191.01976` *m/z* `[M–H]⁻` (anchored, ±7 ppm, raw imzML); DHA (C22:6)
`327.2330` *m/z* `[M–H]⁻`. **Units:** TIC-normalized mean intensity (a.u.) unless a
normalization is named. DHA is a designed negative control (membrane-bound,
organoid-confined) — the citrate/DHA ratio cancels any factor scaling both ions together.

**Biology:** apical-out organoids secrete more citrate into the surrounding gel near-field
(≤50 µm); the report's headline. DHA does not follow, and the outward *slope* is n.s., so
the claim is near-field **level**, not gradient steepness.

---

## 2. Tab → report page map

### Column tables — grouped scatter + box, one column per class (Basolateral-out / Apical-out / Mixed)

| Tab | Page | Value | Meaning |
|---|---|---|---|
| `P1_Citrate_abs`      | 1 | `cit`      | per-organoid mean citrate (absolute TIC) |
| `P1_DHA_abs`          | 1 | `dha`      | per-organoid mean DHA (absolute TIC) |
| `P2_Citrate_metaTIC`  | 2 | `cit_mtic` | citrate / metabolite-TIC pool (% of 341-feature pool) — removes regional ionization bias |
| `P2_DHA_metaTIC`      | 2 | `dha_mtic` | DHA / metabolite-TIC pool (%) |
| `P3_Citrate_secnorm`  | 3 | `cit_rel`  | log₂( organoid / section median ) — removes per-section baseline |
| `P3_DHA_secnorm`      | 3 | `dha_rel`  | same for DHA |
| `P4_CitDHA_ratio`     | 4 | `cit_dha`  | log₂( citrate / DHA ) internal-control ratio |
| `P7_rho_out`          | 7 | `rho_out`  | outward monotonicity (Spearman ρ of zone vs citrate; 0 = flat) |
| `P7_far_index`        | 7 | `far_index`| far-field / surface citrate (≥80 µm) |
| `P8_near50`           | 8 | `near50`   | absolute near-field citrate, gel 0–50 µm |
| `P8_near100`          | 8 | `near100`  | absolute near-field citrate, gel 0–100 µm |

Rows = organoids; columns are unequal length (23/30/31) — **blank cells are padding, not
zeros**.

### Paired tables — page 8b (0–50 µm vs 0–100 µm, same organoid)

| Tab | Rows |
|---|---|
| `P8b_basolateral_out` / `P8b_apical_out` / `P8b_mixed` | organoids of that class |

Columns: `id` (`sid | instance`, a row title), `0-50 um` (= near50), `0-100 um`
(= near100). The two value columns are **paired by row**. Organoids missing either value
are dropped (matches the report).

### XY tables — page 6 outward decay curves

Tabs `P6_surf`, `P6_abs`, `P6_mtic`, `P6_citdha`. Column 1 `distance_um` is the X axis
(7 Voronoi zones: 10/20/50/80/160/250/500 µm; plot on a **log X axis** as the report does).

- **`.xlsx`**: `med_Basolateral-out`, `med_Apical-out` (bold group-median curves), plus
  `med_DHA-control` on `P6_surf` only, then one column per organoid (`<class> | sid inst`)
  = the thin **spaghetti** lines (basolateral-out + apical-out only; mixed omitted, exactly
  as the report plots page 6).
- **`apical_gradient_decay.pzfx`**: the median curves only (bold lines). The per-organoid
  spaghetti lives in the `.xlsx` / long CSV — Prism cannot auto-group 50+ Y columns by
  class, so overlay them manually if you want the spaghetti.

The four metrics: `surf` = citrate / surface-band citrate; `abs` = absolute TIC citrate;
`mtic` = citrate / metabolite-TIC (%); `citdha` = citrate / interior DHA.

### Reference tabs (xlsx only — text, not Prism plots)

- `README` — condensed version of this file.
- `Stats_reference` — the report's pairwise p-values (Mann–Whitney for the class
  comparisons; paired Wilcoxon for page 8b) and which Prism test reproduces each.
  Stars: `*` p<0.05, `**` p<0.01, `***` p<0.001, `ns`.

### Not reproducible in Prism

Pages 9+ (per-section brightfield + MSI ion-overlay images) are spatial **raster images**,
not a Prism graph type — excluded by design.

---

## 3. How to load into Prism

**Column tabs (P1–P4, P7, P8):** New table → **Column** (scatter with bar/box). Paste the
3 class columns; blanks are fine. Analyze → nonparametric → **Mann–Whitney** per pair, or
**Kruskal–Wallis + Dunn's** across all three.

**Paired tabs (P8b):** New table → **Column**, graph *before-after*. Paste `0-50 um` /
`0-100 um` (id = row titles). Analyze → **Wilcoxon matched-pairs signed-rank**.

**XY tabs (P6):** New table → **XY**. Column 1 = X (distance_um); set the X axis to **log**.
Plot `med_*` as connected lines (page-6 bold curves). Or just open
`apical_gradient_decay.pzfx` directly — the XY tables are already built.

> The report treats organoids within a section as pseudo-replicates, so its p-values are
> **descriptive**. Reproduce the same tests in Prism for a like-for-like figure; a fully
> rigorous model would add section as a random effect.

---

## 4. Provenance & regeneration (all inputs read-only)

| Data | File |
|---|---|
| Pages 1–4 | `results/annotation/apical_citrate_dha_per_organoid_normalized.csv` |
| Pages 7, 8, 8b | `results/annotation/apical_gradient_per_organoid.csv` |
| Page 6 | `results/annotation/apical_gradient_profile_long.csv` (from script 11) |
| Statistics | `results/annotation/apical_citrate_dha_stats.csv` |
| Report figure / script | `figures/annotation/apical_citrate_dha_report.pdf` · `R/11_per_organoid_final/03_apical_report.R` |

The **no-suffix** CSVs are used — they match the current PDF (all written Jun 30). Older
`*_consensus` / `*_toppick` variants are earlier annotation runs.

Regenerate everything:
```
Rscript R/11_per_organoid_final/11_export_gradient_profile.R   # page-6 long CSV (loads MSI; ~1-2 min)
Rscript R/11_per_organoid_final/09_prism_export.R              # -> apical_prism_tables.xlsx
Rscript R/11_per_organoid_final/10_pzfx_export.R               # -> *.pzfx (Column + XY)
```
Requires R 4.4.2 + packages `writexl` and `pzfx`.
