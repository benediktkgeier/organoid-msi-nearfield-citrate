# Joy_Tables — Prism-ready data for three report figures

Focused hand-off from `figures/annotation/apical_citrate_dha_report.pdf`. Two files, same
data:

| File | Contents |
|---|---|
| `Joy_Tables.xlsx` | 9 tabs — the 3 figures' data + `README` + `Stats_reference` |
| `Joy_Tables.pzfx` | native GraphPad Prism project, 7 ready-to-graph data tables |

Values are read verbatim from the analysis outputs — nothing recomputed.

---

## The three figures

- **Page 1 — "Per-organoid mean ion intensity by apical orientation"**
  → tabs `P1_Citrate_abs`, `P1_DHA_abs`
- **Page 8 — "Absolute near-field citrate level in the gel, by apical class"**
  → tabs `P8_near50`, `P8_near100`
- **Page 8b — "Near-field citrate: 0-50 µm vs 0-100 µm WITHIN each apical class"**
  → tabs `P8b_basolateral_out`, `P8b_apical_out`, `P8b_mixed`

---

## The data

MALDI-MSI (10 µm pixel) of intestinal organoids in hydrogel. Each organoid is classified
by **apical orientation**, and **each value is one organoid** (per-pixel means aggregated
per organoid — organoids, not pixels, are the plotted points).

| Class key | Column label | Meaning | n |
|---|---|---|---|
| `basolateral_out` | Basolateral-out | apical surface faces inward, basolateral faces gel | 23 |
| `apical_out` | Apical-out | apical surface faces outward, toward the gel | 30 |
| `mixed` | Mixed | mixed / indeterminate | 31 |

- **Citrate** = `191.01976` *m/z* `[M–H]⁻` (anchored, ±7 ppm). **DHA** (C22:6) = `327.2330`
  *m/z* `[M–H]⁻`.
- **Units** = TIC-normalized mean intensity (a.u.).
- **Near-field** = mean citrate in the gel within the given distance band outward from the
  organoid surface (`near50` = 0–50 µm, `near100` = 0–100 µm).

---

## Tab reference

| Tab | Figure | Columns | Value |
|---|---|---|---|
| `P1_Citrate_abs` | pg1 (left) | Basolateral-out / Apical-out / Mixed | per-organoid mean citrate |
| `P1_DHA_abs` | pg1 (right) | Basolateral-out / Apical-out / Mixed | per-organoid mean DHA |
| `P8_near50` | pg8 (left) | Basolateral-out / Apical-out / Mixed | near-field citrate 0–50 µm |
| `P8_near100` | pg8 (right) | Basolateral-out / Apical-out / Mixed | near-field citrate 0–100 µm |
| `P8b_basolateral_out` | pg8b | `0-50 um`, `0-100 um` (+ `id` in xlsx) | paired, one row per organoid |
| `P8b_apical_out` | pg8b | `0-50 um`, `0-100 um` (+ `id` in xlsx) | paired, one row per organoid |
| `P8b_mixed` | pg8b | `0-50 um`, `0-100 um` (+ `id` in xlsx) | paired, one row per organoid |

**Blank cells** in P1/P8 are padding — the class columns have unequal length (23/30/31),
not zeros. In P8b the two value columns are **paired by row** (same organoid); the `id`
column (`sid | instance`, xlsx only) is the row label.

`Stats_reference` (xlsx) lists the report's p-values and which Prism test reproduces each.
Stars: `*` p<0.05, `**` p<0.01, `***` p<0.001, `ns`.

---

## Load into Prism

**P1 & P8 (grouped scatter + box):** New table → **Column** (scatter with bar/box). Paste
the 3 class columns. Analyze → nonparametric → **Mann–Whitney** per pair, or
**Kruskal–Wallis + Dunn's** across all three.

**P8b (paired):** New table → **Column**, graph *before-after*. Paste `0-50 um` /
`0-100 um` (id → row titles). Analyze → **Wilcoxon matched-pairs signed-rank**.

Or just open `Joy_Tables.pzfx` — all seven tables are already built.

> p-values are descriptive (organoids within a section are pseudo-replicates); reproduce the
> same tests in Prism for a like-for-like figure.

---

*Regenerate:* `Rscript R/11_per_organoid_final/12_joy_tables.R` (needs R 4.4.2 + `writexl`,
`pzfx`). Sources: `apical_citrate_dha_per_organoid_normalized.csv`,
`apical_gradient_per_organoid.csv`, `apical_citrate_dha_stats.csv`.