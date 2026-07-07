# Citrate isotopes & adducts report

Self-contained QC report that renders the **full citrate ion family** for the 5
sections with genuine citrate signal, at the project's **normal ±10 ppm**
integration window, and answers two questions:

1. Do the citrate **adducts** (Na, Na₂, Cl, acetate) and the **¹³C isotope** appear
   at the expected masses, and how abundant is each?
2. Is citrate **[M−H]− 191.0197** a genuine ion, or merely the ¹³C (M+1) satellite
   of an upstream monoisotopic ion at **m/z 190.0163** (= 191.0197 − 1.0034)?

Project root: `D:/A_Stanford/HelicobacterMSI/01062026_JoyMetabolGrad/Analysis_R_Final`
R version: **R 4.4.2** (`C:/Program Files/R/R-4.4.2/bin/Rscript.exe`).
Companion to `R/07_metabolite_id/04_citrate_window_images.R` (legacy `R/83`, which instead sweeps
the *special* ±2.5/±3/±5 ppm windows around 191.0197 to study isobar bleed-through).

> **Path note:** this doc uses old flat numbers (`R/82`, `R/83`, `R/83b`). In the fork the
> citrate QC scripts are **Phase 07** (`R/07_metabolite_id/03_citrate_resolution.R`,
> `04_citrate_window_images.R`, `05_citrate_isotopes_adducts.R`). Decode via the legacy→new
> table in `PIPELINE.md`.

## Headline results

**All six citrate-family ions are mass-accurate** (measured intensity-weighted
centroid within ±10 ppm, pooled over the 5 sections):

| Ion | Adduct formula | Theoretical m/z | Measured m/z | Δ ppm | p99.5 (a.u.) |
|---|---|---|---|---|---|
| Citrate [M−H]− (C12) | C₆H₇O₇⁻ | 191.01973 | 191.02016 | +2.29 | 353 |
| ¹³C isotope (C13) | ¹³C¹²C₅H₇O₇⁻ | 192.02308 | 192.02328 | +1.06 | 149 |
| [M−2H+Na]− | C₆H₆O₇Na⁻ | 213.00167 | 213.00151 | −0.76 | 89.6 |
| [M−3H+2Na]− | C₆H₅O₇Na₂⁻ | 234.98361 | 234.98361 | −0.02 | 119 |
| [M+Cl]− | C₆H₈O₇Cl⁻ | 226.99640 | 226.99623 | −0.76 | 104 |
| [M+CH₃COO]− (acetate) | C₈H₁₁O₉⁻ | 251.04086 | 251.04180 | +3.78 | 135 |
| upstream M? (diagnostic) | 190.0163 | 190.01637 | 190.01630 | −0.35 | 173 |

**Isotope cross-check — 191.0197 is NOT the ¹³C of 190.0163.** Two independent lines
of evidence:

- **Intensity ratio** I(191.0197)/I(190.0163) = **2.06**. For 191 to be the ¹³C M+1
  of 190, the precursor would need ~**190 carbons** (ratio × 0.989/0.0107). Citrate
  has 6 → impossible. 191 carries far more signal than any ¹³C tail of 190 could.
- **Spatial Spearman** r(191, 190) = **−0.70** pooled (per-section −0.70, −0.74,
  −0.66, −0.68, −0.67 — all 5 negative). A true isotope satellite is *perfectly
  co-localized* with its monoisotopic peak (r ≈ +1); here the two ions are strongly
  **anti-correlated**, i.e. distinct molecules occupying complementary tissue.

So there *is* a real, well-aligned ion at 190.0163 (comparable abundance), but
citrate [M−H]− 191.0197 is an **independent ion**, not its isotope. The strong
anti-correlation suggests 190.0163 marks a spatially complementary compartment
worth identifying separately.

## Run

```
RS="C:/Program Files/R/R-4.4.2/bin/Rscript.exe"
"$RS" R/83b_citrate_isotopes_adducts.R
```

Single call, no prerequisite steps. ~1 min (reads 5 raw `.ibd` files). Prints the
mass table + isotope-check numbers to the console and writes the PDF below.

## Inputs / outputs

| | path | notes |
|---|---|---|
| script | `R/83b_citrate_isotopes_adducts.R` | self-contained; adapted from `R/83` |
| config | `lib_paths.R` | uses `FIG_DIR`, `CACHE_DIR`, `MSI_PIXEL_UM` (10 µm), `IMG_CLIP_HI` (0.995) |
| data | `cache/imzml/06102026_ao_{0h_sl6a_sec2a, 0h_sl6a_sec5a, 20h_sl4a_2a, 20h_sl4a_3a, 20h_sl4a_3b}.{imzML,ibd}` | the 5 sections with citrate signal; read directly via the low-level binary parser |
| **output** | `figures/metabolites/citrate_isotopes_addocts.pdf` | 2 pages (filename spelling "addocts" kept per request) |

Package dependency: **`viridisLite`** only. (`pdftools`/`png` were used to rasterize
pages for visual QA but are **not** required to produce the report.)

## What the PDF contains

**Page 1 — ion image grid (7 ions × 5 sections + per-row colorbar).**
Render spec (locked v3, identical to `R/83`): **viridis(256), linear, gamma 1.0,
per-ion global p99.5 clip** across the 5 sections, **500 µm** white scale bar,
neutral section box (single condition). Each row's left label cell gives ion name,
adduct formula, theoretical/measured m/z + Δppm, the ±10 ppm window, and the p99.5
clip. The **right-column colorbar** maps summed ion intensity (a.u.) from 0 to that
row's p99.5 clip, so absolute per-ion abundances are directly readable. A purple
footnote states the isotope-check verdict.

**Page 2 — per-section spectral zoom-ins (189.98–192.10).**
Summed ion intensity per **0.0001 Da** bin, **log1p** y-axis (compresses the
dominant 191.046 isobar so the small ¹³C stays visible). Reference lines:
**green dotted = upstream 190.0163, red = citrate C12 191.0197, blue dashed =
¹³C 192.0231**; the 1.0034 Da C12→C13 spacing is annotated. 6th panel carries the
legend + isotope-check verdict box.

## Method notes

- **Masses** derive from the locked citrate value `C12 = 191.019726`. Neutral
  M = C12 + proton = 192.027002; adduct anions add the electron mass.
  ¹³C step = 1.0033548. Constants: H 1.00782503, Na 22.98976928, Cl 34.96885268,
  acetate CH₃COO (C₂H₃O₂) 59.013304.
- **Per-pixel integration** is a summed intensity over each ±10 ppm window, read in
  one pass from the `.ibd` binary; the same pass accumulates the 0.0001-Da spectrum
  histogram and the intensity-weighted m/z (for the measured centroid).
- **Isotope check** uses pooled summed intensities for the ratio and pixel-wise
  Spearman correlation over pixels where either ion > 0.

## Caveats / open items

- **Acetate adduct.** The user's original "[M+CH₃OO]−" notation was confirmed to
  mean the **acetate adduct [M+CH₃COO]−** = neutral M + CH₃COO (C₂H₃O₂, 59.013304),
  theoretical 251.0409. It lands at +3.8 ppm — slightly higher than the other ions
  (which are within ~1 ppm) but inside the ±10 ppm window; a co-eluting species at
  the upper window edge may pull the centroid.
- **Citrate [M−H]− is an unresolved blend** at the instrument's resolving power
  (heavier co-isobar ~+7–10 ppm; see `README.md` / `SESSIONS.md` and
  `R/82_citrate_resolution.R`). Inside the ±10 ppm window the measured centroid
  lands at +2.3 ppm because the window clips before the upper isobar's bulk.
- **190.0163 is unidentified** here — deliberately labelled "upstream M?" only.
  Per project rigour, any ID is tentative and would need an independent match
  (PMID/DOI). The anti-correlation finding is the actionable result.
- Statistics are **descriptive** (n = 1 slide/group); no p-values claimed.

## Reproducibility check (2026-06-23)

- ✅ All 5 `imzML`+`ibd` inputs present in `cache/imzml/` (55–190 MB `.ibd` each).
- ✅ Runs end-to-end in a single `Rscript` call with only `viridisLite`; no
  intermediate caches required (reads raw `.ibd`).
- ✅ Output PDF regenerated (2 pages, ~0.6 MB); console table + isotope-check
  numbers reproduce the headline values above.
- ⚠️ Tooling caveat (project-wide): the Glob tool can't traverse `D:` here — use
  `ls`/Grep with explicit paths to verify files.
