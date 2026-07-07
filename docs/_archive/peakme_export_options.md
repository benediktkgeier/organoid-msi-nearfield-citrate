# PeakMe export — two rendering modes

For multi-sample MSI experiments, there are **two valid ways** to render the per-sample ion-image PNGs that PeakMe consumes. Each optimizes for a different annotation task. Pick the one that matches your goal.

Both modes always produce **one PeakMe zip per sample** (locked rule: never combine samples into one zip — composite PNGs prevent per-sample annotation; see memory entry `feedback_peakme_per_sample.md`).

---

## Option A — Per-PNG normalization (default)

**Script:** `R/02c_peakme_export_per_sample.R` → calls `D:/R/Projects/BIG_MSI/HP_Comparative_2026/peakme/peakme_import.R` as subprocess.

**Output zip filename:** `peakme_upload_<sample_id>.zip` (in `results/peakme/`)

### What it does

Each PNG is rendered with its own dynamic range — the brightest pixel in that ion's image for that sample becomes the brightest viridis color. Per-sample, per-ion, independent normalization.

### When to use

- **Annotating spatial pattern** — you want to see structure in every ion regardless of absolute abundance.
- **Faint ions matter** — low-abundance ions still show clear spatial structure because they're stretched to fill the dynamic range.
- **Single-sample focus** — you're working through each sample's annotations one at a time.

### Trade-off

You **cannot directly compare brightness across samples**. The same ion may look fully saturated in sec1b (low total signal) and equally saturated in sec4b (high total signal), even though the absolute intensities differ by 10×.

### Companion

`R/02d_peakme_sort_all.R` — sorts each metadata.csv by **combined cross-sample mean intensity desc** so the same m/z appears at rank 1 in every zip (easier cross-sample navigation during annotation, even though the visuals don't share a scale).

---

## Option B — Global cross-sample p99.5 clip

**Script:** `R/03_peakme_export_globalclip.R` (custom renderer — does not use `peakme_import.R`)

**Output zip filename:** `peakme_upload_globalclip_<sample_id>.zip` (in `results/peakme/`)

### What it does

For each ion, computes the **p99.5 quantile of positive intensities across ALL samples combined**. Every per-sample PNG for that ion uses the **same** intensity-to-color mapping. Linear, gamma 1.0, viridis 256.

Matches the v3-locked cross-sample comparison rule used in `top120_ion_images.pdf` and the long-term plan for the final report PDFs (`70_report_pdf.R`).

### When to use

- **Calling on-tissue vs off-tissue across samples quickly** — an ion that's faint in sec1b but bright in sec4b will look faint in sec1b's PNG and bright in sec4b's PNG → instant pattern call across samples.
- **Comparing intensity distributions** — brightness now has cross-sample meaning.
- **Detecting matrix / sample-specific artifacts** — an ion that's "matrix only" in 1 sample but real biology in another stands out immediately.

### Trade-off

- **Faint ions look faint** — low-abundance ions don't fill the dynamic range, may look uniformly dim. You may miss spatial structure in low-abundance ions.
- **No TIC spectrum chart per PNG** — Option A's `peakme_import.R` writes a small TIC chart alongside each PNG showing the mean spectrum around that m/z. Option B's custom renderer doesn't (could be added; not currently).

### Companion

metadata.csv is written in combined cross-sample mean intensity desc order out of the box — no separate sort step.

---

## How to choose

| Task | Use |
|---|---|
| First annotation pass — quickly call on-tissue / off-tissue / matrix per ion | **Option B (globalclip)** ✅ |
| Detailed annotation of spatial pattern per ion | Option A (per-PNG norm) |
| Cross-validate annotations across samples | Either — but Option B makes discordance visually obvious |
| When you don't know yet | Generate both. Disk cost ~700 MB × 2 = 1.4 GB. Render time ~25 min × 2 = 50 min. |

Many teams generate **both** and re-upload as needed during annotation rounds.

## File naming convention

To avoid uploading the wrong zip to PeakMe:

- `peakme_upload_<sample_id>.zip` — Option A (per-PNG normalization)
- `peakme_upload_globalclip_<sample_id>.zip` — Option B (global cross-sample p99.5 clip)
- `peakme_upload_COMBINED_DO_NOT_USE.zip` — quarantined output from a long-deprecated combined-zip attempt. Never upload.

## Locked rules these workflows follow

- `feedback_peakme_per_sample.md` — one zip per sample; never combine.
- `feedback_peakme_upload_format.md` — flat zip; bare `<mz>.png` filenames; metadata.csv reordered for display.
- `feedback_ion_image_display.md` — linear scale, p99.5 clip, gamma 1.0, viridis. Both options follow this; they differ only in whether p99.5 is per-PNG or global.

## Future projects

For any new multi-sample MSI project: generate both options as part of phase 1. The first round of PeakMe annotations is typically done with Option B (globalclip) to triage ions fast; a second pass with Option A is sometimes useful for detailed spatial-pattern review.
