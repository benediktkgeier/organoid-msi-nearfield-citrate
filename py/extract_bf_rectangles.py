#!/usr/bin/env python
# Extract the user-drawn brightfield organoid RECTANGLES (one per section, drawn on
# annotate_brightfield_<slide>_guided_rectangle.pdf) -> results/register_if/bf_rectangles.csv.
#
# Handles two gotchas:
#   1. R's pdf() adds a ~3.7% axis-padding margin, so the image does NOT fill the page.
#      We detect the inner image bbox (render w/o annots, find non-white extent) and map
#      each rectangle to IMAGE fractions (not page fractions).
#   2. The user mis-numbered the rectangles relative to the IF sections, so we pair-swap
#      within pairs: 1<->2, 3<->4, 5<->6 (SWAP below).
#
# Usage: python extract_bf_rectangles.py
import fitz, numpy as np, csv
DPI = 120; SWAP = {1: 2, 2: 1, 3: 4, 4: 3, 5: 6, 6: 5}
rows = []
for sl in ["sl6A", "sl4A"]:
    d = fitz.open("figures/register_if/annotate_brightfield_%s_guided_rectangle.pdf" % sl); pg = d[0]
    pix = pg.get_pixmap(dpi=DPI, annots=False)                       # image only -> inner bbox
    arr = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)[:, :, 0]
    nw = arr < 250; cols = np.where(nw.any(axis=0))[0]; rr = np.where(nw.any(axis=1))[0]
    c0, c1, r0, r1 = cols.min(), cols.max(), rr.min(), rr.max(); s = DPI / 72.0
    sq = [a for a in (pg.annots() or []) if a.type[1] == "Square"]
    used = set(int(a.info.get("content", "").strip()) for a in sq if a.info.get("content", "").strip().isdigit())
    missing = [n for n in range(1, 7) if n not in used]; mi = 0
    for a in sq:
        cc = a.info.get("content", "").strip()
        sec = int(cc) if cc.isdigit() else missing[mi]
        if not cc.isdigit(): mi += 1
        sec = SWAP[sec]
        r = a.rect
        fx0 = (r.x0 * s - c0) / (c1 - c0); fx1 = (r.x1 * s - c0) / (c1 - c0)
        fy0 = (r.y0 * s - r0) / (r1 - r0); fy1 = (r.y1 * s - r0) / (r1 - r0)
        rows.append([sl, sec, round(float(fx0), 5), round(float(fy0), 5), round(float(fx1), 5), round(float(fy1), 5)])
with open("results/register_if/bf_rectangles.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["slide", "section", "fx0", "fy0", "fx1", "fy1"])
    for r in sorted(rows): w.writerow(r)
print("wrote", len(rows), "brightfield rectangles -> results/register_if/bf_rectangles.csv")
