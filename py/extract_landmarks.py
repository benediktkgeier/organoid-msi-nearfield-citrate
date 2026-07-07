#!/usr/bin/env python
# Extract numbered organoid-center marks from the annotated landmark sheets.
# Each mark = any annotation whose comment/content is a number (the organoid id);
# the point is the annotation centre. Output: page, number, ptx, pty (PDF pt, top-down).
# Usage: python extract_landmarks.py [annotated.pdf]
#   default annotated.pdf = figures/register_if/landmark_sheets_marked.pdf
import fitz, csv, sys, os, re
src = sys.argv[1] if len(sys.argv) > 1 else "figures/register_if/landmark_sheets_marked.pdf"
out = "results/register_if/landmark_points_raw.csv"
d = fitz.open(src); rows = []
for pi in range(len(d)):
    for a in d[pi].annots() or []:
        c = (a.info.get("content", "") or "").strip()
        m = re.search(r"\d+", c)            # first number in the comment = organoid id
        if m:
            r = a.rect
            rows.append([pi, int(m.group()), round((r.x0 + r.x1) / 2, 2), round((r.y0 + r.y1) / 2, 2)])
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["page", "number", "ptx", "pty"])
    for r in sorted(rows): w.writerow(r)
print("wrote", len(rows), "landmark points ->", out)
for r in sorted(rows): print(r)
