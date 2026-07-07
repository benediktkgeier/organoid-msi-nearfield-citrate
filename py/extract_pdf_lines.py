"""Extract DRAWING annotations (ink / line / polyline) from a marked PDF, with
geometry + colour, so we can recover user-drawn split lines.

This is the geometry-aware sibling of `extract_pdf_annots.py` (which only dumps
/Rect + /Contents for text comments).  Adapted from
D:/R/PeakMe/PeakMe_GCPL/phaseR_extract_seed_annots.py so this project is
self-contained.

Output TSV columns:
  page, subtype, color, x_center, y_center_from_top, x0, y0, x1, y1, path, contents
All coords in PDF points; y is FROM TOP (origin top-left) to match the
R/102 annotation sidecar's `pdf_y_from_top`.  `path` is the FULL polyline as
"x,y;x,y;..." (y from top): the complete drawn stroke for /Ink and /PolyLine,
both endpoints for /Line.  One row per /Ink stroke.

Usage:
    python extract_pdf_lines.py <input.pdf> <output.tsv>
"""
import sys, csv
import pypdf

inp, outp = sys.argv[1], sys.argv[2]
reader = pypdf.PdfReader(inp)
rows = []
for i, page in enumerate(reader.pages, start=1):
    pageH = float(page.mediabox.height)
    annots = page.get("/Annots")
    if annots is None:
        continue
    for a in annots:
        o = a.get_object()
        sub = str(o.get("/Subtype"))
        rect = o.get("/Rect")
        contents = o.get("/Contents") or ""
        if isinstance(contents, bytes):
            try: contents = contents.decode("utf-8", "ignore")
            except Exception: contents = str(contents)
        # colour: interior (/IC) else stroke (/C), as "R,G,B"
        col = o.get("/IC") or o.get("/C") or []
        col = ",".join(str(round(float(c), 3)) for c in col) if col else ""
        x0 = y0 = x1 = y1 = ""
        xc = yc = ""
        path = ""
        if rect:
            rx0, ry0, rx1, ry1 = [float(v) for v in rect]
            xc = round((rx0 + rx1) / 2, 2)
            yc = round(pageH - (ry0 + ry1) / 2, 2)            # from top
        # Line: /L = [x0 y0 x1 y1]
        L = o.get("/L")
        if L:
            x0, y0 = round(float(L[0]), 2), round(pageH - float(L[1]), 2)
            x1, y1 = round(float(L[2]), 2), round(pageH - float(L[3]), 2)
            path = "%.2f,%.2f;%.2f,%.2f" % (x0, y0, x1, y1)
        # PolyLine: /Vertices = [x0 y0 x1 y1 ...]
        V = o.get("/Vertices")
        if V and not path:
            vf = [float(v) for v in V]
            vx = vf[0::2]; vy = vf[1::2]
            path = ";".join("%.2f,%.2f" % (x, pageH - y) for x, y in zip(vx, vy))
        # Ink: /InkList = list of strokes -> ONE ROW PER STROKE
        ink = o.get("/InkList")
        if ink:
            for si, stroke in enumerate(ink):
                fl = [float(v) for v in stroke]
                xs = fl[0::2]; ys = fl[1::2]
                if not xs:
                    continue
                spath = ";".join("%.2f,%.2f" % (x, pageH - y) for x, y in zip(xs, ys))
                rows.append(dict(page=i, subtype=sub, color=col,
                    x_center=round(sum(xs) / len(xs), 2),
                    y_center_from_top=round(pageH - sum(ys) / len(ys), 2),
                    x0=round(min(xs), 2), y0=round(pageH - max(ys), 2),
                    x1=round(max(xs), 2), y1=round(pageH - min(ys), 2),
                    path=spath,
                    contents=("stroke%d/%d " % (si + 1, len(ink)))
                              + str(contents).strip().replace("\n", " ")))
            continue   # strokes emitted above
        rows.append(dict(page=i, subtype=sub, color=col,
                         x_center=xc, y_center_from_top=yc,
                         x0=x0, y0=y0, x1=x1, y1=y1,
                         path=path,
                         contents=str(contents).strip().replace("\n", " ")))

with open(outp, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["page", "subtype", "color", "x_center",
        "y_center_from_top", "x0", "y0", "x1", "y1", "path", "contents"],
        delimiter="\t")
    w.writeheader()
    for r in rows:
        w.writerow(r)
print(f"Wrote {len(rows)} annotation records to {outp}")
