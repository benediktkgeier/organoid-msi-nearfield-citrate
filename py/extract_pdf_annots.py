"""Extract Adobe annotations (FreeText / sticky-note / markup comments) from a
PDF, page by page.  Writes a TSV with columns:
    page, type, rt, x_left, y_bot_pdfcoord, x_right, y_top_pdfcoord, contents

Generic helper (copied verbatim from D:/R/PeakMe/PeakMe_GCPL/phaseZC_extract_annots.py)
so this project is self-contained.  Usage:
    python extract_pdf_annots.py <input.pdf> <output.tsv>
"""
import sys
import pypdf
import csv

inp = sys.argv[1]
outp = sys.argv[2]
reader = pypdf.PdfReader(inp)
rows = []
for i, page in enumerate(reader.pages, start=1):
    annots = page.get("/Annots")
    if annots is None:
        continue
    # Annots may be an indirect reference list
    for a in annots:
        obj = a.get_object()
        subtype = obj.get("/Subtype")
        rect = obj.get("/Rect")
        contents = obj.get("/Contents") or ""
        if isinstance(contents, bytes):
            try:
                contents = contents.decode("utf-8")
            except Exception:
                contents = str(contents)
        rt = obj.get("/RT")        # for reply-type annotations
        # For FreeText / Note / Stamp annotations, /Contents has the text.
        # For markup annotations like /Highlight / /StrikeOut, the comment
        # lives in /Contents of the popup or markup itself.
        rows.append({
            "page": i,
            "type": str(subtype) if subtype else "",
            "rt": str(rt) if rt else "",
            "x_left": float(rect[0]) if rect else "",
            "y_bot_pdfcoord": float(rect[1]) if rect else "",
            "x_right": float(rect[2]) if rect else "",
            "y_top_pdfcoord": float(rect[3]) if rect else "",
            "contents": str(contents).strip().replace("\n", " | "),
        })

with open(outp, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["page","type","rt","x_left","y_bot_pdfcoord",
                                        "x_right","y_top_pdfcoord","contents"],
                        delimiter="\t")
    w.writeheader()
    for r in rows:
        w.writerow(r)
print(f"Wrote {len(rows)} annotation records to {outp}")
