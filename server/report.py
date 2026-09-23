"""PDF report of an item (reportlab, in MassSpec Bench's style): summary, domain + pLDDT track, variant table,
comparisons, methods and citations. Built on demand, so it always reflects the latest variant mapping."""
from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from reportlab.lib import colors  # noqa: E402
from reportlab.lib.enums import TA_LEFT  # noqa: E402
from reportlab.lib.pagesizes import A4  # noqa: E402
from reportlab.lib.styles import ParagraphStyle  # noqa: E402
from reportlab.lib.units import mm  # noqa: E402
from reportlab.platypus import (Image as RLImage, KeepTogether, PageBreak, Paragraph, SimpleDocTemplate,  # noqa: E402
                                Spacer, Table, TableStyle)

from . import core  # noqa: E402
from .common import item_dir, load_result  # noqa: E402

INK = colors.HexColor("#13212b")
MUTED = colors.HexColor("#5b6b75")
ACCENT = colors.HexColor("#0f766e")
WARN = colors.HexColor("#b45309")
LINE = colors.HexColor("#d8e0e4")
AFC = ["#0053D6", "#65CBF3", "#FFDB13", "#FF7D45"]

S = {
    "title": ParagraphStyle("title", fontName="Helvetica-Bold", fontSize=24, leading=28, textColor=INK),
    "kicker": ParagraphStyle("kicker", fontName="Helvetica-Bold", fontSize=8, leading=10, textColor=ACCENT, spaceAfter=2),
    "h1": ParagraphStyle("h1", fontName="Helvetica-Bold", fontSize=16, leading=20, textColor=INK, spaceAfter=4),
    "h2": ParagraphStyle("h2", fontName="Helvetica-Bold", fontSize=11, leading=14, textColor=INK, spaceBefore=4, spaceAfter=2),
    "body": ParagraphStyle("body", fontName="Helvetica", fontSize=9, leading=12.5, textColor=INK, alignment=TA_LEFT),
    "muted": ParagraphStyle("muted", fontName="Helvetica", fontSize=8.5, leading=11.5, textColor=MUTED),
    "sub": ParagraphStyle("sub", fontName="Helvetica-Oblique", fontSize=7.5, leading=10, textColor=MUTED),
    "cell": ParagraphStyle("cell", fontName="Helvetica", fontSize=7.2, leading=9.2, textColor=INK),
    "tile_v": ParagraphStyle("tv", fontName="Helvetica-Bold", fontSize=15, leading=18, textColor=INK),
    "tile_l": ParagraphStyle("tl", fontName="Helvetica", fontSize=7.5, leading=9, textColor=MUTED),
}


def clean(html) -> str:
    """UI text (may hold <b>, <i>, <code> and HTML entities) → reportlab paragraph markup in a WinAnsi font."""
    import html as _h
    t = str(html or "")
    t = re.sub(r"<code>(.*?)</code>", lambda m: "\x01" + m.group(1) + "\x02", t)
    t = re.sub(r"<(/?)(b|i)>", lambda m: f"\x03{m.group(1)}{m.group(2)}\x04", t)
    t = re.sub(r"<[^>]+>", "", t)
    t = _h.unescape(t)
    t = t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    t = t.replace("\x01", "<font face='Courier'>").replace("\x02", "</font>").replace("\x03", "<").replace("\x04", ">")
    for a, b in (("−", "-"), ("→", "->"), ("↔", "<->"), ("≥", "&gt;="), ("≤", "&lt;="), ("≈", "~"), ("β", "beta"),
                 ("α", "alpha"), ("’", "'"), ("‘", "'"), ("“", '"'), ("”", '"'), ("Δ", "delta"), ("·", "-")):
        t = t.replace(a, b)
    return t


def track_png(res: dict, out: Path) -> Path | None:
    """pLDDT bars in the AlphaFold colours, UniProt domain rows and variant ticks on one axis (matplotlib)."""
    af = res.get("alphafold")
    u = res.get("uniprot")
    L = res.get("length") or 0
    if not L:
        return None
    feats = [f for f in (u or {}).get("features", []) if f["type"] in ("Domain", "Coiled coil", "Motif", "Zinc finger", "Repeat",
                                                                         "DNA binding", "Transmembrane", "Signal")]
    regions = [f for f in (u or {}).get("features", []) if f["type"] == "Region"]
    rows = [feats, regions]
    nrow = sum(1 for r in rows if r)
    fig_h = 1.6 + 0.35 * nrow + 0.4
    fig, ax = plt.subplots(figsize=(7.4, fig_h))
    y = 0
    if af:
        pl = af["plddt"]
        for i, v in enumerate(pl):
            if v is None:
                continue
            c = AFC[0] if v > 90 else AFC[1] if v >= 70 else AFC[2] if v >= 50 else AFC[3]
            ax.add_patch(Rectangle((i + 0.5, 0), 1, v / 100, color=c, lw=0))
        ax.text(-L * 0.01, 0.5, "pLDDT", ha="right", va="center", fontsize=7)
    y = -0.25
    for r in rows:
        if not r:
            continue
        # greedy lanes so overlapping features stay readable
        lanes = []
        for f in sorted(r, key=lambda f: (f["start"], -f["end"])):
            for k, end in enumerate(lanes):
                if f["start"] > end:
                    lanes[k] = f["end"]
                    lane = k
                    break
            else:
                lanes.append(f["end"])
                lane = len(lanes) - 1
            yy = y - 0.22 * lane
            col = "#0f766e" if f["type"] == "Domain" else "#7c8b93" if f["type"] == "Region" else "#b45309"
            ax.add_patch(Rectangle((f["start"] - 0.5, yy - 0.16), f["end"] - f["start"] + 1, 0.16, facecolor=col, alpha=0.85,
                                   edgecolor="white", lw=0.6))
            room = int((f["end"] - f["start"] + 1) / L * 120)          # characters that fit at 5.5 pt on a 7.4 in axis
            label = f["description"] or f["type"]
            if room >= 4:
                ax.text((f["start"] + f["end"]) / 2, yy - 0.08, label if len(label) <= room else label[:room - 1] + "…",
                        ha="center", va="center", fontsize=5.5, color="white")
        y -= 0.22 * max(1, len(lanes)) + 0.06
    vs = (u or {}).get("variants", [])
    vres = res.get("variants") or {}
    mine = [r for r in vres.get("rows", []) if r.get("pos") and not r.get("error")]
    if vs:
        for v in vs:
            ax.plot([v["start"], v["start"]], [y - 0.02, y - 0.14], color="#6a7b84", lw=0.5)
        ax.text(-L * 0.01, y - 0.08, "variants", ha="right", va="center", fontsize=7)
    for r in mine:
        ax.plot([r["pos"]], [y - 0.24], marker="v", color="#c9304f", ms=4)
        ax.text(r["pos"], y - 0.36, r["label"], ha="center", va="top", fontsize=5.5, color="#c9304f")
    ax.set_xlim(0, L + 1)
    ax.set_ylim(y - 0.5, 1.05)
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("residue (UniProt numbering)" if u else "residue", fontsize=8)
    ax.tick_params(labelsize=7)
    if af:
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(color=AFC[0], label="> 90"), Patch(color=AFC[1], label="70-90"), Patch(color=AFC[2], label="50-70"),
                           Patch(color=AFC[3], label="< 50")], title="pLDDT", fontsize=6, title_fontsize=6, ncol=4,
                  loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False)
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def build_pdf(item: str) -> Path:
    d = item_dir(item)
    res = load_result(item)
    out = d / "report.pdf"
    W, H = A4
    margin = 16 * mm
    width = W - 2 * margin
    story = [Paragraph("STRUCTURE BENCH · " + ("PROTEIN" if res["kind"] == "protein" else "UPLOADED STRUCTURE") + " REPORT", S["kicker"]),
             Paragraph(clean(res.get("name")), S["title"]), Spacer(1, 2),
             Paragraph(clean(f"{res.get('title', '')}" + (f" — {res['organism']}" if res.get("organism") else "")), S["body"]),
             Paragraph(datetime.now().strftime("Generated %d %B %Y, %H:%M"), S["muted"]), Spacer(1, 10)]
    cells = [[Paragraph(clean(t["value"]), S["tile_v"]), Paragraph(clean(t["label"]), S["tile_l"])] for t in res.get("tiles", [])]
    rows = [[Table([[c[0]], [c[1]]], style=[("LEFTPADDING", (0, 0), (-1, -1), 0)]) for c in cells[i:i + 3]] for i in range(0, len(cells), 3)]
    if rows:
        for r in rows:
            while len(r) < 3:
                r.append("")
        t = Table(rows, colWidths=[width / 3] * 3)
        t.setStyle(TableStyle([("BOX", (0, 0), (-1, -1), 0.6, LINE), ("INNERGRID", (0, 0), (-1, -1), 0.6, LINE),
                               ("TOPPADDING", (0, 0), (-1, -1), 7), ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
                               ("LEFTPADDING", (0, 0), (-1, -1), 9), ("VALIGN", (0, 0), (-1, -1), "TOP")]))
        story += [t, Spacer(1, 12)]
    story.append(Paragraph("Summary", S["h1"]))
    u = res.get("uniprot")
    if u and u.get("function"):
        story += [Paragraph("<b>Function (UniProt).</b> " + clean(re.sub(r"\s*\(PubMed:[^)]*\)", "", u["function"]))[:1800], S["body"]), Spacer(1, 4)]
    if u and u.get("location"):
        story += [Paragraph("<b>Location.</b> " + clean("; ".join(u["location"][:8])), S["body"]), Spacer(1, 4)]
    fl = []
    for f in res.get("flags", []):
        tag = {"ok": ("OK", ACCENT), "warn": ("CHECK", WARN), "info": ("NOTE", MUTED)}[f["level"]]
        fl.append([Paragraph(f"<font color='#{tag[1].hexval()[2:]}'><b>{tag[0]}</b></font>", S["cell"]), Paragraph(clean(f["text"]), S["body"])])
    for k in ("uniprot", "track", "pae", "pdb"):
        y = (res.get("texts") or {}).get(k, {}).get("yours")
        if y:
            fl.append([Paragraph("<font color='#0f766e'><b>READ</b></font>", S["cell"]), Paragraph(clean(y), S["body"])])
    if fl:
        t = Table(fl, colWidths=[16 * mm, width - 16 * mm])
        t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("LINEBELOW", (0, 0), (-1, -2), 0.4, LINE),
                               ("TOPPADDING", (0, 0), (-1, -1), 5), ("BOTTOMPADDING", (0, 0), (-1, -1), 5), ("LEFTPADDING", (0, 0), (-1, -1), 0)]))
        story.append(t)

    png = track_png(res, d / "report_track.png")
    if png:
        from PIL import Image
        with Image.open(png) as im:
            w, h = im.size
        dw = width
        dh = dw * h / w
        block = [PageBreak(), Paragraph("SEQUENCE", S["kicker"]), Paragraph("Domains and AlphaFold confidence", S["h1"]),
                 RLImage(str(png), width=dw, height=dh, hAlign="LEFT"), Spacer(1, 4)]
        tx = (res.get("texts") or {}).get("track", {})
        if tx.get("how"):
            block.append(Paragraph("<b>How to read it.</b> " + clean(tx["how"]), S["muted"]))
        if tx.get("yours"):
            block.append(Paragraph("<b>In your data.</b> " + clean(tx["yours"]), S["body"]))
        story += block

    vres = res.get("variants")
    if vres:
        story += [PageBreak(), Paragraph("VARIANTS", S["kicker"]), Paragraph(f"Variants on {clean(vres['structure_label'])} chain {vres['chain']}", S["h1"]),
                  Paragraph(clean(vres.get("yours", "")), S["body"]), Spacer(1, 6)]
        bcol = "pLDDT" if vres.get("bkind") == "plddt" or res["kind"] == "protein" else "B-factor"
        head = ["Variant", bcol, "Rel. SASA", "SS", "Domain", "AlphaMissense", "UniProt", "Reading (heuristic)"]
        data = [[Paragraph(f"<b>{h}</b>", S["cell"]) for h in head]]
        for r in vres["rows"]:
            if r.get("error"):
                data.append([Paragraph(clean(r.get("label") or r["input"]), S["cell"]), Paragraph(clean(r["error"]), S["cell"])] + [""] * 6)
                continue
            b = r.get("plddt") if bcol == "pLDDT" else r.get("bfactor")
            data.append([Paragraph(f"<b>{clean(r['label'])}</b>", S["cell"]),
                         Paragraph("—" if b is None else f"{b:.0f}", S["cell"]),
                         Paragraph("not modelled" if r.get("unresolved") else ("—" if r.get("rsa") is None else f"{100 * r['rsa']:.0f} % ({r.get('rsa_class')})"), S["cell"]),
                         Paragraph(clean(r.get("ss_word") or "—"), S["cell"]),
                         Paragraph(clean("; ".join((r.get("domain") or [])[:2]) or "—"), S["cell"]),
                         Paragraph("—" if r.get("am_score") is None else f"{r['am_score']:.2f} ({clean(r['am_class'])})", S["cell"]),
                         Paragraph(clean("; ".join(k["description"].split(";")[0] for k in r.get("uniprot_known") or []) or "—"), S["cell"]),
                         Paragraph(clean(" ".join(f["text"] for f in r.get("flags", []))), S["cell"])])
        t = Table(data, repeatRows=1, colWidths=[width * x for x in (0.08, 0.06, 0.1, 0.07, 0.13, 0.1, 0.16, 0.30)])
        t.setStyle(TableStyle([("LINEBELOW", (0, 0), (-1, 0), 0.8, INK), ("LINEBELOW", (0, 1), (-1, -1), 0.3, LINE),
                               ("TOPPADDING", (0, 0), (-1, -1), 2.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5), ("VALIGN", (0, 0), (-1, -1), "TOP"),
                               ("LEFTPADDING", (0, 0), (-1, -1), 2), ("RIGHTPADDING", (0, 0), (-1, -1), 2)]))
        story += [t, Spacer(1, 6), Paragraph("Readings are heuristics (burial, confidence and contacts), not predictions of "
                                             "pathogenicity. AlphaMissense is a prediction (Cheng et al. 2023), not a clinical classification.", S["sub"])]
    if res.get("comparisons"):
        story += [PageBreak(), Paragraph("COMPARISONS", S["kicker"]), Paragraph("Structure comparisons", S["h1"])]
        for c in res["comparisons"]:
            story += [KeepTogether([Paragraph(clean(f"{c['label_b']} chain {c['chain_b']} on {c['label_a']} chain {c['chain_a']}"), S["h2"]),
                                    Paragraph(clean(c["yours"]), S["body"]), Spacer(1, 6)])]
    story += [PageBreak(), Paragraph("METHODS", S["kicker"]), Paragraph("How this was computed", S["h1"])]
    for h, p in res.get("methods", []):
        story += [Paragraph(clean(h), S["h2"]), Paragraph(clean(p), S["body"])]
    if res.get("versions"):
        story += [Spacer(1, 8), Paragraph("Software versions", S["h2"]),
                  Paragraph(clean(", ".join(f"{k} {v}" for k, v in res["versions"].items())), S["muted"])]
    story += [Spacer(1, 8), Paragraph("References", S["h2"])]
    for c in core.citations(res):
        story.append(Paragraph(clean(c), S["muted"]))
    story += [Spacer(1, 6), Paragraph("Data: UniProt and the AlphaFold Protein Structure Database are distributed under CC BY 4.0; "
                                      "PDB entries are in the public domain (CC0).", S["sub"]),
              Spacer(1, 10), Paragraph("Interpretations are generated automatically from the numbers above and are meant to guide, "
                                       "not replace, expert review.", S["sub"])]

    def deco(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 7.5)
        canvas.setFillColor(MUTED)
        canvas.drawString(margin, 10 * mm, "Structure Bench · created by Paul H. Kim, Ph.D.")
        canvas.drawRightString(W - margin, 10 * mm, f"{doc.page}")
        canvas.setStrokeColor(ACCENT)
        canvas.setLineWidth(2)
        canvas.line(margin, H - 10 * mm, margin + 18 * mm, H - 10 * mm)
        canvas.restoreState()

    doc = SimpleDocTemplate(str(out), pagesize=A4, leftMargin=margin, rightMargin=margin, topMargin=16 * mm, bottomMargin=16 * mm,
                            title=f"Structure Bench report — {res.get('name')}", author="Paul H. Kim, Ph.D.")
    doc.build(story, onFirstPage=deco, onLaterPages=deco)
    return out
