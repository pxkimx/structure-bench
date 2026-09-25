"""PDF report of an item (reportlab, in MassSpec Bench's style).

It follows the page section by section — protein (UniProt), predicted structure (confidence, domains, PAE and the
3D model), experimental structures, variants and comparisons — and every figure and table carries the same two
texts the page shows: what it is and how to read it, and what it says for this protein. Built on demand, so it
always reflects the latest variant mapping and comparisons.

3D images: when the report is made from the page, the page sends snapshots of the 3D viewer (the AlphaFold model
from two sides, the mapped variants, each superposition). Without them (the assistant, or a direct download) the
server draws the CA trace instead and says so under the figure.
"""
from __future__ import annotations

import base64
import re
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from reportlab.lib import colors  # noqa: E402
from reportlab.lib.enums import TA_LEFT  # noqa: E402
from reportlab.lib.pagesizes import A4  # noqa: E402
from reportlab.lib.styles import ParagraphStyle  # noqa: E402
from reportlab.lib.units import mm  # noqa: E402
from reportlab.platypus import (CondPageBreak, Image as RLImage, KeepTogether, PageBreak, Paragraph,  # noqa: E402
                                SimpleDocTemplate, Spacer, Table, TableStyle)

from . import core  # noqa: E402
from .common import item_dir, load_result, log_exc  # noqa: E402

INK = colors.HexColor("#13212b")
MUTED = colors.HexColor("#5b6b75")
ACCENT = colors.HexColor("#0f766e")
WARN = colors.HexColor("#b45309")
LINE = colors.HexColor("#d8e0e4")
PANEL = colors.HexColor("#f4f7f8")
AFC = ["#0053D6", "#65CBF3", "#FFDB13", "#FF7D45"]
METHOD_C = {"X-ray": "#2563eb", "EM": "#9333ea", "NMR": "#15803d"}
BURIAL_C = {"buried": "#2a64b0", "intermediate": "#8e9aa3", "exposed": "#d97706"}

S = {
    "title": ParagraphStyle("title", fontName="Helvetica-Bold", fontSize=24, leading=28, textColor=INK),
    "kicker": ParagraphStyle("kicker", fontName="Helvetica-Bold", fontSize=8, leading=10, textColor=ACCENT, spaceAfter=2),
    "h1": ParagraphStyle("h1", fontName="Helvetica-Bold", fontSize=16, leading=20, textColor=INK, spaceAfter=4),
    "h2": ParagraphStyle("h2", fontName="Helvetica-Bold", fontSize=11, leading=14, textColor=INK, spaceBefore=6, spaceAfter=3),
    "h3": ParagraphStyle("h3", fontName="Helvetica-Bold", fontSize=9.5, leading=12, textColor=INK, spaceBefore=4, spaceAfter=1),
    "body": ParagraphStyle("body", fontName="Helvetica", fontSize=9, leading=12.5, textColor=INK, alignment=TA_LEFT),
    "muted": ParagraphStyle("muted", fontName="Helvetica", fontSize=8.5, leading=11.5, textColor=MUTED),
    "how": ParagraphStyle("how", fontName="Helvetica", fontSize=8.3, leading=11.2, textColor=MUTED, spaceBefore=3),
    "yours": ParagraphStyle("yours", fontName="Helvetica", fontSize=8.8, leading=12, textColor=INK),
    "yoursbox": ParagraphStyle("yoursbox", fontName="Helvetica", fontSize=8.8, leading=12, textColor=INK, backColor=PANEL,
                               borderColor=LINE, borderWidth=0.4, borderPadding=(5, 6, 5, 6), leftIndent=6, rightIndent=6),
    "sub": ParagraphStyle("sub", fontName="Helvetica-Oblique", fontSize=7.5, leading=10, textColor=MUTED),
    "cell": ParagraphStyle("cell", fontName="Helvetica", fontSize=7.2, leading=9.2, textColor=INK),
    "cellb": ParagraphStyle("cellb", fontName="Helvetica-Bold", fontSize=7.2, leading=9.2, textColor=INK),
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



# ---------------------------------------------------------------- small helpers
def P(text, style="body"):
    """Paragraph from app text; anything ReportLab cannot parse goes in as plain text rather than losing the report."""
    try:
        return Paragraph(clean(text), S[style])
    except Exception:  # noqa: BLE001
        from xml.sax.saxutils import escape
        return Paragraph(escape(re.sub(r"<[^>]+>", "", str(text or ""))), S[style])


def _fmt(v, nd=1):
    return "—" if v is None else f"{v:.{nd}f}"


def _pct(v):
    return "—" if v is None else f"{100 * v:.0f} %"


class Doc:
    """Collects flowables; knows the text width and where the figure PNGs go."""

    def __init__(self, d: Path, width: float):
        self.d = d
        self.width = width
        self.fig_dir = d / "report_figs"
        self.fig_dir.mkdir(exist_ok=True)
        self.story = []

    def add(self, *xs):
        self.story.extend(xs)

    def section(self, kicker, title, lede=""):
        self.add(PageBreak(), P(kicker.upper(), "kicker"), P(title, "h1"))
        if lede:
            self.add(P(lede, "muted"), Spacer(1, 6))

    def explain(self, how="", yours="", note=""):
        """The two texts every figure and table has on the page: what it shows, and what it says here."""
        out = []
        if how:
            out.append(P("<b>What this shows.</b> " + str(how), "how"))
        if yours:
            # a shaded paragraph, not a one-cell table: a table row cannot break across pages, and a long reading
            # (titin lists hundreds of domains) then stopped the whole PDF with a LayoutError
            out += [Spacer(1, 7), P("<b>In your data.</b> " + str(yours), "yoursbox"), Spacer(1, 3)]
        if note:
            out.append(P(note, "sub"))
        return out

    def image(self, path_or_bytes, max_w=None, max_h=None):
        from PIL import Image
        max_w = max_w or self.width
        max_h = max_h or 190 * mm
        src = path_or_bytes
        if isinstance(src, (bytes, bytearray)):
            import io
            buf = io.BytesIO(src)
            with Image.open(buf) as im:
                w, h = im.size
            buf.seek(0)
            src = buf
        else:
            with Image.open(src) as im:
                w, h = im.size
            src = str(src)
        sc = min(max_w / w, max_h / h)
        return RLImage(src, width=w * sc, height=h * sc, hAlign="LEFT")

    def figure(self, title, img, how="", yours="", note="", sub=""):
        head = [P(title, "h2")] + ([P(sub, "sub")] if sub else [])
        self.add(CondPageBreak(90 * mm), KeepTogether(head + [img]), *self.explain(how, yours, note), Spacer(1, 10))

    def table(self, title, columns, rows, widths, how="", yours="", note="", max_rows=None):
        if title:
            self.add(CondPageBreak(40 * mm), P(title, "h2"))
        self.add(*self.explain(how, yours))
        shown = rows if max_rows is None else rows[:max_rows]
        data = [[P(f"<b>{c}</b>", "cell") for c in columns]] + [
            [c if hasattr(c, "wrap") else P(c, "cell") for c in r] for r in shown]
        t = Table(data, repeatRows=1, colWidths=[self.width * x for x in widths], splitByRow=1)
        t.setStyle(TableStyle([("LINEBELOW", (0, 0), (-1, 0), 0.8, INK), ("LINEBELOW", (0, 1), (-1, -1), 0.3, LINE),
                               ("TOPPADDING", (0, 0), (-1, -1), 2.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
                               ("VALIGN", (0, 0), (-1, -1), "TOP"), ("LEFTPADDING", (0, 0), (-1, -1), 2),
                               ("RIGHTPADDING", (0, 0), (-1, -1), 2)]))
        self.add(Spacer(1, 4), t)
        if max_rows is not None and len(rows) > max_rows:
            self.add(P(f"Showing {max_rows} of {len(rows)} rows; the full list is on the page and in its CSV download.", "sub"))
        if note:
            self.add(P(note, "sub"))
        self.add(Spacer(1, 10))

    def save(self, fig, name):
        out = self.fig_dir / name
        fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        return out


def _snapshot(snaps: dict | None, key: str):
    """A PNG the page sent (data URL), as bytes; None when absent or unreadable."""
    s = (snaps or {}).get(key)
    if not s or "," not in s:
        return None
    try:
        raw = base64.b64decode(s.split(",", 1)[1])
    except Exception:  # noqa: BLE001
        return None
    try:          # trim the white margin the viewer leaves, so the molecule fills the space it is given
        import io
        from PIL import Image, ImageChops
        with Image.open(io.BytesIO(raw)) as im:
            rgb = im.convert("RGB")
            box = ImageChops.difference(rgb, Image.new("RGB", rgb.size, (255, 255, 255))).point(lambda v: 255 if v > 12 else 0).getbbox()
            if box:
                pad = max(8, int(0.02 * max(rgb.size)))
                box = (max(0, box[0] - pad), max(0, box[1] - pad), min(rgb.width, box[2] + pad), min(rgb.height, box[3] + pad))
                out = io.BytesIO()
                rgb.crop(box).save(out, "PNG")
                return out.getvalue()
    except Exception:  # noqa: BLE001
        log_exc("snapshot crop")
    return raw


# ---------------------------------------------------------------- figures (matplotlib)
def bands_png(doc: Doc, af: dict) -> Path:
    fr = [af.get("frac_very_high") or 0, max(0.0, (af.get("frac_confident") or 0) - (af.get("frac_very_high") or 0)),
          af.get("frac_low") or 0, af.get("frac_very_low") or 0]
    lab = ["very high (> 90)", "confident (70–90)", "low (50–70)", "very low (< 50)"]
    fig, ax = plt.subplots(figsize=(7.4, 0.95))
    x = 0
    for f, c, l in zip(fr, AFC, lab):
        ax.barh(0, f, left=x, color=c, height=0.6)
        if f >= 0.06:
            ax.text(x + f / 2, 0, f"{100 * f:.0f}%", ha="center", va="center", fontsize=7.5,
                    color="white" if c in (AFC[0],) else "#13212b")
        x += f
    ax.set_xlim(0, 1)
    ax.axis("off")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=c, label=l) for c, l in zip(AFC, lab)], ncol=4, fontsize=6.5, frameon=False,
              loc="upper center", bbox_to_anchor=(0.5, -0.05))
    return doc.save(fig, "bands.png")


def pae_png(doc: Doc, res: dict) -> Path | None:
    p = doc.d / "pae.json"
    if not p.exists():
        return None
    import json
    P_ = json.loads(p.read_text())
    k, n, b = P_["size"], P_["n"], P_.get("bin", 1)
    m = np.asarray(P_["values"], float).reshape(k, k)
    m[m < 0] = np.nan
    fig, ax = plt.subplots(figsize=(5.2, 4.5))
    im = ax.imshow(m, cmap="Greens_r", vmin=0, vmax=P_.get("max") or 31.75, extent=(0.5, n + 0.5, n + 0.5, 0.5),
                   interpolation="nearest")
    for blk in (res.get("pae") or {}).get("blocks", []):
        for v in (blk["start"] - 0.5, blk["end"] + 0.5):
            ax.axhline(v, color="#c9304f", lw=0.6, alpha=0.7)
            ax.axvline(v, color="#c9304f", lw=0.6, alpha=0.7)
        ax.text((blk["start"] + blk["end"]) / 2, -n * 0.015, blk["name"][:18], ha="center", va="bottom", fontsize=6,
                color="#c9304f")
    ax.set_xlabel("scored residue", fontsize=8)
    ax.set_ylabel("aligned residue", fontsize=8)
    ax.tick_params(labelsize=7)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("expected position error (Å)", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    if b > 1:
        ax.set_title(f"averaged in {b}×{b}-residue blocks for display", fontsize=6.5, color="#5b6b75")
    return doc.save(fig, "pae.png")


def _ca(doc: Doc, meta: dict, chain: str | None = None, file_key: str = "file"):
    """CA coordinates, residue numbers and B-factors of one chain of a registered structure."""
    s, _ = core._parsed(doc.d / meta[file_key])
    ch = chain or meta.get("main_chain")
    model = s[0]
    if ch not in model:
        ch = next(iter(model)).id
    xyz, num, b = [], [], []
    for r in model[ch]:
        if "CA" in r and r.id[0] == " ":
            xyz.append(r["CA"].coord)
            num.append(r.id[1])
            b.append(r["CA"].get_bfactor())
    return np.asarray(xyz, float), num, b


def _pca_frame(xyz):
    """Centre and axes of a point cloud, longest axis first — so a picture fills the page whatever the shape."""
    c = xyz.mean(0)
    _, _, vt = np.linalg.svd(xyz - c, full_matrices=False)
    return c, vt


def _draw_view(ax, chains, c, vt, a, b, depth):
    """Draw CA traces projected on axes a, b; nearer segments (along `depth`) drawn later and thicker."""
    segs = []
    for xyz, cols in chains:
        q = (xyz - c) @ vt.T
        for i in range(len(q) - 1):
            if np.linalg.norm(xyz[i + 1] - xyz[i]) < 4.5:             # no line across chain breaks
                segs.append(((q[i, a], q[i + 1, a]), (q[i, b], q[i + 1, b]), (q[i, depth] + q[i + 1, depth]) / 2, cols[i]))
    if not segs:
        return
    ds = np.array([s_[2] for s_ in segs])
    lo, hi = ds.min(), ds.max() if ds.max() > ds.min() else ds.min() + 1
    for xs, ys, dz, col in sorted(segs, key=lambda s_: s_[2]):
        t = (dz - lo) / (hi - lo)
        ax.plot(xs, ys, color=col, lw=0.9 + 2.2 * t, solid_capstyle="round", alpha=0.55 + 0.45 * t)
    ax.set_aspect("equal")
    ax.axis("off")


def trace_png(doc: Doc, name, xyz, cols, marks=None, overlay=None) -> Path:
    """Fallback 3D picture: the CA trace seen along two of its principal axes (depth shown by line width).
    Used when the page's 3D viewer did not send a snapshot. `overlay` = (xyz, colour) drawn underneath."""
    chains = ([(overlay[0], [overlay[1]] * len(overlay[0]))] if overlay is not None else []) + [(xyz, list(cols))]
    allx = np.vstack([ch[0] for ch in chains])
    c, vt = _pca_frame(allx)
    q = (allx - c) @ vt.T
    ext = q.max(0) - q.min(0)
    long_ = ext[0] > 2.2 * ext[1]
    if long_:                      # a rod: two views stacked, each as wide as the page
        fig, axs = plt.subplots(2, 1, figsize=(7.4, max(2.4, 7.4 * (ext[1] + ext[2]) / ext[0] + 0.9)))
    else:
        fig, axs = plt.subplots(1, 2, figsize=(7.4, max(2.8, 3.7 * max(ext[1], ext[2]) / ext[0] + 0.4)))
    for ax, (a, b, dpt, title) in zip(axs, ((0, 1, 2, "view 1"), (0, 2, 1, "view 2, turned 90° about the long axis"))):
        _draw_view(ax, chains, c, vt, a, b, dpt)
        for (p, lab) in (marks or []):
            pq = (np.asarray(p) - c) @ vt.T
            ax.scatter([pq[a]], [pq[b]], s=30, color="#c026d3", zorder=5, edgecolors="white", linewidths=0.6)
            ax.annotate(lab, (pq[a], pq[b]), xytext=(3, 3), textcoords="offset points", fontsize=6.5, color="#c026d3")
        ax.set_title(title, fontsize=7, color="#5b6b75")
    fig.tight_layout()
    return doc.save(fig, name)


def _afc(v):
    return AFC[0] if v > 90 else AFC[1] if v >= 70 else AFC[2] if v >= 50 else AFC[3]


def coverage_png(doc: Doc, u: dict, af: dict | None, entries: list) -> Path:
    L = u["length"]
    doms = [f for f in u["features"] if f["type"] == "Domain"]
    n = len(entries)
    fig_h = 0.9 + 0.13 * n
    fig, ax = plt.subplots(figsize=(7.4, fig_h))
    y0 = 0
    if af:
        for i, v in enumerate(af["plddt"]):
            if v is not None:
                ax.add_patch(Rectangle((i + 0.5, y0), 1, 0.5, color=_afc(v), lw=0))
    for f in doms:
        ax.add_patch(Rectangle((f["start"] - 0.5, y0 - 0.9), f["end"] - f["start"] + 1, 0.7, color="#0f766e", lw=0))
        if (f["end"] - f["start"]) / L > 0.08:
            ax.text((f["start"] + f["end"]) / 2, y0 - 0.55, (f["description"] or "Domain")[:22], ha="center", va="center",
                    fontsize=5.5, color="white")
    for i, e in enumerate(entries):
        y = y0 - 1.6 - i * 0.9
        ax.add_patch(Rectangle((0.5, y), L, 0.6, color="#eef2f4", lw=0))
        for sgm in e["segments"]:
            s0, s1 = max(1, sgm["start"]), min(L, sgm["end"])
            ax.add_patch(Rectangle((s0 - 0.5, y), s1 - s0 + 1, 0.6, color=METHOD_C.get(e["method"], "#7c8b93"), lw=0))
        ax.text(-L * 0.01, y + 0.3, f"{e['id']} {e['method']} {e.get('resolution_text') or ''}".strip(), ha="right",
                va="center", fontsize=5.5)
    ax.set_xlim(0, L + 1)
    ax.set_ylim(y0 - 1.6 - n * 0.9, 0.7)
    ax.set_yticks([])
    for s_ in ("left", "right", "top"):
        ax.spines[s_].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.set_xlabel("residue (UniProt numbering)", fontsize=8)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=c, label=m) for m, c in METHOD_C.items()], ncol=3, fontsize=6.5, frameon=False,
              loc="lower center", bbox_to_anchor=(0.5, 1.0))
    return doc.save(fig, "coverage.png")


def variants_png(doc: Doc, res: dict, V: dict) -> Path | None:
    rows = [r for r in V["rows"] if not r.get("error") and r.get("pos")]
    if not rows:
        return None
    L = res.get("length") or max(r["pos"] for r in rows)
    use_am = any(r.get("am_score") is not None for r in rows)
    fig, ax = plt.subplots(figsize=(7.4, 2.5))
    for r in rows:
        y = r.get("am_score") if use_am else (r.get("plddt") if r.get("plddt") is not None else r.get("bfactor"))
        if y is None:
            y = 0.0 if use_am else 0
        c = BURIAL_C.get(r.get("rsa_class"), "#c7ccd1")
        ax.plot([r["pos"], r["pos"]], [0, y], color=c, lw=1)
        ax.scatter([r["pos"]], [y], s=26, color=c, zorder=3, edgecolors="white", linewidths=0.6)
        if len(rows) <= 30:
            ax.text(r["pos"], y + (0.03 if use_am else 3), r["label"], fontsize=5.8, ha="center", va="bottom", rotation=90)
    if use_am:
        ax.axhspan(0.564, 1.05, color="#c9304f", alpha=0.06)
        ax.axhspan(0, 0.34, color="#2a64b0", alpha=0.06)
        ax.set_ylim(0, 1.25)
        ax.set_ylabel("AlphaMissense score", fontsize=8)
    else:
        ax.set_ylim(0, 125)
        ax.set_ylabel("pLDDT" if V.get("bkind") == "plddt" else "B-factor", fontsize=8)
    ax.set_xlim(0, L + 1)
    ax.set_xlabel("residue", fontsize=8)
    ax.tick_params(labelsize=7)
    for s_ in ("right", "top"):
        ax.spines[s_].set_visible(False)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=c, label=k) for k, c in BURIAL_C.items()], title="burial",
              fontsize=6.5, title_fontsize=6.5, frameon=False, ncol=3, loc="upper right")
    return doc.save(fig, "variants.png")


def deviation_png(doc: Doc, c: dict, k: int) -> Path | None:
    D = c.get("deviation") or []
    if not D:
        return None
    pos = [d["pos"] for d in D]
    dev = [d["dev"] for d in D]
    fig, ax = plt.subplots(figsize=(7.4, 2.6))
    for sg in c.get("segments") or []:
        ax.axvspan(sg["start"] - 0.5, sg["end"] + 0.5, color="#c9304f", alpha=0.08)
    ax.plot(pos, dev, color="#13212b", lw=0.9)
    ax.scatter([p for p, d in zip(pos, D) if not d.get("core", True)], [d["dev"] for d in D if not d.get("core", True)],
               s=6, color="#d97706", zorder=3, label="outside the fitted core")
    ax.axhline(3, color="#c9304f", lw=0.7, ls="--")
    ax.text(pos[-1], 3.1, "3 Å", fontsize=6.5, color="#c9304f", ha="right", va="bottom")
    ax.set_ylabel("CA deviation (Å)", fontsize=8)
    ax.set_xlabel(f"residue ({c.get('numbering', 'UniProt')} numbering)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_ylim(0, max(4, max(dev) * 1.08))
    pl = [d.get("plddt") for d in D]
    if any(v is not None for v in pl):
        ax2 = ax.twinx()
        ax2.scatter(pos, [v if v is not None else np.nan for v in pl], s=3, c=[_afc(v) if v is not None else "#ccc" for v in pl])
        ax2.set_ylim(0, 105)
        ax2.set_ylabel("pLDDT", fontsize=7, color="#5b6b75")
        ax2.tick_params(labelsize=6.5, colors="#5b6b75")
    for s_ in ("top",):
        ax.spines[s_].set_visible(False)
    if any(not d.get("core", True) for d in D):
        ax.legend(fontsize=6.5, frameon=False, loc="upper left")
    return doc.save(fig, f"deviation_{k}.png")


# ---------------------------------------------------------------- sections
def cover(doc: Doc, res: dict):
    protein = res["kind"] == "protein"
    doc.add(P("STRUCTURE BENCH · " + ("PROTEIN" if protein else "UPLOADED STRUCTURE") + " REPORT", "kicker"),
            P(res.get("name"), "title"), Spacer(1, 2),
            P(f"{res.get('title', '')}" + (f" — <i>{res['organism']}</i>" if res.get("organism") else ""), "body"),
            P(f"Analysis created {res.get('created', '?')} · report generated {datetime.now():%d %B %Y, %H:%M} · "
              f"Structure Bench {res.get('app_version', '')}", "muted"), Spacer(1, 10))
    cells = [[P(str(t["value"]), "tile_v"), P(t["label"], "tile_l")] for t in res.get("tiles", [])]
    rows = [[Table([[c[0]], [c[1]]], style=[("LEFTPADDING", (0, 0), (-1, -1), 0)]) for c in cells[i:i + 3]]
            for i in range(0, len(cells), 3)]
    if rows:
        for r in rows:
            while len(r) < 3:
                r.append("")
        t = Table(rows, colWidths=[doc.width / 3] * 3)
        t.setStyle(TableStyle([("BOX", (0, 0), (-1, -1), 0.6, LINE), ("INNERGRID", (0, 0), (-1, -1), 0.6, LINE),
                               ("TOPPADDING", (0, 0), (-1, -1), 7), ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
                               ("LEFTPADDING", (0, 0), (-1, -1), 9), ("VALIGN", (0, 0), (-1, -1), "TOP")]))
        doc.add(t, Spacer(1, 12))
    doc.add(P("What your data says", "h1"),
            P("The findings the page lists at the top, in the same order. CHECK marks something that changes how the "
              "results should be read; NOTE is context.", "muted"), Spacer(1, 4))
    for f in res.get("flags", []) or [{"level": "ok", "text": "Everything needed was available."}]:
        tag, col = {"ok": ("OK", "#0f766e"), "warn": ("CHECK", "#b45309"), "info": ("NOTE", "#5b6b75")}[f["level"]]
        row = [P(f"<font color='{col}'><b>{tag}</b></font>", "cell"), P(f["text"], "body")]
        if len(str(f["text"])) > 2500:
            doc.add(P(f"<b>{tag}</b> " + str(f["text"]), "body"), Spacer(1, 5))
            continue
        t = Table([row], colWidths=[16 * mm, doc.width - 16 * mm])
        t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("LINEBELOW", (0, 0), (-1, -1), 0.4, LINE),
                               ("TOPPADDING", (0, 0), (-1, -1), 5), ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
                               ("LEFTPADDING", (0, 0), (-1, -1), 0)]))
        doc.add(t)
    doc.add(Spacer(1, 8), P("Contents: " + " · ".join(
        (["Protein (UniProt)", "Predicted structure (AlphaFold)", "Experimental structures"] if protein else
         ["Your structure"]) + (["Variants"] if res.get("variants") else []) +
        (["Structure comparisons"] if res.get("comparisons") else []) + ["Methods and references"]), "muted"))


FEATURE_HOW = {
    "Domains & regions": "Stretches UniProt curators assigned a structural or functional identity: folded domains, coiled "
                         "coils, motifs, repeats and named regions (for example a disordered tail or a binding region). "
                         "Positions are UniProt residue numbers, which are also the AlphaFold model's numbers.",
    "Sites": "Single residues with a known role: active sites, binding sites, cleavage sites and other functional "
             "positions.",
    "PTMs": "Post-translational modifications UniProt lists (phosphorylation, acetylation, lipidation, glycosylation, "
            "cross-links, disulfides), from experiments or by similarity — the evidence is on the UniProt entry.",
    "Natural variants": "Sequence variants found in people (or natural populations) that UniProt curates, with the "
                        "disease or effect noted where known. Positions are in the canonical isoform.",
}


def protein_section(doc: Doc, res: dict):
    u = res.get("uniprot")
    doc.section("UniProt", "Protein", "The curated entry: what the protein does, where it lives, what goes wrong in "
                                      "disease, and every annotated feature.")
    if not u:
        doc.add(*doc.explain(yours="UniProt could not be reached and this entry is not in the local cache, so this "
                                   "section is empty. Refresh the protein on the page when you are back online."))
        return
    kv = [["Accession", f"{u['accession']} · {u.get('id', '')} · {'reviewed (Swiss-Prot)' if u.get('reviewed') else 'unreviewed (TrEMBL)'}"],
          ["Protein", u.get("name")], ["Gene", u.get("gene")],
          ["Organism", f"<i>{u.get('organism')}</i>" + (f" ({u['common']})" if u.get("common") else "") + f" · taxon {u.get('taxon')}"],
          ["Length", f"{u['length']:,} residues" + (f" · signal peptide 1–{u['signal_length']}" if u.get("signal_length") else "")],
          ["Entry", f"version {u.get('version', '?')} · last updated {u.get('last_updated', '?')}"],
          ["Function", re.sub(r"\s*\(PubMed:[^)]*\)", "", u.get("function") or "—")],
          ["Subcellular location", "; ".join(u.get("location") or []) or "—"]]
    if res.get("isoforms"):
        kv.append(["Isoform models", "; ".join(f"{i['entry']} ({i['length']} aa, mean pLDDT {_fmt(i.get('mean_plddt'), 0)})"
                                               for i in res["isoforms"])])
    t = Table([[P(f"<b>{k}</b>", "cell"), P(v, "body")] for k, v in kv], colWidths=[32 * mm, doc.width - 32 * mm],
              splitByRow=1)
    t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("LINEBELOW", (0, 0), (-1, -1), 0.3, LINE),
                           ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
                           ("LEFTPADDING", (0, 0), (-1, -1), 0)]))
    tx = (res.get("texts") or {}).get("uniprot", {})
    doc.add(P(f"{u.get('gene') or u['accession']} — {u.get('name')}", "h2"), t, *doc.explain(tx.get("how"), tx.get("yours")),
            Spacer(1, 8))
    dz = [d for d in u.get("diseases", []) if d.get("name") or d.get("description")]
    if dz:
        doc.table("Disease involvement", ["Disease", "MIM", "UniProt description"],
                  [[f"<b>{d.get('name') or 'Other'}</b>" + (f" ({d['acronym']})" if d.get("acronym") else ""),
                    d.get("mim") or "—", d.get("description") or ""] for d in dz], (0.26, 0.08, 0.66),
                  how="Diseases UniProt curators link to this protein, with the OMIM (MIM) number and UniProt's summary "
                      "of each. The variants behind them are in the Natural variants table below.",
                  yours=f"{len(dz)} disease{'s' if len(dz) != 1 else ''} annotated for {u.get('gene') or u['accession']}.")
    for grp, feats in (u.get("groups") or {}).items():
        if not feats:
            continue
        rows = []
        for f in feats:
            change = f"{f.get('ref') or ''}→{'/'.join(f.get('alt') or []) or '?'}" if f["type"] == "Natural variant" else ""
            rows.append([f["type"], str(f["start"]), str(f["end"]), change, re.sub(r";\s*dbSNP:\S+", "", f.get("description") or "")])
        types = {}
        for f in feats:
            types[f["type"]] = types.get(f["type"], 0) + 1
        doc.table(f"Annotated features — {grp} ({len(feats)})", ["Type", "Start", "End", "Change", "Description"], rows,
                  (0.16, 0.06, 0.06, 0.08, 0.64), how=FEATURE_HOW.get(grp, ""),
                  yours=", ".join(f"{n} {t.lower()}{'s' if n != 1 and not t.endswith('s') else ''}" for t, n in
                                  sorted(types.items(), key=lambda x: -x[1])) + ".", max_rows=400)


def structure_section(doc: Doc, res: dict, snaps: dict | None):
    protein = res["kind"] == "protein"
    af = res.get("alphafold")
    tx = res.get("texts") or {}
    doc.section("AlphaFold" if protein else "Structure", "Predicted structure" if protein else "Your structure",
                "How confident the model is along the sequence, which parts are placed relative to each other, and the "
                "model itself." if protein else "The uploaded file.")
    if protein and not af:
        doc.add(*doc.explain(yours="AlphaFold DB has no model for this protein (or could not be reached), so there is "
                                   "no confidence track, PAE or predicted 3D model. Experimental structures, if any, are "
                                   "in the next section."))
    if af:
        info = [["Model", f"{af['entry']} version {af.get('version')}" + (f" · {af['tool']}" if af.get("tool") else "")],
                ["Created", af.get("created") or "—"],
                ["Covers", f"UniProt residues {af.get('uniprot_start')}–{af.get('uniprot_end')}"],
                ["Mean pLDDT", f"{_fmt(af.get('mean_plddt'))} (AlphaFold DB reports {_fmt(af.get('global_metric'))})"]]
        t = Table([[P(f"<b>{k}</b>", "cell"), P(v, "body")] for k, v in info], colWidths=[32 * mm, doc.width - 32 * mm])
        t.setStyle(TableStyle([("LINEBELOW", (0, 0), (-1, -1), 0.3, LINE), ("LEFTPADDING", (0, 0), (-1, -1), 0),
                               ("TOPPADDING", (0, 0), (-1, -1), 2.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5)]))
        doc.add(P("The AlphaFold model", "h2"), t, Spacer(1, 8))
        n_conf = (af.get("frac_confident") or 0)
        doc.figure("Confidence (pLDDT) across the model", doc.image(bands_png(doc, af)),
                   how="pLDDT is AlphaFold's own per-residue confidence (0–100) that the local structure around a residue "
                       "is right. Above 90 the backbone and most side chains are usually accurate; 70–90 means a correct "
                       "backbone; 50–70 is low confidence; below 50 usually means the region is disordered on its own, or "
                       "folds only with a partner — it is drawn as loose ribbon, not a real extended structure. The bar "
                       "shows what fraction of residues falls in each band.",
                   yours=f"{100 * n_conf:.0f} % of the {len(af['plddt'])} residues are modelled confidently (pLDDT ≥ 70) "
                         f"and {100 * (af.get('frac_very_low') or 0):.0f} % are very low (< 50). Mean pLDDT "
                         f"{_fmt(af.get('mean_plddt'))}.")
    png = track_png(res, doc.fig_dir / "track.png")
    if png:
        doc.figure("Confidence, domains and variants along the sequence", doc.image(png), how=tx.get("track", {}).get("how"),
                   yours=tx.get("track", {}).get("yours"))
    pae = pae_png(doc, res) if af else None
    if pae:
        doc.figure("Predicted aligned error (PAE)", doc.image(pae, max_w=doc.width * 0.75), how=tx.get("pae", {}).get("how"),
                   yours=tx.get("pae", {}).get("yours"))
    elif af:
        doc.add(*doc.explain(yours="The PAE file could not be downloaded, so how the domains pack against each other "
                                   "cannot be judged from this report."))
    # the 3D model
    main_sid = "AF" if "AF" in res["structures"] else next(iter(res["structures"]), None)
    if main_sid:
        meta = res["structures"][main_sid]
        a, b_ = _snapshot(snaps, "model_front"), _snapshot(snaps, "model_back")
        how3 = ("A picture of the 3D model. " + ("It is drawn as a cartoon (helices as spirals, strands as arrows) "
                "coloured by pLDDT with the same four colours as above, laid along its longest axis and seen from two opposite sides. "
                if a else "") + "Low-confidence stretches (yellow, orange) look like long loose loops: that is how "
                "AlphaFold draws regions it cannot place, not a real extended structure. The interactive viewer on the "
                "page can colour the same model by domain, burial or AlphaMissense score and show single residues.")
        if a:
            from PIL import Image
            import io
            with Image.open(io.BytesIO(a)) as im:
                wide = im.size[0] / im.size[1] > 1.6
            if wide:          # a rod: one picture per row, full width
                t = Table([[doc.image(x, max_h=70 * mm)] for x in (a, b_) if x], colWidths=[doc.width])
            else:
                imgs = [doc.image(x, max_w=doc.width / 2 - 2, max_h=80 * mm) for x in (a, b_) if x]
                t = Table([imgs], colWidths=[doc.width / 2] * len(imgs))
            t.setStyle(TableStyle([("LEFTPADDING", (0, 0), (-1, -1), 0), ("RIGHTPADDING", (0, 0), (-1, -1), 2)]))
            doc.figure("The 3D model" if not protein else "The AlphaFold 3D model", t, how=how3,
                       yours=tx.get("viewer", {}).get("yours"), sub=f"{meta['label']} · snapshot from the 3D viewer")
        else:
            try:
                xyz, num, bs = _ca(doc, meta)
                cols = [_afc(v) for v in bs] if meta.get("bkind") == "plddt" else plt.cm.viridis(np.linspace(0, 1, len(xyz)))
                png3 = trace_png(doc, "model_trace.png", xyz, cols)
                doc.figure("The 3D model" if not protein else "The AlphaFold 3D model", doc.image(png3), how=how3,
                           yours=tx.get("viewer", {}).get("yours"),
                           note="Drawn by the server as the CA trace (one line through the backbone), because this report "
                                "was made without the page's 3D viewer. Use the PDF report button on the page for the "
                                "cartoon picture." + ("" if meta.get("bkind") == "plddt" else
                                                     " Coloured from the first residue (purple) to the last (yellow)."),
                           sub=f"{meta['label']} chain {meta.get('main_chain')}")
            except Exception:  # noqa: BLE001
                log_exc("report trace")
    # chains of every structure in the item
    rows = []
    for sid, m in res["structures"].items():
        for c in m.get("chains", []):
            rows.append([m["label"], c["chain"], c.get("molecule") or "—", str(c.get("n", "—")),
                         f"{c['span'][0]}–{c['span'][1]}" if c.get("span") else "—",
                         _pct(c.get("identity")) if c.get("identity") is not None else "—",
                         "yes" if c.get("mapped") else "no",
                         f"{m.get('method') or ''} {_fmt(m.get('resolution'), 2) + ' Å' if m.get('resolution') else ''}".strip() or "—"])
    if rows:
        doc.table("Structures in this analysis", ["Structure", "Chain", "Molecule", "Residues", "UniProt span",
                                                  "Identity", "Mapped", "Method"], rows,
                  (0.2, 0.06, 0.2, 0.08, 0.12, 0.08, 0.07, 0.19),
                  how="Every structure loaded into this analysis and each of its polymer chains. A chain is 'mapped' when "
                      "its sequence was aligned to the UniProt sequence, so its residues carry UniProt numbers; identity "
                      "is over the aligned residues (a mutant or engineered construct shows below 100 %).",
                  yours=f"{len(res['structures'])} structure{'s' if len(res['structures']) != 1 else ''}, "
                        f"{len(rows)} chain{'s' if len(rows) != 1 else ''}.", max_rows=120)


def pdb_section(doc: Doc, res: dict):
    u = res.get("uniprot")
    if not u:
        return
    E = u.get("pdb") or []
    tx = (res.get("texts") or {}).get("pdb", {})
    doc.section("PDB", "Experimental structures", "Every entry UniProt cross-references, placed on the sequence.")
    if not E:
        doc.add(*doc.explain(tx.get("how"), tx.get("yours") or "UniProt lists no experimental structure for this protein."))
        return
    order = sorted(E, key=lambda e: ((e.get("span") or [0])[0], e.get("resolution") or 99))
    shown = order if len(order) <= 60 else sorted(E, key=lambda e: (-(e.get("covered") or 0), e.get("resolution") or 99))[:60]
    shown = sorted(shown, key=lambda e: ((e.get("span") or [0])[0], e.get("resolution") or 99))
    doc.figure(f"PDB entries along the sequence ({len(E)})", doc.image(coverage_png(doc, u, res.get("alphafold"), shown),
                                                                      max_h=230 * mm),
               how=tx.get("how"), yours=tx.get("yours"),
               note="" if len(shown) == len(E) else f"The figure shows the {len(shown)} entries covering the most residues; "
                                                    f"all {len(E)} are in the table.")
    loaded = {m.get("pdb_id") for m in res["structures"].values() if m.get("pdb_id")}
    doc.table("All entries", ["PDB", "Method", "Resolution", "Chains", "UniProt residues", "Covered", "Loaded"],
              [[e["id"], e["method"], e.get("resolution_text") or "—",
                ", ".join(sorted({c for sg in e["segments"] for c in sg["chains"]}))[:40],
                "; ".join(f"{sg['start']}–{sg['end']}" for sg in e["segments"])[:60], str(e.get("covered") or "—"),
                "yes" if e["id"] in loaded else ""] for e in order],
              (0.08, 0.1, 0.1, 0.2, 0.32, 0.1, 0.1),
              how="The same entries as a table, in sequence order. 'UniProt residues' is the span of the deposited "
                  "construct that UniProt maps to this protein; not every residue in it is necessarily resolved in "
                  "the density. 'Loaded' entries were opened in this analysis and are described in the other sections.")


def variant_section(doc: Doc, res: dict, snaps: dict | None):
    V = res.get("variants")
    if not V:
        return
    doc.section("Variants", f"Variants on {V['structure_label']} chain {V['chain']}",
                "Each variant is checked against the reference sequence, placed on the structure and described.")
    how = V.get("how", "")
    if V.get("ss_method") == "phi/psi":
        how += (" Secondary structure here is an approximation from backbone phi/psi angles (the DSSP program is not "
                "installed): 'extended' includes polyproline and disordered stretches, not only beta-strands.")
    doc.add(*doc.explain(how, V.get("yours")), Spacer(1, 8))
    png = variants_png(doc, res, V)
    if png:
        use_am = any(r.get("am_score") is not None for r in V["rows"] if not r.get("error"))
        s = V.get("summary") or {}
        doc.figure("Where the variants sit", doc.image(png),
                   how=("Each variant as a lollipop at its position. Height is the AlphaMissense pathogenicity score "
                        "(0–1; the red band above 0.564 is 'likely pathogenic', the blue band below 0.34 'likely benign', "
                        "between is ambiguous). " if use_am else
                        "Each variant as a lollipop at its position; height is the model's confidence (pLDDT) or the "
                        "B-factor at that residue. ") + "Colour is burial from relative solvent accessibility: blue "
                       "buried (< 20 %), grey intermediate, orange exposed (> 50 %).",
                   yours=f"{s.get('n', len(V['rows']))} variants: {s.get('buried_confident', 0)} buried in a confident "
                         f"region, {s.get('exposed_confident', 0)} exposed in a confident region, "
                         f"{s.get('low_confidence', 0)} in low-confidence regions"
                         + (f", {s.get('am_pathogenic', 0)} likely pathogenic by AlphaMissense" if use_am else "") + ".")
    snap = _snapshot(snaps, "variants")
    how_v = ("The structure the variants were mapped on, with each mapped variant marked in magenta and labelled"
             + (" (sticks and a sphere on the CA atom, as in the page's viewer with 'variants' on)." if snap else "."))
    placed = [r for r in V["rows"] if r.get("key") and not r.get("unresolved")]
    yours_v = (f"{len(placed)} of {len(V['rows'])} variants are resolved in this structure and marked. Buried ones sit "
               "inside the fold; exposed ones are on the surface, where they matter mostly if a partner or ligand binds "
               "there — see the contacts in the table.")
    if snap:
        doc.figure("Variants on the 3D structure", doc.image(snap, max_h=95 * mm), how=how_v, yours=yours_v,
                   sub=V["structure_label"] + " · snapshot from the 3D viewer")
    elif placed:
        try:
            meta = res["structures"][V["structure"]]
            xyz, num, bs = _ca(doc, meta, V["chain"])
            at = {n: i for i, n in enumerate(num)}
            marks = [(xyz[at[r["num"]]], r["label"]) for r in placed if r.get("num") in at]
            cols = [_afc(v) for v in bs] if meta.get("bkind") == "plddt" else ["#9aa6ad"] * len(xyz)
            png3 = trace_png(doc, "variants_trace.png", xyz, cols, marks=marks)
            doc.figure("Variants on the 3D structure", doc.image(png3), how=how_v, yours=yours_v,
                       note="Drawn by the server as the CA trace, because this report was made without the page's 3D "
                            "viewer.", sub=V["structure_label"])
        except Exception:  # noqa: BLE001
            log_exc("report variant trace")
    # which column matches the structure the variants are actually mapped on (variants_png's y-axis label uses the
    # same test) — 'protein' must not force 'pLDDT' here, or a crystal structure's own B-factor never gets shown
    bcol = "pLDDT" if V.get("bkind") == "plddt" else "B-factor"
    rows = []
    for r in V["rows"]:
        if r.get("error"):
            rows.append([f"<b>{r.get('label') or r['input']}</b>", "", "", "", "", "", "", r["error"]])
            continue
        b = r.get("plddt") if bcol == "pLDDT" else r.get("bfactor")
        nb = r.get("neighbours") or {}
        rows.append([f"<b>{r['label']}</b>",
                     str(r["pos"]) + (f" ({r['chain']}:{r['num']})" if r.get("num") not in (None, r["pos"]) else ""),
                     _fmt(b, 0),
                     "not modelled" if r.get("unresolved") else f"{_pct(r.get('rsa'))} {r.get('rsa_class') or ''}",
                     r.get("ss_word") or "—",
                     "—" if r.get("unresolved") or not nb else f"{nb.get('n', 0)} ({nb.get('nonlocal', 0)} non-local)",
                     "; ".join((r.get("domain") or [])[:2]) or "—",
                     "—" if r.get("am_score") is None else f"{r['am_score']:.2f} {r.get('am_class')}"])
    doc.table("Summary table", ["Variant", "Position", bcol, "Rel. SASA", "Sec. str.", "Contacts ≤ 5 Å", "Domain / region",
                                "AlphaMissense"], rows, (0.09, 0.09, 0.06, 0.13, 0.09, 0.12, 0.24, 0.18),
              how="One row per variant. Rel. SASA is the residue's solvent-accessible surface divided by the maximum for "
                  "that amino acid (Shrake–Rupley, Tien et al. 2013 maxima). Contacts are residues with any heavy atom "
                  "within 5 Å; 'non-local' ones are more than 4 residues away in sequence, i.e. packing contacts.")
    doc.add(P("Each variant in detail", "h2"),
            P("Everything the page shows for each variant, with the reading it gives. Readings are heuristics built from "
              "burial, confidence and contacts — not predictions of pathogenicity.", "how"), Spacer(1, 4))
    for r in V["rows"]:
        if r.get("error"):
            doc.add(P(f"<b>{r.get('label') or r['input']}</b> — {r['error']}", "body"), Spacer(1, 4))
            continue
        nb = r.get("neighbours") or {}
        lines = [f"<b>{r['label']}</b> · residue {r['pos']} · {r.get('kind', '')}"]
        if r.get("change"):
            lines.append("Change: " + r["change"])
        if nb and not r.get("unresolved"):
            lines.append("Contacts within 5 Å: " + (", ".join(nb.get("residues", [])[:14]) or "none") +
                         (f" · other chains {', '.join(nb['other_chain'])}" if nb.get("other_chain") else "") +
                         (f" · ligand {', '.join(nb['ligand'])}" if nb.get("ligand") else ""))
        if r.get("uniprot_known"):
            lines.append("UniProt: " + " | ".join(re.sub(r";\s*dbSNP:\S+", "", k["description"]) for k in r["uniprot_known"]))
        elif r.get("uniprot_other_at_pos"):
            lines.append("UniProt lists other changes at this position: " + ", ".join(r["uniprot_other_at_pos"]))
        for f in r.get("flags", []):
            tag = {"ok": "OK", "warn": "CHECK", "info": "NOTE"}[f["level"]]
            lines.append(f"<b>{tag}</b> {f['text']}")
        # one paragraph per line: clean() keeps only <b>/<i>, so a <br/> between them was dropped and the lines ran together
        doc.add(KeepTogether([P(lines[0], "body")] + [P(x, "cell") for x in lines[1:]] + [Spacer(1, 7)]))
    doc.add(P("AlphaMissense (Cheng et al. 2023) is a prediction for the canonical human isoform, not a clinical "
              "classification.", "sub"))


def comparison_section(doc: Doc, res: dict, snaps: dict | None):
    cmps = res.get("comparisons") or []
    if not cmps:
        return
    doc.section("Compare", "Structure comparisons", "Two structures superposed, and where they differ along the sequence.")
    for k, c in enumerate(cmps):
        doc.add(CondPageBreak(120 * mm), P(f"{c['label_b']} chain {c['chain_b']} superposed on {c['label_a']} chain "
                                           f"{c['chain_a']}", "h2"))
        tiles = [["CE RMSD (structural core)", f"{_fmt((c.get('ce') or {}).get('rmsd'), 2)} Å"],
                 ["RMSD on the common core", f"{_fmt(c.get('rmsd_core') if c.get('rmsd_core') is not None else c.get('rmsd_mapped'), 2)} Å "
                                             f"({c.get('n_core') or c.get('n_pairs')} residues)"],
                 ["RMSD on all matched pairs", f"{_fmt(c.get('rmsd_mapped'), 2)} Å ({c.get('n_pairs')} pairs)"],
                 ["Stretches deviating > 3 Å", str(len(c.get("segments") or []))]]
        t = Table([[P(f"<b>{a}</b>", "cell"), P(b, "body")] for a, b in tiles], colWidths=[55 * mm, doc.width - 55 * mm])
        t.setStyle(TableStyle([("LINEBELOW", (0, 0), (-1, -1), 0.3, LINE), ("LEFTPADDING", (0, 0), (-1, -1), 0),
                               ("TOPPADDING", (0, 0), (-1, -1), 2.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5)]))
        doc.add(t, Spacer(1, 4))
        png = deviation_png(doc, c, k)
        if png:
            doc.figure("CA deviation along the sequence", doc.image(png), how=c.get("how"), yours=c.get("yours"))
        else:
            doc.add(*doc.explain(c.get("how"), c.get("yours")))
        snap = _snapshot(snaps, f"cmp:{c['key']}")
        how_s = ("Both structures after superposition: the reference in grey, the moved structure in teal. Where the two "
                 "lines separate is where the deviation plot above rises.")
        if snap:
            doc.figure("The superposition in 3D", doc.image(snap, max_h=90 * mm), how=how_s,
                       yours=f"{c['label_b']} chain {c['chain_b']} (teal) on {c['label_a']} chain {c['chain_a']} (grey), "
                             f"zoomed on the moved chain; {len(c.get('segments') or [])} stretch"
                             f"{'es' if len(c.get('segments') or []) != 1 else ''} deviate by more than 3 Å.",
                       sub="snapshot from the 3D viewer")
        else:
            try:
                ma, mb = res["structures"][c["a"]], dict(res["structures"][c["b"]])
                mb["file"] = c["superposed_file"]
                xa, na, _ = _ca(doc, ma, c["chain_a"])
                xb, _, _ = _ca(doc, mb, c["chain_b"])
                # only the compared stretch of the reference (±8 residues): a small domain superposed on a
                # full-length model is otherwise a speck beside it
                an = [int(re.sub(r"\D", "", d_["a"].split(":")[1]) or 0) for d_ in c.get("deviation") or [] if d_.get("a")]
                crop = ""
                if an:
                    lo, hi = min(an) - 8, max(an) + 8
                    keep = [i for i, n_ in enumerate(na) if lo <= n_ <= hi]
                    if 0 < len(keep) < len(na):
                        xa = xa[keep]
                        crop = f" The reference is cropped to the compared region (its residues {lo}–{hi})."
                png3 = trace_png(doc, f"superposition_{k}.png", xb, ["#0f766e"] * len(xb), overlay=(xa, "#b8c2c8"))
                doc.figure("The superposition in 3D", doc.image(png3), how=how_s,
                           yours=f"{c['label_b']} chain {c['chain_b']} (teal) on {c['label_a']} chain {c['chain_a']} "
                                 f"(grey).{crop}",
                           note="Drawn by the server as CA traces, because this report was made without the page's 3D "
                                "viewer.")
            except Exception:  # noqa: BLE001
                log_exc("report superposition")
        if c.get("segments"):
            doc.table("", ["From", "To", "Residues", "Max deviation (Å)"],
                      [[str(s_["start"]), str(s_["end"]), str(s_["end"] - s_["start"] + 1), _fmt(s_.get("max_dev"), 1)]
                       for s_ in c["segments"]], (0.2, 0.2, 0.2, 0.4),
                      how="Stretches of at least a few consecutive residues deviating by more than 3 Å.")


def methods_section(doc: Doc, res: dict):
    doc.section("Reproducibility", "Methods and references", "What was computed, with which Biopython function, and the "
                                                              "citation for each method.")
    for h, p in res.get("methods", []):
        doc.add(P(h, "h3"), P(p, "body"))
    if res.get("versions"):
        doc.add(Spacer(1, 8), P("Software versions", "h2"),
                P(", ".join(f"{k} {v}" for k, v in res["versions"].items()), "muted"))
    src = res.get("sources") or {}
    if src:
        doc.add(Spacer(1, 6), P("Where the data came from", "h2"),
                P("; ".join(f"{k}: {'downloaded' if v == 'network' else 'local cache' if v == 'cache' else v}"
                            for k, v in src.items()), "muted"))
    doc.add(Spacer(1, 8), P("References", "h2"))
    for c in core.citations(res):
        doc.add(P(c, "muted"))
    doc.add(Spacer(1, 6), P("Data: UniProt and the AlphaFold Protein Structure Database are distributed under CC BY 4.0; "
                            "PDB entries are in the public domain (CC0).", "sub"),
            Spacer(1, 10), P("Interpretations are generated automatically from the numbers above and are meant to guide, "
                             "not replace, expert review.", "sub"))


def build_pdf(item: str, snapshots: dict | None = None) -> Path:
    d = item_dir(item)
    res = load_result(item)
    out = d / "report.pdf"
    W, H = A4
    margin = 16 * mm
    doc = Doc(d, W - 2 * margin)
    protein = res["kind"] == "protein"
    cover(doc, res)
    if protein:
        protein_section(doc, res)
    structure_section(doc, res, snapshots)
    if protein:
        pdb_section(doc, res)
    variant_section(doc, res, snapshots)
    comparison_section(doc, res, snapshots)
    methods_section(doc, res)

    def deco(canvas, d_):
        canvas.saveState()
        canvas.setFont("Helvetica", 7.5)
        canvas.setFillColor(MUTED)
        canvas.drawString(margin, 10 * mm, f"Structure Bench · {res.get('name')} · created by Paul H. Kim, Ph.D.")
        canvas.drawRightString(W - margin, 10 * mm, f"{d_.page}")
        canvas.setStrokeColor(ACCENT)
        canvas.setLineWidth(2)
        canvas.line(margin, H - 10 * mm, margin + 18 * mm, H - 10 * mm)
        canvas.restoreState()

    pdf = SimpleDocTemplate(str(out), pagesize=A4, leftMargin=margin, rightMargin=margin, topMargin=16 * mm,
                            bottomMargin=16 * mm, title=f"Structure Bench report — {res.get('name')}",
                            author="Paul H. Kim, Ph.D.")
    pdf.build(doc.story, onFirstPage=deco, onLaterPages=deco)
    return out
