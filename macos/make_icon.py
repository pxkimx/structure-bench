#!/usr/bin/env python3
"""Draw the app icon from the same geometry as the logo in web/index.html, and build AppIcon.icns.

    python macos/make_icon.py

The coordinates are the SVG's 32-unit viewBox — a rounded panel, a protein chain trace (three cubic Béziers)
and five beads in the AlphaFold confidence colours plus the app's teal — scaled onto Apple's 1024 canvas,
so the icon and the in-app mark cannot drift apart.
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from PIL import Image, ImageDraw

HERE = Path(__file__).resolve().parent
RES = HERE / "Structure Bench.app" / "Contents" / "Resources"
SS = 4
CANVAS, INSET, BOX = 1024, 100, 824
PANEL, TRACE = "#0d1b22", "#6a7b84"
# path d="M5.5 21.5C8 11 11.5 9 14 15.5S20 24 22.5 17 25.5 8.5 27 10" with the S segments expanded
CURVES = [((5.5, 21.5), (8, 11), (11.5, 9), (14, 15.5)),
          ((14, 15.5), (16.5, 22), (20, 24), (22.5, 17)),
          ((22.5, 17), (25, 10), (25.5, 8.5), (27, 10))]
BEADS = [(5.5, 21.5, 2.1, "#FF7D45"), (10.2, 11.6, 2.3, "#65CBF3"), (15, 17.4, 2.5, "#0053D6"),
         (20.4, 20.6, 2.3, "#00e0cf"), (24.4, 11.8, 2.1, "#FFDB13")]


def bezier(p0, p1, p2, p3, n=60):
    out = []
    for i in range(n + 1):
        t = i / n
        a, b, c, d = (1 - t) ** 3, 3 * (1 - t) ** 2 * t, 3 * (1 - t) * t ** 2, t ** 3
        out.append((a * p0[0] + b * p1[0] + c * p2[0] + d * p3[0], a * p0[1] + b * p1[1] + c * p2[1] + d * p3[1]))
    return out


def main():
    N = CANVAS * SS
    im = Image.new("RGBA", (N, N), (0, 0, 0, 0))
    d = ImageDraw.Draw(im)
    u = BOX * SS / 32.0
    o = INSET * SS
    P = lambda x, y: (o + x * u, o + y * u)  # noqa: E731
    d.rounded_rectangle([P(0.5, 0.5), P(31.5, 31.5)], radius=7.5 * u, fill=PANEL)
    pts = [P(*q) for c in CURVES for q in bezier(*c, n=400)]
    r = 1.3 * u / 2
    for x, y in pts:                                       # a stroked path with round joins and caps
        d.ellipse([x - r, y - r, x + r, y + r], fill=TRACE)
    for x, y, r, c in BEADS:
        d.ellipse([P(x - r, y - r), P(x + r, y + r)], fill=c)
    im = im.resize((CANVAS, CANVAS), Image.LANCZOS)
    RES.mkdir(parents=True, exist_ok=True)
    im.save(RES / "AppIcon.png")
    if shutil.which("iconutil"):
        with tempfile.TemporaryDirectory() as td:
            iset = Path(td) / "AppIcon.iconset"
            iset.mkdir()
            for s in (16, 32, 128, 256, 512):
                im.resize((s, s), Image.LANCZOS).save(iset / f"icon_{s}x{s}.png")
                im.resize((2 * s, 2 * s), Image.LANCZOS).save(iset / f"icon_{s}x{s}@2x.png")
            subprocess.run(["iconutil", "-c", "icns", str(iset), "-o", str(RES / "AppIcon.icns")], check=True)
    print("wrote", RES / "AppIcon.png")


if __name__ == "__main__":
    main()
