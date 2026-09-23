# Vendored front-end libraries

| File | Library | Version | Source | Licence |
|---|---|---|---|---|
| `3Dmol-min.js` | 3Dmol.js | **2.5.5** | https://cdn.jsdelivr.net/npm/3dmol@2.5.5/build/3Dmol-min.js (downloaded once, 2026-09-23) | BSD 3-Clause — `3Dmol-LICENSE.txt` |

- sha256 of `3Dmol-min.js`: `f7cc78921ae72e7623e89cdd111434f58c2efddd2ffda1cd212644b406fb8016`
- `3Dmol-min.js.LICENSE.txt` is the licence banner the minified build refers to; `3Dmol-LICENSE.txt` is the
  package's full LICENSE (3Dmol.js incorporates GLmol, Three.js and jQuery code under compatible licences).
- Loaded by `web/index.html` as `vendor/3Dmol-min.js` — never from a CDN at run time, so the viewer works offline.
- Cite: Rego N, Koes D. 3Dmol.js: molecular visualization with WebGL. *Bioinformatics* 31:1322–1324 (2015).

To upgrade: download the new `build/3Dmol-min.js` for a pinned version from jsDelivr, update the table and hash,
and click through the viewer (colour modes, surface, click-to-inspect, PNG export).
