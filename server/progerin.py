"""Progerin (HGPS) mode: what the Hutchinson–Gilford progeria mutation does to lamin A, on the sequence and the model.

Every coordinate is read from the UniProt record of the protein being shown (P02545, human prelamin-A/C) at run
time — the deleted segment from the "Progerin" isoform's alternative-sequence feature, the cleavage site, the
farnesylated cysteine, the mature chain, the lamin C isoform, the domains — so the numbers on the page are the
record's, not this file's. The few facts UniProt does not hold (the cDNA change and its ClinVar record, the mouse
G609G position, the papers) are written here with their source and the date they were read. `EXPECTED` is what the
literature says those record-derived numbers should be; the tests check the two agree.
"""
from __future__ import annotations

import re

HUMAN_LMNA = "P02545"
VERIFIED = "2026-11-22"      # the day every fact below was read from its source (UniProt REST, ClinVar E-utilities, PubMed)

# What the record is expected to say (Eriksson 2003; De Sandre-Giovannoli 2003; UniProt P02545). Tests compare the
# values derived from the UniProt record with these; the page shows the record's values.
EXPECTED = {"prelamin_length": 664, "deletion": (607, 656), "deletion_length": 50, "progerin_length": 614,
            "cleavage": (646, 647), "farnesyl_cys": 661, "caax": "CSIM", "mature_end": 646,
            "lamin_c_shared_end": 566, "lamin_c_tail": "VSGSRR", "lamin_c_length": 572,
            "nls": (417, 422), "ig_fold": (428, 545), "rod": (31, 387)}

# Facts that are not in the UniProt record, each with where it was read.
MUTATION = {
    "cdna": "c.1824C>T", "protein": "p.Gly608= (G608G)", "transcript": "NM_170707.4", "exon": 11,
    "clinvar": "VCV000014500", "clinvar_id": "14500", "dbsnp": "rs58596362",
    "classification": "Pathogenic — criteria provided, multiple submitters, no conflicts (last evaluated 2025-04-30)",
    "conditions": "Hutchinson-Gilford progeria syndrome",
    "effect": ("A silent change (GGC>GGT, glycine either way) that activates a cryptic splice donor site within exon 11. Part of "
               "exon 11 is then spliced out of the mRNA: an in-frame loss of 150 nucleotides, 50 amino acids near the C-terminus "
               "(Eriksson 2003; UniProt isoform 6 note)."),
    "source": f"ClinVar {'VCV000014500'} via NCBI E-utilities, read {VERIFIED}; Eriksson et al. 2003; De Sandre-Giovannoli et al. 2003",
}
OTHER_VARIANTS = [
    {"label": "G608S", "cdna": "c.1822G>A", "clinvar": "VCV000014501", "dbsnp": "rs61064130",
     "classification": "Pathogenic — no assertion criteria provided (2003-05-15)",
     "note": ("Changes the same codon. UniProt VAR_017664: in HGPS; reduced binding to SUN1; may affect splicing by activating "
              "a cryptic splice donor site (Eriksson 2003). It produces progerin by the same mechanism."),
     "source": f"ClinVar VCV000014501 and UniProt P02545 VAR_017664, read {VERIFIED}"},
]
MOUSE = {
    "accession": "P48678", "gene": "Lmna", "organism": "Mus musculus", "length": 665,
    "human_pos": 608, "mouse_pos": 609, "residue": "G", "offset": 1,
    "human_window": "SASGSGAQVGGPISSG", "human_window_range": (599, 614),
    "mouse_window": "AAGGAGAQVGGSISSG", "mouse_window_range": (600, 615),
    "identity": "12 of 16", "farnesyl_cys": 662, "c_terminus": "CSIM",
    "model": "Lmna G609G knock-in: the mouse strain that carries the HGPS mutation, of Osorio et al. 2011",
    "source": f"UniProt P48678 entry version 226 (annotation 2026-11-04), sequence read {VERIFIED}; Osorio et al. 2011 (PMID 22030750)",
}
PAPERS = {
    "eriksson2003": "Eriksson M, et al. Recurrent de novo point mutations in lamin A cause Hutchinson–Gilford progeria syndrome. Nature 423:293–298 (2003). PMID 12714972.",
    "desandre2003": "De Sandre-Giovannoli A, et al. Lamin A truncation in Hutchinson–Gilford progeria. Science 300:2055 (2003). PMID 12702809.",
    "osorio2011": "Osorio FG, et al. Splicing-directed therapy in a new mouse model of human accelerated aging. Sci Transl Med 3:106ra107 (2011). PMID 22030750.",
    "gordon2012": "Gordon LB, et al. Clinical trial of a farnesyltransferase inhibitor in children with Hutchinson–Gilford progeria syndrome. PNAS 109:16666–16671 (2012). PMID 23012407.",
    "dhillon2021": "Dhillon S. Lonafarnib: first approval. Drugs 81:283–289 (2021). PMID 33590450.",
    "koblan2021": "Koblan LW, et al. In vivo base editing rescues Hutchinson–Gilford progeria syndrome in mice. Nature 589:608–614 (2021). PMID 33408413.",
    "clinvar": "Landrum MJ, et al. ClinVar: improvements to accessing data. Nucleic Acids Res 48:D835–D844 (2020). Record VCV000014500 (NM_170707.4(LMNA):c.1824C>T, p.Gly608=).",
}

LANDMARK_NAMES = {"Head": "head", "Coil 1A": "coil", "Linker 1": "linker", "Coil 1B": "coil", "Linker 2": "linker",
                  "Coil 2": "coil", "Tail": "tail"}


def _feat(u: dict, ftype: str, pattern: str):
    rx = re.compile(pattern, re.I)
    return [f for f in u.get("features", []) if f["type"] == ftype and rx.search(f.get("description") or "")]


def _isoform(u: dict, synonym: str):
    for iso in u.get("isoforms") or []:
        if any(synonym.lower() == s.lower() for s in iso.get("synonyms") or []):
            return iso
    return None


def _vsp(u: dict, ids: list) -> list:
    ids = set(ids or [])
    return sorted([f for f in u.get("features", []) if f["type"] == "Alternative sequence" and f.get("id") in ids],
                  key=lambda f: f["start"])


def _mean(vals):
    vals = [v for v in vals if v is not None]
    return round(sum(vals) / len(vals), 1) if vals else None


def block(u: dict | None, af: dict | None = None) -> dict | None:
    """The progerin facts for a UniProt entry, or None when the entry is not human LMNA."""
    if not u or (u.get("accession") or "").upper() != HUMAN_LMNA:
        return None
    seq = u.get("sequence") or ""
    L = len(seq)
    notes = []
    # --- the deleted segment: the "Progerin" isoform's alternative-sequence feature
    iso = _isoform(u, "Progerin")
    vsp = _vsp(u, iso.get("vsp") if iso else [])
    dele = [f for f in vsp if not f.get("alt")]               # a deletion (no replacement sequence)
    if dele:
        ds, de = dele[0]["start"], dele[0]["end"]
        del_src = f"UniProt {u['accession']} feature {dele[0].get('id', '')} — isoform {iso['name']} ({', '.join(iso.get('ids') or [])}, synonym Progerin): residues {ds}–{de} missing"
    else:
        ds, de = EXPECTED["deletion"]
        del_src = "expected from Eriksson 2003 / De Sandre-Giovannoli 2003 — the UniProt record read here has no Progerin isoform feature"
        notes.append("The UniProt record loaded here does not carry the Progerin isoform feature, so the deletion is taken from the literature.")
    dlen = de - ds + 1
    # --- cleavage, farnesylation, propeptide, mature chain
    cl = _feat(u, "Site", r"cleavage")
    cl = [f for f in cl if f["end"] - f["start"] == 1]
    if cl:
        c_after, c_before = cl[0]["start"], cl[0]["end"]
        cl_src = f"UniProt site {c_after}–{c_before} “{cl[0]['description']}”; the entry's PTM note: the last 15 residues are removed by ZMPSTE24"
    else:
        c_after, c_before = EXPECTED["cleavage"]
        cl_src = "expected (UniProt P02545) — no cleavage site feature in the record read here"
        notes.append("No cleavage-site feature in the record loaded here; the ZMPSTE24 site is taken from the literature.")
    lip = _feat(u, "Lipidation", r"farnesyl")
    farn = lip[0]["start"] if lip else EXPECTED["farnesyl_cys"]
    meth = _feat(u, "Modified residue", r"methyl ester")
    meth_pos = meth[0]["start"] if meth else None
    prop = [f for f in u.get("features", []) if f["type"] == "Propeptide" and f["end"] == L]
    prop_rng = (prop[0]["start"], prop[0]["end"]) if prop else (farn + 1, L)
    mature = [f for f in u.get("features", []) if f["type"] == "Chain" and re.search(r"^Lamin-A/C$", f.get("description") or "")]
    m_end = mature[0]["end"] if mature else EXPECTED["mature_end"]
    # --- lamin C
    isoc = _isoform(u, "Lamin-C")
    vc = _vsp(u, isoc.get("vsp") if isoc else [])
    rep = [f for f in vc if f.get("alt")]
    cut = [f for f in vc if not f.get("alt")]
    if rep and cut:
        c_shared, c_tail = rep[0]["start"] - 1, (rep[0]["alt"] or [""])[0]
        lamin_c_src = f"UniProt {rep[0].get('id', '')} ({rep[0]['start']}–{rep[0]['end']} {seq[rep[0]['start'] - 1:rep[0]['end']]}→{c_tail}) and {cut[0].get('id', '')} ({cut[0]['start']}–{cut[0]['end']} missing), isoform {isoc['name']} ({', '.join(isoc.get('ids') or [])})"
    else:
        c_shared, c_tail = EXPECTED["lamin_c_shared_end"], EXPECTED["lamin_c_tail"]
        lamin_c_src = "expected (UniProt isoform C) — not in the record read here"
    # --- landmarks
    land = []
    for f in u.get("features", []):
        d = f.get("description") or ""
        if f["type"] == "Region" and d in LANDMARK_NAMES:
            land.append({"name": d, "start": f["start"], "end": f["end"], "kind": LANDMARK_NAMES[d]})
        elif f["type"] == "Motif" and "nuclear localization" in d.lower():
            land.append({"name": "NLS", "start": f["start"], "end": f["end"], "kind": "motif", "full": d})
        elif f["type"] == "Domain" and d.upper() == "LTD":
            land.append({"name": "Ig fold (LTD)", "start": f["start"], "end": f["end"], "kind": "domain",
                         "full": "Lamin tail domain, the immunoglobulin-like fold"})
        elif f["type"] == "Domain" and "IF rod" in d:
            land.append({"name": "IF rod", "start": f["start"], "end": f["end"], "kind": "rod", "full": d})
    land.sort(key=lambda x: (x["start"], -(x["end"] - x["start"])))
    # --- the progerin sequence
    prog_seq = seq[:ds - 1] + seq[de:]
    caax = seq[farn - 1:farn + 3]
    # --- UniProt natural variants in or near the deletion, and the HGPS ones anywhere
    def lab(v):
        return f"{v.get('ref', '')}{v['start']}{'/'.join(v.get('alt') or []) or '?'}"
    in_del = [{"label": lab(v), "pos": v["start"], "description": v["description"].split("; dbSNP")[0]}
              for v in u.get("variants", []) if ds <= v["start"] <= de]
    hgps = [{"label": lab(v), "pos": v["start"], "description": v["description"].split("; dbSNP")[0],
             "relation": relation(v["start"], ds, de, dlen)}
            for v in u.get("variants", []) if re.search(r"\bHGPS\b", v.get("description") or "")]
    # --- AlphaFold confidence of the parts in question
    conf = None
    if af and af.get("plddt"):
        pl = af["plddt"]
        ig = next((x for x in land if x["kind"] == "domain"), None)
        tail = next((x for x in land if x["kind"] == "tail"), None)
        conf = {"deleted": _mean(pl[ds - 1:de]), "cleavage": _mean(pl[c_after - 1:c_before]),
                "caax": _mean(pl[farn - 1:L]), "ig_fold": _mean(pl[ig["start"] - 1:ig["end"]]) if ig else None,
                "tail": _mean(pl[tail["start"] - 1:tail["end"]]) if tail else None,
                "deleted_confident": sum(1 for v in pl[ds - 1:de] if v is not None and v >= 70),
                "entry": af.get("entry"), "version": af.get("version")}
    return {
        "accession": u["accession"], "gene": u.get("gene"), "prelamin_length": L, "verified": VERIFIED,
        "mutation": MUTATION, "other_variants": OTHER_VARIANTS,
        "deletion": {"start": ds, "end": de, "length": dlen, "sequence": seq[ds - 1:de], "source": del_src,
                     "from_record": bool(dele)},
        "progerin": {"length": len(prog_seq), "sequence": prog_seq, "c_terminus": prog_seq[-4:],
                     "farnesyl_cys": farn - dlen if farn > de else farn, "isoform": ", ".join(iso.get("ids") or []) if iso else ""},
        "mature": {"start": 1, "end": m_end, "source": f"UniProt chain “Lamin-A/C” 1–{m_end}" if mature else "expected"},
        "cleavage": {"after": c_after, "before": c_before, "residues": seq[c_after - 1:c_before],
                     "inside_deletion": ds <= c_after and c_before <= de, "source": cl_src},
        "caax": {"start": farn, "end": farn + 3, "sequence": caax, "farnesyl": farn, "methyl_ester": meth_pos,
                 "propeptide": list(prop_rng), "kept_in_progerin": not (ds <= farn <= de),
                 "source": (f"UniProt lipidation {farn} “{lip[0]['description']}”" if lip else "expected")
                           + (f"; modified residue {meth_pos} “{meth[0]['description']}”" if meth else "")
                           + (f"; propeptide {prop_rng[0]}–{prop_rng[1]} “{prop[0]['description']}”" if prop else "")},
        "lamin_c": {"shared_end": c_shared, "tail": c_tail, "length": c_shared + len(c_tail), "ids": isoc.get("ids") if isoc else [],
                    "source": lamin_c_src},
        "landmarks": land, "mouse": MOUSE,
        "uniprot_variants_in_deletion": in_del, "hgps_variants": hgps,
        "confidence": conf, "notes": notes, "texts": texts(ds, de, dlen, c_after, c_before, farn, caax, conf, c_shared, c_tail),
        "sources": [f"UniProt {u['accession']} entry version {u.get('version', '?')} (updated {u.get('last_updated', '?')}): "
                    "isoform and feature coordinates", MUTATION["source"], MOUSE["source"],
                    *[PAPERS[k] for k in ("eriksson2003", "desandre2003", "osorio2011", "gordon2012", "dhillon2021", "koblan2021")]],
    }


def relation(pos: int, ds: int, de: int, dlen: int) -> dict:
    """Where a prelamin A position stands relative to the HGPS deletion, and its number in progerin."""
    if pos < ds:
        return {"where": "before", "text": f"N-terminal of the deleted segment — present in progerin at the same position ({pos})", "progerin_pos": pos}
    if pos > de:
        return {"where": "after", "text": f"C-terminal of the deleted segment — present in progerin, renumbered {pos - dlen}", "progerin_pos": pos - dlen}
    return {"where": "inside", "text": f"inside the deleted segment {ds}–{de} — this residue does not exist in progerin", "progerin_pos": None}


def texts(ds, de, dlen, c_after, c_before, farn, caax, conf, c_shared, c_tail) -> dict:
    why = (f"Mature lamin A is made from prelamin A in four steps (UniProt's processing note). Farnesyltransferase puts a farnesyl "
           f"lipid on the cysteine of the C-terminal {caax} motif (C{farn}); the last three residues ({caax[1:]}) are removed by RCE1 "
           f"and/or ZMPSTE24; ICMT methylates the new C-terminal cysteine; and finally ZMPSTE24 cuts between {c_after} and "
           f"{c_before}, removing the last 15 residues together with the farnesylated cysteine. That last cut needs the "
           f"farnesylation and methylation to have happened first. In HGPS the silent change c.1824C>T (G608G) activates a cryptic "
           f"splice donor within exon 11, so residues {ds}–{de} ({dlen} amino acids) are missing from the protein — that protein is "
           f"progerin. The {caax} motif is still there, so progerin is farnesylated and methylated like prelamin A; but the "
           f"ZMPSTE24 site ({c_after}|{c_before}) was inside the deleted stretch, so the farnesylated tail is never removed. "
           "Progerin stays permanently farnesylated — farnesylation is what targets prelamin A to the nuclear envelope — and HGPS "
           "cells show visibly abnormal nuclear membranes (Eriksson 2003).")
    lona = ("Lonafarnib is a farnesyltransferase inhibitor: it blocks the first of those steps, preventing the farnesylation and "
            "the accumulation of progerin that follows (Dhillon 2021). In the first clinical trial (Gordon et al. 2012; 25 children "
            "treated for at least two years) 9 increased their rate of weight gain by half or more, and there was preliminary "
            "evidence of improved vascular stiffness, bone structure and hearing. In November 2020 it was approved in the USA "
            "(Zokinvy) to reduce the risk of death in HGPS and for processing-deficient progeroid laminopathies (Dhillon 2021). It "
            "does not touch the mutation or the mis-splicing — the shortened protein is still made — which is why "
            "splice-blocking antisense morpholinos (first shown in the G609G mouse, Osorio et al. 2011) and base editing of the "
            "c.1824C>T change itself (Koblan et al. 2021, in patient fibroblasts and a mouse model) are being pursued.")
    lc = (f"Lamin C, the other main LMNA product, is identical to lamin A up to residue {c_shared}, then ends in {len(c_tail)} "
          f"different residues ({c_tail}) with no CaaX motif. It is never farnesylated and is untouched by the HGPS deletion.")
    if conf and conf.get("deleted") is not None:
        cav = (f"The AlphaFold model places the C-terminal tail, but with low confidence: mean pLDDT {conf['deleted']} over the "
               f"deleted segment {ds}–{de} ({conf['deleted_confident']} of {dlen} residues at pLDDT ≥ 70) and {conf['tail']} over the "
               f"whole tail, against {conf['ig_fold']} in the Ig fold. The highlight therefore shows where the segment is and roughly "
               "where it sits, not a reliable fold (the PAE plot in the Predicted structure section shows how loosely the tail is "
               f"placed against the rest). No structure of progerin itself is predicted here: removing {dlen} residues from a "
               "low-confidence region does not give a model to trust.")
    else:
        cav = ("No AlphaFold model is loaded, so the segment is shown on the sequence only. The prelamin A tail beyond the Ig "
               "fold has low confidence in the AlphaFold model; a progerin structure is not predicted here.")
    return {"why": why, "lonafarnib": lona, "lamin_c": lc, "model_caveat": cav}


def fasta(b: dict) -> str:
    p = b["progerin"]
    hdr = (f">progerin|{b['accession']}-6|Homo sapiens LMNA isoform 6 (progerin): prelamin A with residues "
           f"{b['deletion']['start']}-{b['deletion']['end']} deleted (HGPS {b['mutation']['cdna']}, {b['mutation']['protein']}); "
           f"{p['length']} aa; C-terminal {p['c_terminus']} kept, farnesyl cysteine at progerin position {p['farnesyl_cys']}. "
           f"Structure Bench, from UniProt {b['accession']}")
    body = "\n".join(p["sequence"][i:i + 60] for i in range(0, len(p["sequence"]), 60))
    return hdr + "\n" + body + "\n"


def pymol_script(b: dict, meta: dict, rows: list[dict], name: str, version: str) -> tuple[str, str]:
    """A PyMOL script that colours the deleted segment, the ZMPSTE24 site, the farnesyl cysteine, the NLS and the Ig fold
    on whichever structure is chosen. `rows` map UniProt positions to the structure's residue numbers."""
    from pathlib import Path
    by_pos = {}
    for r in rows:
        if r.get("pos") and r.get("key"):
            by_pos.setdefault(r["pos"], r["key"])

    def sel(s, e):
        keys = [by_pos[p] for p in range(s, e + 1) if p in by_pos]
        if not keys:
            return None
        chains = {}
        for k in keys:
            ch, num = k.split(":")
            chains.setdefault(ch, []).append(num)
        return " or ".join(f"(chain {ch} and resi {'+'.join(nums)})" for ch, nums in chains.items())

    fname = Path(meta["file"]).name
    obj = "sb_" + "".join(ch if ch.isalnum() else "_" for ch in Path(fname).stem)
    d = b["deletion"]
    L = [f"# PyMOL script from Structure Bench {version} — {name} — progerin (HGPS) mode on {meta['label']}",
         f"# Put this file next to {fname} (download it from the app).",
         f"# Residues are UniProt {b['accession']} positions mapped to this structure's numbering.",
         f"load {fname}, {obj}", "hide everything", "show cartoon", "set cartoon_transparency, 0",
         f"color grey80, {obj}", "bg_color white", "set label_size, 18", "set label_color, black"]
    for lm in b["landmarks"]:
        if lm["kind"] in ("domain", "motif"):
            s = sel(lm["start"], lm["end"])
            if s:
                L.append(f"color {'teal' if lm['kind'] == 'domain' else 'purple'}, {obj} and ({s})   # {lm['name']} {lm['start']}–{lm['end']}")
    s = sel(d["start"], d["end"])
    if s:
        L += [f"select hgps_deleted, {obj} and ({s})   # residues {d['start']}–{d['end']} missing in progerin",
              "color orange, hgps_deleted", "show cartoon, hgps_deleted"]
    c = b["cleavage"]
    s = sel(c["after"], c["before"])
    if s:
        L += [f"select zmpste24_site, {obj} and ({s})", "show sticks, zmpste24_site and not name N+C+O", "color red, zmpste24_site and elem C",
              f"label zmpste24_site and name CA and resi {by_pos.get(c['after'], 'x').split(':')[-1]}, \"ZMPSTE24 site {c['after']}|{c['before']}\""]
    s = sel(b["caax"]["farnesyl"], b["caax"]["farnesyl"])
    if s:
        L += [f"select farnesyl_cys, {obj} and ({s})", "show spheres, farnesyl_cys and name SG", "color yellow, farnesyl_cys and elem S",
              f"label farnesyl_cys and name CA, \"C{b['caax']['farnesyl']} farnesyl\""]
    L += ["orient hgps_deleted" if sel(d["start"], d["end"]) else "orient", "set sphere_scale, 0.6"]
    return "\n".join(L) + "\n", f"{obj}_progerin.pml"
