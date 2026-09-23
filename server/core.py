"""Structure Bench analyses: gene → protein → structure.

Every public function takes and returns JSON-serialisable dicts and is used unchanged by the HTTP API, the
assistant and the self-test. An *item* (see common.py) is one protein lookup or one uploaded structure;
variant mappings and comparisons are stored inside it, so reloading the page shows the same thing.
"""
from __future__ import annotations

import json
import shutil
import threading
import time
from collections import OrderedDict
from pathlib import Path

import numpy as np

from . import alphafold as afdb
from . import net
from . import structure as st
from . import uniprot as up
from . import variants as vr
from .common import (UserFacingError, VERSION, WORK, dump, item_dir, load_result, log_exc, new_item_dir, now_iso,
                     save_result, touch, versions)

SCHEMA = 1            # bump when result.json changes shape, so an old saved lookup is rebuilt instead of reused
_LOCKS: dict[str, threading.Lock] = {}
_GLOBAL = threading.Lock()


def _lock(item: str) -> threading.Lock:
    with _GLOBAL:
        return _LOCKS.setdefault(item, threading.Lock())


# ---------------------------------------------------------------- parsed-structure cache
_PARSED: OrderedDict = OrderedDict()


def _parsed(path: Path):
    k = (str(path), path.stat().st_mtime)
    if k in _PARSED:
        _PARSED.move_to_end(k)
        return _PARSED[k]
    s, fmt = st.load(path, path.stem)
    _PARSED[k] = (s, fmt)
    while len(_PARSED) > 6:
        _PARSED.popitem(last=False)
    return s, fmt


def _model(item_path: Path, meta: dict):
    s, _ = _parsed(item_path / meta["file"])
    return s[0], s


def rows_for(item_path: Path, sid: str, chain: str) -> list[dict]:
    p = item_path / "structures" / f"{sid}_{chain}.residues.json"
    return json.loads(p.read_text())["rows"] if p.exists() else []


def residues_doc(item_path: Path, sid: str, chain: str) -> dict:
    p = item_path / "structures" / f"{sid}_{chain}.residues.json"
    return json.loads(p.read_text()) if p.exists() else {}


# ---------------------------------------------------------------- small text helpers
def _runs(idx: list[int], gap: int = 1, minlen: int = 1) -> list[tuple[int, int]]:
    """Consecutive integers → [(start, end)], bridging gaps up to `gap`."""
    out = []
    for i in sorted(idx):
        if out and i - out[-1][1] <= gap:
            out[-1][1] = i
        else:
            out.append([i, i])
    return [(a, b) for a, b in out if b - a + 1 >= minlen]


def _rng(runs) -> str:
    return ", ".join(f"{a}–{b}" if a != b else f"{a}" for a, b in runs)


def _pct(x) -> str:
    return f"{100 * x:.0f} %"


def _mean(v):
    v = [x for x in v if x is not None]
    return float(np.mean(v)) if v else None


# ---------------------------------------------------------------- look up
def lookup(gene: str | None = None, species="human", accession: str | None = None, refresh: bool = False,
           reuse: bool = True) -> dict:
    """Gene symbol + organism, or a UniProt accession → an item (or a list of choices when several reviewed
    entries match the gene). Returns {"item": id} | {"choices": [...], "message": ...}."""
    query = {"gene": gene, "species": species, "accession": accession}
    flags = []
    acc = (accession or "").strip().upper()
    if not acc and gene and up.is_accession(gene) and not gene.isalpha():
        acc = gene.strip().upper()
    if acc:
        if not up.is_accession(acc):
            raise UserFacingError(f"'{accession}' is not a UniProt accession (they look like P02545 or Q9Y6K9).")
        acc = acc.split("-")[0]
    else:
        if not gene:
            raise UserFacingError("Type a gene symbol (for example LMNA) or a UniProt accession (P02545).")
        try:
            res = up.search_gene(gene, species, refresh)
        except ValueError as e:
            raise UserFacingError(str(e)) from None
        except net.NetError as e:
            raise UserFacingError(f"UniProt could not be reached, so the gene symbol cannot be turned into a protein. "
                                  f"If you know the UniProt accession, type that instead — earlier look-ups work offline. ({e})") from None
        hits = res["hits"]
        if not hits:
            raise UserFacingError(f"UniProt has no entry with gene name {gene} in {res['species']}. Check the spelling, the "
                                  "organism, or search by accession.")
        if not res["reviewed"]:
            flags.append({"level": "warn", "text": f"No reviewed (Swiss-Prot) entry for {gene} in {res['species']}; using the "
                                                   "unreviewed (TrEMBL) entry — its annotation is automatic."})
        if len(hits) > 1:
            return {"choices": hits, "message": f"{len(hits)} {'reviewed ' if res['reviewed'] else ''}UniProt entries have gene "
                                                 f"name {gene} in {res['species']}. Pick one."}
        acc = hits[0]["accession"]
        query["taxon"] = res["taxon"]
    if reuse and not refresh:
        for d in sorted(WORK.iterdir(), key=lambda x: x.stat().st_mtime, reverse=True):
            try:
                r = json.loads((d / "result.json").read_text())
                if r.get("kind") == "protein" and r.get("accession") == acc and r.get("schema") == SCHEMA:
                    touch(d.name)
                    return {"item": d.name, "reused": True}
            except Exception:  # noqa: BLE001
                continue
    d = new_item_dir()
    try:
        build_protein(d, acc, query, refresh, flags)
    except Exception:
        shutil.rmtree(d, ignore_errors=True)
        raise
    return {"item": d.name}


def build_protein(d: Path, acc: str, query: dict, refresh: bool = False, flags: list | None = None) -> dict:
    t0 = time.time()
    flags = list(flags or [])
    res = {"kind": "protein", "item": d.name, "accession": acc, "created": now_iso(), "app_version": VERSION, "schema": SCHEMA,
           "query": query, "flags": flags, "tiles": [], "structures": {}, "texts": {}, "sources": {},
           "variants": None, "comparisons": [], "timings": {}}
    dump({"kind": "protein", "query": query, "accession": acc}, d / "request.json")
    # --- UniProt
    u = None
    try:
        raw, src = up.fetch_entry(acc, refresh)
        if not raw.get("primaryAccession"):
            raise UserFacingError(f"UniProt has no entry {acc}.")
        if raw.get("entryType", "").lower().startswith("inactive"):
            raise UserFacingError(f"UniProt entry {acc} is inactive (merged or deleted); search by gene instead.")
        u = up.parse_entry(raw)
        res["sources"]["uniprot"] = src
        if u["accession"] != acc:
            flags.append({"level": "info", "text": f"{acc} is now a secondary accession; UniProt's primary accession is <b>{u['accession']}</b>."})
            acc = res["accession"] = u["accession"]
    except net.NetError as e:
        if "does not exist" in str(e):
            raise UserFacingError(f"UniProt has no entry {acc}. Check the accession.") from None
        flags.append({"level": "info", "text": f"UniProt could not be reached — showing the structure only. ({e})"})
    res["timings"]["uniprot"] = round(time.time() - t0, 2)
    res["uniprot"] = u
    # --- AlphaFold
    t1 = time.time()
    pred, iso = None, []
    try:
        preds, src = afdb.predictions(acc, refresh)
        pred, iso = afdb.pick(preds, acc)
        res["sources"]["alphafold"] = src
        if pred is None:
            flags.append({"level": "info", "text": f"AlphaFold DB has no canonical model for {acc} (only isoform models), so the "
                                                   "confidence track and PAE are not available."})
    except net.NetError as e:
        msg = str(e)
        if "does not exist" in msg:
            flags.append({"level": "info", "text": f"AlphaFold DB has no model for {acc} — typical for very long proteins "
                                                   "(> 2,700 residues outside the human proteome), viral proteins and some newer entries."})
        else:
            flags.append({"level": "info", "text": f"AlphaFold DB could not be reached — the predicted structure, PAE and "
                                                   f"AlphaMissense panels are empty for now. ({msg})"})
    res["isoforms"] = [{k: v for k, v in i.items() if k != "sequence"} for i in iso]
    (d / "isoforms.json").write_text(json.dumps(iso))
    af = None
    if pred:
        try:
            af = _add_alphafold(d, res, pred, u, refresh, flags)
        except net.NetError as e:
            flags.append({"level": "info", "text": f"The AlphaFold model file could not be downloaded — {e}"})
        except Exception:  # noqa: BLE001
            flags.append({"level": "warn", "text": "The AlphaFold model could not be read: " + log_exc("alphafold model")})
    res["timings"]["alphafold"] = round(time.time() - t1, 2)
    if u is None and af is None:
        raise UserFacingError(f"Neither UniProt nor AlphaFold DB could be reached for {acc}, and nothing about it is in the "
                              "local cache. Connect to the internet once; after that this protein works offline.")
    # --- PAE, AlphaMissense
    t2 = time.time()
    if pred and af:
        _add_pae(d, res, pred, u, refresh, flags)
        _add_am(d, res, pred, u, refresh, flags)
    res["timings"]["pae_am"] = round(time.time() - t2, 2)
    srcs = set(res["sources"].values())
    if "example" in srcs:
        flags.append({"level": "info", "text": "Loaded from the bundled offline example (UniProt, AlphaFold, PAE and AlphaMissense files "
                                               "shipped with the app) because the internet was not reachable. The data are as of the date they were bundled."})
    elif srcs and srcs <= {"cache"}:
        flags.append({"level": "info", "text": "Loaded from the local cache — nothing was downloaded. Use <b>Refresh</b> to fetch the current versions."})
    # --- names, tiles, texts
    name = f"{u['gene'] or acc} · {acc}" if u else acc
    res["name"] = name
    res["title"] = (u["name"] if u else (pred or {}).get("uniprotDescription", "")) or acc
    res["organism"] = u["organism"] if u else (pred or {}).get("organismScientificName", "")
    L = u["length"] if u else len((pred or {}).get("sequence") or "")
    res["length"] = L
    res["sequence"] = u["sequence"] if u else (pred or {}).get("sequence", "")
    _texts_protein(res)
    _findings(res, flags)
    res["tiles"] = [{"value": acc, "label": "UniProt accession"}, {"value": f"{L:,}", "label": "residues"}]
    if af:
        res["tiles"] += [{"value": f"{af['mean_plddt']:.1f}", "label": "mean pLDDT (AlphaFold)"},
                         {"value": _pct(af["frac_confident"]), "label": "residues confident (pLDDT ≥ 70)"}]
    if u:
        res["tiles"] += [{"value": str(len(u["pdb"])), "label": "PDB entries"},
                         {"value": str(len(u["variants"])), "label": "UniProt natural variants"}]
    res["summary"] = {"accession": acc, "gene": u["gene"] if u else "", "length": L,
                      "plddt": af["mean_plddt"] if af else None, "pdb": len(u["pdb"]) if u else 0}
    res["methods"] = methods(res)
    res["versions"] = versions("biopython", "numpy")
    res["versions"]["Structure Bench"] = VERSION
    res["timings"]["total"] = round(time.time() - t0, 2)
    dump(res, d / "result.json")
    return res


def _add_alphafold(d: Path, res: dict, pred: dict, u: dict | None, refresh: bool, flags: list) -> dict:
    data, src, fname = afdb.model_file(pred, refresh)
    res["sources"]["model"] = src
    (d / "structures" / fname).write_bytes(data)
    s, fmt = _parsed(d / "structures" / fname)
    model = s[0]
    chains = st.polymer_chains(model)
    if not chains:
        raise ValueError("no protein chain in the AlphaFold file")
    ch = chains[0]
    ref = u["sequence"] if u else (pred.get("sequence") or "")
    start = int(pred.get("uniprotStart") or 1)
    if ref and ch["sequence"] == ref[start - 1:start - 1 + len(ch["sequence"])]:
        pos = [start + i for i in range(len(ch["residues"]))]
        ident = 1.0
    else:
        m = st.map_sequence(ch["sequence"], ref) if ref else {"pos": [r.id[1] for r in ch["residues"]], "identity": 1.0}
        pos, ident = m["pos"], m["identity"]
        if ref:
            flags.append({"level": "warn", "text": f"The AlphaFold model's sequence differs from the current UniProt sequence "
                                                   f"(identity {ident:.0%} after alignment) — UniProt may have revised the sequence "
                                                   "since the model was built. Residues were mapped through the alignment."})
    tab = st.residue_table(model, ch["chain"], ch["residues"], pos, d / "structures" / fname, "plddt")
    tab["sid"] = "AF"
    tab["chain"] = ch["chain"]
    dump(tab, d / "structures" / f"AF_{ch['chain']}.residues.json")
    pl = [r["b"] for r in tab["rows"]]
    by_pos = [None] * (len(ref) if ref else len(pl))
    for r in tab["rows"]:
        if r["pos"] and 1 <= r["pos"] <= len(by_pos):
            by_pos[r["pos"] - 1] = r["b"]
    arr = np.array([x for x in pl if x is not None])
    af = {"entry": pred.get("entryId"), "version": pred.get("latestVersion"), "file": f"structures/{fname}",
          "mean_plddt": round(float(arr.mean()), 2), "global_metric": pred.get("globalMetricValue"),
          "frac_very_high": round(float((arr > 90).mean()), 4), "frac_confident": round(float((arr >= 70).mean()), 4),
          "frac_low": round(float(((arr >= 50) & (arr < 70)).mean()), 4), "frac_very_low": round(float((arr < 50).mean()), 4),
          "db_fractions": {k: pred.get(k) for k in ("fractionPlddtVeryHigh", "fractionPlddtConfident", "fractionPlddtLow", "fractionPlddtVeryLow")},
          "plddt": by_pos, "uniprot_start": start, "uniprot_end": pred.get("uniprotEnd"), "source": src,
          "created": pred.get("modelCreatedDate"), "tool": pred.get("toolUsed"),
          "page": f"https://alphafold.ebi.ac.uk/entry/{res['accession']}", "cif_url": pred.get("cifUrl"),
          "pae_image": pred.get("paeImageUrl"), "ss_method": tab["ss_method"], "ss_label": tab["ss_label"]}
    res["alphafold"] = af
    res["structures"]["AF"] = {"sid": "AF", "label": f"AlphaFold {pred.get('entryId')} v{pred.get('latestVersion')}",
                               "kind": "alphafold", "file": f"structures/{fname}", "fmt": fmt, "bkind": "plddt",
                               "chains": [{"chain": ch["chain"], "n": len(ch["residues"]), "first": ch["first"], "last": ch["last"],
                                           "identity": ident, "mapped": True, "molecule": res.get("title", "")}],
                               "main_chain": ch["chain"], "method": "predicted (AlphaFold)", "resolution": None}
    return af


def _add_pae(d: Path, res: dict, pred: dict, u: dict | None, refresh: bool, flags: list):
    try:
        obj, res["sources"]["pae"] = afdb.pae_json(pred, refresh)
        m, mx = afdb.parse_pae(obj)
    except net.NetError as e:
        flags.append({"level": "info", "text": f"The PAE file could not be downloaded, so domain packing cannot be judged — {e}"})
        return
    except Exception:  # noqa: BLE001
        flags.append({"level": "info", "text": "The PAE file could not be read: " + log_exc("pae")})
        return
    blocks = up.main_domains(u)
    basis = "UniProt domains"
    if len(blocks) < 2:
        blocks = afdb.confident_segments(res["alphafold"]["plddt"])
        basis = "confident segments (pLDDT ≥ 70, ≥ 30 residues)"
    summ = afdb.summarize_pae(m, blocks)
    summ["basis"] = basis
    summ["max"] = mx
    summ["n"] = int(m.shape[0])
    disp = afdb.pae_display(m)
    disp["max"] = mx
    dump(disp, d / "pae.json")
    res["pae"] = summ


def _add_am(d: Path, res: dict, pred: dict, u: dict | None, refresh: bool, flags: list):
    if not pred.get("amAnnotationsUrl"):
        res["alphamissense"] = {"available": False, "why": "AlphaFold DB publishes AlphaMissense for canonical human isoforms only."}
        return
    try:
        am, src = afdb.am_table(pred, refresh)
    except net.NetError as e:
        res["alphamissense"] = {"available": False, "why": str(e)}
        flags.append({"level": "info", "text": f"AlphaMissense scores could not be downloaded — {e}"})
        return
    L = res["uniprot"]["length"] if res.get("uniprot") else len(pred.get("sequence") or "")
    means = afdb.am_residue_means(am, L)
    (d / "am.json").write_text(json.dumps({k: list(v) for k, v in am.items()}))
    n_path = sum(1 for s, c in am.values() if c == "likely pathogenic")
    res["sources"]["alphamissense"] = src
    res["alphamissense"] = {"available": True, "n": len(am), "residue_mean": means, "source": src,
                            "frac_pathogenic": round(n_path / len(am), 4) if am else None}


def _am(d: Path) -> dict:
    p = d / "am.json"
    return {k: tuple(v) for k, v in json.loads(p.read_text()).items()} if p.exists() else {}


# ---------------------------------------------------------------- findings box
def _findings(res: dict, flags: list):
    """The 'What your data says' lines for a protein: identity, model confidence, domain packing, PDB coverage."""
    u, af, pae = res.get("uniprot"), res.get("alphafold"), res.get("pae")
    top = []
    if u:
        top.append({"level": "ok" if u["reviewed"] else "warn",
                    "text": f"UniProt <b>{u['accession']}</b> — {u['name']} ({u['gene'] or '—'}, <i>{u['organism']}</i>), {u['length']:,} residues, "
                            + ("reviewed (Swiss-Prot)." if u["reviewed"] else "unreviewed (TrEMBL): annotation is automatic, treat features with care.")})
    if af:
        lvl = "ok" if af["frac_confident"] >= 0.7 else "info"
        top.append({"level": lvl, "text": f"AlphaFold model {af['entry']} (v{af['version']}): mean pLDDT <b>{af['mean_plddt']:.1f}</b>; "
                                          f"{_pct(af['frac_confident'])} of residues confident (≥ 70), {_pct(af['frac_very_low'])} very low (&lt; 50)"
                                          + (" — a large part of this protein is probably disordered on its own, so the model shows it as loose ribbon."
                                             if af["frac_very_low"] >= 0.2 else ".")})
    if pae and pae.get("pairs"):
        bad = [p for p in pae["pairs"] if p["call"] != "confident"]
        if bad:
            what = "domains" if pae.get("basis", "").startswith("UniProt") else "parts"
            top.append({"level": "warn", "text": "The relative placement of " + "; ".join(f"{p['a']} and {p['b']} (mean PAE {p['mean']:.0f} Å)" for p in bad[:3])
                                                 + f" is not confident — do not read anything into how these {what} sit against each other in the model."})
        else:
            top.append({"level": "ok", "text": "All annotated domains have confident relative positions (mean PAE &lt; 10 Å between them)."})
    if u:
        n = len(u["pdb"])
        top.append({"level": "info", "text": f"{n} experimental structure{'s' if n != 1 else ''} listed by UniProt" + (" — see Experimental structures for which parts they cover." if n else " — the AlphaFold model is the only 3D view.")})
    if af and af.get("ss_method") == "phi/psi":
        top.append({"level": "info", "text": "Secondary structure is approximated from backbone phi/psi angles because the DSSP program is not installed; "
                                             "'helix-like' and 'extended' are heuristics, and 'extended' also covers disordered stretches."})
    flags[:0] = top


# ---------------------------------------------------------------- texts for the protein panels
def _texts_protein(res: dict):
    u, af, pae = res.get("uniprot"), res.get("alphafold"), res.get("pae")
    T = res["texts"]
    T["uniprot"] = {"how": "The UniProt entry: curated function, where the protein lives, diseases linked to it, and every "
                           "annotated sequence feature. Positions use the canonical isoform's numbering, as do AlphaFold and "
                           "AlphaMissense — so all panels line up."}
    if u:
        g = u["groups"]
        dz = [x["acronym"] or x["name"] for x in u["diseases"] if x["acronym"] or x["name"]]
        T["uniprot"]["yours"] = (f"{u['name']} ({u['gene']}, {u['organism']}), {u['length']:,} residues, "
                                 f"{'reviewed (Swiss-Prot)' if u['reviewed'] else 'unreviewed (TrEMBL)'}. "
                                 f"{len(g['Domains & regions'])} domains/regions, {len(g['PTMs'])} PTM sites, "
                                 f"{len(u['variants'])} natural variants"
                                 + (f"; linked to {len(dz)} disease{'s' if len(dz) != 1 else ''} ({', '.join(dz[:6])}{'…' if len(dz) > 6 else ''})." if dz else "."))
    T["track"] = {"how": "Each residue is a thin bar coloured by AlphaFold's per-residue confidence (pLDDT, read from the CA "
                         "B-factor column of the model): dark blue > 90, light blue 70–90, yellow 50–70, orange < 50. pLDDT is "
                         "the model's own estimate of how well it placed the residue relative to its neighbours. Below, on the "
                         "same axis, are the UniProt domains and regions and a tick for every natural variant. Hover for details."}
    if af:
        pl = af["plddt"]
        L = len(pl)
        verylow = [i + 1 for i, v in enumerate(pl) if v is not None and v < 50]
        lowruns = _runs(verylow, gap=3, minlen=10)
        parts = [f"{_pct(af['frac_confident'])} of residues are confident (pLDDT ≥ 70) and {_pct(af['frac_very_low'])} "
                 f"very low (< 50); mean pLDDT {af['mean_plddt']:.1f}."]
        if af["frac_low"] >= 0.3:
            parts.append(f"{_pct(af['frac_low'])} sit in the low band (50–70): the backbone there is only roughly placed — "
                         "treat side-chain positions and contacts in those stretches as guesses.")
        rows = json.loads((WORK / res["item"] / "structures" / "AF_A.residues.json").read_text())["rows"] \
            if (WORK / res["item"] / "structures" / "AF_A.residues.json").exists() else []
        ssp = {r["pos"]: r["ss"] for r in rows if r["pos"]}
        word = st.SS_WORD.get(af.get("ss_method"), st.SS_WORD["phi/psi"])
        within = {b["name"]: b["within"] for b in (pae or {}).get("blocks", [])}
        for dom in up.main_domains(u):
            seg = [pl[i - 1] for i in range(dom["start"], min(dom["end"], L) + 1) if pl[i - 1] is not None]
            if not seg:
                continue
            m = float(np.mean(seg))
            hfrac = np.mean([ssp.get(i) == "H" for i in range(dom["start"], dom["end"] + 1)])
            efrac = np.mean([ssp.get(i) == "E" for i in range(dom["start"], dom["end"] + 1)])
            shape = (f"modelled as mostly {word['H']} ({_pct(hfrac)} of residues)" if hfrac >= 0.6 else
                     f"modelled mostly as {word['E']} (β-strand-like) residues ({_pct(efrac)}; {_pct(hfrac)} {word['H']})" if efrac >= 0.4 else
                     f"modelled with {_pct(hfrac)} {word['H']} and {_pct(efrac)} {word['E']} residues")
            conf = "high" if m >= 80 else "moderate" if m >= 70 else "low"
            s = f"The {dom['name']} {dom['start']}–{dom['end']} is {shape}, with {conf} local confidence (mean pLDDT {m:.0f})"
            w = within.get(dom["name"])
            if w is not None and w >= 10 and m >= 70:
                s += f", but see PAE for how its parts pack: the mean PAE within it is {w:.0f} Å"
            parts.append(s + ".")
        if lowruns:
            doms = up.main_domains(u)
            outside = [(a, b) for a, b in lowruns if not any(a >= x["start"] and b <= x["end"] for x in doms)]
            txt = _rng(lowruns[:8]) + ("…" if len(lowruns) > 8 else "")
            parts.append(f"Very low confidence at {txt}"
                         + (" (outside the annotated domains)" if doms and len(outside) == len(lowruns) else "")
                         + " — likely disordered in isolation, not a structure to interpret.")
        if u and u["variants"]:
            vpos = sorted({v["start"] for v in u["variants"]})
            inlow = sum(1 for p in vpos if p <= L and pl[p - 1] is not None and pl[p - 1] < 50)
            parts.append(f"{len(vpos)} positions carry UniProt natural variants; {inlow} of them fall where pLDDT < 50.")
        T["track"]["yours"] = " ".join(parts)
    T["pae"] = {"how": "Predicted aligned error (PAE): for each pair of residues (x, y), AlphaFold's expected error, in Å, in the "
                       "position of residue y if the model were superposed on residue x. Dark (low PAE) means their relative "
                       "position is confident. So a dark square on the diagonal is a confidently folded unit; a dark "
                       "off-diagonal block means two units are confidently placed relative to each other, and a light one "
                       "means the model's arrangement of those units is essentially arbitrary — even if each unit alone has "
                       "high pLDDT. Domain boundaries are marked on both axes; hover for the value."}
    if pae:
        pairs = pae.get("pairs", [])
        s = []
        if pairs:
            for p in pairs[:6]:
                call = {"confident": "confident relative placement", "uncertain": "uncertain relative placement",
                        "not determined": "relative placement not determined — their arrangement in the model is arbitrary"}[p["call"]]
                s.append(f"{p['a']} ↔ {p['b']}: mean PAE {p['mean']:.1f} Å — {call}.")
            conf = [p for p in pairs if p["call"] == "confident"]
            s.insert(0, ("The domains of this protein have confident relative positions." if len(conf) == len(pairs) else
                         "None of the domain pairs have a confident relative position." if not conf else
                         f"{len(conf)} of {len(pairs)} domain pairs have a confident relative position."))
        for b in pae.get("blocks", []):
            if b["within"] >= 10:
                s.append(f"Within {b['name']} ({b['start']}–{b['end']}) the mean PAE is {b['within']:.0f} Å: its parts are confident "
                         "locally but not rigidly placed against each other.")
            else:
                s.append(f"{b['name']} ({b['start']}–{b['end']}) is a compact, confidently folded unit (mean PAE within {b['within']:.0f} Å).")
        if not pae.get("blocks"):
            s.append("No domain of ≥ 20 residues is annotated and no confident segment of ≥ 30 residues was found, so there is nothing to compare.")
        s.append(f"(Blocks used: {pae.get('basis')}; thresholds &lt; 10 Å confident, 10–20 Å uncertain, &gt; 20 Å not determined — "
                 "rules of thumb, not AlphaFold-defined cut-offs.)")
        T["pae"]["yours"] = " ".join(s)
    T["viewer"] = {"how": "The 3D model, drawn with 3Dmol.js. Colour by AlphaFold confidence (same colours as the track), by chain, "
                          "by UniProt domain, or by a per-residue value (relative solvent accessibility, or the mean AlphaMissense "
                          "score over all 19 substitutions). Click a residue for its name, number, pLDDT/B-factor and domain. "
                          "Mapped variants are shown as sticks with a sphere on the CA."}
    if af:
        T["viewer"]["yours"] = (f"Showing {af['entry']} (version {af['version']}), {len(af['plddt']):,} residues, coloured by pLDDT. "
                                + ("Orange and yellow stretches are low-confidence: in the model they often look like long loose "
                                   "ribbons, which is how AlphaFold draws regions it cannot place — not a real extended structure."
                                   if af["frac_very_low"] > 0.05 else "Almost the whole chain is modelled with confidence."))
    T["pdb"] = {"how": "Experimental structures of this protein listed by UniProt: one bar per PDB entry along the sequence, "
                       "showing which residues the deposited construct spans (not every one of them is necessarily resolved). "
                       "Colour shows the method; the label gives the resolution. Load puts the entry in the viewer and makes it "
                       "available for variant mapping and comparison."}
    if u:
        T["pdb"]["yours"] = pdb_yours(u, af, (pae or {}).get("blocks") or up.main_domains(u))


def pdb_yours(u: dict, af: dict | None, blocks: list | None = None) -> str:
    pdb = u["pdb"]
    L = u["length"]
    if not pdb:
        return ("UniProt lists no experimental structure for this protein, so the AlphaFold model is the only 3D view — "
                "read it together with its pLDDT and PAE.")
    cov = np.zeros(L + 2, bool)
    for e in pdb:
        for s in e["segments"]:
            cov[max(1, s["start"]):min(L, s["end"]) + 1] = True
    frac = cov[1:L + 1].mean()
    best = min((e for e in pdb if e["resolution"]), key=lambda e: e["resolution"], default=None)
    meth = {}
    for e in pdb:
        meth[e["method"]] = meth.get(e["method"], 0) + 1
    parts = [f"{len(pdb)} entries ({', '.join(f'{v} {k}' for k, v in meth.items())}) together span {_pct(frac)} of the sequence."]
    if best:
        parts.append(f"Best resolution {best['resolution']:.2f} Å ({best['id']}, residues {best['span'][0]}–{best['span'][1]}).")
        big = min((e for e in pdb if e["resolution"] and e["span"] and e["covered"] >= 50), key=lambda e: e["resolution"], default=None)
        if big and big is not best:
            parts.append(f"Best among entries spanning ≥ 50 residues: {big['resolution']:.2f} Å ({big['id']}, {big['span'][0]}–{big['span'][1]}).")
    for dom in (blocks if blocks is not None else up.main_domains(u)):
        n = [e["id"] for e in pdb if any(s["start"] <= dom["end"] and s["end"] >= dom["start"] and
                                           (min(s["end"], dom["end"]) - max(s["start"], dom["start"]) + 1) >= 0.5 * (dom["end"] - dom["start"] + 1)
                                           for s in e["segments"])]
        whole = [e["id"] for e in pdb if any(s["start"] <= dom["start"] and s["end"] >= dom["end"] for s in e["segments"])]
        parts.append(f"{dom['name']} {dom['start']}–{dom['end']}: " +
                     (f"{len(whole)} entr{'y' if len(whole) == 1 else 'ies'} cover it whole ({', '.join(whole[:5])}{'…' if len(whole) > 5 else ''})" if whole
                      else f"no entry covers it whole; {len(n)} cover at least half of it" + (f" ({', '.join(n[:5])})" if n else "")) + ".")
    gaps = _runs([i for i in range(1, L + 1) if not cov[i]], gap=0, minlen=20)
    if gaps:
        extra = ""
        if af:
            low = [g for g in gaps if (_mean(af["plddt"][g[0] - 1:g[1]]) or 0) < 50]
            if low:
                extra = f" ({len(low)} of these {len(gaps)} are also very low confidence in AlphaFold — probably disordered, which is often why they were left out of crystallised constructs)"
        parts.append(f"No experimental structure spans {_rng(gaps[:6])}{'…' if len(gaps) > 6 else ''}{extra}.")
    return " ".join(parts)


# ---------------------------------------------------------------- loading more structures
def list_pdb_entries(item: str) -> dict:
    res = load_result(item)
    u = res.get("uniprot")
    if not u:
        return {"entries": [], "note": "No UniProt data for this item, so no list of PDB entries."}
    loaded = {s["pdb_id"] for s in res["structures"].values() if s.get("pdb_id")}
    return {"entries": [{**e, "loaded": e["id"] in loaded} for e in u["pdb"]], "yours": res["texts"].get("pdb", {}).get("yours")}


def load_pdb(item: str, pdb_id: str, refresh: bool = False) -> dict:
    """Download https://files.rcsb.org/download/{id}.cif (cached), add it to the item, map its chains to UniProt."""
    pdb_id = (pdb_id or "").strip().upper()
    if not (len(pdb_id) == 4 and pdb_id[0].isdigit() and pdb_id.isalnum()):
        raise UserFacingError(f"'{pdb_id}' is not a PDB ID (four characters, starting with a digit, like 1IFR).")
    with _lock(item):
        d = item_dir(item)
        res = load_result(item)
        sid = f"PDB-{pdb_id}"
        if sid in res["structures"] and not refresh:
            return {"sid": sid, "structure": res["structures"][sid], "already": True}
        try:
            data, src = net.fetch(f"https://files.rcsb.org/download/{pdb_id}.cif", f"pdb/{pdb_id}.cif",
                                  f"PDB entry {pdb_id}", refresh)
        except net.NetError as e:
            raise UserFacingError(f"PDB entry {pdb_id} could not be downloaded from RCSB and is not in the local cache — {e}") from None
        fname = f"{pdb_id}.cif"
        (d / "structures" / fname).write_bytes(data)
        entry = next((e for e in (res.get("uniprot") or {}).get("pdb", []) if e["id"] == pdb_id), None)
        try:
            meta = _register(d, res, sid, fname, "pdb", f"PDB {pdb_id}",
                             hint_chains=[c for s in (entry or {}).get("segments", []) for c in s["chains"]])
        except Exception:
            (d / "structures" / fname).unlink(missing_ok=True)
            raise
        meta["pdb_id"] = pdb_id
        if entry:
            meta["method"], meta["resolution"] = entry["method"], entry["resolution"]
        meta["source"] = src
        save_result(item, res)
        return {"sid": sid, "structure": meta}


def load_isoform(item: str, entry: str) -> dict:
    """Add an AlphaFold isoform model (listed under the canonical one) as another structure."""
    with _lock(item):
        d = item_dir(item)
        res = load_result(item)
        iso = next((i for i in res.get("isoforms", []) if i["entry"] == entry), None)
        if not iso:
            raise UserFacingError(f"{entry} is not one of this protein's AlphaFold isoform models.")
        try:
            data, src = net.fetch(iso["cif"], f"alphafold/{iso['cif'].split('/')[-1]}", f"the AlphaFold model {entry}")
        except net.NetError as e:
            raise UserFacingError(str(e)) from None
        fname = iso["cif"].split("/")[-1]
        (d / "structures" / fname).write_bytes(data)
        sid = entry
        meta = _register(d, res, sid, fname, "alphafold-isoform", f"AlphaFold {entry} (isoform)")
        meta["bkind"] = "plddt"
        save_result(item, res)
        return {"sid": sid, "structure": meta}


def _register(d: Path, res: dict, sid: str, fname: str, kind: str, label: str, hint_chains=None) -> dict:
    """Parse a structure file in the item, map every protein chain (to UniProt when there is one), and compute the
    residue table (pLDDT/B-factor, SASA alone and in the complex, secondary structure) for its main chain."""
    path = d / "structures" / fname
    try:
        s, fmt = _parsed(path)
    except KeyError as e:
        log_exc("parse " + fname)
        raise UserFacingError(f"{fname} could not be read as a structure: the {e} field is missing — the file has no "
                              "usable atom records.") from None
    except Exception as e:  # noqa: BLE001
        log_exc("parse " + fname)
        raise UserFacingError(f"{fname} could not be read as a structure: {e}") from None
    model = s[0]
    chains = st.polymer_chains(model)
    if not chains:
        raise UserFacingError(f"{fname} contains no protein chain (only nucleic acid, ligands or CA-less residues?).")
    names = st.entity_names(path) if fmt == "cif" else {}
    u = res.get("uniprot")
    ref = (u or {}).get("sequence") or (res.get("sequence") if res.get("kind") == "protein" else "")
    info = []
    for ch in chains:
        c = {"chain": ch["chain"], "n": len(ch["residues"]), "first": ch["first"], "last": ch["last"],
             "molecule": names.get(ch["chain"], ""), "fragments": ch["fragments"]}
        if ref:
            m = st.map_sequence(ch["sequence"], ref)
            need = min(15, int(0.8 * len(ch["residues"])))
            c.update({"identity": m["identity"], "aligned": m["aligned"], "mapped": m["identity"] >= 0.85 and m["aligned"] >= need,
                      "span": [min(p for p in m["pos"] if p), max(p for p in m["pos"] if p)] if m["aligned"] else None})
            c["_pos"] = m["pos"]
        else:
            c.update({"identity": None, "mapped": False})
        info.append(c)
    if ref:
        cand = [c for c in info if c["mapped"]]
        if not cand:
            best = max(info, key=lambda c: c.get("aligned") or 0)
            what = ", ".join(f"{c['chain']}" + (f" ({c['molecule']})" if c["molecule"] else "") for c in info[:6])
            near = (f" The closest is chain {best['chain']}: {best['identity']:.0%} identity over {best['aligned']} residues."
                    if best.get("aligned") else "")
            raise UserFacingError(f"No chain in {fname} matches this protein's sequence (its protein chains: {what}).{near} "
                                  "It may be an entry of a different protein or a homologue from another species; UniProt "
                                  "lists this protein's own entries under Experimental structures.")
        hints = set(hint_chains or [])
        best = max(c["aligned"] for c in cand)
        near_full = [c for c in cand if c["aligned"] >= 0.9 * best]
        copies = {c["chain"] for c in cand}
        if len(near_full) > 1:
            # several copies of the protein: prefer the one touching the most of something else (DNA, a partner,
            # a ligand) — in a p53–DNA crystal that is the copy actually bound to the DNA
            for c in near_full:
                c["_contacts"] = st.foreign_contacts(model, c["chain"], copies)
        main = max(near_full, key=lambda c: (c["chain"] in hints, c.get("_contacts", 0), c["aligned"]))
        pos = main["_pos"]
    else:
        main = max(info, key=lambda c: c["n"])
        pos = [None] * main["n"]
    ch = next(x for x in chains if x["chain"] == main["chain"])
    bkind = "plddt" if kind.startswith("alphafold") or _looks_like_af(model) else "bfactor"
    tab = st.residue_table(model, ch["chain"], ch["residues"], pos, path, bkind)
    tab["sid"], tab["chain"] = sid, ch["chain"]
    if not ref:          # upload item: the structure's own numbering is the reference
        for r in tab["rows"]:
            r["pos"] = r["num"] if not r["icode"] else None
    dump(tab, d / "structures" / f"{sid}_{ch['chain']}.residues.json")
    for c in info:
        c.pop("_pos", None)
        c.pop("_contacts", None)
    view = fname
    if fmt == "bcif":                     # 3Dmol.js reads mmCIF; keep the original and add a text copy for the viewer
        view = Path(fname).stem + ".view.cif"
        st.write(s, d / "structures" / view, "cif")
    nmodels = len(s)
    have = sorted({r["pos"] for r in tab["rows"] if r["pos"]})
    gaps = _runs([p for p in range(have[0], have[-1] + 1) if p not in set(have)], gap=0) if have else []
    meta = {"sid": sid, "label": label, "kind": kind, "file": f"structures/{fname}", "view_file": f"structures/{view}",
            "fmt": fmt, "chains": info, "main_chain": ch["chain"], "bkind": bkind, "names": names,
            "partners": tab["partners"], "ss_method": tab["ss_method"], "models": nmodels,
            "method": "", "resolution": None, "span": [have[0], have[-1]] if have else None,
            "unresolved": [list(g) for g in gaps]}
    if nmodels > 1:
        meta["note"] = f"{nmodels} models in the file (NMR ensemble?) — only the first is used."
    res["structures"][sid] = meta
    return meta


def _looks_like_af(model) -> bool:
    """An uploaded AlphaFold/ColabFold model: B-factors all in 0–100 and the same value on every atom of a residue."""
    try:
        bs = []
        for r in list(model.get_residues())[:80]:
            v = {round(a.get_bfactor(), 2) for a in r}
            if len(v) != 1:
                return False
            bs.append(v.pop())
        return bool(bs) and 0 <= min(bs) and max(bs) <= 100 and max(bs) > 1
    except Exception:  # noqa: BLE001
        return False


# ---------------------------------------------------------------- uploads
def upload(filename: str, data: bytes, item: str | None = None) -> dict:
    """A structure file from the user. With `item`, it is added to that protein; otherwise a new offline item."""
    name = Path(filename or "structure").name
    if not data:
        raise UserFacingError(f"{name} is empty.")
    safe = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in name)[:80] or "structure"
    if item:
        with _lock(item):
            d = item_dir(item)
            res = load_result(item)
            (d / "structures" / safe).write_bytes(data)
            if st.sniff(d / "structures" / safe) == "unknown":
                (d / "structures" / safe).unlink()
                raise UserFacingError(f"{name} is not a PDB, mmCIF or BinaryCIF file (.pdb, .cif, .bcif, optionally .gz).")
            sid = "UP-" + Path(safe).stem.split(".")[0][:24]
            try:
                meta = _register(d, res, sid, safe, "upload", f"Upload {name}")
            except Exception:
                (d / "structures" / safe).unlink(missing_ok=True)
                raise
            save_result(item, res)
            return {"item": item, "sid": sid, "structure": meta}
    d = new_item_dir()
    try:
        (d / "structures" / safe).write_bytes(data)
        if st.sniff(d / "structures" / safe) == "unknown":
            raise UserFacingError(f"{name} is not a PDB, mmCIF or BinaryCIF file (.pdb, .cif, .bcif, optionally .gz).")
        res = {"kind": "upload", "item": d.name, "name": name, "title": name, "created": now_iso(), "app_version": VERSION,
               "query": {"upload": name}, "flags": [], "tiles": [], "structures": {}, "texts": {}, "sources": {},
               "variants": None, "comparisons": [], "uniprot": None, "alphafold": None, "pae": None, "timings": {}}
        dump({"kind": "upload", "file": safe}, d / "request.json")
        t0 = time.time()
        meta = _register(d, res, "UP-1", safe, "upload", name)
        res["timings"]["register"] = round(time.time() - t0, 2)
        rows = rows_for(d, "UP-1", meta["main_chain"])
        doc = residues_doc(d, "UP-1", meta["main_chain"])
        res["length"] = len(rows)
        res["sequence"] = "".join(r["aa"] for r in rows)
        bvals = [r["b"] for r in rows if r["b"] is not None]
        res["flags"] = [
            {"level": "info", "text": "Uploaded structure, analysed offline: no UniProt features, AlphaFold confidence, PAE or "
                                      "AlphaMissense (they need a UniProt accession). Residues use the file's own numbering."},
            {"level": "ok", "text": f"{len(meta['chains'])} protein chain{'s' if len(meta['chains']) != 1 else ''} "
                                    f"({', '.join(c['chain'] + ' ' + str(c['n']) + ' aa' for c in meta['chains'][:8])}); main chain "
                                    f"{meta['main_chain']} numbered {rows[0]['num']}–{rows[-1]['num']}."},
        ]
        if meta["bkind"] == "plddt":
            res["flags"].append({"level": "info", "text": f"The B-factor column looks like a predicted-model confidence (pLDDT, mean "
                                                          f"{np.mean(bvals):.1f}), so it is coloured and read as pLDDT."})
        if meta.get("note"):
            res["flags"].append({"level": "info", "text": meta["note"]})
        if doc.get("ss_method") == "phi/psi":
            res["flags"].append({"level": "info", "text": "Secondary structure is an approximation from backbone phi/psi (the DSSP "
                                                          "program is not installed)."})
        res["tiles"] = [{"value": str(len(meta["chains"])), "label": "protein chains"}, {"value": str(len(rows)), "label": f"residues in chain {meta['main_chain']}"},
                        {"value": f"{np.mean(bvals):.1f}" if bvals else "—", "label": "mean pLDDT" if meta["bkind"] == "plddt" else "mean CA B-factor"}]
        res["summary"] = {"chains": len(meta["chains"]), "length": len(rows)}
        res["texts"] = {"viewer": {"how": "Your structure, drawn with 3Dmol.js. Colour by chain, by B-factor (pLDDT if the file is a "
                                          "predicted model) or by relative solvent accessibility; click a residue for its details.",
                                   "yours": f"{name}: {len(meta['chains'])} protein chain(s), {sum(c['n'] for c in meta['chains'])} residues in total."}}
        res["methods"] = methods(res)
        res["versions"] = versions("biopython", "numpy")
        res["versions"]["Structure Bench"] = VERSION
        dump(res, d / "result.json")
        return {"item": d.name, "sid": "UP-1", "structure": meta}
    except Exception:
        shutil.rmtree(d, ignore_errors=True)
        raise


# ---------------------------------------------------------------- variants
def _structure_meta(res: dict, sid: str | None) -> dict:
    sid = sid or ("AF" if "AF" in res["structures"] else next(iter(res["structures"]), None))
    if not sid or sid not in res["structures"]:
        have = ", ".join(res["structures"]) or "none"
        raise UserFacingError(f"There is no structure '{sid}' in this item (loaded: {have}). Load it first.")
    return res["structures"][sid]


def _known(u: dict | None, pos: int, alt: str) -> list[dict]:
    if not u:
        return []
    out = []
    for v in u["variants"]:
        if v["start"] == pos and v["end"] == pos and alt in (v.get("alt") or []):
            out.append({"description": v["description"], "id": v.get("id", "")})
    return out


def _nearby_known(u, pos):
    if not u:
        return []
    return sorted({f"{v.get('ref', '')}{v['start']}{'/'.join(v.get('alt') or [])}" for v in u["variants"] if v["start"] == pos})


def map_variants(item: str, text: str, structure: str | None = None, chain: str | None = None, save: bool = True) -> dict:
    """Parse, validate and place variants on a structure of the item. See the module docstring of variants.py."""
    with _lock(item):
        d = item_dir(item)
        res = load_result(item)
        meta = _structure_meta(res, structure)
        sid = meta["sid"]
        ch = chain or meta["main_chain"]
        if ch not in [c["chain"] for c in meta["chains"]]:
            raise UserFacingError(f"Chain {ch} is not a protein chain of {meta['label']} ({', '.join(c['chain'] for c in meta['chains'])}).")
        rows = rows_for(d, sid, ch)
        if not rows:                          # a secondary chain: compute its table on demand
            rows = _chain_table(d, res, meta, ch)
        doc = residues_doc(d, sid, ch)
        by_pos = {r["pos"]: r for r in rows if r["pos"]}
        u = res.get("uniprot")
        protein = res["kind"] == "protein"
        seq = res.get("sequence") or ""
        af_pl = (res.get("alphafold") or {}).get("plddt") or []
        am = _am(d) if protein else {}
        am_ok = bool(am)
        isoforms = json.loads((d / "isoforms.json").read_text()) if (d / "isoforms.json").exists() else []
        chain_offsets = []
        for f in (u or {}).get("features", []):
            if f["type"] in ("Chain", "Propeptide") and f["start"] > 1:
                chain_offsets.append((f"numbering from the start of the mature {f['type'].lower()} at residue {f['start']}", f["start"] - 1))
        model, _ = _model(d, meta)
        parsed = vr.parse_many(text)
        if not parsed:
            raise UserFacingError("No variants given. Type one per line, e.g. R482W or p.Arg482Trp.")
        if len(parsed) > 200:
            raise UserFacingError(f"{len(parsed)} variants — map at most 200 at a time.")
        out = []
        for v in parsed:
            row = {"input": v["input"]}
            if "error" in v:
                row.update({"error": v["error"], "flags": [{"level": "warn", "text": v["error"]}]})
                out.append(row)
                continue
            row.update({k: v[k] for k in ("ref", "pos", "alt", "kind", "label")})
            if protein and seq:
                chk = vr.check_reference(v, seq, (u or {}).get("signal_length", 0), isoforms, chain_offsets,
                                         f"{res['accession']} (canonical isoform)")
            else:
                r0 = by_pos.get(v["pos"])
                got = r0["aa"] if r0 else None
                chk = {"ok": got == v["ref"], "message": "" if got == v["ref"] else
                       (f"Residue {v['pos']} of chain {ch} is {vr.NAME.get(got, got)}, not {vr.NAME.get(v['ref'])}." if got
                        else f"Chain {ch} has no residue numbered {v['pos']}."),
                       "suggest": [] if got == v["ref"] else ["Upload items use the file's own residue numbering; check which numbering the variant uses."]}
            row["ref_ok"] = chk["ok"]
            if not chk["ok"]:
                row["error"] = chk["message"] + " " + " ".join(chk.get("suggest", []))
                row["flags"] = [{"level": "warn", "text": row["error"]}]
                out.append(row)
                continue
            p = v["pos"]
            row["domain"] = up.domain_at(u, p)
            row["coiled_coil"] = any(f["start"] <= p <= f["end"] and (f["type"] == "Coiled coil" or
                                     (f["type"] == "Region" and "coil" in f["description"].lower()))
                                     for f in (u or {}).get("features", []))
            row["uniprot_known"] = _known(u, p, v["alt"])
            row["uniprot_other_at_pos"] = _nearby_known(u, p)
            row["plddt"] = af_pl[p - 1] if protein and 0 < p <= len(af_pl) else (None if meta["bkind"] != "plddt" else None)
            if am_ok and v["kind"] == "missense":
                s_ = am.get(v["label"])
                if s_:
                    row["am_score"], row["am_class"] = s_
            r = by_pos.get(p)
            row["structure"] = sid
            row["chain"] = ch
            if r is None:
                row["unresolved"] = True
            else:
                row.update({"key": r["key"], "num": r["num"], "rsa": r["rsa"], "rsa_alone": r["rsa_alone"],
                            "rsa_class": st.rsa_class(r["rsa"]), "rsa_class_alone": st.rsa_class(r["rsa_alone"]), "ss": r["ss"],
                            "ss_word": st.SS_WORD.get(doc.get("ss_method", "phi/psi"), st.SS_WORD["phi/psi"])[r["ss"]]})
                if meta["bkind"] == "plddt":
                    if not protein or sid != "AF":
                        row["model_plddt"] = r["b"]
                    if not protein:
                        row["plddt"] = r["b"]
                else:
                    row["bfactor"] = r["b"]
                try:
                    nb = st.neighbours(model, r["key"], 5.0, meta.get("names"))
                    row["neighbours"] = {"n": nb["n_total"], "nonlocal": nb["n_nonlocal"],
                                         "residues": [f"{x['resname']}{x['num']}" for x in nb["same_chain"] if x["seq_sep"] > 2][:14],
                                         "other_chain": [f"{x['chain']}:{x['resname']}{x['num']}" + (f" ({x['molecule']})" if x.get("molecule") else "") for x in nb["other_chain"]][:10],
                                         "ligand": [f"{x['resname']} {x['chain']}{x['num']}" for x in nb["ligand"]][:10],
                                         "water": nb["water"]}
                    if nb["ligand"]:
                        row["ligand_contacts"] = ", ".join(sorted({x["resname"] for x in nb["ligand"]}))
                    lost = (r["sasa_alone"] or 0) - (r["sasa"] or 0)
                    if nb["other_chain"] and lost > 5:
                        partners = sorted({x["chain"] for x in nb["other_chain"]})
                        nm = meta.get("names", {})
                        same = [c for c in partners if nm.get(c) and nm.get(c) == nm.get(ch)]
                        row["interface"] = ", ".join(f"chain {c}" + (f" ({nm[c]})" if nm.get(c) else "") for c in partners)
                        row["interface_lost"] = round(lost, 1)
                        if same:
                            row["interface_note"] = "another copy of the same protein — a real dimer contact or a crystal contact"
                except Exception:  # noqa: BLE001
                    row["neighbours"] = None
                    log_exc("neighbours")
            row["change"] = vr.change_words(v["ref"], v["alt"])
            # the confidence the reading uses is the chosen structure's own: a crystal that resolved a residue is
            # informative even where the AlphaFold model is not (pLDDT is still shown, for reference)
            row["reading_plddt"] = (row.get("model_plddt", row.get("plddt")) if meta["bkind"] == "plddt" else None)
            row["flags"] = vr.interpret({**row, "plddt": row["reading_plddt"]})
            out.append(row)
        result = {"structure": sid, "structure_label": meta["label"], "chain": ch, "rows": out, "text": text,
                  "am_available": am_ok, "ss_method": doc.get("ss_method"), "bkind": meta["bkind"], "time": now_iso()}
        result["summary"], result["yours"] = _variant_summary(result, res)
        result["how"] = ("Each variant is checked against the reference sequence, placed on the chosen structure (AlphaFold "
                         "positions are UniProt positions; a PDB chain is mapped to UniProt through a global BLOSUM62 alignment), "
                         "and described: AlphaFold confidence (pLDDT) or crystallographic B-factor at that residue; relative solvent "
                         "accessibility (Shrake–Rupley SASA ÷ the Tien 2013 maximum; buried < 20 %, exposed > 50 %); residues "
                         "within 5 Å; secondary structure; the UniProt domain and whether UniProt lists this exact change; and the "
                         "AlphaMissense score. The flags are heuristics — rules of thumb for where to look, not predictions of "
                         "pathogenicity.")
        if save:
            res["variants"] = result
            dump(result, d / "variants.json")
            save_result(item, res)
        return result


def _chain_table(d: Path, res: dict, meta: dict, chain: str) -> list[dict]:
    model, _ = _model(d, meta)
    chs = {c["chain"]: c for c in st.polymer_chains(model)}
    c = chs[chain]
    ref = res.get("sequence") if res["kind"] == "protein" else None
    pos = st.map_sequence(c["sequence"], ref)["pos"] if ref else [r.id[1] for r in c["residues"]]
    tab = st.residue_table(model, chain, c["residues"], pos, d / meta["file"], meta["bkind"])
    tab["sid"], tab["chain"] = meta["sid"], chain
    dump(tab, d / "structures" / f"{meta['sid']}_{chain}.residues.json")
    return tab["rows"]


def _variant_summary(vres: dict, res: dict) -> tuple[dict, str]:
    rows = [r for r in vres["rows"] if not r.get("error")]
    bad = [r for r in vres["rows"] if r.get("error")]
    syn = [r for r in rows if r["kind"] == "synonymous"]
    mis = [r for r in rows if r["kind"] == "missense"]
    unres = [r for r in mis if r.get("unresolved")]
    placed = [r for r in mis if not r.get("unresolved")]
    low = [r for r in placed if r.get("reading_plddt") is not None and r["reading_plddt"] < 50]
    conf = [r for r in placed if r.get("reading_plddt") is None or r["reading_plddt"] >= 70]
    fold = lambda r: r.get("rsa_alone") if r.get("rsa_alone") is not None else r.get("rsa")  # noqa: E731
    buried = [r for r in conf if fold(r) is not None and fold(r) < 0.2]
    exposed = [r for r in conf if fold(r) is not None and fold(r) > 0.5]
    iface = [r for r in placed if r.get("interface")]
    known = [r for r in rows if r.get("uniprot_known")]
    amp = [r for r in mis if r.get("am_class") == "likely pathogenic"]
    lab = lambda xs: ", ".join(r["label"] for r in xs[:12]) + (f" … and {len(xs) - 12} more" if len(xs) > 12 else "")  # noqa: E731
    parts = [f"{len(vres['rows'])} variant{'s' if len(vres['rows']) != 1 else ''} on {vres['structure_label']} chain {vres['chain']}."]
    if bad:
        parts.append(f"{len(bad)} could not be placed ({', '.join(r.get('label') or r['input'] for r in bad)}) — see the error in its row.")
    where = "a confident region" if vres.get("bkind") == "plddt" else "the resolved structure"
    if buried:
        parts.append(f"Buried in {where}: {lab(buried)} — the ones most likely to destabilise the fold.")
    if exposed:
        parts.append(f"Exposed in {where}: {lab(exposed)} — possible binding-surface positions.")
    partly = [r for r in conf if fold(r) is not None and 0.2 <= fold(r) <= 0.5]
    if partly:
        parts.append(f"Partly buried in {where}: {lab(partly)} — the effect depends on how the side chain packs.")
    if low:
        parts.append(f"In low-confidence regions (pLDDT < 50), where structure says little: {lab(low)}.")
    if unres:
        parts.append(f"Not modelled in this structure: {lab(unres)}.")
    if iface:
        parts.append(f"At an interface in this entry: {lab(iface)}.")
    if syn:
        parts.append(f"Synonymous (no amino-acid change, not assessable structurally): {lab(syn)}.")
    if known:
        parts.append(f"Listed by UniProt as natural variants: {lab(known)}.")
    if vres["am_available"] and mis:
        parts.append(f"AlphaMissense calls {len(amp)} of {len(mis)} missense changes likely pathogenic" + (f" ({lab(amp)})" if amp else "") + " — a prediction, not a clinical classification.")
    elif mis and res["kind"] == "protein":
        parts.append("AlphaMissense is not available for this protein (canonical human isoforms only).")
    meta = res["structures"].get(vres["structure"], {})
    copies = [c["chain"] for c in meta.get("chains", []) if c.get("mapped")]
    if len(copies) > 1:
        parts.append(f"This entry holds {len(copies)} copies of the protein (chains {', '.join(copies)}); this reading is for chain "
                     f"{vres['chain']}. Copies can differ in their contacts (one may sit on DNA or a partner, another not) — "
                     "choose another chain to compare.")
    s = {"n": len(vres["rows"]), "errors": len(bad), "buried_confident": len(buried), "exposed_confident": len(exposed),
         "low_confidence": len(low), "unresolved": len(unres), "synonymous": len(syn), "interface": len(iface),
         "am_pathogenic": len(amp)}
    return s, " ".join(parts)


def alphamissense(item: str, variants: str) -> dict:
    res = load_result(item)
    d = item_dir(item)
    am = _am(d)
    if not am:
        why = (res.get("alphamissense") or {}).get("why") or "AlphaMissense needs a canonical human UniProt protein."
        return {"available": False, "why": why, "rows": []}
    rows = []
    seq = res.get("sequence") or ""
    for v in vr.parse_many(variants):
        if "error" in v:
            rows.append({"input": v["input"], "error": v["error"]})
            continue
        if v["kind"] != "missense":
            rows.append({"input": v["input"], "label": v["label"], "note": "AlphaMissense scores missense changes only."})
            continue
        ok = 0 < v["pos"] <= len(seq) and seq[v["pos"] - 1] == v["ref"]
        s = am.get(v["label"])
        rows.append({"input": v["input"], "label": v["label"], "ref_ok": ok, "score": s[0] if s else None,
                     "class": s[1] if s else None, **({} if ok else {"error": f"position {v['pos']} is not {v['ref']} in the canonical sequence"})})
    return {"available": True, "rows": rows,
            "note": "AlphaMissense (Cheng et al. 2023, Science): a prediction of pathogenicity from sequence and structure context. "
                    "Scores > 0.564 are called likely pathogenic, < 0.34 likely benign, in between ambiguous. Canonical isoform only."}


# ---------------------------------------------------------------- residue environment
def residue_environment(item: str, residue: str, radius: float = 5.0, structure: str | None = None) -> dict:
    """Everything about one residue: its row (pLDDT/B, SASA, SS), domain, and neighbours within `radius` Å.
    `residue` is a UniProt position ('482'), a variant ('R482W'), or a structure key ('A:482')."""
    d = item_dir(item)
    res = load_result(item)
    meta = _structure_meta(res, structure)
    ch = meta["main_chain"]
    r = str(residue).strip()
    rows = rows_for(d, meta["sid"], ch)
    if ":" in r:
        key = r
        row = next((x for x in rows if x["key"] == key), None)
    else:
        v = vr.parse(r)
        p = v["pos"] if "pos" in v else int("".join(c for c in r if c.isdigit()) or 0)
        row = next((x for x in rows if x["pos"] == p), None)
        if row is None:
            raise UserFacingError(f"Residue {p} is not modelled in {meta['label']} chain {ch}.")
        key = row["key"]
    radius = max(2.0, min(float(radius or 5.0), 12.0))
    model, _ = _model(d, meta)
    try:
        nb = st.neighbours(model, key, radius, meta.get("names"))
    except KeyError:
        raise UserFacingError(f"{key} is not a residue of {meta['label']}.") from None
    pos = row["pos"] if row else None
    out = {"structure": meta["sid"], "label": meta["label"], "residue": row, "radius": radius, "neighbours": nb,
           "domain": up.domain_at(res.get("uniprot"), pos) if pos else [],
           "plddt_af": ((res.get("alphafold") or {}).get("plddt") or [None] * (pos or 0))[pos - 1] if pos and res.get("alphafold") else None,
           "rsa_class": st.rsa_class(row["rsa"]) if row else None}
    return out


# ---------------------------------------------------------------- comparison
def compare_structures(item: str, a: str, b: str, chain_a: str | None = None, chain_b: str | None = None) -> dict:
    with _lock(item):
        d = item_dir(item)
        res = load_result(item)
        ma, mb = _structure_meta(res, a), _structure_meta(res, b)
        ca, cb = chain_a or ma["main_chain"], chain_b or mb["main_chain"]
        if ma["sid"] == mb["sid"] and ca == cb and a != b:
            pass
        ra = rows_for(d, ma["sid"], ca) or _chain_table(d, res, ma, ca)
        rb = rows_for(d, mb["sid"], cb) or _chain_table(d, res, mb, cb)
        model_a, sa = _model(d, ma)
        model_b, sb = _model(d, mb)
        t0 = time.time()
        try:
            ce = st.ce_rmsd(model_a, ca, model_b, cb)
        except Exception as e:  # noqa: BLE001
            log_exc("cealign")
            ce = {"rmsd": None, "error": f"CE alignment failed: {e}"}
        # residue correspondence: by UniProt position when both are mapped; otherwise by aligning b's chain to a's
        if res["kind"] == "protein":
            pa = {r["pos"]: r for r in ra if r["pos"]}
            pairs = [(p, pa[p], r) for r in rb if r["pos"] and r["pos"] in pa for p in [r["pos"]]]
            numbering = "UniProt"
        else:
            m = st.map_sequence("".join(r["aa"] for r in rb), "".join(r["aa"] for r in ra))
            pairs = [(ra[p - 1]["num"], ra[p - 1], r) for r, p in zip(rb, m["pos"]) if p]
            numbering = f"{ma['label']} chain {ca}"
        fixed, moving, keep = [], [], []
        chain_a_obj, chain_b_obj = model_a[ca], model_b[cb]
        for p, x, y in pairs:
            try:
                fa = st.find_residue(model_a, x["key"])["CA"]
                fb = st.find_residue(model_b, y["key"])["CA"]
            except KeyError:
                continue
            fixed.append(fa)
            moving.append(fb)
            keep.append((p, x, y))
        if len(fixed) < 3:
            raise UserFacingError(f"{ma['label']} and {mb['label']} share fewer than 3 matched residues — they do not overlap "
                                  "in sequence, so they cannot be superposed residue by residue.")
        moved, _ = st.load(d / mb["file"], "moved")          # fresh parse — see structure.superpose_pairs
        mv = [st.find_residue(moved[0], y["key"])["CA"] for _, _, y in keep]
        rms, (rot, tran) = st.superpose_pairs(fixed, mv, moved)
        dev = []
        af_pl = (res.get("alphafold") or {}).get("plddt") or []
        for (p, x, y), fa, mb_ca in zip(keep, fixed, mv):
            dist = float(np.linalg.norm(fa.coord - mb_ca.coord))
            pl = af_pl[p - 1] if res["kind"] == "protein" and 0 < p <= len(af_pl) else (x["b"] if ma["bkind"] == "plddt" else None)
            dev.append({"pos": p, "a": x["key"], "b": y["key"], "dev": round(dist, 2), "plddt": pl})
        fname = f"superposed_{mb['sid']}_{cb}_on_{ma['sid']}_{ca}.cif"
        st.write(moved, d / "structures" / fname, "cif")
        big = [x["pos"] for x in dev if x["dev"] > 3.0]
        segs = _runs(big, gap=3, minlen=1)
        seg_rows = []
        ends = {dev[0]["pos"], dev[-1]["pos"]} if dev else set()
        for s0, s1 in segs:
            xs = [x for x in dev if s0 <= x["pos"] <= s1]
            pls = [x["plddt"] for x in xs if x["plddt"] is not None]
            seg_rows.append({"start": s0, "end": s1, "n": len(xs), "max_dev": max(x["dev"] for x in xs),
                             "mean_plddt": round(float(np.mean(pls)), 1) if pls else None,
                             "terminus": s1 - s0 <= 5 and any(abs(s0 - e) <= 3 or abs(s1 - e) <= 3 for e in ends)})
        cmp = {"key": f"{ma['sid']}:{ca}|{mb['sid']}:{cb}", "a": ma["sid"], "b": mb["sid"], "chain_a": ca, "chain_b": cb,
               "label_a": ma["label"], "label_b": mb["label"], "ce": ce, "rmsd_mapped": round(rms, 3), "n_pairs": len(keep),
               "numbering": numbering, "deviation": dev, "segments": seg_rows, "superposed_file": f"structures/{fname}",
               "time": round(time.time() - t0, 2)}
        cmp["how"] = ("CE (combinatorial extension, Shindyalov & Bourne 1998; Bio.PDB.cealign) finds the best structural "
                      "alignment of the two chains from their CA atoms alone and reports its RMSD. Then residues are paired by "
                      f"{'UniProt position' if numbering == 'UniProt' else 'sequence alignment'}, the second structure is "
                      "superposed on those pairs (Bio.PDB.Superimposer, least squares on CA), and the distance between paired CA "
                      "atoms is plotted along the sequence with AlphaFold confidence on top. Deviations > 3 Å are marked.")
        cmp["yours"] = _compare_yours(cmp, res)
        res["comparisons"] = [c for c in res.get("comparisons", []) if c["key"] != cmp["key"]] + [cmp]
        save_result(item, res)
        return cmp


def _compare_yours(c: dict, res: dict) -> str:
    parts = []
    if c["ce"].get("rmsd") is not None:
        parts.append(f"CE overall RMSD {c['ce']['rmsd']:.2f} Å.")
    else:
        parts.append(c["ce"].get("error", "CE did not run."))
    parts.append(f"Superposed on {c['n_pairs']} residue pairs matched by {'UniProt position' if c['numbering'] == 'UniProt' else 'sequence'}: "
                 f"RMSD {c['rmsd_mapped']:.2f} Å.")
    segs = c["segments"]
    nbig = sum(1 for x in c["deviation"] if x["dev"] > 3)
    if c["n_pairs"] and nbig > 0.5 * c["n_pairs"]:
        parts.append(f"{nbig} of {c['n_pairs']} matched residues deviate by more than 3 Å: one rigid superposition cannot fit "
                     "both structures, which usually means a hinge or a different packing of segments (for a long helix, a small "
                     "kink swings the far end a long way). CE's RMSD, over the part that does align, is the better measure of "
                     "fold similarity; the per-residue curve here mostly shows that global mismatch rather than local "
                     "differences. To compare domains, load entries that cover one domain each.")
        return " ".join(parts)
    if c["ce"].get("rmsd") is not None and c["rmsd_mapped"] > c["ce"]["rmsd"] + 1.5:
        parts.append("The sequence-matched RMSD is much larger than CE's: part of the chain moves as a block (a hinge or a "
                     "differently packed domain), which CE leaves out of its core alignment.")
    if not segs:
        parts.append("No stretch deviates by more than 3 Å — the two agree closely over the shared residues.")
    exp = [s for s in segs if s["mean_plddt"] is not None and s["mean_plddt"] < 70]
    real = [s for s in segs if s["mean_plddt"] is not None and s["mean_plddt"] >= 70]
    unk = [s for s in segs if s["mean_plddt"] is None]
    fmt = lambda xs: "; ".join((f"{s['start']}–{s['end']}" if s["start"] != s["end"] else f"{s['start']}") + f" (up to {s['max_dev']:.1f} Å"
                               + (f", mean pLDDT {s['mean_plddt']:.0f}" if s["mean_plddt"] is not None else "")
                               + (", at the end of the shared region — chain ends are often frayed" if s.get("terminus") else "") + ")"
                               for s in xs[:6])  # noqa: E731
    if exp:
        parts.append(f"Large deviation where AlphaFold is not confident — expected: {fmt(exp)}.")
    if real:
        parts.append(f"Large deviation in confident regions — a real conformational difference, worth explaining (a bound ligand, "
                     f"a partner protein, or crystal contacts in the experimental entry): {fmt(real)}.")
    if unk:
        parts.append(f"Deviating stretches (no confidence data): {fmt(unk)}.")
    return " ".join(parts)


# ---------------------------------------------------------------- exports
def export_structure(item: str, sid: str, bmode: str = "flag", fmt: str = "cif") -> tuple[Path, str]:
    """The structure with its B-factor column replaced: 'am' = mean AlphaMissense score per residue (0–1),
    'flag' = 1 at mapped variant residues and 0 elsewhere. Returns (path, download name)."""
    d = item_dir(item)
    res = load_result(item)
    meta = _structure_meta(res, sid)
    s2, _ = st.load(d / meta["file"], "export")              # fresh parse: the cached one must not be modified
    vals = {}
    ch = meta["main_chain"]
    rows = rows_for(d, meta["sid"], ch)
    if bmode == "am":
        means = (res.get("alphamissense") or {}).get("residue_mean")
        if not means:
            raise UserFacingError("AlphaMissense is not available for this protein, so there is nothing to write into the B-factor column.")
        for r in rows:
            if r["pos"] and r["pos"] <= len(means) and means[r["pos"] - 1] is not None:
                vals[r["key"]] = means[r["pos"] - 1]
    else:
        vres = res.get("variants") or {}
        if not vres or vres.get("structure") != meta["sid"]:
            raise UserFacingError(f"Map variants on {meta['label']} first; the flag column marks the mapped variant residues.")
        for r in vres["rows"]:
            if r.get("key"):
                vals[r["key"]] = 1.0
    st.replace_bfactor(s2, vals, 0.0)
    for m in list(s2)[1:]:
        s2.detach_child(m.id)
    fmt = "pdb" if fmt == "pdb" else "cif"
    name = f"{res.get('accession') or Path(meta['file']).stem}_{meta['sid']}_{'alphamissense' if bmode == 'am' else 'variants'}.{fmt}"
    out = d / "exports" / name
    out.parent.mkdir(exist_ok=True)
    try:
        st.write(s2, out, fmt)
    except Exception as e:  # noqa: BLE001 - PDB format limits (atoms > 99,999, two-letter chains)
        if fmt == "pdb":
            raise UserFacingError(f"This structure does not fit the old PDB format ({e}); download mmCIF instead.") from None
        raise
    return out, name


def _script_residues(res: dict, meta: dict) -> list[dict]:
    vres = res.get("variants") or {}
    if not vres or vres.get("structure") != meta["sid"]:
        return []
    return [{"label": r["label"], "chain": r["chain"], "num": r["num"]} for r in vres["rows"] if r.get("key")]


def viewer_script(item: str, program: str = "pymol", sid: str | None = None, variants: str | None = None) -> tuple[str, str]:
    """A PyMOL .pml or ChimeraX .cxc that loads the structure, colours it and shows the variants as labelled sticks."""
    res = load_result(item)
    meta = _structure_meta(res, sid)
    if variants:
        vres = map_variants(item, variants, meta["sid"], save=False)
        vs = [{"label": r["label"], "chain": r["chain"], "num": r["num"]} for r in vres["rows"] if r.get("key")]
    else:
        vs = _script_residues(res, meta)
    merged = {}                                   # one label per residue: R482W and R482Q → "R482W/Q"
    for v in vs:
        k = (v["chain"], v["num"])
        if k in merged:
            merged[k]["label"] += "/" + v["label"][-1]
        else:
            merged[k] = dict(v)
    vs = list(merged.values())
    fname = Path(meta["file"]).name
    url = (res.get("alphafold") or {}).get("cif_url") if meta["sid"] == "AF" else None
    pid = meta.get("pdb_id")
    obj = "sb_" + "".join(ch if ch.isalnum() else "_" for ch in Path(fname).stem)   # PyMOL names must not start with a digit
    plddt = meta["bkind"] == "plddt"
    if program == "chimerax":
        L = [f"# ChimeraX script from Structure Bench {VERSION} — {res.get('name')} — {meta['label']}",
             f"# Put this file next to {fname} (download it from the app), or replace the open line with:",
             f"#   open {pid}" if pid else (f"#   open {url}" if url else "#   open <your file>"),
             f"open {fname}", "cartoon", "hide atoms",
             "color bfactor palette alphafold" if plddt else "color bychain"]
        for v in vs:
            spec = f"/{v['chain']}:{v['num']}"
            L += [f"show {spec} atoms", f"style {spec} stick", f"color {spec} magenta target a",
                  f"label {spec} text \"{v['label']}\" height 1.2"]
        if vs:
            L.append("view " + " ".join(f"/{v['chain']}:{v['num']}" for v in vs))
        return "\n".join(L) + "\n", f"{obj}.cxc"
    L = [f"# PyMOL script from Structure Bench {VERSION} — {res.get('name')} — {meta['label']}",
         f"# Put this file next to {fname} (download it from the app), or replace the load line with:",
         f"#   fetch {pid}, {obj}" if pid else (f"#   load {url}, {obj}" if url else "#   load <your file>"),
         f"load {fname}, {obj}", "hide everything", "show cartoon", "set cartoon_transparency, 0"]
    if plddt:
        L += ["# AlphaFold confidence colours (pLDDT in the B-factor column)",
              "set_color af_vhigh, [0.000, 0.325, 0.839]", "set_color af_conf, [0.396, 0.796, 0.953]",
              "set_color af_low, [1.000, 0.859, 0.075]", "set_color af_vlow, [1.000, 0.490, 0.271]",
              f"color af_vlow, {obj}", f"color af_low, {obj} and b > 50", f"color af_conf, {obj} and b > 70",
              f"color af_vhigh, {obj} and b > 90"]
    else:
        L.append(f"util.cbc {obj}")
    for v in vs:
        sel = f"{obj} and chain {v['chain']} and resi {v['num']}"
        L += [f"show sticks, {sel} and not name N+C+O", f"color magenta, {sel} and elem C",
              f"show spheres, {sel} and name CA", "set sphere_scale, 0.4",
              f"label {sel} and name CA, \"{v['label']}\""]
    L += ["set label_size, 18", "set label_color, black", "bg_color white"]
    if vs:
        L.append("orient " + " or ".join(f"({obj} and chain {v['chain']} and resi {v['num']})" for v in vs))
    return "\n".join(L) + "\n", f"{obj}.pml"


# ---------------------------------------------------------------- for the assistant
def item_summary(item: str) -> dict:
    """A compact view of an item for the assistant: no per-residue arrays."""
    res = load_result(item)
    u = res.get("uniprot") or {}
    out = {k: res.get(k) for k in ("item", "kind", "name", "title", "accession", "organism", "length", "tiles", "flags")}
    if u:
        out["uniprot"] = {"function": u.get("function", "")[:1500], "location": u.get("location"),
                          "diseases": [{k: x[k] for k in ("name", "acronym")} for x in u.get("diseases", [])],
                          "domains": [f"{f['type']}: {f['description']} {f['start']}-{f['end']}" for f in u["groups"]["Domains & regions"]][:40],
                          "n_variants": len(u.get("variants", [])), "n_pdb": len(u.get("pdb", []))}
    af = res.get("alphafold")
    if af:
        out["alphafold"] = {k: af[k] for k in ("entry", "version", "mean_plddt", "frac_very_high", "frac_confident", "frac_low", "frac_very_low", "ss_method")}
    if res.get("pae"):
        out["pae"] = res["pae"]
    out["texts"] = {k: v.get("yours") for k, v in (res.get("texts") or {}).items() if v.get("yours")}
    out["structures"] = [{k: s.get(k) for k in ("sid", "label", "kind", "main_chain", "method", "resolution", "partners")}
                         for s in res["structures"].values()]
    if res.get("variants"):
        v = res["variants"]
        out["variants"] = {"structure": v["structure"], "yours": v["yours"], "rows": [
            {k: r.get(k) for k in ("label", "error", "kind", "plddt", "bfactor", "rsa", "rsa_class", "ss_word", "domain", "am_score",
                                   "am_class", "uniprot_known", "interface", "ligand_contacts", "unresolved")}
            | {"flags": [f["text"] for f in r.get("flags", [])]} for r in v["rows"]]}
    if res.get("comparisons"):
        out["comparisons"] = [{k: c[k] for k in ("a", "b", "chain_a", "chain_b", "ce", "rmsd_mapped", "n_pairs", "segments", "yours")}
                              for c in res["comparisons"]]
    return out


# ---------------------------------------------------------------- methods
CITATIONS = {
    "biopython": "Cock PJA, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25:1422–1423 (2009).",
    "alphafold": "Jumper J, et al. Highly accurate protein structure prediction with AlphaFold. Nature 596:583–589 (2021).",
    "afdb": "Varadi M, et al. AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. Nucleic Acids Res 52:D368–D375 (2024).",
    "alphamissense": "Cheng J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. Science 381:eadg7492 (2023).",
    "uniprot": "The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. Nucleic Acids Res 53:D609–D617 (2025).",
    "shrake": "Shrake A, Rupley JA. Environment and exposure to solvent of protein atoms. Lysozyme and insulin. J Mol Biol 79:351–371 (1973).",
    "tien": "Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO. Maximum allowed solvent accessibilites of residues in proteins. PLoS One 8:e80635 (2013).",
    "ce": "Shindyalov IN, Bourne PE. Protein structure alignment by incremental combinatorial extension (CE) of the optimal path. Protein Eng 11:739–747 (1998).",
    "blosum": "Henikoff S, Henikoff JG. Amino acid substitution matrices from protein blocks. PNAS 89:10915–10919 (1992).",
    "pdb": "Berman HM, et al. The Protein Data Bank. Nucleic Acids Res 28:235–242 (2000).",
    "3dmol": "Rego N, Koes D. 3Dmol.js: molecular visualization with WebGL. Bioinformatics 31:1322–1324 (2015).",
    "dssp": "Kabsch W, Sander C. Dictionary of protein secondary structure. Biopolymers 22:2577–2637 (1983).",
}


def methods(res: dict) -> list[list[str]]:
    bv = versions("biopython").get("biopython", "?")
    M = []
    if res["kind"] == "protein":
        M.append(["Protein annotation", f"The UniProtKB entry {res.get('accession')} was retrieved from the UniProt REST API "
                  "(rest.uniprot.org, JSON); gene symbols were resolved with Bio.UniProt.search "
                  f"(gene_exact, organism_id, reviewed) in Biopython {bv}. Feature positions refer to the canonical isoform. "
                  "[UniProt Consortium 2025]"])
        af = res.get("alphafold")
        if af:
            M.append(["Predicted structure", f"The AlphaFold model {af['entry']} (database version {af['version']}) and its predicted "
                      "aligned error were retrieved from the AlphaFold Protein Structure Database (entries listed with "
                      "Bio.PDB.alphafold_db.get_predictions) and parsed with Bio.PDB.MMCIFParser. Per-residue pLDDT was read from "
                      "the CA B-factor. Mean PAE between UniProt domains was computed from the symmetrised PAE matrix. "
                      "[Jumper 2021; Varadi 2024]"])
        if (res.get("alphamissense") or {}).get("available"):
            M.append(["Variant effect prediction", "AlphaMissense pathogenicity scores and classes for every possible missense "
                      "substitution of the canonical isoform were taken from the AlphaFold DB substitution table "
                      "(likely benign < 0.34, ambiguous 0.34–0.564, likely pathogenic > 0.564). [Cheng 2023]"])
    M.append(["Solvent accessibility", "Per-residue solvent-accessible surface area was computed with the Shrake–Rupley "
              "algorithm (Bio.PDB.SASA.ShrakeRupley, probe 1.4 Å, 100 points per atom; hydrogens and waters removed) for each "
              "chain alone and together with every chain and ligand within 8 Å, and divided by the theoretical maximum of Tien "
              "et al. 2013 (Bio.PDB.DSSP.residue_max_acc['Wilke']). Residues with relative SASA < 20 % were called buried and "
              "> 50 % exposed. [Shrake & Rupley 1973; Tien 2013]"])
    M.append(["Contacts and secondary structure", "Neighbours were all residues with a heavy atom within 5 Å of any heavy atom of "
              "the residue (Bio.PDB.NeighborSearch). "
              + ("Secondary structure was assigned by DSSP (Bio.PDB.DSSP). [Kabsch & Sander 1983]" if st.dssp_available() else
                 "The DSSP program was not installed, so secondary structure was approximated from backbone phi/psi angles "
                 "(Bio.PDB.Polypeptide.get_phi_psi_list): helix-like if phi −160…−20° and psi −120…50° for ≥ 4 consecutive "
                 "residues, extended if phi −180…−45° and psi ≥ 90° or ≤ −150° for ≥ 3 residues. This is a heuristic.")])
    M.append(["Structure mapping and comparison", "Chains of experimental or uploaded structures were extracted with "
              "Bio.PDB.PPBuilder and mapped to UniProt numbering by global alignment (Bio.Align.PairwiseAligner, BLOSUM62, gap "
              "open −10, extend −0.5, end gaps free). Structures were compared with Bio.PDB.cealign.CEAligner (CE) and, on "
              "residue pairs matched by UniProt position, superposed with Bio.PDB.Superimposer to give per-residue CA deviations. "
              "[Shindyalov & Bourne 1998; Henikoff 1992]"])
    M.append(["Software", f"Structure Bench {VERSION} (Paul H. Kim, Ph.D.), Biopython {bv} [Cock 2009], NumPy; 3D view with "
              "3Dmol.js 2.5.5 [Rego & Koes 2015]."])
    return M


def citations(res: dict) -> list[str]:
    keys = ["uniprot", "alphafold", "afdb", "alphamissense", "shrake", "tien", "biopython", "ce", "blosum", "pdb", "3dmol"]
    if res["kind"] != "protein":
        keys = ["shrake", "tien", "biopython", "ce", "blosum", "3dmol"]
    if st.dssp_available():
        keys.append("dssp")
    return [CITATIONS[k] for k in keys]
