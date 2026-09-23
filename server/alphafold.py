"""AlphaFold Database: predictions for an accession, the model file, PAE and AlphaMissense substitutions."""
from __future__ import annotations

import csv
import io
import json

import numpy as np

from . import net

AF_COLORS = [(90, "#0053D6", "Very high (pLDDT > 90)"), (70, "#65CBF3", "Confident (70–90)"),
             (50, "#FFDB13", "Low (50–70)"), (0, "#FF7D45", "Very low (< 50)")]
AM_CLASS = {"LBen": "likely benign", "Amb": "ambiguous", "LPath": "likely pathogenic",
            "likely_benign": "likely benign", "ambiguous": "ambiguous", "likely_pathogenic": "likely pathogenic"}


def plddt_band(v: float | None) -> str:
    if v is None:
        return "—"
    return "very high" if v > 90 else "confident" if v >= 70 else "low" if v >= 50 else "very low"


def predictions(acc: str, refresh=False) -> tuple[list[dict], str]:
    """All AlphaFold DB entries for an accession (canonical + isoforms), via Bio.PDB.alphafold_db.get_predictions."""
    from Bio.PDB import alphafold_db
    acc = acc.strip().upper()
    return net.cached_call(f"alphafold/{acc}.predictions.json", lambda: list(alphafold_db.get_predictions(acc)),
                           f"the AlphaFold DB entry for {acc}", refresh)


def pick(preds: list[dict], acc: str) -> tuple[dict | None, list[dict]]:
    """The canonical model (entryId AF-<acc>-F1) and the isoform models, listed as options."""
    canon = next((p for p in preds if p.get("entryId") == f"AF-{acc}-F1"), None)
    if canon is None:          # some entries are keyed by modelEntityId only
        canon = next((p for p in preds if p.get("modelEntityId") == f"AF-{acc}-F1"), None)
    iso = []
    for p in preds:
        if p is canon:
            continue
        iso.append({"entry": p.get("entryId") or p.get("modelEntityId"), "accession": p.get("uniprotAccession"),
                    "length": len(p.get("sequence") or p.get("uniprotSequence") or ""), "mean_plddt": p.get("globalMetricValue"),
                    "description": p.get("uniprotDescription", ""), "cif": p.get("cifUrl"),
                    "sequence": p.get("sequence") or p.get("uniprotSequence") or ""})
    return canon, iso


def _name(url: str) -> str:
    return url.rstrip("/").split("/")[-1]


def model_file(pred: dict, refresh=False) -> tuple[bytes, str, str]:
    url = pred["cifUrl"]
    data, src = net.fetch(url, f"alphafold/{_name(url)}", f"the AlphaFold model {pred.get('entryId')}", refresh)
    return data, src, _name(url)


def pae_json(pred: dict, refresh=False):
    url = pred.get("paeDocUrl")
    if not url:
        raise net.NetError("This AlphaFold entry lists no PAE file.")
    return net.fetch_json(url, f"alphafold/{_name(url)}", "the PAE (predicted aligned error) file", refresh)


def parse_pae(obj) -> tuple[np.ndarray, float]:
    """PAE matrix (N×N, Å) and the maximum possible value, from either AlphaFold DB layout.

    New (v3+):  [{"predicted_aligned_error": [[...], ...], "max_predicted_aligned_error": 31.75}]
    Old (v1/2): [{"residue1": [1,1,...], "residue2": [1,2,...], "distance": [...], "max_predicted_aligned_error": 31.75}]
    A bare dict (not wrapped in a list) is accepted for both.
    """
    if isinstance(obj, list):
        if not obj:
            raise ValueError("The PAE file is empty.")
        obj = obj[0]
    if not isinstance(obj, dict):
        raise ValueError("The PAE file is not in a layout this app recognises.")
    mx = float(obj.get("max_predicted_aligned_error") or 31.75)
    if "predicted_aligned_error" in obj:
        m = np.asarray(obj["predicted_aligned_error"], dtype=float)
    elif "distance" in obj and "residue1" in obj and "residue2" in obj:
        r1 = np.asarray(obj["residue1"], dtype=int)
        r2 = np.asarray(obj["residue2"], dtype=int)
        n = int(max(r1.max(), r2.max()))
        m = np.full((n, n), np.nan)
        m[r1 - 1, r2 - 1] = np.asarray(obj["distance"], dtype=float)
    else:
        raise ValueError("The PAE file has neither 'predicted_aligned_error' nor 'residue1/residue2/distance'.")
    if m.ndim != 2 or m.shape[0] != m.shape[1]:
        raise ValueError(f"The PAE matrix is not square ({m.shape}).")
    return m, mx


def pae_display(m: np.ndarray, maxn: int = 700) -> dict:
    """The matrix for the canvas: rounded to 0.1 Å, block-averaged when larger than maxn so the page stays light."""
    n = m.shape[0]
    b = max(1, int(np.ceil(n / maxn)))
    if b > 1:
        k = int(np.ceil(n / b))
        pad = np.full((k * b, k * b), np.nan)
        pad[:n, :n] = m
        d = np.nanmean(pad.reshape(k, b, k, b), axis=(1, 3))
    else:
        d = m
    return {"n": n, "bin": b, "size": d.shape[0], "values": np.round(np.nan_to_num(d, nan=-1), 1).ravel().tolist()}


def am_table(pred: dict, refresh=False) -> tuple[dict, str]:
    """AlphaMissense for every possible missense change: {'R482W': (0.93, 'likely pathogenic'), ...}."""
    url = pred.get("amAnnotationsUrl")
    if not url:
        raise net.NetError("AlphaFold DB has no AlphaMissense table for this entry (it exists for canonical human "
                           "UniProt isoforms only).")
    data, src = net.fetch(url, f"alphafold/{_name(url)}", "the AlphaMissense table", refresh)
    return parse_am(data.decode("utf-8", "replace")), src


def parse_am(text: str) -> dict:
    out = {}
    rd = csv.reader(io.StringIO(text))
    head = next(rd, None) or []
    try:
        iv, isc, icl = head.index("protein_variant"), head.index("am_pathogenicity"), head.index("am_class")
    except ValueError:
        iv, isc, icl = 0, 1, 2
    for row in rd:
        if len(row) <= max(iv, isc, icl):
            continue
        try:
            out[row[iv].strip()] = (float(row[isc]), AM_CLASS.get(row[icl].strip(), row[icl].strip()))
        except ValueError:
            continue
    return out


def am_residue_means(am: dict, length: int) -> list[float | None]:
    """Mean AlphaMissense score over the 19 substitutions at each position (a per-residue sensitivity summary)."""
    acc = [[] for _ in range(length + 1)]
    for k, (s, _) in am.items():
        try:
            p = int(k[1:-1])
        except ValueError:
            continue
        if 1 <= p <= length:
            acc[p].append(s)
    return [round(float(np.mean(v)), 4) if v else None for v in acc[1:]]


def summarize_pae(m: np.ndarray, blocks: list[dict]) -> dict:
    """Mean PAE within and between blocks (1-based inclusive ranges). Symmetrised: PAE(i,j) and PAE(j,i) differ."""
    n = m.shape[0]
    rows = []
    for b in blocks:
        s, e = max(1, b["start"]), min(n, b["end"])
        if e - s < 4:
            continue
        rows.append({**b, "start": s, "end": e, "within": round(float(np.nanmean(m[s - 1:e, s - 1:e])), 1)})
    pairs = []
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            a, b = rows[i], rows[j]
            v = np.nanmean(np.concatenate([m[a["start"] - 1:a["end"], b["start"] - 1:b["end"]].ravel(),
                                           m[b["start"] - 1:b["end"], a["start"] - 1:a["end"]].ravel()]))
            pairs.append({"a": a["name"], "b": b["name"], "mean": round(float(v), 1),
                          "call": "confident" if v < 10 else "uncertain" if v < 20 else "not determined"})
    return {"blocks": rows, "pairs": pairs}


def confident_segments(plddt: list[float | None], minlen: int = 30, thr: float = 70) -> list[dict]:
    """Runs of pLDDT ≥ thr at least minlen long (gaps ≤ 3 bridged) — used as PAE blocks when UniProt has < 2 domains."""
    segs, s, gap = [], None, 0
    for i, v in enumerate(plddt + [None] * 5):
        ok = v is not None and v >= thr
        if ok:
            if s is None:
                s = i
            gap, last = 0, i
        elif s is not None:
            gap += 1
            if gap > 3:
                if last - s + 1 >= minlen:
                    segs.append({"name": f"Confident segment {s+1}–{last+1}", "start": s + 1, "end": last + 1})
                s, gap = None, 0
    return segs


def load_json_file(p):
    return json.loads(open(p).read())
