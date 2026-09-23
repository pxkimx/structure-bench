"""Protein variant notation: parsing, reference-residue checks, and the plain-language interpretation flags."""
from __future__ import annotations

import re

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H",
       "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
       "Tyr": "Y", "Val": "V", "Sec": "U", "Pyl": "O", "Ter": "*"}
AA1 = set("ACDEFGHIKLMNPQRSTVWYUO")
NAME = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His",
        "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
        "Y": "Tyr", "V": "Val", "U": "Sec", "O": "Pyl", "*": "stop"}
# residue volume (Å³, Zamyatnin 1972) and charge at pH 7, for the size/charge wording of the heuristics
VOLUME = {"G": 60.1, "A": 88.6, "S": 89.0, "C": 108.5, "D": 111.1, "P": 112.7, "N": 114.1, "T": 116.1, "E": 138.4,
          "V": 140.0, "Q": 143.8, "H": 153.2, "M": 162.9, "I": 166.7, "L": 166.7, "K": 168.6, "R": 173.4, "F": 189.9,
          "Y": 193.6, "W": 227.8}
CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1}
HYDROPHOBIC = set("AVILMFWC")

_TOKEN = re.compile(r"[\s,;]+")
_ONE = re.compile(r"^\(?([A-Z])(\d+)([A-Z*=]|Ter|TER)\)?$")
_THREE = re.compile(r"^(?:p\.)?\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=|\*)\)?$")


def split(text: str) -> list[str]:
    return [t for t in _TOKEN.split(text or "") if t.strip()]


def parse(token: str) -> dict:
    """One variant → {input, ref, pos, alt, kind} or {input, error}.

    Accepted: R482W, p.R482W, p.Arg482Trp, Arg482Trp, p.(Arg482Trp), G608G, p.Gly608=, p.Gly608Gly, R482*, p.Arg482Ter.
    """
    t = token.strip()
    m = _THREE.match(t)
    if m:
        ref = AA3.get(m.group(1).capitalize())
        alt_raw = m.group(3)
        alt = ref if alt_raw == "=" else "*" if alt_raw == "*" else AA3.get(alt_raw.capitalize())
        if ref and alt:
            return _mk(t, ref, int(m.group(2)), alt)
    body = t[2:] if t[:2].lower() == "p." else t
    for cand in (body, body.upper()):
        m = _ONE.match(cand)
        if m and m.group(1) in AA1:
            a = m.group(3)
            alt = m.group(1) if a == "=" else "*" if a in ("*", "X", "Ter", "TER") else a
            if alt in AA1 or alt == "*":
                return _mk(t, m.group(1), int(m.group(2)), alt)
    hint = ""
    if re.search(r"(fs|del|ins|dup|ext)", t, re.I):
        hint = " Frameshifts, insertions and deletions cannot be placed on a single residue here."
    elif re.match(r"^c\.", t, re.I):
        hint = " That is a DNA (c.) change; give the protein change (p.) instead."
    return {"input": t, "error": f"'{t}' is not a protein substitution this app can read (examples: R482W, p.Arg482Trp, G608G).{hint}"}


def _mk(t, ref, pos, alt):
    kind = "synonymous" if ref == alt else "nonsense" if alt == "*" else "missense"
    return {"input": t, "ref": ref, "pos": pos, "alt": alt, "kind": kind, "label": f"{ref}{pos}{alt}"}


def parse_many(text: str) -> list[dict]:
    seen, out = set(), []
    for tok in split(text):
        v = parse(tok)
        k = v.get("label") or v["input"]
        if k in seen:
            continue
        seen.add(k)
        out.append(v)
    return out


def check_reference(v: dict, seq: str, signal_len: int = 0, isoforms: list[dict] | None = None,
                    chain_offsets: list[tuple[str, int]] | None = None, what: str = "the UniProt sequence") -> dict:
    """Does the reference residue match the sequence at that position? If not, say what would explain it.

    Tries: the mature-chain numbering (signal peptide / propeptide removed), numbering without the initiator
    methionine, and every isoform whose sequence is known. Returns {"ok": bool, "message": str, "suggest": [...]}.
    """
    pos, ref = v["pos"], v["ref"]
    n = len(seq)
    if 1 <= pos <= n and seq[pos - 1] == ref:
        return {"ok": True, "message": ""}
    got = seq[pos - 1] if 1 <= pos <= n else None
    msg = (f"Position {pos} in {what} is {NAME.get(got, got)} ({got}), not {NAME.get(ref, ref)} ({ref})." if got
           else f"Position {pos} is beyond the end of {what} ({n} residues).")
    sugg = []
    offs = []
    if signal_len:
        offs.append((f"mature-chain numbering (after removing the {signal_len}-residue signal peptide)", signal_len))
    for lab, o in chain_offsets or []:
        if o and o != signal_len:
            offs.append((lab, o))
    offs.append(("numbering that skips the initiator methionine", 1))
    for lab, o in offs:
        p2 = pos + o
        if 1 <= p2 <= n and seq[p2 - 1] == ref:
            sugg.append(f"It matches if the variant uses {lab}: in UniProt numbering it would be "
                        f"{ref}{p2}{v['alt'] if v['kind'] != 'synonymous' else ref}.")
    for iso in isoforms or []:
        s = iso.get("sequence") or ""
        if 1 <= pos <= len(s) and s[pos - 1] == ref:
            sugg.append(f"It matches isoform {iso.get('accession') or iso.get('entry')} at position {pos} — the variant "
                        "may be numbered on that isoform rather than the canonical sequence.")
    if not sugg:
        sugg.append("Check the transcript/isoform the variant was called on (the canonical UniProt isoform is used here), "
                    "and whether the numbering is for the precursor or the mature protein.")
    return {"ok": False, "message": msg, "suggest": sugg}


def change_words(ref: str, alt: str) -> str:
    """'larger, positive → neutral' style description of a substitution (heuristic, from volume and charge)."""
    if alt == "*" or ref == alt:
        return ""
    bits = []
    dv = VOLUME.get(alt, 0) - VOLUME.get(ref, 0)
    if abs(dv) >= 30:
        bits.append(f"{'larger' if dv > 0 else 'smaller'} side chain ({dv:+.0f} Å³)")
    cr, ca = CHARGE.get(ref, 0), CHARGE.get(alt, 0)
    if cr != ca:
        w = {1: "positive", -1: "negative", 0: "neutral"}
        bits.append(f"charge {w[cr]} → {w[ca]}")
    if (ref in HYDROPHOBIC) != (alt in HYDROPHOBIC):
        bits.append("hydrophobic → polar" if ref in HYDROPHOBIC else "polar → hydrophobic")
    if alt == "P":
        bits.append("proline introduced")
    if ref == "G":
        bits.append("glycine lost")
    if alt == "G":
        bits.append("glycine introduced")
    if ref == "C" or alt == "C":
        bits.append("cysteine " + ("lost" if ref == "C" else "introduced"))
    return "; ".join(bits)


def interpret(row: dict) -> list[dict]:
    """Heuristic reading of one mapped variant → list of {level, text}. Deliberately simple and labelled as such."""
    out = []
    kind = row.get("kind")
    if kind == "synonymous":
        return [{"level": "info", "text": "No amino-acid change — structure cannot assess this. Synonymous variants can "
                                          "still act on the RNA, for example by creating or destroying a splice site."}]
    if row.get("unresolved"):
        return [{"level": "info", "text": "This residue is not modelled in the chosen structure, so nothing structural "
                                          "can be said from it; try the AlphaFold model or another PDB entry."}]
    if kind == "nonsense":
        out.append({"level": "warn", "text": f"Stop codon: everything after residue {row['pos']} is missing from the protein "
                                             "(if the transcript escapes nonsense-mediated decay). The structure shows what would be lost."})
        return out
    pl = row.get("plddt")
    rsa = row.get("rsa")
    ss = row.get("ss_word", "")
    chg = row.get("change", "")
    if pl is not None and pl < 50:
        out.append({"level": "info", "text": "In a low-confidence region (pLDDT < 50) — structural reading not meaningful "
                                             "here; the region is probably disordered in isolation."})
        return out
    confident = pl is None or pl >= 70           # experimental structure (pl None) counts as a resolved fold
    if pl is not None and pl < 70:
        out.append({"level": "info", "text": f"pLDDT {pl:.0f} (low, 50–70): the local backbone is only roughly placed; read the "
                                             "burial and contacts below with caution."})
    fold = row.get("rsa_alone") if row.get("rsa_alone") is not None else rsa    # burial by the chain's own fold
    if fold is not None and fold < 0.20:
        if confident:
            where = "a confident fold" if pl is not None else "the folded, resolved structure"
            out.append({"level": "warn", "text": f"Buried residue in {where}; a charge/size change here is likely to be "
                                                 "destabilising." + (f" This change: {chg}." if chg else "")})
        else:
            out.append({"level": "info", "text": "Buried in the model, but the model is not confident here."})
    elif fold is not None and fold > 0.50:
        if confident:
            where = "a confident region" if pl is not None else "the resolved structure"
            out.append({"level": "info", "text": f"Exposed residue in {where}: less likely to destabilise the fold, "
                                                 "but it may sit on a binding surface — check partner and ligand contacts."})
    elif fold is not None:
        out.append({"level": "info", "text": "Partly buried (relative SASA 20–50 %): an intermediate position; the effect "
                                             "depends on the side chain's packing." + (f" This change: {chg}." if chg else "")})
    if row.get("coiled_coil") and row.get("structure") == "AF":
        out.append({"level": "info", "text": "UniProt annotates a coiled coil here. Coiled coils form only as dimers or larger "
                                             "assemblies, and this is a single-chain model — residues that pack against the partner "
                                             "helix look exposed here, so the burial reading underestimates it."})
    if row.get("interface"):
        extra = ""
        if row.get("rsa_alone") is not None and rsa is not None:
            extra = f" (relative SASA {row['rsa_alone']:.0%} for this chain alone, {rsa:.0%} with its partners)"
        note = f" The partner is {row['interface_note']}." if row.get("interface_note") else \
            " A PDB entry is one crystal's contents — an interface can also be a crystal contact."
        out.append({"level": "warn", "text": f"At an interface: it loses surface when {row['interface']} is present{extra}, so it may "
                                             "affect that interaction." + note})
    if row.get("ligand_contacts"):
        out.append({"level": "warn", "text": f"Within 5 Å of ligand {row['ligand_contacts']} in this structure."})
    if row.get("alt") == "P" and row.get("ss") == "H":
        out.append({"level": "warn", "text": f"Proline introduced into a {ss} segment — proline usually breaks helices."})
    return out
