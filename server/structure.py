"""Structure computations with Bio.PDB: reading, chain sequences, sequence→UniProt mapping, solvent
accessibility, neighbours, secondary structure, superposition, and writing files.

Pure Biopython + numpy + stdlib. Residues are addressed by a key string "A:482" (chain:number, plus an
insertion code when there is one, "A:52A") so every result is JSON-serialisable.
"""
from __future__ import annotations

import gzip
import shutil
import subprocess
import warnings
from pathlib import Path

import numpy as np

MAX_ACC_SOURCE = "Wilke"          # Tien et al. 2013 theoretical maxima, as tabulated in Bio.PDB.DSSP
WATER = {"HOH", "WAT", "DOD", "H2O"}
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
         "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
         "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O"}


def max_acc() -> dict:
    from Bio.PDB.DSSP import residue_max_acc
    return residue_max_acc[MAX_ACC_SOURCE]


def rkey(res) -> str:
    chain = res.get_parent().id
    _, num, ic = res.id
    return f"{chain}:{num}{ic.strip()}"


def one(res) -> str:
    return THREE.get(res.get_resname(), "X")


# ---------------------------------------------------------------- reading
def _patch_binary_cif():
    """Biopython 1.88's BinaryCIF string decoder indexes the unique-string table with the lookup value -1 that
    marks a missing value ('?'). When a column is missing in every row — pdbx_PDB_ins_code in every AlphaFold DB
    .bcif — the table is empty and parsing fails with 'index -1 is out of bounds'. This replacement maps -1 to ''.
    """
    try:
        import Bio.PDB.binary_cif as bc
    except Exception:  # noqa: BLE001 - msgpack missing: .bcif simply is not offered
        return

    def string_array(column):
        encoding = column["data"]["encoding"][-1]
        offsets = bc._decode({"data": {"data": encoding["offsets"], "encoding": encoding["offsetEncoding"]}})
        text = encoding["stringData"]
        uniq = np.array([text[offsets[i]:offsets[i + 1]] for i in range(len(offsets) - 1)] + [""], dtype=object)
        lookups = np.asarray(bc._decode({"data": {"data": column["data"]["data"], "encoding": encoding["dataEncoding"]}}))
        lookups = np.where(lookups < 0, len(uniq) - 1, lookups)
        column["data"]["data"] = uniq[lookups]
        column["data"]["encoding"].pop()

    bc._decoders["StringArray"] = string_array


_patch_binary_cif()


def sniff(path: Path) -> str:
    """'cif' | 'pdb' | 'bcif' from the content (extensions lie, and .gz hides them)."""
    with open(path, "rb") as f:
        raw = f.read(4096)
    if raw[:2] == b"\x1f\x8b":
        try:
            with gzip.open(path, "rb") as f:
                raw = f.read(4096)
        except (OSError, EOFError):
            return "unknown"
    head = raw.lstrip()
    if head.startswith(b"data_") or b"\n_atom_site." in raw or b"loop_" in raw[:400]:
        return "cif"
    try:
        txt = head[:200].decode("ascii")
        if any(txt.startswith(k) for k in ("HEADER", "ATOM", "HETATM", "REMARK", "CRYST1", "MODEL", "TITLE", "COMPND", "EXPDTA")):
            return "pdb"
    except UnicodeDecodeError:
        pass
    if head[:1] in (b"\x80", b"\x81", b"\x82", b"\x83", b"\x84", b"\x85", b"\x86", b"\x87", b"\x88", b"\x89", b"\x8a",
                    b"\x8b", b"\x8c", b"\x8d", b"\x8e", b"\x8f", b"\xde", b"\xdf"):
        return "bcif"
    return "unknown"


def is_gzip(path: Path) -> bool:
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def text_copy(path: Path) -> Path:
    """The file itself, or for a gzipped one a decompressed copy next to it ('x.cif.gz' → 'x.unzipped.cif'). The copy
    gets its own name so it can never overwrite another upload called 'x.cif'; the viewer and the mmCIF name reader
    use it, since neither reads gzip."""
    path = Path(path)
    if not is_gzip(path):
        return path
    out = path.with_name(text_copy_name(path.name))
    if not out.exists() or out.stat().st_mtime < path.stat().st_mtime:
        with gzip.open(path, "rb") as fi, open(out, "wb") as fo:
            shutil.copyfileobj(fi, fo, 1 << 20)
    return out


def text_copy_name(name: str) -> str:
    base = name[:-3] if name.endswith(".gz") else name
    stem, dot, ext = base.rpartition(".")
    return f"{stem}.unzipped.{ext}" if dot and stem else f"{base}.unzipped"


def gunzipped_size(data: bytes, limit: int) -> int:
    """Size of gzip-compressed bytes once decompressed, counting only up to just past `limit` (no huge allocation)."""
    import zlib
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    n = 0
    for i in range(0, len(data), 1 << 20):
        n += len(d.decompress(data[i:i + (1 << 20)], limit + 1 - n if limit + 1 > n else 1))
        while d.unconsumed_tail and n <= limit:
            n += len(d.decompress(d.unconsumed_tail, limit + 1 - n))
        if n > limit:
            return n
    return n


def load(path: Path, sid: str = "s"):
    """Parse a .cif / .pdb / .bcif (optionally .gz) with the matching Bio.PDB parser, quietly."""
    from Bio.PDB import MMCIFParser, PDBParser
    path = Path(path)
    fmt = sniff(path)
    src = text_copy(path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if fmt == "cif":
            return MMCIFParser(QUIET=True).get_structure(sid, str(src)), "cif"
        if fmt == "pdb":
            return PDBParser(QUIET=True).get_structure(sid, str(src)), "pdb"
        if fmt == "bcif":
            from Bio.PDB.binary_cif import BinaryCIFParser
            s = BinaryCIFParser().get_structure(sid, str(src))
            # Some encoders (AlphaFold DB's python-mmcif) store coordinates and B-factors as strings, which
            # BinaryCIFParser passes through unchanged; everything downstream needs numbers.
            for a in all_atoms(s):
                a.coord = np.asarray(a.coord, dtype=np.float32)
                a.bfactor = float(a.bfactor) if a.bfactor not in (None, "", "?", ".") else 0.0
                a.occupancy = float(a.occupancy) if a.occupancy not in (None, "", "?", ".") else 1.0
            return s, "bcif"
    raise ValueError(f"{path.name} is not a PDB, mmCIF or BinaryCIF file this app can read.")


def write(structure, path: Path, fmt: str = "cif"):
    from Bio.PDB import MMCIFIO, PDBIO
    io = PDBIO() if fmt == "pdb" else MMCIFIO()
    io.set_structure(structure)
    io.save(str(path))


def entity_names(path: Path) -> dict:
    """Author chain id → molecule name from an mmCIF (_entity.pdbx_description). Empty for PDB/bcif files."""
    try:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        d = MMCIF2Dict(str(path))
        ids = d.get("_entity.id", [])
        desc = d.get("_entity.pdbx_description", [])
        names = dict(zip(ids, desc))
        out = {}
        for eid, strands in zip(d.get("_entity_poly.entity_id", []), d.get("_entity_poly.pdbx_strand_id", [])):
            for c in strands.split(","):
                out[c.strip()] = short_name(names.get(eid, ""))
        return out
    except Exception:  # noqa: BLE001 - names are a nicety
        return {}


def short_name(n: str) -> str:
    """'DNA (5'-D(*TP*TP*TP...)-3')' → 'DNA'; long names cut to 48 characters."""
    n = (n or "").strip().strip("'\" ")
    for k in ("DNA", "RNA"):
        if n.upper().startswith(k + " (") or n.upper().startswith(k + "("):
            return k
    return n if len(n) <= 48 else n[:47] + "…"


def foreign_contacts(model, chain_id: str, copies: set, cutoff: float = 4.5) -> int:
    """Heavy atoms of other molecules (not copies of this protein, not water) within `cutoff` Å of the chain."""
    from Bio.PDB import NeighborSearch
    other = [a for a in model.get_atoms() if a.get_parent().get_parent().id not in copies | {chain_id}
             and a.get_parent().get_resname() not in WATER and (a.element or "").upper() not in ("H", "D")]
    if not other:
        return 0
    ns = NeighborSearch(other)
    hit = set()
    for a in model[chain_id].get_atoms():
        for b in ns.search(a.coord, cutoff):
            hit.add(id(b))
    return len(hit)


def polymer_chains(model) -> list[dict]:
    """Protein chains of a model as PPBuilder sees them: residues in order, and their one-letter sequence.

    PPBuilder breaks a chain wherever consecutive residues are not peptide-bonded (missing density), so the
    pieces are concatenated per chain. Non-standard amino acids (MSE…) are kept, as X or their parent letter.
    """
    from Bio.PDB import PPBuilder
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pps = PPBuilder().build_peptides(model, aa_only=False)
    by = {}
    for pp in pps:
        cid = pp[0].get_parent().id
        by.setdefault(cid, []).extend(list(pp))
    out = []
    for chain in model:
        res = by.get(chain.id, [])
        if len(res) < 3:
            continue
        out.append({"chain": chain.id, "residues": res, "sequence": "".join(one(r) for r in res),
                    "first": int(res[0].id[1]), "last": int(res[-1].id[1]),
                    "fragments": sum(1 for pp in pps if pp[0].get_parent().id == chain.id)})
    return out


# ---------------------------------------------------------------- sequence mapping
def aligner():
    """Global alignment, BLOSUM62, gap open −10 / extend −0.5; end gaps free, since a crystallised construct is
    usually a fragment of the UniProt sequence (and may carry a tag)."""
    from Bio.Align import PairwiseAligner, substitution_matrices
    a = PairwiseAligner()
    a.mode = "global"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -10
    a.extend_gap_score = -0.5
    try:                                  # Biopython ≥ 1.86 names
        a.end_insertion_score = 0.0
        a.end_deletion_score = 0.0
    except AttributeError:
        a.target_end_gap_score = 0.0
        a.query_end_gap_score = 0.0
    return a


def _clean(s: str) -> str:
    ok = set("ARNDCQEGHILKMFPSTWYVBZX")
    return "".join(c if c in ok else "X" for c in s.upper())


MAX_PLACEHOLDER = 300     # longest run of unresolved residues (from the author numbering) filled in before aligning


def map_sequence(query: str, ref: str, window: int = 15, min_local: float = 0.6, nums: list | None = None) -> dict:
    """Map every position of `query` (structure chain) onto `ref` (UniProt). Returns 1-based ref positions
    (None where the chain residue has no counterpart), identity over mapped pairs, and ref coverage.

    `nums` (the chain's author residue numbers) fills unresolved stretches with placeholders before aligning, so
    a residue next to a gap lands where its numbering says rather than on an equally good match a few residues
    away (2GS2: E722 before three missing residues is EGFR E746, not E749). Jumps of more than MAX_PLACEHOLDER or
    backwards (fusion partners numbered 1002…) get no placeholders.

    Aligned pairs are kept only where some window of `window` consecutive alignment events is at least
    `min_local` identical, counting chain residues that the alignment inserts (no UniProt counterpart) as
    mismatches — that is how a fusion partner, an expression tag or a linker smeared across the sequence by the
    global alignment is kept from being mapped onto the protein (2RH1's T4 lysozyme, 3SN6's). Then every kept run
    is trimmed at both ends until it starts with an identical pair followed by at least 4 identical of 5.
    """
    if not query or not ref:
        return {"pos": [None] * len(query), "identity": 0.0, "aligned": 0, "coverage": 0.0, "dropped": 0}
    q2, real = [], []                       # the query with placeholders; real[j] = query index, or None
    for i, c in enumerate(query):
        if nums is not None and i > 0:
            try:
                gap = int(nums[i]) - int(nums[i - 1]) - 1
            except (TypeError, ValueError):
                gap = 0
            if 1 <= gap <= MAX_PLACEHOLDER:
                q2 += ["X"] * gap
                real += [None] * gap
        q2.append(c)
        real.append(i)
    aln = aligner().align(_clean(ref), _clean("".join(q2)))[0]
    ev = []                                  # alignment events in order: (query index, ref index | None, identical)
    qprev = 0
    for (t0, t1), (q0, q1) in zip(*aln.aligned):
        for j in range(qprev, int(q0)):      # chain residues the alignment inserts: no UniProt counterpart
            if real[j] is not None:
                ev.append((real[j], None, False))
        for k in range(int(t1 - t0)):
            ti, j = int(t0) + k, int(q0) + k
            if real[j] is not None:
                qi = real[j]
                ev.append((qi, ti, ref[ti].upper() == query[qi].upper()))
        qprev = int(q1)
    for j in range(qprev, len(q2)):
        if real[j] is not None:
            ev.append((real[j], None, False))
    match = np.array([m for _, _, m in ev], dtype=float)
    n_e = len(ev)
    keep = np.array([t is not None for _, t, _ in ev], bool)
    if n_e >= window:
        cs = np.concatenate([[0.0], np.cumsum(match)])
        good = (cs[window:] - cs[:-window]) / window >= min_local          # window starting at j is good
        win = np.zeros(n_e, bool)
        for j in np.flatnonzero(good):
            win[j:j + window] = True
        keep &= win
        # trim both ends of every kept run: an end must be an identical pair with ≥ 4 of the 5 events next to it
        # (inwards, itself included) identical, spaced in UniProt as they are in the chain's numbering — a lone
        # chance match in a fusion partner or tag (1IFR's GSH: H434 is not LMNA H433) is not an anchor
        num = [int(nums[qi]) if nums is not None and isinstance(nums[qi], (int, np.integer)) else qi for qi, _, _ in ev]

        def anchor(i, step, lo, hi):
            ks = [k for k in range(i, i + 5 * step, step) if lo <= k <= hi]
            if ev[i][1] is None:
                return False
            if nums is not None and all(ev[k][1] is not None and num[k] == ev[k][1] + 1 for k in ks):
                return True        # numbered exactly as UniProt: the end belongs even where it differs (4OBE's KRAS4B tail)
            if not match[i] or sum(match[k] for k in ks) < min(4, len(ks)):
                return False
            pr = sorted(k for k in ks if ev[k][1] is not None)
            return all(ev[b][1] - ev[a][1] == num[b] - num[a] for a, b in zip(pr, pr[1:]))
        i = 0
        while i < n_e:
            if not win[i]:
                i += 1
                continue
            j = i
            while j < n_e and win[j]:
                j += 1
            a, b = i, j - 1
            while a <= b and not anchor(a, 1, a, b):
                keep[a] = False
                a += 1
            while b >= a and not anchor(b, -1, a, b):
                keep[b] = False
                b -= 1
            i = j
    pos = [None] * len(query)
    same = n = 0
    for (qi, ti, m), k in zip(ev, keep):
        if k and ti is not None:
            pos[qi] = ti + 1
            n += 1
            same += m
    n_pairs = sum(1 for _, t, _ in ev if t is not None)
    return {"pos": pos, "identity": round(same / n, 4) if n else 0.0, "aligned": n,
            "coverage": round(n / len(ref), 4), "score": float(aln.score), "dropped": int(n_pairs - n)}


# ---------------------------------------------------------------- building sub-models
def is_ligand(r) -> bool:
    """A HETATM group that is not an amino acid (MSE and other modified residues stay part of the chain)."""
    return r.id[0].startswith("H_") and r.get_resname() not in THREE and "CA" not in r


def _new_model(chains, keep_het=True):
    """A detached Structure/Model holding copies of these chains, without hydrogens or water
    (and without ligands when keep_het is False)."""
    from Bio.PDB.Model import Model
    from Bio.PDB.Structure import Structure
    s = Structure("sub")
    m = Model(0)
    s.add(m)
    for ch in chains:
        c = ch.copy()
        for r in list(c):
            if r.get_resname() in WATER or (not keep_het and is_ligand(r)):
                c.detach_child(r.id)
                continue
            for a in list(r):
                if (a.element or "").upper() in ("H", "D"):
                    r.detach_child(a.id)
            if len(r) == 0:
                c.detach_child(r.id)
        if len(c):
            m.add(c)
    return s, m


def near_chains(model, chain_id: str, cutoff: float = 8.0) -> list[str]:
    """Other chains with any heavy atom within `cutoff` Å of this chain (the only ones that can bury its surface)."""
    from Bio.PDB import NeighborSearch
    atoms = [a for a in model.get_atoms() if (a.element or "").upper() not in ("H", "D")
             and a.get_parent().get_resname() not in WATER]
    ns = NeighborSearch(atoms)
    near = set()
    for a in model[chain_id].get_atoms():
        if (a.element or "").upper() in ("H", "D"):
            continue
        for b in ns.search(a.coord, cutoff):
            c = b.get_parent().get_parent().id
            if c != chain_id:
                near.add(c)
    return sorted(near)


def sasa(model, chain_id: str) -> dict:
    """Per-residue SASA (Å², Shrake–Rupley) for one chain: alone, and in the context of the whole model.

    'Whole model' means the chain plus every chain/ligand within 8 Å of it — farther ones cannot touch
    its surface, and leaving them out keeps big complexes fast. Waters and hydrogens are ignored.
    """
    from Bio.PDB.SASA import ShrakeRupley
    sr = ShrakeRupley()
    chain = model[chain_id]
    _, alone = _new_model([chain], keep_het=False)
    sr.compute(alone, level="R")
    out = {rkey(r): {"alone": float(r.sasa)} for r in alone[chain_id]}
    partners = near_chains(model, chain_id)
    ligands = sorted({r.get_resname() for r in chain if is_ligand(r)})
    if partners or ligands:
        _, cplx = _new_model([chain] + [model[c] for c in partners], keep_het=True)
        sr.compute(cplx, level="R")
        for r in cplx[chain_id]:
            out.setdefault(rkey(r), {})["complex"] = float(r.sasa)
    else:
        for k in out:
            out[k]["complex"] = out[k]["alone"]
    return {"per_residue": out, "partners": partners, "ligands": ligands}


def rel(sasa_value, resname) -> float | None:
    m = max_acc().get(resname)
    if not m or sasa_value is None:
        return None
    return round(float(sasa_value) / m, 3)


def rsa_class(v) -> str:
    if v is None:
        return "—"
    return "buried" if v < 0.20 else "exposed" if v > 0.50 else "intermediate"


# ---------------------------------------------------------------- secondary structure
def dssp_available() -> str | None:
    return shutil.which("mkdssp") or shutil.which("dssp")


def secondary_structure(model, chain_id: str, path: Path | None = None) -> dict:
    """{'method': ..., 'label': {...}, 'ss': {key: 'H'|'E'|'C'}}.

    With mkdssp on PATH, Bio.PDB.DSSP (Kabsch & Sander 1983). Without it — the usual case, since the binary
    is not a Python package — a phi/psi approximation: each residue's backbone dihedrals are placed in a
    Ramachandran region, and a run of ≥ 4 helical residues is called helix-like, ≥ 3 extended residues
    extended. It cannot see hydrogen bonds, so 'extended' also covers polyproline and disordered stretches,
    not only paired β-strands.
    """
    exe = dssp_available()
    if exe and path is not None:
        try:
            from Bio.PDB.DSSP import DSSP
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                d = DSSP(model, str(path), dssp=exe)
            ss = {}
            for (cid, rid) in d.keys():
                if cid != chain_id:
                    continue
                code = d[(cid, rid)][2]
                ss[f"{cid}:{rid[1]}{rid[2].strip()}"] = "H" if code in "HGI" else "E" if code in "EB" else "C"
            if ss:
                return {"method": "dssp", "label": "DSSP (Kabsch & Sander 1983)", "ss": ss}
        except Exception:  # noqa: BLE001 - fall through to the approximation
            from .common import log_exc
            log_exc("dssp")
    from Bio.PDB import PPBuilder
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pps = [pp for pp in PPBuilder().build_peptides(model[chain_id], aa_only=False)]
    ss = {}
    for pp in pps:
        raw = []
        for res, (phi, psi) in zip(pp, pp.get_phi_psi_list()):
            if phi is None or psi is None:
                raw.append("C")
                continue
            f, s = np.degrees(phi), np.degrees(psi)
            if -160 <= f <= -20 and -120 <= s <= 50:
                raw.append("H")
            elif -180 <= f <= -45 and (s >= 90 or s <= -150):
                raw.append("E")
            else:
                raw.append("C")
        lab = _smooth(raw)
        for res, c in zip(pp, lab):
            ss[rkey(res)] = c
    return {"method": "phi/psi", "label": "approximation from backbone phi/psi (no DSSP binary installed)", "ss": ss}


def _smooth(raw: list[str]) -> list[str]:
    out = ["C"] * len(raw)
    i = 0
    while i < len(raw):
        j = i
        while j < len(raw) and raw[j] == raw[i]:
            j += 1
        if (raw[i] == "H" and j - i >= 4) or (raw[i] == "E" and j - i >= 3):
            out[i:j] = [raw[i]] * (j - i)
        i = j
    return out


SS_WORD = {"dssp": {"H": "helix", "E": "strand", "C": "coil/turn"},
           "phi/psi": {"H": "helix-like", "E": "extended", "C": "other"}}


# ---------------------------------------------------------------- per-residue table
def residue_table(model, chain_id: str, residues: list, ref_pos: list, path: Path | None = None,
                  bkind: str = "plddt") -> dict:
    """One row per chain residue: reference (UniProt) position, CA B-factor/pLDDT, SASA (alone and in the
    model), relative SASA, secondary structure. Also returns which chains/ligands are nearby."""
    sa = sasa(model, chain_id)
    ss = secondary_structure(model, chain_id, path)
    rows = []
    for res, p in zip(residues, ref_pos):
        k = rkey(res)
        s = sa["per_residue"].get(k, {})
        ca = res["CA"] if "CA" in res else None
        b = float(ca.get_bfactor()) if ca is not None else None
        rows.append({"key": k, "num": int(res.id[1]), "icode": res.id[2].strip(), "aa": one(res), "resname": res.get_resname(),
                     "pos": p, "b": None if b is None else round(b, 2),
                     "sasa": round(s["complex"], 1) if "complex" in s else None,
                     "sasa_alone": round(s["alone"], 1) if "alone" in s else None,
                     "rsa": rel(s.get("complex"), res.get_resname()), "rsa_alone": rel(s.get("alone"), res.get_resname()),
                     "ss": ss["ss"].get(k, "C")})
    return {"rows": rows, "partners": sa["partners"], "ss_method": ss["method"], "ss_label": ss["label"], "bkind": bkind}


# ---------------------------------------------------------------- neighbours
def find_residue(model, key: str):
    cid, rest = key.split(":", 1)
    num = int("".join(ch for ch in rest if ch.isdigit() or ch == "-"))
    ic = "".join(ch for ch in rest if ch.isalpha()) or " "
    chain = model[cid]
    for r in chain:
        if r.id[1] == num and r.id[2] == ic:
            return r
    raise KeyError(key)


def neighbours(model, key: str, radius: float = 5.0, names: dict | None = None) -> dict:
    """Residues with any heavy atom within `radius` Å of any heavy atom of the residue (NeighborSearch)."""
    from Bio.PDB import NeighborSearch
    target = find_residue(model, key)
    heavy = [a for a in model.get_atoms() if (a.element or "").upper() not in ("H", "D")]
    ns = NeighborSearch(heavy)
    found = {}
    for a in target:
        if (a.element or "").upper() in ("H", "D"):
            continue
        for b in ns.search(a.coord, radius):
            r = b.get_parent()
            if r is target:
                continue
            d = float(np.linalg.norm(a.coord - b.coord))
            k = rkey(r)
            if k not in found or d < found[k][1]:
                found[k] = (r, d)
    tc = target.get_parent().id
    out = {"same_chain": [], "other_chain": [], "ligand": [], "water": 0}
    for k, (r, d) in sorted(found.items(), key=lambda kv: kv[1][1]):
        c = r.get_parent().id
        rn = r.get_resname()
        if rn in WATER:
            out["water"] += 1
            continue
        row = {"key": k, "resname": rn, "num": int(r.id[1]), "chain": c, "min_dist": round(d, 2)}
        if is_ligand(r):
            out["ligand"].append(row)
        elif c != tc:
            row["molecule"] = (names or {}).get(c, "")
            out["other_chain"].append(row)
        else:
            row["seq_sep"] = int(abs(int(r.id[1]) - int(target.id[1])))
            out["same_chain"].append(row)
    out["n_total"] = len(out["same_chain"]) + len(out["other_chain"]) + len(out["ligand"])
    out["n_nonlocal"] = sum(1 for x in out["same_chain"] if x["seq_sep"] > 2) + len(out["other_chain"]) + len(out["ligand"])
    out["target"] = {"key": key, "resname": target.get_resname()}
    return out


# ---------------------------------------------------------------- comparison
def ce_rmsd(model_a, chain_a: str, model_b, chain_b: str) -> dict:
    """Bio.PDB.cealign.CEAligner on the CA atoms of two single chains. Returns RMSD and the aligned length."""
    from Bio.PDB.cealign import CEAligner
    _, ma = _new_model([model_a[chain_a]], keep_het=False)
    _, mb = _new_model([model_b[chain_b]], keep_het=False)
    ce = CEAligner()
    ce.set_reference(ma)
    ce.align(mb, transform=False)
    n = None
    try:
        n = len(ce._coord)          # guide atoms in the mobile chain (for reporting only)
    except Exception:  # noqa: BLE001
        pass
    return {"rmsd": round(float(ce.rms), 3), "n_guide": n}


def all_atoms(structure):
    """Every atom object, including every alternate location of a disordered atom."""
    for a in structure.get_atoms():
        if a.is_disordered():
            yield from a.disordered_get_list()
        else:
            yield a


def _kabsch(X: np.ndarray, Y: np.ndarray):
    """(rot, tran) that best puts Y on X, in Bio.PDB's convention: Y @ rot + tran."""
    from Bio.SVDSuperimposer import SVDSuperimposer
    sup = SVDSuperimposer()
    sup.set(X, Y)
    sup.run()
    return sup.get_rotran()


def superpose_pairs(fixed_ca: list, moving_ca: list, moving_structure, cutoff: float = 3.0, cycles: int = 5):
    """Superpose matched CA pairs on their largest common rigid core, then move the whole moving structure.

    A least-squares fit over every pair is dragged by any part that moved as a block — a flexible tail, a hinged
    domain — and then shows the whole chain as 'deviating' (EGFR kinase, AlphaFold vs 1M17: well under 1 Å over the
    kinase core but 12 Å over all pairs, because a C-terminal tail sits elsewhere). So the fit is seeded on every
    stretch of 20 consecutive pairs, the seed that brings the most pairs within `cutoff` Å wins, and the fit is
    refined on the pairs within `cutoff` until they stop changing (≤ `cycles` times) — the idea behind LGA/TM-align
    style superposition. Falls back to all pairs when no core of ≥ 3 pairs exists. Returns (rms over all pairs,
    (rot, tran), core mask, rms over the core).

    The moving structure must be a freshly parsed one, not Entity.copy(): a copied DisorderedAtom keeps pointing
    at the original's selected altloc, so transforming the copy silently leaves those atoms where they were.
    """
    X = np.array([a.coord for a in fixed_ca], dtype=float)
    Y = np.array([a.coord for a in moving_ca], dtype=float)
    n = len(X)
    dev = lambda rt: np.linalg.norm(X - (Y @ rt[0] + rt[1]), axis=1)  # noqa: E731
    best = None
    L = min(n, 20)
    for s0 in range(0, n - L + 1, max(1, L // 4)):
        rt = _kabsch(X[s0:s0 + L], Y[s0:s0 + L])
        k = int((dev(rt) <= cutoff).sum())
        if best is None or k > best[0]:
            best = (k, rt)
    core = np.ones(n, bool)
    if best and best[0] >= 3:
        core = dev(best[1]) <= cutoff
        for _ in range(cycles):
            new = dev(_kabsch(X[core], Y[core])) <= cutoff
            if new.sum() < 3 or (new == core).all():
                break
            core = new
    rot, tran = _kabsch(X[core], Y[core])
    d = dev((rot, tran))
    for a in all_atoms(moving_structure):
        a.transform(rot, tran)
    return float(np.sqrt(np.mean(d ** 2))), (rot, tran), core, float(np.sqrt(np.mean(d[core] ** 2)))


def replace_bfactor(structure, values: dict, default: float = 0.0):
    """Set every atom's B-factor to values[residue key] (or default). Used for 'colour by value' exports."""
    for r in structure[0].get_residues():
        v = values.get(rkey(r), default)
        for a in r.get_atoms():
            for x in (a.disordered_get_list() if a.is_disordered() else [a]):
                x.set_bfactor(float(v))


def run_ok(cmd) -> bool:
    try:
        subprocess.run(cmd, capture_output=True, timeout=5)
        return True
    except Exception:  # noqa: BLE001
        return False
