"""UniProt: gene → accession search (Bio.UniProt.search) and the full entry (REST JSON), parsed into plain dicts."""
from __future__ import annotations

import re

from . import net

SPECIES = {
    "human": (9606, "Homo sapiens"), "mouse": (10090, "Mus musculus"), "rat": (10116, "Rattus norvegicus"),
    "zebrafish": (7955, "Danio rerio"), "fly": (7227, "Drosophila melanogaster"),
    "worm": (6239, "Caenorhabditis elegans"), "yeast": (559292, "Saccharomyces cerevisiae S288C"),
    "ecoli": (83333, "Escherichia coli K-12"),
}
ALIASES = {"homo sapiens": "human", "hs": "human", "mus musculus": "mouse", "mm": "mouse", "rattus norvegicus": "rat",
           "danio rerio": "zebrafish", "drosophila": "fly", "drosophila melanogaster": "fly", "d. melanogaster": "fly",
           "c. elegans": "worm", "caenorhabditis elegans": "worm", "elegans": "worm", "s. cerevisiae": "yeast",
           "saccharomyces cerevisiae": "yeast", "budding yeast": "yeast", "e. coli": "ecoli", "e.coli": "ecoli",
           "escherichia coli": "ecoli", "coli": "ecoli"}
ACC_RE = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-(\d+))?$")


def taxon(species) -> tuple[int, str]:
    """'human' / 'Mouse' / '9606' / 10090 → (taxon id, label). A number not in the list is accepted as-is."""
    s = str(species if species not in (None, "") else "human").strip().lower()
    s = ALIASES.get(s, s)
    if s in SPECIES:
        return SPECIES[s]
    if s.isdigit():
        tid = int(s)
        for t, lab in SPECIES.values():
            if t == tid:
                return t, lab
        return tid, f"taxon {tid}"
    raise ValueError(f"'{species}' is not a species this app knows. Use human, mouse, rat, zebrafish, fly, worm, "
                     "yeast, ecoli, or type the NCBI taxon ID (for example 9606).")


def is_accession(s: str) -> bool:
    return bool(ACC_RE.match((s or "").strip().upper()))


FIELDS = ["accession", "id", "protein_name", "gene_names", "organism_name", "organism_id", "length", "reviewed"]


def _hit(r: dict) -> dict:
    pd = r.get("proteinDescription", {})
    name = (pd.get("recommendedName") or (pd.get("submissionNames") or [{}])[0]).get("fullName", {}).get("value", "")
    genes = r.get("genes") or []
    return {"accession": r.get("primaryAccession"), "id": r.get("uniProtkbId"), "name": name,
            "gene": (genes[0].get("geneName") or {}).get("value", "") if genes else "",
            "organism": r.get("organism", {}).get("scientificName", ""), "taxon": r.get("organism", {}).get("taxonId"),
            "length": r.get("sequence", {}).get("length"),
            "reviewed": "reviewed" in (r.get("entryType") or "").lower() and "unreviewed" not in (r.get("entryType") or "").lower()}


def search_gene(gene: str, species="human", refresh=False) -> dict:
    """Reviewed (Swiss-Prot) entries whose gene name is exactly `gene` in that organism.

    Falls back to unreviewed (TrEMBL) entries when there is no reviewed one, and says so.
    """
    from Bio import UniProt
    gene = (gene or "").strip()
    if not re.match(r"^[A-Za-z0-9][A-Za-z0-9_.\-/]{0,39}$", gene):
        raise ValueError(f"'{gene}' does not look like a gene symbol (letters, digits, - . _ only).")
    tid, lab = taxon(species)
    out = {"gene": gene, "taxon": tid, "species": lab, "reviewed": True, "hits": [], "source": None}
    for reviewed in ("true", "false"):
        q = f"gene_exact:{gene} AND organism_id:{tid} AND reviewed:{reviewed}"

        def run(q=q):
            res = UniProt.search(q, fields=FIELDS, batch_size=25)
            hits = []
            for i, r in enumerate(res):
                if i >= 25:
                    break
                hits.append(_hit(r))
            return hits

        hits, src = net.cached_call(f"uniprot/search/{net.key(q)}.json", run, f"the UniProt search for {gene} ({lab})", refresh)
        out["source"] = src
        if hits:
            out["hits"] = hits
            out["reviewed"] = reviewed == "true"
            break
    return out


def fetch_entry(acc: str, refresh=False) -> tuple[dict, str]:
    acc = acc.strip().upper().split("-")[0]
    return net.fetch_json(f"https://rest.uniprot.org/uniprotkb/{acc}.json", f"uniprot/{acc}.json",
                          f"the UniProt entry for {acc}", refresh)


# ---------------------------------------------------------------- parsing
FEATURE_GROUPS = {
    "Domains & regions": ["Domain", "Region", "Coiled coil", "Motif", "Zinc finger", "Repeat", "DNA binding",
                          "Transmembrane", "Intramembrane", "Topological domain", "Signal", "Transit peptide",
                          "Propeptide", "Chain", "Peptide", "Compositional bias"],
    "Sites": ["Active site", "Binding site", "Site"],
    "PTMs": ["Modified residue", "Lipidation", "Glycosylation", "Disulfide bond", "Cross-link"],
    "Natural variants": ["Natural variant"],
}
TRACK_TYPES = ["Domain", "Region", "Coiled coil", "Motif", "Zinc finger", "Repeat", "DNA binding", "Transmembrane",
               "Intramembrane", "Signal", "Transit peptide", "Propeptide"]


def _union_len(segs: list[dict]) -> int:
    """Residues covered by the union of a PDB entry's (start, end) ranges. UniProt's 'Chains' property lists each
    chain's range separately when they are not identical enough to join with '/' (LMNA 6YSH: 'A=25-70, B=26-70',
    both chains covering the same 46-residue stretch) — summing the per-chain ranges would count that stretch
    twice, making a homodimer resolved once look like it covers twice as much sequence as it does."""
    if not segs:
        return 0
    merged = []
    for s, e in sorted({(x["start"], x["end"]) for x in segs}):
        if merged and s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return sum(e - s + 1 for s, e in merged)


def _loc(f):
    s = f.get("location", {}).get("start", {}).get("value")
    e = f.get("location", {}).get("end", {}).get("value")
    return s, e


def _txt(c):
    return " ".join(t.get("value", "") for t in (c.get("texts") or []))


def parse_entry(d: dict) -> dict:
    acc = d.get("primaryAccession")
    pd = d.get("proteinDescription", {})
    name = (pd.get("recommendedName") or (pd.get("submissionNames") or [{}])[0]).get("fullName", {}).get("value", "")
    genes = d.get("genes") or []
    gene = (genes[0].get("geneName") or {}).get("value", "") if genes else ""
    org = d.get("organism", {})
    seq = d.get("sequence", {}).get("value", "")
    comments = d.get("comments") or []
    function = " ".join(_txt(c) for c in comments if c.get("commentType") == "FUNCTION" and not c.get("molecule")) or \
        " ".join(_txt(c) for c in comments if c.get("commentType") == "FUNCTION")
    locs = []
    for c in comments:
        if c.get("commentType") == "SUBCELLULAR LOCATION":
            for s in c.get("subcellularLocations") or []:
                v = s.get("location", {}).get("value")
                if v and v not in locs:
                    locs.append(v + (f" ({c['molecule']})" if c.get("molecule") else ""))
    diseases = []
    for c in comments:
        if c.get("commentType") != "DISEASE":
            continue
        dz = c.get("disease") or {}
        if dz:
            diseases.append({"name": dz.get("diseaseId", ""), "acronym": dz.get("acronym", ""),
                             "description": dz.get("description", ""),
                             "mim": (dz.get("diseaseCrossReference") or {}).get("id", ""),
                             "note": _txt(c.get("note") or {})})
        else:
            diseases.append({"name": "", "acronym": "", "description": _txt(c.get("note") or {}), "mim": "", "note": ""})
    isoforms = []
    for c in comments:
        if c.get("commentType") == "ALTERNATIVE PRODUCTS":
            for iso in c.get("isoforms") or []:
                isoforms.append({"name": (iso.get("name") or {}).get("value", ""), "ids": iso.get("isoformIds") or [],
                                 "status": iso.get("isoformSequenceStatus", "")})
    feats = []
    for f in d.get("features") or []:
        s, e = _loc(f)
        if s is None or e is None:
            continue
        item = {"type": f.get("type"), "start": int(s), "end": int(e), "description": f.get("description", "") or ""}
        alt = f.get("alternativeSequence") or {}
        if f.get("type") in ("Natural variant", "Mutagenesis") and alt:
            item["ref"] = alt.get("originalSequence", "")
            item["alt"] = alt.get("alternativeSequences") or []
        if f.get("featureId"):
            item["id"] = f["featureId"]
        feats.append(item)
    groups = {g: [f for f in feats if f["type"] in types] for g, types in FEATURE_GROUPS.items()}
    variants = [f for f in feats if f["type"] == "Natural variant"]
    signal = next((f for f in feats if f["type"] == "Signal"), None)
    pdb = []
    for x in d.get("uniProtKBCrossReferences") or []:
        if x.get("database") != "PDB":
            continue
        props = {p["key"]: p["value"] for p in x.get("properties") or []}
        res = props.get("Resolution", "-")
        try:
            resv = float(res.split()[0])
        except Exception:  # noqa: BLE001
            resv = None
        segs = []
        for part in (props.get("Chains") or "").split(","):
            part = part.strip()
            m = re.match(r"^([\w/]+)=(-?\d+)-(-?\d+)$", part)
            if m:
                segs.append({"chains": m.group(1).split("/"), "start": int(m.group(2)), "end": int(m.group(3))})
        cov = _union_len(segs)
        pdb.append({"id": x["id"], "method": props.get("Method", ""), "resolution": resv, "resolution_text": res,
                    "segments": segs, "span": [min(s["start"] for s in segs), max(s["end"] for s in segs)] if segs else None,
                    "covered": cov})
    return {"accession": acc, "id": d.get("uniProtkbId"), "name": name, "gene": gene,
            "organism": org.get("scientificName", ""), "common": org.get("commonName", ""), "taxon": org.get("taxonId"),
            "reviewed": "unreviewed" not in (d.get("entryType") or "").lower(),
            "length": len(seq), "sequence": seq, "function": function, "location": locs, "diseases": diseases,
            "isoforms": isoforms, "features": feats, "groups": groups, "variants": variants,
            "signal_length": signal["end"] if signal and signal["start"] == 1 else 0,
            "pdb": pdb, "version": (d.get("entryAudit") or {}).get("entryVersion"),
            "last_updated": (d.get("entryAudit") or {}).get("lastAnnotationUpdateDate")}


def track_features(u: dict) -> list[dict]:
    """Features worth a row on the sequence track: domains first, then the rest by type."""
    order = {t: i for i, t in enumerate(TRACK_TYPES)}
    fs = [f for f in u["features"] if f["type"] in order]
    return sorted(fs, key=lambda f: (order[f["type"]], f["start"]))


def domain_at(u: dict | None, pos: int) -> list[str]:
    """Names of the UniProt domains/regions covering a position, most specific (shortest) first."""
    if not u:
        return []
    hits = [f for f in u["features"] if f["type"] in ("Domain", "Region", "Coiled coil", "Motif", "Zinc finger", "Repeat",
                                                     "DNA binding", "Transmembrane", "Signal", "Propeptide", "Topological domain")
            and f["start"] <= pos <= f["end"]]
    rank = {"Domain": 0, "Region": 2, "Propeptide": 3, "Topological domain": 3}
    hits.sort(key=lambda f: (rank.get(f["type"], 1), f["end"] - f["start"]))   # domains, then specific features, then regions
    return [f"{f['description'] or f['type']}" + ("" if f["type"] in ("Domain", "Region") else f" ({f['type'].lower()})")
            + f" {f['start']}–{f['end']}" for f in hits]


def main_domains(u: dict | None) -> list[dict]:
    """Blocks for PAE summaries: UniProt 'Domain' features (non-overlapping, ≥ 20 residues)."""
    if not u:
        return []
    ds = sorted([f for f in u["features"] if f["type"] == "Domain" and f["end"] - f["start"] >= 19], key=lambda f: f["start"])
    out = []
    for f in ds:
        if out and f["start"] <= out[-1]["end"]:
            continue
        out.append({"name": f["description"] or "Domain", "start": f["start"], "end": f["end"]})
    return out
