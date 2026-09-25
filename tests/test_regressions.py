"""Regression tests for the bugs found in the 0.1.1 stress test (lookups, mapping, comparison, uploads, reports).
Offline, on the bundled LMNA example and small synthetic inputs.

    .venv/bin/python -m pytest tests -q
"""
import gzip
import json
import threading

import numpy as np
import pytest
from conftest import MODEL, UNIPROT

from server import agent, core, net, report
from server import alphafold as afdb
from server import structure as st
from server import uniprot as up
from server import variants as vr
from server.common import WORK, UserFacingError


def _seq():
    return up.parse_entry(json.loads(UNIPROT.read_text()))["sequence"]


def _item():
    return core.lookup(accession="P02545")["item"]


# ---------------------------------------------------------------- lookup
def test_primary_gene_name_beats_a_synonym(monkeypatch):
    """ALB is also a synonym of FBF1: the entry whose primary gene name is the query must win, not the first hit."""
    hits = [{"accession": "Q8TES7", "gene": "FBF1", "name": "Fas-binding factor 1", "reviewed": True},
            {"accession": "P02545", "gene": "LMNA", "name": "Prelamin-A/C", "reviewed": True}]
    monkeypatch.setattr(up, "search_gene", lambda g, s, r=False: {"hits": hits, "reviewed": True, "species": "Homo sapiens", "taxon": 9606})
    r = core.lookup(gene="LMNA", species="human", reuse=False)
    res = core.load_result(r["item"])
    assert res["accession"] == "P02545"
    assert any("synonym of FBF1" in f["text"] for f in res["flags"])


def test_several_primary_matches_still_ask(monkeypatch):
    hits = [{"accession": "P42771", "gene": "CDKN2A"}, {"accession": "Q8N726", "gene": "CDKN2A"}, {"accession": "X1", "gene": "OTHER"}]
    monkeypatch.setattr(up, "search_gene", lambda g, s, r=False: {"hits": hits, "reviewed": True, "species": "Homo sapiens", "taxon": 9606})
    r = core.lookup(gene="CDKN2A")
    assert [h["accession"] for h in r["choices"]] == ["P42771", "Q8N726", "X1"]      # synonym-only match listed last
    assert "only as a synonym" in r["message"]


def test_isoform_accession_says_canonical_is_shown():
    r = core.lookup(accession="P02545-2")
    assert core.load_result(r["item"])["accession"] == "P02545"
    assert "isoform" in r["notice"] and "canonical" in r["notice"]


def test_opaque_alphafold_entry_id_is_the_canonical_model():
    """SARS-CoV-2 spike's only AFDB entry is AF-0000000365840314 (ColabFold), keyed to P0DTC2 from residue 1."""
    preds = [{"entryId": "AF-0000000365840314", "uniprotAccession": "P0DTC2", "uniprotStart": 1, "isComplex": False,
              "cifUrl": "x", "sequence": "MFVF"}]
    canon, iso = afdb.pick(preds, "P0DTC2")
    assert canon is preds[0] and iso == []
    canon, iso = afdb.pick([{"entryId": "AF-P02545-2-F1", "uniprotAccession": "P02545-2", "uniprotStart": 1}], "P02545")
    assert canon is None and len(iso) == 1                      # an isoform model is never taken as the canonical one


def test_concurrent_identical_lookups_make_one_item():
    for d in list(WORK.iterdir()):
        try:
            if json.loads((d / "result.json").read_text()).get("accession") == "P02545":
                import shutil
                shutil.rmtree(d)
        except Exception:  # noqa: BLE001
            pass
    out = []
    ts = [threading.Thread(target=lambda: out.append(core.lookup(accession="P02545")["item"])) for _ in range(6)]
    [t.start() for t in ts]
    [t.join() for t in ts]
    assert len(out) == 6 and len(set(out)) == 1


def test_refresh_that_cannot_reach_the_internet_says_so():
    net.fetch("https://example.invalid/x", "alphafold/AF-P02545-F1-model_v6.cif", "the model")      # seeds the cache
    data, src = net.fetch("https://example.invalid/x", "alphafold/AF-P02545-F1-model_v6.cif", "the model", refresh=True)
    assert src == "stale" and data
    r = core.lookup(accession="P02545", refresh=True)
    flags = " ".join(f["text"] for f in core.load_result(r["item"])["flags"])
    assert "Refresh was asked for" in flags and "Loaded from the local cache" not in flags


# ---------------------------------------------------------------- variants
def test_position_zero_is_not_beyond_the_end():
    chk = vr.check_reference(vr.parse("A0G"), _seq())
    assert "numbering starts at 1" in chk["message"] and "beyond" not in chk["message"]


# ---------------------------------------------------------------- chain → UniProt mapping
def test_mapping_drops_a_fusion_partner_inside_the_chain():
    """2RH1/3SN6-style construct: receptor, then a T4-lysozyme-like insert numbered 1001…, then receptor again."""
    ref = _seq()
    junk = "NIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMG"
    q = ref[99:230] + junk + ref[262:360]
    nums = list(range(100, 231)) + list(range(1001, 1001 + len(junk))) + list(range(263, 361))
    m = st.map_sequence(q, ref, nums=nums)
    a, b = 131, 131 + len(junk)
    assert all(p is None for p in m["pos"][a:b]), [p for p in m["pos"][a:b] if p]
    assert m["pos"][:a] == list(range(100, 231)) and m["pos"][b:] == list(range(263, 361))


def test_mapping_places_a_residue_before_a_gap_by_its_number():
    """2GS2: 'IKE' + 3 missing + 'ATS' — E could align to E746 or E749; the numbering says 746."""
    ref = "MSTQLAEGKIPVAIKELREATSPKANKEILDEAYVMASVDNPHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDN"
    i = ref.index("IKELREATS")
    q = ref[:i + 3] + ref[i + 6:]                      # 'LRE' unresolved
    nums = list(range(1, i + 4)) + list(range(i + 7, len(ref) + 1))
    m = st.map_sequence(q, ref, nums=nums)
    assert m["pos"][i + 2] == i + 3                    # E at its own number, not three further on
    assert m["pos"] == nums


# ---------------------------------------------------------------- comparison
def test_core_superposition_is_not_dragged_by_a_moved_tail():
    s1, _ = st.load(MODEL, "a")
    s2, _ = st.load(MODEL, "b")
    for r in s2[0]["A"]:
        if 360 <= r.id[1] <= 387:                      # swing the end of the rod away by 15 Å
            for a in r:
                a.coord = a.coord + np.array([15.0, 0, 0], dtype=np.float32)
    keys = range(300, 388)
    fixed = [s1[0]["A"][k]["CA"] for k in keys]
    moving = [s2[0]["A"][k]["CA"] for k in keys]
    rms, _, core_mask, rms_core = st.superpose_pairs(fixed, moving, s2)
    assert rms_core < 0.01 and core_mask.sum() == 60 and rms > 5            # the 60 unmoved pairs are the core
    assert s2[0]["A"][320]["CA"].coord == pytest.approx(s1[0]["A"][320]["CA"].coord, abs=1e-3)


def test_compare_rejects_bad_or_foreign_chains():
    item = _item()
    with pytest.raises(UserFacingError, match="not a protein chain"):
        core.compare_structures(item, "AF", "AF", "Q", None)
    with pytest.raises(UserFacingError, match="Choose two structures"):
        core.compare_structures(item, "", "")


def test_residue_environment_bad_input_is_a_message():
    item = _item()
    for res, rad in (("A:xyz", 5), ("482", "big")):
        with pytest.raises(UserFacingError):
            core.residue_environment(item, res, rad)


# ---------------------------------------------------------------- uploads
def test_gzipped_upload_gets_a_readable_viewer_file_and_its_own_id():
    item = _item()
    data = gzip.compress(MODEL.read_bytes())
    a = core.upload("model.cif", MODEL.read_bytes(), item)
    b = core.upload("model.cif.gz", data, item)
    assert a["sid"] != b["sid"]
    d = WORK / item
    assert (d / a["structure"]["file"]).read_bytes() == MODEL.read_bytes()          # the first upload is intact
    view = d / b["structure"]["view_file"]
    assert view.suffix == ".cif" and view.read_bytes()[:5] == b"data_"


def test_gzipped_upload_is_limited_by_its_unpacked_size(monkeypatch):
    monkeypatch.setattr(core, "MAX_STRUCTURE_BYTES", 100_000)
    with pytest.raises(UserFacingError, match="unpacks to more than"):
        core.upload("big.cif.gz", gzip.compress(MODEL.read_bytes()))


def test_file_without_atoms_is_a_message():
    with pytest.raises(UserFacingError, match="no atom records|not a PDB"):
        core.upload("h.pdb", b"HEADER    junk\n")


def test_structure_residue_differing_from_uniprot_is_flagged():
    item = _item()
    lines = []
    for line in MODEL.read_text().splitlines():
        if line.startswith("ATOM") and " ARG A 1 482 " in line:           # label and author fields of R482
            line = line.replace("ARG", "LYS")
        lines.append(line)
    r = core.upload("r482k.cif", ("\n".join(lines) + "\n").encode(), item)
    ch = next(c for c in r["structure"]["chains"] if c["chain"] == "A")
    assert ch["differences"] == ["R482K"] and "R482K" in r["structure"]["note"]
    v = core.map_variants(item, "R482W", r["sid"], save=False)
    assert v["rows"][0]["flags"][0]["level"] == "warn" and "is Lys, not Arg" in v["rows"][0]["flags"][0]["text"]


# ---------------------------------------------------------------- reports
def test_report_file_name_format():
    assert core.report_filename({"kind": "protein", "created": "2026-09-23 10:11", "uniprot": {"gene": "LMNA"},
                                 "accession": "P02545"}) == "2026-09-23_StructureBench-LMNA_report.pdf"
    assert core.report_filename({"kind": "upload", "created": "2026-01-02 00:00", "name": "my model (v2).pdb.gz"}) \
        == "2026-01-02_StructureBench-my-model-v2_report.pdf"


def test_pdf_survives_a_reading_longer_than_a_page():
    item = _item()
    res = core.load_result(item)
    res["texts"]["pdb"]["yours"] = "Ig-like 1 6–96: no entry covers it whole. " * 400      # titin-sized reading
    core.save_result(item, res)
    p = report.build_pdf(item)
    assert p.read_bytes()[:4] == b"%PDF"


# ---------------------------------------------------------------- assistant tools
def test_assistant_tool_errors_are_sentences():
    item = _item()
    out = agent.run_tool("load_pdb_entry", {"item": item})
    assert out.startswith("error:") and "pdb_id" in out and "KeyError" not in out
    out = agent.run_tool("residue_environment", {"item": item, "residue": 482})
    assert '"key": "A:482"' in out
    out = agent.run_tool("lookup_protein", {"accession": 12345})
    assert out.startswith("error:") and "AttributeError" not in out


def test_report_has_every_section_and_explains_each_figure():
    """The 0.1.1 report left out the UniProt entry, the PAE, the 3D model and the PDB coverage (Paul, 2026-09-23).
    Every section of the page must now be in it, and every figure must carry its 'What this shows' text."""
    import base64
    import io

    pymupdf = pytest.importorskip("pymupdf")          # dev-only: reads the PDF back
    from PIL import Image

    r = core.lookup(accession="P02545")
    core.map_variants(r["item"], "E145K\nR482W", "AF")
    core.compare_structures(r["item"], "AF", "AF")
    # without the page: server-drawn CA traces stand in for the 3D pictures
    text = "".join(pg.get_text() for pg in pymupdf.open(str(report.build_pdf(r["item"]))))
    flat = " ".join(text.split())
    for must in ("What your data says", "Function", "Disease involvement", "Annotated features", "Confidence (pLDDT)",
                 "Predicted aligned error", "The AlphaFold 3D model", "Structures in this analysis",
                 "PDB entries along the sequence", "Where the variants sit", "Variants on the 3D structure",
                 "Each variant in detail", "CA deviation along the sequence", "The superposition in 3D",
                 "Methods and references", "Drawn by the server"):
        assert must in flat, must
    assert flat.count("What this shows.") >= 10 and flat.count("In your data.") >= 10
    # with the page's snapshots: they are used (and trimmed), and the fallback note disappears
    buf = io.BytesIO()
    im = Image.new("RGB", (900, 500), "white")
    im.paste((0, 83, 214), (300, 200, 600, 300))
    im.save(buf, "PNG")
    url = "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()
    text = "".join(pg.get_text() for pg in pymupdf.open(str(report.build_pdf(r["item"], {"model_front": url, "model_back": url,
                                                                                       "variants": url}))))
    assert "snapshot from the 3D viewer" in text and "cartoon" in text
    assert report._snapshot({"x": url}, "x") and Image.open(io.BytesIO(report._snapshot({"x": url}, "x"))).size[0] < 400


def test_report_is_reachable_by_both_routes():
    """The page POSTs the 3D snapshots to /report; the GET /report.pdf route is what it falls back to when the
    running server is older than the page it serves (Paul, 2026-09-23: rebuilding the app under a running copy left
    the two out of step and the PDF button just failed). Both must exist, and GET must need no body."""
    from server.app import app
    routes = {(r.path, m) for r in app.router.routes for m in getattr(r, "methods", ())}
    assert ("/api/items/{item}/report", "POST") in routes
    assert ("/api/items/{item}/report.pdf", "GET") in routes
    import inspect

    from server import app as appmod
    sig = inspect.signature(appmod.report_pdf)
    assert sig.parameters["snapshots"].default is None      # GET works without snapshots
