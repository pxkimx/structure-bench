"""Progerin (HGPS) mode: every number derived from the UniProt record must match the literature, and the mode must
appear only for human LMNA. Runs offline on the bundled P02545 example."""
import copy

import pytest

from server import core
from server import progerin as pg
from server.common import UserFacingError


@pytest.fixture(scope="module")
def item():
    return core.lookup(accession="P02545")["item"]


@pytest.fixture(scope="module")
def res(item):
    return core.load_result(item)


@pytest.fixture(scope="module")
def b(res):
    assert res.get("progerin"), "an LMNA lookup must carry the progerin block"
    return res["progerin"]


def test_deletion_matches_uniprot_and_the_literature(b, res):
    u = res["uniprot"]
    assert b["deletion"]["from_record"], "the deletion must come from the record's Progerin isoform feature, not from memory"
    assert (b["deletion"]["start"], b["deletion"]["end"]) == pg.EXPECTED["deletion"] == (607, 656)
    assert b["deletion"]["length"] == pg.EXPECTED["deletion_length"] == 50
    vsp = [f for f in u["features"] if f["type"] == "Alternative sequence" and f.get("id") == "VSP_053505"]
    assert vsp and (vsp[0]["start"], vsp[0]["end"]) == (607, 656)
    iso = next(i for i in u["isoforms"] if "Progerin" in i.get("synonyms", []))
    assert "VSP_053505" in iso["vsp"] and iso["ids"] == ["P02545-6"]
    assert b["deletion"]["sequence"] == u["sequence"][606:656]


def test_progerin_length_and_c_terminus(b, res):
    assert b["prelamin_length"] == pg.EXPECTED["prelamin_length"] == 664
    assert b["progerin"]["length"] == pg.EXPECTED["progerin_length"] == 614
    seq = res["uniprot"]["sequence"]
    assert b["progerin"]["sequence"] == seq[:606] + seq[656:]
    assert b["progerin"]["c_terminus"] == b["caax"]["sequence"] == pg.EXPECTED["caax"] == "CSIM"
    assert b["caax"]["farnesyl"] == pg.EXPECTED["farnesyl_cys"] == 661 and b["caax"]["kept_in_progerin"]
    assert b["progerin"]["farnesyl_cys"] == 611 and b["progerin"]["sequence"][610] == "C"
    assert b["caax"]["methyl_ester"] == 661 and b["caax"]["propeptide"] == [662, 664]


def test_zmpste24_site_is_inside_the_deletion(b):
    assert (b["cleavage"]["after"], b["cleavage"]["before"]) == pg.EXPECTED["cleavage"] == (646, 647)
    assert b["cleavage"]["residues"] == "YL"
    assert b["cleavage"]["inside_deletion"]
    assert b["mature"]["end"] == pg.EXPECTED["mature_end"] == 646


def test_lamin_c_and_landmarks(b):
    assert (b["lamin_c"]["shared_end"], b["lamin_c"]["tail"], b["lamin_c"]["length"]) == (566, "VSGSRR", 572)
    assert b["lamin_c"]["ids"] == ["P02545-2"]
    names = {x["name"]: (x["start"], x["end"]) for x in b["landmarks"]}
    assert names["NLS"] == pg.EXPECTED["nls"] == (417, 422)
    assert names["Ig fold (LTD)"] == pg.EXPECTED["ig_fold"] == (428, 545)
    assert names["IF rod"] == pg.EXPECTED["rod"] == (31, 387)
    assert names["Coil 1A"] == (34, 70) and names["Coil 1B"] == (81, 218) and names["Coil 2"] == (243, 383)
    assert names["Tail"][1] == 664


def test_mouse_mapping(b):
    m = b["mouse"]
    assert m["accession"] == "P48678" and m["length"] == 665
    assert m["human_pos"] == 608 and m["mouse_pos"] == 609 and m["offset"] == 1
    # the alignment that justifies the offset: position 608 of the human window is G, and so is mouse 609
    hw, (hs, _) = m["human_window"], m["human_window_range"]
    mw, (ms, _) = m["mouse_window"], m["mouse_window_range"]
    assert hw[608 - hs] == "G" and mw[609 - ms] == "G"
    assert sum(a == c for a, c in zip(hw, mw)) == 12 and m["identity"] == "12 of 16"
    assert m["farnesyl_cys"] == b["caax"]["farnesyl"] + m["offset"]


def test_mutation_record_and_relations(b):
    assert b["mutation"]["cdna"] == "c.1824C>T" and b["mutation"]["clinvar"] == "VCV000014500" and b["mutation"]["dbsnp"] == "rs58596362"
    assert b["mutation"]["classification"].startswith("Pathogenic")
    assert pg.relation(482, 607, 656, 50) == {"where": "before", "text": pg.relation(482, 607, 656, 50)["text"], "progerin_pos": 482}
    assert pg.relation(644, 607, 656, 50)["where"] == "inside" and pg.relation(644, 607, 656, 50)["progerin_pos"] is None
    assert pg.relation(661, 607, 656, 50)["progerin_pos"] == 611
    assert any(v["label"].startswith("G608") for v in b["hgps_variants"])
    assert any(v["label"].startswith("R644") for v in b["uniprot_variants_in_deletion"])


def test_confidence_is_reported_honestly(b):
    c = b["confidence"]
    assert c and c["deleted"] is not None and c["ig_fold"] is not None
    assert c["deleted"] < c["ig_fold"], "the deleted tail segment must be less confident than the Ig fold in the AlphaFold model"
    assert "not a reliable fold" in b["texts"]["model_caveat"]


def test_exports(item, b):
    text, name = core.progerin_fasta(item)
    lines = text.strip().split("\n")
    assert lines[0].startswith(">progerin|P02545-6") and name.endswith(".fasta")
    assert "".join(lines[1:]) == b["progerin"]["sequence"] and len("".join(lines[1:])) == 614
    script, sname = core.progerin_script(item, "AF")
    assert sname.endswith("_progerin.pml")
    assert "hgps_deleted" in script and "zmpste24_site" in script and "farnesyl_cys" in script
    assert "resi 607" in script.replace("+", " ").replace("resi 607+", "resi 607 ") or "607+" in script


def test_not_lmna_has_no_progerin_mode(res):
    u = copy.deepcopy(res["uniprot"])
    u["accession"] = "P04637"           # the same features under another accession: still not LMNA
    assert pg.block(u, res.get("alphafold")) is None
    assert pg.block(None) is None
    assert pg.block({"accession": "P48678", "sequence": "M" * 665, "features": [], "isoforms": [], "variants": []}) is None


def test_exports_refuse_a_non_lmna_item(tmp_path, monkeypatch):
    monkeypatch.setattr(core, "load_result", lambda item: {"kind": "protein", "accession": "P04637"})
    with pytest.raises(UserFacingError):
        core.progerin_fasta("x")


def test_methods_and_citations_mention_the_sources(res):
    heads = [h for h, _ in core.methods(res)]
    assert "Progerin (HGPS) mode" in heads
    cites = core.citations(res)
    assert any("Eriksson" in c for c in cites) and any("Osorio" in c for c in cites) and any("VCV000014500" in c for c in cites)
