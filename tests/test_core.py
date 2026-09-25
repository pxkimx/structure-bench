"""Offline checks of the core computations on the bundled LMNA example (examples/lmna).

    .venv/bin/python -m pytest tests -q
"""
import json

import numpy as np
import pytest
from conftest import AM, MODEL, PAE, UNIPROT

from server import alphafold as afdb
from server import core, net
from server import structure as st
from server import uniprot as up
from server import variants as vr


# ---------------------------------------------------------------- variant notation
@pytest.mark.parametrize("text,ref,pos,alt,kind", [
    ("R482W", "R", 482, "W", "missense"),
    ("p.R482W", "R", 482, "W", "missense"),
    ("p.Arg482Trp", "R", 482, "W", "missense"),
    ("Arg482Trp", "R", 482, "W", "missense"),
    ("p.(Arg482Trp)", "R", 482, "W", "missense"),
    ("r482w", "R", 482, "W", "missense"),
    ("G608G", "G", 608, "G", "synonymous"),
    ("p.Gly608=", "G", 608, "G", "synonymous"),
    ("p.Gly608Gly", "G", 608, "G", "synonymous"),
    ("R482*", "R", 482, "*", "nonsense"),
    ("p.Arg482Ter", "R", 482, "*", "nonsense"),
])
def test_parse_formats(text, ref, pos, alt, kind):
    v = vr.parse(text)
    assert "error" not in v, v
    assert (v["ref"], v["pos"], v["alt"], v["kind"]) == (ref, pos, alt, kind)


@pytest.mark.parametrize("bad", ["c.1824C>T", "p.R482fs", "482W", "RW", "X12Y", "p.Arg482Xyz"])
def test_parse_rejects(bad):
    assert "error" in vr.parse(bad)


def test_parse_many_splits_and_dedupes():
    vs = vr.parse_many("R482W, p.Arg482Trp\nE145K;G608G  R527P")
    assert [v["label"] for v in vs] == ["R482W", "E145K", "G608G", "R527P"]


# ---------------------------------------------------------------- reference check
def _u():
    return up.parse_entry(json.loads(UNIPROT.read_text()))


def test_reference_ok_and_wrong():
    u = _u()
    assert vr.check_reference(vr.parse("R482W"), u["sequence"])["ok"]
    bad = vr.check_reference(vr.parse("A482W"), u["sequence"])
    assert not bad["ok"]
    assert "Arg (R), not Ala (A)" in bad["message"]


def test_reference_signal_peptide_offset():
    # a toy protein with a 20-residue signal peptide: mature-chain numbering is 20 lower than UniProt's
    seq = "M" + "A" * 19 + "KLVRSTYW" * 10
    v = vr.parse("R4W")                                 # residue 4 of the mature chain = UniProt 24 = R
    assert seq[23] == "R"
    chk = vr.check_reference(v, seq, signal_len=20)
    assert not chk["ok"]
    assert any("signal peptide" in s and "R24W" in s for s in chk["suggest"]), chk


def test_reference_beyond_end():
    chk = vr.check_reference(vr.parse("R999W"), _u()["sequence"])
    assert not chk["ok"] and "beyond the end" in chk["message"]


# ---------------------------------------------------------------- structure
@pytest.fixture(scope="module")
def af():
    s, fmt = st.load(MODEL)
    assert fmt == "cif"
    return s


def test_plddt_equals_ca_bfactor_and_r482w_maps_to_482(af):
    item = core.lookup(accession="P02545")["item"]
    v = core.map_variants(item, "R482W", "AF")
    row = v["rows"][0]
    assert row["pos"] == 482 and row["num"] == 482 and row["key"] == "A:482"
    ca = af[0]["A"][482]["CA"]
    assert af[0]["A"][482].get_resname() == "ARG"
    assert row["plddt"] == pytest.approx(ca.get_bfactor(), abs=1e-6)
    assert row["am_score"] == pytest.approx(0.9191, abs=1e-4)      # AlphaMissense table value for R482W


def test_relative_sasa_range(af):
    ch = st.polymer_chains(af[0])[0]
    tab = st.residue_table(af[0], "A", ch["residues"], list(range(1, 665)))
    rsa = np.array([r["rsa"] for r in tab["rows"]], dtype=float)
    assert len(rsa) == 664 and np.isfinite(rsa).all()
    assert rsa.min() >= 0 and rsa.max() <= 1.3
    assert np.median(rsa) < 1.0
    # a single chain: alone and 'in the model' must be the same
    assert all(r["rsa"] == r["rsa_alone"] for r in tab["rows"])


def test_neighbour_counts(af):
    nb = st.neighbours(af[0], "A:482", 5.0)
    assert 4 <= nb["n_total"] <= 30                 # a surface arginine in a β-sandwich
    assert all(x["min_dist"] <= 5.0 for x in nb["same_chain"])
    assert {"A:481", "A:483"} <= {x["key"] for x in nb["same_chain"]}   # its sequence neighbours are always there
    assert nb["n_nonlocal"] <= nb["n_total"]
    big = st.neighbours(af[0], "A:482", 10.0)
    assert big["n_total"] > nb["n_total"]


def test_calcium_ion_is_a_ligand_not_a_residue_with_a_ca_atom():
    """is_ligand() tells a modified amino acid that stays in the chain (MSE: a full N-CA-C backbone) from a small
    HETATM group. A bound Ca2+ ion's PDB atom name is also 'CA' (same string as the alpha-carbon atom), so a check
    of 'CA' in r alone calls it part of the chain instead of a ligand — hiding calcium-site contacts (EF-hands,
    C2 domains) behind a same-chain or interface reading instead of 'near ligand Ca'."""
    import io
    import warnings

    from Bio.PDB import PDBParser
    pdb = ("ATOM      1  N   ALA A   1      11.104  13.207   2.100  1.00 20.00           N\n"
           "ATOM      2  CA  ALA A   1      12.560  13.207   2.100  1.00 20.00           C\n"
           "ATOM      3  C   ALA A   1      13.000  14.600   2.100  1.00 20.00           C\n"
           "ATOM      4  O   ALA A   1      12.300  15.600   2.100  1.00 20.00           O\n"
           "HETATM    5 CA    CA A 101      20.000  20.000  20.000  1.00 30.00          CA\n")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = PDBParser(QUIET=True).get_structure("x", io.StringIO(pdb))
    ala, ca_ion = list(s[0]["A"])
    assert st.is_ligand(ala) is False                  # a real residue, not a ligand
    assert st.is_ligand(ca_ion) is True                 # the calcium ion, not part of the chain


def test_ce_self_alignment_rmsd_zero(af):
    out = st.ce_rmsd(af[0], "A", af[0], "A")
    assert out["rmsd"] == pytest.approx(0.0, abs=1e-3)


def test_mapping_undoes_offset(af, tmp_path):
    """Renumber a copy of the model by +10 and check that the sequence mapping gives back UniProt numbering."""
    s2, _ = st.load(MODEL, "shifted")
    chain = s2[0]["A"]
    for r in sorted(list(chain), key=lambda r: -r.id[1]):          # highest first so new ids never collide
        r.id = (r.id[0], r.id[1] + 10, r.id[2])
    out = tmp_path / "shifted.cif"
    st.write(s2, out, "cif")
    s3, _ = st.load(out, "reloaded")
    ch = st.polymer_chains(s3[0])[0]
    assert ch["first"] == 11 and ch["last"] == 674
    m = st.map_sequence(ch["sequence"], _u()["sequence"])
    assert m["identity"] == 1.0 and m["aligned"] == 664
    assert all(p == r.id[1] - 10 for p, r in zip(m["pos"], ch["residues"]))


def test_mapping_fragment_and_tag():
    u = _u()["sequence"]
    frag = "MGSSHHHHHH" + u[430:545]                     # a His-tagged LTD construct
    m = st.map_sequence(frag, u)
    assert m["pos"][:10] == [None] * 10                 # the tag maps nowhere
    assert m["pos"][10] == 431 and m["pos"][-1] == 545


# ---------------------------------------------------------------- PAE
def test_pae_new_layout():
    obj = json.loads(PAE.read_text())
    m, mx = afdb.parse_pae(obj)
    assert m.shape == (664, 664) and mx == pytest.approx(31.75)


def test_pae_old_layout():
    n = 4
    r1 = [i for i in range(1, n + 1) for _ in range(n)]
    r2 = [j for _ in range(n) for j in range(1, n + 1)]
    dist = [float(i * 10 + j) for i, j in zip(r1, r2)]
    m, mx = afdb.parse_pae([{"residue1": r1, "residue2": r2, "distance": dist, "max_predicted_aligned_error": 31.75}])
    assert m.shape == (4, 4) and m[1, 2] == 23.0 and m[3, 0] == 41.0
    m2, _ = afdb.parse_pae({"predicted_aligned_error": [[0, 1], [2, 0]]})       # bare dict, no list
    assert m2.shape == (2, 2)
    with pytest.raises(ValueError):
        afdb.parse_pae([{"something": 1}])


def test_pae_domain_summary_lmna():
    m, _ = afdb.parse_pae(json.loads(PAE.read_text()))
    s = afdb.summarize_pae(m, up.main_domains(_u()))
    names = [b["name"] for b in s["blocks"]]
    assert names == ["IF rod", "LTD"]
    ltd = s["blocks"][1]
    assert ltd["within"] < 10                       # a compact domain
    assert s["pairs"][0]["call"] == "not determined"  # rod ↔ LTD placement is arbitrary in the model


# ---------------------------------------------------------------- AlphaMissense
def test_alphamissense_parse():
    am = afdb.parse_am(AM.read_text())
    assert len(am) == 664 * 19
    assert am["M1A"] == (pytest.approx(0.7048), "likely pathogenic")
    means = afdb.am_residue_means(am, 664)
    assert len(means) == 664 and all(0 <= x <= 1 for x in means)


# ---------------------------------------------------------------- the item as a whole
def test_offline_lookup_and_example_variants():
    assert net.offline()
    item = core.lookup(gene="LMNA", species="human")["item"]
    res = core.load_result(item)
    assert res["accession"] == "P02545" and res["alphafold"]["entry"] == "AF-P02545-F1"
    v = core.map_variants(item, "E145K\nR482W\nR527P\nG608G", "AF")
    rows = {r["label"]: r for r in v["rows"]}
    assert set(rows) == {"E145K", "R482W", "R527P", "G608G"}
    assert not any(r.get("error") for r in rows.values())
    assert rows["G608G"]["kind"] == "synonymous"
    assert "splice" in rows["G608G"]["flags"][0]["text"]
    assert rows["E145K"]["uniprot_known"]                         # UniProt lists E145K (HGPS, atypical)
    assert any("coiled coil" in f["text"] for f in rows["E145K"]["flags"])
    assert rows["R527P"]["am_class"] == "likely pathogenic"


def test_wrong_reference_is_reported_in_mapping():
    item = core.lookup(accession="P02545")["item"]
    v = core.map_variants(item, "A482W", "AF", save=False)
    assert v["rows"][0]["error"].startswith("Position 482 in P02545")


def test_scripts_and_exports():
    item = core.lookup(accession="P02545")["item"]
    core.map_variants(item, "R482W\nR482Q", "AF")
    pml, name = core.viewer_script(item, "pymol", "AF")
    assert name.endswith(".pml") and 'label sb_AF_P02545_F1_model_v6 and chain A and resi 482 and name CA, "R482W/Q"' in pml
    cxc, _ = core.viewer_script(item, "chimerax", "AF")
    assert "color bfactor palette alphafold" in cxc and "/A:482" in cxc
    p, _ = core.export_structure(item, "AF", "flag", "pdb")
    flagged = {int(l[22:26]) for l in p.read_text().splitlines() if l.startswith("ATOM") and float(l[60:66]) == 1.0}
    assert flagged == {482}
    p2, _ = core.export_structure(item, "AF", "am", "cif")
    s, _ = st.load(p2)
    assert 0 <= s[0]["A"][482]["CA"].get_bfactor() <= 1
