"""Offline self-test: every analysis and every assistant tool on the bundled LMNA example.

    .venv/bin/python -m server.selftest

Runs in a temporary SB_HOME with the network switched off (SB_OFFLINE=1), so it checks the installation and the
offline path at the same time. Takes about 10–20 seconds.
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
import time
import traceback

os.environ.setdefault("SB_OFFLINE", "1")
os.environ["SB_HOME"] = tempfile.mkdtemp(prefix="sb_selftest_")

from . import agent, core, report  # noqa: E402
from .common import EXAMPLES  # noqa: E402

RESULTS = []


def check(name, fn):
    t0 = time.time()
    try:
        msg = fn() or ""
        RESULTS.append((True, name))
        print(f"  PASS  {name:<48} {time.time() - t0:5.2f}s  {msg}", flush=True)
    except Exception as e:  # noqa: BLE001
        RESULTS.append((False, name))
        print(f"  FAIL  {name:<48} {type(e).__name__}: {e}", flush=True)
        traceback.print_exc()


def main():
    print(f"Structure Bench self-test — offline={os.environ['SB_OFFLINE']} home={os.environ['SB_HOME']}")
    st = {}

    def lookup():
        r = core.lookup(gene="LMNA", species="human")
        st["item"] = r["item"]
        res = core.load_result(r["item"])
        assert res["accession"] == "P02545"
        assert res["alphafold"] and res["pae"] and res["alphamissense"]["available"]
        return f"{res['name']}, mean pLDDT {res['alphafold']['mean_plddt']}"
    check("lookup LMNA (gene → UniProt → AlphaFold, PAE, AM)", lookup)

    def by_acc():
        r = core.lookup(accession="P02545")
        assert r["item"] == st["item"] and r.get("reused")
        return "reused the saved lookup"
    check("lookup by accession reuses the item", by_acc)

    def variants():
        v = core.map_variants(st["item"], "E145K\nR482W\nR527P\nG608G\nA482W", "AF")
        bad = [r for r in v["rows"] if r.get("error")]
        assert len(bad) == 1 and bad[0]["input"] == "A482W"
        return v["yours"][:90] + "…"
    check("map_variants (E145K, R482W, R527P, G608G + a wrong one)", variants)

    def env():
        e = core.residue_environment(st["item"], "R482W", 5)
        assert e["residue"]["pos"] == 482 and e["neighbours"]["n_total"] > 3
        return f"{e['neighbours']['n_total']} residues within 5 Å of R482"
    check("residue_environment", env)

    def am():
        a = core.alphamissense(st["item"], "R482W, G608G")
        assert a["available"] and a["rows"][0]["class"] == "likely pathogenic"
        return f"R482W {a['rows'][0]['score']} {a['rows'][0]['class']}"
    check("alphamissense", am)

    def pdb_list():
        p = core.list_pdb_entries(st["item"])
        assert len(p["entries"]) >= 20
        return f"{len(p['entries'])} PDB entries listed"
    check("list_pdb_entries", pdb_list)

    def pdb_offline():
        try:
            core.load_pdb(st["item"], "1IFR")
        except core.UserFacingError as e:
            assert "could not be downloaded" in str(e)
            return "offline → a plain message, as it should"
        return "loaded (was cached)"
    check("load_pdb while offline degrades to a message", pdb_offline)

    def upload_into():
        data = (EXAMPLES / "lmna" / "AF-P02545-F1-model_v6.cif").read_bytes()
        r = core.upload("my_model.cif", data, st["item"])
        st["up"] = r["sid"]
        return r["sid"]
    check("upload a structure into the item", upload_into)

    def compare():
        c = core.compare_structures(st["item"], "AF", st["up"])
        assert c["ce"]["rmsd"] < 0.01 and c["rmsd_mapped"] < 0.01 and c["n_pairs"] == 664
        return f"CE RMSD {c['ce']['rmsd']} Å on itself"
    check("compare_structures (CE + per-residue deviation)", compare)

    def upload_new():
        data = (EXAMPLES / "lmna" / "AF-P02545-F1-model_v6.cif").read_bytes()
        r = core.upload("offline_upload.cif", data)
        v = core.map_variants(r["item"], "R482W", None)
        assert v["rows"][0]["num"] == 482
        return f"new item {r['item']}, residues numbered as in the file"
    check("upload as a new offline item + mapping", upload_new)

    def exports():
        p, _ = core.export_structure(st["item"], "AF", "am", "cif")
        q, _ = core.export_structure(st["item"], "AF", "flag", "pdb")
        pml, _ = core.viewer_script(st["item"], "pymol", "AF")
        cxc, _ = core.viewer_script(st["item"], "chimerax", "AF")
        assert "R482W" in pml and "R482W" in cxc
        return f"{p.name}, {q.name}, .pml, .cxc"
    check("exports (B-factor files, PyMOL, ChimeraX)", exports)

    def pdf():
        p = report.build_pdf(st["item"])
        assert p.stat().st_size > 20000
        return f"{p.stat().st_size // 1024} KB"
    check("PDF report", pdf)

    def tools():
        outs = {}
        for name, inp in [("list_items", {}), ("get_item_result", {"item": st["item"]}),
                          ("list_pdb_entries", {"item": st["item"]}),
                          ("map_variants", {"item": st["item"], "variants": "R482W", "structure": "AF"}),
                          ("residue_environment", {"item": st["item"], "residue": "482"}),
                          ("compare_structures", {"item": st["item"], "a": "AF", "b": st["up"]}),
                          ("alphamissense", {"item": st["item"], "variants": "R482W"}),
                          ("pymol_script", {"item": st["item"], "variants": "R482W"}),
                          ("lookup_protein", {"gene": "LMNA", "species": "human"})]:
            out = agent.run_tool(name, inp)
            assert not out.startswith("error"), f"{name}: {out[:200]}"
            outs[name] = len(out)
        json.loads(agent.run_tool("get_item_result", {"item": st["item"]}))
        return f"{len(outs)} tools answered"
    check("assistant tools (without calling Claude)", tools)

    ok = sum(1 for r in RESULTS if r[0])
    print(f"\n{ok}/{len(RESULTS)} checks passed.")
    return 0 if ok == len(RESULTS) else 1


if __name__ == "__main__":
    sys.exit(main())
