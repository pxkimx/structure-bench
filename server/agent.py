"""In-app Claude assistant. Streams a tool-using conversation over SSE (same protocol as RNAseq Bench:
{type: text|tool|tool_result|done|error}). The tools call the same core functions the UI uses."""
from __future__ import annotations

import json
import os
import threading
from typing import Iterator

from . import core
from .common import UserFacingError, list_items, load_settings

DEFAULT_MODEL = "claude-sonnet-5"
FALLBACK_MODELS = ["claude-sonnet-5", "claude-opus-5-5", "claude-haiku-4-5-20251001"]


def available_models() -> dict:
    cfg = load_settings()
    key = cfg.get("api_key") or os.environ.get("ANTHROPIC_API_KEY")
    if not key:
        return {"ok": False, "error": "no key", "models": FALLBACK_MODELS}
    try:
        import anthropic
        client = anthropic.Anthropic(api_key=key)
        ids = [m.id for m in client.models.list(limit=100)]
        ids = [i for i in ids if i.startswith("claude")]
        ids = sorted(ids, key=lambda i: (0 if "sonnet" in i else 1 if "opus" in i else 2))
        return {"ok": True, "models": ids or FALLBACK_MODELS}
    except Exception as e:  # noqa: BLE001
        return {"ok": False, "error": f"{type(e).__name__}: {str(e)[:200]}", "models": FALLBACK_MODELS}


def _resolve_model(client, wanted: str) -> str:
    """Use the configured model if the key can see it; otherwise the newest Sonnet (models get retired)."""
    try:
        ids = [m.id for m in client.models.list(limit=100)]
    except Exception:  # noqa: BLE001
        return wanted
    if wanted in ids:
        return wanted
    for pat in ("sonnet", "opus", "haiku"):
        for i in ids:
            if pat in i:
                return i
    return wanted


SYSTEM = """You are the assistant built into Structure Bench, a local app that takes a biologist from gene to protein to
structure: UniProt annotation, the AlphaFold model with its per-residue confidence (pLDDT) and predicted aligned error
(PAE), experimental PDB entries, AlphaMissense scores, and variant mapping (solvent accessibility by Shrake–Rupley,
contacts within 5 Å, secondary structure, domain, interfaces), plus structure comparison (CE alignment, per-residue CA
deviation). You talk to a bench biologist looking at a results page. Explain in plain language, be concise, and quote
the numbers from tool results — never invent a value, a residue number or a PDB ID; read them with the tools first.
Be honest about limits: pLDDT < 50 means the structure says nothing there; high pLDDT with high PAE between two
domains means their arrangement is unknown; relative SASA on a single-chain model ignores partners (coiled coils,
complexes); secondary structure without the DSSP program is a phi/psi approximation; AlphaMissense is a prediction,
not a clinical classification; an interface in a PDB entry can be a crystal contact. When the user asks for a script,
call pymol_script and give the script in a code block. Keep answers short unless asked for depth."""

TOOLS = [
    {"name": "lookup_protein", "description": "Look up a protein by gene symbol + species (human, mouse, rat, zebrafish, fly, worm, yeast, ecoli or an NCBI taxon id) or by UniProt accession. Fetches UniProt, the AlphaFold model, PAE and AlphaMissense (cached; works offline for anything looked up before). Returns the item id, or a list of choices when several reviewed entries match.",
     "input_schema": {"type": "object", "properties": {"gene": {"type": "string"}, "accession": {"type": "string"}, "species": {"type": "string"}}}},
    {"name": "list_items", "description": "Saved analyses (items), newest first: id, kind (protein/upload), name.",
     "input_schema": {"type": "object", "properties": {}}},
    {"name": "get_item_result", "description": "Summary of an item: tiles, findings (flags), UniProt function/domains/diseases, AlphaFold confidence, PAE between domains, the per-panel 'yours' readings, loaded structures, the last variant mapping and comparisons.",
     "input_schema": {"type": "object", "properties": {"item": {"type": "string"}}, "required": ["item"]}},
    {"name": "list_pdb_entries", "description": "Experimental structures UniProt lists for this protein: PDB id, method, resolution, chains and the residue ranges they span, and whether each is loaded.",
     "input_schema": {"type": "object", "properties": {"item": {"type": "string"}}, "required": ["item"]}},
    {"name": "load_pdb_entry", "description": "Download a PDB entry (mmCIF from RCSB, cached) into the item so it can be used for mapping and comparison. Returns its structure id (PDB-XXXX) and chains.",
     "input_schema": {"type": "object", "properties": {"item": {"type": "string"}, "pdb_id": {"type": "string"}}, "required": ["item", "pdb_id"]}},
    {"name": "map_variants", "description": "Map protein variants (R482W, p.Arg482Trp, G608G…, newline- or comma-separated) onto a structure of the item ('AF' = AlphaFold model, default; 'PDB-1IFR' = a loaded PDB entry). Checks the reference residue, then gives pLDDT/B-factor, relative SASA and burial class, neighbours within 5 Å, secondary structure, domain, UniProt annotation, AlphaMissense, interface, and heuristic flags. The result also appears on the page.",
     "input_schema": {"type": "object", "properties": {"item": {"type": "string"}, "variants": {"type": "string"}, "structure": {"type": "string"}}, "required": ["item", "variants"]}},
    {"name": "residue_environment", "description": "One residue's environment in a structure: its pLDDT/B-factor, SASA, secondary structure, domain and every residue, other chain or ligand within `radius` Å (default 5). `residue` is a UniProt position ('482'), a variant ('R482W') or a structure key ('A:482').",
     "input_schema": {"type": "object", "properties": {"item": {"type": "string"}, "residue": {"type": "string"}, "radius": {"type": "number"}, "structure": {"type": "string"}}, "required": ["item", "residue"]}},
    {"name": "compare_structures", "description": "Compare two structures of the item (structure ids such as 'AF', 'PDB-1IFR', 'UP-1'): CE RMSD, RMSD after superposing residue pairs matched by UniProt position, per-residue CA deviation, and where large deviations coincide with low pLDDT. Load PDB entries first.",
     "input_schema": {"type": "object", "properties": {"item": {"type": "string"}, "a": {"type": "string"}, "b": {"type": "string"}}, "required": ["item", "a", "b"]}},
    {"name": "alphamissense", "description": "AlphaMissense pathogenicity score and class for missense variants of the canonical human isoform (Cheng et al. 2023).",
     "input_schema": {"type": "object", "properties": {"item": {"type": "string"}, "variants": {"type": "string"}}, "required": ["item", "variants"]}},
    {"name": "pymol_script", "description": "A PyMOL (.pml) or ChimeraX (.cxc) script that loads the structure, colours it (AlphaFold confidence or by chain) and shows the given variants as labelled sticks.",
     "input_schema": {"type": "object", "properties": {"item": {"type": "string"}, "variants": {"type": "string"}, "structure": {"type": "string"}, "program": {"type": "string", "enum": ["pymol", "chimerax"]}}, "required": ["item"]}},
]


def _slim_variants(v: dict) -> dict:
    keep = ("label", "input", "error", "kind", "pos", "num", "plddt", "bfactor", "rsa", "rsa_alone", "rsa_class", "ss_word",
            "domain", "am_score", "am_class", "uniprot_known", "interface", "ligand_contacts", "unresolved", "change", "neighbours")
    return {"structure": v["structure"], "structure_label": v["structure_label"], "chain": v["chain"], "yours": v["yours"],
            "ss_method": v.get("ss_method"),
            "rows": [{k: r.get(k) for k in keep if r.get(k) not in (None, [], "")} | {"flags": [f["text"] for f in r.get("flags", [])]}
                     for r in v["rows"]]}


def run_tool(name: str, inp: dict) -> str:
    try:
        if name == "lookup_protein":
            r = core.lookup(inp.get("gene"), inp.get("species") or "human", inp.get("accession"))
            if "item" in r:
                s = core.item_summary(r["item"])
                return json.dumps({"item": r["item"], "name": s.get("name"), "flags": s.get("flags"), "tiles": s.get("tiles")})
            return json.dumps(r)
        if name == "list_items":
            return json.dumps(list_items(20))
        item = inp.get("item")
        if name == "get_item_result":
            return json.dumps(core.item_summary(item))[:60000]
        if name == "list_pdb_entries":
            r = core.list_pdb_entries(item)
            r["entries"] = [{k: e.get(k) for k in ("id", "method", "resolution", "segments", "loaded")} for e in r["entries"]]
            return json.dumps(r)[:40000]
        if name == "load_pdb_entry":
            r = core.load_pdb(item, inp["pdb_id"])
            m = r["structure"]
            return json.dumps({"sid": r["sid"], "main_chain": m["main_chain"], "chains": m["chains"], "partners": m.get("partners"),
                               "method": m.get("method"), "resolution": m.get("resolution")})
        if name == "map_variants":
            return json.dumps(_slim_variants(core.map_variants(item, inp["variants"], inp.get("structure"))))[:60000]
        if name == "residue_environment":
            r = core.residue_environment(item, str(inp["residue"]), float(inp.get("radius") or 5), inp.get("structure"))
            return json.dumps(r)[:40000]
        if name == "compare_structures":
            c = core.compare_structures(item, inp["a"], inp["b"])
            c = {k: v for k, v in c.items() if k not in ("deviation", "how")}
            return json.dumps(c)
        if name == "alphamissense":
            return json.dumps(core.alphamissense(item, inp["variants"]))
        if name == "pymol_script":
            text, fname = core.viewer_script(item, inp.get("program") or "pymol", inp.get("structure"), inp.get("variants"))
            return f"file name: {fname}\n\n{text}"
        return f"unknown tool {name}"
    except UserFacingError as e:
        return f"error: {e}"
    except Exception as e:  # noqa: BLE001
        from .common import log_exc
        log_exc(f"agent tool {name}")
        return f"error: {type(e).__name__}: {e}"


# ---------------------------------------------------------------- conversations
CONV: dict[str, list] = {}
LOCK = threading.Lock()


def chat(conv_id: str, user_text: str, item: str | None) -> Iterator[str]:
    """Yields SSE 'data:' lines: {type: text|tool|tool_result|done|error}."""
    try:
        import anthropic
    except ImportError:
        yield _sse({"type": "error", "text": "The anthropic package is missing: pip install anthropic"}); return
    cfg = load_settings()
    key = cfg.get("api_key") or os.environ.get("ANTHROPIC_API_KEY")
    if not key:
        yield _sse({"type": "error", "text": "No API key. Open Settings (gear icon) and paste an Anthropic API key."}); return
    client = anthropic.Anthropic(api_key=key)
    model = _resolve_model(client, cfg.get("model") or DEFAULT_MODEL)
    with LOCK:
        hist = CONV.setdefault(conv_id, [])
    ctx = ""
    if item:
        try:
            s = core.item_summary(item)
            ctx = (f"\n\nThe user is currently viewing item '{item}' ({s['kind']}): {s.get('name')} — {s.get('title')}. "
                   "Key findings:\n" + "\n".join(f"- [{f['level']}] {f['text']}" for f in s.get("flags") or [])
                   + "\nLoaded structures: " + ", ".join(f"{x['sid']} ({x['label']})" for x in s.get("structures") or [])
                   + ("\nA variant mapping is on the page (use get_item_result or map_variants to read it)." if s.get("variants") else ""))
        except Exception:  # noqa: BLE001
            ctx = f"\n\nThe user is viewing item '{item}'."
    hist.append({"role": "user", "content": user_text})
    for _ in range(12):  # tool-use rounds
        try:
            with client.messages.stream(model=model, max_tokens=4000, system=SYSTEM + ctx, tools=TOOLS, messages=hist) as stream:
                for ev in stream:
                    if ev.type == "content_block_delta" and getattr(ev.delta, "type", "") == "text_delta":
                        yield _sse({"type": "text", "text": ev.delta.text})
                msg = stream.get_final_message()
        except Exception as e:  # noqa: BLE001
            m = str(e)
            if "authentication" in m.lower() or "invalid x-api-key" in m.lower():
                m = "The API key was rejected. Open Settings & API key and paste a valid key from console.anthropic.com."
            elif "credit" in m.lower() or "billing" in m.lower():
                m = "Anthropic says this key has no credit. Add a prepaid balance at console.anthropic.com → Billing."
            yield _sse({"type": "error", "text": f"{type(e).__name__}: {m}"}); hist.pop(); return
        hist.append({"role": "assistant", "content": [b.model_dump() for b in msg.content]})
        tool_uses = [b for b in msg.content if b.type == "tool_use"]
        if not tool_uses:
            break
        results = []
        for tu in tool_uses:
            yield _sse({"type": "tool", "name": tu.name, "input": tu.input})
            out = run_tool(tu.name, tu.input)
            yield _sse({"type": "tool_result", "name": tu.name, "text": out[:1500]})
            results.append({"type": "tool_result", "tool_use_id": tu.id, "content": out})
        hist.append({"role": "user", "content": results})
    yield _sse({"type": "done"})


def _sse(d: dict) -> str:
    return f"data: {json.dumps(d)}\n\n"


def reset(conv_id: str):
    with LOCK:
        CONV.pop(conv_id, None)
