"""Structure Bench web server.  Run:  uvicorn server.app:app --port 8767"""
from __future__ import annotations

import json
import os
import threading
import traceback

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse, Response, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from . import agent, core, report, uniprot
from .common import (HOME, ROOT, VERSION, UserFacingError, item_dir, json_default, list_items, load_result, load_settings,
                     save_settings)

app = FastAPI(title="Structure Bench", version=VERSION)


def J(x) -> Response:
    """JSON with numpy scalars converted (FastAPI's encoder cannot serialise numpy.int32 from Bio.PDB residue ids)."""
    return Response(json.dumps(x, default=json_default), media_type="application/json")
print(f"Structure Bench {VERSION} — data in {HOME} — python {__import__('platform').python_version()}", flush=True)


@app.exception_handler(UserFacingError)
async def _user_error(_: Request, exc: UserFacingError):
    return JSONResponse({"detail": str(exc)}, status_code=400)


@app.exception_handler(Exception)
async def _any_error(_: Request, exc: Exception):
    tb = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
    print(tb, flush=True)
    return JSONResponse({"detail": f"Something went wrong inside the app ({type(exc).__name__}: {str(exc)[:200]}). "
                                   "The details are in the server log.", "trace": tb[-3000:]}, status_code=500)


# ---------------------------------------------------------------- look up
class LookupIn(BaseModel):
    gene: str | None = None
    species: str | None = "human"
    accession: str | None = None
    refresh: bool = False
    reuse: bool = True


@app.post("/api/lookup")
def lookup(body: LookupIn):
    return J(core.lookup(body.gene, body.species or "human", body.accession, body.refresh, body.reuse and not body.refresh))


@app.get("/api/species")
def species():
    return [{"key": k, "taxon": v[0], "name": v[1]} for k, v in uniprot.SPECIES.items()]


@app.post("/api/example")
def example():
    """The bundled LMNA example (works offline) with four variants mapped on the AlphaFold model."""
    r = core.lookup(accession="P02545")
    core.map_variants(r["item"], "E145K\nR482W\nR527P\nG608G", "AF")
    return J(r)


# ---------------------------------------------------------------- items
@app.get("/api/items")
def items(limit: int = 60):
    return list_items(limit)


@app.get("/api/items/{item}")
def get_item(item: str):
    r = load_result(item)
    r["item"] = item
    r["citations"] = core.citations(r)
    return J(r)


@app.get("/api/items/{item}/residues/{sid}")
def residues(item: str, sid: str, chain: str | None = None):
    res = load_result(item)
    meta = res["structures"].get(sid)
    if not meta:
        raise HTTPException(404, f"No structure {sid} in this item.")
    d = core.residues_doc(item_dir(item), sid, chain or meta["main_chain"])
    if not d:
        core._chain_table(item_dir(item), res, meta, chain or meta["main_chain"])
        d = core.residues_doc(item_dir(item), sid, chain or meta["main_chain"])
    return J(d)


@app.get("/api/items/{item}/pae")
def pae(item: str):
    p = item_dir(item) / "pae.json"
    if not p.exists():
        raise HTTPException(404, "No PAE for this item.")
    return FileResponse(p, media_type="application/json")


@app.get("/api/items/{item}/files/{path:path}")
def files(item: str, path: str):
    d = item_dir(item)
    p = (d / path).resolve()
    if not str(p).startswith(str(d.resolve())) or not p.is_file():
        raise HTTPException(404)
    return FileResponse(p, filename=p.name if p.suffix in (".cif", ".pdb", ".bcif", ".pdf", ".csv", ".gz") else None)


class PdbIn(BaseModel):
    pdb_id: str


@app.post("/api/items/{item}/pdb")
def load_pdb(item: str, body: PdbIn):
    return J(core.load_pdb(item, body.pdb_id))


class IsoIn(BaseModel):
    entry: str


@app.post("/api/items/{item}/isoform")
def load_isoform(item: str, body: IsoIn):
    return J(core.load_isoform(item, body.entry))


class VarIn(BaseModel):
    text: str
    structure: str | None = None
    chain: str | None = None


@app.post("/api/items/{item}/variants")
def variants(item: str, body: VarIn):
    return J(core.map_variants(item, body.text, body.structure, body.chain))


class EnvIn(BaseModel):
    residue: str
    radius: float = 5.0
    structure: str | None = None


@app.post("/api/items/{item}/environment")
def environment(item: str, body: EnvIn):
    return J(core.residue_environment(item, body.residue, body.radius, body.structure))


class CmpIn(BaseModel):
    a: str
    b: str
    chain_a: str | None = None
    chain_b: str | None = None


@app.post("/api/items/{item}/compare")
def compare(item: str, body: CmpIn):
    return J(core.compare_structures(item, body.a, body.b, body.chain_a, body.chain_b))


@app.get("/api/items/{item}/export")
def export(item: str, structure: str = "AF", bfactor: str = "flag", fmt: str = "cif"):
    p, name = core.export_structure(item, structure, bfactor, fmt)
    return FileResponse(p, filename=name, media_type="chemical/x-mmcif" if fmt != "pdb" else "chemical/x-pdb")


@app.get("/api/items/{item}/script")
def script(item: str, program: str = "pymol", structure: str | None = None):
    text, name = core.viewer_script(item, program, structure)
    return Response(text, media_type="text/plain", headers={"Content-Disposition": f'attachment; filename="{name}"'})


@app.get("/api/items/{item}/report.pdf")
def report_pdf(item: str):
    p = report.build_pdf(item)
    res = load_result(item)
    return FileResponse(p, filename=f"StructureBench_{(res.get('name') or item).replace(' · ', '_').replace(' ', '_')}.pdf",
                        media_type="application/pdf")


@app.post("/api/upload")
async def upload(file: UploadFile = File(...), item: str = Form("")):
    data = await file.read()
    if len(data) > 200 * 1024 * 1024:
        raise UserFacingError("That file is larger than 200 MB — too big for the viewer. Upload one model or a subset of chains.")
    return J(core.upload(file.filename, data, item or None))


# ---------------------------------------------------------------- settings, assistant
class Settings(BaseModel):
    api_key: str | None = None
    model: str | None = None


@app.get("/api/settings")
def get_settings():
    cfg = load_settings()
    key = cfg.get("api_key") or os.environ.get("ANTHROPIC_API_KEY") or ""
    return {"version": VERSION, "has_key": bool(key), "key_hint": (key[:7] + "…" + key[-4:]) if key else "",
            "model": cfg.get("model") or agent.DEFAULT_MODEL, "home": str(HOME)}


@app.post("/api/settings")
def set_settings(body: Settings):
    save_settings({k: v for k, v in body.model_dump().items() if v})
    return get_settings()


@app.get("/api/models")
def list_models():
    return agent.available_models()


@app.post("/api/quit")
def quit_server():
    threading.Timer(0.5, lambda: os._exit(0)).start()
    return {"ok": True}


class ChatIn(BaseModel):
    conv: str
    message: str
    item: str | None = None


@app.post("/api/agent/chat")
def agent_chat(body: ChatIn):
    return StreamingResponse(agent.chat(body.conv, body.message, body.item), media_type="text/event-stream",
                             headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"})


@app.post("/api/agent/reset")
def agent_reset(body: ChatIn):
    agent.reset(body.conv)
    return {"ok": True}


app.mount("/", StaticFiles(directory=ROOT / "web", html=True), name="web")
