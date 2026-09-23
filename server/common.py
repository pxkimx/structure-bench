"""Shared paths, workspace items, settings and error helpers.

Every analysis the user runs is an *item*: a folder under SB_HOME/work/<id>/ holding request.json (what
was asked), result.json (what the UI, the PDF report and the assistant all read), the structure files it
used and any variant mapping or comparison done on it. Nothing in here imports FastAPI.
"""
from __future__ import annotations

import json
import os
import time
import traceback
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HOME = Path(os.environ.get("SB_HOME", Path.home() / "Library" / "Application Support" / "StructureBench"))
WORK = HOME / "work"
CACHE = HOME / "cache"          # every download, so a lookup done once also works offline
EXAMPLES = ROOT / "examples"
for _d in (WORK, CACHE):
    _d.mkdir(parents=True, exist_ok=True)
VERSION = (ROOT / "VERSION").read_text().strip() if (ROOT / "VERSION").exists() else "dev"


class UserFacingError(Exception):
    """An error whose message is written for the user (bad input, unknown gene) — shown without a traceback."""


def log_exc(where: str) -> str:
    """Print the current exception's traceback to the server log and return a short 'Type: msg' string."""
    import sys
    et, ev, _ = sys.exc_info()
    print(f"[{where}] non-fatal error:\n" + traceback.format_exc(), flush=True)
    return f"{et.__name__ if et else 'Error'}: {str(ev)[:160]}"


def versions(*pkgs) -> dict:
    from importlib.metadata import version
    out = {}
    for k in pkgs:
        try:
            out[k] = version(k)
        except Exception:  # noqa: BLE001
            pass
    return out


def json_default(o):
    try:
        import numpy as np
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return None if not np.isfinite(o) else float(o)
        if isinstance(o, np.bool_):
            return bool(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
    except Exception:  # noqa: BLE001
        pass
    if isinstance(o, Path):
        return str(o)
    return str(o)


def dump(obj, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, default=json_default))
    tmp.replace(path)          # never leave a half-written result.json for the UI to trip over


# ---------------------------------------------------------------- workspace items
def new_item_dir() -> Path:
    d = WORK / uuid.uuid4().hex[:10]
    (d / "structures").mkdir(parents=True)
    return d


def item_dir(item: str) -> Path:
    item = str(item or "").strip()
    p = (WORK / item).resolve()
    if not item or not str(p).startswith(str(WORK.resolve())) or not (p / "result.json").exists():
        raise UserFacingError(f"There is no saved analysis called '{item}'. Open one from Recent on the home page.")
    return p


def load_result(item: str) -> dict:
    return json.loads((item_dir(item) / "result.json").read_text())


def save_result(item: str, res: dict):
    dump(res, item_dir(item) / "result.json")


def list_items(limit: int = 60) -> list[dict]:
    rows = []
    if not WORK.exists():
        return rows
    for d in sorted(WORK.iterdir(), key=lambda x: x.stat().st_mtime, reverse=True):
        if not d.is_dir() or not (d / "result.json").exists():
            continue
        try:
            r = json.loads((d / "result.json").read_text())
            rows.append({"item": d.name, "kind": r.get("kind"), "name": r.get("name"), "t": d.stat().st_mtime * 1000,
                         "summary": r.get("summary") or {}})
        except Exception:  # noqa: BLE001 - one broken folder must not hide the rest
            continue
        if len(rows) >= limit:
            break
    return rows


def touch(item: str):
    try:
        os.utime(item_dir(item), None)
    except Exception:  # noqa: BLE001
        pass


def now_iso() -> str:
    return time.strftime("%Y-%m-%d %H:%M")


# ---------------------------------------------------------------- settings (API key, model)
CONFIG = Path.home() / ".structure-bench" / "config.json"


def load_settings() -> dict:
    """Settings saved from the UI. The API key falls back to RNAseq Bench's, so one key serves every Bench app."""
    cfg = {}
    try:
        cfg = json.loads(CONFIG.read_text())
    except Exception:  # noqa: BLE001
        pass
    if not cfg.get("api_key"):
        try:
            k = json.loads((Path.home() / ".rnaseq-bench" / "config.json").read_text()).get("api_key")
            if k:
                cfg["api_key"] = k
        except Exception:  # noqa: BLE001
            pass
    return cfg


def save_settings(d: dict):
    CONFIG.parent.mkdir(parents=True, exist_ok=True)
    cur = {}
    try:
        cur = json.loads(CONFIG.read_text())
    except Exception:  # noqa: BLE001
        pass
    cur.update(d)
    CONFIG.write_text(json.dumps(cur))
    try:
        CONFIG.chmod(0o600)
    except Exception:  # noqa: BLE001
        pass
