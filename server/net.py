"""Every download goes through here: a 20 s timeout, an on-disk cache, and a bundled-example fallback.

Rules (the app's internet rule):
- cache first — a lookup done once works offline afterwards, and never waits on the network again;
- a failed download raises NetError with a sentence a biologist can act on, which callers turn into an
  `info` flag rather than an exception on screen;
- when the network fails and the file is one of the bundled examples (LMNA), the example copy is used.

Biopython's own fetchers (Bio.UniProt.search, Bio.PDB.alphafold_db.get_predictions) call urlopen without a
timeout. They look `urlopen` up in their own module namespace, so it is replaced there with one that has
the 20 s timeout — the Biopython code path is kept, it just cannot hang the app.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import time
import urllib.error
import urllib.request
from pathlib import Path

from .common import CACHE, EXAMPLES

TIMEOUT = 20
UA = "StructureBench/0.1 (local desktop app; Biopython)"


class NetError(Exception):
    """A download that failed. The message is written for the user."""


def offline() -> bool:
    """SB_OFFLINE=1 forces the no-network path (tests, self-test)."""
    return os.environ.get("SB_OFFLINE", "") not in ("", "0", "false")


def _urlopen_with_timeout(url, *a, **kw):
    if offline():
        raise urllib.error.URLError("offline mode (SB_OFFLINE=1)")
    kw.setdefault("timeout", TIMEOUT)
    if isinstance(url, str):
        url = urllib.request.Request(url, headers={"User-Agent": UA})
    return urllib.request.urlopen(url, *a, **kw)


def patch_biopython():
    try:
        import Bio.UniProt as _u
        _u.urlopen = _urlopen_with_timeout
    except Exception:  # noqa: BLE001
        pass
    try:
        import Bio.PDB.alphafold_db as _a
        _a.urlopen = _urlopen_with_timeout
    except Exception:  # noqa: BLE001
        pass


patch_biopython()


def _seed_path(rel: str) -> Path | None:
    """A bundled example file standing in for this cache entry, if there is one (examples/*/manifest.json)."""
    for man in EXAMPLES.glob("*/manifest.json"):
        try:
            m = json.loads(man.read_text())
        except Exception:  # noqa: BLE001
            continue
        f = m.get("cache", {}).get(rel)
        if f and (man.parent / f).exists():
            return man.parent / f
    return None


def cache_path(rel: str) -> Path:
    p = (CACHE / rel).resolve()
    if not str(p).startswith(str(CACHE.resolve())):
        raise ValueError("bad cache path")
    return p


def fetch(url: str, rel: str, what: str, refresh: bool = False) -> tuple[bytes, str]:
    """Return (content, source) where source is 'cache' | 'network' | 'example'.

    what — a short phrase for messages ("the UniProt entry for P02545").
    """
    p = cache_path(rel)
    if p.exists() and p.stat().st_size > 0 and not refresh:
        return p.read_bytes(), "cache"
    err = None
    if not offline():
        t0 = time.time()
        try:
            with _urlopen_with_timeout(url) as r:
                data = r.read()
            p.parent.mkdir(parents=True, exist_ok=True)
            tmp = p.with_suffix(p.suffix + ".part")
            tmp.write_bytes(data)
            tmp.replace(p)
            print(f"[net] {url} {len(data)//1024} KB in {time.time()-t0:.1f}s", flush=True)
            return data, "network"
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise NetError(f"{what[0].upper() + what[1:]} does not exist (HTTP 404 from {url.split('/')[2]}).") from None
            err = f"HTTP {e.code}"
        except Exception as e:  # noqa: BLE001 - timeouts, DNS, TLS, reset: all mean "not reachable now"
            err = f"{type(e).__name__}: {str(e)[:120]}"
    else:
        err = "offline mode"
    if p.exists() and p.stat().st_size > 0:          # refresh failed: the cached copy is still better than nothing
        return p.read_bytes(), "cache"
    seed = _seed_path(rel)
    if seed:
        p.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(seed, p)
        return p.read_bytes(), "example"
    host = url.split("/")[2] if "//" in url else url
    raise NetError(f"{host} could not be reached ({err}), and {what} is not in the local cache yet.")


def fetch_json(url: str, rel: str, what: str, refresh: bool = False):
    data, src = fetch(url, rel, what, refresh)
    try:
        return json.loads(data), src
    except Exception:
        # a truncated or HTML error page in the cache must not poison every later lookup
        cache_path(rel).unlink(missing_ok=True)
        raise NetError(f"{what[0].upper() + what[1:]} came back unreadable; try again.") from None


def cached_call(rel: str, fn, what: str, refresh: bool = False):
    """Cache the JSON result of a Biopython fetcher (e.g. a UniProt search) the same way as a download."""
    p = cache_path(rel)
    if p.exists() and not refresh:
        return json.loads(p.read_text()), "cache"
    err = None
    if not offline():
        try:
            out = fn()
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(json.dumps(out))
            return out, "network"
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise NetError(f"{what[0].upper() + what[1:]} does not exist (HTTP 404).") from None
            err = f"HTTP {e.code}"
        except Exception as e:  # noqa: BLE001
            err = f"{type(e).__name__}: {str(e)[:120]}"
    else:
        err = "offline mode"
    if p.exists():
        return json.loads(p.read_text()), "cache"
    seed = _seed_path(rel)
    if seed:
        p.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(seed, p)
        return json.loads(p.read_text()), "example"
    raise NetError(f"UniProt could not be reached ({err}), and {what} is not in the local cache yet.")


def key(*parts) -> str:
    return hashlib.sha1("|".join(map(str, parts)).encode()).hexdigest()[:16]
