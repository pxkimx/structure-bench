"""Incomplete API requests get the app's plain 400, not FastAPI's raw 422. Runs a real server on a free port."""
import json
import os
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


@pytest.fixture(scope="module")
def port():
    port = free_port()
    env = {**os.environ, "SB_OFFLINE": "1"}
    p = subprocess.Popen([sys.executable, "-m", "uvicorn", "server.app:app", "--port", str(port)], cwd=ROOT, env=env,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        for _ in range(200):
            try:
                urllib.request.urlopen(f"http://127.0.0.1:{port}/api/settings", timeout=1)
                break
            except Exception:  # noqa: BLE001
                time.sleep(0.1)
        else:
            raise RuntimeError("the server did not start")
        yield port
    finally:
        p.kill()
        p.wait(5)


def post(port, path, body):
    req = urllib.request.Request(f"http://127.0.0.1:{port}{path}", data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=5) as r:
            return r.status, json.loads(r.read())
    except urllib.error.HTTPError as e:
        return e.code, json.loads(e.read())


def test_missing_field_is_a_plain_400(port):
    code, body = post(port, "/api/items/none/environment", {})
    assert code == 400
    assert "residue" in body["detail"] and body["detail"].startswith("The request is incomplete")
    code, body = post(port, "/api/items/none/isoform", {})
    assert code == 400 and "entry" in body["detail"]
    code, body = post(port, "/api/items/none/environment", {"residue": "A:482", "radius": "wide"})
    assert code == 400 and "radius" in body["detail"]
