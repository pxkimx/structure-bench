"""Every test runs offline, on the bundled LMNA example, in a throw-away SB_HOME."""
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
os.environ["SB_OFFLINE"] = "1"
os.environ["SB_HOME"] = tempfile.mkdtemp(prefix="sb_test_")

EX = ROOT / "examples" / "lmna"
MODEL = EX / "AF-P02545-F1-model_v6.cif"
PAE = EX / "AF-P02545-F1-predicted_aligned_error_v6.json"
AM = EX / "AF-P02545-F1-aa-substitutions.csv"
UNIPROT = EX / "uniprot_P02545.json"
