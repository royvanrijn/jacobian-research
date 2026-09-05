"""Source bytes for tests of frozen builders, distinct from active runtimes."""
from hashlib import sha256
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"elliptic-curves/cas"))
from pointed_quartic_migration import MIGRATED, REGRESSION_REVISION


def historical_digest(path):
    name = str(path.resolve().relative_to(ROOT))
    if name in MIGRATED:
        data = subprocess.check_output(["git", "show", REGRESSION_REVISION+":"+name], cwd=ROOT)
    else:
        data = path.read_bytes()
    return sha256(data).hexdigest()
