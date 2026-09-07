#!/usr/bin/env python3
"""Hash-pinned theorem and independent replay for the F13-to-F12 reduction."""

from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PINNED_FILES = {
    "scripts/verify_macfarlane_f12_reduction.py": (
        "60e1cfdf3132763d4285a296873bf82f6deacc807316600c4a60597c12a1c6ff"
    ),
    "scripts/audit_macfarlane_f12_independent.py": (
        "3a2e85cc222c5d3318b5479a7c400c4f3f298460d0ea6e916e1de592d778abf2"
    ),
    "verified/TWELVE_VARIABLE_DEGREE_THREE_KELLER_COUNTEREXAMPLE.md": (
        "66ca76270d8bfba8cc0279eaf2411958151027ed0f37f5e3bc4124fd627c1c71"
    ),
    "scripts/audit_macfarlane_g20_dimension_reduction.py": (
        "d29ad30d8f4d14890dc21f24c70ced9a4cee2735cb05d8902abe53370518a404"
    ),
}
CHECKERS = (
    "scripts/verify_macfarlane_f12_reduction.py",
    "scripts/audit_macfarlane_f12_independent.py",
)


def main() -> None:
    assert PINNED_FILES, "F12 suite pins have not been finalized"
    for filename, expected in PINNED_FILES.items():
        actual = hashlib.sha256((ROOT / filename).read_bytes()).hexdigest()
        if actual != expected:
            raise AssertionError(f"stale F12 suite pin for {filename}: {actual}")
    for checker in CHECKERS:
        subprocess.run([sys.executable, checker], cwd=ROOT, check=True)
    print("PASS F12 suite: exact reduction and independent determinant replay")


if __name__ == "__main__":
    main()
