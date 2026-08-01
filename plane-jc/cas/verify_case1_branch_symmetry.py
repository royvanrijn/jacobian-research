#!/usr/bin/env python3
"""Run the pinned exact replay for the Case-1 sign involution.

The raw degree-35 coefficient serialization is not itself equivariant under the
visible sign change; equivariance appears after the archive's exact descent to
the intrinsic quintic field.  Keep that descent in one authoritative place and
turn its four markers into a stable repository-level regression marker.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
EXACT_REPLAY = REPO / (
    "plane-jc/external/zenodo-21479814/"
    "bilLkarkariy-jc2-72-108-exact-certificates-d9ea4fd/"
    "release_bundle/exact_replay"
)
CHECKER = EXACT_REPLAY / "verify_hne0_branch_symmetry.py"
EXPECTED = (
    "SYSTEM_SYMMETRY_PASS",
    "BRANCH1_EXACT_IDENTITY_PASS",
    "BRANCH2_EXACT_IDENTITY_PASS",
    "GMPY2_EXACT_PASS",
)

assert CHECKER.is_file(), f"missing pinned symmetry checker: {CHECKER}"
completed = subprocess.run(
    [sys.executable, CHECKER.name],
    cwd=EXACT_REPLAY,
    text=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    check=False,
)
print(completed.stdout, end="")
assert completed.returncode == 0, (
    f"pinned symmetry replay exited with status {completed.returncode}"
)
lines = set(completed.stdout.splitlines())
missing = [marker for marker in EXPECTED if marker not in lines]
assert not missing, f"missing exact symmetry markers: {missing}"
print("CASE1_BRANCH_SYMMETRY_PASS")
