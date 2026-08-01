#!/usr/bin/env python3
"""Verify the terminal exact data for the two Case-1 branches.

The raw degree-35 coefficient serialization is not directly equivariant under
the visible sign change.  Decode through the pinned archive's exact quintic
field implementation first, then apply the involution fixing h and negating
(u1,u2), together with the recorded row units.

The new adjacent-minor proof, replayed by the following CI step, forces every
hard-ideal solution onto h=N=0.  This checker therefore also replays the pinned
serialized Nullstellensatz certificates excluding that special fibre in both
sign branches.  It intentionally does not replay the superseded 89 MB general
h-membership certificate.
"""
from __future__ import annotations

import pickle
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
EXACT_REPLAY = REPO / (
    "plane-jc/external/zenodo-21479814/"
    "bilLkarkariy-jc2-72-108-exact-certificates-d9ea4fd/"
    "release_bundle/exact_replay"
)
sys.path.insert(0, str(EXACT_REPLAY))

from degree5_core import decode_poly, sign_substitution  # noqa: E402

BRANCH1 = EXACT_REPLAY / "hne0_polred.pkl"
BRANCH2 = EXACT_REPLAY / "hne0_branch2_polred.pkl"
ROW_SCALES = (1, 1, -1, -1, -1, -1)
SERIALIZED_CHECKER = EXACT_REPLAY / "verify_serialized_certificates.py"
SERIALIZED_MARKERS = (
    "CASE2_SERIALIZED_EXACT_PASS",
    "s=c H0_SERIALIZED_EXACT_PASS",
    "s=-c H0_SERIALIZED_EXACT_PASS",
    "ALL_SERIALIZED_EXACT_CERTIFICATES_PASS",
)

branch1 = [decode_poly(item) for item in pickle.loads(BRANCH1.read_bytes())]
branch2 = [decode_poly(item) for item in pickle.loads(BRANCH2.read_bytes())]
assert len(branch1) == len(branch2) == len(ROW_SCALES)

for index, (left, right, scale) in enumerate(
    zip(branch1, branch2, ROW_SCALES), start=1
):
    transported = sign_substitution(left, (1, -1, -1)) * scale
    assert transported.terms == right.terms, (
        f"degree-five Case-1 branch symmetry fails on residual {index}"
    )

print("SYSTEM_SYMMETRY_PASS")
print("CASE1_HARD_RESIDUAL_SYMMETRY_PASS")

completed = subprocess.run(
    [sys.executable, SERIALIZED_CHECKER.name],
    cwd=EXACT_REPLAY,
    text=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    check=False,
)
print(completed.stdout, end="")
assert completed.returncode == 0, (
    f"serialized exact certificate replay exited with {completed.returncode}"
)
lines = set(completed.stdout.splitlines())
missing = [marker for marker in SERIALIZED_MARKERS if marker not in lines]
assert not missing, f"missing serialized certificate markers: {missing}"

print("CASE1_SPECIAL_FIBRE_CERTIFICATES_PASS")
print("CASE1_BRANCH_SYMMETRY_PASS")
