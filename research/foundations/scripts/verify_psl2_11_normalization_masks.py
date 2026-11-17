#!/usr/bin/env python3
"""Robust driver for the exact PSL2(11) normalization-mask replay."""

from hashlib import sha256
from pathlib import Path
import shutil
import subprocess


MASK_SOURCE_SHA256 = (
    "18db2e56f09d9d0d38ca027ec1d88511855b783eb7ce25bcf94c5fc3db5d0695"
)
REQUIRED_MARKERS = (
    "PASS C5 explicit normalization masks M1--M4",
    "PASS C6 filtered Riemann--Roch dimensions 23, 29, and 41",
    "PASS C6 uniform-pole masks N2--N5",
    "PASS C6 exchanged masks N6--N7 with forced pole pair (6,7)",
    "PASS C6 asymmetric mask N1 with pole pair (7,5)",
    "PASS explicit normalization representatives for all C5/C6 masks",
)


def run_normalization_mask_check(singular=None):
    """Run Singular and reject silent procedure failures or stale source."""

    if singular is None:
        singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for the normalization-mask replay")

    mask_script = Path(__file__).with_suffix(".sing")
    actual_hash = sha256(mask_script.read_bytes()).hexdigest()
    if actual_hash != MASK_SOURCE_SHA256:
        raise RuntimeError(
            "stale normalization-mask source hash: "
            f"expected {MASK_SOURCE_SHA256}, got {actual_hash}"
        )

    result = subprocess.run(
        [singular, "-q", str(mask_script)],
        text=True,
        capture_output=True,
        check=False,
    )
    output = result.stdout + result.stderr
    if (
        result.returncode != 0
        or "FAIL" in output
        or "   ?" in output
        or any(marker not in output for marker in REQUIRED_MARKERS)
    ):
        raise RuntimeError("Singular normalization-mask replay failed:\n" + output)

    for line in result.stdout.splitlines():
        if line.startswith("PASS"):
            print(line)
    print("PASS robust Singular failure and source-hash audit")


if __name__ == "__main__":
    run_normalization_mask_check()
