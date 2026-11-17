#!/usr/bin/env python3
"""Verify the complete compact Case-1 determinantal closure.

The default mode combines the new adjacent-minor characteristic-zero
standard-basis decision with the two archived special-fibre certificates and
the exact sign-branch transport.  ``--audit-existing-only`` performs only
committed-byte and compact reconstruction checks; it never starts Singular or
multiplies the large archived certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPLAY = (
    ROOT
    / "plane-jc/external/zenodo-21479814/"
    "bilLkarkariy-jc2-72-108-exact-certificates-d9ea4fd/"
    "release_bundle/exact_replay"
)
sys.path.insert(0, str(ROOT / "scripts"))

import research_jc72_108_case1_rankdrop as rankdrop  # noqa: E402


ADDITIONAL_INPUT_SHA256 = {
    "hne0_polred.pkl": (
        "5a6e423d74ef09fc9c7a7282c500bda566018d7e56a93124665796bbe417cedf"
    ),
    "hne0_branch2_polred.pkl": (
        "f108089ec3abc4714bba5986e8c881f4a83f72f38480e9e1548fd803aafaf717"
    ),
    "degree5_core.py": (
        "90bc904196e98052cc40fe22f92ef596db8a6954233e3b52fe28cf195d090aa9"
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def verify_manifest(filename: str) -> int:
    """Validate one archive checksum manifest and reject ambiguous paths."""

    seen: set[str] = set()
    count = 0
    manifest = REPLAY / filename
    for line_number, line in enumerate(manifest.read_text().splitlines(), 1):
        if not line:
            continue
        try:
            expected, relative = line.split(maxsplit=1)
        except ValueError as error:
            raise RuntimeError(f"malformed {filename} line {line_number}") from error
        relative = relative.removeprefix("*")
        candidate = Path(relative)
        if candidate.is_absolute() or ".." in candidate.parts:
            raise RuntimeError(f"unsafe path in {filename}: {relative}")
        if relative in seen:
            raise RuntimeError(f"duplicate path in {filename}: {relative}")
        seen.add(relative)
        actual = sha256(REPLAY / candidate)
        if actual != expected:
            raise RuntimeError(
                f"{filename}: {relative} changed: expected {expected}, got {actual}"
            )
        count += 1
    if not count:
        raise RuntimeError(f"empty checksum manifest: {filename}")
    return count


def audit_existing_only() -> None:
    exact_count = verify_manifest("EXACT_SHA256SUMS.txt")
    reconstructed_count = verify_manifest("RECONSTRUCTED_CERTIFICATES.sha256")
    for relative, expected in ADDITIONAL_INPUT_SHA256.items():
        actual = sha256(REPLAY / relative)
        if actual != expected:
            raise RuntimeError(
                f"unmanifested determinantal input changed: {relative}: "
                f"expected {expected}, got {actual}"
            )
    rankdrop.audit_existing_only()
    print(
        "JC2_72_108_CASE1_DETERMINANTAL_COMMITTED_AUDIT_PASS "
        f"(manifested files {exact_count}+{reconstructed_count}; "
        "3 additional sign-transport inputs; no Singular or large multiplication)"
    )


def run_checked(command: list[str], *, cwd: Path, marker: str) -> None:
    completed = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode or marker not in completed.stdout:
        raise RuntimeError(
            f"replay failed or omitted {marker}:\n"
            + completed.stdout[-4000:]
            + completed.stderr[-4000:]
        )
    print(completed.stdout, end="")


def full_replay(singular_timeout: int) -> None:
    audit_existing_only()
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for the exact determinantal replay")
    with tempfile.TemporaryDirectory(prefix="jc72-108-rankdrop-") as directory:
        singular_input = Path(directory) / "case1-rankdrop.sing"
        rankdrop.write_singular(singular_input)
        completed = subprocess.run(
            [singular, "-q", str(singular_input)],
            text=True,
            capture_output=True,
            timeout=singular_timeout,
            check=False,
        )
        required = (
            "CASE1_SPECIAL_FIBRE_BEZOUT_PASS",
            "CASE1_N_MOD_H_PASS",
            "RANKDROP_EXACT_UNIT_PASS",
        )
        if completed.returncode or any(
            marker not in completed.stdout for marker in required
        ):
            raise RuntimeError(
                "exact adjacent-minor replay failed:\n"
                + completed.stdout[-4000:]
                + completed.stderr[-4000:]
            )
        print(completed.stdout, end="")

    run_checked(
        [sys.executable, str(REPLAY / "verify_serialized_certificates.py")],
        cwd=REPLAY,
        marker="ALL_SERIALIZED_EXACT_CERTIFICATES_PASS",
    )
    run_checked(
        [sys.executable, str(REPLAY / "verify_hne0_branch_symmetry.py")],
        cwd=REPLAY,
        marker="BRANCH2_EXACT_IDENTITY_PASS",
    )
    print("JC2_72_108_CASE1_DETERMINANTAL_CLOSURE_PASS")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-existing-only", action="store_true")
    parser.add_argument("--singular-timeout", type=int, default=600)
    args = parser.parse_args()
    if args.singular_timeout <= 0:
        parser.error("--singular-timeout must be positive")
    if args.audit_existing_only:
        audit_existing_only()
    else:
        full_replay(args.singular_timeout)


if __name__ == "__main__":
    main()
