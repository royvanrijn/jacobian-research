#!/usr/bin/env python3
"""Verify exact exceptional slices inside the HC4NHM16 squarefree-line row."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[1]
BASE_DRIVER = ROOT / "scripts" / "research_hc4_smooth_quartic_simple_line.py"
STRATA_DRIVER = ROOT / "scripts" / "research_hc4_smooth_quartic_simple_line_strata.py"
BASE_SHA256 = "99f3f94c9ee4bac0a489f25916ff290b076d33e7165e88b0a952754548c419ec"
STRATA_SHA256 = "73ddbc0416bb4a480831fb319fda5bfb9e98c336e7205c17c489568d8b3add6c"


EXPECTED: dict[str, tuple[list[str], list[str]]] = {
    "h2-zero": (
        ["w", "v", "u", "b14", "b11", "b10", "b9", "b8", "b5", "b2"],
        ["b13,", "b12,", "b7,", "b6,", "b3+3*b10-2*u", "b1+b9"],
    ),
    "h2-zero-tau-cubic": (
        ["w", "v", "u", "b14", "b11", "b10", "b9", "b8", "b5", "b4", "b2"],
        ["b13,", "b12,", "b7,", "b6,", "b3+3*b10-2*u", "3*b0+b4+2*v"],
    ),
    "tau0-delta-p-nonzero": (
        ["w", "v", "u", "b14", "b13", "b11", "b8", "b5", "b2"],
        ["b4+2*v", "b0", "b12+(2*c*m)*u+(-2*c*m^2)*v"],
    ),
    "tau0-delta-p-zero": (
        ["w", "v", "u", "b14", "b13", "b11", "b8", "b5", "b2"],
        ["b10,", "b4+2*v", "b3-2*u", "b0"],
    ),
    "tau0-delta-m3": (
        ["w", "b14", "b11", "b8", "b5", "b2"],
        ["b4+2*v", "m*b10+b9", "3*m^2*b12-80*b13+(-288*c)*u"],
    ),
    "taum1-p-equals-r": (
        ["w", "v", "b14", "b13", "b11", "b8", "b5", "b2"],
        ["u-v", "b6-b13"],
    ),
    "taum1-p-equals-r-quadratic": (
        ["w", "v", "b14", "b11", "b10", "b9", "b8", "b5", "b2", "b12^2"],
        ["u-v", "b7-b12", "b1+b9", "b0+b10"],
    ),
    "taum1-linear-factor": (
        ["w", "v", "b14", "b13", "b11", "b8", "b5", "b2"],
        ["u-v"],
    ),
    "taum1-linear-pivot-zero": (
        ["w", "v", "b14", "b12", "b11", "b8", "b5", "b2"],
        ["u-v", "11*b13+(18*p)*v"],
    ),
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_basis(output: str, begin: str, end: str) -> list[str]:
    body = output.split(begin + "\n", 1)[1].split("\n" + end, 1)[0]
    return [line.rstrip(",") for line in body.splitlines() if line]


def verify_group(group: str) -> None:
    from research_hc4_smooth_quartic_simple_line_strata import build_program

    program, equation_count, _ = build_program(group, "full")
    assert equation_count == (67 if group == "tau0-delta-m3" else 81)
    result = subprocess.run(
        ["Singular", "--no-tty", "--quiet"],
        input=program,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=300,
        check=False,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "DONE" in output and "? error" not in output, output
    expected_basis, linear_markers = EXPECTED[group]
    assert parse_basis(output, "REDUCED_BEGIN", "REDUCED_END") == expected_basis
    if group == "tau0-delta-m3":
        assert parse_basis(output, "CORE_BEGIN", "CORE_END") == ["v", "u", "b10"]
        assert parse_basis(
            output,
            "SUBSTITUTION_REMAINDER_BEGIN",
            "SUBSTITUTION_REMAINDER_END",
        ) == ["0"] * 9
    assert parse_basis(
        output, "SUPPORT_REMAINDER_BEGIN", "SUPPORT_REMAINDER_END"
    ) == ["0"] * 18
    linear = output.split("LINEAR_BEGIN\n", 1)[1].split("\nLINEAR_END", 1)[0]
    for marker in linear_markers:
        assert marker in linear, (group, marker)
    print(f"PASS {group}: exact staged basis has only determinant-zero support")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--group", choices=("all", *EXPECTED), default="all")
    parser.add_argument(
        "--audit-existing-only",
        action="store_true",
        help="verify committed equation-builder provenance without Singular replay",
    )
    args = parser.parse_args()
    assert digest(BASE_DRIVER) == BASE_SHA256
    assert digest(STRATA_DRIVER) == STRATA_SHA256
    if args.audit_existing_only:
        print(
            "PASS committed HC4 smooth-quartic exceptional-slice provenance "
            "is intact; no symbolic or Singular replay"
        )
        return
    assert shutil.which("Singular") is not None
    groups = EXPECTED if args.group == "all" else (args.group,)
    for group in groups:
        verify_group(group)
    print("THEOREM: registered squarefree-line exceptional slices are empty")


if __name__ == "__main__":
    main()
