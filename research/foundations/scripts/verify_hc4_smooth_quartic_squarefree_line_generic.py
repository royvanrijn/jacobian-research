#!/usr/bin/env python3
"""Verify the generic squarefree-line gate following HC4NHM14."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "research_hc4_smooth_quartic_simple_line.py"
DRIVER_SHA256 = "99f3f94c9ee4bac0a489f25916ff290b076d33e7165e88b0a952754548c419ec"


def verify_driver_hash() -> None:
    digest = hashlib.sha256(DRIVER.read_bytes()).hexdigest()
    assert digest == DRIVER_SHA256, (digest, DRIVER_SHA256)
    print(f"PASS pinned equation driver sha256:{digest}")


def verify_generic_standard_basis() -> None:
    assert shutil.which("Singular") is not None
    result = subprocess.run(
        [
            sys.executable,
            str(DRIVER),
            "--case",
            "squarefree-line",
            "--triangular",
        ],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=240,
        check=False,
    )
    if result.returncode != 0:
        raise AssertionError(result.stdout + result.stderr)
    output = result.stdout
    for marker in (
        "equation_count=81",
        "linear_basis_size=10",
        "reduced_generator_count=71",
        "nonlinear_basis_size=8",
        "basis_size=18",
        "dimension=0",
    ):
        assert marker in output, marker

    basis_text = output.split("BASIS_BEGIN\n", 1)[1].split("\nBASIS_END", 1)[0]
    basis = [entry.rstrip(",") for entry in basis_text.splitlines() if entry]
    assert basis == [
        "w",
        "v",
        "u",
        "b14",
        "b13",
        "b12",
        "b11",
        "b10",
        "b9",
        "b8",
        "b7",
        "b6",
        "b5",
        "b4",
        "b3",
        "b2",
        "b1",
        "b0",
    ]
    print("PASS generic coefficient-field basis is the 18-variable maximal ideal")


def verify_degenerate_remainder() -> None:
    x, y, z = sp.symbols("x y z")
    p, q, r, c0, c1, c2 = sp.symbols("p q r c0 c1 c2")
    matrix = sp.Matrix(
        [
            [0, 0, -y**2],
            [0, 0, x**2],
            [
                -y**2,
                x**2,
                p * x**2 + q * x * y + r * y**2 + z * (c0 * x + c1 * y + c2 * z),
            ],
        ]
    )
    assert sp.expand(matrix.det()) == 0
    print("PASS the unique generic relaxed fiber has det(A)=0")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--audit-existing-only",
        action="store_true",
        help="verify committed equation-builder provenance without Singular replay",
    )
    arguments = parser.parse_args()

    verify_driver_hash()
    if arguments.audit_existing_only:
        print(
            "PASS committed HC4 smooth-quartic generic-line provenance is intact; "
            "no symbolic or Singular replay"
        )
        return
    verify_generic_standard_basis()
    verify_degenerate_remainder()
    print("THEOREM: generic squarefree-line simple-residual-line stratum is empty")


if __name__ == "__main__":
    main()
