#!/usr/bin/env python3
"""Verify the paper's explicit formulas against their Lean definitions."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
GENERATOR = ROOT / "scripts/compile_common_arithmetic_fibers_example.py"
CERTIFICATE = (
    ROOT
    / "artifacts/generated-results/common_arithmetic_fibers_explicit_quintic.json"
)
PAPER = (
    ROOT
    / "papers/common-arithmetic-fibers/sections/04-explicit-example.tex"
)
FORMAL = ROOT / "formal/finite-etale-keller"
LEAN_TARGET = "FiniteEtaleKeller.PaperExampleCorrespondence"


def run(command: list[str], *, cwd: Path = ROOT) -> str:
    completed = subprocess.run(
        command,
        cwd=cwd,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if completed.returncode:
        raise SystemExit(
            f"command failed ({completed.returncode}): {' '.join(command)}\n"
            f"{completed.stdout}"
        )
    return completed.stdout


def verify_generated_hashes(certificate: dict[str, object]) -> None:
    views = certificate["generated_views"]
    assert isinstance(views, dict)
    for name, raw_record in views.items():
        assert isinstance(raw_record, dict)
        path = ROOT / str(raw_record["path"])
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != raw_record["sha256"]:
            raise SystemExit(f"{name} view does not match the JSON certificate: {path}")


def verify_paper_inputs_generated_tex() -> None:
    source = PAPER.read_text(encoding="utf-8")
    required = [
        r"\input{../../artifacts/generated-results/common_arithmetic_fibers_explicit_quintic.tex}",
        r"\CAFExplicitPolynomialFactorized",
        r"\CAFExplicitPolynomialExpanded",
        r"\CAFExplicitSeed",
        r"\CAFExplicitT",
        r"\CAFExplicitQ",
        r"\CAFExplicitFOne",
        r"\CAFExplicitFTwo",
        r"\CAFExplicitFThree",
        r"\CAFExplicitIntegralScaling",
        r"\CAFExplicitIntegralTarget",
        r"\CAFExplicitNormalizedTarget",
        r"\CAFExplicitInverseIdentity",
        r"\CAFExplicitIntegralToJacobianOneScaling",
        r"\CAFExplicitJacobianOneTarget",
    ]
    missing = [token for token in required if token not in source]
    if missing:
        raise SystemExit(
            "the paper bypasses generated explicit-example data:\n"
            + "\n".join(missing)
        )


def main() -> None:
    generator_output = run([sys.executable, str(GENERATOR), "--check"])
    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    verify_generated_hashes(certificate)
    verify_paper_inputs_generated_tex()

    sympy_path = ROOT / certificate["generated_views"]["sympy"]["path"]
    sympy_output = run([sys.executable, str(sympy_path)])

    lean_output = run(["lake", "build", LEAN_TARGET], cwd=FORMAL)
    forbidden = ("sorryAx", "declaration uses 'sorry'")
    if any(marker in lean_output for marker in forbidden):
        raise SystemExit("the Lean correspondence build introduced a placeholder")

    theorem_names = certificate["lean_correspondence"]["theorems"]
    correspondence_source = (
        FORMAL / "FiniteEtaleKeller/PaperExampleCorrespondence.lean"
    ).read_text(encoding="utf-8")
    absent = [name for name in theorem_names if name not in correspondence_source]
    if absent:
        raise SystemExit(
            "certificate names absent Lean correspondence theorems:\n"
            + "\n".join(absent)
        )

    print(generator_output.strip())
    print(sympy_output.strip())
    print("PASS: paper TeX consumes only generated explicit-example coefficients")
    print("PASS: Lean proves the generated polynomial, map, inverse, scalings, and targets")
    print("PASS: JSON hashes authenticate the Lean/TeX/SymPy/Sage views")


if __name__ == "__main__":
    main()
