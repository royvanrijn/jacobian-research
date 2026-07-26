#!/usr/bin/env python3
"""Reject `sorryAx` and require axiom reports for the public Lean certificates."""

from __future__ import annotations

import sys
from pathlib import Path

REQUIRED = (
    "FiniteEtaleKeller.jacobianDet_generalGaugeMap",
    "FiniteEtaleKeller.jacobianDet_generalGaugeJacobianOneMap",
    "FiniteEtaleKeller.generalGaugeInversePolynomial_derivative",
    "FiniteEtaleKeller.LocalizedPolynomialRoot.localizedAlgHomEquiv",
    "FiniteEtaleKeller.generalGaugeMap_announcedSeed",
    "FiniteEtaleKeller.automaticRealizationMap_certificate",
    "FiniteEtaleKeller.automaticJacobianOneFiberRepresentingEquiv_natural",
)


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: check_lean_axiom_log.py BUILD_LOG", file=sys.stderr)
        return 2

    path = Path(sys.argv[1])
    text = path.read_text(encoding="utf-8", errors="replace")

    failures: list[str] = []
    if "sorryAx" in text:
        failures.append("Lean build log contains `sorryAx`")

    for declaration in REQUIRED:
        if declaration not in text:
            failures.append(f"missing #print axioms report for {declaration}")

    if failures:
        print("Lean axiom audit failed:", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1

    print("PASS: public certificate axiom reports are present and contain no sorryAx")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
