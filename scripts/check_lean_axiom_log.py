#!/usr/bin/env python3
"""Reject `sorryAx` and require reports for each publication certificate."""

from __future__ import annotations

import sys
from pathlib import Path

FIBERS_REQUIRED = (
    "FiniteEtaleKeller.generalGaugeFunctionFieldComparison",
    "FiniteEtaleKeller.generalGaugeGeometricDegree_eq",
    "FiniteEtaleKeller.automaticRealizationFunctionFieldComparison",
    "FiniteEtaleKeller.automaticRealizationFunctionField_finrank",
    "FiniteEtaleKeller.automaticRealization_pageOne",
    "FiniteEtaleKeller.paperPolynomialPresentation_pageOne",
    "FiniteEtaleKeller.finiteEtalePresentation",
    "FiniteEtaleKeller.abstractFiniteEtale_pageOne",
    "FiniteEtaleKeller.paperAbstractFiniteEtale_pageOne",
)

HASSE_REQUIRED = (
    "FiniteEtaleKeller.FixedHasseFamily.paperMap_normalization_inverse",
    "FiniteEtaleKeller.FixedHasseFamily.paperMap_geometricDegree",
    "FiniteEtaleKeller.FixedHasseFamily.jacobianDet_paperMap",
    "FiniteEtaleKeller.FixedHasseFamily.paperFiberRepresentingEquiv",
    "FiniteEtaleKeller.FixedHasseFamily.paperFiberPoint_hasse_certificate",
    "FiniteEtaleKeller.FixedHasseFamily.targetProjectiveContent_eq_one",
    "FiniteEtaleKeller.FixedHasseFamily.targetProjectiveHeight_eq",
    "FiniteEtaleKeller.FixedHasseFamily.paperParameter_certificate",
    "FiniteEtaleKeller.FixedHasseFamily.fixedHassePaper_certificate",
    "FiniteEtaleKeller.FixedHasseFamily.HasseCoreCondition.mul",
    "FiniteEtaleKeller.FixedHasseFamily.prime_not_rational_cube",
    "FiniteEtaleKeller.FixedHasseFamily.primeParameter_certificate",
)

CERTIFICATES = {
    "fibers": FIBERS_REQUIRED,
    "fixed-hasse": HASSE_REQUIRED,
}


def main() -> int:
    if len(sys.argv) == 2:
        certificate = "fibers"
        path_arg = sys.argv[1]
    elif len(sys.argv) == 3 and sys.argv[1] in CERTIFICATES:
        certificate = sys.argv[1]
        path_arg = sys.argv[2]
    else:
        print(
            "usage: check_lean_axiom_log.py "
            "[fibers|fixed-hasse] BUILD_LOG",
            file=sys.stderr,
        )
        return 2

    path = Path(path_arg)
    text = path.read_text(encoding="utf-8", errors="replace")

    failures: list[str] = []
    if "sorryAx" in text:
        failures.append("Lean build log contains `sorryAx`")

    for declaration in CERTIFICATES[certificate]:
        if declaration not in text:
            failures.append(f"missing #print axioms report for {declaration}")

    if failures:
        print("Lean axiom audit failed:", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1

    print(
        f"PASS: {certificate} certificate axiom reports are present "
        "and contain no sorryAx"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
