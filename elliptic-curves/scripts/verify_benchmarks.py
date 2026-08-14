#!/usr/bin/env python3
"""Verify exact numeric baselines; rank certificates have a separate replay."""

from __future__ import annotations

import json
import ast
import shutil
import subprocess
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
import sys


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
DATA = PROGRAM_ROOT / "data" / "benchmarks.json"
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.fermigier import (  # noqa: E402
    FERMIGIER_E22_RECONSTRUCTION_SHIFT,
    FERMIGIER_REPORTED_PARAMETER,
    fermigier_canonical_coefficients,
    fermigier_discriminant_factor,
    fermigier_quartic,
    twelve_visible_points,
    weierstrass_c_invariants,
    weierstrass_discriminant,
)


def gp_conductor(coefficients: list[int]) -> int:
    gp = shutil.which("gp")
    if gp is None:
        raise SystemExit("PARI/GP executable 'gp' is required")
    model = ",".join(str(coefficient) for coefficient in coefficients)
    completed = subprocess.run(
        [gp, "-q", "-f"],
        input=f"E=ellinit([{model}]);print(ellglobalred(E)[1]);\n",
        text=True,
        capture_output=True,
        check=True,
    )
    if "***" in completed.stdout + completed.stderr:
        raise SystemExit(completed.stdout + completed.stderr)
    return int(completed.stdout.strip())


def gp_minimal_model_from_quartic(
    coefficients: tuple[Fraction, ...],
) -> tuple[list[int], int, Fraction]:
    gp = shutil.which("gp")
    if gp is None:
        raise SystemExit("PARI/GP executable 'gp' is required")
    terms = []
    for degree, coefficient in enumerate(coefficients):
        terms.append(
            f"({coefficient.numerator}/{coefficient.denominator})*x^{degree}"
        )
    quartic = "+".join(terms)
    program = f"""
x='x;y='y;
r={quartic};
E=ellinit(ellfromeqn(y^2-r));
M=ellminimalmodel(E);
print("MODEL=",[M.a1,M.a2,M.a3,M.a4,M.a6]);
print("CONDUCTOR=",ellglobalred(E)[1]);
print("J=",E.j);
"""
    completed = subprocess.run(
        [gp, "-q", "-f"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    if "***" in completed.stdout + completed.stderr:
        raise SystemExit(completed.stdout + completed.stderr)
    values = dict(
        line.split("=", 1)
        for line in completed.stdout.splitlines()
        if "=" in line
    )
    return (
        ast.literal_eval(values["MODEL"]),
        int(values["CONDUCTOR"]),
        Fraction(values["J"]),
    )


def gp_minimal_model_from_weierstrass(
    coefficients: tuple[Fraction, ...],
) -> tuple[list[int], int]:
    gp = shutil.which("gp")
    if gp is None:
        raise SystemExit("PARI/GP executable 'gp' is required")
    model = ",".join(
        f"({coefficient.numerator}/{coefficient.denominator})"
        for coefficient in coefficients
    )
    program = f"""
E=ellinit([{model}]);
M=ellminimalmodel(E);
print("MODEL=",[M.a1,M.a2,M.a3,M.a4,M.a6]);
print("CONDUCTOR=",ellglobalred(E)[1]);
"""
    completed = subprocess.run(
        [gp, "-q", "-f"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    if "***" in completed.stdout + completed.stderr:
        raise SystemExit(completed.stdout + completed.stderr)
    values = dict(
        line.split("=", 1)
        for line in completed.stdout.splitlines()
        if "=" in line
    )
    return ast.literal_eval(values["MODEL"]), int(values["CONDUCTOR"])


def main() -> None:
    data = json.loads(DATA.read_text())
    for adapter_parameter in range(-12, 13):
        canonical_coefficients = fermigier_canonical_coefficients(
            adapter_parameter
        )
        assert weierstrass_discriminant(canonical_coefficients) == (
            fermigier_discriminant_factor(adapter_parameter)
        )
        if adapter_parameter == 0:
            continue
        e, d, c, b, a = fermigier_quartic(2 * adapter_parameter).quartic
        quartic_i = 12 * a * e - 3 * b * d + c * c
        quartic_j = (
            72 * a * c * e
            + 9 * b * c * d
            - 27 * a * d * d
            - 27 * b * b * e
            - 2 * c**3
        )
        c4, c6 = weierstrass_c_invariants(canonical_coefficients)
        assert quartic_i == 101232**4 * adapter_parameter**4 * c4
        assert quartic_j == 2 * 101232**6 * adapter_parameter**6 * c6
    target = data["targets"]["rank_conductor"]
    with localcontext() as context:
        context.prec = 120
        bound = Decimal(target["strict_log_conductor_bound"])
        maximum_conductor = int(bound.exp())
        assert maximum_conductor == int(target["maximum_integer_conductor"])
        e22 = data["curves"]["fermigier_e22"]
        assert e22["rank_lower_bound_reproduced_here"] is True
        assert e22["rank_lower_bound_checker"] == (
            "elliptic-curves/scripts/verify_fermigier_rank_certificates.py"
        )
        conductor = int(e22["conductor"])
        coefficients = [int(value) for value in e22["weierstrass_coefficients"]]
        assert gp_conductor(coefficients) == conductor
        literal_quartic = fermigier_quartic(FERMIGIER_REPORTED_PARAMETER)
        assert len(set(twelve_visible_points(literal_quartic))) == 12
        literal_model, literal_conductor, literal_j = gp_minimal_model_from_quartic(
            literal_quartic.quartic
        )
        assert literal_model == [
            1,
            0,
            1,
            -1223348097402005168062873899944,
            -213263015130965060475376699543914227367884158,
        ]
        assert literal_conductor == int(
            "3336936695055698757544757721801363002721636124100955091377369358240007970"
        )
        assert literal_model != coefficients

        reconstruction = fermigier_quartic(
            FERMIGIER_E22_RECONSTRUCTION_SHIFT
        )
        assert len(set(twelve_visible_points(reconstruction))) == 12
        derived_model, derived_conductor, reconstruction_j = gp_minimal_model_from_quartic(
            reconstruction.quartic
        )
        assert derived_model == coefficients
        assert derived_conductor == conductor
        assert literal_j != reconstruction_j
        canonical_model, canonical_conductor = gp_minimal_model_from_weierstrass(
            fermigier_canonical_coefficients(FERMIGIER_REPORTED_PARAMETER)
        )
        assert canonical_model == coefficients
        assert canonical_conductor == conductor
        computed_log = Decimal(conductor).ln()
        recorded_log = Decimal(e22["log_conductor"])
        assert abs(computed_log - recorded_log) < Decimal("1e-35")
        assert conductor > maximum_conductor
    # This only checks that the record model is an elliptic curve; its 29-point
    # independence certificate is external and deliberately not inferred here.
    e29_coefficients = [
        int(value)
        for value in data["curves"]["elkies_klagsbrun_e29"][
            "weierstrass_coefficients"
        ]
    ]
    a1, a2, a3, a4, a6 = e29_coefficients
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    discriminant = -b2 * b2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    assert discriminant != 0
    print(
        "PASS benchmarks: literal Fermigier shift discrepancy is pinned; the "
        "doubled-shift/canonical adapter reconstructs exact E22; rank replay "
        "is delegated to the pinned independence checker"
    )


if __name__ == "__main__":
    main()
