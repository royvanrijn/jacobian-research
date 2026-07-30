#!/usr/bin/env python3
"""Fast exact replay of the cubic-component hbar^7 unit certificate.

The generated Singular program contains the 401 consistency polynomials in
the ten effective coordinates of the reduced fifth-order lift component
over

    K = Q[a]/(94*a^3 + 335*a^2 + 400*a + 160).

It verifies that 27 residuals are already nonzero constants.  The selected
first residual is the X^18 coefficient from the earlier branchwise
certificate.  This checker pins the generated program, independently checks
the displayed Bezout inverse over Q[a], and replays both the direct one-row
identity and the full standard-basis/lift audit in Singular.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
import shutil
import subprocess

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree_five_cubic_h7_unit_certificate.sing"
)
EXPECTED_SHA256 = (
    "86eeadee714614dba8794eb392d087e3bcedcb53ce51517b0d83acda8200e980"
)


def main() -> None:
    digest = hashlib.sha256(CERTIFICATE.read_bytes()).hexdigest()
    if digest != EXPECTED_SHA256:
        raise AssertionError(
            f"stale hbar^7 certificate: {digest} != {EXPECTED_SHA256}"
        )

    a = sp.symbols("a")
    cubic = 94 * a**3 + 335 * a**2 + 400 * a + 160
    quadratic = 587583566 * a**2 + 1388701707 * a + 831388850
    bezout_quadratic = (
        32313555201 * a**2 + 79786133680 * a + 49319661920
    ) / sp.Integer(1637349242961920)
    bezout_cubic = -(
        201988446756823689 * a + 256263622855091438
    ) / sp.Integer(1637349242961920)
    assert sp.expand(
        bezout_quadratic * quadratic + bezout_cubic * cubic
    ) == 1

    residual = sp.Rational(2189187, 83886080) * quadratic
    expected_residual = (
        sp.Rational(643165152050421, 41943040) * a**2
        + sp.Rational(3040127723842209, 83886080) * a
        + sp.Rational(182006566236495, 8388608)
    )
    assert sp.expand(residual - expected_residual) == 0
    inverse = (
        sp.Rational(10771185067, 14243378945796) * a**2
        + sp.Rational(1813321220, 971139473577) * a
        + sp.Rational(12329915480, 10682534209347)
    )
    assert sp.rem(
        sp.together(residual * inverse - 1).as_numer_denom()[0],
        cubic,
        domain=sp.QQ,
    ) == 0

    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required for the hbar^7 certificate")
    result = subprocess.run(
        [singular, "-q", str(CERTIFICATE)],
        check=True,
        capture_output=True,
        text=True,
    )
    required = {
        "H7_CONSISTENCY_FIELD": "Q(a)",
        "H7_CONSISTENCY_GENERATORS": "401",
        "H7_CONSISTENCY_ACTIVE_PARAMETERS": "10/27",
        "H7_CONSTANT_RESIDUALS": "27",
        "H7_DIRECT_CERTIFICATE_GENERATOR": "1",
        "H7_DIRECT_CERTIFICATE_VERIFIED": "1",
        "H7_CONSISTENCY_GB_SIZE": "1",
        "H7_CONSISTENCY_UNIT": "1",
        "H7_NULLSTELLENSATZ_NONZERO": "1",
        "H7_NULLSTELLENSATZ_MAX_DEGREE": "0",
        "H7_NULLSTELLENSATZ_TERMS": "1",
        "H7_NULLSTELLENSATZ_VERIFIED": "1",
    }
    values = dict(
        line.split("=", 1)
        for line in result.stdout.splitlines()
        if "=" in line
    )
    for key, expected in required.items():
        if values.get(key) != expected:
            raise AssertionError(
                f"{key}: {values.get(key)!r} != {expected!r}"
            )
    selected = values.get("H7_SELECTED_CONSTANT_RESIDUAL", "")
    if "generator=1,output=(18, 0, 0)" not in selected:
        raise AssertionError(selected)

    print("PASS: 401 exact hbar^7 consistency equations use ten parameters")
    print("PASS: 27 residuals are nonzero constants on the full A^27")
    print("PASS: the selected X^18 residual has an exact Bezout inverse in K")
    print("PASS: direct and Singular-lift unit identities verify exactly")
    print("CONCLUSION: the reduced cubic fifth-order component is obstructed")


if __name__ == "__main__":
    main()
