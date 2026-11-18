#!/usr/bin/env python3
"""Exact certificate for the degree-eight alternating Laurent pencil.

Run from the repository root:
    .venv/bin/python research/arithmetic/verify_degree_eight_alternating_pencil.py
Optional: --write-receipt PATH

Requires sympy==1.14.0 (research/requirements.txt). All algebra is exact over
Q[e]/(4148928 e^2+791056 e+145375) and F_11. The script checks the algebraic
certificate; the monodromy and descent arguments are in the accompanying note.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import sympy as sp


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-receipt", type=Path)
    args = parser.parse_args()
    if sp.__version__ != "1.14.0":
        raise RuntimeError("Replay requires the pinned sympy==1.14.0")

    z, s, e, u = sp.symbols("z s e u")
    p = 4148928 * e**2 + 791056 * e + 145375
    a = -sp.Rational(625, 21952)
    h = z**8 + z**6 + e * z**2
    f = h - s * z + a
    q = (
        s**4
        + (3246152 * e - 22375) * s**2 / 19765032
        - sp.Rational(625, 3515908409649216) * (6143242504 * e + 991745875)
    )
    critical = 7 * u**4 + 5 * u**3 + e * u - a

    def red(expr: sp.Expr) -> sp.Expr:
        return sp.factor(sp.rem(sp.expand(expr), p, e))

    def norm(expr: sp.Expr) -> sp.Expr:
        value = red(expr)
        degree = sp.degree(value, e)
        return sp.factor(sp.resultant(p, value, e) / sp.LC(sp.Poly(p, e)) ** degree)

    assert sp.Poly(p, e).is_irreducible
    explicit_e = (-1009 + 94 * sp.sqrt(-329)) / 10584
    assert sp.simplify(p.subs(e, explicit_e)) == 0

    discriminant = sp.discriminant(f, z)
    assert red(discriminant + 7**7 * q**2) == 0
    critical_discriminant = red(sp.discriminant(critical, u))
    expected_critical_discriminant = -(32616905272 * e + 10887485125) / 784147392
    assert red(critical_discriminant - expected_critical_discriminant) == 0
    branch_discriminant = red(sp.discriminant(q, s))
    assert norm(a) != 0
    assert norm(critical_discriminant) != 0
    assert norm(branch_discriminant) != 0

    # The place of Q(sqrt(-329),sqrt(-7)) above 11 is w -> 1, v -> 2:
    # e=(-1009+94*w)/10584 -> 10, a -> 5, s -> 2.
    prime = 11
    e_residue, a_residue, s_residue = 10, 5, 2

    def residue(rational: sp.Expr) -> int:
        num, den = sp.fraction(rational)
        return int(num) * pow(int(den), -1, prime) % prime

    assert (-329) % prime == 1
    assert (-7) % prime == 2**2
    assert residue(sp.Rational(-1009 + 94, 10584)) == e_residue
    assert residue(a) == a_residue
    assert int(p.subs(e, e_residue)) % prime == 0
    cubic = sp.Poly(z**3 + 4*z**2 - 4*z - 3, z, modulus=prime)
    quintic = sp.Poly(z**5 - 4*z**4 - z**3 + 2*z**2 - 2*z + 2,
                      z, modulus=prime)
    reduced_f = sp.Poly(z**8 + z**6 + 10*z**2 + 5 - 2*z, z, modulus=prime)
    assert cubic * quintic == reduced_f
    assert cubic.is_irreducible
    assert quintic.is_irreducible
    assert sp.gcd(reduced_f, reduced_f.diff()).degree() == 0

    # Negative controls catch a changed discriminant certificate, a repeated
    # critical polynomial, and a reducible candidate Frobenius factor.
    assert red(discriminant + 7**7 * (q + 1)**2) != 0
    assert sp.discriminant(u**2 * (u**2 + 1), u) == 0
    assert not sp.Poly(z**3-z, z, modulus=prime).is_irreducible

    report = {
        "schema_version": 1,
        "kind": "internal_exact_replay",
        "command": ".venv/bin/python research/arithmetic/verify_degree_eight_alternating_pencil.py",
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "arithmetic": "Q[e]/(4148928e^2+791056e+145375) and F_11",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input": {
            "H": "z^8+z^6+e*z^2",
            "a": str(a),
            "coefficient_relation": str(p),
            "Q": str(q),
        },
        "checks": {
            "coefficient_polynomial_irreducible": True,
            "discriminant_equals_negative_7_power_7_times_Q_squared": True,
            "critical_polynomial_degree": 4,
            "critical_discriminant": str(critical_discriminant),
            "critical_discriminant_norm": str(norm(critical_discriminant)),
            "branch_polynomial_degree": 4,
            "branch_discriminant_norm": str(norm(branch_discriminant)),
            "a_nonzero": True,
            "mod_11_e_a_s": [e_residue, a_residue, s_residue],
            "mod_11_factor_degrees": [3, 5],
            "mod_11_factors_irreducible_and_distinct": True,
            "negative_controls": 3,
        },
        "scope": (
            "Exact algebraic identities, critical and branch squarefreeness, "
            "and the Frobenius factor certificate. Monodromy uses Jordan's "
            "theorem, simplicity of A8, and the proof in the accompanying note. "
            "No independent human review or formal verification is claimed."
        ),
        "result": "PASS",
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.write_receipt:
        args.write_receipt.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
