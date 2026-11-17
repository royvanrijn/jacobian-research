#!/usr/bin/env python3
"""Exact replays for the GVC(3) shifted-power and minimum-frontier results."""

from __future__ import annotations

import importlib.util
import json
from math import comb, factorial
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "scripts" / "verify_gvc3_homogeneous_counterexample.py"
OUTPUT = ROOT / "artifacts" / "generated-results" / "gvc3_power_tail_and_minimum.json"

spec = importlib.util.spec_from_file_location("gvc3_base", BASE)
assert spec and spec.loader
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)


def double_factorial(value: int) -> int:
    result = 1
    for entry in range(value, 0, -2):
        result *= entry
    return result


def build_base():
    rho = base.add(base.monomial((0, 0, 2)), base.monomial((1, 1, 0)))
    a = base.add(rho, base.monomial((2, 0, 0)))
    c = base.add(
        base.multiply(base.Y, base.power(rho, 2)),
        base.scale(
            base.multiply(base.multiply(base.X, base.power(base.T, 2)), rho),
            -2,
        ),
        base.scale(base.monomial((3, 0, 2)), -1),
    )
    p6 = base.multiply(a, base.power(c, 2))
    delta = base.add(
        base.monomial((1, 1, 0), 4),
        base.monomial((0, 0, 2)),
    )
    return rho, p6, delta


def expected_tail(k: int, m: int, d: int) -> int:
    return (
        2 ** ((k + 2) * m + (k - 1) * d)
        * factorial(k * m + (k - 1) * d)
        * factorial(2 * m + 2 * d)
        * double_factorial(2 * k * m + 2 * (k + 1) * d + 1)
        // double_factorial(4 * m + 4 * d + 1)
        * comb(m + d - 1, d - 1)
    )


def expected_terminal_trace(k: int, n: int) -> int:
    return (
        2 ** ((k + 2) * n - 3)
        * factorial(k * n - 1)
        * factorial(2 * n)
        * double_factorial(2 * k * n + 3)
        // double_factorial(4 * n + 1)
    )


def verify_power_tails():
    rho, p6, delta = build_base()
    records = []
    for k in (6, 7):
        pk = base.multiply(base.power(rho, k - 6), p6)
        for d in range(1, 4 if k == 6 else 3):
            for m in range(1, 4 if k == 6 else 3):
                n = m + d
                output = base.apply_operator(
                    base.power(delta, k * m),
                    base.power(pk, n),
                )
                detector = base.multiply(
                    base.power(delta, (k - 1) * d),
                    base.monomial((0, 2 * d, 0)),
                )
                value = base.apolar_scalar(detector, output)
                expected = expected_tail(k, m, d)
                assert value == expected
                assert value != 0
                records.append(
                    {"k": k, "m": m, "d": d, "detector": str(value)}
                )
    return records


def verify_exact_trace_depth():
    _, p6, delta = build_base()
    records = []
    for n in range(1, 5):
        power = base.power(p6, n)
        pure = base.apply_operator(base.power(delta, 6 * n), power)
        assert pure == {}
        terminal = base.apply_operator(base.power(delta, 6 * n - 1), power)
        assert terminal
        value = base.apolar_scalar(base.monomial((0, 2, 0)), terminal)
        expected = expected_terminal_trace(6, n)
        assert value == expected
        records.append(
            {
                "n": n,
                "terminal_terms": len(terminal),
                "dy2_value": str(value),
            }
        )
    return records


def sphere_moment(expr, order, x, t):
    coefficient = sp.expand(expr**order).coeff(x, 0)
    return sp.factor(sp.integrate(coefficient, (t, -1, 1)) / 2)


def verify_quartic_parity_obstruction():
    x, t, a, b, c = sp.symbols("x t a b c")
    y = (1 - t**2) / x
    even = x * y - 2 * t**2 - x**2 * t**2
    odd = y - 3 * x * t**2
    h = a * x + b * y + c * t

    # Affine projective chart alpha=1: E+H*O.
    f = sp.expand(even + h * odd)
    moments = [sphere_moment(f, order, x, t) for order in range(1, 5)]
    assert sp.simplify(
        moments[0] - 2 * (5 * a - 3 * b) / 15
    ) == 0
    substituted = [
        sp.factor(value.subs(b, sp.Rational(5, 3) * a))
        for value in moments[1:]
    ]
    z = sp.symbols("z")

    def c2_to_z(expr):
        numerator, denominator = sp.together(expr).as_numer_denom()
        poly = sp.Poly(numerator, c)
        result = 0
        for (degree,), coefficient in poly.terms():
            assert degree % 2 == 0
            result += coefficient * z ** (degree // 2)
        return sp.factor(result / denominator)

    e2, e3, e4 = [c2_to_z(value) for value in substituted]
    eq2 = 28 * a**2 - 52 * a + 27 * z - 63
    assert sp.rem(
        sp.together(e2).as_numer_denom()[0], eq2, z
    ) == 0
    z_value = sp.solve(eq2, z)[0]

    g3 = sp.factor(
        sp.together(e3.subs(z, z_value)).as_numer_denom()[0]
    )
    g4 = sp.factor(
        sp.together(e4.subs(z, z_value)).as_numer_denom()[0]
    )
    g3 = sp.primitive(g3, a)[1]
    g4 = sp.primitive(g4, a)[1]
    if sp.LC(sp.Poly(g3, a)) < 0:
        g3 = -g3
    if sp.LC(sp.Poly(g4, a)) < 0:
        g4 = -g4

    expected_g3 = (
        7164 * a**3 - 36868 * a**2 - 81341 * a - 24453
    )
    expected_g4 = (
        77776 * a**4
        + 137224 * a**3
        - 745076 * a**2
        - 1119246 * a
        - 198747
    )
    assert g3 == expected_g3
    assert g4 == expected_g4
    resultant = sp.resultant(g3, g4, a)
    expected_resultant = -6466167050191094761727778002592000
    assert resultant == expected_resultant
    assert resultant != 0

    # Projective boundary alpha=0: H*O.
    boundary = sp.expand(h * odd)
    boundary_moments = [
        sphere_moment(boundary, order, x, t)
        for order in range(1, 4)
    ]
    assert sp.simplify(
        boundary_moments[0] - 2 * (5 * a - 3 * b) / 15
    ) == 0
    boundary_substituted = [
        sp.factor(value.subs(b, sp.Rational(5, 3) * a))
        for value in boundary_moments[1:]
    ]
    expected_boundary_2 = -sp.Rational(4, 315) * (
        28 * a**2 + 27 * c**2
    )
    expected_boundary_3 = sp.Rational(16, 15015) * a * (
        584 * a**2 + 819 * c**2
    )
    assert sp.simplify(
        boundary_substituted[0] - expected_boundary_2
    ) == 0
    assert sp.simplify(
        boundary_substituted[1] - expected_boundary_3
    ) == 0
    boundary_residual = sp.factor(
        (584 * a**2 + 819 * c**2).subs(
            c**2, -sp.Rational(28, 27) * a**2
        )
    )
    assert boundary_residual == -sp.Rational(796, 3) * a**2

    return {
        "affine_chart_moments": [str(value) for value in moments],
        "affine_chart_g3": str(g3),
        "affine_chart_g4": str(g4),
        "affine_chart_resultant": str(resultant),
        "projective_boundary_moments": [
            str(value) for value in boundary_moments
        ],
        "projective_boundary_after_m1": [
            str(value) for value in boundary_substituted
        ],
        "projective_boundary_residual": str(boundary_residual),
    }


def main() -> None:
    artifact = {
        "format": "gvc3-power-tail-and-minimum-v1",
        "shifted_power_checks": verify_power_tails(),
        "exact_polyharmonic_depth": verify_exact_trace_depth(),
        "quartic_linear_parity_obstruction": (
            verify_quartic_parity_obstruction()
        ),
        "conclusion": {
            "all_shifted_powers_fail": True,
            "self_multiplier_rank_one_sic": True,
            "one_profile_minimum_laplacian_power": 6,
            "common_linear_parity_completion_k2": (
                "empty through moment four over characteristic zero"
            ),
        },
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS maximal shifted-power detector")
    print("PASS exact polyharmonic depth")
    print("PASS complete common-linear parity quartic obstruction")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
