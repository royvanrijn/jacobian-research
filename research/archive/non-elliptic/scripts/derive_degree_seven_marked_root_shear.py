#!/usr/bin/env python3
"""Derive the degree-seven monic-slice completing shear by exact interpolation.

For each rational parameter point this replays the normalized Hamiltonian
homotopy used by ``explore_all_degree_fixed_gamma.py`` and solves the scalar
four-residue equation for the Q^2-shear.  The unsheared residue is quadratic
in the seed coefficients because the defining relative vector field is
bilinear in the two classical symbols.  Six samples therefore determine the
answer; a seventh sample is retained as an independent prediction check.

This calculation is intentionally slower than the sparse PBW verifier.  It
certifies the rank-two completion formula, not a quantum statement.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


w, X, Q, Z, s2 = sp.symbols("w X Q Z s2")
V, rho = sp.symbols("V rho")
A = -sp.Rational(8, 7)
R = 2 * X * (1 - 3 * X * Q / 2)


def hamiltonian(function):
    return sp.Matrix(
        [
            3 * X**2 * sp.diff(function, Z),
            (2 - 6 * X * Q) * sp.diff(function, Z),
            -3 * X**2 * sp.diff(function, X)
            + (6 * X * Q - 2) * sp.diff(function, Q),
        ]
    )


def completing_shear(factor):
    """Return the exact completing shear for one marked-root factor."""

    seed = sp.expand(w**2 * (w - 1) * factor)
    derivative = sp.diff(seed, w)
    complement = sp.expand(w * derivative - seed)
    assert derivative.subs(w, 1) == -1
    assert sp.diff(seed, w, 2).subs(w, 1) == -9

    W = Z + s2 * Q**2
    Y = Q - X * W / 3
    source_u = 1 - 3 * X * Y / (2 * A)
    source_gamma = 1 - 3 * X * Q / 2
    marked = sp.expand(source_u * source_gamma)
    S = sp.cancel(
        -2
        * A
        * (
            source_u
            + complement.subs(w, marked) / source_gamma**2
        )
        / (3 * X**2)
    )
    T = sp.cancel((1 + derivative.subs(w, marked) / source_gamma) / X)

    SX, SQ, SZ = (sp.diff(S, variable) for variable in (X, Q, Z))
    TX, TQ, TZ = (sp.diff(T, variable) for variable in (X, Q, Z))
    family_vector = sp.Matrix(
        [
            SZ * TQ - SQ * TZ,
            SX * TZ - SZ * TX,
            SQ * TX - SX * TQ,
        ]
    ).applyfunc(sp.cancel)
    reference_vector = sp.Matrix(
        [(1 + 3 * X * Q) / 2, -3 * Q**2, 9 * Q * Z / 2]
    )
    difference = family_vector - reference_vector

    integrand_z = sp.cancel(difference[0] / (3 * X**2))
    numerator_z, denominator_z = sp.together(integrand_z).as_numer_denom()
    zero_primitive = sum(
        coefficient
        * Z ** (monomial[0] + 1)
        / (denominator_z * (monomial[0] + 1))
        for monomial, coefficient in sp.Poly(numerator_z, Z).terms()
    )
    residual = (difference - hamiltonian(zero_primitive)).applyfunc(sp.cancel)
    assert residual[0] == residual[1] == 0

    integrand_v = sp.cancel(
        residual[2]
        .subs(X, 1 / V)
        .subs(Q, V * (2 - rho * V) / 3)
        / 3
    )
    numerator_v, denominator_v = sp.together(integrand_v).as_numer_denom()
    fixed_r_primitive = sum(
        coefficient
        * V ** (monomial[0] + 1)
        / (denominator_v * (monomial[0] + 1))
        for monomial, coefficient in sp.Poly(numerator_v, V).terms()
    )
    completion = sp.cancel(
        fixed_r_primitive.subs({V: 1 / X, rho: R}) + zero_primitive
    )
    numerator, denominator = sp.together(completion).as_numer_denom()
    pole_order = sp.Poly(denominator, X).degree()
    residue = sp.factor(sp.Poly(numerator, X).coeff_monomial(X**0))
    solutions = sp.solve(residue, s2)
    assert len(solutions) == 1
    candidate = sp.factor(solutions[0])
    polynomial_x = sp.Poly(numerator, X)
    assert all(
        sp.factor(polynomial_x.coeff_monomial(X**exponent).subs(s2, candidate))
        == 0
        for exponent in range(pole_order)
    )
    return candidate


def shear_at(sigma_value: int, tau_value: int):
    sigma = sp.Rational(sigma_value)
    tau = sp.Rational(tau_value)
    factor = (
        w**4
        + sigma * w**3
        + tau * w**2
        + (-sp.Rational(13, 2) - 3 * sigma - 2 * tau) * w
        + sp.Rational(9, 2)
        + 2 * sigma
        + tau
    )
    return completing_shear(factor)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    interpolation_points = ((0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2))
    test_point = (2, 3)
    records = {
        point: shear_at(*point) for point in (*interpolation_points, test_point)
    }

    sigma, tau = sp.symbols("sigma tau")
    coefficients = sp.symbols("c0:6")
    quadratic = (
        coefficients[0]
        + coefficients[1] * sigma
        + coefficients[2] * tau
        + coefficients[3] * sigma**2
        + coefficients[4] * sigma * tau
        + coefficients[5] * tau**2
    )
    solution = sp.solve(
        [
            quadratic.subs({sigma: point[0], tau: point[1]}) - records[point]
            for point in interpolation_points
        ],
        coefficients,
        dict=True,
    )
    assert len(solution) == 1
    candidate = sp.factor(quadratic.subs(solution[0]))
    expected = sp.factor(
        3
        * (
            26104 * sigma**2
            + 21736 * sigma * tau
            + 134160 * sigma
            + 4576 * tau**2
            + 56160 * tau
            + 75285
        )
        / 28028
    )
    assert sp.factor(candidate - expected) == 0
    assert candidate.subs({sigma: 2, tau: 3}) == records[test_point]

    certificate = {
        "scope": "exact classical rank-two shear interpolation",
        "interpolation_points": [
            {"sigma": x, "tau": y, "shear": str(records[(x, y)])}
            for x, y in interpolation_points
        ],
        "independent_test_point": {
            "sigma": test_point[0],
            "tau": test_point[1],
            "shear": str(records[test_point]),
        },
        "quadratic_shear": str(candidate),
    }
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.output:
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered)

    print("PASS: six exact homotopy samples determine the quadratic shear")
    print("PASS: the seventh exact sample is predicted correctly")
    print(f"shear = {candidate}")


if __name__ == "__main__":
    main()
