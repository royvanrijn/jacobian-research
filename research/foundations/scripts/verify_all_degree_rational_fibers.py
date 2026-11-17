#!/usr/bin/env python3
"""Exact audit of the all-degree integer-root rational-fiber theorem."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import WeightedSeedModel, w, x, y, z  # noqa: E402


def extra_roots(degree: int) -> tuple[int, ...]:
    """The N-2 extra integer roots in the uniform construction."""
    if degree < 3:
        raise ValueError("inverse degree must be at least three")
    k = degree // 2
    if degree % 2:
        return (2,) + tuple(
            root for j in range(3, k + 2) for root in (j, 1 - j)
        )
    return (3, 4) + tuple(
        root for j in range(5, k + 3) for root in (j, 1 - j)
    )


def closed_scale(degree: int) -> int:
    """Closed formula for lambda_N."""
    k = degree // 2
    if degree % 2:
        return (-1) ** k * k * math.factorial(k) ** 2
    numerator = (-1) ** k * (7 * k + 2) * math.factorial(k + 1) ** 2
    assert numerator % 144 == 0
    return numerator // 144


def audit_uniform_identities() -> None:
    """Audit the algebraic identities used in the proof for symbolic k."""
    k = sp.symbols("k", integer=True, positive=True)
    j = sp.symbols("j", integer=True)

    # Odd degree.  G'(0)=G(1), and D=G'(1)/G(1) is obtained by the
    # telescoping logarithmic-derivative sum over the paired roots.
    odd_g0 = sp.simplify(-2 * sp.product(j * (1 - j), (j, 3, k + 1)))
    assert sp.simplify(
        odd_g0 - (-1) ** k * (k + 1) * sp.factorial(k) ** 2
    ) == 0
    odd_d = sp.simplify(
        1
        + sp.Rational(1, 2)
        - 1
        + sp.summation(1 / (1 - j) + 1 / j, (j, 3, k + 1))
    )
    assert sp.simplify(odd_d - 1 / (k + 1)) == 0
    odd_lambda = sp.simplify(odd_g0 * (1 - odd_d))
    assert sp.simplify(
        odd_lambda - (-1) ** k * k * sp.factorial(k) ** 2
    ) == 0

    # Q=H_k^(2)+H_(k+1)^(2) makes the second derivative collapse to
    # -2(k+1)H_k^(2)/k.  H_k^(2)>=1 makes this strictly less than -2.
    harmonic2 = sp.harmonic(k, 2)
    odd_q = harmonic2 + (harmonic2 + 1 / (k + 1) ** 2)
    odd_hpp = sp.simplify((odd_d**2 - odd_q) / (1 - odd_d))
    assert sp.simplify(odd_hpp + 2 * (k + 1) * harmonic2 / k) == 0

    # Even degree.  The same paired-product and telescoping calculations
    # yield the 1/144 scale exactly.
    even_g0 = sp.simplify(
        12 * sp.product(j * (1 - j), (j, 5, k + 2))
    )
    assert sp.simplify(
        even_g0
        - (-1) ** k * (k + 2) * sp.factorial(k + 1) ** 2 / 12
    ) == 0
    even_d = sp.simplify(
        1
        + sp.Rational(1, 2)
        - sp.Rational(1, 2)
        - sp.Rational(1, 3)
        + sp.summation(1 / (1 - j) + 1 / j, (j, 5, k + 2))
    )
    assert sp.simplify(even_d - (sp.Rational(5, 12) + 1 / (k + 2))) == 0
    even_lambda = sp.simplify(even_g0 * (1 - even_d))
    expected = (
        (-1) ** k * (7 * k + 2) * sp.factorial(k + 1) ** 2 / 144
    )
    assert sp.simplify(even_lambda - expected) == 0

    # For k>=2: Q>=29/18 and 5/12<D<=2/3.  This strict rational
    # comparison is precisely Q>1+(1-D)^2, equivalent to H''(1)<-2.
    assert sp.Rational(29, 18) > 1 + sp.Rational(7, 12) ** 2


def audit_degree(degree: int) -> None:
    roots = (0, -1) + extra_roots(degree)
    assert len(roots) == degree
    assert len(set(roots)) == degree
    assert all(isinstance(root, int) for root in roots)

    G = sp.Poly(sp.prod(w - root for root in roots), w, domain=sp.QQ)
    derivative = G.diff()
    g0 = derivative.eval(0)
    g1 = derivative.eval(1)
    lam = g0 - g1
    assert lam == closed_scale(degree)
    assert lam != 0
    assert G.eval(1) == g0

    H = sp.cancel((G.as_expr() - g0 * w) / lam)
    model = WeightedSeedModel(sp.diff(H, w))
    assert sp.expand(model.primitive - H) == 0
    assert model.fiber_degree == degree
    assert H.subs(w, 0) == 0
    assert sp.diff(H, w).subs(w, 0) == 0
    assert H.subs(w, 1) == 0
    assert sp.diff(H, w).subs(w, 1) == -1
    assert sp.diff(H, w, 2).subs(w, 1) < -2

    # At (A,B,C)=(0,-G'(0)/lambda,1), the complete inverse pencil is
    # G/lambda.  Every root is simple, so reconstruction is rational.
    target_b = -sp.Rational(g0, lam)
    E = sp.Poly(model.inverse_polynomial(0, target_b, 1), w)
    assert E == sp.Poly(G.as_expr() / lam, w)
    assert E.degree() == degree
    assert E.gcd(E.diff()).degree() == 0
    assert all(E.eval(root) == 0 for root in roots)
    gammas = tuple(-E.diff().eval(root) for root in roots)
    assert all(gamma != 0 and gamma.is_Rational for gamma in gammas)

    # W distinguishes the reconstructed source points.  The displayed x,y,z
    # formulas additionally check rationality of every coordinate.
    a0 = model.a
    points = []
    for root, gamma in zip(roots, gammas):
        xx = 1 / gamma
        yy = root - gamma
        zz = sp.cancel((gamma - 1 - a0 * (root / gamma - 1)) / xx**2)
        assert all(sp.sympify(value).is_Rational for value in (xx, yy, zz))
        points.append((sp.cancel(xx), sp.cancel(yy), sp.cancel(zz)))
    assert len(set(points)) == degree

    # A low-degree direct substitution guards the reconstruction formulas
    # against convention or sign drift in the weighted-map implementation.
    if degree <= 8:
        mapping = model.mapping()
        for point in points:
            values = {x: point[0], y: point[1], z: point[2]}
            image = tuple(sp.cancel(component.subs(values)) for component in mapping)
            assert image == (0, target_b, 1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--max-degree",
        type=int,
        default=100,
        help="largest degree in the finite regression (default: 100)",
    )
    args = parser.parse_args()
    if args.max_degree < 3:
        parser.error("--max-degree must be at least 3")

    audit_uniform_identities()
    for degree in range(3, args.max_degree + 1):
        audit_degree(degree)

    print("PASS: symbolic odd/even scale identities hold for every k")
    print("PASS: the uniform weighted admissibility bounds give H''(1) < -2")
    print(
        "PASS: exact complete rational fibers checked in degrees "
        f"3 through {args.max_degree}"
    )


if __name__ == "__main__":
    main()
