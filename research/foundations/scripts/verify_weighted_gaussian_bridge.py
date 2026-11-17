#!/usr/bin/env python3
"""Exact audit of the weighted-seed to Gaussian-moment bridge.

For h(0)=1, the construction in
extended-geometry/WEIGHTED_GAUSSIAN_BRIDGE.md produces a polynomial map
Phi=(phi_1,phi_2).  Its fixed branch has determinant correction one, so
P=W_1 phi_1(Z_1,Z_2)+W_2 phi_2(Z_1,Z_2) has all pure Gaussian moments zero,
while Q=Z_1 reads the nonpolynomial branch g=u h(g).

The general proof is written in the companion note.  This script checks the
universal algebra symbolically and performs bounded Wick/coefficient audits
for canonical and split weighted seeds.
"""

from __future__ import annotations

from fractions import Fraction
from math import factorial

import sympy as sp


z, y, u = sp.symbols("z y u")


def bridge(seed: sp.Expr, scale: sp.Expr = sp.Integer(1)) -> tuple[sp.Expr, ...]:
    """Return h,D,R,L,B,phi_1,phi_2 for h=1+scale*seed."""
    h = sp.expand(1 + scale * seed)
    assert sp.expand(h.subs(z, 0)) == 1
    hp = sp.diff(h, z)
    D = sp.expand(h - z * hp)
    quotient, remainder = sp.div(sp.Poly(D - 1, z), sp.Poly(z, z))
    assert remainder.is_zero
    R = quotient.as_expr()
    L = sp.expand(D * y + z * R)
    B = sp.expand(h * R * (1 + y) - h * hp * (1 + y) ** 2)
    phi_1 = h
    phi_2 = sp.expand(-h * R * (1 + y) + L * B)
    return h, D, R, L, B, phi_1, phi_2


Sparse = dict[tuple[int, int], Fraction]


def sparse(expr: sp.Expr) -> Sparse:
    out: Sparse = {}
    for (z_degree, y_degree), coefficient in sp.Poly(expr, z, y).terms():
        rational = sp.Rational(coefficient)
        out[(z_degree, y_degree)] = Fraction(int(rational.p), int(rational.q))
    return out


def multiply(left: Sparse, right: Sparse, cutoff: int) -> Sparse:
    out: Sparse = {}
    for (az, ay), ac in left.items():
        for (bz, by), bc in right.items():
            exponent = az + bz, ay + by
            if exponent[0] > cutoff or exponent[1] > cutoff:
                continue
            out[exponent] = out.get(exponent, Fraction(0)) + ac * bc
    return {exponent: coefficient for exponent, coefficient in out.items() if coefficient}


def power(poly: Sparse, exponent: int, cutoff: int) -> Sparse:
    result: Sparse = {(0, 0): Fraction(1)}
    base = poly
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply(result, base, cutoff)
        base = multiply(base, base, cutoff)
        remaining //= 2
    return result


def gaussian_moment(phi_1: sp.Expr, phi_2: sp.Expr, observable: sp.Expr, m: int) -> Fraction:
    """Evaluate E[A(Z_1,Z_2) P^m] by exact circular Wick contraction."""
    p1 = sparse(phi_1)
    p2 = sparse(phi_2)
    obs = sparse(observable)
    normalized = Fraction(0)
    for a in range(m + 1):
        b = m - a
        product = multiply(power(p1, a, m), power(p2, b, m), m)
        product = multiply(product, obs, m)
        normalized += product.get((a, b), Fraction(0))
    return factorial(m) * normalized


def coefficient(expr: sp.Expr, variable: sp.Symbol, degree: int) -> Fraction:
    value = sp.Poly(expr, variable).coeff_monomial(variable**degree)
    rational = sp.Rational(value)
    return Fraction(int(rational.p), int(rational.q))


def main() -> None:
    # Symbolic fixed-branch and determinant cancellation for a seed with a
    # parameter.  This also checks the exact pencil slice
    # H(g)-s*g+t_0=0, s=(lambda*u)^(-1), t_0=lambda^(-1).
    lam = sp.symbols("lambda", nonzero=True)
    seed = z**3 * (1 - z)
    h, D, R, L, B, phi_1, phi_2 = bridge(seed, lam)
    branch = {u: z / h, y: 1 / D - 1}
    assert sp.factor((u * phi_1 - z).subs(branch)) == 0
    assert sp.factor((u * phi_2 - y).subs(branch)) == 0
    determinant = (1 - u * sp.diff(phi_1, z)) * (1 - u * sp.diff(phi_2, y))
    assert sp.factor(determinant.subs(branch)) == 1
    pencil = seed - z / (lam * u) + 1 / lam
    assert sp.factor(pencil.subs(u, z / h)) == 0

    # Check the two identities that make the correction transparent on the
    # graph L=0: 1+y=D^(-1) and the y-derivative has its prescribed value.
    assert sp.factor(L.subs(y, 1 / D - 1)) == 0
    derivative_on_branch = sp.factor(sp.diff(phi_2, y).subs(y, 1 / D - 1))
    assert sp.factor(derivative_on_branch + h * sp.diff(h, z) / D) == 0
    print("PASS weighted Gaussian bridge: universal fixed branch and determinant cancellation")

    examples = [
        ("canonical-cubic", z**2 * (1 - z), sp.Integer(1)),
        ("canonical-quartic", z**3 * (1 - z), sp.Integer(1)),
        ("canonical-quintic", z**4 * (1 - z), sp.Rational(2, 3)),
        ("split-quartic", z**2 * (z - 1) * (z - 3), sp.Integer(1)),
    ]
    for name, example_seed, example_scale in examples:
        h, D, R, L, B, phi_1, phi_2 = bridge(example_seed, example_scale)
        assert sp.Poly(phi_2, z, y).total_degree() >= 1
        for m in range(1, 13):
            pure = gaussian_moment(phi_1, phi_2, sp.Integer(1), m)
            mixed = gaussian_moment(phi_1, phi_2, z, m)
            expected = factorial(m - 1) * coefficient(h**m, z, m - 1)
            assert pure == 0
            assert mixed == expected
        print(f"PASS weighted Gaussian bridge: {name} Wick moments through m=12")

    # The most direct three-real-Gaussian half-pair ansatz would require this
    # rational correction v.  It is polynomial for Long's affine h=1+z, but
    # not for the audited nonlinear weighted seeds.  This is an obstruction
    # only to that ansatz, not to every possible dimension-three bridge.
    affine_h = 1 + z
    affine_D = sp.expand(affine_h - z * sp.diff(affine_h, z))
    affine_v = sp.cancel(
        -affine_h
        * sp.diff(affine_h, z)
        * (2 * affine_h - z * sp.diff(affine_h, z))
        / (2 * affine_D**2)
    )
    assert sp.factor(affine_v + (1 + z) * (2 + z) / 2) == 0
    for _, example_seed, example_scale in examples:
        h = sp.expand(1 + example_scale * example_seed)
        D = sp.expand(h - z * sp.diff(h, z))
        numerator = sp.expand(-h * sp.diff(h, z) * (2 * h - z * sp.diff(h, z)))
        _, remainder = sp.div(sp.Poly(numerator, z), sp.Poly(D**2, z))
        assert not remainder.is_zero
    print("PASS weighted Gaussian bridge: exact half-pair obstruction on audited nonlinear seeds")


if __name__ == "__main__":
    main()
