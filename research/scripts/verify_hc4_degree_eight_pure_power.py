#!/usr/bin/env python3
"""Exact degree-eight HC4 closure of the pure eighth-power packet.

Scoped to the synchronized scalar reverse-Schur packet."""
from __future__ import annotations

import itertools
import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_DIR = ROOT / "artifacts" / "generated-results"
ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)

x, y, z = sp.symbols("x y z")


def binary_border(f: sp.Expr) -> sp.Matrix:
    grad = sp.Matrix([sp.diff(f, x), sp.diff(f, y)])
    hess = sp.hessian(f, (x, y))
    return sp.Matrix.vstack(
        sp.Matrix.hstack(sp.zeros(1, 1), grad.T),
        sp.Matrix.hstack(grad, hess),
    )


def schur_face(f: sp.Expr, g: sp.Expr, q: sp.Expr) -> sp.Expr:
    border = binary_border(f)
    column = sp.Matrix([g, sp.diff(g, x), sp.diff(g, y)])
    return sp.expand(
        q * border.det() - (column.T * border.adjugate() * column)[0]
    )


def bordered(c: sp.Expr) -> sp.Expr:
    grad = sp.Matrix([sp.diff(c, v) for v in (x, y, z)])
    hess = sp.hessian(c, (x, y, z))
    return sp.expand((grad.T * hess.adjugate() * grad)[0])


def homogeneous(symbols: tuple[sp.Symbol, ...], degree: int) -> sp.Expr:
    if degree < 0:
        return sp.Integer(0)
    return sum(symbols[i] * x ** (degree - i) * y**i for i in range(degree + 1))


def radical_power(
    basis: sp.GroebnerBasis, polynomial: sp.Expr, maximum: int = 12
) -> int | None:
    for power in range(1, maximum + 1):
        if basis.reduce(sp.expand(polynomial**power))[1] == 0:
            return power
    return None


def weighted_face(
    f: sp.Expr,
    g: sp.Expr,
    q: sp.Expr,
    e: int,
    prefix: str,
) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    """Complete weight-eight potential with wt(z)=8-e."""
    r_degree = 3 * e - 16
    s_degree = 4 * e - 24
    r_symbols = (
        sp.symbols(f"{prefix}_r0:{r_degree + 1}") if r_degree >= 0 else ()
    )
    s_symbols = (
        sp.symbols(f"{prefix}_s0:{s_degree + 1}") if s_degree >= 0 else ()
    )
    r = homogeneous(tuple(r_symbols), r_degree)
    s = homogeneous(tuple(s_symbols), s_degree)
    potential = f + z * g + z**2 * q / 2 + z**3 * r / 6 + z**4 * s / 24
    return potential, tuple(r_symbols) + tuple(s_symbols)


def resonance(d: int, m: int, e: int, n: int) -> int:
    return (
        d * d * m
        + d * d * n * n
        - 2 * d * e * m * n
        - 2 * d * e * m
        - d * m * m
        - d * m
        + e * e * m * m
        + 2 * e * m * m
        + m * m
    )


def root_weight(d: int, multiplicity: int, transverse_degree: int) -> int:
    ceiling = (multiplicity + 1) // 2
    for order in range(ceiling):
        if resonance(d, multiplicity, transverse_degree, order) == 0:
            return order
    return ceiling


# ---------------------------------------------------------------------------
# 7. Pure eighth power
# ---------------------------------------------------------------------------
# e<=4 is already a fixed cylinder: every same-weight term uses only x,z.
for e in range(5):
    assert root_weight(8, 8, e) == e
    assert 2 * e - 8 <= 0

# e=5: the x^4*y*z direction and then every y-dependent q coefficient die.
a5, b5 = sp.symbols("a5 b5")
Q0, Q1, Q2 = sp.symbols("Q0 Q1 Q2")
pure5 = x**8 + z * x**4 * (a5 * x + b5 * y) + z**2 * (
    Q0 * x**2 + Q1 * x * y + Q2 * y**2
) / 2
Jpure5 = sp.Poly(bordered(pure5), x, y, z)
assert Jpure5.coeff_monomial(x**22) == -64 * b5**2
Jpure5_b = sp.Poly(Jpure5.as_expr().subs(b5, 0), x, y, z)
assert Jpure5_b.coeff_monomial(x**14 * y**2 * z**2) == -192 * Q2**2
assert sp.factor(
    Jpure5_b.coeff_monomial(x**16 * z**2).subs(Q2, 0) + 64 * Q1**2
) == 0

# e=6: the constant Schur face kills the two y-dependent coefficients of g;
# four further square coefficients kill q1..q4 and then r1,r2.
a6, b6, c6 = sp.symbols("a6 b6 c6")
qq = tuple(sp.symbols("pure6_q0:5"))
rr = tuple(sp.symbols("pure6_r0:3"))
ss = sp.symbols("pure6_s")
gpure6 = x**4 * (a6 * x**2 + b6 * x * y + c6 * y**2)
qpure6 = homogeneous(qq, 4)
constant_pure6 = sp.Poly(schur_face(x**8, gpure6, qpure6), x, y)
assert constant_pure6.coeff_monomial(x**24) == 64 * b6**2
assert constant_pure6.coeff_monomial(x**23 * y) == 256 * b6 * c6
assert constant_pure6.coeff_monomial(x**22 * y**2) == 256 * c6**2
pure6 = (
    x**8
    + a6 * x**6 * z
    + z**2 * qpure6 / 2
    + z**3 * homogeneous(rr, 2) / 6
    + z**4 * ss / 24
)
Jpure6 = sp.Poly(bordered(pure6), x, y, z)
assert Jpure6.coeff_monomial(x**14 * y**6 * z**2) == -640 * qq[4] ** 2
assert sp.factor(
    Jpure6.coeff_monomial(x**16 * y**4 * z**2).subs(qq[4], 0)
    + 384 * qq[3] ** 2
) == 0
assert sp.factor(
    Jpure6.coeff_monomial(x**18 * y**2 * z**2).subs({qq[4]: 0, qq[3]: 0})
    + 192 * qq[2] ** 2
) == 0
assert sp.factor(
    Jpure6.coeff_monomial(x**20 * z**2).subs(
        {qq[4]: 0, qq[3]: 0, qq[2]: 0}
    )
    + 64 * qq[1] ** 2
) == 0
q_tail_zero = {qq[i]: 0 for i in range(1, 5)}
assert sp.factor(
    Jpure6.coeff_monomial(x**14 * y**2 * z**4).subs(q_tail_zero)
    + sp.Rational(128, 3) * rr[2] ** 2
) == 0
assert sp.factor(
    Jpure6.coeff_monomial(x**16 * z**4).subs({**q_tail_zero, rr[2]: 0})
    + 16 * rr[1] ** 2
) == 0



result = {
    "scope": "pure eighth-power synchronized scalar packet",
    "status": "fixed cylinder",
}
output = ARTIFACT_DIR / "hc4_degree8_pure_power.json"
output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
