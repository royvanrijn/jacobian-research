#!/usr/bin/env python3
"""Close degree-nine scalar HC4 packets with five odd root multiplicities."""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_DIR = ROOT / "artifacts" / "generated-results"
ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)

x, y, z = sp.symbols("x y z")
lam, mu, a = sp.symbols("lam mu a")
r = sp.symbols("r0:4")
s = sp.symbols("s0:2")


def homogeneous(symbols: tuple[sp.Symbol, ...], degree: int) -> sp.Expr:
    return sum(
        symbols[index] * x ** (degree - index) * y**index
        for index in range(degree + 1)
    )


def bordered(c: sp.Expr) -> sp.Expr:
    gradient = sp.Matrix([sp.diff(c, variable) for variable in (x, y, z)])
    hessian = sp.hessian(c, (x, y, z))
    return sp.expand((gradient.T * hessian.adjugate() * gradient)[0])


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


def local_value(m: int) -> sp.Rational:
    n = (m + 1) // 2
    return -sp.Rational(resonance(9, m, 7, n), 9 * m * (9 - m))


values = {m: local_value(m) for m in range(1, 9)}
assert values[1] == values[3] == sp.Rational(11, 18)

# With five odd multiplicities, H and C are constants at e=7.  Every
# partition except (3,3,1,1,1) contains unequal local values and dies before
# the complete face.
partitions = (
    (5, 1, 1, 1, 1),
    (4, 1, 1, 1, 1, 1),
    (3, 3, 1, 1, 1),
    (3, 2, 1, 1, 1, 1),
    (2, 2, 1, 1, 1, 1, 1),
)
for partition in partitions:
    value_set = {values[multiplicity] for multiplicity in set(partition)}
    if partition == (3, 3, 1, 1, 1):
        assert value_set == {sp.Rational(11, 18)}
    else:
        assert len(value_set) > 1

B = x * y * (x - y) * (x - lam * y) * (x - mu * y)
potential = (
    x**2 * y**2 * B
    + a * z * x * y * B
    + sp.Rational(11, 36) * a**2 * z**2 * B
    + z**3 * homogeneous(r, 3) / 6
    + z**4 * homogeneous(s, 1) / 24
)
field = sp.QQ.frac_field(lam, mu)
coefficients = sp.Poly(
    bordered(potential), x, y, z, domain=field[r + s + (a,)]
).coeffs()
basis = sp.groebner(
    coefficients, *r, *s, a, order="grevlex", domain=field
)
expressions = [sp.factor(polynomial.as_expr()) for polynomial in basis.polys]
expected = [a**2, *r, *s]
assert expressions == expected

result = {
    "scope": "degree-nine synchronized scalar packets with five odd multiplicities",
    "status": "closed",
    "unique_local_survivor": "(3,3,1,1,1)",
    "complete_face_basis": [str(value) for value in expressions],
}
output = ARTIFACT_DIR / "hc4_degree9_five_odd.json"
output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
