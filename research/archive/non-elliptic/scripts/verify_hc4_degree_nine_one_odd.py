#!/usr/bin/env python3
"""Close nonresonant degree-nine scalar HC4 packets with one odd root."""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_DIR = ROOT / "artifacts" / "generated-results"
ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)

x, y, t = sp.symbols("x y t")
lam = sp.symbols("lam")


def binary_border(f: sp.Expr) -> sp.Matrix:
    gradient = sp.Matrix([sp.diff(f, x), sp.diff(f, y)])
    hessian = sp.hessian(f, (x, y))
    return sp.Matrix.vstack(
        sp.Matrix.hstack(sp.zeros(1, 1), gradient.T),
        sp.Matrix.hstack(gradient, hessian),
    )


def schur_face(f: sp.Expr, g: sp.Expr, q: sp.Expr) -> sp.Expr:
    border = binary_border(f)
    column = sp.Matrix([g, sp.diff(g, x), sp.diff(g, y)])
    return sp.expand(
        q * border.det() - (column.T * border.adjugate() * column)[0]
    )


def homogeneous(symbols: tuple[sp.Symbol, ...], degree: int) -> sp.Expr:
    return sum(
        symbols[index] * x ** (degree - index) * y**index
        for index in range(degree + 1)
    )


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


def local_value(m: int, e: int) -> sp.Rational:
    n = (m + 1) // 2
    return -sp.Rational(resonance(9, m, e, n), 9 * m * (9 - m))


# At e=5 the nonresonant (4,2,2,1) and (2,2,2,2,1) packets have constant
# H and C, but incompatible local values.
assert len({local_value(m, 5) for m in (4, 2, 1)}) == 3
assert local_value(2, 5) != local_value(1, 5)

# ---------------------------------------------------------------------------
# (4,2,2,1): complete projective Schur cover at e=6 and e=7
# ---------------------------------------------------------------------------

f_4221 = x**4 * y**2 * (x - y) ** 2 * (x - lam * y)
root_4221 = x**2 * y * (x - y) * (x - lam * y)
ratio = sp.symbols("ratio")

q6 = sp.symbols("q6_0:4")
for H, variables in (
    (x + ratio * y, (*q6, ratio, lam)),
    (y, (*q6, lam)),
):
    equations = sp.Poly(
        schur_face(f_4221, root_4221 * H, homogeneous(q6, 3)), x, y
    ).coeffs()
    basis = sp.groebner(equations, *variables, order="lex")
    assert len(basis.polys) == 1 and basis.polys[0].as_expr() == 1

q7 = sp.symbols("q7_0:6")
h1, h2 = sp.symbols("h1 h2")
charts = (
    (x**2 + h1 * x * y + h2 * y**2, (*q7, h1, h2, lam), "leading"),
    (x * y + h2 * y**2, (*q7, h2, lam), "middle"),
    (y**2, (*q7, lam), "trailing"),
)
for H, variables, chart in charts:
    equations = sp.Poly(
        schur_face(f_4221, root_4221 * H, homogeneous(q7, 5)), x, y
    ).coeffs()
    basis = sp.groebner(equations, *variables, order="lex")
    if chart == "leading":
        # Only lambda=0 or 1 survive; both are lower-root collisions.
        assert sp.factor(basis.polys[-1].as_expr()) == lam**2 * (lam - 1) ** 2
    else:
        assert len(basis.polys) == 1 and basis.polys[0].as_expr() == 1

# ---------------------------------------------------------------------------
# Odd root at infinity: ODE closure for (3,2,2,2) and (2,2,2,2,1)
# ---------------------------------------------------------------------------

# At e=6 interpolation gives C=47/126 H^2 and
# 21 A H'^2 + 3 H^2 A'' - 14 H A'H'=0.
# An affine normalization gives H=1 or H=t.  The first forces A''=0; the
# second is diagonal on t^j with coefficient 3j^2-17j+21.
for degree in (3, 4):
    assert all(3 * j * j - 17 * j + 21 != 0 for j in range(degree + 1))


def classify_e7_ode(degree: int, linear_multiplier: bool) -> dict[str, list[str]]:
    coefficients = sp.symbols(f"a{degree}_0:{degree}")
    A = t**degree + sum(
        coefficients[index] * t**index for index in range(degree)
    )
    if linear_multiplier:
        l0, l1 = sp.symbols(f"l{degree}_0 l{degree}_1")
        L = l1 * t + l0
        variables = (*coefficients, l0, l1)
    else:
        k = sp.symbols(f"k{degree}")
        L = k
        variables = (*coefficients, k)
    normal_forms = {
        "constant": sp.Integer(1),
        "linear": t,
        "double": t**2,
        "distinct": t * (t - 1),
    }
    record: dict[str, list[str]] = {}
    for name, H in normal_forms.items():
        equation = sp.expand(
            126 * A * L * sp.diff(A, t, 2)
            - 98 * L * sp.diff(A, t) ** 2
            - 63 * A * sp.diff(H, t) ** 2
            - 25 * H**2 * sp.diff(A, t, 2)
            + 70 * H * sp.diff(A, t) * sp.diff(H, t)
        )
        basis = sp.groebner(
            sp.Poly(equation, t).all_coeffs(), *variables, order="lex"
        )
        expressions = [sp.factor(polynomial.as_expr()) for polynomial in basis.polys]
        record[name] = [str(expression) for expression in expressions]
        if name != "double":
            assert expressions == [1]
        elif degree == 3:
            l0, l1 = variables[-2:]
            assert expressions == [
                1070696 * coefficients[0] + 1972593 * l0**3,
                11638 * coefficients[1] - 41013 * l0**2,
                23 * coefficients[2] + 63 * l0,
                l0**4,
                7 * l1 - 1,
            ]
        else:
            assert expressions == [*coefficients, 7 * variables[-1] - 1]
    return record


cubic_record = classify_e7_ode(3, True)
quartic_record = classify_e7_ode(4, False)
# In the double-root H chart these bases leave A=t^3 or A=t^4, hence no
# squarefree finite A and no packet with the declared distinct roots.

result = {
    "scope": "nonresonant degree-nine synchronized scalar packets with one odd multiplicity",
    "status": "closed",
    "cubic_ode": cubic_record,
    "quartic_ode": quartic_record,
}
output = ARTIFACT_DIR / "hc4_degree9_one_odd.json"
output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
