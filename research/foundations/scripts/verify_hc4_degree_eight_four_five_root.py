#!/usr/bin/env python3
"""Exact degree-eight HC4 closure of the four-/five-root and universal weight-six packets.

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
# 3. Four-root W=5 charts: (4,2,1,1) and (3,2,2,1)
# ---------------------------------------------------------------------------
lam = sp.symbols("lam")


def solve_linear_q(
    f: sp.Expr, g: sp.Expr, q_symbols: tuple[sp.Symbol, ...], q: sp.Expr
) -> tuple[dict[sp.Symbol, sp.Expr], list[sp.Expr]]:
    polynomial = sp.Poly(schur_face(f, g, q), x, y)
    equations = [coefficient for _, coefficient in sorted(
        polynomial.terms(), key=lambda term: -term[0][0]
    )]
    matrix = sp.Matrix(
        [[sp.diff(equation, variable) for variable in q_symbols] for equation in equations]
    )
    pivot_rows = None
    for rows in itertools.combinations(range(len(equations)), len(q_symbols)):
        determinant = sp.factor(matrix[list(rows), :].det())
        if determinant != 0 and not determinant.has(lam):
            pivot_rows = rows
            break
    assert pivot_rows is not None
    solution = sp.solve(
        [equations[row] for row in pivot_rows],
        q_symbols,
        dict=True,
        simplify=False,
    )[0]
    residuals = [
        sp.factor(equation.subs(solution))
        for index, equation in enumerate(equations)
        if index not in pivot_rows
    ]
    return solution, [residual for residual in residuals if residual != 0]


for label, multiplicities in {
    "4211": (4, 2, 1, 1),
    "3221": (3, 2, 2, 1),
}.items():
    roots = (x, y, x - y, x - lam * y)
    f = sp.prod(root**m for root, m in zip(roots, multiplicities))
    root_factor = sp.prod(
        root ** ((m + 1) // 2) for root, m in zip(roots, multiplicities)
    )
    for e in (5, 6):
        h_degree = e - 5
        hs = tuple(sp.symbols(f"p{label}_e{e}_h0:{h_degree + 1}"))
        g = root_factor * homogeneous(hs, h_degree)
        q_degree = 2 * e - 8
        qs = tuple(sp.symbols(f"p{label}_e{e}_q0:{q_degree + 1}"))
        q = homogeneous(qs, q_degree)
        _, residuals = solve_linear_q(f, g, qs, q)
        if e == 5:
            normalized = [
                sp.together(residual / hs[0] ** 2).as_numer_denom()[0]
                for residual in residuals
            ]
            assert sp.Poly(sp.gcd_list(normalized), lam).degree() == 0
        else:
            ratio = sp.symbols(f"p{label}_ratio")
            finite = [
                sp.together(residual.subs({hs[0]: 1, hs[1]: ratio}))
                .as_numer_denom()[0]
                for residual in residuals
            ]
            finite_basis = sp.groebner(finite, ratio, lam, order="lex")
            assert len(finite_basis.polys) == 1
            assert finite_basis.polys[0].as_expr() == 1
            infinity = [
                sp.together(residual.subs({hs[0]: 0, hs[1]: 1}))
                .as_numer_denom()[0]
                for residual in residuals
            ]
            assert sp.Poly(sp.gcd_list(infinity), lam).degree() == 0


# ---------------------------------------------------------------------------
# 4. Universal W=6 obstruction: f=A^2 B, g=A B, deg A=2, deg B=4
# ---------------------------------------------------------------------------
A, A1, A2, B, B1, B2, C = sp.symbols("A A1 A2 B B1 B2 C")


def affine_schur_data(
    d: int,
    e: int,
    F: sp.Expr,
    F1: sp.Expr,
    F2: sp.Expr,
    G: sp.Expr,
    G1: sp.Expr,
) -> tuple[sp.Expr, sp.Expr]:
    denominator = -d * F * (d * F * F2 - (d - 1) * F1**2)
    numerator = (
        -d**2 * F**2 * G1**2
        + 2 * d * e * F * F1 * G * G1
        + d**2 * F * F2 * G**2
        - 2 * d * e * F * F2 * G**2
        - d * F * F2 * G**2
        - d**2 * F1**2 * G**2
        + 2 * d * e * F1**2 * G**2
        + 2 * d * F1**2 * G**2
        - e**2 * F1**2 * G**2
        - 2 * e * F1**2 * G**2
        - F1**2 * G**2
    )
    return sp.factor(denominator), sp.factor(numerator)


F = A**2 * B
F1 = 2 * A * A1 * B + A**2 * B1
F2 = 2 * (A1**2 + A * A2) * B + 4 * A * A1 * B1 + A**2 * B2
G = A * B
G1 = A1 * B + A * B1
denominator, numerator = affine_schur_data(8, 6, F, F1, F2, G, G1)
core_w6 = sp.factor((numerator - C * B * denominator) / (A**4 * B**2))
assert sp.factor(core_w6.subs(B, 0) - A**2 * B1**2 * (31 - 56 * C)) == 0
core_w6_special = sp.factor(core_w6.subs(C, sp.Rational(31, 56)))
expected_w6 = -sp.Rational(32, 7) * B * (
    A**2 * B2 - 3 * A * A1 * B1 + 2 * A * A2 * B + 2 * A1**2 * B
)
assert sp.factor(core_w6_special - expected_w6) == 0

# Quadratic normal forms.  A=t gives no nonzero quartic solution.  A=t^2
# leaves only t^3,t^4, which is never squarefree.
t = sp.symbols("t")
betas = sp.symbols("beta0:5")
quartic = sum(betas[i] * t**i for i in range(5))
operator_distinct = sp.expand(t**2 * sp.diff(quartic, t, 2) - 3 * t * sp.diff(quartic, t) + 2 * quartic)
assert [sp.expand(operator_distinct).coeff(t, i) for i in range(5)] == [
    2 * betas[0],
    -betas[1],
    -2 * betas[2],
    -betas[3],
    2 * betas[4],
]
operator_double = sp.expand(t**2 * sp.diff(quartic, t, 2) - 6 * t * sp.diff(quartic, t) + 12 * quartic)
assert [sp.expand(operator_double).coeff(t, i) for i in range(5)] == [
    12 * betas[0],
    6 * betas[1],
    2 * betas[2],
    0,
    0,
]
remaining_quartic = betas[3] * t**3 + betas[4] * t**4
assert sp.resultant(remaining_quartic, sp.diff(remaining_quartic, t), t) == 0


# ---------------------------------------------------------------------------
# 5. The only five-root W=5 partition: (2,2,2,1,1)
# ---------------------------------------------------------------------------
# e=5: q=C B.  A B are coprime, A has three simple roots and B two.
denominator5, numerator5 = affine_schur_data(8, 5, F, F1, F2, G, G1)
core_w5_e5 = sp.factor((numerator5 - C * B * denominator5) / (4 * A**4 * B**2))
assert sp.factor(core_w5_e5.subs(B, 0) - A**2 * B1**2 * (3 - 14 * C)) == 0
core_w5_e5_special = sp.factor(core_w5_e5.subs(C, sp.Rational(3, 14)))
expected_w5_e5 = -sp.Rational(2, 7) * B * (
    9 * A**2 * B2
    - 6 * A * A1 * B1
    + 18 * A * A2 * B
    - 10 * A1**2 * B
)
assert sp.factor(core_w5_e5_special - expected_w5_e5) == 0
assert sp.factor(
    (expected_w5_e5 / (-sp.Rational(2, 7) * B)).subs(A, 0)
    + 10 * A1**2 * B
) == 0

# e=6: g=A B L and q=B C.  Values of C at the three A-roots and the two
# B-roots are incompatible for a nonzero linear L.
L, L1 = sp.symbols("L L1")
G_linear = A * B * L
G_linear_1 = A1 * B * L + A * B1 * L + A * B * L1
denominator6, numerator6 = affine_schur_data(
    8, 6, F, F1, F2, G_linear, G_linear_1
)
core_w5_e6 = sp.factor(
    (numerator6 - B * C * denominator6) / (A**4 * B**2)
)
assert sp.factor(
    core_w5_e6.subs(A, 0) - A1**2 * B**2 * (44 * L**2 - 96 * C)
) == 0
assert sp.factor(
    core_w5_e6.subs(B, 0) - A**2 * B1**2 * (31 * L**2 - 56 * C)
) == 0
assert sp.Rational(11, 24) - sp.Rational(31, 56) == -sp.Rational(2, 21)



result = {
    "scope": "degree-eight four-/five-root and universal weight-six packets",
    "status": "closed",
}
output = ARTIFACT_DIR / "hc4_degree8_four_five_root.json"
output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
