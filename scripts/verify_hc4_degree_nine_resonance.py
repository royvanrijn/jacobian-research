#!/usr/bin/env python3
from __future__ import annotations
import itertools
import json
from pathlib import Path
import sympy as sp
ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_DIR = ROOT / 'artifacts' / 'generated-results'
ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
x, y, z = sp.symbols('x y z')
lam, mu, ratio, a = sp.symbols('lam mu ratio a')

def binary_border(f: sp.Expr) -> sp.Matrix:
    gradient = sp.Matrix([sp.diff(f, x), sp.diff(f, y)])
    hessian = sp.hessian(f, (x, y))
    return sp.Matrix.vstack(sp.Matrix.hstack(sp.zeros(1, 1), gradient.T), sp.Matrix.hstack(gradient, hessian))

def schur_face(f: sp.Expr, g: sp.Expr, q: sp.Expr) -> sp.Expr:
    border = binary_border(f)
    column = sp.Matrix([g, sp.diff(g, x), sp.diff(g, y)])
    return sp.expand(q * border.det() - (column.T * border.adjugate() * column)[0])

def bordered(c: sp.Expr) -> sp.Expr:
    gradient = sp.Matrix([sp.diff(c, variable) for variable in (x, y, z)])
    hessian = sp.hessian(c, (x, y, z))
    return sp.expand((gradient.T * hessian.adjugate() * gradient)[0])

def solve_line_q(f: sp.Expr, g: sp.Expr, prefix: str) -> tuple[sp.Expr, list[sp.Expr]]:
    q0, q1 = sp.symbols(f'{prefix}_q0 {prefix}_q1')
    q = q0 * x + q1 * y
    equations = [coefficient for _, coefficient in sp.Poly(schur_face(f, g, q), x, y).terms()]
    matrix = sp.Matrix([[sp.diff(equation, variable) for variable in (q0, q1)] for equation in equations])
    pivot_rows = None
    for rows in itertools.combinations(range(len(equations)), 2):
        if sp.factor(matrix[list(rows), :].det()) != 0:
            pivot_rows = rows
            break
    assert pivot_rows is not None
    solution = sp.solve([equations[row] for row in pivot_rows], (q0, q1), dict=True, simplify=False)[0]
    residuals = []
    for index, equation in enumerate(equations):
        if index in pivot_rows:
            continue
        numerator = sp.factor(sp.together(equation.subs(solution)).as_numer_denom()[0])
        if numerator != 0:
            residuals.append(numerator)
    return (sp.expand(q.subs(solution)), residuals)

def resonance(d: int, m: int, e: int, n: int) -> int:
    return d * d * m + d * d * n * n - 2 * d * e * m * n - 2 * d * e * m - d * m * m - d * m + e * e * m * m + 2 * e * m * m + m * m
assert resonance(9, 3, 5, 1) == 0
direct_empty = ((x ** 4 * y ** 3 * (x - y) * (x - lam * y), x ** 2 * y * (x - y) * (x - lam * y), (lam,), 'p4311'), (x ** 3 * y ** 2 * (x - y) ** 2 * (x - lam * y) * (x - mu * y), x * y * (x - y) * (x - lam * y) * (x - mu * y), (lam, mu), 'p32211'), (x ** 3 * y ** 2 * (x - y) ** 2 * (x - lam * y) ** 2, x * y * (x - y) * (x - lam * y) * (x + ratio * y), (ratio, lam), 'p3222finite'), (x ** 3 * y ** 2 * (x - y) ** 2 * (x - lam * y) ** 2, x * y ** 2 * (x - y) * (x - lam * y), (lam,), 'p3222infinity'))
for f, g, variables, prefix in direct_empty:
    _, residuals = solve_line_q(f, g, prefix)
    basis = sp.groebner(residuals, *variables, order='lex')
    assert len(basis.polys) == 1 and basis.polys[0].as_expr() == 1
f_3321 = x ** 3 * y ** 3 * (x - y) ** 2 * (x - lam * y)
root_3321 = x * y * (x - y) * (x - lam * y)
_, residuals = solve_line_q(f_3321, root_3321 * (x + ratio * y), 'p3321finite')
basis = sp.groebner(residuals, ratio, lam, order='lex')
assert [sp.factor(polynomial.as_expr()) for polynomial in basis.polys] == [(3 * ratio + 2 * lam - 5) ** 2, (ratio - 1) * (lam - 1) ** 2, (lam - 1) ** 3]
_, residuals = solve_line_q(f_3321, root_3321 * y, 'p3321infinity')
assert sp.groebner(residuals, lam, order='lex').polys[0].as_expr() == 1
f_33111 = x ** 3 * y ** 3 * (x - y) * (x - lam * y) * (x - mu * y)
root_33111 = x * y * (x - y) * (x - lam * y) * (x - mu * y)
q_33111, residuals = solve_line_q(f_33111, root_33111, 'p33111')
basis = sp.groebner(residuals, lam, mu, order='lex')
assert [sp.factor(polynomial.as_expr()) for polynomial in basis.polys] == [lam + mu + 1, mu ** 2 + mu + 1]
lam_value = -mu - 1
potential = f_33111.subs(lam, lam_value) + a * z * root_33111.subs(lam, lam_value) + a ** 2 * z ** 2 * q_33111.subs(lam, lam_value) / 2
polynomial = sp.Poly(bordered(potential), x, y, z)
coefficient = polynomial.coeff_monomial(x ** 16 * y ** 4 * z)
remainder = sp.rem(sp.Poly(coefficient / a ** 3, mu, domain=sp.QQ), sp.Poly(mu ** 2 + mu + 1, mu, domain=sp.QQ)).as_expr()
assert remainder == 36
result = {'scope': 'degree-nine synchronized scalar m=3,e=5 resonance', 'status': 'closed', 'only_schur_survivor': ['lambda+mu+1', 'mu^2+mu+1'], 'terminal_coefficient': '[x^16 y^4 z] = 36 a^3'}
output = ARTIFACT_DIR / 'hc4_degree9_resonance.json'
output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print(json.dumps(result, indent=2))
