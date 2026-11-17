#!/usr/bin/env python3
from __future__ import annotations
import gc
import itertools
import json
from pathlib import Path
import sympy as sp
ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_DIR = ROOT / 'artifacts' / 'generated-results'
ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
x, y, z = sp.symbols('x y z')

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

def homogeneous(symbols: tuple[sp.Symbol, ...], degree: int) -> sp.Expr:
    return sum((symbols[index] * x ** (degree - index) * y ** index for index in range(degree + 1)))

def resonance(d: int, m: int, e: int, n: int) -> int:
    return d * d * m + d * d * n * n - 2 * d * e * m * n - 2 * d * e * m - d * m * m - d * m + e * e * m * m + 2 * e * m * m + m * m

def local_value(m: int, e: int) -> sp.Rational:
    n = (m + 1) // 2
    return -sp.Rational(resonance(9, m, e, n), 9 * m * (9 - m))

def complete_e7_basis(f: sp.Expr, g: sp.Expr, q: sp.Expr, lam: sp.Symbol, a: sp.Symbol) -> list[sp.Expr]:
    r = sp.symbols('r0:4')
    s = sp.symbols('s0:2')
    potential = f + z * g + z ** 2 * q / 2 + z ** 3 * homogeneous(r, 3) / 6 + z ** 4 * homogeneous(s, 1) / 24
    field = sp.QQ.frac_field(lam)
    coefficients = sp.Poly(bordered(potential), x, y, z, domain=field[r + s + (a,)]).coeffs()
    basis = sp.groebner(coefficients, *r, *s, a, order='grevlex', domain=field)
    result = [sp.factor(polynomial.as_expr()) for polynomial in basis.polys]
    del potential, coefficients, basis
    sp.core.cache.clear_cache()
    gc.collect()
    return result

def solve_linear_q(f: sp.Expr, g: sp.Expr, q_symbols: tuple[sp.Symbol, ...], q: sp.Expr) -> list[sp.Expr]:
    equations = [coefficient for _, coefficient in sp.Poly(schur_face(f, g, q), x, y).terms()]
    matrix = sp.Matrix([[sp.diff(equation, variable) for variable in q_symbols] for equation in equations])
    pivot_rows = None
    for rows in itertools.combinations(range(len(equations)), len(q_symbols)):
        if sp.factor(matrix[list(rows), :].det()) != 0:
            pivot_rows = rows
            break
    assert pivot_rows is not None
    solution = sp.solve([equations[row] for row in pivot_rows], q_symbols, dict=True, simplify=False)[0]
    residuals = []
    for index, equation in enumerate(equations):
        if index in pivot_rows:
            continue
        numerator = sp.factor(sp.together(equation.subs(solution)).as_numer_denom()[0])
        if numerator != 0:
            residuals.append(numerator)
    return residuals
values_e6 = {m: local_value(m, 6) for m in range(1, 9)}
values_e7 = {m: local_value(m, 7) for m in range(1, 9)}
assert values_e7[1] == values_e7[3] == sp.Rational(11, 18)
partitions = ((6, 1, 1, 1), (5, 2, 1, 1), (4, 3, 1, 1), (4, 2, 1, 1, 1), (3, 3, 2, 1), (3, 2, 2, 1, 1), (2, 2, 2, 1, 1, 1))
for partition in partitions:
    assert len({values_e6[multiplicity] for multiplicity in set(partition)}) > 1
lam, a = sp.symbols('lam a')
expected = [a ** 2, *sp.symbols('r0:4'), *sp.symbols('s0:2')]
for f, g, q in ((x ** 6 * y * (x - y) * (x - lam * y), a * x ** 4 * y * (x - y) * (x - lam * y), sp.Rational(11, 18) * a ** 2 * x ** 2 * y * (x - y) * (x - lam * y)), (x ** 4 * y ** 3 * (x - y) * (x - lam * y), a * x ** 3 * y ** 2 * (x - y) * (x - lam * y), sp.Rational(11, 18) * a ** 2 * x ** 2 * y * (x - y) * (x - lam * y)), (x ** 3 * y ** 3 * (x - y) ** 2 * (x - lam * y), a * x ** 2 * y ** 2 * (x - y) ** 2 * (x - lam * y), sp.Rational(11, 18) * a ** 2 * x * y * (x - y) ** 2 * (x - lam * y))):
    assert complete_e7_basis(f, g, q, lam, a) == expected
u, v, ratio = sp.symbols('u v ratio')
q_symbols = sp.symbols('q0:6')
f = x ** 5 * y ** 2 * (x - y) * (x - lam * y)
g = x ** 3 * y * (x - y) * (x - lam * y) * (u * x + v * y)
residuals = solve_linear_q(f, g, q_symbols, homogeneous(q_symbols, 5))
finite = [residual.subs({u: 1, v: ratio}) for residual in residuals]
assert sp.groebner(finite, ratio, lam, order='lex').polys[0].as_expr() == 1
infinity = [residual.subs({u: 0, v: 1}) for residual in residuals]
assert sp.Poly(sp.gcd_list(infinity), lam).degree() == 0
assert values_e7[1] == values_e7[3]
assert values_e7[4] != values_e7[1]
assert values_e7[2] != values_e7[1]
result = {'scope': 'degree-nine synchronized scalar packets with three odd multiplicities', 'status': 'closed', 'e6_root_values': {str(key): str(value) for key, value in values_e6.items()}, 'e7_root_values': {str(key): str(value) for key, value in values_e7.items()}, 'complete_face_basis': [str(value) for value in expected]}
output = ARTIFACT_DIR / 'hc4_degree9_three_odd.json'
output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print(json.dumps(result, indent=2))
