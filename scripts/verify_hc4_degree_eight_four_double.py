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
lam = sp.symbols('lam')

def binary_border(f: sp.Expr) -> sp.Matrix:
    grad = sp.Matrix([sp.diff(f, x), sp.diff(f, y)])
    hess = sp.hessian(f, (x, y))
    return sp.Matrix.vstack(sp.Matrix.hstack(sp.zeros(1, 1), grad.T), sp.Matrix.hstack(grad, hess))

def schur_face(f: sp.Expr, g: sp.Expr, q: sp.Expr) -> sp.Expr:
    border = binary_border(f)
    column = sp.Matrix([g, sp.diff(g, x), sp.diff(g, y)])
    return sp.expand(q * border.det() - (column.T * border.adjugate() * column)[0])

def bordered(c: sp.Expr) -> sp.Expr:
    grad = sp.Matrix([sp.diff(c, v) for v in (x, y, z)])
    hess = sp.hessian(c, (x, y, z))
    return sp.expand((grad.T * hess.adjugate() * grad)[0])

def homogeneous(symbols: tuple[sp.Symbol, ...], degree: int) -> sp.Expr:
    if degree < 0:
        return sp.Integer(0)
    return sum((symbols[i] * x ** (degree - i) * y ** i for i in range(degree + 1)))

def radical_power(basis: sp.GroebnerBasis, polynomial: sp.Expr, maximum: int=12) -> int | None:
    for power in range(1, maximum + 1):
        if basis.reduce(sp.expand(polynomial ** power))[1] == 0:
            return power
    return None

def weighted_face(f: sp.Expr, g: sp.Expr, q: sp.Expr, e: int, prefix: str) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    r_degree = 3 * e - 16
    s_degree = 4 * e - 24
    r_symbols = sp.symbols(f'{prefix}_r0:{r_degree + 1}') if r_degree >= 0 else ()
    s_symbols = sp.symbols(f'{prefix}_s0:{s_degree + 1}') if s_degree >= 0 else ()
    r = homogeneous(tuple(r_symbols), r_degree)
    s = homogeneous(tuple(s_symbols), s_degree)
    potential = f + z * g + z ** 2 * q / 2 + z ** 3 * r / 6 + z ** 4 * s / 24
    return (potential, tuple(r_symbols) + tuple(s_symbols))

def solve_linear_q(f: sp.Expr, g: sp.Expr, q_symbols: tuple[sp.Symbol, ...], q: sp.Expr) -> tuple[dict[sp.Symbol, sp.Expr], list[sp.Expr]]:
    polynomial = sp.Poly(schur_face(f, g, q), x, y)
    equations = [coefficient for _, coefficient in sorted(polynomial.terms(), key=lambda term: -term[0][0])]
    matrix = sp.Matrix([[sp.diff(equation, variable) for variable in q_symbols] for equation in equations])
    pivot_rows = None
    for rows in itertools.combinations(range(len(equations)), len(q_symbols)):
        determinant = sp.factor(matrix[list(rows), :].det())
        if determinant != 0 and (not determinant.has(lam)):
            pivot_rows = rows
            break
    assert pivot_rows is not None
    solution = sp.solve([equations[row] for row in pivot_rows], q_symbols, dict=True, simplify=False)[0]
    residuals = [sp.factor(equation.subs(solution)) for index, equation in enumerate(equations) if index not in pivot_rows]
    return (solution, [residual for residual in residuals if residual != 0])
R4 = x * y * (x - y) * (x - lam * y)
f4 = R4 ** 2
a4 = sp.symbols('a4')
potential4 = f4 + a4 * z * R4 + a4 ** 2 * z ** 2 / 16
J4 = sp.Poly(bordered(potential4), x, y, z)
assert J4.coeff_monomial(x ** 13 * y ** 3 * z) == 9 * a4 ** 3
h50, h51 = sp.symbols('h50 h51')
g5 = R4 * (h50 * x + h51 * y)
q5s = sp.symbols('q5_0:3')
q5 = homogeneous(tuple(q5s), 2)
_, residual5 = solve_linear_q(f4, g5, tuple(q5s), q5)
ratio5 = sp.symbols('ratio5')
finite5 = [sp.together(residual.subs({h50: 1, h51: ratio5})).as_numer_denom()[0] for residual in residual5]
assert sp.groebner(finite5, ratio5, lam, order='lex').polys[0].as_expr() == 1
infinity5 = [sp.together(residual.subs({h50: 0, h51: 1})).as_numer_denom()[0] for residual in residual5]
assert sp.Poly(sp.gcd_list(infinity5), lam).degree() == 0
h60, h61, h62 = sp.symbols('h60 h61 h62')
g6 = R4 * (h60 * x ** 2 + h61 * x * y + h62 * y ** 2)
q6s = sp.symbols('q6_0:5')
q6 = homogeneous(tuple(q6s), 4)
solution6, residual6 = solve_linear_q(f4, g6, tuple(q6s), q6)
ratio61, ratio62 = sp.symbols('ratio61 ratio62')
finite6 = [sp.together(residual.subs({h60: 1, h61: ratio61, h62: ratio62})).as_numer_denom()[0] for residual in residual6]
basis6 = sp.groebner(finite6, ratio61, ratio62, lam, order='lex')
exception = (lam - 2) * (lam + 1) * (2 * lam - 1)
assert basis6.reduce(lam ** 2 * (lam - 1) ** 2 * exception)[1] == 0
ratio_inf = sp.symbols('ratio_inf')
infinity6 = [sp.together(residual.subs({h60: 0, h61: 1, h62: ratio_inf})).as_numer_denom()[0] for residual in residual6]
basis_inf = sp.groebner(infinity6, ratio_inf, lam, order='lex')
assert basis_inf.reduce(exception)[1] == 0
deep_infinity = [sp.together(residual.subs({h60: 0, h61: 0, h62: 1})).as_numer_denom()[0] for residual in residual6]
assert sp.Poly(sp.gcd_list(deep_infinity), lam).degree() == 0
q6_solution = q6.subs(solution6)
ap, bp = sp.symbols('ap bp')
harmonic_orbit = {sp.Integer(-1), 1 - sp.Integer(-1), 1 / sp.Integer(-1), 1 / (1 - sp.Integer(-1)), sp.Integer(-1) / (sp.Integer(-1) - 1), (sp.Integer(-1) - 1) / sp.Integer(-1)}
assert harmonic_orbit == {sp.Integer(-1), sp.Integer(2), sp.Rational(1, 2)}
exception_charts = ((sp.Integer(-1), {h60: ap, h61: bp, h62: -ap}),)
for index, (value, h_substitution) in enumerate(exception_charts):
    rr = tuple(sp.symbols(f'harm{index}_r0:3'))
    ss = sp.symbols(f'harm{index}_s')
    specialized = f4.subs(lam, value) + z * g6.subs(lam, value).subs(h_substitution) + z ** 2 * q6_solution.subs(lam, value).subs(h_substitution) / 2 + z ** 3 * homogeneous(rr, 2) / 6 + z ** 4 * ss / 24
    harmonic_basis = sp.groebner(sp.Poly(bordered(specialized), x, y, z).coeffs(), *rr, ss, ap, bp, order='grevlex')
    assert [sp.factor(p.as_expr()) for p in harmonic_basis.polys] == [ap ** 3, ap ** 2 * bp, ap * bp ** 2, bp ** 3, rr[0], rr[1], rr[2], ss]
result = {'scope': 'degree-eight four-double-root packet', 'status': 'closed', 'harmonic_cross_ratios': ['-1', '2', '1/2']}
output = ARTIFACT_DIR / 'hc4_degree8_four_double.json'
output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print(json.dumps(result, indent=2))
