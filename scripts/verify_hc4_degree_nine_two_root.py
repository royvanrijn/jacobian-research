#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
import sympy as sp
ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_DIR = ROOT / 'artifacts' / 'generated-results'
ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
x, y, z = sp.symbols('x y z')

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
    return sum((symbols[i] * x ** (degree - i) * y ** i for i in range(degree + 1)))

def resonance(d: int, m: int, e: int, n: int) -> int:
    return d * d * m + d * d * n * n - 2 * d * e * m * n - 2 * d * e * m - d * m * m - d * m + e * e * m * m + 2 * e * m * m + m * m

def root_weight(d: int, multiplicity: int, transverse_degree: int) -> int:
    for order in range((multiplicity + 1) // 2):
        if resonance(d, multiplicity, transverse_degree, order) == 0:
            return order
    return (multiplicity + 1) // 2

def radical_power(basis: sp.GroebnerBasis, p: sp.Expr, limit: int=15) -> int | None:
    for power in range(1, limit + 1):
        if basis.reduce(sp.expand(p ** power))[1] == 0:
            return power
    return None

def solve_q(f: sp.Expr, g: sp.Expr, degree: int) -> sp.Expr:
    q_symbols = tuple(sp.symbols(f'q0:{degree + 1}'))
    q = homogeneous(q_symbols, degree)
    equations = sp.Poly(schur_face(f, g, q), x, y).coeffs()
    solutions = sp.solve(equations, q_symbols, dict=True, simplify=False)
    assert len(solutions) == 1
    return sp.expand(q.subs(solutions[0]))

def complete_face(f: sp.Expr, g: sp.Expr, q: sp.Expr, e: int, prefix: str):
    potential = f + z * g + z ** 2 * q / 2
    tails: list[sp.Symbol] = []
    r_degree = 3 * e - 18
    if r_degree >= 0:
        r = tuple(sp.symbols(f'{prefix}_r0:{r_degree + 1}'))
        potential += z ** 3 * homogeneous(r, r_degree) / 6
        tails.extend(r)
    s_degree = 4 * e - 27
    if s_degree >= 0:
        s = tuple(sp.symbols(f'{prefix}_s0:{s_degree + 1}'))
        potential += z ** 4 * homogeneous(s, s_degree) / 24
        tails.extend(s)
    return (potential, tuple(tails))
roots = (x, y, x - y)
closed: dict[str, dict[str, list[int]]] = {}
for label, multiplicities in {'8+1': (8, 1), '7+2': (7, 2), '6+3': (6, 3), '5+4': (5, 4)}.items():
    f = x ** multiplicities[0] * y ** multiplicities[1]
    closed[label] = {}
    for e in range(5, 8):
        orders = [root_weight(9, m, e) for m in multiplicities]
        weight = sum(orders)
        if weight > e:
            continue
        root_factor = x ** orders[0] * y ** orders[1]
        h_degree = e - weight
        h = tuple(sp.symbols(f"p{label.replace('+', '_')}_e{e}_h0:{h_degree + 1}"))
        g = sp.expand(root_factor * homogeneous(h, h_degree))
        q = solve_q(f, g, 2 * e - 9)
        potential, tails = complete_face(f, g, q, e, f"p{label.replace('+', '_')}_e{e}")
        basis = sp.groebner(sp.Poly(bordered(potential), x, y, z).coeffs(), *tails, *h, order='grevlex')
        powers = [radical_power(basis, coefficient) for coefficient in h]
        assert all((power is not None for power in powers))
        closed[label][str(e)] = [int(power) for power in powers]
result = {'scope': 'degree-nine synchronized scalar packet with exactly two leading roots', 'status': 'closed to fixed cylinders', 'partitions': closed}
(ARTIFACT_DIR / 'hc4_degree9_two_root.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print(json.dumps(result, indent=2))
