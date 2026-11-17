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

def resonance(d: int, m: int, e: int, n: int) -> int:
    return d * d * m + d * d * n * n - 2 * d * e * m * n - 2 * d * e * m - d * m * m - d * m + e * e * m * m + 2 * e * m * m + m * m

def root_weight(d: int, multiplicity: int, transverse_degree: int) -> int:
    ceiling = (multiplicity + 1) // 2
    for order in range(ceiling):
        if resonance(d, multiplicity, transverse_degree, order) == 0:
            return order
    return ceiling
for multiplicity in range(1, 8):
    for transverse_degree in range(7):
        assert root_weight(8, multiplicity, transverse_degree) == (multiplicity + 1) // 2

def integer_partitions(total: int, maximum: int | None=None):
    if total == 0:
        yield ()
        return
    maximum = total if maximum is None else min(total, maximum)
    for first in range(maximum, 0, -1):
        for tail in integer_partitions(total - first, first):
            yield ((first,) + tail)
all_degree_eight = set(integer_partitions(8))
inherited = {(1, 1, 1, 1, 1, 1, 1, 1), (2, 1, 1, 1, 1, 1, 1)}
fixed_ledger = {(7, 1), (6, 2), (6, 1, 1), (5, 3), (5, 2, 1), (4, 4), (4, 3, 1), (4, 2, 2), (3, 3, 2)}
four_root_w5 = {(4, 2, 1, 1), (3, 2, 2, 1)}
universal_w6 = {(5, 1, 1, 1), (4, 1, 1, 1, 1), (3, 3, 1, 1), (3, 2, 1, 1, 1), (2, 2, 1, 1, 1, 1)}
five_root_w5 = {(2, 2, 2, 1, 1)}
four_double = {(2, 2, 2, 2)}
weight_excess = {(3, 1, 1, 1, 1, 1)}
pure_power = {(8,)}
assert all_degree_eight == inherited | fixed_ledger | four_root_w5 | universal_w6 | five_root_w5 | four_double | weight_excess | pure_power
assert sum(((m + 1) // 2 for m in (3, 1, 1, 1, 1, 1))) == 7
fixed_partitions: dict[str, tuple[int, ...]] = {'71': (7, 1), '62': (6, 2), '53': (5, 3), '44': (4, 4), '611': (6, 1, 1), '521': (5, 2, 1), '431': (4, 3, 1), '422': (4, 2, 2), '332': (3, 3, 2)}

def check_fixed_partition(name: str, multiplicities: tuple[int, ...]) -> None:
    roots = (x, y, x - y)
    f = sp.prod((root ** m for root, m in zip(roots, multiplicities)))
    root_factor = sp.prod((root ** ((m + 1) // 2) for root, m in zip(roots, multiplicities)))
    weight = sp.Poly(root_factor, x, y).total_degree()
    for e in range(weight, 7):
        h_degree = e - weight
        hs = tuple(sp.symbols(f'p{name}_e{e}_h0:{h_degree + 1}'))
        g = sp.expand(root_factor * homogeneous(hs, h_degree))
        q_degree = 2 * e - 8
        qs = tuple(sp.symbols(f'p{name}_e{e}_q0:{q_degree + 1}'))
        q = homogeneous(qs, q_degree)
        coefficients = sp.Poly(schur_face(f, g, q), x, y).coeffs()
        constant_basis = sp.groebner(coefficients, *qs, *hs, order='grevlex')
        h_powers = [radical_power(constant_basis, h, 10) for h in hs]
        if all((power is not None for power in h_powers)):
            continue
        if name == '422' and e == 6:
            continue
        potential, higher_symbols = weighted_face(f, g, q, e, f'p{name}_e{e}')
        full_basis = sp.groebner(sp.Poly(bordered(potential), x, y, z).coeffs(), *qs, *higher_symbols, *hs, order='grevlex')
        assert all((radical_power(full_basis, h, 12) is not None for h in hs))
for partition_name, partition in fixed_partitions.items():
    check_fixed_partition(partition_name, partition)
f422 = x ** 4 * y ** 2 * (x - y) ** 2
root422 = x ** 2 * y * (x - y)
a422 = sp.symbols('a422')
g422 = -a422 * x ** 2 * y * (x - y) * (x + 2 * y) * (3 * x - 2 * y)
q422 = a422 ** 2 * (sp.Rational(33, 8) * x ** 4 + 23 * x ** 3 * y - 21 * x ** 2 * y ** 2 - 4 * x * y ** 3 + 2 * y ** 4)
r422 = sp.symbols('r422_0:3')
s422 = sp.symbols('s422')
potential422 = f422 + z * g422 + z ** 2 * q422 / 2 + z ** 3 * homogeneous(tuple(r422), 2) / 6 + z ** 4 * s422 / 24
basis422 = sp.groebner(sp.Poly(bordered(potential422), x, y, z).coeffs(), *r422, s422, a422, order='grevlex')
assert [sp.factor(p.as_expr()) for p in basis422.polys] == [a422 ** 3, r422[0], r422[1], r422[2], s422]
result = {'scope': 'degree-eight scalar root census and two-/three-root charts', 'status': 'closed', 'partition_count': len(all_degree_eight)}
output = ARTIFACT_DIR / 'hc4_degree8_partition_census.json'
output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print(json.dumps(result, indent=2))
