#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
import sympy as sp
ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_DIR = ROOT / 'artifacts' / 'generated-results'
ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
x, y, z = sp.symbols('x y z')

def border(f):
    g = sp.Matrix([sp.diff(f, x), sp.diff(f, y)])
    H = sp.hessian(f, (x, y))
    return sp.Matrix.vstack(sp.Matrix.hstack(sp.zeros(1, 1), g.T), sp.Matrix.hstack(g, H))

def schur(f, g, q):
    B = border(f)
    c = sp.Matrix([g, sp.diff(g, x), sp.diff(g, y)])
    return sp.expand(q * B.det() - (c.T * B.adjugate() * c)[0])

def bordered(c):
    g = sp.Matrix([sp.diff(c, v) for v in (x, y, z)])
    H = sp.hessian(c, (x, y, z))
    return sp.expand((g.T * H.adjugate() * g)[0])

def hom(v, d):
    return sum((v[i] * x ** (d - i) * y ** i for i in range(d + 1)))

def solve_q(f, g, degree):
    qs = tuple(sp.symbols(f'q0:{degree + 1}'))
    q = hom(qs, degree)
    solutions = sp.solve(sp.Poly(schur(f, g, q), x, y).coeffs(), qs, dict=True, simplify=False)
    assert len(solutions) == 1
    return sp.expand(q.subs(solutions[0]))

def radical_power(G, p, limit=15):
    for k in range(1, limit + 1):
        if G.reduce(sp.expand(p ** k))[1] == 0:
            return k
    return None
roots = (x, y, x - y)
a = sp.symbols('a')
record = {}
ray_data = {'6+2+1': ((6, 2, 1), (sp.Rational(49, 9), sp.Rational(14, 3), 1)), '5+2+2': ((5, 2, 2), (sp.Rational(-14, 25), -1, 1)), '4+4+1': ((4, 4, 1), (1, 2, 1)), '4+3+2': ((4, 3, 2), (sp.Rational(-9, 8), sp.Rational(-3, 4), 1))}
for name, (multiplicities, vector) in ray_data.items():
    f = sp.prod((root ** m for root, m in zip(roots, multiplicities)))
    root_factor = sp.prod((root ** ((m + 1) // 2) for root, m in zip(roots, multiplicities)))
    H = vector[0] * x ** 2 + vector[1] * x * y + vector[2] * y ** 2
    g = sp.expand(a * root_factor * H)
    q = solve_q(f, g, 5)
    r = tuple(sp.symbols(f"{name.replace('+', '_')}_r0:4"))
    s = tuple(sp.symbols(f"{name.replace('+', '_')}_s0:2"))
    potential = f + z * g + z ** 2 * q / 2 + z ** 3 * hom(r, 3) / 6 + z ** 4 * hom(s, 1) / 24
    G = sp.groebner(sp.Poly(bordered(potential), x, y, z).coeffs(), *r, *s, a, order='grevlex')
    expected = [a ** 3, *r, *s]
    actual = [sp.factor(p.as_expr()) for p in G.polys]
    assert actual == expected
    record[name] = {'basis': [str(v) for v in actual], 'transverse_radical': 'a=0'}
f = x ** 3 * y ** 3 * (x - y) ** 3
u, v = sp.symbols('u v')
H = u * x ** 2 - 2 * (u + v) * x * y + v * y ** 2
g = x * y * (x - y) * H
q = solve_q(f, g, 1)
G = sp.groebner(sp.Poly(bordered(f + z * g + z ** 2 * q / 2), x, y, z).coeffs(), u, v, order='grevlex')
actual = [sp.factor(p.as_expr()) for p in G.polys]
assert actual == [u ** 3, u * v * (u + v), v ** 3]
record['3+3+3@e=5'] = {'basis': [str(w) for w in actual], 'transverse_radical': 'u=v=0'}
a6, r6 = sp.symbols('a6 r6')
g = a6 * x ** 2 * y ** 2 * (x - y) ** 2
q = solve_q(f, g, 3)
G = sp.groebner(sp.Poly(bordered(f + z * g + z ** 2 * q / 2 + z ** 3 * r6 / 6), x, y, z).coeffs(), r6, a6, order='grevlex')
actual = [sp.factor(p.as_expr()) for p in G.polys]
expected = [a6 ** 2 * r6, a6 ** 3 - sp.Rational(81, 2) * r6, r6 ** 2]
assert all((sp.expand(left - right) == 0 for left, right in zip(actual, expected)))
assert radical_power(G, a6) == 5 and radical_power(G, r6) == 2
record['3+3+3@e=6'] = {'basis': [str(w) for w in actual], 'transverse_radical': 'a6=r6=0'}
result = {'status': 'all six three-root survivors closed', 'record': record}
(ARTIFACT_DIR / 'hc4_degree9_three_root_survivors.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
