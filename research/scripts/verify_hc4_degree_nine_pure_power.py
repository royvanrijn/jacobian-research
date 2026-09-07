#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
import sympy as sp
ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_DIR = ROOT / 'artifacts' / 'generated-results'
ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
x, y, z = sp.symbols('x y z')

def hom(v, d):
    return sum((v[i] * x ** (d - i) * y ** i for i in range(d + 1)))

def bordered(c):
    g = sp.Matrix([sp.diff(c, v) for v in (x, y, z)])
    H = sp.hessian(c, (x, y, z))
    return sp.expand((g.T * H.adjugate() * g)[0])

def binary_border(f):
    g = sp.Matrix([sp.diff(f, x), sp.diff(f, y)])
    H = sp.hessian(f, (x, y))
    return sp.Matrix.vstack(sp.Matrix.hstack(sp.zeros(1, 1), g.T), sp.Matrix.hstack(g, H))

def schur(f, g, q):
    B = binary_border(f)
    c = sp.Matrix([g, sp.diff(g, x), sp.diff(g, y)])
    return sp.expand(q * B.det() - (c.T * B.adjugate() * c)[0])

def resonance(d, m, e, n):
    return d * d * m + d * d * n * n - 2 * d * e * m * n - 2 * d * e * m - d * m * m - d * m + e * e * m * m + 2 * e * m * m + m * m

def root_weight(d, m, e):
    for n in range((m + 1) // 2):
        if resonance(d, m, e, n) == 0:
            return n
    return (m + 1) // 2
for e in range(5):
    assert root_weight(9, 9, e) == e
a5 = sp.symbols('a5')
q5 = tuple(sp.symbols('e5_q0:2'))
potential = x ** 9 + a5 * x ** 5 * z + z ** 2 * hom(q5, 1) / 2
J = sp.Poly(bordered(potential), x, y, z)
assert J.coeff_monomial(x ** 16 * z ** 2) == -81 * q5[1] ** 2
a6, b6 = sp.symbols('a6 b6')
q6 = tuple(sp.symbols('e6_q0:4'))
r6 = sp.symbols('e6_r')
g6 = x ** 5 * (a6 * x + b6 * y)
constant = sp.Poly(schur(x ** 9, g6, hom(q6, 3)), x, y)
assert constant.coeff_monomial(x ** 26) == 81 * b6 ** 2
potential = x ** 9 + a6 * x ** 6 * z + z ** 2 * hom(q6, 3) / 2 + z ** 3 * r6 / 6
J = sp.Poly(bordered(potential), x, y, z)
assert J.coeff_monomial(x ** 16 * y ** 4 * z ** 2) == -486 * q6[3] ** 2
assert sp.factor(J.coeff_monomial(x ** 18 * y ** 2 * z ** 2).subs(q6[3], 0) + 243 * q6[2] ** 2) == 0
assert sp.factor(J.coeff_monomial(x ** 20 * z ** 2).subs({q6[3]: 0, q6[2]: 0}) + 81 * q6[1] ** 2) == 0
a7, b7, c7 = sp.symbols('a7 b7 c7')
q7 = tuple(sp.symbols('e7_q0:6'))
r7 = tuple(sp.symbols('e7_r0:4'))
s7 = tuple(sp.symbols('e7_s0:2'))
g7 = x ** 5 * (a7 * x ** 2 + b7 * x * y + c7 * y ** 2)
constant = sp.Poly(schur(x ** 9, g7, hom(q7, 5)), x, y)
assert constant.coeff_monomial(x ** 28) == 81 * b7 ** 2
assert constant.coeff_monomial(x ** 27 * y) == 324 * b7 * c7
assert constant.coeff_monomial(x ** 26 * y ** 2) == 324 * c7 ** 2
potential = x ** 9 + a7 * x ** 7 * z + z ** 2 * hom(q7, 5) / 2 + z ** 3 * hom(r7, 3) / 6 + z ** 4 * hom(s7, 1) / 24
J = sp.Poly(bordered(potential), x, y, z)
zero = {}
checks = [(q7[5], x ** 16 * y ** 8 * z ** 2, -1215), (q7[4], x ** 18 * y ** 6 * z ** 2, -810), (q7[3], x ** 20 * y ** 4 * z ** 2, -486), (q7[2], x ** 22 * y ** 2 * z ** 2, -243), (q7[1], x ** 24 * z ** 2, -81), (r7[3], x ** 16 * y ** 4 * z ** 4, -sp.Rational(405, 4)), (r7[2], x ** 18 * y ** 2 * z ** 4, -54), (r7[1], x ** 20 * z ** 4, -sp.Rational(81, 4)), (s7[1], x ** 16 * z ** 6, -sp.Rational(9, 4))]
for variable, monomial, coefficient in checks:
    actual = sp.factor(J.coeff_monomial(monomial).subs(zero))
    assert sp.factor(actual - coefficient * variable ** 2) == 0
    zero[variable] = 0
result = {'scope': 'pure ninth-power synchronized scalar packet', 'status': 'all e=0,...,7 weighted faces are fixed cylinders', 'surviving_terms': 'depend only on x,z'}
(ARTIFACT_DIR / 'hc4_degree9_pure_power.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
