#!/usr/bin/env python3
"""Finite exact algebra checks for the written odd Hopf/Image transfer.

Does not reverify the retained all-order Gaussian cubic classification.
"""
from itertools import product
from pathlib import Path
import hashlib
import json
import math

import sympy as S


W, V, Z, Y = S.symbols('W V Z Y')
x, y, z = S.symbols('x y z')
s = W*Z + V*Y
h = (W*Y + V*Z, (W*Y - V*Z)/S.I, W*Z - V*Y)


def scalar_contraction(expr):
    value = S.Integer(0)
    for (a, b, c, d), coefficient in S.Poly(S.expand(expr), W, V, Z, Y).terms():
        assert a+b == c+d, 'scalar input must be balanced'
        if (a, b) == (c, d):
            value += coefficient*math.factorial(c)*math.factorial(d)
    return S.expand(value)


def sphere_monomial(a, b, c):
    if any(n % 2 for n in (a, b, c)):
        return S.Integer(0)
    return S.prod(S.factorial2(n-1) for n in (a, b, c))/S.factorial2(a+b+c+1)


def main():
    root = Path(__file__).resolve().parent
    assert S.expand(sum(v*v for v in h)-s*s) == 0
    powers = [(a, b, 3-a-b) for a in range(4) for b in range(4-a)]
    biforms = [S.Poly(S.expand(h[0]**a*h[1]**b*h[2]**c), W, V, Z, Y)
               for a, b, c in powers]
    basis = [W**(3-i)*V**i*Z**(3-j)*Y**j for i in range(4) for j in range(4)]
    substitution = S.Matrix([[p.coeff_monomial(v) for p in biforms] for v in basis])
    assert substitution.shape == (16, 10) and substitution.rank() == 10

    hopf_checks = 0
    for degree in range(5):
        for a in range(degree+1):
            for b in range(degree-a+1):
                for c in range(degree-a-b+1):
                    expr = h[0]**a*h[1]**b*h[2]**c*s**(degree-a-b-c)
                    actual = scalar_contraction(expr)/math.factorial(degree+1)
                    assert S.simplify(actual-sphere_monomial(a, b, c)) == 0
                    hopf_checks += 1

    polar_checks = 0
    for a, b, c in product(range(9), repeat=3):
        total = a+b+c
        if total > 8:
            continue
        gaussian = (S.Integer(0) if any(n % 2 for n in (a, b, c))
                    else S.prod(S.factorial2(n-1) for n in (a, b, c)))
        radial = 2**S.Rational(total, 2)*S.gamma(S.Rational(total+3, 2))/S.gamma(S.Rational(3, 2))
        assert S.simplify(gaussian-radial*sphere_monomial(a, b, c)) == 0
        polar_checks += 1

    # For a monomial input, compare every possible output coefficient
    # with contraction after multiplication by the detecting dual monomial.
    # Linearity makes this a check of the coefficient identity on a basis.
    coefficient_checks = 0
    negative_degree_checks = 0
    for a in range(6):
        for b in range(6):
            for r in range(a+1):
                for p in range(b+1):
                    dual = (r, a-r)
                    coordinate = (p, b-p)
                    survives = all(u <= v for u, v in zip(dual, coordinate))
                    if a > b:
                        assert not survives
                        negative_degree_checks += 1
                        continue
                    for g0 in range(b-a+1):
                        gamma = (g0, b-a-g0)
                        matches = all(v == u+g for u, v, g in zip(dual, coordinate, gamma))
                        lhs = (math.prod(math.factorial(v)//math.factorial(v-u)
                                         for u, v in zip(dual, coordinate))
                               if survives and matches else 0)
                        rhs = (S.Rational(math.prod(math.factorial(v) for v in coordinate),
                                          math.prod(math.factorial(g) for g in gamma))
                               if matches else 0)
                        assert lhs == rhs
                        coefficient_checks += 1

    A = x+S.I*y
    norm = x*x+y*y+z*z
    G = A*z*z
    H3, H1 = G-norm*A/5, A/5
    laplace = lambda f: S.expand(sum(S.diff(f, v, 2) for v in (x, y, z)))
    assert laplace(H3) == laplace(H1) == 0
    assert S.expand(H3+norm*H1-G) == 0 and H3 != 0 and H1 != 0
    P = 2*W*Y*h[2]**2
    assert S.expand(G.subs(dict(zip((x, y, z), h)), simultaneous=True)-P) == 0
    for m in range(1, 9):
        assert scalar_contraction(P**m) == 0
    sharp_values = {}
    for b in range(1, 7):
        actual = scalar_contraction(P**b*(2*V*Z)**b)
        angular = sum((-1)**j*S.binomial(b, j)/S.Integer(2*b+2*j+1)
                      for j in range(b+1))
        assert actual == math.factorial(4*b+1)*angular and actual > 0
        sharp_values[str(b)] = str(actual)
        assert scalar_contraction(P**(b+1)*(2*V*Z)**b) == 0

    receipt = {
        'claim': 'SIC2B33ODD20261116',
        'status': 'PASS',
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'sympy_version': S.__version__,
        'cubic_hopf_substitution_rank': 10,
        'hopf_sphere_basis_checks': hopf_checks,
        'gaussian_polar_monomial_checks': polar_checks,
        'coefficient_extraction_checks': coefficient_checks,
        'negative_output_degree_checks': negative_degree_checks,
        'sharp_mixed_slice_values_at_m_equal_b': sharp_values,
        'boundary': 'Finite exact algebra only. The written all-order transfer imports '
                    'G3HOMOGENEOUS20261114; no independent specialist, external or formal review.',
    }
    (root/'odd_hopf_cubic_transfer_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
