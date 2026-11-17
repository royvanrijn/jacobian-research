#!/usr/bin/env python3
"""Exact bounded controls for FIXED_MAP_PLANE_SECTIONS.md; SymPy required."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import platform

import sympy as sp


def fixed_map(x, y, z):
    x = sp.cancel(x/2)
    t = sp.cancel(1+x*y)
    q = sp.cancel(t*t*z + sp.Rational(3, 8)*y*y*(1+3*t))
    first = -t*q
    second = (y+8*x*q-sp.Rational(40, 9)*t*q
              -sp.Rational(32, 3)*t*t*x*x*q**4
              +sp.Rational(80, 9)*t*t*x**3*q**5)
    third = (x*(5-3*t)-sp.Rational(8, 3)*x**3*z
             +sp.Rational(16, 3)*(x*q)**4-sp.Rational(16, 3)*(x*q)**5)
    return tuple(sp.cancel(v) for v in (first, second, third))


def mobius(n):
    answer, p = 1, 2
    while p*p <= n:
        if n % p == 0:
            n //= p
            answer = -answer
            if n % p == 0:
                return 0
        p += 1
    if n > 1:
        answer = -answer
    return answer


def run():
    pi, b, c, s, u, alpha, beta, delta = sp.symbols('pi b c s u alpha beta delta')
    q = sp.Rational
    controls = {}
    ys = b-q(40, 9)*u
    section_c0 = fixed_map(sp.Integer(0), ys, -u-q(3, 2)*ys**2)
    assert section_c0 == (u, b, 0)
    controls['c_zero_section_coordinates'] = 3
    zs = q(3, 64)*b**3*(10/b+q(16, 3)*(3*b/16)**4-q(16, 3)*(3*b/16)**5-c)
    section_u0 = fixed_map(4/b, -b/2, zs)
    assert section_u0 == (0, b, c)
    assert fixed_map(c, sp.Integer(0), sp.Integer(0)) == (0, 0, c)
    controls['u_zero_section_and_remaining_line_coordinates'] = 6

    E = (pi**5*s**5-q(3, 2)*pi**4*s**4+q(3, 2)*pi*s**3
         +(-q(5, 4)*pi-q(9, 32)*b)*s**2+q(9, 16)*s-q(9, 32)*c)
    restriction = sp.expand(E.subs(c, alpha*pi+beta*b+delta))
    N = sp.expand(restriction.subs(b, 0))
    assert sp.expand(restriction-N+q(9, 32)*b*(s*s+beta)) == 0
    resultant = sp.Poly(sp.resultant(N, s*s+beta, s), pi)
    assert resultant.degree() == 10 and resultant.LC() == beta**5
    at_beta0 = sp.factor(resultant.as_expr().subs(beta, 0))
    assert sp.expand(at_beta0-q(81, 1024)*(alpha*pi+delta)**2) == 0
    controls['structural_generic_plane_resultant_controls'] = 3
    finite_planes = 0
    for a in (-1, 0, 1):
        for v in (-1, 0, 1):
            for d in (-1, 0, 1):
                Ns = N.subs({alpha: a, beta: v, delta: d})
                common = sp.Poly(sp.gcd(Ns, s*s+v), s, domain=sp.QQ.frac_field(pi))
                if a == v == d == 0:
                    assert common.degree() == 1
                else:
                    assert common.degree() == 0
                finite_planes += 1
    controls['finite_nontautological_plane_gcd_controls'] = finite_planes
    assert E.subs({c: 0, s: 0}) == 0
    assert sp.diff(E, s).subs(s, 0) == q(9, 16)
    assert sp.diff(E, c) == -q(9, 32)
    assert sp.Poly(E, s).degree() == 5
    for a in (-2, -1, 1, 2):
        assert sp.Poly(E.subs(pi, a), s).degree() == 5
    controls['special_plane_and_scalar_parameter_controls'] = 8

    C = (q(32, 9)*pi**5*s**5-q(16, 3)*pi**4*s**4
         +q(16, 3)*pi*s**3-q(40, 9)*pi*s*s+2*s-b*s*s)
    assert sp.expand(E.subs(c, C)) == 0
    assert sp.expand(sp.diff(C, s)-q(32, 9)*sp.diff(E, s)) == 0
    controls['root_incidence_identities'] = 2

    primitive = []
    for B in (1, 2, 3, 5, 10, 20):
        plane = {(den, x, y) for den in range(1, B+1)
                 for x in range(-B, B+1) for y in range(-B, B+1)
                 if gcd(den, x, y) == 1}
        line = {(den, x) for den in range(1, B+1) for x in range(-B, B+1)
                if gcd(den, x) == 1}
        plane_formula = sum(mobius(d)*(B//d)*(2*(B//d)+1)**2 for d in range(1, B+1))
        line_formula = sum(mobius(d)*(B//d)*(2*(B//d)+1) for d in range(1, B+1))
        assert plane_formula == len(plane) and line_formula == len(line)
        union = {(den, 0, x, y) for den, x, y in plane}
        union |= {(den, x, y, 0) for den, x, y in plane}
        assert len(union) == 2*len(plane)-len(line)
        primitive.append({'height': B, 'one_plane': len(plane),
                          'intersection_line': len(line), 'two_plane_union': len(union)})
    controls['primitive_plane_union_height_controls'] = len(primitive)
    here = Path(__file__).resolve().parent
    root = here.parents[1]
    proof = here/'FIXED_MAP_PLANE_SECTIONS.md'
    inputs = ['unpublished/fixed-map-hasse-failures/sections/02-fixed-map.tex',
              'research/arithmetic/sources/buggenhout-et-al-thin-projective-2506.13471v2.pdf']
    return {'status': 'pass', 'python': platform.python_version(), 'sympy': sp.__version__,
            'command': '.venv/bin/python research/arithmetic/verify_fixed_map_plane_sections.py --write-receipt',
            'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
            'proof_sha256': sha256(proof.read_bytes()).hexdigest(),
            'input_sha256': {p: sha256((root/p).read_bytes()).hexdigest() for p in inputs},
            'checks': controls, 'primitive_counts': primitive,
            'exact_structural_resultant': {'leading_pi_degree': 10,
                                          'leading_pi_coefficient': 'beta^5',
                                          'beta_zero_resultant': '81*(alpha*pi+delta)^2/1024'},
            'ambient_lower_coefficient': '8/zeta(3)',
            'scope': ('Exact symbolic section and resultant identities, 27 bounded plane gcd controls, '
                      'and primitive plane-union enumeration through height20. These do not prove '
                      'the unrestricted irreducibility, projective thin-set theorem, or asymptotic.'),
            'formal_verification': False, 'independent_whole_proof_replay': False,
            'external_review': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-receipt', action='store_true')
    args = parser.parse_args()
    result = run()
    output = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.write_receipt:
        Path(__file__).resolve().with_name('fixed_map_plane_sections_checks.json').write_text(output, encoding='utf-8')
    print(output, end='')


if __name__ == '__main__':
    main()
