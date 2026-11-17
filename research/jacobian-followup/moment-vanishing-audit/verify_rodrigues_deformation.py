#!/usr/bin/env python3
"""Exact finite replay of Rodrigues jets and the local rigidity certificate.

The all-order deductions use the written jet and filtered-ring arguments;
this is not formal or independent specialist verification.
"""
from pathlib import Path
import hashlib
import json
import math

from conformal_euler_cubic import sphere_average
from rodrigues_deformation import (
    S, z, w, t, a, b, F0, G1, G2, BASIS,
    marked_monomial, marked, linear_data, implicit_jets,
)


def main():
    root = Path(__file__).resolve().parent
    direct_checks = 0
    for power in range(9):
        for i in range(4):
            for j in range(4-i):
                for k in range(4-i-j):
                    actual = sphere_average(z**i*w**j*t**k*F0**power)
                    assert actual == marked_monomial(power, i, j, k)
                    direct_checks += 1

    jacobian, pivots, inverse, coordinates = linear_data()
    assert jacobian.rank() == 10
    assert coordinates[:, 10:14].rank() == 4
    assert coordinates[:, 10:].rank() == 6
    assert all(marked(f, n-1) == 0 for n in range(11, 16) for f in BASIS)
    for r in range(1, 4):
        actual = math.factorial(10*r)//math.factorial(9*r)*marked(w**(3*r), 9*r)
        expected = (-S.Rational(1, 2))**(3*r)*S.factorial(10*r)/S.factorial2(18*r+1)
        assert actual == expected and actual != 0

    def checkpoint(prefix):
        def save(degree, jet, tail):
            path = root/(prefix+str(degree)+'_checkpoint.json')
            data = {'order': degree, 'jet': str(jet),
                    'obstructions': {str(n): str(f) for n, f in tail.items()}}
            if path.exists():
                previous = json.loads(path.read_text())
                assert previous == data, 'changed exact checkpoint '+path.name
            else:
                path.write_text(json.dumps(data, indent=2)+'\n')
            print(prefix+str(degree)+': '+str(len(tail))+' nonzero tail residuals', flush=True)
        return save

    jets, tails = implicit_jets(5, checkpoint('rodrigues_order'))
    b_jets, b_tails = implicit_jets(6, checkpoint('rodrigues_b_axis_order'), initial=b*G2)
    assert tails[2] == tails[3] == {}
    assert set(tails[4]) == {11, 12}
    assert set(tails[5]) == {11, 12, 13, 14, 15}
    assert all(not b_tails[r] for r in range(2, 6))
    assert b_tails[6] == {12: -S.Rational(27910144, 52441021269)*b**6}
    assert S.expand(jets[2]-(-S.Rational(400, 5733)*a*a*t*w
                             +S.Rational(25, 9)*a*a*w**3
                             -S.Rational(190, 1911)*a*a*z
                             +S.Rational(40, 231)*a*b*(t*t-w*z)
                             -S.Rational(25, 3)*a*b*t*w*w
                             +5*b*b*w*(t*t-w*z/2)+b*b*w/3)) == 0

    def weight_part(poly, weight):
        return S.Add(*(coefficient*a**i*b**j for (i, j), coefficient
                       in S.Poly(S.expand(poly), a, b).terms() if 3*i+2*j == weight))

    actual_p = weight_part(tails[4][11]+tails[5][11], 11)
    actual_q = weight_part(tails[4][12]+tails[5][12]+b_tails[6][12], 12)
    P = a*b*(1265477*a*a-4866582*b**3)
    Q = 5242521625*a**4-174930621075*a*a*b**3+112837775358*b**6
    assert S.expand(actual_p+S.Rational(64000, 123934268957499)*P) == 0
    assert S.expand(actual_q+S.Rational(4096, 868407422585195493)*Q) == 0
    ratio = S.Rational(4866582, 1265477)
    separation = 5242521625*ratio**2-174930621075*ratio+112837775358
    assert separation != 0
    basis = S.groebner([P, Q], a, b, order='lex')
    leading = [g.LM(order='lex').exponents for g in basis.polys]
    standard = [(i, j) for i in range(10) for j in range(12)
                if not any(i >= u and j >= v for u, v in leading)]
    assert leading == [(4, 0), (3, 1), (2, 4), (1, 7), (0, 10)]
    assert len(standard) == 22
    assert max(i+j for i, j in standard) == 9
    assert max(3*i+2*j for i, j in standard) == 18
    assert (0, 9) in standard

    report = {
        'claims': ['SIC2B33JET20261116', 'SIC2B33LOCAL20261116'],
        'status': 'PASS',
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'arithmetic_sha256': hashlib.sha256((root/'rodrigues_deformation.py').read_bytes()).hexdigest(),
        'sympy_version': S.__version__,
        'direct_sphere_monomial_checks': direct_checks,
        'first_ten_jacobian_rank': 10,
        'pivot_basis_indices_zero_based': pivots,
        'orbit_tangent_dimension': 4, 'all_moment_tangent_dimension': 6,
        'third_order_obstruction': 'zero through moment30; higher moments vanish by the finite-jet lemma',
        'fourth_order_residuals': {str(n): str(f) for n, f in tails[4].items()},
        'fifth_order_residuals': {str(n): str(f) for n, f in tails[5].items()},
        'sixth_order_b_axis_residuals': {str(n): str(f) for n, f in b_tails[6].items()},
        'parameter_weights': {'a': 3, 'b': 2},
        'primitive_weighted_initial_R11': str(P),
        'primitive_weighted_initial_R12': str(Q),
        'ratio_separation': str(separation),
        'leading_monomial_exponents_lex': leading,
        'standard_monomial_exponents': standard,
        'first_twelve_transverse_length': 22,
        'first_twelve_transverse_nilpotence_index': 10,
        'full_moment_transverse_length_upper_bound': 22,
        'full_moment_local_generating_prefix_upper_bound': 90,
        'conclusion': 'The reduced local pure locus is exactly the scalar/orthogonal Rodrigues orbit.',
        'boundary': 'The first-twelve local transverse algebra has length22; the full moment algebra '
                    'is a quotient, with exact length not computed. Other semistable components '
                    'away from this orbit and the full (3,3) Image theorem remain open. '
                    'No independent specialist, external or formal verification.',
    }
    (root/'rodrigues_deformation_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in
                      {'standard_monomial_exponents', 'fifth_order_residuals'}}, indent=2))


if __name__ == '__main__':
    main()
