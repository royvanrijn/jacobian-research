#!/usr/bin/env python3
"""Finite exact Haar source/period controls; not a general holonomic solver."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import sys

import verify_haar_shear_witnesses as S

ROOT = Path(__file__).resolve().parent
ONE = {(0, 0, 0, 0): Q(1)}
X, Y, Z, W = S.VARS


def weight(poly):
    return {m: c*(m[0]+m[1]-m[2]-m[3]) for m, c in poly.items()
            if m[0]+m[1]-m[2]-m[3]}


def monomials(degree, variables=4):
    return [a for a in product(range(degree+1), repeat=variables)
            if sum(a) <= degree]


def source_controls():
    cases = {
        'linear': X,
        'mixed': S.add(X, S.power(Z, 2)),
        'octahedral_highest': {(8, 0, 0, 0): Q(1), (4, 4, 0, 0): Q(14),
                               (0, 8, 0, 0): Q(1)},
        'icosahedral_highest': {(11, 1, 0, 0): Q(1), (6, 6, 0, 0): Q(11),
                                (1, 11, 0, 0): Q(-1)},
    }
    count = 0
    for name, p in cases.items():
        pm = ONE
        for m in range(4):
            for alpha in monomials(2):
                q = {alpha: Q(1)}
                qp = S.mul(q, pm)
                assert S.haar(S.mul(S.DET, qp)) == S.haar(qp)
                count += 1
                for derivation in [S.raise_, S.lower, weight]:
                    value = S.haar(derivation(qp))
                    assert value == 0, (name, m, alpha)
                    rhs = S.haar(S.mul(derivation(q), pm))
                    if m:
                        rhs += m*S.haar(S.mul(q, S.mul(derivation(p), S.power(p, m-1))))
                    assert rhs == value
                    count += 1
            pm = S.mul(pm, p)
    # In the lowering equation at P=x, Q=y, m=1 the plus sign is essential.
    lhs = S.haar(S.mul(S.lower(Y), X))
    rhs = S.haar(S.mul(Y, S.lower(X)))
    assert lhs == Q(1, 2) and rhs == Q(-1, 2)
    assert lhs+rhs == 0 and lhs-rhs == 1
    return {'exact_source_coefficient_identities': count,
            'incorrect_tangent_sign_detected': True,
            'wrong_sign_residual': str(lhs-rhs)}


def add3(*polys):
    out = defaultdict(Q)
    for p in polys:
        for a, c in p.items():
            out[a] += c
    return {a: c for a, c in out.items() if c}


def mul3(p, q):
    out = defaultdict(Q)
    for a, c in p.items():
        for b, d in q.items():
            out[tuple(i+j for i, j in zip(a, b))] += c*d
    return {a: c for a, c in out.items() if c}


def scale3(p, c):
    return {a: c*d for a, d in p.items() if c*d}


def sphere3(p):
    # u=conjugate(v), tau real; uv+tau^2=1 on the unit two-sphere.
    out = Q(0)
    for (a, b, c), coefficient in p.items():
        if a == b and c % 2 == 0:
            beta = sum((Q((-1)**j*comb(a, j), c+2*j+1)
                        for j in range(a+1)), Q(0))
            out += coefficient*beta
    return out


def normalizer_control():
    one = {(0, 0, 0): Q(1)}
    u = {(1, 0, 0): Q(1)}
    U, V, T2 = {(2, 0, 0): Q(1)}, {(0, 2, 0): Q(1)}, {(0, 0, 2): Q(1)}
    U2 = mul3(U, U)
    inside = add3(scale3(one, 2), scale3(mul3(T2, add3(U2, scale3(U, 2), scale3(one, 2))), -1))
    correction = mul3(mul3(T2, add3(U, scale3(one, 2))), inside)
    pn = mul3(add3(one, U), add3(V, scale3(correction, -1)))
    ppower = one
    checked = 0
    for m in range(7):
        assert sphere3(ppower) == int(m == 0)
        for s in range(1, 4):
            marker = {(2*s, 0, 0): Q(1)}
            expected = (Q(4**(2*m)*factorial(2*m)**2, factorial(4*m+1))
                        *comb(m-1, s-1) if 1 <= s <= m else Q(0))
            assert sphere3(mul3(marker, ppower)) == expected, (m, s)
            checked += 1
        ppower = mul3(ppower, pn)
    # Exponential coefficients of G=1+F_U, checked beyond the direct integrals.
    a = [Q(4**(2*m)*factorial(2*m)**2,
           factorial(4*m+1)*factorial(m)) for m in range(17)]
    for m in range(16):
        assert (Q(m)+Q(3, 4))*(Q(m)+Q(5, 4))*a[m+1] == (Q(m)+Q(1, 2))*a[m]
    f = a.copy()
    f[0] = Q(0)
    # M=1/16+(-t/16-3t^2/2)D+(2t^2-t^3)D^2+t^3D^3.
    for m in range(16):
        diagonal = (Q(1, 16)-Q(m, 16)+2*m*(m-1)+m*(m-1)*(m-2))*f[m]
        previous = (Q(-3, 2)*(m-1)-(m-1)*(m-2))*f[m-1] if m else Q(0)
        assert diagonal+previous == 0, m
    for v in range(9):
        infinity = -Q(3, 2)*v-v*(v-1)
        origin = Q(1, 16)-Q(v, 16)+2*v*(v-1)+v*(v-1)*(v-2)
        assert infinity == -v*(Q(v)+Q(1, 2))
        assert origin == (v-1)*(Q(v)-Q(1, 4))*(Q(v)+Q(1, 4))
        assert (infinity == 0) == (v == 0)
        assert (origin == 0) == (v == 1)
    assert f[0] == 0 and f[1] == Q(8, 15)
    return {'direct_marked_sphere_integrals': checked,
            'direct_pure_sphere_integrals': 7,
            'scalar_operator_coefficient_checks': 16,
            'infinity_indicial': '-v*(v+1/2)',
            'origin_indicial': '(v-1)*(v-1/4)*(v+1/4)',
            'polynomial_degree_bound': 0, 'first_nonzero_jet_index': 1,
            'first_nonzero_coefficient': str(f[1]),
            'finite_nonpolynomiality_certificate': True,
            'entireness_proof': 'compactness in the written note',
            'actual_period_rank_lower_bound': 2,
            'exceptional_marker_constructed': False}


def exceptional_families():
    cases = {}
    for name, n, coeffs in [('2O', 8, {0: Q(1), 4: Q(14), 8: Q(1)}),
                             ('2I', 12, {1: Q(1), 6: Q(11), 11: Q(-1)})]:
        p = {(n-i, i, 0, 0): c for i, c in coeffs.items()}
        q = {(0, 0, n-i, i): c for i, c in coeffs.items()}
        assert S.compact_star(p) == q and not S.lower(q)
        ppowers, qpowers = [ONE], [ONE]
        for _ in range(6):
            ppowers.append(S.mul(ppowers[-1], p))
            qpowers.append(S.mul(qpowers[-1], q))
        coefficients = []
        checked = 0
        for m in range(1, 7):
            row = []
            for j in range(m+1):
                exact = comb(m, j)*S.haar(S.mul(ppowers[j], qpowers[m-j]))
                expected = comb(m, j)*S.norm_moment(coeffs, n, j) if 2*j == m else Q(0)
                assert exact == expected
                row.append(str(exact))
                checked += 1
            coefficients.append(row)
        second = Q(coefficients[1][1])
        assert second == (Q(16, 15) if name == '2O' else Q(25, 546))
        cases[name] = {'n': n, 'second_moment_coefficient_of_ab': str(second),
                       'direct_parameter_moment_coefficients_checked': checked,
                       'coefficient_rows_a_power_0_to_m': coefficients,
                       'all_order_ordinary_ideal': '(a*b), by written proof',
                       'all_order_pure_support': 'a*b=0, all points marked-safe',
                       'q_axis_sharp_coefficient': str(Q(1, n)),
                       'sheared_axis_sharp_coefficient_t_nonzero': str(Q(n+1, n))}
    return cases


def main():
    receipt = {'completed': True, 'source_controls': source_controls(),
               'normalizer_entire_control': normalizer_control(),
               'exceptional_two_axis_families': exceptional_families(),
               'general_b_function_solver_implemented': False,
               'full_exceptional_coefficient_spaces_decided': False,
               'exceptional_Mathieu_Zhao_failure_constructed': False,
               'universal_proofs_verified_by_checker': False,
               'independent_replay': False, 'external_review': False,
               'formal_verification': False,
               'python_version': sys.version.split()[0],
               'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'imported_arithmetic_sha256': hashlib.sha256((ROOT/'verify_haar_shear_witnesses.py').read_bytes()).hexdigest()}
    (ROOT/'haar_period_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
