#!/usr/bin/env python3
"""Exact identities and bounded controls for ETF1--ETF5; no all-primes infinitude claim."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
import math
from pathlib import Path
import platform

from verify_elliptic_exception_explicit_local_criterion import valuation, residue, lower_hull

A, B = -10656, 27091730
AP, BP = 21312, 5183416


def on_curve(P, a=A, b=B):
    if P is None:
        return True
    x, y = P
    return y*y == x*x*x+a*x*x+b*x


def add(P, Q, a=A, b=B):
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P
    u, v = Q
    if x == u and y == -v:
        return None
    slope = (3*x*x+2*a*x+b)/(2*y) if P == Q else (v-y)/(u-x)
    xx = slope*slope-a-x-u
    return xx, slope*(x-xx)-y


def negative(P):
    return None if P is None else (P[0], -P[1])


def phi(P):
    if P is None or P[0] == 0:
        return None
    x, y = P
    return y*y/(x*x), y*(x*x-B)/(x*x)


def psi(P):
    if P is None or P[0] == 0:
        return None
    x, y = P
    return y*y/(4*x*x), y*(x*x-BP)/(8*x*x)


def cubic(x):
    z = F(3, 819200)*(128*x-25865)
    c = (x-5200)/1920
    return [z*c, -z, F(0), F(1)]


def normalized(poly, root_scale):
    return [v*root_scale**i/root_scale**3 for i, v in enumerate(poly)]


def ppow(p, n):
    return F(p**n) if n >= 0 else F(1, p**(-n))


def simple_roots(poly, p):
    values = [residue(v, p) for v in poly]
    return [w for w in range(p)
            if sum(v*w**i for i, v in enumerate(values)) % p == 0
            and sum(i*v*w**(i-1) for i, v in enumerate(values) if i) % p != 0]


def sf_product(a, b):
    return a*b//math.gcd(a, b)**2


def square_Q5(v):
    k = valuation(v, 5)
    return k is None or (k % 2 == 0 and pow(residue(v/ppow(5, k), 5), 2, 5) == 1)


def reduce3(P):
    if P is None or (P[0] != 0 and valuation(P[0], 3) < 0):
        return None
    return residue(P[0]-3552, 3), residue(P[1], 3)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-receipt', action='store_true')
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    checks = []

    def require(name, condition, **details):
        if not condition:
            raise AssertionError(name)
        checks.append({'name': name, 'passed': True, **details})

    P = (F(373248, 49), F(59459616, 343))
    T = (F(0), F(0))
    Q = (F(3584, 25), F(4301248, 125))
    S = (F(92179201, 6400), F(553221611649, 512000))
    R = add(P, T)
    q0 = 2709173
    require('rational_points_on_both_curves', all(on_curve(V) for V in [P, T, S, R])
            and on_curve(Q, AP, BP))
    require('explicit_dual_isogeny_image_with_sign', negative(psi(Q)) == S)
    require('dual_composition_is_doubling', psi(phi(P)) == add(P, P))
    require('torsion_translation_formula',
            R == (B/P[0], -B*P[1]/P[0]**2)
            and R == (F(663747385, 186624), F(-6525490184045, 80621568)))
    require('descent_witness_squareclasses',
            P[0]/2 == F(432, 7)**2 and B == 10*q0
            and Q[0]/14 == F(16, 5)**2 and F(BP, 26446) == F(14)**2
            and S[0] == F(9601, 80)**2)
    image1 = {1, 2, 5*q0, 10*q0}
    generated1 = {1, 2, 10*q0, sf_product(2, 10*q0)}
    image2 = {1, 14, 1889, 26446}
    generated2 = {1, 14, 26446, sf_product(14, 26446)}
    require('both_known_descent_images_generated', generated1 == image1 and generated2 == image2)
    require('mod_two_quotient_size_from_descent_kernel',
            len(image1)*len(image2)//2 == 8 == 2**(2+1),
            scope='Finite quotient arithmetic; the kernel and odd-index arguments are written in ETF1.')
    require('global_descent_classes_at_two', (5*q0) % 8 == 1
            and all((d == 1 or d % 8 == 1) for d in [1, 5*q0]))

    unit_classes, valuation_one_classes = set(), set()
    for x in range(1, 32):
        vx = valuation(x, 2)
        if vx not in (0, 1):
            continue
        if any((y*y-x*x*x-A*x*x-B*x) % 32 == 0 for y in range(32)):
            (unit_classes if vx == 0 else valuation_one_classes).add(
                x % 8 if vx == 0 else x//2 % 8)
    require('complete_mod32_forbidden_integral_unit_classes',
            unit_classes == {3} and valuation_one_classes == {3},
            scope='All x,y mod32 with v2(x)=0 or1; local descent classes are3 or6, absent from the global image.')

    require('formal_depth_four_witnesses',
            valuation(R[0], 2) == valuation(S[0], 2) == -8
            and valuation(R[1], 2) == valuation(S[1], 2) == -12)
    tR = -(R[0]-3552)/R[1]
    tS = -(S[0]-3552)/S[1]
    require('exact_formal_parameters',
            valuation(tR, 2) == valuation(tS, 2) == 4
            and residue(tR/16, 256) == 47 and residue(tS/16, 256) == 123)

    local_two_cases = 0
    all_two = True
    for n in range(4, 31):
        for unit in range(1, 16, 2):
            x = ppow(2, -2*n)*unit**2
            poly = cubic(x)
            pts = [(i, valuation(v, 2)) for i, v in enumerate(poly) if v]
            all_two &= lower_hull(pts) == [(0, -15-4*n), (3, 0)]
            all_two &= valuation(poly[1], 2) == -8-2*n
            rootval = F(-15-4*n, 3)
            all_two &= (rootval.denominator == 1) == (n % 3 == 0)
            if n % 3 == 0:
                norm = normalized(poly, ppow(2, int(rootval)))
                all_two &= all(v == 0 or valuation(v, 2) >= 0 for v in norm)
                all_two &= simple_roots(norm, 2) == [1]
            d = -F(3, 1600)*x
            all_two &= valuation(d, 2) % 2 == 0
            all_two &= residue(d/ppow(2, valuation(d, 2)), 8) == 5
            translated_x = B/x
            translated_poly = [v*ppow(2, 18-3*i) for i, v in enumerate(cubic(translated_x))]
            all_two &= valuation(translated_x, 2) == 1+2*n
            all_two &= all(v == 0 or valuation(v, 2) >= 0 for v in translated_poly)
            all_two &= simple_roots(translated_poly, 2) == [1]
            local_two_cases += 2
    require('bounded_two_adic_Newton_and_simple_root_controls', all_two, count=local_two_cases,
            scope='4<=formal depth<=30; odd units1..15; square formal x and its T0-translate.')
    P2 = (F(3442688, 529), F(304670016, 12167))
    twoP, twoP2 = add(P, P), add(P2, P2)
    require('rational_positive_and_negative_two_adic_controls',
            on_curve(P2) and valuation(twoP[0], 2) == -10
            and valuation(twoP2[0], 2) == -12,
            scope='2P has obstructed depth5;2P2 has soluble depth6. No universal obstruction for K=Q(sqrt(-3)).')

    five_count = 0
    all_five = True
    for x in range(1, 50):
        if x % 5 == 1:
            poly = [v*ppow(5, 3-i) for i, v in enumerate(cubic(F(x)))]
            all_five &= [residue(v, 5) for v in poly] == [2, 2, 0, 1]
            all_five &= 3 in simple_roots(poly, 5)
            five_count += 1
        elif x % 5 == 2:
            d = -F(3, 1600)*x
            all_five &= valuation(d, 5) % 2 == 0
            all_five &= pow(residue(d/ppow(5, valuation(d, 5)), 5), 2, 5) == 1
            five_count += 1
    for k in [1, 3, 5, 7]:
        for unit in range(1, 126):
            if unit % 5 == 0:
                continue
            x = F(5**k*unit)
            if not square_Q5(x*x*x+A*x*x+B*x):
                continue
            fv = valuation(128*x-25865, 5)
            poly = cubic(x)
            if k == 1:
                all_five &= fv in [1, 2]
                testpoly = [5*v for v in poly] if fv == 1 else poly
            else:
                all_five &= valuation(poly[0], 5) == 0 and valuation(poly[1], 5) == -1
                testpoly = [5*v for v in poly]
            all_five &= bool(simple_roots(testpoly, 5))
            five_count += 1
    for n in range(1, 13):
        for unit in range(1, 11):
            if unit % 5 == 0:
                continue
            x = ppow(5, -2*n)*unit**2
            poly = cubic(x)
            rootval = F(-4*n-3, 3)
            pts = [(i, valuation(v, 5)) for i, v in enumerate(poly) if v]
            all_five &= lower_hull(pts) == [(0, -4*n-3), (3, 0)]
            all_five &= (rootval.denominator == 1) == (n % 3 == 0)
            if n % 3 == 0:
                norm = normalized(poly, ppow(5, int(rootval)))
                all_five &= len(simple_roots(norm, 5)) == 1
            five_count += 1
    require('bounded_five_adic_controls', all_five, count=five_count,
            scope='Unit x1..49; positive valuations1,3,5,7 with units1..125 filtered only by exact local Q5-square condition for E; formal depths1..12 with units1..10.')

    x0 = F(25865, 128)
    value0 = x0*(F(639735, 128))**2
    deriv0 = 3*x0*x0+2*A*x0+B
    require('five_adic_Taylor_obstruction_constants',
            value0 == x0**3+A*x0*x0+B*x0
            and valuation(value0, 5) == 3 and valuation(deriv0, 5) >= 1)

    require('three_adic_reduction_coordinates',
            reduce3(T) == reduce3(P) == (0, 0) and reduce3(S) == (1, 0) and reduce3(R) is None)
    joint_count = 0
    for N in range(5, 13):
        modulus = 2**(N-4)
        r, s = residue(tR/16, modulus), residue(tS/16, modulus)
        generated = {(e, (k*r+d*s) % modulus, e, d)
                     for e in range(2) for d in range(2) for k in range(modulus)}
        expected = {(e, u, e, d) for e in range(2)
                    for d in range(2) for u in range(modulus)}
        require('joint_two_primary_model_N'+str(N), generated == expected
                and len(expected) == 4*modulus,
                scope='Exact finite cyclic model; formal addition is additive through these precisions by Lemma5.9 at depth4.')
        joint_count += len(expected)

    a2 = F(1, 7)
    require('two_adic_local_proportion', F(1, 2)*(1+a2) == F(4, 7))
    require('joint_two_three_local_proportion',
            F(1, 4)*(1+a2*(F(1, 13)+F(1, 3))) == F(289, 1092))
    retained = 0
    total = 0
    for comp in range(2):
        for three_B in range(2):
            for formal_two in range(4):
                for formal_three in range(3):
                    pass2 = comp == 1 or formal_two == 0
                    if comp:
                        pass3 = three_B == 1
                    else:
                        pass3 = formal_three == 0
                    retained += pass2 and pass3
                    total += 1
    require('fixed_finite_joint_quotient_bound',
            F(retained, total) == F(7, 24), retained=retained, total=total,
            scope='E6(2),E3(3): untested two-adic tail1/4; three-adic O/B tails1/3, A forbidden,A+B retained.')
    tail_controls = []
    for N in range(6, 16):
        modulus = 2**(N-4)
        soluble = tail = 0
        for u in range(modulus):
            if u == 0:
                tail += 1
            elif (4+valuation(u, 2)) % 3 == 0:
                soluble += 1
        low, high = F(soluble, modulus), F(soluble+tail, modulus)
        tail_controls.append(low <= a2 <= high)
    require('finite_formal_mass_brackets', all(tail_controls),
            scope='E4/E_N,6<=N<=15; finite lower/upper annulus masses bracket1/7.')

    sources = [
        'ELLIPTIC_EXCEPTION_TWO_FIVE_ADIC_CLASSIFICATION.md',
        'ELLIPTIC_EXCEPTION_ARITHMETIC.md',
        'ELLIPTIC_EXCEPTION_COUNT.md',
        'ELLIPTIC_EXCEPTION_EXPLICIT_LOCAL_CRITERION.md',
        'ELLIPTIC_EXCEPTION_COMPLETE_THREE_ADIC_CLASSIFICATION.md',
        'verify_elliptic_exception_explicit_local_criterion.py',
        'sources/dujella-number-theory-2021-section15-5.pdf',
        'sources/milne-elliptic-curves-2006.pdf',
        'sources/milne-ant300.pdf',
    ]
    report = {
        'status': 'passed',
        'command': 'python3 research/arithmetic/verify_elliptic_exception_two_five_adic_classification.py --write-receipt',
        'python': platform.python_version(),
        'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'sources_sha256': {p: sha256((here/p).read_bytes()).hexdigest() for p in sources},
        'check_count': len(checks), 'checks': checks,
        'scope': (
            'Exact point/isogeny/descent arithmetic, bounded two/five-adic Newton and Hensel controls, '
            'joint finite-group image models and exact local-count factors. Written descent kernels, '
            'local formal groups, the complete three-adic classification and height equidistribution '
            'establish the unrestricted statements. No basis assertion or elliptic Hasse infinitude.'),
        'formal_verification': False, 'independent_whole_proof_replay': False, 'external_review': False,
    }
    if args.write_receipt:
        (here/'elliptic_exception_two_five_adic_classification_checks.json').write_text(
            json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': 'passed', 'check_count': len(checks),
                      'two_adic_cases': local_two_cases, 'five_adic_cases': five_count,
                      'joint_model_states': joint_count, 'joint_proportion': '289/1092'}))


if __name__ == '__main__':
    main()
