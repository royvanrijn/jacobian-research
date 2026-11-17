#!/usr/bin/env python3
"""Exact identities and bounded controls for ELV1--ELV4; not an infinitude test."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
import math
from pathlib import Path
import platform


def trim(p):
    p = list(map(F, p))
    while len(p) > 1 and not p[-1]:
        p.pop()
    return p


def add(a, b):
    out = [F(0)]*max(len(a), len(b))
    for i, value in enumerate(a):
        out[i] += value
    for i, value in enumerate(b):
        out[i] += value
    return trim(out)


def scale(a, value):
    return trim([F(value)*v for v in a])


def mul(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out)


def power(a, exponent):
    out = [F(1)]
    for _ in range(exponent):
        out = mul(out, a)
    return out


def compose(a, b):
    out = [F(0)]
    for c in reversed(a):
        out = add(mul(out, b), [c])
    return out


def derivative(a):
    return trim([i*a[i] for i in range(1, len(a))] or [0])


def remainder(a, b):
    a, b = trim(a), trim(b)
    if b == [0]:
        raise ZeroDivisionError
    while a != [0] and len(a) >= len(b):
        shift = len(a)-len(b)
        factor = a[-1]/b[-1]
        for i, value in enumerate(b):
            a[shift+i] -= factor*value
        a = trim(a)
    return a


def gcd_poly(a, b):
    a, b = trim(a), trim(b)
    while b != [0]:
        a, b = b, remainder(a, b)
    return scale(a, 1/a[-1])


def evaluate(a, value):
    out = F(0)
    for c in reversed(a):
        out = out*value+c
    return out


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, math.isqrt(n)+1, 2))


def valuation(value, p):
    value = F(value)
    if not value:
        return None
    result = 0
    numerator, denominator = abs(value.numerator), value.denominator
    while numerator % p == 0:
        numerator //= p
        result += 1
    while denominator % p == 0:
        denominator //= p
        result -= 1
    return result


def residue(value, p):
    value = F(value)
    return value.numerator % p * pow(value.denominator % p, -1, p) % p


def ppow(p, exponent):
    return F(p**exponent) if exponent >= 0 else F(1, p**(-exponent))


def lower_hull(points):
    hull = []
    for point in sorted(points):
        while len(hull) >= 2:
            a, b = hull[-2:]
            cross = (b[0]-a[0])*(point[1]-b[1])-(b[1]-a[1])*(point[0]-b[0])
            if cross > 0:
                break
            hull.pop()
        hull.append(point)
    return hull


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-receipt', action='store_true')
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    checks = []

    def require(name, value, **details):
        if not value:
            raise AssertionError(name)
        checks.append({'name': name, 'passed': True, **details})

    x = [F(0), F(1)]
    v = scale(add(x, [432]), F(3, 6400))
    vL, lam, kap = F(243483, 819200), F(10, 9), F(42649, 16384)
    z = add(v, [-vL])
    c = add(scale(z, lam), [-kap])
    Fpoly = [-25865, 128]
    require('coordinate_and_rational_function_identity', z == scale(Fpoly, F(3, 819200)))
    require('centered_constant_coefficient', c == scale(add(x, [-5200]), F(1, 1920)))

    s0, r, shift = F(1816347, 819200), F(-12, 5), F(4, 5)
    sv = add([s0], scale(v, -1))
    tv = add(add(scale(power(v, 2), lam), scale(v, -lam*s0)), [F(125, 1024)])
    shifted_constant = add(add([shift**3+r*shift**2], scale(sv, shift)), tv)
    shifted_linear = add([3*shift**2+2*r*shift], sv)
    require('exact_cubic_translation',
            shifted_constant == mul(z, c) and shifted_linear == scale(z, -1)
            and 3*shift+r == 0)
    disc = add(scale(power(z, 3), 4), scale(mul(power(z, 2), power(c, 2)), -27))
    Qv = compose([143467857, -110592000, 20480000], v)
    old_disc = scale(mul(power(scale(z, 819200), 2), Qv), F(-1, 412316860416000000))
    require('centered_cubic_discriminant_identity', disc == old_disc)
    require('quadratic_discriminant_coordinate_identity',
            add([F(81, 100)], scale(v, -4)) == scale(x, F(-3, 1600)))
    require('exceptional_prime_is_prime_and_split_mod_three',
            is_prime(42649) and 42649 % 3 == 1)
    require('unit_constants_outside_exceptional_set',
            lam == F(2*5, 3**2) and kap == F(42649, 2**14)
            and F(3, 819200) == F(3, 2**15*5**2))
    require('prime_1889_is_inert_and_not_exceptional',
            is_prime(1889) and 1889 % 3 == 2 and 1889 not in {2, 3, 5, 42649})

    primes = [p for p in range(7, 102) if is_prime(p)]
    newton_count = normalization_count = squareclass_count = 0
    all_newton = all_normalization = all_squareclass = True
    for p in primes:
        for unit in range(1, min(p, 9)):
            for k in range(-6, 7):
                if k == 0:
                    continue
                zv = ppow(p, k)*unit
                cv = lam*zv-kap
                coeff = [zv*cv, -zv, F(0), F(1)]
                points = [(i, valuation(value, p)) for i, value in enumerate(coeff) if value]
                expected = [(0, k if k > 0 else 2*k), (3, 0)]
                all_newton &= lower_hull(points) == expected
                all_newton &= valuation(cv, p) == (0 if k > 0 else k)
                root_value = F(k if k > 0 else 2*k, 3)
                all_newton &= (root_value.denominator == 1) == (k % 3 == 0)
                newton_count += 1
                if k % 3 == 0:
                    wscale = ppow(p, k//3 if k > 0 else 2*k//3)
                    norm = [value*wscale**i/wscale**3 for i, value in enumerate(coeff)]
                    reductions = [residue(value, p) for value in norm]
                    all_normalization &= all(valuation(value, p) is None
                                             or valuation(value, p) >= 0 for value in norm)
                    all_normalization &= (reductions[1:3] == [0, 0]
                                          and reductions[-1] == 1 and reductions[0] != 0)
                    all_normalization &= (-27*reductions[0]**2) % p != 0
                    if p % 3 == 2:
                        roots = [u for u in range(p) if (u**3+reductions[0]) % p == 0]
                        all_normalization &= (len(roots) == 1 and roots[0] != 0
                                              and (3*roots[0]**2) % p != 0)
                    normalization_count += 1
                dv = 4*zv**3-27*(zv*cv)**2
                expected_val = 2*k if k > 0 else 4*k
                all_squareclass &= valuation(dv, p) == expected_val
                dunit = residue(dv/ppow(p, expected_val), p)
                all_squareclass &= pow(dunit, (p-1)//2, p) == pow((-3) % p, (p-1)//2, p)
                squareclass_count += 1
    require('bounded_Newton_polygons', all_newton, count=newton_count,
            scope='All primes 7<=p<=101; 1<=unit<=min(p-1,8); -6<=k<=6, k!=0.')
    require('bounded_normalized_unramified_reductions', all_normalization,
            count=normalization_count,
            scope='Same domain with 3|k; separability for every p, unique simple root when p=2 mod3.')
    require('bounded_discriminant_minus_three_squareclasses', all_squareclass,
            count=squareclass_count, scope='Same nonzero-k domain; exact valuations and unit Legendre symbols.')

    unit_count = 0
    all_unit_cases = True
    prime_1889_passed = True
    for p in primes+[1889]:
        l, hkap = residue(lam, p), residue(kap, p)
        for u in range(1, p):
            constant = u*(l*u-hkap) % p
            dmod = (4*u**3-27*constant**2) % p
            split_resolvent = dmod != 0 and pow(dmod, (p-1)//2, p) == 1
            simple_root = any((w**3-u*w+constant) % p == 0
                              and (3*w*w-u) % p != 0 for w in range(p))
            passed = split_resolvent or simple_root
            all_unit_cases &= passed
            if p == 1889:
                prime_1889_passed &= passed
            unit_count += 1
    require('bounded_integral_unit_case', all_unit_cases, count=unit_count,
            scope='Every nonzero z residue for every prime 7<=p<=101, and for p=1889.')
    require('complete_unit_residue_control_at_1889', prime_1889_passed, count=1888)

    old_path = here/'elliptic_exception_hasse_certificate.json'
    field_path = here/'elliptic_exception_field_discriminant_certificate.json'
    old = json.loads(old_path.read_text(encoding='utf-8'))
    field = json.loads(field_path.read_text(encoding='utf-8'))
    seed_v = F(old['v'])
    seed_x = F(6400, 3)*seed_v-432
    seed_z = seed_v-vL
    seed_F = 128*seed_x-25865
    require('exact_seed_x_coordinate', seed_x == F(173387072000, 9601**2))
    require('exact_seed_rational_function_factorization',
            seed_F == F(5*13*1321*230703199, 9601**2)
            and seed_z == F(3, 819200)*seed_F)
    require('seed_factors_are_prime',
            all(is_prime(p) for p in [5, 13, 1321, 9601, 230703199]))
    predicted_support = {p for p in [13, 1321, 9601, 230703199] if valuation(seed_F, p) % 3}
    actual_support = {int(p) for p, exponent in field['conductor_factorization'] if int(exponent)}
    require('seed_conductor_support_matches_maximal_order_certificate',
            predicted_support == actual_support
            and all(p % 3 == 1 for p in predicted_support),
            support=sorted(predicted_support))

    integer_count = 0
    all_integer_cases = True
    bad_primes = [p for p in range(11, 60) if is_prime(p) and p % 3 == 2]
    for a in range(-40, 41):
        for denom in range(1, 17):
            if math.gcd(a, denom) != 1:
                continue
            N = 128*a-25865*denom**2
            if N == 0:
                continue
            all_integer_cases &= 128 % math.gcd(N, denom) == 0
            for p in bad_primes:
                vn, vc = valuation(N, p), valuation(denom, p)
                vf = valuation(F(N, denom**2), p)
                all_integer_cases &= vf == vn-2*vc
                all_integer_cases &= not (vn > 0 and vc > 0)
                all_integer_cases &= ((vf % 3 == 0) == (vn % 3 == 0 and vc % 3 == 0))
                integer_count += 1
    require('bounded_integer_numerator_denominator_equivalence', all_integer_cases,
            count=integer_count,
            scope='-40<=a<=40, 1<=c<=16, gcd(a,c)=1; all primes 11<=p<=59, p=2 mod3.')

    E = [F(0), F(27091730), F(-10656), F(1)]
    x0 = F(25865, 128)
    require('elliptic_cubic_is_squarefree', gcd_poly(E, derivative(E)) == [1])
    require('cubeclass_branch_value_nonzero',
            evaluate(E, x0) == F(10585532408369625, 2097152)
            and evaluate(E, x0) == x0*F(639735, 128)**2)
    cubeclass_controls = [F(1), F(-1), F(19), F(5, 7)]
    all_genus = True
    for rclass in cubeclass_controls:
        curve = compose(E, [x0, F(0), F(0), rclass/128])
        all_genus &= len(curve)-1 == 9 and gcd_poly(curve, derivative(curve)) == [1]
    require('bounded_hyperelliptic_polynomial_controls', all_genus,
            representatives=[str(v) for v in cubeclass_controls],
            scope='Exact degree9 and squarefreeness for four representatives; written preimage argument proves all nonzero r.')
    require('genus_four_arithmetic', (9-1)//2 == 4 and 2*4+2 == 10)

    files = [
        'ELLIPTIC_EXCEPTION_EXPLICIT_LOCAL_CRITERION.md',
        'ELLIPTIC_EXCEPTION_ARITHMETIC.md',
        'ELLIPTIC_EXCEPTION_FIELD_DISCRIMINANT.md',
        'RANK_FIVE_DISCRIMINANT_CLASSIFICATION.md',
        'elliptic_exception_hasse_certificate.json',
        'elliptic_exception_field_discriminant_certificate.json',
        'sources/milne-ant300.pdf',
        'sources/faltings-1983-endlichkeitssaetze.pdf',
    ]
    report = {
        'status': 'passed',
        'command': 'python3 research/arithmetic/verify_elliptic_exception_explicit_local_criterion.py --write-receipt',
        'python': platform.python_version(),
        'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'sources_sha256': {p: sha256((here/p).read_bytes()).hexdigest() for p in files},
        'check_count': len(checks), 'checks': checks,
        'scope': (
            'Exact rational and polynomial identities, exceptional-prime arithmetic, certified-seed factorization, '
            'bounded Newton/residue/valuation controls, and genus-polynomial identities. The written proofs '
            'establish the all-prime reduction and fixed-cubeclass finiteness using local-field theorems, RDF3 '
            'and Faltings. No all-places elliptic infinitude follows from these controls.'),
        'formal_verification': False,
        'independent_whole_proof_replay': False,
        'external_review': False,
    }
    if args.write_receipt:
        (here/'elliptic_exception_explicit_local_criterion_checks.json').write_text(
            json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': 'passed', 'check_count': len(checks),
                      'newton_controls': newton_count, 'unit_residue_controls': unit_count,
                      'integer_controls': integer_count}))


if __name__ == '__main__':
    main()
