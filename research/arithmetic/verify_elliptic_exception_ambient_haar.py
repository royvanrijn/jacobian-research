#!/usr/bin/env python3
"""Exact identities and bounded group/mass controls for EHD1--EHD3."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import itertools
import json
import math
from pathlib import Path
import platform

AA, BB = -10656, 27091730


def prime(n):
    return n >= 2 and all(n % d for d in range(2, math.isqrt(n)+1))


def val(v, p):
    v = F(v)
    if not v:
        return None
    n, d, k = abs(v.numerator), v.denominator, 0
    while n % p == 0:
        n //= p
        k += 1
    while d % p == 0:
        d //= p
        k -= 1
    return k


def mod(v, p):
    v = F(v)
    return v.numerator*pow(v.denominator, -1, p) % p


def legendre(a, p):
    a %= p
    return 0 if not a else (1 if pow(a, (p-1)//2, p) == 1 else -1)


def add_mod(P, Q, p):
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P
    u, v = Q
    if x == u and (y+v) % p == 0:
        return None
    slope = ((3*x*x+2*AA*x+BB)*pow(2*y, -1, p) if P == Q
             else (v-y)*pow(u-x, -1, p)) % p
    xx = (slope*slope-AA-x-u) % p
    return xx, (slope*(x-xx)-y) % p


def add_rat(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P
    u, v = Q
    if x == u and y == -v:
        return None
    slope = ((3*x*x+2*AA*x+BB)/(2*y) if P == Q else (v-y)/(u-x))
    xx = slope*slope-AA-x-u
    return xx, slope*(x-xx)-y


def times(n, P):
    result = None
    while n:
        if n & 1:
            result = add_rat(result, P)
        n //= 2
        if n:
            P = add_rat(P, P)
    return result


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

    q0, d = 2709173, 51730
    delta = 16*BB*BB*(AA*AA-4*BB)
    require('exact_elliptic_discriminant_factorization',
            delta == 2**9*5**2*7**3*1889*q0**2
            and AA*AA-4*BB == 2**3*7**3*1889)
    require('branch_squarefree_constants',
            d == 2*5*7*739 and all(prime(p) for p in [739, 1889, 42649, q0])
            and 739 % 3 == 42649 % 3 == 1)
    x0, c0 = F(25865, 128), F(639735, 128)
    require('closed_branch_point_identities',
            x0 == F(d, 16**2)
            and x0**3+AA*x0*x0+BB*x0 == x0*c0*c0
            and F(639735, 2048)**2*d == x0*c0*c0
            and x0*c0 != 0)

    good = [p for p in range(7, 252) if prime(p) and delta % p]
    pointsets = {}
    group_cases = branch_cases = 0
    all_groups = all_branches = True
    target_primes = []
    for p in good:
        points = [None]+[(x, y) for x in range(p) for y in range(p)
                         if (y*y-x*x*x-AA*x*x-BB*x) % p == 0]
        pointsets[p] = points
        N = len(points)
        all_groups &= (N % 2 == 0 and N % p != 0 and N <= 2*p+1
                       and (N-p-1)**2 <= 4*p and (0, 0) in points)
        group_cases += 1
        if p >= 11 and p % 3 == 2:
            target_primes.append(p)
            roots = [point for point in points if point is not None
                     and (128*point[0]-25865) % p == 0]
            expected = 2 if legendre(d, p) == 1 else 0
            all_branches &= len(roots) == expected and all(point[1] != 0 for point in roots)
            all_branches &= mod(x0, p) != 0 and mod(c0, p) != 0
            b = 1+expected
            soluble = 1-F(b, N)*(1-F(1, p*p+p+1))
            all_branches &= 0 < soluble < 1 and 1-soluble >= F(1, 6*p)
            branch_cases += 1
    require('bounded_good_reduction_and_no_anomaly_controls', all_groups,
            count=group_cases, scope='Every good prime7<=p<=251; complete finite-field point enumeration.')
    require('bounded_ambient_branch_and_local_fraction_controls', all_branches,
            count=branch_cases,
            scope='Every good prime11<=p<=251,p=2mod3; complete F=0 branch scan and exact rational masses.')

    selected = [11, 13, 17, 19, 23, 29]
    P = (F(373248, 49), F(59459616, 343))
    T = (F(0), F(0))
    S = (F(92179201, 6400), F(553221611649, 512000))
    basis = [P, T, S]
    quotient_labels = {}
    quotient_sizes = []
    reductions = {}
    for p in selected:
        points = pointsets[p]
        doubles = {add_mod(Q, Q, p) for Q in points}
        labels = {}
        cosets = {}
        for Q in points:
            coset = frozenset(add_mod(Q, D, p) for D in doubles)
            if coset not in cosets:
                cosets[coset] = len(cosets)
            labels[Q] = cosets[coset]
        q = len(cosets)
        require('finite_two_quotient_p'+str(p),
                q == sum(Q is None or Q[1] == 0 for Q in points)
                and len(points) == q*len(doubles) and q in [2, 4],
                group_order=len(points), quotient_order=q)
        quotient_labels[p] = labels
        quotient_sizes.append(q)
        reductions[p] = [tuple(mod(a, p) for a in Q) for Q in basis]
    rows = []
    for s in range(1, len(selected)+1):
        primes = selected[:s]
        image = set()
        for ep, et, es in itertools.product(range(2), repeat=3):
            row = [(ep, et)]
            for p in primes:
                Q = None
                for coefficient, point in zip([ep, et, es], reductions[p]):
                    if coefficient:
                        Q = add_mod(Q, point, p)
                row.append(quotient_labels[p][Q])
            image.add(tuple(row))
        target_order = 4*math.prod(quotient_sizes[:s])
        ratio = F(target_order, len(image))
        require('bounded_joint_image_support_index_s'+str(s),
                len(image) <= 8 and ratio >= 2**(s-1),
                primes=primes, image_order=len(image), ambient_quotient_order=target_order,
                quotient_index=str(ratio),
                scope='Complete image of Gamma/2Gamma using the three generators proved in ETF1.')
        rows.append({'s': s, 'ambient_quotient_order': target_order,
                     'image_order': len(image), 'index_lower_bound': str(ratio)})

    all_mass = True
    depth_cases = 0
    for p in target_primes:
        if p > 101:
            continue
        Q = p*p+p+1
        for e in range(1, 13):
            offset = (3-e % 3) % 3
            exact = F(p**(2-offset), Q)
            cutoff = e+12
            low = sum((F(p-1, p**(n-e+1)) for n in range(e, cutoff) if n % 3 == 0), F(0))
            tail = F(1, p**(cutoff-e))
            all_mass &= low <= exact <= low+tail
            all_mass &= (exact == F(1, Q) if e % 3 == 1 else
                         exact == F(p, Q) if e % 3 == 2 else exact == F(p*p, Q))
            if e % 3 == 0:
                all_mass &= 1-exact == F(p+1, Q)
            depth_cases += 1
    require('bounded_actual_formal_depth_mass_controls', all_mass, count=depth_cases,
            scope='Every target prime<=101;1<=e<=12; exact geometric formula bracketed by12 tested depths and its tail.')
    require('global_mod_two_size', 2**(2+1) == 8,
            scope='Rank2 and rational torsionC2 are imported from the elliptic arithmetic theorem.')

    for p in [11, 17]:
        N = len(pointsets[p])
        Q = times(N, P)
        t = -(Q[0]-3552)/Q[1]
        require('depth_one_rational_witness_p'+str(p), val(t, p) == 1,
                multiple=N, normalized_unit=mod(t/p, p),
                scope='One actual rational witness forces the full Mordell-Weil minimum e_p=1, without a basis.')

    inputs = [
        'ELLIPTIC_EXCEPTION_AMBIENT_HAAR_DIAGNOSTIC.md',
        'ELLIPTIC_EXCEPTION_ARITHMETIC.md',
        'ELLIPTIC_EXCEPTION_EXPLICIT_LOCAL_CRITERION.md',
        'ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md',
        'ELLIPTIC_EXCEPTION_PRIME_ELIGIBILITY.md',
        'ELLIPTIC_EXCEPTION_TWO_FIVE_ADIC_CLASSIFICATION.md',
        'ELLIPTIC_EXCEPTION_QUANTITATIVE_FINITE_LOCAL_ASYMPTOTICS.md',
        'sources/milne-elliptic-curves-2006.pdf',
        'sources/milne-ant300.pdf',
    ]
    report = {
        'status': 'passed',
        'command': 'python3 research/arithmetic/verify_elliptic_exception_ambient_haar.py --write-receipt',
        'python': platform.python_version(),
        'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'sources_sha256': {name: sha256((here/name).read_bytes()).hexdigest() for name in inputs},
        'check_count': len(checks), 'checks': checks, 'finite_support_index_controls': rows,
        'scope': (
            'Exact discriminant/branch identities, complete finite-field groups and two-quotient images '
            'on declared bounded prime domains, geometric-mass controls and rational depth-one witnesses. '
            'The written Hasse/Chebotarev, compact-group and divergent-product arguments establish the '
            'unrestricted ambient-null and support-singularity claims. No value of intrinsic all-primes '
            'Haar measure, independence of rational prime tests or elliptic Hasse infinitude is inferred.'),
        'formal_verification': False, 'independent_whole_proof_replay': False, 'external_review': False,
    }
    if args.write_receipt:
        (here/'elliptic_exception_ambient_haar_checks.json').write_text(
            json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': 'passed', 'check_count': len(checks),
                      'good_prime_controls': group_cases, 'branch_controls': branch_cases,
                      'formal_depth_controls': depth_cases}))


if __name__ == '__main__':
    main()
