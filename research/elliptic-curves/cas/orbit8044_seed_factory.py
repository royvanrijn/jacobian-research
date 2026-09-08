#!/usr/bin/env python3
"""Exact, Sage-free contracts for the height-ordered orbit8044 seed factory."""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import gcd, isqrt

import det1092_funnel as f
from alternate_quartic_covers import short_add
from v3_warm_support import atomic, point_tuple, read, require, sha


def parameters(max_height, maximum):
    """Prefix in (max(|a|,b),b,a) order; b>0, gcd(a,b)=1, including0/1."""
    require(type(max_height) is int and 1 <= max_height <= 10000, 'height cap outside bounds')
    require(type(maximum) is int and 1 <= maximum <= 100000, 'address cap outside bounds')
    count = 0
    for h in range(1, max_height+1):
        for b in range(1, h+1):
            for a in range(-h, h+1) if b == h else (-h, h):
                if gcd(a, b) == 1:
                    yield dict(id=f'orbit8044-{count:06d}', index=count,
                               u=str(F(a, b)), a=a, b=b, height=h)
                    count += 1
                    if count == maximum:
                        return


def image(chart, a, b):
    n = f.homogeneous(list(map(int, chart['s_numerator'])), a, b)
    d = f.homogeneous(list(map(int, chart['s_denominator'])), a, b)
    require(n or d, 'parametrization has a base point')
    return F(n, d) if d else None


def rational_j(model):
    A, B = map(F, model[3:])
    require(4*A**3+27*B**2 != 0, 'singular equation')
    return F(1728*4*A**3, 4*A**3+27*B**2)


def rational_power(value, degree):
    value = F(value)
    if value <= 0:
        return False
    def root(n):
        lo, hi = 0, 1 << ((n.bit_length()+degree-1)//degree)
        while lo+1 < hi:
            mid = (lo+hi)//2
            if mid**degree <= n:
                lo = mid
            else:
                hi = mid
        return hi if hi**degree == n else lo
    return all(root(n)**degree == n for n in (value.numerator, value.denominator))


def rational_isomorphic(left, right):
    """Exact Q-isomorphism for short nonsingular equations, including special j."""
    require(len(left) == len(right) == 5 and not any(map(F, (*left[:3], *right[:3]))), 'short equations required')
    if rational_j(left) != rational_j(right):
        return False
    A, B = map(F, left[3:]); C, D = map(F, right[3:])
    if not A:
        return not C and rational_power(B/D, 6)
    if not B:
        return not D and rational_power(A/C, 4)
    # If A=u^4 C and B=u^6 D, then u^2=BC/(AD).
    u2 = B*C/(A*D)
    return rational_power(u2, 2) and A == u2**2*C and B == u2**3*D


def word_point(model, base, word):
    require(len(word) == len(base) and all(type(n) is int and abs(n) <= 8 for n in word), 'bounded exact word required')
    result = None
    for n, P in zip(word, base):
        Q = P if n >= 0 else (P[0], -P[1])
        for _ in range(abs(n)):
            result = short_add(model, result, Q)
    return result


def obvious_dependencies(model, base, points, trace_word):
    """Exact signed singles/pairs, plus the exact two-branch trace identity."""
    n = len(base); found = {}; targets = set(points)
    def consider(P, word):
        if P in targets:
            found.setdefault(P, word)
        if P is not None and (P[0], -P[1]) in targets:
            found.setdefault((P[0], -P[1]), [-v for v in word])
    consider(None, [0]*n)
    for i, P in enumerate(base):
        word = [0]*n; word[i] = 1; consider(P, word)
    for i, j in combinations_with_replacement(range(n), 2):
        if all(P in found for P in points):
            break
        word = [0]*n; word[i] += 1; word[j] += 1
        consider(short_add(model, base[i], base[j]), word)
        if i != j:
            word = [0]*n; word[i] = 1; word[j] = -1
            consider(short_add(model, base[i], (base[j][0], -base[j][1])), word)
    trace_checked = False
    if len(points) == 2:
        trace_checked = short_add(model, *points) == word_point(model, base, trace_word)
        require(trace_checked, 'conic trace identity differs')
        for i in range(2):
            if points[i] in found:
                word = [a-b for a, b in zip(trace_word, found[points[i]])]
                require(word_point(model, base, word) == points[1-i], 'inferred generic word failed exact addition')
                found.setdefault(points[1-i], word)
    return dict(trace_identity_checked=trace_checked,
        points=[dict(index=i, status='EXACT_INHERITED_WORD' if P in found else 'UNKNOWN',
                     word=found.get(P)) for i, P in enumerate(points)])


def queue_admission(folder, case, registry):
    """Idempotent queue publication before its admission record; no seed is dropped."""
    seed = folder/'seeds'/case
    packet = read(seed/'m18.json'); checked = read(seed/'standalone-replay.json')
    require(packet['status'] == 'CERTIFIED_M18' and checked['status'] == 'PASS_STANDALONE_FUNNEL_M18', 'unverified seed')
    require(checked['rank_lower_bound'] == checked['matrix_rank'] == 18, 'rank gate differs')
    require(checked['inputs'][str((seed/'m18.json').relative_to(f.ROOT))] == sha(seed/'m18.json'), 'seed changed after proof')
    require(packet['intake_sha256'] == sha(seed/'intake.json'), 'intake binding differs')
    model = packet['curve']; j = str(rational_j(model))
    matches = [r for r in registry if r['j'] == j and rational_isomorphic(model, r['curve'])]
    result = dict(case=case, parameter=packet['parameter'], curve=model, j=j,
        seed_sha256=sha(seed/'m18.json'), proof_sha256=sha(seed/'standalone-replay.json'),
        protocol_sha256=sha(folder/'protocol.json'),
        status='DUPLICATE_RATIONAL_ISOMORPHISM' if matches else 'QUEUED_M18',
        duplicate_of=matches[0]['case'] if matches else None)
    if not matches:
        atomic(folder/'v3-queue'/f'{case}.json', dict(status='READY_FOR_UNCHANGED_V3', **{k:v for k,v in result.items() if k != 'status'}), immutable=True)
    atomic(folder/'admissions'/f'{case}.json', result, immutable=True)
    return result
