#!/usr/bin/env python3
"""Exact bounded audit of the classical C3 primitive individual-norm criterion.
The unrestricted criterion requires the written ideal-class and signed-class
proof. This checker makes no additive sumset or finite-exception assertion.
"""


if not __debug__:
    raise RuntimeError(
        "This certificate checker requires active Python assertions; "
        "run without -O, -OO, or PYTHONOPTIMIZE."
    )
import argparse
from array import array
from hashlib import sha256
import json
import sys
from math import gcd, isqrt
from pathlib import Path

QR23 = {pow(a, 2, 23) for a in range(1, 23)}


def require_count(actual, predicted):
    if actual != predicted:
        raise AssertionError((actual, predicted))


def count_hash(counts):
    """Use fixed-width little-endian bytes for a portable receipt hash."""
    encoded = array('I', counts)
    assert encoded.itemsize == 4
    if sys.byteorder != 'little':
        encoded.byteswap()
    return sha256(encoded.tobytes()).hexdigest()


def coordinate_counts(bound, a, b, c):
    """Complete ellipse enumeration: 4aQ=(2ax+by)^2+23y^2."""
    assert b*b-4*a*c == -23
    counts = array('I', [0]) * (bound+1)
    for y in range(isqrt(4*a*bound//23)+1):
        root = isqrt(4*a*bound-23*y*y)
        lo = -((root+b*y)//(2*a))
        hi = (root-b*y)//(2*a)
        for x in range(lo, hi+1):
            n = a*x*x+b*x*y+c*y*y
            if n and gcd(x, y) == 1:
                assert n <= bound
                counts[n] += 1 if y == 0 else 2
    return counts


def sieve(bound):
    spf = array('I', range(bound+1))
    for p in range(2, isqrt(bound)+1):
        if spf[p] == p:
            for k in range(p*p, bound+1, p):
                if spf[k] == k:
                    spf[k] = p
    return spf, [p for p in range(2,bound+1) if spf[p] == p]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bound', type=int, default=10**6)
    parser.add_argument('--output', default=str(Path(__file__).with_name('primitive_norm_minus23_checks.json')))
    args = parser.parse_args()
    bound = args.bound
    if not 529 <= bound < 2**32:
        parser.error('--bound must be at least 529 and below 2**32')
    q0 = coordinate_counts(bound, 1, 1, 6)
    q1 = coordinate_counts(bound, 2, 1, 3)
    spf, primes = sieve(bound)
    kind = bytearray(bound+1) # 1 principal split, 2 nonprincipal split, 3 inert, 4 ramified
    prime_counts = dict(principal_split=0, nonprincipal_split=0, inert=0, ramified=0)
    for p in primes:
        if p == 23:
            kind[p] = 4
            prime_counts['ramified'] += 1
            assert q0[p] and not q1[p]
        elif p % 23 in QR23:
            assert bool(q0[p]) != bool(q1[p]), (p, q0[p], q1[p])
            kind[p] = 1 if q0[p] else 2
            prime_counts['principal_split' if q0[p] else 'nonprincipal_split'] += 1
        else:
            kind[p] = 3
            prime_counts['inert'] += 1
            assert not q0[p] and not q1[p]
    predicted = array('I', [0]) * (bound+1)
    local_count = 0
    class_failures = []
    r_histogram = {}
    for n in range(1,bound+1):
        rest, s, r, allowed = n, 0, 0, True
        while rest > 1:
            p = spf[rest]
            e = 0
            while rest % p == 0:
                rest //= p
                e += 1
            if kind[p] == 3 or (kind[p] == 4 and e >= 2):
                allowed = False
            elif kind[p] in (1,2):
                s += 1
                r += kind[p] == 2 and e % 3 != 0
        if allowed:
            local_count += 1
            r_histogram[str(r)] = r_histogram.get(str(r),0)+1
            orientation_count_numerator = (2**r+2*(-1)**r)
            assert orientation_count_numerator % 3 == 0
            predicted[n] = 2**(s-r+1) * (orientation_count_numerator//3)
            if r == 1:
                class_failures.append(n)
        require_count(q0[n], predicted[n])
    assert q0[2] == q0[3] == q0[4] == 0
    assert all(q0[n] for n in (1,6,8,12,18,23,36))
    assert 2 in class_failures and 8 not in class_failures
    # Signed-class criterion and exact roots-of-unity count, independent finite enumeration.
    for r in range(17):
        residues = [0]
        for _ in range(r):
            residues = [(z+eps)%3 for z in residues for eps in (-1,1)]
        count = residues.count(0)
        assert count == (2**r+2*(-1)**r)//3
        assert (count > 0) == (r != 1)
    # Independent square-box regression challenges the ellipse endpoints and signs.
    small = min(bound, 1000)
    square = array('I', [0])*(small+1)
    ycap = isqrt(4*small//23)+1
    xcap = isqrt(small)+ycap+1
    for x in range(-xcap,xcap+1):
        for y in range(-ycap,ycap+1):
            n=x*x+x*y+6*y*y
            if 1 <= n <= small and gcd(x,y)==1:
                square[n] += 1
    assert list(square) == list(q0[:small+1])
    # Challenge the same count validator with wrong membership predictions.
    assert kind[2] == kind[3] == 2
    negative_controls = []
    for n, wrong, label in (
        (2, 2, 'injected nonprincipal split prime 2'),
        (4, 2, 'injected imprimitive per-copy orientation at 4'),
        (8, 0, 'deleted allowed nonprincipal cubic power 8'),
        (5, 2, 'injected inert norm 5'),
        (529, 2, 'injected forbidden ramified square 23^2'),
    ):
        try:
            require_count(q0[n], wrong)
        except AssertionError:
            negative_controls.append(label)
        else:
            raise AssertionError(('negative control accepted', n, wrong))
    receipt = dict(
        checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        bound=bound,
        coordinate_principal_form=[1,1,6],
        coordinate_nonprincipal_form=[2,1,3],
        primitive_principal_norm_count=sum(bool(v) for v in q0),
        complete_local_prime_factor_candidate_count=local_count,
        class_obstruction_count=len(class_failures),
        first_class_obstructions=class_failures[:40],
        prime_class_counts=prime_counts,
        active_nonprincipal_prime_count_histogram=r_histogram,
        full_coordinate_vs_class_formula_equality=True,
        formula='2^(s-r+1)*(2^r+2*(-1)^r)/3 on inert-free n with v_23(n)<=1, else 0',
        signed_class_regression_through_r=16,
        independent_square_box_regression_bound=small,
        negative_controls=negative_controls,
        count_hash_encoding='unsigned 32-bit little endian',
        primitive_counts_sha256=count_hash(q0),
        predicted_counts_sha256=count_hash(predicted),
        scope='Exact coordinate-count audit only for 1<=n<=bound. The unrestricted individual-norm criterion depends on the written C3 ideal-class proof. No additive sumset or finite-exception classification is asserted.'
    )
    Path(args.output).write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__ == '__main__':
    main()
