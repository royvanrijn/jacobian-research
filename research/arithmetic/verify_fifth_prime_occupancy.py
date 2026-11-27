"""Exact occupancy controls and bounded census; Python standard library only.

This checks reconstruction and rational exponents, not an asymptotic saving.
The census uses positive normalized coefficients <= V and ordered distinct
odd primes in [Z,2Z). Its field window is X=Z^10/V^2, with exact comparisons.
"""

import argparse
import json
from collections import Counter
from fractions import Fraction as F
from functools import cache
from itertools import permutations, product
from math import gcd, isqrt


def primes_between(z):
    return [p for p in range(max(3, z), 2*z)
            if all(p % r for r in range(2, isqrt(p)+1))]


def recover(primes, coefficients):
    p1, p2, p3 = primes
    x, y, z = coefficients
    a, b, c = (p**5 for p in primes)
    aa, bb, cc = a*y*y, b*z*z, c*x*x
    assert (aa-bb) % x == (cc-bb) % y == (aa-cc) % z == 0
    traces = [F(aa+bb-cc, y*z), F(bb+cc-aa, x*z),
              F(bb-aa-cc, x*y)]
    assert all(t.denominator == 1 for t in traces)
    u, v, w = map(int, traces)
    defect = F(4*a*b-u*u, x*x)
    assert defect.denominator == 1 and defect != 0
    e = int(defect)
    assert v*v+e*y*y == 4*b*c
    assert w*w+e*z*z == 4*a*c
    assert u*v-e*x*y == 2*b*w
    assert u*y+v*x == 2*b*z
    return (u, v, w), e


@cache
def squarefree_part(n):
    """Return d,delta with n=d*delta^2, for positive n."""
    d = delta = 1
    p = 2
    while p*p <= n:
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        d *= p**(exponent % 2)
        delta *= p**(exponent // 2)
        p = 3 if p == 2 else p+2
    return d*n, delta


def in_window(primes, coefficients, z_scale, v_cap):
    traces, e = recover(primes, coefficients)
    if max(map(abs, traces)) > 128*z_scale**5:
        return None
    if abs(e)*max(q*q for q in coefficients) > 128**2*z_scale**10:
        return None
    # Enumerate the short *complete* scale interval instead of factoring a
    # possibly huge defect. Only a candidate d in [X,2X) needs factoring.
    n = abs(e)
    low = isqrt(n*v_cap*v_cap // (2*z_scale**10))+1
    high = isqrt(n*v_cap*v_cap // z_scale**10)
    for delta in range(low, high+1):
        if n % (delta*delta) or gcd(delta, primes[0]*primes[1]*primes[2]) != 1:
            continue
        d = n // (delta*delta)
        if e < 0 and d == 1:
            continue  # The rational field.
        assert z_scale**10 <= d*v_cap*v_cap < 2*z_scale**10
        if squarefree_part(d)[1] == 1:
            return "imaginary" if e > 0 else "real"
    return None


def incident_units(primes, coefficients):
    p1, p2, p3 = primes
    x, y, z = coefficients
    return (gcd(x, p1*p2) == gcd(y, p2*p3) == gcd(z, p1*p3) == 1)


def census(z_scale, v_cap, keep_tuples=False):
    """Enumerate divisor records, CRT roots, then their actual small occupants."""
    primes = primes_between(z_scale)
    prime_set = set(primes)
    roots = [{} for _ in range(v_cap+1)]
    roots[1] = {0: [0]}
    for modulus in range(2, v_cap+1):
        for r in range(modulus):
            if gcd(r, modulus) == 1:
                roots[modulus].setdefault(pow(r, 5, modulus), []).append(r)
    counts = Counter()
    found = set()
    for p1, p2 in permutations(primes, 2):
        a, b = p1**5, p2**5
        for y, z in product(range(1, v_cap+1), repeat=2):
            if gcd(y, z) != 1 or gcd(y, p2) != 1 or gcd(z, p1) != 1:
                continue
            difference = a*y*y-b*z*z
            assert difference != 0
            q = y*z
            zy = 0 if y == 1 else z*pow(z, -1, y)
            yz = 0 if z == 1 else y*pow(y, -1, z)
            for x in range(1, v_cap+1):
                if difference % x or gcd(x, y*z*p1*p2) != 1:
                    continue
                counts["divisor_records"] += 1
                ry = (b*z*z*pow(x*x, -1, y)) % y if y > 1 else 0
                rz = (a*y*y*pow(x*x, -1, z)) % z if z > 1 else 0
                left, right = roots[y].get(ry, []), roots[z].get(rz, [])
                if left and right:
                    counts["records_with_compatible_classes"] += 1
                counts["compatible_classes"] += len(left)*len(right)
                for r1, r2 in product(left, right):
                    r = (r1*zy+r2*yz) % q
                    first = z_scale+(r-z_scale) % q
                    for p3 in range(first, 2*z_scale, q):
                        counts["integer_occupants"] += 1
                        if p3 not in prime_set or p3 in (p1, p2):
                            continue
                        pp, qq = (p1, p2, p3), (x, y, z)
                        assert incident_units(pp, qq)
                        counts["distinct_prime_occupants"] += 1
                        if keep_tuples:
                            found.add((pp, qq))
                        # Strict threshold min(q_i)>V/Z^(1/3), without rounding.
                        residual = min(qq)**3*z_scale > v_cap**3
                        if residual:
                            counts["residual_prime_occupants"] += 1
                        sign = in_window(pp, qq, z_scale, v_cap)
                        if sign:
                            counts["window_"+sign] += 1
                            if residual:
                                counts["residual_window_"+sign] += 1
    return dict(sorted(counts.items())), found


def controls():
    # Both signed fixtures are independent powers of explicit prime generators.
    for pp, qq, expected in (
        ((3, 11, 19), (-41, -7516, 281), 2*116**2),
        ((23, 71, 127), (-28594789, -479683349, -347203498), -2*28**2),
    ):
        assert recover(pp, qq)[1] == expected
        bad = (pp[0], pp[1], pp[2]+2)
        try:
            recover(bad, qq)
        except AssertionError:
            pass
        else:
            raise AssertionError("corrupted prime tuple was accepted")

    # A positive control for the field test itself (not a dyadic prime-box
    # fixture): fifth powers in Q(sqrt(-3)), from norms 7,13,19.
    field_pp, field_qq = (7, 13, 19), (11787, 527737, -117658)
    assert recover(field_pp, field_qq)[1] == 12
    assert in_window(field_pp, field_qq, 10, 65000) == "imaginary"
    assert in_window(field_pp, field_qq, 10, 100) is None

    # The residual coefficient sector is nonempty in the critical Z=16 box,
    # but this actual small-prime occupant fails its prescribed field window.
    pp, qq = (17, 23, 29), (57, 113, 53)
    assert min(qq)**3*16 > 128**3
    assert recover(pp, qq) == ((-5081129, 22042455, -10354167), 3304670587)
    assert squarefree_part(3304670587) == (3304670587, 1)
    assert in_window(pp, qq, 16, 128) is None

    # An independent direct six-variable enumeration checks CRT coverage.
    _, by_crt = census(10, 7, keep_tuples=True)
    direct = set()
    for pp in permutations(primes_between(10), 3):
        a, b, c = (p**5 for p in pp)
        for qq in product(range(1, 8), repeat=3):
            x, y, z = qq
            if gcd(x, y) != 1 or gcd(x*y, z) != 1 or not incident_units(pp, qq):
                continue
            if (a*y*y-b*z*z) % x == (c*x*x-b*z*z) % y == (a*y*y-c*x*x) % z == 0:
                direct.add((pp, qq))
    assert by_crt == direct and direct

    # Exact intersections and endpoint envelopes for the conditional fourth moment.
    theta = F(1, 3)
    b_new = (7-theta)/(20-2*theta)
    assert b_new == F(10, 29)
    assert (F(9)-theta)/2 == F(13, 3)
    assert 11-theta == F(32, 3)
    assert (2-5*b_new)/(F(1, 2)-b_new) == F(16, 9)
    low, high = F(11, 32), F(23, 66)
    assert low < b_new < F(9, 26) < high
    assert (2-5*low)/(F(1, 2)-low) == F(9, 5)
    assert (2-5*high)/(F(1, 2)-high) == F(17, 10)
    target = 1+3*b_new
    assert target == F(59, 29)
    assert F(53, 26)-target == F(3, 754)
    assert F(11, 2)-10*high == F(133, 66) < target
    assert 4-F(23, 4)*high == F(527, 264) < target
    assert 2 < target
    return len(direct)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Z", type=int, default=10)
    parser.add_argument("--V", type=int, default=56)
    args = parser.parse_args()
    if args.Z < 2 or args.V < 1:
        parser.error("Z>=2 and V>=1 are required")
    matched = controls()
    counts, _ = census(args.Z, args.V)
    print(json.dumps({
        "scope": "bounded census and algebra/exponent controls; no asymptotic saving",
        "Z": args.Z, "V": args.V, "X": str(F(args.Z**10, args.V**2)),
        "coefficient_signs": "positive; no multiplicative factor for other signs",
        "direct_control_tuples": matched, "counts": counts,
    }, indent=2))


if __name__ == "__main__":
    main()
