#!/usr/bin/env python3
"""Exact finite support for CN2EffectiveUnit5; not a spectral proof."""

if not __debug__:
    raise RuntimeError("This certificate requires Python assertions enabled.")

import argparse
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from math import isqrt, prod
from pathlib import Path


def primes_through(bound):
    flags = bytearray(b"\x01") * (bound + 1)
    flags[:2] = b"\x00\x00"
    for p in range(2, isqrt(bound) + 1):
        if flags[p]:
            flags[p*p:bound+1:p] = b"\x00" * ((bound-p*p)//p+1)
    return [p for p in range(2, bound+1) if flags[p]]


def hmul(x, y):
    a,b,c,d = x
    e,f,g,h = y
    return (a*e-b*f-c*g-d*h, a*f+b*e+c*h-d*g,
            a*g-b*h+c*e+d*f, a*h+b*g-c*f+d*e)


def omul5(x, y):
    """O/5O in the basis 1,iota,j,iota*j, with iota^2=-5."""
    a,b,c,d = x
    e,f,g,h = y
    return ((a*e-c*g)%5, (a*f+b*e-d*g+c*h)%5,
            (a*g+c*e)%5, (a*h+b*g+d*e-c*f)%5)


def norm5(x):
    return (x[0]*x[0]+x[2]*x[2])%5


def legendre(a, p):
    r = pow(a%p, (p-1)//2, p)
    return -1 if r == p-1 else r


def sigma(n):
    return sum(d+n//d if d*d != n else d
               for d in range(1, isqrt(n)+1) if n%d == 0)


def tau(n):
    return sum(2 if d*d != n else 1
               for d in range(1, isqrt(n)+1) if n%d == 0)


def lattice_counts(bound):
    r = [0]*(bound+1)
    even_pair = [0]*(bound+1)
    odd_pair = [0]*(bound+1)
    for y in range(-isqrt(bound//5), isqrt(bound//5)+1):
        xmax = isqrt(bound-5*y*y)
        for x in range(-xmax, xmax+1):
            n = x*x+5*y*y
            r[n] += 1
            if x%2 and y%2:
                odd_pair[n] += 1
            if x%2 == 0 and y%2 == 0:
                even_pair[n] += 1
    return r, even_pair, odd_pair


def exact_constants():
    small = primes_through(63)
    euler_E = prod(max(F((a+1)**6 * p**(5*a),
                         sum(p**j for j in range(a+1))**6)
                      for a in range(7))
                   for p in small if p not in (2,5))
    euler_M = prod(max(F((a+1)**6, p**a) for a in range(7))
                   for p in small if p%20 in (1,3,7,9))
    assert euler_E < F(15,2)**6
    assert euler_M < 7**6
    assert F(8,7)**6 < 3
    assert 4*7**6 < 9**6
    B = sum((F(2,p*(p-1)) for p in small
             if p%20 in (11,13,17,19)), F(0))
    tail_starts = [71,73,77,79]
    assert {a%20 for a in tail_starts} == {11,13,17,19}
    B_upper = B + F(5,2)*sum(
        (F(1,a*a)+F(1,20*a) for a in tail_starts), F(0))
    assert B_upper < F(1,16)
    # The progression bounds used for every real P, with generous endpoint loss.
    for P in range(1,1000):
        count = sum(n%20 in (11,13,17,19) for n in range(1,P+1))
        assert count <= F(P,5)+2
    assert 40**12 < 2**64
    sufficient_upper = F(21,80)+F(261,400)+F(513,2**21)
    assert sufficient_upper == F(47985177,52428800)
    assert sufficient_upper < F(15,16)
    # d_e/2^(2e/3)>3/4, certified without irrational arithmetic.
    degrees = [1,2,2]
    assert all(64*d**3 > 27*2**(2*e) for e,d in enumerate(degrees))
    assert F(3,4*3888) > F(1,6000)
    return {"E_sixth_product": str(euler_E),
            "M_sixth_product": str(euler_M),
            "inert_mass_upper": str(B_upper),
            "infinite_sieve_upper_at_m_2_64": str(sufficient_upper),
            "margin_below_15_over_16": str(F(15,16)-sufficient_upper),
            "dyadic_degrees": degrees}


def check_local_five():
    elements = list(product(range(5), repeat=4))
    units = [q for q in elements if norm5(q)]
    axis = [q for q in units if q[:2] == (0,0) or q[2:] == (0,0)]
    assert len(units) == 400 and len(axis) == 40
    identity = (1,0,0,0)
    for q in units:
        inv_n = pow(norm5(q), -1, 5)
        conj = (q[0], -q[1], -q[2], -q[3])
        inv = tuple(t*inv_n%5 for t in conj)
        assert omul5(q, inv) == omul5(inv,q) == identity
    axis_set = set(axis)
    for q in axis:
        assert legendre(norm5(q),5) == 1
        for r in axis:
            assert omul5(q,r) in axis_set
    # Conditional norm fibers independently authenticate the quadratic bias.
    densities = {}
    for n in range(1,5):
        whole = [q for q in units if norm5(q) == n]
        bad = [q for q in axis if norm5(q) == n]
        density = F(len(bad), len(whole))
        assert density == F(1+legendre(n,5), 10)
        densities[str(n)] = str(density)
    # Separate Borel reduction computes the local maximal-order index six.
    gl2_units = sum((a*d-b*c)%5 != 0
                    for a,b,c,d in product(range(5), repeat=4))
    borel_units = sum((a*d)%5 != 0
                      for a,b,d in product(range(5), repeat=3))
    assert (gl2_units, borel_units) == (480,80)
    total_mass = F(1,12)*3*F(gl2_units,borel_units)
    principal_probability = F(1,2)/total_mass
    assert total_mass == F(3,2) and principal_probability == F(1,3)
    return {"ring_units": len(units), "axis_units": len(axis),
            "projective_quotient_order": len(units)//4,
            "projective_axis_order": len(axis)//4,
            "conditional_axis_densities": densities,
            "maximal_mod5_units": gl2_units,
            "borel_mod5_units": borel_units,
            "projective_total_mass": str(total_mass),
            "principal_probability": str(principal_probability)}


def check_dyadic():
    units = [q for q in product(range(2),repeat=4) if sum(q)%2]
    patterns = [q for q in product(range(2),repeat=4) if sum(q) == 2]
    orbits = []
    while patterns:
        q = patterns[0]
        orbit = {tuple(t%2 for t in hmul(q,u)) for u in units}
        assert len(orbit) == 2
        assert all(sum(v) == 2 for v in orbit)
        orbits.append(orbit)
        patterns = [p for p in patterns if p not in orbit]
    assert len(orbits) == 3
    good = {q for orbit in orbits for q in orbit
            if (q[0]+q[1])%2 and (q[2]+q[3])%2}
    assert len(good) == 4
    assert sum(orbit <= good for orbit in orbits) == 2
    for q,i in [((1,1,0,0),(0,1,0,0)),
                ((1,0,1,0),(0,0,1,0)),
                ((1,0,0,1),(0,0,0,1))]:
        assert hmul(q,q) == tuple(2*t for t in i)
    return {"norm_two_parity_orbits":
            [sorted(list(v)) for v in orbits],
            "cross_pair_orbits_selected": 2}


def check_shells(bound=512):
    r, even_pair, odd_pair = lattice_counts(bound)
    cumulative = 0
    for n in range(1,25):
        cumulative += r[n]
        assert cumulative <= 2*n
    # eta(2z)^2 eta(10z)^2 coefficients; this does not prove the modular identity.
    c = [0]*(bound+1)
    c[1] = 1
    for step in list(range(2,bound+1,2))+list(range(10,bound+1,10)):
        previous = c[:]
        for n in range(bound+1):
            c[n] = previous[n]
            if n >= step:
                c[n] -= 2*previous[n-step]
            if n >= 2*step:
                c[n] += previous[n-2*step]
    full = [sum(r[a]*r[n-a] for a in range(n+1))
            for n in range(bound+1)]
    for n in range(1,bound+1):
        sig = lambda d: sigma(n//d) if n%d == 0 else 0
        assert 3*full[n] == 4*sig(1)-16*sig(4)+20*sig(5)-80*sig(20)+8*c[n]
        m,e = n,0
        while m%2 == 0:
            m//=2
            e+=1
        if n%5 == 0 or n%8 not in (2,3,4,6,7):
            continue
        if e == 0:
            selected = full[n]
            # Neither a nonprimitive pair nor a zero component occurs.
            assert all(not (r[a]*r[n-a]) or
                       (a and n-a and not even_pair[a] and not even_pair[n-a])
                       for a in range(n+1))
        elif e == 1:
            selected = sum(r[a]*r[n-a] for a in range(1,n,2))
        else:
            selected = sum(odd_pair[a]*odd_pair[n-a] for a in range(n+1))
            assert selected == full[n]-full[n//4]
        d = [1,2,2][e]
        gap = 4*d*sigma(m)-3*selected
        if gap > 0:
            assert gap*gap <= 64*d*d*tau(m)**2*m
    # Negative control for the tempting but false full degree-three e=1 choice.
    assert full[6] == 16
    assert sum(r[a]*r[6-a] for a in range(1,6,2)) == 8
    return {"full_shell_formula_comparison_through": bound,
            "selected_shell_spectral_normalization_regression_through": bound,
            "lattice_initial_bound": 24,
            "negative_full_shell_at_six": 16,
            "correct_odd_odd_selection_at_six": 8}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("primitive_norm_minus20_effective_checks.json"))
    args = parser.parse_args()
    receipt = {"checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "local_five": check_local_five(),
               "dyadic": check_dyadic(),
               "shells": check_shells(),
               "constants": exact_constants(),
               "scope": "Finite local, normalization, Euler and rational certificates only. The written Hecke transfer, ideal-class mass argument and infinite prime-tail proof are not established by execution.",
               "claimed_theorem_scope": "N>=2^66, 5 does not divide N, N mod8 in {2,3,4,6,7}; ordered primitive numerical-value count >=N^(2/3)/6000. No unrestricted exception-list completion."}
    args.output.write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
