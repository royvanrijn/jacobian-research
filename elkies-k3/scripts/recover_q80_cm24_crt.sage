#!/usr/bin/env sage
"""Reproduce the retracted standard-P1 q=80 CRT experiment.

This is an exact CRT/LLL reconstruction over Q(sqrt(-6)).  Candidate
projective triples (a,b,c) represent (a+b*sqrt(-6))/c.  The residues at
7, 11, and 29 are now known to have the wrong P2 component at the residual
I2, so this is a negative historical diagnostic, not target recognition.
The partial seeds at 31 fail the P3 gate and are excluded.  A split-prime 97 filter
removes false joint reconstructions before exact discriminants are computed.
"""

from sage.all import *
from itertools import product as cartesian_product


K.<s> = QuadraticField(-6)
KT.<T> = PolynomialRing(K)
KU.<z> = PolynomialRing(K)
primes = (7, 11, 29)
crt_modulus = prod(primes)
candidate_cache = {}

filter_field = GF(97)
filter_s = filter_field(-6).sqrt(all=True)[0]
FT.<t> = PolynomialRing(filter_field)
FU.<u> = PolynomialRing(filter_field)
filter_jet_matrix_inverse = matrix(
    filter_field, 4, 4,
    lambda row, column: FU((1+u) ** (4+column))[row],
).inverse()
reduction_cache = {}


def candidates(residues, roots, count):
    key = (tuple(residues), tuple(roots), count)
    if key in candidate_cache:
        return candidate_cache[key]
    root_crt = ZZ(CRT_list(list(roots), list(primes)))
    residue_crt = ZZ(CRT_list(list(residues), list(primes)))
    lattice = matrix(ZZ, [
        [crt_modulus, 0, 0],
        [-root_crt, 1, 0],
        [residue_crt, 0, 1],
    ]).LLL()
    values = {}
    for coefficients in cartesian_product(range(-10, 11), repeat=3):
        if coefficients == (0, 0, 0):
            continue
        a, b, c = map(ZZ, vector(ZZ, coefficients) * lattice)
        if c == 0:
            continue
        common = gcd((a, b, c))
        a, b, c = a // common, b // common, c // common
        if any(c % prime == 0 for prime in primes):
            continue
        if c < 0:
            a, b, c = -a, -b, -c
        value = K(a+b*s) / c
        score = a*a + 6*b*b + c*c
        values[value] = min(score, values.get(value, infinity))
    answer = tuple(
        value for value, _ in sorted(values.items(), key=lambda item: item[1])[:count]
    )
    assert len(answer) == count
    candidate_cache[key] = answer
    return answer


def reduce97(value):
    if value not in reduction_cache:
        rational, radical = K(value).list()
        reduction_cache[value] = (
            filter_field(rational.numerator()) / filter_field(rational.denominator())
            + filter_s * filter_field(radical.numerator()) / filter_field(radical.denominator())
        )
    return reduction_cache[value]


def ambient(base_t, base_u, d, p, q, e, inverse_jet_matrix=None):
    base = base_t.base_ring()
    r = -3*d**2 + 3-p-q
    A = base_t**2 * (-3+p*base_t+q*base_t**2+r*base_t**3)
    A_one = base_u(A(base_t=1+base_u))
    v = (A_one+3*d**2) / (-3*d**2)
    branch = 2*d**3 * (
        1 + base(3)/2*v + base(3)/8*v**2 - base(1)/16*v**3
    )
    jet_matrix = inverse_jet_matrix
    if jet_matrix is None:
        jet_matrix = matrix(
            base, 4, 4,
            lambda row, column: base_u((1+base_u)**(4+column))[row],
        ).inverse()
    fixed = vector(base, [
        base_u(2*(1+base_u)**3+e*(1+base_u)**8)[row]
        for row in range(4)
    ])
    b = jet_matrix * (vector(base, [branch[row] for row in range(4)])-fixed)
    B = base_t**3 * (
        2 + sum(b[row]*base_t**(row+1) for row in range(4)) + e*base_t**5
    )
    return A, B


def cm_discriminant(A, B, variable):
    discriminant = 4*A**3 + 27*B**2
    fixed = variable**7 * (variable-1)**4
    residual, remainder = discriminant.quo_rem(fixed)
    return not remainder and residual.discriminant() == 0


hits = []
tested = 0
filter_passed = 0
for root29 in (9, 20):
    roots = (1, 4, root29)
    d_values = candidates((3, 3, 3), roots, 10)
    p_values = candidates((4, 2, 5), roots, 25)
    q_values = candidates((3, 3, 8), roots, 25)
    e_values = candidates((2, 9, 5), roots, 25)
    for d, p, q, e in cartesian_product(
        d_values, p_values, q_values, e_values
    ):
        tested += 1
        try:
            A97, B97 = ambient(
                t, u,
                reduce97(d), reduce97(p), reduce97(q), reduce97(e),
                filter_jet_matrix_inverse,
            )
        except ZeroDivisionError:
            continue
        if not cm_discriminant(A97, B97, t):
            continue
        filter_passed += 1
        A, B = ambient(T, z, d, p, q, e)
        if cm_discriminant(A, B, T):
            hits.append((roots, d, p, q, e, A, B))
            print(
                f"Q80CRT|HIT|roots={roots}|"
                f"d={d}|p={p}|q={q}|e={e}",
                flush=True,
            )

print(
    f"Q80CRT|SUMMARY|tested={tested}|mod97={filter_passed}|hits={len(hits)}|"
    "primes=7,11,29|status=RETRACTED_PROFILE_BRANCH_DIAGNOSTIC",
    flush=True,
)
