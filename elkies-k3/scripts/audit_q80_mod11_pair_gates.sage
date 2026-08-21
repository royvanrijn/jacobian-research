#!/usr/bin/env sage
"""Audit affine pair-equality multiplicities for the two full p=11 q=80 hits.

Finite-field square-root signs are arbitrary, so all relative signs of P2 and
P3 are reported after fixing P1.  GCD degrees are necessary affine pieces of
the transported pair intersections (2,3,1); local contributions at singular
fibers or infinity still require a separate intersection calculation.
"""

from itertools import product as itertools_product
from sage.all import *


field = GF(11)
polynomials.<T> = PolynomialRing(field)


def polynomial(coefficients):
    return sum(field(value)*T**index for index, value in enumerate(coefficients))


def pair_intersection(P, Q):
    """Return (P-Q).O from the reduced chord denominator."""
    XP, YP, ZP = P
    XQ, YQ, ZQ = Q
    D = XP*ZQ**2-XQ*ZP**2
    S = YP*ZQ**3+YQ*ZP**3
    H = D*ZP*ZQ
    N = S**2-D**2*(XP*ZQ**2+XQ*ZP**2)
    first = gcd(H, N)
    cancellation = gcd(H, N//first)
    H_reduced = H//cancellation
    N_reduced = N//cancellation**2
    finite = H_reduced.degree()
    excess = N_reduced.degree()-2*finite-4
    assert excess <= 0 or excess % 2 == 0
    at_infinity = max(0, excess//2)
    return finite+at_infinity


seeds = {
    "standard": {
        "P1": ((0, 1, 2, 0), (0, 0, 5, 3, 3, 0)),
        "P2": ((9, 6, 4, 7, 1), (6, 6, 0, 3, 1, 5, 1)),
        "P3": (5, (4, 7, 5, 10, 8, 0, 3), (3, 1, 8, 10, 5, 3, 5, 5, 0, 4)),
    },
    "second": {
        "P1": ((0, 9, 5, 0), (0, 0, 4, 6, 1, 0)),
        "P2": ((1, 2, 10, 2, 5), (10, 8, 7, 10, 7, 10, 2)),
        "P3": (9, (4, 7, 6, 10, 3, 7, 1), (3, 1, 2, 4, 6, 1, 5, 5, 5, 1)),
    },
}


for label, data in seeds.items():
    x1, y1 = map(polynomial, data["P1"])
    x2, y2 = map(polynomial, data["P2"])
    pole, x3, y3 = data["P3"]
    Z = T-field(pole)
    x3, y3 = polynomial(x3), polynomial(y3)
    for sign2, sign3 in itertools_product((1, -1), repeat=2):
        gcds = (
            gcd(x1-x2, y1-sign2*y2),
            gcd(x1*Z**2-x3, y1*Z**3-sign3*y3),
            gcd(x2*Z**2-x3, sign2*y2*Z**3-sign3*y3),
        )
        gates = tuple(value.degree() for value in gcds)
        intersections = (
            pair_intersection((x1, y1, polynomials.one()), (x2, sign2*y2, polynomials.one())),
            pair_intersection((x1, y1, polynomials.one()), (x3, sign3*y3, Z)),
            pair_intersection((x2, sign2*y2, polynomials.one()), (x3, sign3*y3, Z)),
        )
        print(
            f"Q80PAIR11|seed={label}|signs=1,{sign2},{sign3}|"
            f"affine_gcd_degrees={gates}|resolved_intersections={intersections}|gcds={gcds}",
            flush=True,
        )
