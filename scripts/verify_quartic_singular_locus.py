#!/usr/bin/env python3
"""Radical certificate for the quartic nonproperness hypersurface."""

import sympy as sp


A, B, C, U = sp.symbols("A B C U")
Q4 = (
    A-B**2
    + C**2*(27*B**4-36*A*B**2+8*A**2)
    + 16*A**3*C**4
)
N = C*Q4
jacobian_ideal = [N, sp.diff(N, A), sp.diff(N, B), sp.diff(N, C)]
jacobian_gb = sp.groebner(jacobian_ideal, A, B, C, order="lex")

# Three geometric strata: the component intersection, the lifted node, and the
# two conjugate lifted cusps.
plane_intersection = [C, A-B**2]
node_lift = [B, 4*A*C**2+1]
cusp_lifts = [12*A*C**2-1, 27*B**2*C**2-2]


def intersection(left, right):
    gb = sp.groebner(
        [U*polynomial for polynomial in left]
        + [(1-U)*polynomial for polynomial in right],
        U, A, B, C, order="lex",
    )
    return [
        sp.factor(polynomial.as_expr())
        for polynomial in gb.polys
        if not polynomial.as_expr().has(U)
    ]


radical_generators = intersection(
    intersection(plane_intersection, node_lift), cusp_lifts
)
assert radical_generators == [
    16*A**2*C**2+4*A-27*B**4*C**2-4*B**2,
    B*(16*A-27*B**4*C**2-16*B**2),
    C*(4*A*C**2-18*B**2*C**2+1),
    B*C*(27*B**2*C**2-2),
]
radical_gb = sp.groebner(radical_generators, A, B, C, order="lex")

# J is contained in the proposed radical, and the square of every radical
# generator lies in J. Hence sqrt(J) equals the displayed intersection.
assert all(radical_gb.reduce(polynomial)[1] == 0 for polynomial in jacobian_ideal)
assert all(
    jacobian_gb.reduce(sp.expand(generator**2))[1] == 0
    for generator in radical_generators
)

# Direct substitution checks each component. The cusp polynomial is squarefree
# in BC and splits into the two complex cusp lifts.
assert all(sp.factor(polynomial.subs({C: 0, A: B**2})) == 0
           for polynomial in jacobian_ideal)
assert all(sp.factor(polynomial.subs({B: 0, A: -1/(4*C**2)})) == 0
           for polynomial in jacobian_ideal)
for sign in (1, -1):
    cusp_substitution = {
        A: 1/(12*C**2),
        B: sign*sp.sqrt(6)/(9*C),
    }
    assert all(sp.factor(polynomial.subs(cusp_substitution)) == 0
               for polynomial in jacobian_ideal)
X = sp.symbols("X")
assert sp.gcd(sp.Poly(27*X**2-2, X), sp.Poly(54*X, X)).degree() == 0

print("PASS: radical(Jac(C*Q4)) is certified by square membership")
print("PASS: singular locus contains exactly the plane intersection, node lift, and two cusp lifts")
print("PASS: the cusp pair is reduced over C")
