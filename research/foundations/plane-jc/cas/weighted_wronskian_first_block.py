#!/usr/bin/env python3
"""Structural checks for the weighted-Wronskian first block.

This is deliberately much smaller than the archived Groebner replay.  It
checks the operator reduction, the six-dimensional cokernel, the binary
transvectant identity, the mu_7 grading of the archived lex basis, and the
descent of its degree-35 eliminant to a quotient quintic.
"""

from pathlib import Path
import re

import sympy as sp


t, x, y = sp.symbols("t x y")

# A=tU and D=t^2 V turn 2AD'-3A'D=t^2 into T_U(V)=1.
u = sp.symbols("u0:8")
v = sp.symbols("v0:11")
U = sum(u[i] * t**i for i in range(8))
V = sum(v[i] * t**i for i in range(11))
A = t * U
D = t**2 * V
weighted_wronskian = sp.expand(2 * A * sp.diff(D, t) - 3 * sp.diff(A, t) * D)
reduced_operator = sp.expand(U * V + 2 * t * U * sp.diff(V, t) - 3 * t * sp.diff(U, t) * V)
assert sp.expand(weighted_wronskian - t**2 * reduced_operator) == 0

# If V=V0+US, the correction operator maps P_3 to P_9.  The apparent
# degree-ten coefficient is 1+2*3-7=0, and a generic specialization has rank 4.
c = sp.symbols("c0:4")
S = sum(c[i] * t**i for i in range(4))
correction = sp.expand(U * S + 2 * t * U * sp.diff(S, t) - t * sp.diff(U, t) * S)
assert sp.Poly(correction, t).degree() <= 9
matrix = sp.Matrix(
    [[sp.Poly(correction, t).nth(i).coeff(cj) for cj in c] for i in range(10)]
)
assert matrix.subs({u[i]: i + 1 for i in range(8)}).rank() == 4

# For a monic degree-eight A, the affine basis t^i dt/y (0<=i<=6) has
# one logarithmic direction.  Its residue at infinity gives the displayed
# linear relation and leaves the six compact-curve de Rham directions.
s = sp.symbols("s")
r = sp.symbols("r0:7")
a5_inf, a6_inf, a7_inf = sp.symbols("a5_inf a6_inf a7_inf")
inverse_sqrt = sp.series(
    (1 + a7_inf * s + a6_inf * s**2 + a5_inf * s**3) ** sp.Rational(-1, 2),
    s,
    0,
    4,
).removeO()
residue = sp.expand(
    sum(r[i] * s ** (2 - i) for i in range(7)) * inverse_sqrt
).coeff(s, -1)
expected_residue = (
    r[3]
    - a7_inf * r[4] / 2
    + (3 * a7_inf**2 / 8 - a6_inf / 2) * r[5]
    + (-a5_inf / 2 + 3 * a7_inf * a6_inf / 4 - 5 * a7_inf**3 / 16) * r[6]
)
assert sp.expand(residue - expected_residue) == 0

# Homogenization identifies the operator with a binary first transvectant.
a = sp.symbols("a0:9")
d = sp.symbols("d0:13")
binary_A = sum(a[i] * x**i * y ** (8 - i) for i in range(9))
binary_D = sum(d[i] * x**i * y ** (12 - i) for i in range(13))
binary_jacobian = sp.diff(binary_A, y) * sp.diff(binary_D, x) - sp.diff(
    binary_A, x
) * sp.diff(binary_D, y)
affine_A = sum(a[i] * t**i for i in range(9))
affine_D = sum(d[i] * t**i for i in range(13))
assert sp.expand(
    binary_jacobian.subs({x: t, y: 1})
    - (8 * affine_A * sp.diff(affine_D, t) - 12 * sp.diff(affine_A, t) * affine_D)
) == 0

# Read the preserved lex basis rather than recomputing its Groebner basis.
root = Path(__file__).resolve().parents[1]
lex_output = (
    root
    / "external/zenodo-21479814/"
    "bilLkarkariy-jc2-72-108-exact-certificates-d9ea4fd/"
    "release_bundle/exact_replay/firstblock_Q_exact.out"
)
lines = lex_output.read_text().splitlines()
lex_lines = [line for line in lines if re.match(r"L\[[1-6]\]=", line)]
assert len(lex_lines) == 6

a2, a3, a4, a5, a6, a7, q = sp.symbols("a2 a3 a4 a5 a6 a7 q")
avars = (a2, a3, a4, a5, a6, a7)
local_dict = {str(z): z for z in avars}
lex = [
    sp.Poly(sp.sympify(line.split("=", 1)[1].replace("^", "**"), locals=local_dict), *avars)
    for line in lex_lines
]

# The residual scaling is a_k -> zeta^(k-1) a_k.  Every lex relation is
# homogeneous for this Z/7 grading.
weights = {avars[i]: i + 1 for i in range(6)}
for relation in lex:
    term_weights = {
        sum(weights[var] * exponent for var, exponent in zip(avars, monomial)) % 7
        for monomial, _ in relation.terms()
    }
    assert len(term_weights) == 1

# L[1] is h(a7^7), not an intrinsically degree-35 equation.
H = lex[0].as_expr()
H_terms = sp.Poly(H, a7).terms()
assert {exponent[0] % 7 for exponent, _ in H_terms} == {0}
h = sum(coefficient * q ** (exponent[0] // 7) for exponent, coefficient in H_terms)
assert sp.expand(h.subs(q, a7**7) - H) == 0
assert sp.degree(h, q) == 5
assert len(sp.factor_list(h, q)[1]) == 1
galois_group, contained_in_alternating = sp.polys.numberfields.galois_group(
    sp.Poly(h, q, domain=sp.QQ)
)
assert galois_group.order() == 120
assert not contained_in_alternating

# Multiplying the a_k relation by a7^(k-1) expresses the invariant
# x_k=a_k*a7^(k-1) as a polynomial in q, for k=2,...,6.
quotient_degrees = []
for k, relation in zip(range(6, 1, -1), lex[1:]):
    ak = avars[k - 2]
    expression = relation.as_expr()
    coefficient = sp.diff(expression, ak)
    assert not coefficient.has(*avars)
    rest = sp.expand(expression - coefficient * ak)
    invariant_numerator = sp.Poly(-rest * a7 ** (k - 1), a7)
    assert {exponent[0] % 7 for exponent, _ in invariant_numerator.terms()} == {0}
    invariant_polynomial = sum(
        value * q ** (exponent[0] // 7)
        for exponent, value in invariant_numerator.terms()
    ) / coefficient
    quotient_degrees.append(sp.degree(invariant_polynomial, q))

print("PASS: A=tU, D=t^2V gives the weighted Euler operator T_U(V)=1")
print("PASS: correction P3 -> P9 has generic rank 4 and cokernel dimension 6")
print("PASS: the infinity-residue kernel gives six compact de Rham directions")
print("PASS: the binary Jacobian is 4(2AD'-3A'D)")
print("PASS: all six archived lex relations respect the residual mu_7 grading")
print("PASS: H(a7)=h(a7^7) with quotient degree", sp.degree(h, q))
print("PASS: the quotient quintic is irreducible with Galois group S5")
print("PASS: invariant-coordinate relation degrees", quotient_degrees)
