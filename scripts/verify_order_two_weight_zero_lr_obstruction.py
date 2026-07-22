#!/usr/bin/env python3
"""Exact weight-zero order-two audit for the filtered LR obstruction."""

import sympy as sp


u, gamma, w = sp.symbols("u gamma w")
v = u - 1
S = gamma - 1 + sp.Rational(8, 7) * v
W = u * gamma


def H(root):
    return sp.expand(root**2 * (root - 1) * (3 * root**2 - 11 * root + 6) / 60)


# The invariant presentation F_2=(x^-2 a,x^-1 b,x gamma) and its
# logarithmic differential matrix in the source basis
# (delta x/x,delta u,delta gamma).
H_2 = H(W)
p_2 = sp.diff(H(w), w).subs(w, W)
c_2 = sp.Rational(1, 30)
q_2 = sp.cancel((W * p_2 - H_2) / c_2)
a = sp.expand(sp.cancel(u + q_2 / gamma**2))
b = sp.expand(sp.cancel(c_2 + p_2 / gamma))
J = sp.Matrix(
    [
        [-2 * a, sp.diff(a, u), sp.diff(a, gamma)],
        [-b, sp.diff(b, u), sp.diff(b, gamma)],
        [gamma, 0, 1],
    ]
)
assert sp.factor(J.det()) == c_2

# The saturated target quotient is
# R/(a,b^2) + R/(b,a gamma) + R/(gamma).
G_1 = sp.groebner([a, b**2], gamma, u, order="grevlex")
G_2 = sp.groebner([b, a * gamma], gamma, u, order="grevlex")


def quotient_remainders(f, g, h):
    """Class of (x*f,y*g,z*h) in the three saturated quotient summands."""
    source = sp.Matrix(
        [
            f,
            (f + g) * v,
            -sp.Rational(8, 7) * (f + g) * v + (2 * f + h) * S,
        ]
    )
    target = J * source
    return (
        G_1.reduce(sp.expand(target[0]))[1],
        G_2.reduce(sp.expand(target[1]))[1],
        sp.expand(target[2].subs(gamma, 0)),
    )


# Ordinary-degree 25, torus-weight-zero diagonal fields have invariant
# coefficients of weighted degree 24 for deg(v)=2 and deg(S)=3.
monomials = [v ** (12 - 3 * index) * S ** (2 * index) for index in range(5)]

# Compute the exact kernel of the quotient map on the 15-dimensional space
# (x*f,y*g,z*h), with f,g,h in the span of the five monomials.
columns = []
row_keys = set()
for component in range(3):
    for monomial in monomials:
        f = monomial if component == 0 else 0
        g = monomial if component == 1 else 0
        h = monomial if component == 2 else 0
        data = {}
        for quotient_index, remainder in enumerate(quotient_remainders(f, g, h)):
            for exponent, coefficient in sp.Poly(
                sp.expand(remainder), gamma, u
            ).terms():
                key = (quotient_index,) + exponent
                data[key] = coefficient
                row_keys.add(key)
        columns.append(data)

row_keys = sorted(row_keys)
quotient_matrix = sp.Matrix(
    [[column.get(key, 0) for column in columns] for key in row_keys]
)
assert quotient_matrix.rank() == 10

# The five-dimensional kernel is exactly P(v,S)*(x,-y,-2z).
expected_kernel = []
for index in range(5):
    vector = sp.zeros(15, 1)
    vector[index] = 1
    vector[5 + index] = -1
    vector[10 + index] = -2
    assert quotient_matrix * vector == sp.zeros(quotient_matrix.rows, 1)
    expected_kernel.append(vector)
assert sp.Matrix.hstack(*expected_kernel).rank() == 5

# Normalize the nonzero torus residue to N*(x,0,-3z), where N=v^6*S^4.
# Every representative in this leading sector is therefore
#
#   Z=N*(x,0,-3z)+P*(x,-y,-2z).
#
# The restriction of Z to the invariant plane is
# N*(v d/dv-S d/dS), independently of P.  If E(P) denotes this Euler
# derivative, the diagonal coefficients of (DZ)Z/2 are A,B,C below.
parameters = sp.symbols("p0:5")
N = v**6 * S**4
P = sum(parameter * monomial for parameter, monomial in zip(parameters, monomials))
E_P = sum(
    parameters[index] * (12 - 5 * index) * monomials[index]
    for index in range(5)
)
A = sp.expand(((N + P) ** 2 + 2 * N**2 + N * E_P) / 2)
B = sp.expand((P**2 - N * E_P) / 2)
C = sp.expand((3 * N**2 + 12 * N * P + 4 * P**2 - 2 * N * E_P) / 2)

# The last quotient summand alone detects every quadratic reconstruction.
# Its residue is the gamma=0 specialization of delta(gamma).
normal_residue = sp.expand(
    (
        -sp.Rational(8, 7) * (A + B) * v
        + (2 * A + C) * S
    ).subs(gamma, 0)
)
residue_coefficients = sp.Poly(normal_residue, u).coeffs()
parameter_ideal = sp.groebner(residue_coefficients, *parameters, order="grevlex")
assert len(parameter_ideal.polys) == 1
assert parameter_ideal.polys[0].as_expr() == 1

print("PASS: the degree-25 weight-zero saturated kernel is P(v,S)*(x,-y,-2z)")
print("PASS: no weight-zero kernel choice kills the order-two quadratic residue")
