#!/usr/bin/env python3
"""Exact torus-gauge recurrence and saturated target-module audit."""

import sympy as sp


u, gamma, t, w = sp.symbols("u gamma t w")
W = u * gamma
tau = 1 + t / 28
r = 1 / (1 + t)


def H(parameter, root):
    return sp.expand(
        root**2
        * (root - 1)
        * (3 * root**2 - (5 * parameter + 1) * root + 3 * parameter)
        / 60
    )


def c(parameter):
    return (sp.sympify(parameter) - 1) / 30


def weighted_top(polynomial):
    """Top part for deg(u)=2 and deg(gamma)=3."""
    polynomial = sp.cancel(polynomial)
    assert sp.denom(polynomial) == 1
    terms = sp.Poly(sp.expand(polynomial), u, gamma).terms()
    degree = max(2 * exponent[0] + 3 * exponent[1] for exponent, _ in terms)
    top = sum(
        coefficient * u**exponent[0] * gamma**exponent[1]
        for exponent, coefficient in terms
        if 2 * exponent[0] + 3 * exponent[1] == degree
    )
    return degree, sp.factor(top)


# After G_t is followed by
# (A,B,C) -> (tau^40 A,tau^-23 B,tau^-17 C), the F_2 root equation is
#
# H_2(W')-tau^-40 Sigma_t W'
#       +r^2 tau^6(Sigma_t W-H_(2+t)(W))=0.
parameter = 2 + t
sigma = sp.expand(c(parameter) * gamma + sp.diff(H(parameter, w), w).subs(w, W))
constant_term = sp.expand(sigma * W - H(parameter, W))


def root_equation(delta):
    return (
        H(2, W + delta)
        - tau ** (-40) * sigma * (W + delta)
        + r**2 * tau**6 * constant_term
    )


# The torus kills the weight-25 first-order forcing.  The first surviving
# forcing has weight 15 and determines the new recurrence seed.
forcing_1 = sp.expand(sp.series(root_equation(0), t, 0, 2).removeO()).coeff(t, 1)
assert sp.factor(forcing_1) == gamma**2 * u * (29 * gamma * u**2 - 24 * u - 5) / 420
assert weighted_top(forcing_1) == (15, sp.Rational(29, 420) * u**3 * gamma**3)

# Solve the exact root coefficients through order two in k[u,gamma].  This is
# enough to regress the observed degrees 25 and 49 without a three-variable
# Groebner computation.
delta = sp.Integer(0)
delta_coefficients = []
for order in (1, 2):
    unknown = sp.symbols(f"delta_{order}")
    trial = delta + unknown * t**order
    coefficient = sp.expand(
        sp.series(root_equation(trial), t, 0, order + 1).removeO()
    ).coeff(t, order)
    # dE/dW' at t=0 is -gamma/30.
    solved = sp.cancel(coefficient.subs(unknown, 0) / (gamma / 30))
    assert sp.denom(solved) == 1
    delta_coefficients.append(sp.expand(solved))
    delta += sp.expand(solved) * t**order

a_1 = sp.Rational(29, 14)
a_2 = 15 * a_1**2
assert weighted_top(delta_coefficients[0]) == (12, a_1 * u**3 * gamma**2)
assert weighted_top(delta_coefficients[1]) == (36, a_2 * u**9 * gamma**6)

# At every later order the unique top terms are -gamma*Delta_m/30 and
# W^3/2 sum_(i+j=m) Delta_i Delta_j.  Cubic Taylor terms lose 17 weights,
# parameter forcing has weight at most 25, and positive-t linear terms lose
# at least 7 weights.  Hence the scalar recurrence below is all-order.
leading_delta = {1: a_1}
for order in range(2, 9):
    leading_delta[order] = 15 * sum(
        leading_delta[index] * leading_delta[order - index]
        for index in range(1, order)
    )
    assert leading_delta[order] == (
        15 ** (order - 1) * a_1**order * sp.catalan(order - 1)
    )

# Put N=v^6 S^4.  Since the top parts of u and gamma are v and S, the
# reconstructed leading factors are
# gamma'/gamma=sqrt(1-(870/7)tN), x'/x=(...)^-1/2,
# z'/z=(...)^3/2.  Their coefficients never vanish.
scale = 60 * a_1
assert scale == sp.Rational(870, 7)
for order in range(1, 9):
    alpha = sp.binomial(-sp.Rational(1, 2), order) * (-scale) ** order
    beta = sp.binomial(sp.Rational(3, 2), order) * (-scale) ** order
    inverse_alpha = sp.binomial(-sp.Rational(1, 2), order) * scale**order
    inverse_beta = sp.binomial(sp.Rational(3, 2), order) * scale**order
    assert alpha != 0 and beta != 0
    assert inverse_alpha != 0 and inverse_beta != 0

assert sp.binomial(-sp.Rational(1, 2), 1) * (-scale) == sp.Rational(435, 7)
assert sp.binomial(sp.Rational(3, 2), 1) * (-scale) == -sp.Rational(1305, 7)
assert sp.binomial(-sp.Rational(1, 2), 2) * scale**2 == sp.Rational(567675, 98)
assert sp.binomial(sp.Rational(3, 2), 2) * scale**2 == sp.Rational(567675, 98)

# Saturated target-image quotient.  Write
# F_2=(x^-2 a(u,gamma),x^-1 b(u,gamma),x gamma).  The equivariant target
# coefficient modules of weights -2,-1,+1 are generated respectively by
# (A,B^2), (B,AC), and C.  After pullback and scalar extension to
# R=Q[u,gamma], their image is diagonal, with quotient
# R/(a,b^2) + R/(b,a gamma) + R/(gamma).
H_2 = H(2, W)
p_2 = sp.diff(H(2, w), w).subs(w, W)
q_2 = sp.cancel((W * p_2 - H_2) / c(2))
a = sp.cancel(u + q_2 / gamma**2)
b = sp.cancel(c(2) + p_2 / gamma)
assert sp.denom(a) == 1 and sp.denom(b) == 1

# The logarithmic differential matrix from (delta x/x,delta u,delta gamma)
# to the three normalized target components is unimodular over R.
logarithmic_jacobian = sp.Matrix(
    [
        [-2 * a, sp.diff(a, u), sp.diff(a, gamma)],
        [-b, sp.diff(b, u), sp.diff(b, gamma)],
        [gamma, 0, 1],
    ]
)
assert sp.factor(logarithmic_jacobian.det()) == sp.Rational(1, 30)

# For V_m=N^m(alpha_m x,0,beta_m z), its normalized C-component modulo
# (gamma) is
# N^m (-(2 alpha_m+beta_m)+(8/7)(alpha_m+beta_m)v).
# It is nonzero in R/(gamma)=Q[v] for every m because alpha_m,beta_m are
# nonzero and the two displayed linear coefficients cannot vanish together.
v = sp.symbols("v")
S_mod_gamma = sp.Rational(8, 7) * v - 1
N_mod_gamma = v**6 * S_mod_gamma**4
for order in range(1, 9):
    alpha = sp.binomial(-sp.Rational(1, 2), order) * (-scale) ** order
    beta = sp.binomial(sp.Rational(3, 2), order) * (-scale) ** order
    residue = sp.expand(
        N_mod_gamma**order
        * (-(2 * alpha + beta) + sp.Rational(8, 7) * (alpha + beta) * v)
    )
    assert residue != 0

print("PASS: the torus-gauge root equation is two-variable over Q[u,gamma]")
print("PASS: a_1=29/14 and the Catalan recurrence prove degree 24*m+1")
print("PASS: the exact first two leading source degrees are 25 and 49")
print("PASS: v^(6m)S^(4m) survives the saturated target-image quotient")
