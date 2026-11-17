#!/usr/bin/env python3
"""Verify the one-pair Schur-completion obstruction for the Dvorsky cubic.

Pairing the unpaired variable t with one new variable s gives the
nondegenerate quadratic operator

    L = d_a d_d - d_b d_c + d_t d_s.

This checker first parametrizes every homogeneous cubic F with

    F(a,b,c,d,t,0) = (t+c)(a*d+b*t)

that satisfies L(F)=0.  It then computes the unrestricted two-jet of an
arbitrary polynomial or formal harmonic lift along s=0 and verifies

    L^2(F^2)|_(a=b=c=d=s=0) = 12*t^2 - 8*r_c*t.

Thus the coefficient of t^2 is always 12; no hyperplane lift for this
canonical six-variable quadratic operator satisfies the first two moments.
"""

from __future__ import annotations

import sympy as sp


a, b, c, d, t, s = sp.symbols("a b c d t s")
qaa, qab, qac, qad, qbb, qbc, qbd, qcc, qcd, qdd = sp.symbols(
    "qaa qab qac qad qbb qbc qbd qcc qcd qdd"
)
la, lb, lc, ld, r = sp.symbols("la lb lc ld r")


def laplacian(polynomial: sp.Expr) -> sp.Expr:
    """Apply the split nondegenerate quadratic operator once."""

    return sp.expand(
        sp.diff(polynomial, a, d)
        - sp.diff(polynomial, b, c)
        + sp.diff(polynomial, t, s)
    )


P = (t + c) * (a * d + b * t)

# Once L(F)=0 is imposed, every homogeneous cubic lift has this normal
# form.  R0 is an arbitrary quadratic in a,b,c,d.  Its D-Laplacian is the
# scalar delta=qad-qbc, which forces the displayed t-term in R1.
R0 = (
    qaa * a**2
    + qab * a * b
    + qac * a * c
    + qad * a * d
    + qbb * b**2
    + qbc * b * c
    + qbd * b * d
    + qcc * c**2
    + qcd * c * d
    + qdd * d**2
)
delta = qad - qbc
R2 = -c * t + R0
R1 = la * a + lb * b + lc * c + ld * d - sp.Rational(1, 2) * delta * t
F = sp.expand(P + s * R2 + s**2 * R1 + r * s**3)

assert sp.expand(F.subs(s, 0) - P) == 0
assert laplacian(F) == 0

second_moment = laplacian(laplacian(sp.expand(F**2)))
second_axis = sp.expand(second_moment.subs({a: 0, b: 0, c: 0, d: 0, s: 0}))
assert second_axis == 12 * t**2

coefficient = sp.Poly(second_moment, a, b, c, d, t, s).coeff_monomial(t**2)
assert coefficient == 12

# An arbitrary lift has expansion
#
#   F=P+s*(-c*t+R0)+s^2*(-t*D(R0)/2+S0)+O(s^3).
#
# Only the two-jets of R0 and S0 in a,b,c,d can affect the restriction of
# L^2(F^2) to the t-axis.  Parametrize those jets completely.
coordinates = (a, b, c, d)
rho0 = sp.symbols("rho0")
rho_linear = sp.symbols("rho_a rho_b rho_c rho_d")
rho_quadratic = sp.symbols(
    "rho_aa rho_ab rho_ac rho_ad rho_bb "
    "rho_bc rho_bd rho_cc rho_cd rho_dd"
)
sigma0 = sp.symbols("sigma0")
sigma_linear = sp.symbols("sigma_a sigma_b sigma_c sigma_d")
sigma_quadratic = sp.symbols(
    "sigma_aa sigma_ab sigma_ac sigma_ad sigma_bb "
    "sigma_bc sigma_bd sigma_cc sigma_cd sigma_dd"
)


def general_two_jet(
    constant: sp.Symbol,
    linear: tuple[sp.Symbol, ...],
    quadratic: tuple[sp.Symbol, ...],
) -> sp.Expr:
    answer = constant + sum(
        coefficient * variable
        for coefficient, variable in zip(linear, coordinates)
    )
    position = 0
    for left in range(len(coordinates)):
        for right in range(left, len(coordinates)):
            answer += (
                quadratic[position]
                * coordinates[left]
                * coordinates[right]
            )
            position += 1
    assert position == len(quadratic)
    return answer


R0_jet = general_two_jet(rho0, rho_linear, rho_quadratic)
S0_jet = general_two_jet(sigma0, sigma_linear, sigma_quadratic)
directional_laplacian_R0 = (
    sp.diff(R0_jet, a, d) - sp.diff(R0_jet, b, c)
)
jet_lift = sp.expand(
    P
    + s * (-c * t + R0_jet)
    + s**2 * (
        -sp.Rational(1, 2) * t * directional_laplacian_R0
        + S0_jet
    )
)
jet_first_moment = laplacian(jet_lift)
assert sp.expand(jet_first_moment.subs(s, 0)) == 0
assert sp.expand(sp.diff(jet_first_moment, s).subs(s, 0)) == 0

jet_second_moment = laplacian(laplacian(sp.expand(jet_lift**2)))
jet_second_axis = sp.factor(
    jet_second_moment.subs({a: 0, b: 0, c: 0, d: 0, s: 0})
)
assert jet_second_axis == 4 * t * (3 * t - 2 * rho_linear[2])
assert sp.expand(jet_second_axis).coeff(t, 2) == 12

print("PASS: the harmonic normal form covers every homogeneous one-pair lift")
print("PASS: the second pure moment has unavoidable t^2 coefficient 12")
print("PASS: every polynomial/formal hyperplane lift for this operator is excluded")
