#!/usr/bin/env python3
"""Explore the missing kappa=-1 chart one shear monomial at a time."""

import argparse
import sympy as sp


parser = argparse.ArgumentParser()
parser.add_argument("--shear-degree", type=int, default=2)
parser.add_argument("--x-degree", type=int, default=0)
args = parser.parse_args()

w, X, Q, Z, tau, coefficient = sp.symbols("w X Q Z tau coefficient")
rho = sp.symbols("rho")

H = sp.expand(
    w**2
    * (w - 1)
    * (
        tau * w**2
        + (sp.Rational(3, 2) - 2 * tau) * w
        + tau
        - sp.Rational(5, 2)
    )
)
p = sp.diff(H, w)
q = sp.expand(w * p - H)

psi = coefficient * X**args.x_degree * Q ** args.shear_degree
W = Z + psi
source_u = 1 + X * W
source_gamma = 1 + X**2 * Q
marked = sp.expand(source_u * source_gamma)
S = sp.cancel(
    (source_u + q.subs(w, marked) / source_gamma**2) / (2 * X**2)
)
T = sp.cancel((1 + p.subs(w, marked) / source_gamma) / X)
R = 2 * X * source_gamma

base_jacobian = sp.Matrix([S, T, R]).jacobian((X, Q, Z))
assert sp.factor(base_jacobian.det()) == -1


def quotient_bracket(left, right):
    return sp.cancel(
        2
        * X**3
        * (
            sp.diff(left, X) * sp.diff(right, Z)
            - sp.diff(left, Z) * sp.diff(right, X)
        )
        - (2 + 6 * X**2 * Q)
        * (
            sp.diff(left, Q) * sp.diff(right, Z)
            - sp.diff(left, Z) * sp.diff(right, Q)
        )
    )


assert quotient_bracket(S, T) == 1
assert quotient_bracket(R, S) == 0
assert quotient_bracket(R, T) == 0

SX, SQ, SZ = (sp.diff(S, variable) for variable in (X, Q, Z))
TX, TQ, TZ = (sp.diff(T, variable) for variable in (X, Q, Z))
w_family = sp.Matrix(
    [SZ * TQ - SQ * TZ, SX * TZ - SZ * TX, SQ * TX - SX * TQ]
).applyfunc(sp.cancel)
w_E = sp.Matrix(
    [(1 - 3 * X**2 * Q) / 2, 9 * X * Q**2 / 2, -6 * X * Q * Z]
)
difference = (w_family - w_E).applyfunc(sp.cancel)


def hamiltonian(f):
    return sp.Matrix(
        [
            -2 * X**3 * sp.diff(f, Z),
            (2 + 6 * X**2 * Q) * sp.diff(f, Z),
            2 * X**3 * sp.diff(f, X)
            - (2 + 6 * X**2 * Q) * sp.diff(f, Q),
        ]
    )


f_zero = sp.integrate(sp.cancel(difference[0] / (-2 * X**3)), Z)
residual = (difference - hamiltonian(f_zero)).applyfunc(sp.cancel)
assert residual[0] == residual[1] == 0

# At fixed rho=R, the remaining Hamiltonian vector field is
# 2*X^3*d/dX.  The transformed integrand is a finite Laurent polynomial, so
# integrate it termwise rather than invoking a general rational integrator.
Q_at_fixed_R = (rho - 2 * X) / (2 * X**3)
integrand = sp.cancel(
    residual[2].subs(Q, Q_at_fixed_R) / (2 * X**3)
)
numerator, denominator = sp.together(integrand).as_numer_denom()
denominator_degree = sp.Poly(denominator, X).degree()
denominator_unit = sp.factor(denominator / X**denominator_degree)
numerator_in_X = sp.Poly(numerator, X)
assert numerator_in_X.coeff_monomial(X ** (denominator_degree - 1)) == 0
h_fixed_R = sum(
    term_coefficient
    * X ** (monomial[0] - denominator_degree + 1)
    / (
        denominator_unit
        * (monomial[0] - denominator_degree + 1)
    )
    for monomial, term_coefficient in numerator_in_X.terms()
)
h = sp.cancel(h_fixed_R.subs(rho, R))
assert sp.cancel(residual[2] - hamiltonian(h)[2]) == 0
f = sp.cancel(f_zero + h)

f_numerator, f_denominator = sp.together(f).as_numer_denom()
negative_order = sp.Poly(f_denominator, X).degree()
f_numerator_in_X = sp.Poly(f_numerator, X)

A_tau = 2 * tau**2 - 15 * tau - 18
expected_by_degree = {
    0: {
        0: -2 * A_tau,
        1: -420 * coefficient,
        2: 8 * Q * A_tau,
        3: 1260 * Q * coefficient,
    },
    1: {
        0: 84 * coefficient,
        1: -2 * A_tau,
        2: -420 * Q * coefficient,
        3: 8 * Q * A_tau,
        4: 1260 * Q**2 * coefficient,
    },
    2: {
        0: -15 * coefficient,
        2: 105 * Q * coefficient,
        3: -2 * A_tau,
        4: -420 * Q**2 * coefficient,
        5: 8 * Q * A_tau,
        6: 1260 * Q**3 * coefficient,
    },
    3: {
        0: 28 * coefficient,
        2: -252 * Q * coefficient,
        4: 1260 * Q**2 * coefficient,
        5: -22 * A_tau,
        6: -4620 * Q**3 * coefficient,
        7: 88 * Q * A_tau,
        8: 13860 * Q**4 * coefficient,
    },
}
if args.x_degree == 0 and args.shear_degree in expected_by_degree:
    expected = expected_by_degree[args.shear_degree]
    for exponent in range(negative_order):
        actual = f_numerator_in_X.coeff_monomial(X**exponent)
        assert sp.factor(actual - expected.get(exponent, 0)) == 0

if (args.x_degree, args.shear_degree) == (1, 1):
    exceptional_obstruction = 105 * coefficient - 2 * A_tau
    assert sp.factor(f_denominator - 1680 * X**4) == 0
    assert sp.factor(
        f_numerator_in_X.coeff_monomial(X**0) - exceptional_obstruction
    ) == 0
    assert f_numerator_in_X.coeff_monomial(X**1) == 0
    assert sp.factor(
        f_numerator_in_X.coeff_monomial(X**2)
        + 4 * Q * exceptional_obstruction
    ) == 0
    assert f_numerator_in_X.coeff_monomial(X**3) == 0
    exceptional_candidate = 2 * A_tau / 105
    assert sp.factor(
        exceptional_obstruction.subs(coefficient, exceptional_candidate)
    ) == 0
    W_back = sp.symbols("W_back")
    Z_back = W_back - exceptional_candidate * X * Q
    assert sp.expand(
        (Z + exceptional_candidate * X * Q).subs(Z, Z_back) - W_back
    ) == 0

print("shear monomial = X^%d Q^%d" % (args.x_degree, args.shear_degree))
print("denominator =", sp.factor(f_denominator))
for exponent in range(negative_order):
    residue = sp.factor(f_numerator_in_X.coeff_monomial(X**exponent))
    if residue != 0:
        print(f"residue numerator [X^{exponent}] =", residue)
if args.x_degree == 0 and args.shear_degree in expected_by_degree:
    print("PASS: exact kappa=-1 monomial residue response")
if (args.x_degree, args.shear_degree) == (1, 1):
    print("PASS: X*Q shear completes every kappa=-1 seed")
