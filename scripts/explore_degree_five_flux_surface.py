#!/usr/bin/env python3
"""Explore rank-two flux on the full normalized degree-five seed surface.

This is deliberately an exploratory calculation rather than a theorem
certificate.  The normalized seeds are parametrized by

    H=w^2(w-1)(t*w^2+(kappa/2-2*t+2)*w-kappa/2+t-3),

so that H'(1)=-1 and H''(1)=kappa.  The corresponding weighted coefficient
is a=-(1+kappa)/(2+kappa).  We use the seed-dependent diagonal source change
which sends 2*C to the *fixed* quotient coordinate R=2*X-3*X^2*Q, and then
compute the localized relative Hamiltonian and its negative-X residues.
"""

import argparse
import sympy as sp


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--psi-degree", type=int, default=2)
    parser.add_argument("--kappa-value")
    parser.add_argument("--tau-value")
    args = parser.parse_args()

    w = sp.symbols("w")
    x, y, z = sp.symbols("x y z")
    X, Q, Z = sp.symbols("X Q Z")
    kappa, tau = sp.symbols("kappa tau")
    shear = {
        degree: sp.symbols(f"s{degree}")
        for degree in range(2, args.psi_degree + 1)
    }

    H = sp.expand(
        w**2
        * (w - 1)
        * (
            tau * w**2
            + (kappa / 2 - 2 * tau + 2) * w
            - kappa / 2
            + tau
            - 3
        )
    )
    p = sp.diff(H, w)
    q = sp.expand(w * p - H)
    assert sp.expand(H.subs(w, 1)) == 0
    assert sp.expand(p.subs(w, 1)) == -1
    assert sp.expand(sp.diff(H, w, 2).subs(w, 1)) == kappa

    parameter_values = {}
    if args.kappa_value is not None:
        parameter_values[kappa] = sp.sympify(args.kappa_value)
    if args.tau_value is not None:
        parameter_values[tau] = sp.sympify(args.tau_value)
    if parameter_values:
        H = sp.expand(H.subs(parameter_values))
        p = sp.diff(H, w)
        q = sp.expand(w * p - H)

    a = sp.cancel((-(1 + kappa) / (2 + kappa)).subs(parameter_values))
    # Constant and linear terms are suppressed: they belong to the elementary
    # R-preserving base normalization that the residue is meant to quotient.
    psi = sum(shear[degree] * Q**degree for degree in shear)
    W = Z + psi
    Y = Q - X * W / 3
    # Work directly with the weighted invariants after substitution.  This
    # avoids expanding the source map before the decisive cancellations.
    source_v = -3 * X * Y / (2 * a)
    source_gamma = 1 - 3 * X * Q / 2
    source_u = 1 + source_v
    marked = sp.expand(source_u * source_gamma)
    A_sub = (source_u + q.subs(w, marked) / source_gamma**2) / X**2
    B_sub = (1 + p.subs(w, marked) / source_gamma) / X
    S = sp.cancel(-2 * a * A_sub / 3)
    T = sp.cancel(B_sub)
    R = sp.expand(2 * X * source_gamma)
    assert sp.expand(R - (2 * X - 3 * X**2 * Q)) == 0

    SX, SQ, SZ = (sp.diff(S, variable) for variable in (X, Q, Z))
    TX, TQ, TZ = (sp.diff(T, variable) for variable in (X, Q, Z))
    # The determinant is -1 by the weighted Jacobian theorem and the chosen
    # diagonal scalings.  Writing the relevant inverse-Jacobian column as
    # cofactors is dramatically cheaper than a generic symbolic determinant.
    w_family = sp.Matrix(
        [SZ * TQ - SQ * TZ, SX * TZ - SZ * TX, SQ * TX - SX * TQ]
    ).applyfunc(sp.cancel)
    w_E = sp.Matrix([(1 + 3 * X * Q) / 2, -3 * Q**2, 9 * Q * Z / 2])
    difference = (w_family - w_E).applyfunc(sp.cancel)

    def hamiltonian(f):
        return sp.Matrix(
            [
                3 * X**2 * sp.diff(f, Z),
                (2 - 6 * X * Q) * sp.diff(f, Z),
                -3 * X**2 * sp.diff(f, X)
                + (6 * X * Q - 2) * sp.diff(f, Q),
            ]
        )

    f_zero = sp.integrate(sp.cancel(difference[0] / (3 * X**2)), Z)
    residual = (difference - hamiltonian(f_zero)).applyfunc(sp.cancel)
    assert residual[0] == residual[1] == 0

    inv_x, rho = sp.symbols("inv_x rho")
    residual_vrho = sp.cancel(
        residual[2].subs({X: 1 / inv_x, Q: inv_x * (2 - rho * inv_x) / 3})
    )
    h_vrho = sp.integrate(residual_vrho / 3, inv_x)
    h = sp.cancel(
        h_vrho.subs({inv_x: 1 / X, rho: 2 * X - 3 * X**2 * Q})
    )
    f = sp.cancel(f_zero + h)
    assert all(
        sp.cancel(difference[index] - hamiltonian(f)[index]) == 0
        for index in range(3)
    )

    numerator, denominator = sp.together(f).as_numer_denom()
    denominator_x_degree = sp.Poly(denominator, X).degree()
    numerator_x = sp.Poly(numerator, X)
    residues = [
        sp.factor(numerator_x.coeff_monomial(X**index))
        for index in range(denominator_x_degree)
    ]
    denominator_unit = sp.factor(denominator / X**denominator_x_degree)

    print("normalized seed H =", sp.factor(H))
    print("weighted a =", a)
    print("psi =", psi)
    print("denominator = X^%d * (%s)" % (denominator_x_degree, denominator_unit))
    for index, residue in enumerate(residues):
        print(f"residue numerator [X^{index}] =", residue)
    print("nonzero residue count =", sum(residue != 0 for residue in residues))


if __name__ == "__main__":
    main()
