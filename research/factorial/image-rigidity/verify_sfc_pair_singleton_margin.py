#!/usr/bin/env python3
"""Exact controls for the confluent Gaussian pair-plus-singleton margin.

No sampling supplies the uniform lower bound: its proof is in the note.
The controls derive the rows from Gaussian integrals and independently
check the polynomial certificate, including its coalescent endpoint.
"""

from math import comb

import sympy as sp


def main():
    A, B, c, delta, weight, q, kappa = sp.symbols(
        "A B c delta lambda q kappa")
    x, y = sp.symbols("x y")
    rows = [A + weight*c,
            A*A + q*B*B - q*delta*A*B + weight*c*c,
            A**3 + kappa*q*A*B*B - kappa*q*delta*A*A*B + weight*c**3]

    # E[phi_delta(Z_m)^j] = exp(-j(m-j) delta^2/(2m)).
    # For m=2 or 3 the mixed terms have just one value, x or y.
    for m in (1, 2, 3):
        off_diagonal = {1: 1, 2: x, 3: y}[m]
        raw = sum(comb(m, j) * (A-B/delta)**(m-j) * (B/delta)**j
                  * (1 if j in (0, m) else off_diagonal)
                  for j in range(m+1)) + weight*c**m
        exact_q = 2*(1-x)/delta**2
        exact_kappa = 3*(1-y)/(2*(1-x))
        expected = rows[m-1].subs({q: exact_q, kappa: exact_kappa})
        assert sp.cancel(raw-expected) == 0

    J = (1-kappa)*(A*A-weight*A*c+weight**2*c*c)-kappa*weight*c*c
    rho = weight*(weight+1)*(1+(kappa-1)*weight)
    assert sp.expand(rows[2]-kappa*A*rows[1]-J*rows[0]-rho*c**3) == 0

    # Direct confluent integration of (A+B Z_m)^m, Var Z_m=1/m.
    for m, row in enumerate(rows, 1):
        integral = sum(comb(m, j)*A**(m-j)*B**j
                       * (sp.factorial2(j-1)/sp.Integer(m)**(j//2) if j else 1)
                       for j in range(0, m+1, 2)) + weight*c**m
        endpoint = row.subs({delta: 0, q: sp.Rational(1, 2), kappa: 2})
        assert sp.expand(integral-endpoint) == 0

    # Dropping the derivative direction must destroy the endpoint quadratic.
    assert sp.expand(rows[1].subs({delta: 0, q: 0})
                     -(A*A+B*B/2+weight*c*c)) != 0
    print("Passed Gaussian-kernel row derivation, symbolic residual identity,")
    print("all three coalescent integrals, and a lost-derivative negative control.")
    print("Uniform margin: q_R^(3/2)/(512 (Lambda+1)^3 Lambda (Lambda^2+3 Lambda+4)).")


if __name__ == "__main__":
    main()
