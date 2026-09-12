#!/usr/bin/env python3
"""Exact principal-part obstructions for natural nonlinear graph shears.

Three all-degree nonlinear potential classes are tested:

1. B=M+grad V(Q) in the natural chart Q=(x,q,R,T);
2. chart 1110 with V depending on the exceptional variables (D,T); and
3. chart 1111 with V depending on the exceptional variables (D,S).

In the first class the D^2 coefficient forces V_TT to have a pole
-3/(16X).  In the two exceptional classes, restriction to X=0 forces a
one-variable Hessian datum, and comparison away from X=0 contradicts that
dependence.  No polynomial potential in any of these classes yields constant
nonzero determinant.
"""

from __future__ import annotations

import runpy

import sympy as sp


namespace = runpy.run_path("scripts/search_hc4_graph_polarizations.py")
X, Y, W, D = namespace["h_variables"]
T = namespace["Th"]
S = namespace["Sh"]

# The two recurring projection minors.
A0 = (
    3 * W * X**4 * Y**2
    + 6 * W * X**3 * Y
    + 3 * W * X**2
    + 9 * X**3 * Y**3
    + 12 * X**2 * Y**2
    - X * Y
    - 3
)
C0 = (
    3 * W * X**3 * Y
    + 3 * W * X**2
    + 9 * X**2 * Y**2
    + 3 * X * Y
    - 4
)
rational_forced_term = sp.cancel((A0 + 3) / (4 * X * C0))


def main() -> None:
    # Natural chart: the coefficient of D^2 in det d(M+grad V(Q)) is
    # -6 X^2 (A0+4 X C0 V_TT).  Constancy forces the bracket to vanish.
    v_tt = sp.symbols("v_tt")
    d_squared_bracket = sp.expand(A0 + 4 * X * C0 * v_tt)
    assert d_squared_bracket.subs({X: 0}) == -3
    forced_v_tt = sp.cancel(-A0 / (4 * X * C0))
    assert sp.limit(X * forced_v_tt, X, 0) == -sp.Rational(3, 16)
    print(
        "PASS natural chart: constancy forces "
        "V_TT=-A0/(4*X*C0) with principal part -3/(16X)"
    )

    # Exceptional chart 1111, V=V(D,S).  If A=V_DD and Delta=det Hess V,
    # the determinant is (Delta*A0-4*A*X*C0)/6.  At X=0, constant nonzero
    # determinant forces Delta=kappa != 0.  Then A=kappa*rational_forced_term.
    # Its X=0 boundary value is kappa*Y/16, but S=(W+4Y^2)/2 there; fixing S
    # and varying Y contradicts A being a function of (D,S).
    assert sp.factor(S.subs({X: 0}) - (W + 4 * Y**2) / 2) == 0
    assert sp.limit(rational_forced_term, X, 0) == Y / 16
    print(
        "PASS chart 1111, V(D,S): the forced V_DD boundary value Y/16 "
        "does not factor through S=(W+4Y^2)/2"
    )

    # Exceptional chart 1110, V=V(D,T).  Constant nonzero determinant first
    # forces A=V_DD=kappa != 0, then Delta=-kappa*rational_forced_term.
    # If the latter factored through T, its X=0 value would force it to equal
    # -kappa*T/16 identically.  The following nonzero numerator disproves that.
    assert T.subs({X: 0}) == Y
    discrepancy_numerator = sp.factor(
        sp.together(rational_forced_term - T / 16).as_numer_denom()[0]
    )
    assert discrepancy_numerator != 0
    assert discrepancy_numerator.free_symbols
    print(
        "PASS chart 1110, V(D,T): the forced Hessian determinant does not "
        "factor through T"
    )

    print(
        "PASS: all-degree nonlinear shear obstructions hold in the natural "
        "chart and the two exceptional two-variable classes"
    )


if __name__ == "__main__":
    main()
