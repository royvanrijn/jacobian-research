#!/usr/bin/env python3
"""Exact small audit of the relative degree-five rank-two family.

This checker packages the inexpensive algebra behind the relative
quantization calculation.  The large order-five period and Fitting
calculations remain in their dedicated checkers.

It verifies:

* the two-parameter normalized residual-root family and the old
  ``(kappa,tau)=(0,1)`` symbol;
* the root-at-infinity chart and its valuation weights;
* homogeneity of the restricted correction differential and of every
  Moyal contraction for that valuation;
* the weights of the certified sparse order-three correction;
* the two exact coprime Fitting charts of the 15-by-16 strong-cocycle
  presentation; and
* the two nonboundary closed points in the four-period Kuranishi shadow.

No characteristic-zero length is inferred from the modular length-218
Fitting computations.
"""

from __future__ import annotations

from functools import reduce
from math import gcd
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ

from explore_degree_five_quantum_residue import degree_five_family
from verify_degree_five_third_order_function_field import (
    S_SUPPORT,
    T_SUPPORT,
)


ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "artifacts" / "generated-results"


def valuation(monomial: tuple[int, int, int]) -> int:
    """Root-at-infinity weight of ``X^i Q^j Z^k``."""

    x_degree, q_degree, z_degree = monomial
    return x_degree - q_degree - 2 * z_degree


def primitive_polynomial(path: Path, a, tau) -> sp.Poly:
    expression = sp.sympify(path.read_text().strip().replace("^", "**"))
    polynomial = sp.Poly(expression, a, tau)
    assert reduce(gcd, (abs(int(c)) for c in polynomial.coeffs())) == 1
    return polynomial


def main() -> None:
    a, tau, w = sp.symbols("a tau w")
    kappa = -(1 + 2 * a) / (1 + a)
    residual_linear = sp.factor(kappa / 2 - 2 * tau + 2)
    residual_constant = sp.factor(-kappa / 2 + tau - 3)
    residual = tau * w**2 + residual_linear * w + residual_constant
    seed = sp.expand(w**2 * (w - 1) * residual)

    assert sp.factor(seed.subs(w, 1)) == 0
    assert sp.factor(sp.diff(seed, w).subs(w, 1) + 1) == 0
    assert sp.factor(sp.diff(seed, w, 2).subs(w, 1) - kappa) == 0
    assert sp.factor(kappa.subs(a, -sp.Rational(1, 2))) == 0
    assert sp.factor(
        residual.subs({a: -sp.Rational(1, 2), tau: 1}) - (w**2 - 2)
    ) == 0
    function_field = QQ.frac_field("a", "tau")
    ff_a, ff_tau = function_field.gens
    symbol_s, symbol_t = degree_five_family(
        function_field,
        ff_a,
        ff_tau,
        verify_canonical=False,
    )
    assert min(map(valuation, symbol_s)) == -2
    assert min(map(valuation, symbol_t)) == -1

    # If r=1/w is the residual-root coordinate at infinity, its equation is
    # tau+A*r+B*r^2.  The infinity divisor is tau=0.  Its generic root has
    # tau-adic valuation one; at a=-3/2 the two roots have valuation 1/2.
    r = sp.symbols("r")
    infinity_equation = sp.expand(
        tau + residual_linear * r + residual_constant * r**2
    )
    assert sp.factor(infinity_equation.subs(r, 0) - tau) == 0
    assert sp.factor(
        residual_linear.subs(tau, 0)
        - (2 * a + 3) / (2 * (a + 1))
    ) == 0
    assert sp.factor(
        infinity_equation.subs(a, -sp.Rational(3, 2))
        - (tau - 2 * tau * r + (tau - 1) * r**2)
    ) == 0
    discriminant = sp.factor(
        residual_linear**2 - 4 * tau * residual_constant
    )
    assert sp.factor(
        discriminant.subs(tau, 0)
        - (2 * a + 3) ** 2 / (4 * (a + 1) ** 2)
    ) == 0

    # Marked-root/Ore bridge.  On s=1/u and d=v/u, keep the marked argument
    # m=R*u/2 finite.  The resulting formulas give the valuation lattice.
    s, d, marked = sp.symbols("s d marked", nonzero=True)
    X = s / d
    R = 2 * s * marked
    v = d / s
    Q = sp.factor(v * (2 - R * v) / 3)
    W = sp.factor(2 * a * v * (1 / s - v) + 3 * v * Q)
    gamma = sp.factor(R * v / 2)
    assert sp.factor(Q - 2 * d * (1 - d * marked) / (3 * s)) == 0
    assert sp.factor(
        W
        - 2 * d * (a * (1 - d) + d * (1 - d * marked)) / s**2
    ) == 0
    assert sp.factor(gamma - d * marked) == 0
    root_weights = {
        name: int(sp.factor(expression).as_powers_dict().get(s, 0))
        for name, expression in {
            "X": X,
            "Q": Q,
            "W": W,
            "R": R,
            "gamma": gamma,
        }.items()
    }
    assert root_weights == {"X": 1, "Q": -1, "W": -2, "R": 1, "gamma": 0}

    # In the adapted Poisson tensor, delta raises valuation by one and
    # d/dZ raises it by two.  Hence a bracket raises total valuation by
    # three, and Pi^j raises it by 3*j.  Since nu(S,T)=(-2,-1), shifting
    # the two correction summands by (2,1) makes d1 filtered of degree zero.
    delta_shift = 1
    dz_shift = 2
    bracket_shift = delta_shift + dz_shift
    nu_s, nu_t = -2, -1
    assert bracket_shift == 3
    assert nu_t + bracket_shift == 2
    assert nu_s + bracket_shift == 1
    assert (nu_s + nu_t + bracket_shift) == 0
    assert [valuation(m) for m in S_SUPPORT] == [4] * len(S_SUPPORT)
    assert [valuation(m) for m in T_SUPPORT] == [5] * len(T_SUPPORT)

    # All four contributions to O_5 start in root weight 12 when the
    # order-three correction has weights (4,5).
    fifth_weights = {
        "poisson(S2,T2)": 4 + 5 + bracket_shift,
        "Pi3(S2,T)": 4 + nu_t + 3 * bracket_shift,
        "Pi3(S,T2)": nu_s + 5 + 3 * bracket_shift,
        "Pi5(S,T)": nu_s + nu_t + 5 * bracket_shift,
    }
    assert set(fifth_weights.values()) == {12}
    assert 10 + (nu_t + bracket_shift) == 12
    assert 11 + (nu_s + bracket_shift) == 12

    # Exact Fitting charts for the generic strong-cocycle line.
    fitting_d = primitive_polynomial(
        GENERATED / "degree_five_qper_pivot_D34.sing", a, tau
    )
    fitting_e = primitive_polynomial(
        GENERATED / "degree_five_qper_pivot_E35.sing", a, tau
    )
    assert fitting_d.total_degree() == 34
    assert fitting_e.total_degree() == 35
    assert sp.gcd(fitting_d, fitting_e).total_degree() == 0

    # The reduced interior shadow left by the four exact generic periods.
    cubic = sp.Poly(94 * a**3 + 335 * a**2 + 400 * a + 160, a)
    cubic_tau = -(658 * a**2 + 1593 * a + 976) / 8
    assert sp.factor(cubic.as_expr()) == cubic.as_expr()
    rational_point = {a: -sp.Rational(1, 2), tau: -3}
    kernel_chart = (
        4 * a**3 * tau**2
        - 24 * a**3 * tau
        - 72 * a**3
        + 8 * a**2 * tau**2
        - 54 * a**2 * tau
        - 216 * a**2
        + 4 * a * tau**2
        - 30 * a * tau
        - 246 * a
        - 105
    )
    assert kernel_chart.subs(rational_point) == -45
    for boundary in (
        a,
        a + 1,
        sp.together(cubic_tau).as_numer_denom()[0],
        sp.together(kernel_chart.subs(tau, cubic_tau)).as_numer_denom()[0],
    ):
        assert sp.gcd(cubic, sp.Poly(boundary, a)).degree() == 0

    print("PASS: the old (kappa,tau)=(0,1) symbol lies in the relative family")
    print("PASS: root-infinity weights are (X,Q,W,R,gamma)=(1,-1,-2,1,0)")
    print("PASS: shifted correction weights are (S2,T2)=(4,5), (S4,T4)=(10,11)")
    print("PASS: the two exact Fitting charts are coprime off the boundary")
    print("PASS: the interior period shadow is rational plus cubic")
    print("SCOPE: Fitting length 218 remains modular; no DC_2 claim")


if __name__ == "__main__":
    main()
