"""Exact identity checks and bounded controls for the affine isotropic chart.

The written proof supplies completeness and unrestricted degree coverage.
This checker supplies symbolic identity checks and explicitly bounded controls.
Run from the repository root with .venv/bin/python and --output PATH to retain
a receipt. No receipt is written unless that option is supplied.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from fractions import Fraction
from pathlib import Path

import sympy as sp


Z, W, Y, T = sp.symbols("z w y t")


def zero(expr):
    assert sp.cancel(sp.expand(expr)) == 0, expr


def expectation(expr):
    result = 0
    for (p, q, r), coefficient in sp.Poly(sp.expand(expr), Z, W, Y).terms():
        if p == q and r % 2 == 0:
            result += coefficient * math.factorial(p) * math.factorial(r) / (
                2 ** (r // 2) * math.factorial(r // 2)
            )
    return sp.expand(result)


def sparse(expr):
    return {powers: Fraction(coefficient) for powers, coefficient in
            sp.Poly(expr, Z, W, Y).terms()}


def multiply(left, right):
    result = {}
    for alpha, a in left.items():
        for beta, b in right.items():
            gamma = tuple(x + y for x, y in zip(alpha, beta))
            result[gamma] = result.get(gamma, Fraction(0)) + a * b
    return {gamma: value for gamma, value in result.items() if value}


def wick(poly, marker=(0, 0, 0)):
    result = Fraction(0)
    for powers, coefficient in poly.items():
        p, q, r = (a + b for a, b in zip(powers, marker))
        if p == q and r % 2 == 0:
            result += coefficient * math.factorial(p) * math.factorial(r) / (
                2 ** (r // 2) * math.factorial(r // 2)
            )
    return result


def radial_controls():
    """Compare the composed divergent-series formula with direct factorial Wick."""
    radial, y, t = sp.symbols("R y t")
    cases = [(1, 2, 3, -1, 2, -4), (0, 1, 0, 0, 0, 0),
             (1, 0, 0, 0, 0, -2), (0, 0, 2, -2, 1, 0)]
    order = 8
    for params in cases:
        alpha, beta, a, gamma, b, k = map(sp.Rational, params)
        phase = alpha * radial**2 + beta * radial * y + a * y**2 + gamma * radial + b * y + k
        d = 1 - 2 * a * t
        quadratic = alpha * t + beta**2 * t**2 / (2 * d)
        linear = gamma * t + beta * b * t**2 / d
        scalar = k * t + b**2 * t**2 / (2 * d)
        argument = sp.series(quadratic / (1 - linear)**2, t, 0, order + 1).removeO()
        composition = sp.Integer(1)
        argument_power = sp.Integer(1)
        for n in range(1, order + 1):
            argument_power = sp.series(argument_power * argument, t, 0, order + 1).removeO()
            composition += sp.Rational(math.factorial(2*n), math.factorial(n)) * argument_power
        prefactor = sp.series(sp.exp(scalar) / (sp.sqrt(d) * (1-linear)), t, 0, order+1).removeO()
        formula = sp.series(prefactor * composition, t, 0, order+1).removeO().expand()
        power = sp.Integer(1)
        for m in range(order+1):
            contraction = sp.Integer(0)
            for (i, j), coefficient in sp.Poly(power, radial, y).terms():
                if j % 2 == 0:
                    contraction += coefficient * math.factorial(i) * sp.Rational(
                        math.factorial(j), 2**(j//2) * math.factorial(j//2))
            zero(formula.coeff(t, m) * math.factorial(m) - contraction)
            power = sp.expand(power * phase)
    return {"parameter_cases": cases, "moment_orders": [0, order]}


def checks():
    c, d = sp.symbols("c d", nonzero=True)
    ell = sp.Function("ell")(Z)
    h = c + d * Z
    v = -d - 3 * d**2 * Z / (2 * c) - d**3 * Z**2 / (2 * c**2)
    b = h**2 * ell
    k = -c**2 * Z * h * ell**2 / 2
    denominator = h - 2 * Z * v
    zero(denominator - h**3 / c**2)
    zero((h - Z * sp.diff(h, Z))**2 * denominator - h**3)
    zero(2 * denominator * k + Z * b**2)
    root = c * T / (1 - d * T)
    zero(root - T * h.subs(Z, root))
    j = 1 - d * T
    precision = 1 - 2 * T * v.subs(Z, root)
    zero(j**2 * precision - 1)

    alpha, beta, gamma = sp.symbols("alpha beta gamma")
    phase = W * h + v * Y**2 + b * Y + k
    def direction(expr):
        return (alpha * sp.diff(expr, Z) + beta * sp.diff(expr, W)
                + gamma * sp.diff(expr, Y))
    second = sp.expand(direction(direction(phase)))
    zero(second.coeff(Y, 2) + alpha**2 * d**3 / c**2)
    zero(second.subs(alpha, 0) - 2 * gamma**2 * v)

    abstract_b = sp.Function("b")(Z)
    forward_y = Y - Z * abstract_b / c
    forward_w = W + abstract_b * Y / c - Z * abstract_b**2 / (2 * c**2)
    inverse_y = Y + Z * abstract_b / c
    inverse_w = W - abstract_b * Y / c - Z * abstract_b**2 / (2 * c**2)
    zero(2 * Z * forward_w + forward_y**2 - 2 * Z * W - Y**2)
    zero(forward_y.subs({Y: inverse_y, W: inverse_w}, simultaneous=True) - Y)
    zero(forward_w.subs({Y: inverse_y, W: inverse_w}, simultaneous=True) - W)

    # Mixed constant source: phase expansion and the exact moving-coordinate bridge.
    lam = sp.symbols("lambda", nonzero=True)
    a = sp.symbols("a")
    abstract_h, abstract_u, abstract_k = (
        sp.Function(name)(Z) for name in ("h", "u", "k")
    )
    s = (Z / T - abstract_h) / lam
    phi = -s**2 / 2 + T * (a * s**2 + abstract_u * s + abstract_k)
    expanded_phi = (
        -Z**2 / (2 * lam**2 * T**2)
        + (Z * abstract_h + a * Z**2) / (lam**2 * T)
        - abstract_h**2 / (2 * lam**2)
        - 2 * a * Z * abstract_h / lam**2 + Z * abstract_u / lam
        + T * (a * abstract_h**2 / lam**2 - abstract_u * abstract_h / lam + abstract_k)
    )
    zero(phi - expanded_phi)
    j_z = 1 - T * sp.diff(abstract_h, Z)
    r_y = T * lam / j_z
    r_t = Z / (T * j_z)  # root equation: h(r)+lambda Y=r/t
    # F_Y = -Y + phi_z r_Y + Y; check the spatial Stein identity at s=Y.
    f_z = T * (sp.diff(abstract_u, Z) * Y + sp.diff(abstract_k, Z))
    f_y = T * (2 * a * Y + abstract_u) + f_z * r_y
    zero(sp.diff(phi, Z).subs(Z / T - abstract_h, lam * Y) * r_y - (f_y - Y))
    # Check r_Yt=r_tY along the root and H'/H=-1/t.
    total_t_r_y = sp.diff(r_y, T) + sp.diff(r_y, Z) * r_t
    total_y_r_t = sp.diff(r_t, Z) * r_y
    zero(total_t_r_y - total_y_r_t)

    # General mixed source: spatial relation before logarithmic specialization.
    abstract_c = sp.Function("c")(Z)
    general_s = (Z / T - abstract_h) / abstract_c
    general_phi = -general_s**2 / 2 + T * (
        a * general_s**2 + abstract_u * general_s + abstract_k
    )
    general_j = 1 - T * (sp.diff(abstract_h, Z) + sp.diff(abstract_c, Z) * Y)
    general_r_y = T * abstract_c / general_j
    general_f_y = T * (2 * a * Y + abstract_u) + f_z * general_r_y
    zero(sp.diff(general_phi, Z).subs(
        Z / T - abstract_h, abstract_c * Y
    ) * general_r_y - (general_f_y - Y))
    # Explicit time-volume identity for a general amplitude q(z).
    amp = sp.Function("Q")(Z)
    general_r_t = Z / (T * general_j)
    volume = 1 / (T * abstract_c)
    # With F_t along the root = phi_t + phi_z*r_t, cancellation reduces to this.
    lhs_remainder = (
        (sp.diff(volume, T) + sp.diff(volume, Z) * general_r_t)
        * general_r_y * amp
        + volume * (sp.diff(general_r_y, T)
                    + sp.diff(general_r_y, Z) * general_r_t) * amp
        + volume * general_r_y * sp.diff(amp, Z) * general_r_t
        + volume * general_r_y * amp / T
    )
    rhs_remainder = sp.diff(volume * general_r_t * amp, Z) * general_r_y + (
        sp.diff(volume * general_r_t * amp, Y)
    )
    zero(lhs_remainder - rhs_remainder)

    # Root-zero tilt: exact reduction of the transformed denominator.
    eta, gamma, u0 = sp.symbols("eta gamma u0")
    delta_root = sp.symbols("sqrt_delta", nonzero=True)
    mu = T * u0 / delta_root**2
    denom = 1 - T * eta - T * gamma * mu
    qtilt = T * gamma / (delta_root * denom)
    zero(1 - T * (eta + gamma * (mu + Y / delta_root)) - denom * (1 - qtilt * Y))

    # Fixed logarithmic volume: partial fractions at several simple punctures.
    squarefree_cases = [Z*(Z-1), (Z+2)*(Z-1)*(Z-3),
                        (Z+3)*(Z+1)*(Z-2)*(Z-4), Z*(Z+1)*(Z+2)*(Z-2)*(Z-3)]
    for polynomial in squarefree_cases:
        roots = sp.polys.polytools.ground_roots(polynomial, Z)
        assert all(multiplicity == 1 for multiplicity in roots.values())
        zero(1/polynomial - sum(1/(sp.diff(polynomial, Z).subs(Z, root) * (Z-root))
                               for root in roots))
    # A repeated-root rejection control: this nonzero phase is pure by weight.
    assert sp.discriminant(Z**2, Z) == 0
    repeated_phase = sparse(Z**2 * W * Y)
    repeated_power = {(0, 0, 0): Fraction(1)}
    for m in range(1, 9):
        repeated_power = multiply(repeated_power, repeated_phase)
        assert wick(repeated_power) == 0

    # Repeated-square reciprocal phase: exact top coefficient and t-infinity powers.
    xr, reciprocal_r = sp.symbols("x r", nonzero=True)
    kappa = sp.symbols("kappa", nonzero=True)
    root_alpha, h0, h1, h2, h3, h4 = sp.symbols("root_alpha h0 h1 h2 h3 h4")
    u0x, u1x, u2x, k0x, k1x = sp.symbols("u0x u1x u2x k0x k1x")
    reciprocal_s = ((root_alpha*reciprocal_r-h0)*xr**2
                    + (reciprocal_r-h1)*xr - h2 - h3/xr - h4/xr**2) / kappa
    reciprocal_phi = sp.expand(-reciprocal_s**2/2 + (
        a*reciprocal_s**2 + (u0x+u1x/xr+u2x/xr**2)*reciprocal_s + k0x+k1x/xr
    )/reciprocal_r)
    top_reciprocal = -(1-2*a/reciprocal_r)*(root_alpha*reciprocal_r-h0)**2/(2*kappa**2)
    zero(reciprocal_phi.coeff(xr, 4) - top_reciprocal)
    # Multiplication by r clears every parameter pole, irrespective of x powers.
    sp.Poly(sp.expand(reciprocal_r*xr**4*reciprocal_phi), reciprocal_r, xr)
    nu, vv = sp.symbols("u v", nonzero=True)
    moving_centre = sp.Function("xi")(nu*vv)
    zero(nu*sp.diff(moving_centre, nu) - vv*sp.diff(moving_centre, vv))

    # Algebraic coefficients used in the polynomial quadratic-form orbit exclusion.
    fcoord, gcoord, leading, constant = sp.symbols("F G A B")
    quadform = Y**2 + 2*Z*W
    transformed = (1+fcoord) * (quadform-gcoord**2)/(2*fcoord) - (1+fcoord)*(2+fcoord)*gcoord**2/2
    zero(transformed - (1+fcoord)/(2*fcoord) * (quadform-(1+fcoord)**2*gcoord**2))
    coeff_y2 = sp.expand(transformed.subs(gcoord, leading*Y+constant)).coeff(Y, 2)
    zero(coeff_y2 - (1+fcoord)*(1-(1+fcoord)**2*leading**2)/(2*fcoord))
    orbit_f = leading*Y + constant
    orbit_h = Y/(2*leading) + sp.symbols("E")
    orbit_polynomial = (1+orbit_f)*orbit_h - (1+orbit_f)*(2+orbit_f)*gcoord**2/2
    zero(sp.expand(orbit_polynomial).coeff(Y, 2) - (1-leading**2*gcoord**2)/2)

    alpha, beta, radial_a = sp.symbols("alpha beta radial_a")
    a_numerator = sp.expand(alpha*T*(1-2*radial_a*T) + beta**2*T**2/2)
    zero(a_numerator.coeff(T, 1) - alpha)
    zero(a_numerator.subs(alpha, 0).coeff(T, 2) - beta**2/2)
    radial_result = radial_controls()

    # Exact all-order polynomials from the shear; direct Wick controls are bounded.
    shear_phase = W + Z**2 * Y - Z**5 / 2
    transformed_w = W - Z**2 * Y - Z**5 / 2
    marker_periods = {
        q: expectation((transformed_w**q).subs(Z, Z + T))
        for q in range(1, 5)
    }
    expected = {
        1: -T**5 / 2,
        2: T**10 / 4 - 4 * T**4,
        3: -T**15 / 8 + 6 * T**9 - 18 * T**3,
        4: T**20 / 16 - 6 * T**14 + 84 * T**8 - 48 * T**2,
    }
    for q in marker_periods:
        zero(marker_periods[q] - expected[q])
        assert sp.degree(marker_periods[q], T) == 5 * q
    power = {(0, 0, 0): Fraction(1)}
    base = sparse(shear_phase)
    for m in range(33):
        if m:
            assert wick(power) == 0
        for q, period in marker_periods.items():
            assert wick(power, (0, q, 0)) == (
                sp.expand(period).coeff(T, m) * math.factorial(m)
            )
        power = multiply(power, base)

    # Long's phase: pure through the stated control range, one fixed unsafe marker.
    long_phase = (1 + Z) * (W - (2 + Z) * Y**2 / 2)
    power = {(0, 0, 0): Fraction(1)}
    base = sparse(long_phase)
    for m in range(9):
        assert wick(power) == (1 if m == 0 else 0)
        assert wick(power, (1, 0, 0)) == (0 if m == 0 else math.factorial(m))
        power = multiply(power, base)

    # A positive-weight example with nonconstant v.
    weight_phase = W * (Z**2 + Z**4) + Z**3 * Y**2 + (Z + Z**2) * Y + Z**5
    power = {(0, 0, 0): Fraction(1)}
    base = sparse(weight_phase)
    for m in range(8):
        assert wick(power) == (1 if m == 0 else 0)
        if m > 2:
            assert wick(power, (0, 2, 0)) == 0
        power = multiply(power, base)

    # Rejection controls: changing the shear correction or adding constant Y².
    assert expectation((W + Z**2 * Y - Z**5 / 3)**6) != 0
    assert expectation(Y**2 + Z * W) != 0
    return {
        "symbolic": [
            "arbitrary-ell unsafe-branch rational identities",
            "root equation and determinant cancellation",
            "directional second derivative and isotropic-axis obstruction",
            "shear inverse and Gaussian quadratic-form preservation",
            "mixed constant-source phase expansion and moving-coordinate identities",
            "general mixed-source spatial and time-volume identities",
            "root-zero Gaussian tilt denominator identity",
            "multipuncture squarefree partial-fraction identities",
            "repeated-square reciprocal quartic and moving-centre derivative identities",
            "Long polynomial-orbit coefficient identities",
            "radial quartic quadratic-coefficient numerator",
            "four exact all-order marker polynomials from the shear",
        ],
        "bounded_controls": {
            "shear_B_z_squared": {"moment_orders": [0, 32], "W_marker_powers": [1, 4]},
            "long_phase": {"moment_orders": [0, 8], "markers": ["1", "Z"]},
            "positive_weight_example": {"moment_orders": [0, 7], "markers": ["1", "W^2"]},
            "rejection_controls": 2,
            "squarefree_partial_fractions": {"degrees": [2, 3, 4, 5]},
            "repeated_root_pure_rejection_control": {"phase": "Z^2 W Y", "moment_orders": [1, 8]},
            "radial_quartic_formula": radial_result,
        },
        "marker_periods_B_z_squared": {str(q): str(p) for q, p in marker_periods.items()},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = checks()
    here = Path(__file__).resolve()
    note = here.with_name("GAUSSIAN_AFFINE_ISOTROPIC_CHART.md")
    receipt = {
        "status": "passed",
        "scope": "Exact symbolic identities and bounded Wick controls; completeness and arbitrary-degree classification are proved in the written note, not by these controls.",
        "checked_on": "2026-11-17",
        "command": ".venv/bin/python research/moments/vanishing/verify_gaussian_affine_isotropic_chart.py --output research/moments/vanishing/gaussian_affine_isotropic_chart_verification.json",
        "environment": {"python": platform.python_version(), "sympy": sp.__version__},
        "checker_sha256": hashlib.sha256(here.read_bytes()).hexdigest(),
        "proof_sha256": hashlib.sha256(note.read_bytes()).hexdigest(),
        "checks": result,
        "independent_replay": False,
        "formal_verification": False,
        "external_review": False,
    }
    if args.output:
        args.output.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
