#!/usr/bin/env python3
"""Finite diagnostics for DG7--DG12; no global partition inequality is tested."""

from fractions import Fraction
import json
import math
import platform


def quadrant_formula(a, c):
    d = a - c * c
    return 1 / d - c * math.acos(c / math.sqrt(a)) / d**1.5


def simpson(f, start, stop, intervals=8192):
    assert intervals > 0 and intervals % 2 == 0
    step = (stop - start) / intervals
    total = f(start) + f(stop)
    total += 4 * math.fsum(f(start + j * step) for j in range(1, intervals, 2))
    total += 2 * math.fsum(f(start + j * step) for j in range(2, intervals, 2))
    return total * step / 3


def quadrant_integral(a, c):
    # Transform u = R cos(theta)/sqrt(a), v = R sin(theta), then
    # integrate R^3 exp(-q R^2/2) exactly. This does not use arccos.
    def integrand(theta):
        cs = math.cos(theta) * math.sin(theta)
        return 2 * cs / (a * (1 + 2 * c * cs / math.sqrt(a)) ** 2)

    return simpson(integrand, 0, math.pi / 2)


def relaxation_coefficient(rho, quadrant):
    sigma2 = 1 - rho * rho
    multiplier = sigma2 * (
        quadrant(1 + 2 * sigma2, rho / 2)
        - quadrant(1 + sigma2, rho / 2)
        + quadrant(1, rho / 2) / 4
    )
    radial_dirichlet = sigma2 / (3 * (4 - rho * rho))
    return multiplier - radial_dirichlet


def main():
    # With x = r^2/2 the radial law is Exp(1). b = exp(-x)-1/2.
    mean = Fraction(1, 2) - Fraction(1, 2)
    norm2 = Fraction(1, 3) - Fraction(1, 2) + Fraction(1, 4)
    # Integral |grad b|^2 d gamma_2 = integral 2x exp(-3x) dx.
    dirichlet = Fraction(2, 9)
    assert mean == 0 and norm2 == Fraction(1, 12)
    assert dirichlet == Fraction(2, 9)
    for rho in (Fraction(0), Fraction(1, 2), Fraction(9, 10)):
        pairing = 1 / (4 - rho * rho) - Fraction(1, 4)
        assert norm2 - pairing == (1 - rho * rho) / (3 * (4 - rho * rho))

    cases = [(1, 0), (1, 0.5), (1.5, 0.25), (2.5, 0.49)]
    discrepancies = [abs(quadrant_formula(a, c) - quadrant_integral(a, c))
                     for a, c in cases]
    assert max(discrepancies) < 2e-12
    assert abs(quadrant_formula(1, 0) - 1) < 1e-15

    rows = []
    for rho in (0.1, 0.5, 0.9, 0.95, 0.99, 0.999, 0.9999):
        exact_expression = relaxation_coefficient(rho, quadrant_formula)
        quadrature = relaxation_coefficient(rho, quadrant_integral)
        assert abs(exact_expression - quadrature) < 2e-12
        if rho >= 0.99:
            assert exact_expression > 0
        if rho <= 0.95:
            assert exact_expression < 0
        rows.append({"rho": rho, "E_rho": exact_expression,
                     "E_rho_over_one_minus_rho": exact_expression / (1 - rho)})

    limit = (4 - 2 * math.pi / math.sqrt(3)) / 9
    assert 0.04 < limit < 0.042
    error = abs(rows[-1]["E_rho_over_one_minus_rho"] - limit)
    assert error < 0.00013
    # Reject a missing factor two when converting sigma^2 to 1-rho.
    assert abs(rows[-1]["E_rho_over_one_minus_rho"] - limit / 2) > 0.02

    print(json.dumps({
        "status": "pass",
        "python": platform.python_version(),
        "arithmetic": "exact fractions and IEEE-754 binary64 diagnostics",
        "quadrature": "composite Simpson on [0,pi/2], 8192 subintervals",
        "quadrant_max_absolute_discrepancy": max(discrepancies),
        "radial_mean": str(mean), "radial_norm_squared": str(norm2),
        "radial_dirichlet_energy": str(dirichlet),
        "near_one_limit": limit, "finite_correlation_controls": rows,
        "scope": "Finite normalization diagnostics only; written DG1-DG12 supply analytic claims. No unrestricted Gaussian partition search or interval certification.",
    }, indent=2))


if __name__ == "__main__":
    main()
