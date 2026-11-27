#!/usr/bin/env python3
"""Exact scalar audit for the stationary Gaussian boundary theorem (QB).

Run from the repository root with Python 3.12.3; no third-party packages:
    python3 research/plurality-stability/check_stationary_boundary_margin.py

This verifies constant arithmetic, not the geometric proof, Gaussian balance,
or a cover of the remaining shape domain. All Gaussian estimates used by the
theorem have analytic proofs in BALANCED_AFFINE_CENTROID_OPTIMIZATION.md.
"""

import json
from fractions import Fraction as F
from math import factorial


def audit(eta=F(1, 10**11), hull_margin=F(1, 25), radius=4):
    # Rational bounds imported from the note's analytic non-extreme argument.
    moment_square_upper = F(324, 3125)
    three_score_constant_upper = F(9, 25)
    high_value_lower = F(7, 20)
    tetrahedron_lower = F(1, 3) + F(23, 1331)
    value_upper = 4 * moment_square_upper
    singular_cap = F(1, 400000)

    assert tetrahedron_lower > high_value_lower
    assert value_upper < F(1, 2)
    reduced_score_square_upper = three_score_constant_upper * (
        value_upper - high_value_lower / 3
    )
    value_square_lower = (high_value_lower - hull_margin / 2)**2
    robust_margin = value_square_lower - reduced_score_square_upper
    assert high_value_lower > hull_margin / 2
    assert robust_margin > 0

    # lambda <= eta*W < eta/2 < singular_cap**2.
    assert eta / 2 < singular_cap**2
    projected_hull_margin = F(1, 50)
    assert hull_margin - 2 * singular_cap > projected_hull_margin
    dependence_lower = projected_hull_margin / 4
    assert dependence_lower == F(1, 200)

    # The tangent strip bound uses pi > 3; the tail uses phi(0) < 1/2.
    strip_upper = radius * singular_cap / (
        6 * dependence_lower * hull_margin
    )
    exponential_lower = sum(
        (F(radius * radius, 2)**j / factorial(j) for j in range(6)), F(0)
    )
    assert exponential_lower > 500
    tail_upper = F(1, 1000)
    stationary_area_lower = hull_margin / 4
    contradiction_margin = stationary_area_lower - strip_upper - tail_upper
    assert contradiction_margin > 0
    return {
        "eta": str(eta),
        "robust_hull_margin": str(hull_margin),
        "robust_value_squared_margin": str(robust_margin),
        "smallest_singular_value_upper": str(singular_cap),
        "strip_area_upper": str(strip_upper),
        "exponential_series_lower": str(exponential_lower),
        "plane_tail_area_upper": str(tail_upper),
        "stationary_area_lower": str(stationary_area_lower),
        "area_contradiction_margin": str(contradiction_margin),
    }


if __name__ == "__main__":
    certificate = audit()
    # Detect changes that would invalidate distinct parts of the certificate.
    for invalid in ({"eta": F(1, 10**9)},
                    {"hull_margin": F(1, 10)},
                    {"radius": 1}):
        try:
            audit(**invalid)
        except AssertionError:
            continue
        raise AssertionError(f"accepted invalid constants: {invalid}")
    print(json.dumps({"status": "PASS", "constants": certificate,
                      "invalid_constant_controls": 3,
                      "scope": "exact scalar arithmetic; geometric proof is separate"},
                     indent=2))
