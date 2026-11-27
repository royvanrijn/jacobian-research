#!/usr/bin/env python3
"""Controls for balanced variations and the non-extreme first-moment bound.

Run from the repository root with Python 3.12.3 and NumPy 1.26.4:
    python3 research/plurality-stability/check_balanced_variations.py

Finite-difference and quadrature diagnostics are not interval certificates.
The rational checks verify the displayed constant arithmetic, not the
general geometric arguments in the accompanying note.
"""

import json
from fractions import Fraction
from math import exp, pi, sqrt
from statistics import NormalDist

import numpy as np

import search_balanced_affine as affine


def nonextreme_bound():
    F = Fraction
    endpoint = F(2, 3)
    integral_upper = endpoint - endpoint**3 / 6 + endpoint**5 / 40
    quartile_margin = F(1, 4) - F(2, 5) * integral_upper
    exponential_lower = 1 + F(4, 9) + F(4, 9)**2 / 2
    root_margin = F(25, 24)**2 + F(3, 8) * F(25, 24) - F(729, 500)
    tetra_margin = F(7, 22) * F(3, 4) * (1 + F(7, 33))**2 - F(1, 3)
    assert integral_upper == F(754, 1215)
    assert quartile_margin == F(43, 24300) > 0
    assert exponential_lower == F(125, 81)
    assert root_margin == F(637, 36000) > 0
    assert tetra_margin == F(23, 1331) > 0
    # The decimals are illustrative; the strict comparison uses the fractions.
    quartile = NormalDist().inv_cdf(0.75)
    bound = 3 / (16 * pi) * (sqrt(1 + 64 * exp(-quartile**2)) - 1)
    return {
        "bound_decimal": bound,
        "quartile_decimal": quartile,
        "exact_quartile_margin": str(quartile_margin),
        "exact_exponential_lower": str(exponential_lower),
        "exact_root_margin": str(root_margin),
        "exact_tetrahedron_margin_above_one_third": str(tetra_margin),
    }


def square_lifts(rule):
    sites = np.array([[1, 1], [-1, 1], [-1, -1], [1, -1]]) / (2 * sqrt(2 * pi))
    direction = np.array([1, -1, 1, -1]) / 2
    rows = []
    for step in (0.01, 0.003, 0.001):
        lifted = np.column_stack((sites, step * direction))
        cells = affine.AffineCells(lifted, rule)
        masses, areas = cells.masses_and_facets(np.zeros(4))
        moments = affine.DIFFERENCE.T @ (areas[:, None] * cells.normal)
        value = 2 * np.sum(lifted * moments) - np.sum(lifted**2)
        quotient = float((value - 1 / pi) / step**2)
        error = float(np.max(np.abs(masses - 0.25)))
        assert error < 1e-12
        assert quotient > 0
        rows.append({"step": step, "quadratic_quotient": quotient, "balance_error": error})
    assert abs(rows[-1]["quadratic_quotient"] - 1) < 0.003
    return rows


def balanced_hessian(rule):
    sites = affine.TETRAHEDRON @ np.array(
        [[1.2, 0.1, -0.2], [0.1, 0.75, 0.2], [-0.2, 0.2, 1.1]]
    )
    direction = np.array(
        [[0.2, -0.7, 0.6], [0.5, 0.1, -0.8], [-0.3, 0.2, 0.1], [-0.4, 0.4, 0.1]]
    )
    base = affine.AffineCells(sites, rule)
    bias, _, _ = base.balance(np.zeros(4))
    _, areas = base.masses_and_facets(bias)
    laplacian = affine.DIFFERENCE.T @ ((areas / base.length)[:, None] * affine.DIFFERENCE)
    step = 2e-5
    fixed_masses, fixed_moments, balanced_moments = [], [], []
    for sign in (-1, 1):
        cells = affine.AffineCells(sites + sign * step * direction, rule)
        masses, areas = cells.masses_and_facets(bias)
        fixed_masses.append(masses)
        fixed_moments.append(affine.DIFFERENCE.T @ (areas[:, None] * cells.normal))
        balanced_moments.append(cells.balance(bias.copy())[1])
    mixed = (fixed_masses[1] - fixed_masses[0]) / (2 * step)
    fixed_second = np.sum((fixed_moments[1] - fixed_moments[0]) * direction) / (2 * step)
    correction = float(mixed @ np.linalg.solve(laplacian + affine.MEAN_PROJECTION, mixed))
    predicted = 2 * (fixed_second - correction - np.sum(direction**2))
    observed = 2 * (
        np.sum((balanced_moments[1] - balanced_moments[0]) * direction) / (2 * step)
        - np.sum(direction**2)
    )
    discrepancy = float(abs(predicted - observed))
    assert correction > 0.01  # This control must detect omission of the bias correction.
    assert discrepancy < 1e-7
    return {
        "step": step,
        "balancing_correction": correction,
        "hessian_via_schur_complement": float(predicted),
        "hessian_via_balanced_gradient": float(observed),
        "absolute_discrepancy": discrepancy,
    }


if __name__ == "__main__":
    rule = affine.rule(80)
    affine.check_facet_calculus(rule)
    print(json.dumps({
        "nonextreme_bound": nonextreme_bound(),
        "square_lifts": square_lifts(rule),
        "balanced_hessian": balanced_hessian(rule),
    }, indent=2))
