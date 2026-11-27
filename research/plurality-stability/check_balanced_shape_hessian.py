#!/usr/bin/env python3
"""Diagnostics for the five balanced shape directions and exact trace identity.

Python 3.12.3, NumPy 1.26.4; run from the repository root:
    OPENBLAS_NUM_THREADS=1 python3 research/plurality-stability/check_balanced_shape_hessian.py

The face formulas and quotient identities have written analytic proofs.
Quadrature and finite differences here are controls, not interval certificates.
The explicit local-radius arithmetic is checked separately with fractions.
"""

import json
import platform
from fractions import Fraction as F
from math import atan, exp, pi, sin, sqrt

import numpy as np

import search_balanced_affine as affine


LABEL_BASIS = sqrt(3) * affine.TETRAHEDRON


def density(x):
    return exp(-x * x / 2) / sqrt(2 * pi)


def face_moments(cells, bias):
    """Ambient unnormalized Gaussian area, first and second moments per face."""
    _, area = cells.masses_and_facets(bias)
    first = np.zeros((6, 3))
    second = np.zeros((6, 3, 3))
    for e, (i, j) in enumerate(affine.PAIRS):
        n, d = cells.normal[e], cells.length[e]
        h = (bias[j] - bias[i]) / d
        rivals = cells.sites[i] - cells.sites[affine.OTHERS[e]]
        normal_parts = rivals @ n
        tangent = rivals - normal_parts[:, None] * n
        lengths = np.linalg.norm(tangent, axis=1)
        u, v = tangent / lengths[:, None]
        thresholds = (bias[affine.OTHERS[e]] - bias[i] - h * normal_parts) / lengths
        a, b = thresholds
        rho = float(u @ v)
        sigma = sqrt(1 - rho * rho)
        probability = area[e] / density(h)
        f = density(a) * float(affine.normal_cdf(np.array((rho * a - b) / sigma)))
        g = density(b) * float(affine.normal_cdf(np.array((rho * b - a) / sigma)))
        tangent_first = f * u + g * v
        eta = density(a) * density((b - rho * a) / sigma)
        tangent_second = probability * (np.eye(3) - np.outer(n, n))
        tangent_second += a * f * np.outer(u, u) + b * g * np.outer(v, v)
        tangent_second += eta / sigma * (
            np.outer(u, v) + np.outer(v, u) - rho * (np.outer(u, u) + np.outer(v, v))
        )
        first[e] = density(h) * (h * n * probability + tangent_first)
        second[e] = density(h) * (
            h * h * probability * np.outer(n, n) + tangent_second
            + h * (np.outer(n, tangent_first) + np.outer(tangent_first, n))
        )
    return area, first, second


def hessian_data(sites, rule):
    cells = affine.AffineCells(sites, rule)
    bias, moments, masses = cells.balance(np.zeros(4))
    area, first, second = face_moments(cells, bias)
    laplacian = affine.DIFFERENCE.T @ ((area / cells.length)[:, None] * affine.DIFFERENCE)
    q = np.zeros((9, 9))
    mixed = np.zeros((4, 9))
    cell_second = np.zeros((4, 3, 3))
    for e, difference in enumerate(affine.DIFFERENCE):
        label = LABEL_BASIS.T @ difference
        q += np.kron(np.outer(label, label), second[e] / cells.length[e])
        mixed += np.outer(difference, np.kron(label, first[e] / cells.length[e]))
        edge = difference @ sites
        cell_second += difference[:, None, None] * np.outer(edge, first[e])[None, :, :] / cells.length[e]
    correction = mixed.T @ np.linalg.solve(laplacian + affine.MEAN_PROJECTION, mixed)
    hessian = 2 * (q - correction - np.eye(9))
    assert np.max(np.abs(cell_second - cell_second.transpose(0, 2, 1))) < 2e-12
    assert np.max(np.abs(masses - 0.25)) < 6e-12
    assert np.max(np.abs(hessian - hessian.T)) < 2e-12
    return dict(cells=cells, bias=bias, moments=moments, area=area, first=first,
                second=second, laplacian=laplacian, mixed=mixed,
                cell_second=cell_second, q=q, correction=correction, hessian=hessian)


def shape_frame(sites):
    """Five orthonormal site directions orthogonal to rotation and scale."""
    values, axes = np.linalg.eigh(sites.T @ sites)
    assert values[0] > 0
    rotated = sites @ axes
    _, _, vectors = np.linalg.svd(np.sqrt(values).reshape(1, 3))
    directions = [rotated @ np.diag(d / np.sqrt(values)) @ axes.T for d in vectors[1:]]
    rotations = []
    for a, b in ((0, 1), (0, 2), (1, 2)):
        t = np.zeros((3, 3))
        t[a, b] = sqrt(values[b] / (values[a] * (values[a] + values[b])))
        t[b, a] = sqrt(values[a] / (values[b] * (values[a] + values[b])))
        directions.append(rotated @ t @ axes.T)
        skew = np.zeros((3, 3))
        skew[a, b], skew[b, a] = 1, -1
        rotations.append(rotated @ skew @ axes.T / sqrt(values[a] + values[b]))
    radial = sites / np.linalg.norm(sites)
    all_directions = directions + rotations + [radial]
    coordinates = np.array([(LABEL_BASIS.T @ v).reshape(9) for v in all_directions]).T
    assert np.max(np.abs(coordinates.T @ coordinates - np.eye(9))) < 3e-14
    return coordinates[:, :5], coordinates[:, 5:8], coordinates[:, 8]


def radius_control(sites, data, tetrahedral=False):
    """Check the radius metric, full rank, face restriction and tetrahedral sign.

    The asymmetric input is exactly rebalanced but is not stationary. Use its
    actual Laplacian inverse, rather than the stationary identity L^+=P.
    """
    g = np.sum(sites**2, axis=1)
    gram = sites @ sites.T
    metric = np.diag(g) - gram / 4
    assert np.linalg.eigvalsh(metric)[0] > 0
    assert np.max(np.abs(metric @ np.ones(4) - g)) < 2e-14
    _, _, vectors = np.linalg.svd(g.reshape(1, 4))
    frame = vectors[1:].T
    restricted_metric = frame.T @ metric @ frame
    values, axes = np.linalg.eigh(restricted_metric)
    frame = frame @ axes @ np.diag(1 / np.sqrt(values))
    directions = []
    for a in frame.T:
        v = a[:, None] * sites
        v -= v.mean(axis=0)
        assert np.max(np.abs(v.sum(axis=0))) < 2e-14
        assert abs(np.sum(sites * v)) < 2e-14
        assert np.max(np.abs(sites.T @ v - v.T @ sites)) < 2e-14
        assert np.max(np.abs(2 * np.sum(sites * v, axis=1) - 2 * metric @ a)) < 2e-14
        directions.append((LABEL_BASIS.T @ v).reshape(9))
    coordinates = np.array(directions).T
    assert np.max(np.abs(coordinates.T @ coordinates - np.eye(3))) < 2e-14
    assert np.linalg.matrix_rank(metric @ frame) == 3
    q = np.zeros((4, 4))
    mixed = np.zeros((4, 4))
    for e, (i, j) in enumerate(affine.PAIRS):
        edge = np.zeros((4, 3))
        edge[i], edge[j] = sites[i], -sites[j]
        q += edge @ data['second'][e] @ edge.T / data['cells'].length[e]
        mixed += np.outer(affine.DIFFERENCE[e], edge @ data['first'][e]) / data['cells'].length[e]
    correction = mixed.T @ np.linalg.solve(data['laplacian'] + affine.MEAN_PROJECTION, mixed)
    h_radius = 2 * (q - correction - metric)
    direct = frame.T @ h_radius @ frame
    compressed = coordinates.T @ data['hessian'] @ coordinates
    discrepancy = float(np.max(np.abs(direct - compressed)))
    assert discrepancy < 3e-13
    eigenvalues = np.linalg.eigvalsh(direct)
    if tetrahedral:
        theta = atan(sqrt(2))
        expected = -1 + sin(2 * theta) / (2 * theta) - pi / (3 * theta**2)
        assert np.max(np.abs(eigenvalues - expected)) < 3e-12
    return dict(restriction_dimension=3, direct_face_restriction_error=discrepancy,
                eigenvalues=eigenvalues.tolist(), stationary_input=tetrahedral,
                scope='metric and compression controls; no unequal-radius stationary sign certificate')


def finite_difference_control(sites, rule):
    data = hessian_data(sites, rule)
    observed = np.zeros((9, 9))
    step = 2e-5
    for column in range(9):
        direction = LABEL_BASIS @ np.eye(9)[column].reshape(3, 3)
        gradients = []
        for sign in (-1, 1):
            changed = sites + sign * step * direction
            _, moments, _ = affine.AffineCells(changed, rule).balance(data["bias"].copy())
            gradients.append((LABEL_BASIS.T @ (2 * (moments - changed))).reshape(9))
        observed[:, column] = (gradients[1] - gradients[0]) / (2 * step)
    error = float(np.max(np.abs(observed - data["hessian"])))
    assert error < 2e-7
    # The omitted-correction matrix must be detectably wrong at this input.
    omitted_error = float(np.linalg.norm(2 * data["correction"]))
    assert omitted_error > 0.01
    shape_frame(sites)  # Also checks the anisotropic horizontal metric.
    c = sites.T @ sites
    w = np.trace(c)
    whitening = np.linalg.inv(c) - 3 / w * np.eye(3)
    direction = sites @ whitening
    coordinates = (LABEL_BASIS.T @ direction).reshape(9)
    mixed = data["mixed"] @ coordinates
    via_cells = np.array([np.trace(d @ whitening) for d in data["cell_second"]])
    assert np.max(np.abs(mixed - via_cells)) < 3e-13
    anisotropy_norm = float(np.sum(direction**2))
    assert abs(anisotropy_norm - (np.trace(np.linalg.inv(c)) - 9 / w)) < 3e-13
    assert abs(np.sum(direction * sites)) < 3e-14
    assert np.max(np.abs(c @ whitening - whitening @ c)) < 3e-14
    return dict(step=step, maximum_entry_discrepancy=error,
                omitted_correction_matrix_norm=omitted_error,
                anisotropy_direction_squared_norm=anisotropy_norm,
                cell_second_moment_mixed_error=float(np.max(np.abs(mixed - via_cells))),
                radius_restriction=radius_control(sites, data))


def tetrahedral_control(rule):
    sites = sqrt(affine.TETRAHEDRON_MASS) * affine.TETRAHEDRON
    data = hessian_data(sites, rule)
    shape, rotation, radial = shape_frame(sites)
    h = data["hessian"]
    assert np.linalg.norm(data["moments"] - sites) < 2e-13
    assert np.max(np.abs(h @ rotation)) < 2e-12
    assert np.max(np.abs(h @ radial + 2 * radial)) < 2e-12
    theta = atan(sqrt(2))
    t = sin(2 * theta) / (2 * theta)
    off = -1 + t - pi / (3 * theta * theta)
    diagonal = -0.5 - 1.5 * t
    eigenvalues = np.linalg.eigvalsh(shape.T @ h @ shape)
    assert np.max(np.abs(eigenvalues - np.sort([off] * 3 + [diagonal] * 2))) < 3e-12
    c = sites.T @ sites
    balancing_trace = sum(np.trace(d @ np.linalg.solve(c, d)) for d in data["cell_second"])
    assert abs(balancing_trace - pi / (2 * theta * theta)) < 2e-12
    tau = np.trace(data["second"], axis1=1, axis2=2) / data["area"]
    assert np.max(np.abs(tau - 2)) < 2e-12
    formula = float(tau.sum() - 2 * balancing_trace - 16)
    expected = -4 - pi / (theta * theta)
    assert abs(formula - expected) < 3e-12
    assert abs(np.trace(shape.T @ h @ shape) - formula) < 3e-12
    return dict(eigenvalues=eigenvalues.tolist(), shape_trace=formula,
                balancing_trace=float(balancing_trace), stationarity_residual=float(np.linalg.norm(data["moments"] - sites)),
                radius_restriction=radius_control(sites, data, tetrahedral=True))


def local_radius_arithmetic():
    # Density derivatives: E[R^(2j+1)] / E[R] = 1,4,24,192 in dimension 3.
    radial_multipliers = [1]
    for j in range(1, 4):
        radial_multipliers.append(radial_multipliers[-1] * (2 * j + 2))
    assert radial_multipliers == [1, 4, 24, 192]
    assert sum(a * b for a, b in zip([60, 75, 18, 1], radial_multipliers)) == 984
    assert (1 + F(1, 10)) / (1 - F(1, 10)) ** 3 < 2
    radius = F(1, 100000)
    margin = F(1, 54) - 656 * radius - F(1, 100)
    assert margin > 0
    return dict(polar_frobenius_radius=str(radius), value_gap_coefficient="1/100",
                exact_taylor_margin=str(margin))


def main():
    rule = affine.rule(80)
    sites = affine.TETRAHEDRON @ np.array([[1.2, 0.1, -0.2], [0.1, 0.75, 0.2], [-0.2, 0.2, 1.1]])
    result = dict(status="PASS", python=platform.python_version(), numpy=np.__version__,
                  scope="analytic identities with quadrature/finite-difference controls, not a global or interval certificate",
                  tetrahedron=tetrahedral_control(rule),
                  asymmetric_hessian=finite_difference_control(sites, rule),
                  local_radius=local_radius_arithmetic(),
                  assurance="internal controls; no independent review or formal verification")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
