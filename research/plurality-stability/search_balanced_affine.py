#!/usr/bin/env python3
"""Explore degree-one Gaussian stability via balanced affine partitions."""

import json
from itertools import combinations
from math import asin, atan, erfc, pi, sin, sqrt

import numpy as np
from numpy.polynomial.legendre import leggauss


PAIRS = np.array(list(combinations(range(4), 2)))
LEFT, RIGHT = PAIRS.T
OTHERS = np.array([[k for k in range(4) if k not in pair] for pair in PAIRS])
DIFFERENCE = np.eye(4)[LEFT] - np.eye(4)[RIGHT]
MEAN_PROJECTION = np.full((4, 4), 0.25)
TETRAHEDRON = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / sqrt(12)
TETRAHEDRON_MASS = 12 / pi * (1 / 4 + asin(1 / 3) / (2 * pi)) ** 2


def normal_cdf(x):
    return np.fromiter((erfc(-a / sqrt(2)) / 2 for a in x.flat), float, x.size).reshape(x.shape)


def centered_unit(sites):
    sites = sites - sites.mean(axis=0)
    return sites / np.linalg.norm(sites)


class AffineCells:
    def __init__(self, sites, rule):
        self.sites = sites
        self.nodes, self.weights = rule
        differences = sites[:, None] - sites[None, :]
        self.length = np.linalg.norm(differences[LEFT, RIGHT], axis=1)
        self.normal = differences[LEFT, RIGHT] / self.length[:, None]
        rivals = differences[LEFT[:, None], OTHERS]
        self.normal_rival = np.sum(rivals * self.normal[:, None, :], axis=2)
        tangent = rivals - self.normal_rival[:, :, None] * self.normal[:, None, :]
        self.tangent_length = np.linalg.norm(tangent, axis=2)
        tangent_unit = tangent / self.tangent_length[:, :, None]
        correlation = np.sum(tangent_unit[:, 0] * tangent_unit[:, 1], axis=1)
        angle = np.arcsin(np.clip(correlation, -1, 1))
        self.theta = angle[:, None] * self.nodes
        self.coefficient = angle / (2 * pi)

        self.conical_mass = np.empty(4)
        for i in range(4):
            rows = np.delete(differences[i], i, axis=0)
            rows = rows / np.linalg.norm(rows, axis=1)[:, None]
            correlations = (rows @ rows.T)[np.triu_indices(3, 1)]
            self.conical_mass[i] = 1 / 8 + np.arcsin(np.clip(correlations, -1, 1)).sum() / (4 * pi)

    def masses_and_facets(self, bias):
        offset = (bias[RIGHT] - bias[LEFT]) / self.length
        thresholds = (bias[LEFT, None] - bias[OTHERS] + offset[:, None] * self.normal_rival)
        thresholds = thresholds / self.tangent_length
        h, k = thresholds.T
        # Plackett's bivariate normal integral with correlation = sin(theta).
        exponent = -0.5 * ((h[:, None] - np.sin(self.theta) * k[:, None]) / np.cos(self.theta)) ** 2
        exponent -= 0.5 * k[:, None] ** 2
        # The final node evaluates the facet of the requested affine cells.
        # Earlier nodes integrate along the line from zero to the given bias.
        t = np.append(self.nodes, 1.0)
        probability = np.prod(normal_cdf(t[:, None, None] * thresholds), axis=2)
        probability += self.coefficient * (np.exp(t[:, None, None] ** 2 * exponent) @ self.weights)
        facets = probability * np.exp(-0.5 * (t[:, None] * offset) ** 2) / sqrt(2 * pi)
        displacement = (self.weights @ facets[:-1]) * (bias[LEFT] - bias[RIGHT]) / self.length
        return self.conical_mass + DIFFERENCE.T @ displacement, facets[-1]

    def balance(self, bias):
        for _ in range(24):
            mass, facets = self.masses_and_facets(bias)
            residual = mass - 0.25
            if np.max(np.abs(residual)) < 5e-12:
                centroids = DIFFERENCE.T @ (facets[:, None] * self.normal)
                return bias, centroids, mass
            laplacian = DIFFERENCE.T @ ((facets / self.length)[:, None] * DIFFERENCE)
            step = np.linalg.solve(laplacian + MEAN_PROJECTION, residual)
            for j in range(12):
                trial = bias - np.ldexp(step, -j)
                trial -= trial.mean()
                new_mass, _ = self.masses_and_facets(trial)
                if np.linalg.norm(new_mass - 0.25) < np.linalg.norm(residual):
                    bias = trial
                    break
            else:
                raise ArithmeticError("bias Newton step did not reduce estimated imbalance")
        raise ArithmeticError("bias solve exhausted its iteration limit")


def rule(degree):
    nodes, weights = leggauss(degree)
    return (nodes + 1) / 2, weights / 2


def follow(sites, coarse, fine):
    bias = np.zeros(4)
    for iteration in range(100):
        bias, moments, _ = AffineCells(sites, coarse).balance(bias)
        if iteration == 0:
            initial_objective = float(np.sum(moments**2))
            initial_centroid_norm_range = float(np.ptp(np.linalg.norm(moments, axis=1)))
        successor = centered_unit(moments)
        if np.linalg.norm(successor - sites) < 1e-7:
            break
        sites = successor
    bias, moments, masses = AffineCells(sites, fine).balance(bias)
    moment_mass = float(np.sum(moments**2))
    return {
        "initial_objective": initial_objective,
        "initial_centroid_norm_range": initial_centroid_norm_range,
        "objective": moment_mass,
        "gap": moment_mass - TETRAHEDRON_MASS,
        "stationarity_residual": float(np.linalg.norm(centered_unit(moments) - sites)),
        "isotropic_gram_error": float(np.linalg.norm(sites.T @ sites - np.eye(3) / 3)),
        "centroid_norms": np.linalg.norm(moments, axis=1).tolist(),
        "masses": masses.tolist(),
        "iterations": iteration + 1,
        "sites": sites.tolist(),
        "bias": bias.tolist(),
    }


def check_facet_calculus(quadrature_rule):
    sites = sqrt(TETRAHEDRON_MASS) * TETRAHEDRON
    theta = atan(sqrt(2))
    angular_term = sin(2 * theta) / (2 * theta)
    direction = np.array([[1.0, 0.8, -0.2], [0.8, -0.3, 0.4], [-0.2, 0.4, 0.6]])
    perturbation = sites @ direction
    step = 1e-5
    expected = TETRAHEDRON_MASS / 3 * (
        (1 + angular_term) / 2 * np.sum((direction - np.diag(direction.diagonal())) ** 2)
        + 3 * (1 - angular_term) / 4 * np.sum((direction.diagonal() - np.trace(direction) / 3) ** 2)
    )
    directional_moments = []
    for sign in (-1, 1):
        cells = AffineCells(sites + sign * step * perturbation, quadrature_rule)
        _, areas = cells.masses_and_facets(np.zeros(4))
        moments = DIFFERENCE.T @ (areas[:, None] * cells.normal)
        directional_moments.append(np.sum(moments * perturbation))
    assert abs((directional_moments[1] - directional_moments[0]) / (2 * step) - expected) < 2e-9

    cells = AffineCells(sites + 0.2 * perturbation, quadrature_rule)
    bias = np.array([-0.09, 0.15, 0.08, -0.14])
    direction = np.array([1.0, 0.3, -0.7, -0.6])
    _, areas = cells.masses_and_facets(bias)
    matrix = DIFFERENCE.T @ ((areas / cells.length)[:, None] * DIFFERENCE)
    lower, _ = cells.masses_and_facets(bias - step * direction)
    upper, _ = cells.masses_and_facets(bias + step * direction)
    assert np.max(np.abs((upper - lower) / (2 * step) - matrix @ direction)) < 2e-9


def initial_sites(rng):
    yield "simplex", "regular", TETRAHEDRON
    for j in range(32):
        yield f"random-{j}", "random", centered_unit(rng.normal(size=(4, 3)))
    for compression in (0.07, 0.2, 0.5):
        for j in range(8):
            sites = rng.normal(size=(4, 3)) * [1, 1, compression]
            yield f"flattened-{compression}-{j}", "flattened", centered_unit(sites)
    propeller = np.array([[1, 0, 0], [-0.5, sqrt(3) / 2, 0], [-0.5, -sqrt(3) / 2, 0], [0, 0, 1]])
    for height in (0.04, 0.15, 0.5, 1, 3):
        for j in range(4):
            sites = propeller * [1, 1, height] + 0.03 * rng.normal(size=(4, 3))
            yield f"axial-{height}-{j}", "axial", centered_unit(sites)


def main():
    coarse, fine = rule(40), rule(80)
    check_facet_calculus(fine)
    control = follow(TETRAHEDRON, coarse, fine)
    assert abs(control["gap"]) < 1e-13, control
    rows = {}
    for name, family, sites in initial_sites(np.random.default_rng(468)):
        try:
            result = follow(sites, coarse, fine)
        except (ArithmeticError, np.linalg.LinAlgError) as error:
            result = {"error": str(error)}
        rows[name] = {"family": family, **result}
    print(json.dumps({"tetrahedron": TETRAHEDRON_MASS, "runs": rows}, indent=2))


if __name__ == "__main__":
    main()
