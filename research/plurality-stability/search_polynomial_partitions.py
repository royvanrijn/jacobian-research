#!/usr/bin/env python3
"""Numerically search balanced four-cell Gaussian partitions in three dimensions."""

import json
import math

import numpy as np


CORRELATIONS = (0.1, 0.35, 0.65, 0.9)
SEED = 468
SCREEN_PAIRS = 2**14
REFINEMENT_PAIRS = 2**15
HOLDOUT_PAIRS = 2**20
RANDOM_DIRECTIONS = 80
STEPS = (0.25, 0.0625, 0.015625, 0.00390625)
BASELINE = np.tile([1.0, 0.0, 0.0, 0.0, 0.0], (3, 1))
ORTHANTS = BASELINE.copy()
ORTHANTS[2, 0] = 0.0
QUADRATIC = np.tile([0.0, 1.0, 0.0, 0.0, 0.0], (3, 1))


def polynomial_features(points):
    features = []
    for axis in range(3):
        x, y, z = np.roll(points, -axis, axis=1).T
        features.append([x, y * z, x * (x * x - 3), x * (y * y - 1), x * (z * z - 1)])
    return np.asarray(features)


def labels(features, coefficients):
    first, second, third = np.einsum("jkn,jk->jn", features, coefficients)
    # Rows encode the four sign vectors (+++), (+--), (-+-), (--+).
    return np.argmax(
        [first + second + third, first - second - third,
         -first + second - third, -first - second + third],
        axis=0,
    )


def agreement(features, coefficients):
    left, right = features
    return labels(left, coefficients) == labels(right, coefficients)


def gaussian_pair(sample_count, correlation, rng):
    x = rng.standard_normal((sample_count, 3))
    z = rng.standard_normal(x.shape)
    y = correlation * x + math.sqrt((1 - correlation) * (1 + correlation)) * z
    return polynomial_features(x), polynomial_features(y)


def distinct_scores(coefficients):
    # Two rows with nonzero polynomials suffice, since each pair of sign
    # vectors differs in two coordinates with disjoint monomial supports.
    return np.count_nonzero(np.any(coefficients != 0, axis=1)) >= 2


def grid(rng):
    yield BASELINE
    yield ORTHANTS
    yield QUADRATIC
    for axis in range(1, 5):
        for size in (-2.0, -1.0, -0.5, -0.25, -0.0625, 0.0625, 0.25, 0.5, 1.0, 2.0):
            coefficients = BASELINE.copy()
            coefficients[:, axis] = size
            yield coefficients
    for size in (0.0625, 0.25, 1.0, 2.0):
        for _ in range(RANDOM_DIRECTIONS):
            direction = rng.integers(-4, 5, BASELINE.shape) / 4
            yield BASELINE + size * direction


def refine(start, features):
    current = start.copy()
    score = np.count_nonzero(agreement(features, current))
    for step in STEPS:
        for row, column in np.ndindex(current.shape):
            best = current
            for direction in (-1, 1):
                candidate = current.copy()
                candidate[row, column] += direction * step
                if distinct_scores(candidate):
                    candidate_score = np.count_nonzero(agreement(features, candidate))
                    if candidate_score > score:
                        best, score = candidate, candidate_score
            current = best
    return current, int(score)


def summarize(scores, baseline):
    sample_count = len(scores)
    wins = int(np.count_nonzero(scores & ~baseline))
    losses = int(np.count_nonzero(baseline & ~scores))
    delta = (wins - losses) / sample_count
    standard_error = math.sqrt(max(0.0, (wins + losses) / sample_count - delta**2)
                               / (sample_count - 1))
    return dict(agreement=int(np.count_nonzero(scores)), wins=wins, losses=losses,
                delta=delta, paired_standard_error=standard_error)


def search(correlation, seeds):
    grid_seed, screening_seed, refinement_seed, holdout_seed = seeds.spawn(4)
    screening = gaussian_pair(SCREEN_PAIRS, correlation, np.random.default_rng(screening_seed))
    screened = sorted(
        ((int(np.count_nonzero(agreement(screening, candidate))), candidate)
         for candidate in grid(np.random.default_rng(grid_seed)) if distinct_scores(candidate)),
        key=lambda item: item[0],
        reverse=True,
    )
    refinement = gaussian_pair(REFINEMENT_PAIRS, correlation, np.random.default_rng(refinement_seed))
    optimized = []
    seen = set()
    for start in (BASELINE, ORTHANTS, QUADRATIC, *(item[1] for item in screened[:3])):
        candidate, score = refine(start, refinement)
        key = tuple(candidate.flat)
        if key not in seen:
            seen.add(key)
            optimized.append(("refined", candidate, score))

    holdout = gaussian_pair(HOLDOUT_PAIRS, correlation, np.random.default_rng(holdout_seed))
    baseline = agreement(holdout, BASELINE)
    orthants = agreement(holdout, ORTHANTS)
    quadratic = agreement(holdout, QUADRATIC)
    measured = [
        dict(kind=kind, coefficients=coefficients.tolist(), refinement_agreement=score,
             holdout=summarize(agreement(holdout, coefficients), baseline))
        for kind, coefficients, score in sorted(optimized, key=lambda item: -item[2])
    ]
    return dict(
        rho=correlation,
        screening_candidates=len(screened),
        refinement_baseline_agreement=int(np.count_nonzero(agreement(refinement, BASELINE))),
        holdout_baseline_agreement=int(np.count_nonzero(baseline)),
        holdout_orthants_agreement=int(np.count_nonzero(orthants)),
        holdout_quadratic_agreement=int(np.count_nonzero(quadratic)),
        candidates=measured,
    )


def main():
    seeds = np.random.SeedSequence(SEED).spawn(len(CORRELATIONS))
    result = dict(
        numpy=np.__version__,
        seed=SEED,
        screening_pairs=SCREEN_PAIRS,
        refinement_pairs=REFINEMENT_PAIRS,
        holdout_pairs=HOLDOUT_PAIRS,
        correlations=[],
    )
    for correlation, seed in zip(CORRELATIONS, seeds):
        result["correlations"].append(search(correlation, seed))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
