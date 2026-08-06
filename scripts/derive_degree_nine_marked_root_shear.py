#!/usr/bin/env python3
"""Derive the degree-nine monic marked-root completing shear exactly."""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
from pathlib import Path

import sympy as sp

from derive_degree_seven_marked_root_shear import completing_shear, w


def shear_at(point: tuple[int, int]):
    sigma = sp.Rational(point[0])
    tau = sp.Rational(point[1])
    factor = (
        w**6
        + sigma * w**5
        + tau * w**4
        + (-sp.Rational(17, 2) - 5 * sigma - 4 * tau) * w
        + sp.Rational(13, 2)
        + 4 * sigma
        + 3 * tau
    )
    return point, completing_shear(factor)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    interpolation_points = ((0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2))
    test_point = (2, 3)
    points = (*interpolation_points, test_point)
    if args.jobs == 1:
        pairs = [shear_at(point) for point in points]
    else:
        with mp.get_context("spawn").Pool(args.jobs) as pool:
            pairs = pool.map(shear_at, points)
    records = dict(pairs)

    sigma, tau = sp.symbols("sigma tau")
    coefficients = sp.symbols("c0:6")
    quadratic = (
        coefficients[0]
        + coefficients[1] * sigma
        + coefficients[2] * tau
        + coefficients[3] * sigma**2
        + coefficients[4] * sigma * tau
        + coefficients[5] * tau**2
    )
    solutions = sp.solve(
        [
            quadratic.subs({sigma: point[0], tau: point[1]})
            - records[point]
            for point in interpolation_points
        ],
        coefficients,
        dict=True,
    )
    assert len(solutions) == 1
    candidate = sp.factor(quadratic.subs(solutions[0]))
    assert candidate.subs({sigma: 2, tau: 3}) == records[test_point]

    certificate = {
        "scope": "exact degree-nine classical rank-two shear interpolation",
        "family_factor": (
            "w^6+sigma*w^5+tau*w^4"
            "+(-17/2-5*sigma-4*tau)*w+13/2+4*sigma+3*tau"
        ),
        "interpolation_points": [
            {"sigma": x, "tau": y, "shear": str(records[(x, y)])}
            for x, y in interpolation_points
        ],
        "independent_test_point": {
            "sigma": test_point[0],
            "tau": test_point[1],
            "shear": str(records[test_point]),
        },
        "quadratic_shear": str(candidate),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print("PASS: six exact homotopy samples determine the quadratic shear")
    print("PASS: the seventh exact sample is predicted correctly")
    print(f"shear = {candidate}")


if __name__ == "__main__":
    main()
