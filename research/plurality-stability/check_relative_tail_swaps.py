#!/usr/bin/env python3
"""Finite tensor and geometry controls for the analytic tail-swap proof.

Run with Python 3.12.3, standard library only, from the repository root:
    python3 research/plurality-stability/check_relative_tail_swaps.py

These are floating-point controls, not interval certificates or a proof for
all R. Product truncated-normal integrals assemble the moment tensors before
comparison with the scalar formulas in RELATIVE_HERMITE_TAIL_OBSTRUCTION.md.
"""

import itertools
import json
import math
import platform


ROOT2 = math.sqrt(2)
C = (0.25 + math.asin(1 / 3) / (2 * math.pi)) / math.sqrt(math.pi)
H = 1 / (2 * math.pi * math.sqrt(3))
SIGNS = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
# Express the original even sign vectors in the U,N,V orthonormal frame.
SITES = [(x, (y + z) / ROOT2, (y - z) / ROOT2) for x, y, z in SIGNS]
BASE_M = [[C * x for x in site] for site in SITES]
BASE_Q = [
    [[H * (site[j] * site[k] - (j == k)) for k in range(3)] for j in range(3)]
    for site in SITES
]


def phi(x):
    return math.exp(-x * x / 2) / math.sqrt(2 * math.pi)


def interval_moments(left, right):
    """Integrals of 1,x,x^2 against the standard Gaussian density."""
    assert left < right
    if right <= 0:
        mass = (math.erfc(-right / ROOT2) - math.erfc(-left / ROOT2)) / 2
    elif left >= 0:
        mass = (math.erfc(left / ROOT2) - math.erfc(right / ROOT2)) / 2
    else:
        mass = (math.erf(right / ROOT2) - math.erf(left / ROOT2)) / 2
    return mass, phi(left) - phi(right), mass + left * phi(left) - right * phi(right)


def box_tensors(intervals):
    moments = [interval_moments(*interval) for interval in intervals]
    mass = math.prod(row[0] for row in moments)
    first = [math.prod(moments[k][int(k == j)] for k in range(3)) for j in range(3)]
    second = [
        [math.prod(moments[k][int(k == j) + int(k == ell)] for k in range(3))
         for ell in range(3)]
        for j in range(3)
    ]
    q_tensor = [[second[j][k] - mass * (j == k) for k in range(3)] for j in range(3)]
    return mass, first, q_tensor


def dot(a, b):
    return math.fsum(x * y for x, y in zip(a, b))


def frobenius(a, b):
    return math.fsum(dot(row_a, row_b) for row_a, row_b in zip(a, b))


def tensor_changes(first, q_tensor):
    reflection = (1, -1, 1)
    delta_m = [first[j] * (1 - reflection[j]) for j in range(3)]
    delta_q = [[q_tensor[j][k] * (1 - reflection[j] * reflection[k])
                for k in range(3)] for j in range(3)]
    # Expand the increments, avoiding cancellation in W(T)+tiny_change-W(T).
    dw1 = math.fsum(
        2 * dot(BASE_M[i], [sign * x for x in delta_m]) + dot(delta_m, delta_m)
        for i, sign in ((0, 1), (1, -1))
    )
    dw2 = math.fsum(
        sign * frobenius(BASE_Q[i], delta_q) + frobenius(delta_q, delta_q) / 2
        for i, sign in ((0, 1), (1, -1))
    )
    return -dw1, dw2


def check_case(radius, positive_u=False):
    u_interval = (radius, radius + 1) if positive_u else (-radius - 1, -radius)
    intervals = (u_interval, (-4 * radius - 1, -4 * radius), (-1, 1))
    for corner in itertools.product(*intervals):
        scores = [dot(site, corner) for site in SITES]
        assert all(scores[1] > scores[k] for k in (0, 2, 3))
        reflected = (corner[0], -corner[1], corner[2])
        reflected_scores = [dot(site, reflected) for site in SITES]
        assert all(reflected_scores[0] > reflected_scores[k] for k in (1, 2, 3))

    mass, first, q_tensor = box_tensors(intervals)
    reflected_mass, reflected_first, reflected_q = box_tensors(
        (u_interval, (4 * radius, 4 * radius + 1), (-1, 1))
    )
    assert mass > 0 and mass == reflected_mass  # Exact symmetry of the implementation.
    assert first[1] == -reflected_first[1]
    assert math.isclose(q_tensor[0][1], -reflected_q[0][1], rel_tol=2e-14)
    a_value, b_value = -first[1], q_tensor[0][1]
    deficit, second_gain = tensor_changes(first, q_tensor)
    expected_deficit = 8 * ROOT2 * C * a_value - 8 * a_value**2
    expected_gain = 8 * ROOT2 * H * b_value + 8 * b_value**2
    assert math.isclose(deficit, expected_deficit, rel_tol=3e-14)
    assert math.isclose(second_gain, expected_gain, rel_tol=3e-14)
    assert 0 < a_value < C / ROOT2
    assert deficit > 0 and a_value > 4 * radius * mass
    assert radius < abs(b_value / a_value) < radius + 1
    ratio = second_gain / deficit
    if positive_u:
        assert second_gain < 0  # Detect reversing the mixed-moment sign.
    else:
        assert ratio > H / C * radius
    return {
        "R": radius,
        "positive_U_control": positive_u,
        "log10_box_mass": math.log10(mass),
        "delta": deficit,
        "W2_gain_over_delta": ratio,
        "proved_lower_bound_for_negative_U": H / C * radius if not positive_u else None,
    }


if __name__ == "__main__":
    assert math.isclose(sum(dot(m, m) for m in BASE_M), 12 * C**2, rel_tol=1e-14)
    assert math.isclose(sum(frobenius(q, q) for q in BASE_Q) / 2, 1 / math.pi**2,
                        rel_tol=1e-14)
    rows = [check_case(radius, positive_u) for radius in (1, 2, 4, 8)
            for positive_u in (False, True)]
    print(json.dumps({
        "python": platform.python_version(),
        "scope": "finite floating-point geometry and tensor controls; unbounded proof is written",
        "cases": rows,
        "outcome": "PASS",
    }, indent=2))
