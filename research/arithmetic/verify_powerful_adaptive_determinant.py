#!/usr/bin/env python3
"""Replay the rational parameter checks for the adaptive determinant."""

from fractions import Fraction
from typing import NamedTuple


class Parameters(NamedTuple):
    depth: int
    amplification: int
    max_v_degree: int
    t_degree_rate: Fraction
    q_lower: Fraction
    q_upper: Fraction


def audit(saving, t_height, residual_height, parameters):
    depth, amplification, max_v_degree, t_rate, q_lower, q_upper = parameters
    dimension_rate = (max_v_degree + 1) * t_rate
    precision_loss = dimension_rate * (depth + 4 * amplification * t_height - 2 * q_lower)
    height = amplification * (t_height * t_rate + max_v_degree * residual_height)
    divisibility = q_lower * dimension_rate / 2
    if not (precision_loss < 1 and height < divisibility
            and depth * dimension_rate < 1
            and 2 * amplification * residual_height < 1
            and q_upper / amplification <= 1 - saving):
        raise ValueError(f"Invalid parameters at saving={saving}: {parameters}")
    return precision_loss, divisibility - height


def high_q_search(saving, t_height, residual_height):
    """Find feasible high-q choices within the prime-power search window."""
    for depth in range(1, 7):
        amplification = depth + 1
        if Fraction(depth, amplification) > 1 - saving:
            continue
        q_lower = depth * (1 - saving)
        max_dimension = Fraction(1, depth + 4 * amplification * t_height - 2 * q_lower)
        for max_v_degree in range(1, 51):
            margin = q_lower / 2 - amplification * t_height / (max_v_degree + 1)
            if margin > 0 and max_dimension * margin > amplification * max_v_degree * residual_height:
                yield depth, amplification, max_v_degree


def prime_power_ceiling():
    """Maximize the two limiting inequalities when n is p**k."""
    # For fixed depth, the smallest possible order is depth + 1 and
    # increasing it decreases the upper bound. At depth >= 7, use
    # (1/2 - 1/(m+1))/m <= 1/12 to bound the answer by 1/300.
    # At degree >= 51, use 1/(2*m*(3*depth+4)) < 1/300.
    candidates = []
    for depth in range(1, 7):
        amplification = depth + 1
        for max_v_degree in range(1, 51):
            coefficient = Fraction(depth, 2 * amplification) - Fraction(1, max_v_degree + 1)
            saving = coefficient / (max_v_degree * (4 * amplification - depth))
            candidates.append((saving, depth, amplification, max_v_degree))
    return max(candidates)


def main():
    for saving, t_height, residual_height, high in [
        (Fraction(1, 320), Fraction(1), Fraction(1, 320),
         Parameters(2, 3, 5, Fraction(2, 125), Fraction(319, 160), Fraction(2))),
        (Fraction(1, 300), Fraction(1), Fraction(1, 300),
         Parameters(3, 4, 4, Fraction(48, 3125), Fraction(299, 100), Fraction(3))),
        (Fraction(1, 153), Fraction(306, 305), Fraction(1, 305),
         Parameters(3, 4, 4, Fraction(3817, 250000), Fraction(152, 51), Fraction(3))),
    ]:
        low = Parameters(1, 1, 7, Fraction(3, 100), Fraction(1, 2), 1 - saving)
        if saving in (Fraction(1, 320), Fraction(1, 300)):
            if not (t_height >= 1 and residual_height == saving):
                raise ValueError("The dyadic c-based height bounds do not fit.")
        else:
            if not (t_height * (1 - saving / 2) >= 1
                    and residual_height * (1 - saving / 2) >= saving / 2):
                raise ValueError("The adaptive target height bounds do not fit.")
        for case, parameters in [("low", low), ("high", high)]:
            loss, gap = audit(saving, t_height, residual_height, parameters)
            print(f"saving={saving}, {case}: precision loss <= {loss}, height gap >= {gap}")
        if saving != Fraction(1, 320):
            optimized = list(high_q_search(saving, t_height, residual_height))
            if optimized != [(3, 4, 4)]:
                raise ValueError(f"Unexpected high-q search result: {optimized}")
            print(f"Feasible high-q choices at {saving}: {optimized}")
    optimal = prime_power_ceiling()
    if optimal != (Fraction(7, 2080), 3, 4, 4):
        raise ValueError(f"Unexpected prime-power ceiling: {optimal}")
    print(f"Prime-power ceiling in this criterion: {optimal[0]} at {optimal[1:]}")


if __name__ == "__main__":
    main()
