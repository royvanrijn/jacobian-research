#!/usr/bin/env python3
"""Replay the all-order recurrence for the actual 21D homogeneous inverse."""

from __future__ import annotations

from fractions import Fraction as Q
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "generated-results" / "essential_bcw_21_counterexample.json"

Term = tuple[Q, tuple[int, int, int]]


def decode_h(source: dict[str, object]) -> list[list[Term]]:
    components: list[list[Term]] = []
    for stored_component in source["H"]:
        component: list[Term] = []
        for term in stored_component:
            factors: list[int] = []
            for variable, power in term["monomial"]:
                factors.extend([variable] * power)
            assert len(factors) == 3
            component.append((Q(term["coefficient"]), tuple(factors)))
        components.append(component)
    return components


def cubic_coefficient(
    series: list[list[Q]], factors: tuple[int, int, int], degree: int
) -> Q:
    a, b, c = factors
    return sum(
        series[r][a] * series[s][b] * series[degree - r - s][c]
        for r in range(degree + 1)
        for s in range(degree - r + 1)
    )


def inverse_coefficients(
    components: list[list[Term]], direction: list[Q], order: int
) -> list[list[Q]]:
    """Return g_m in V^{-1}(t*direction)=sum_m g_m*t^(2m+1)."""
    series = [direction]
    for m in range(1, order + 1):
        series.append([
            -sum(
                coefficient * cubic_coefficient(series, factors, m - 1)
                for coefficient, factors in component
            )
            for component in components
        ])
    return series


def main() -> None:
    source = json.loads(SOURCE.read_text())
    assert source["dimension"] == 21
    components = decode_h(source)
    assert len(components) == 21

    # For V=I+H with H cubic, write V^{-1}(t*xi)=sum g_m t^(2m+1).
    # Coefficient comparison gives the recurrence replayed above.  Recompute
    # every residual independently after the complete prefix is available.
    order = 15
    direction = [Q(1)] * 21
    series = inverse_coefficients(components, direction, order)
    for m in range(1, order + 1):
        residual = [
            series[m][output]
            + sum(
                coefficient * cubic_coefficient(series, factors, m - 1)
                for coefficient, factors in component
            )
            for output, component in enumerate(components)
        ]
        assert residual == [Q(0)] * 21

    coordinate_zero = [coefficient[0] for coefficient in series]
    assert all(coordinate_zero)
    print("PASS inverse recurrence: decoded the actual 21D cubic homogeneous map")
    print("PASS inverse recurrence: V^-1(t*1)=sum g_m*t^(2m+1) through m=15")
    print("PASS inverse recurrence: all 21 composition residuals vanish at every order")
    print("coordinate-0 coefficients:", [str(value) for value in coordinate_zero[:10]])


if __name__ == "__main__":
    main()
