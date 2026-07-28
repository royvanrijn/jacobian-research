#!/usr/bin/env python3
"""Audit two finite provenance-compression routes for the stored collision.

This is deliberately not a minimum-pair proof.  It checks:

1. literal inverse-recurrence dependency pruning cannot remove any of the
   twenty active coordinates of the identity slice; and
2. every currently stored cubicization circuit with an explicit collision
   has at least twenty active nonidentity outputs after all collision-constant
   identity outputs are specialized.
"""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts" / "generated-results"
SOURCE = ARTIFACTS / "essential_bcw_21_counterexample.json"

STORED_CIRCUITS = (
    "essential_bcw_21_counterexample.json",
    "low_complexity_bcw_21_counterexample.json",
    "constant_kernel_bcw_22_counterexample.json",
    "index_reduced_bcw_22_counterexample.json",
    "rank_compressed_bcw_24_counterexample.json",
    "rank_reduced_bcw_24_counterexample.json",
    "shared_bcw_33_counterexample.json",
    "long_bcw_79_counterexample.json",
    "cubic_homogeneous_counterexample.json",
)


def dependency_closure(dependencies: list[set[int]], seed: int) -> set[int]:
    closure = {seed}
    frontier = [seed]
    while frontier:
        output = frontier.pop()
        for variable in dependencies[output]:
            if variable not in closure:
                closure.add(variable)
                frontier.append(variable)
    return closure


def active_identity_slice_dimension(source: dict[str, object]) -> tuple[int, list[int]]:
    components = source["H"]
    points = source["collision_points"]
    identity_outputs = []
    for index, component in enumerate(components):
        if component:
            continue
        values = {point[index] for point in points}
        if len(values) == 1:
            identity_outputs.append(index)
    return len(components) - len(identity_outputs), identity_outputs


def foundational_divergence() -> dict[tuple[int, int, int], Fraction]:
    """Return div(N) for the normalized foundational map K=I+N."""
    components = (
        {
            (3, 0, 1): Fraction(-1, 2),
            (2, 1, 0): Fraction(-3, 2),
        },
        {
            (3, 2, 1): Fraction(3),
            (2, 3, 0): Fraction(9),
            (2, 1, 1): Fraction(6),
            (1, 2, 0): Fraction(12),
            (1, 0, 1): Fraction(3),
        },
        {
            (3, 3, 1): Fraction(1),
            (2, 4, 0): Fraction(3),
            (2, 2, 1): Fraction(3),
            (1, 3, 0): Fraction(7),
            (1, 1, 1): Fraction(3),
            (0, 2, 0): Fraction(4),
        },
    )
    answer: dict[tuple[int, int, int], Fraction] = {}
    for variable, component in enumerate(components):
        for exponent, coefficient in component.items():
            power = exponent[variable]
            if not power:
                continue
            derivative = list(exponent)
            derivative[variable] -= 1
            monomial = tuple(derivative)
            answer[monomial] = answer.get(monomial, Fraction(0)) + power * coefficient
    return {monomial: coefficient for monomial, coefficient in answer.items() if coefficient}


def main() -> None:
    divergence = foundational_divergence()
    assert divergence == {
        (3, 3, 0): Fraction(1),
        (3, 1, 1): Fraction(6),
        (2, 2, 0): Fraction(30),
        (2, 0, 1): Fraction(9, 2),
        (1, 1, 0): Fraction(24),
    }

    source = json.loads(SOURCE.read_text())
    assert source["dimension"] == 21
    dependencies = [
        {
            variable
            for term in component
            for variable, _power in term["monomial"]
        }
        for component in source["H"]
    ]

    full_closure = dependency_closure(dependencies, 0)
    assert full_closure == set(range(21))

    # Coordinate 20 is the known identity output and is treated as a scalar
    # after specialization.  Every other coordinate is still required by the
    # literal recurrence for inverse coordinate zero.
    sliced_dependencies = [
        {variable for variable in dependencies[index] if variable != 20}
        for index in range(20)
    ]
    sliced_closure = dependency_closure(sliced_dependencies, 0)
    assert sliced_closure == set(range(20))

    profiles = []
    for filename in STORED_CIRCUITS:
        artifact = json.loads((ARTIFACTS / filename).read_text())
        dimension, identity_outputs = active_identity_slice_dimension(artifact)
        profiles.append((filename, dimension, identity_outputs))
    minimum = min(dimension for _filename, dimension, _outputs in profiles)
    assert minimum == 20

    print("PASS provenance compression: normalized 3D canonical contraction fails at its first pure moment")
    print("  div(N)=x^3*y^3+6*x^3*y*z+30*x^2*y^2+(9/2)*x^2*z+24*x*y")
    print("PASS provenance compression: inverse-coordinate-0 dependency closure is all 21 variables")
    print("PASS provenance compression: after X_20=1, literal recurrence closure is all 20 active variables")
    print("PASS provenance compression: stored circuit identity-slice census")
    for filename, dimension, identity_outputs in profiles:
        print(f"  {filename}: active={dimension}, specialized_identity_outputs={identity_outputs}")
    print("COMPUTATION ONLY: the stored-circuit census has minimum active dimension 20")


if __name__ == "__main__":
    main()
