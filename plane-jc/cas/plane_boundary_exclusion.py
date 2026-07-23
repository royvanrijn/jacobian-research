#!/usr/bin/env python3
"""Exact arithmetic checks for the plane boundary-exclusion theorem.

This module checks the Riemann--Hurwitz budget after the geometric
residue-immersion lemma has removed all ramification on the affine
dicritical curve.  It also replays the repository's first numerically
admissible one-dicritical boundary package.

It is not a coefficient search for a plane Keller counterexample.
"""

from dataclasses import asdict, dataclass
from typing import Sequence

from boundary_lattice_prefilter import (
    boundary_intersection_matrix,
    standard_completion,
)
from intrinsic_a2_boundary import (
    IntrinsicA2Boundary,
    audit_a2_boundary,
    audit_keller_pole_vector,
)


@dataclass(frozen=True)
class HurwitzBudget:
    degree: int
    pole_orders: tuple[int, ...]
    required_total_ramification: int
    puncture_ramification: int
    forced_affine_ramification: int

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class ConductorCollisionBudget:
    transverse_index: int
    affine_degree: int
    generic_degree: int
    minimum_collision_length: int
    length_deficit: int
    verdict: str

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def hurwitz_budget(degree: int, pole_orders: Sequence[int]) -> HurwitzBudget:
    """Return the ramification forced away from the punctures."""

    poles = tuple(pole_orders)
    if degree <= 0:
        raise ValueError("degree must be positive")
    if not poles or any(order <= 0 for order in poles):
        raise ValueError("every puncture must have positive pole order")
    if sum(poles) != degree:
        raise ValueError("pole orders must sum to the cover degree")

    total = 2 * degree - 2
    at_punctures = sum(order - 1 for order in poles)
    return HurwitzBudget(
        degree=degree,
        pole_orders=poles,
        required_total_ramification=total,
        puncture_ramification=at_punctures,
        forced_affine_ramification=total - at_punctures,
    )


def one_puncture_budget(degree: int) -> HurwitzBudget:
    """Budget for a finite A1 -> A1 residue map."""

    return hurwitz_budget(degree, (degree,))


def two_puncture_budgets(degree: int) -> tuple[HurwitzBudget, ...]:
    """All pole splits for a finite Gm -> A1 residue map."""

    if degree < 2:
        return ()
    return tuple(
        hurwitz_budget(degree, (left, degree - left))
        for left in range(1, degree)
    )


def conductor_collision_budget(
    transverse_index: int,
    affine_degree: int = 1,
) -> ConductorCollisionBudget:
    """Audit a two-point conductor collision.

    The generic decomposition has degree ``e + affine_degree``.  Two
    distinct boundary points in one fiber contribute at least ``2e``.
    """

    if transverse_index <= 0:
        raise ValueError("transverse index must be positive")
    if affine_degree <= 0:
        raise ValueError("affine degree must be positive")
    degree = transverse_index + affine_degree
    collision_length = 2 * transverse_index
    deficit = collision_length - degree
    if transverse_index == 1:
        verdict = "transverse index one would make the cover finite etale"
    elif deficit > 0:
        verdict = "excluded by finite-flat fiber length"
    else:
        verdict = "fiber length permits gluing; further data are required"
    return ConductorCollisionBudget(
        transverse_index=transverse_index,
        affine_degree=affine_degree,
        generic_degree=degree,
        minimum_collision_length=collision_length,
        length_deficit=deficit,
        verdict=verdict,
    )


def first_free_depth_package() -> dict[str, object]:
    """Replay the accepted degree-six source-boundary package."""

    configuration, initial_form, _ = standard_completion("A2")
    for parent in ("L", "E1", "E2"):
        configuration = configuration.blow_up((parent,))

    boundary = IntrinsicA2Boundary(
        names=configuration.names,
        intersection_matrix=boundary_intersection_matrix(
            configuration, initial_form
        ),
        genera=(0,) * len(configuration.names),
    )
    boundary_audit = audit_a2_boundary(boundary)
    pole_audit = audit_keller_pole_vector(
        boundary_audit, (3, 2, 1, 0), require_nonproper=True
    )

    return {
        "passes_intrinsic_gates": pole_audit.passes,
        "canonical_coefficients": tuple(
            int(value) for value in boundary_audit.canonical_coefficients
        ),
        "pole_vector": pole_audit.pole_vector,
        "hyperplane_intersections": pole_audit.hyperplane_intersections,
        "geometric_degree": pole_audit.geometric_degree,
        "dicritical_candidates": pole_audit.dicritical_candidates,
        "dicritical_target_line_degree": (
            pole_audit.hyperplane_intersections[-1]
        ),
        "residue_conclusion": (
            "degree-one residue map onto a line; excluded for a Keller map"
        ),
    }


if __name__ == "__main__":
    import json

    report = {
        "one_puncture_degrees_1_to_8": [
            one_puncture_budget(degree).as_dict()
            for degree in range(1, 9)
        ],
        "two_puncture_degrees_2_to_8": [
            budget.as_dict()
            for degree in range(2, 9)
            for budget in two_puncture_budgets(degree)
        ],
        "conductor_collision_degrees_2_to_8": [
            conductor_collision_budget(
                transverse_index=degree - 1,
                affine_degree=1,
            ).as_dict()
            for degree in range(2, 9)
        ],
        "first_free_depth_package": first_free_depth_package(),
    }
    print(json.dumps(report, indent=2, sort_keys=True))
