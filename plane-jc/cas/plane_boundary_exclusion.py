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


@dataclass(frozen=True)
class OneDicriticalNormalizationCertificate:
    """Finite-normalization data needed by the conductor theorem.

    This is deliberately not constructible from a source boundary tree
    alone.  ``target_transfer_certified`` records that the source dicritical,
    its punctures, and its local ramification data have been transported to
    the finite Zariski--Main normalization over the original affine target.
    """

    name: str
    generic_degree: int
    transverse_index: int
    residue_degree: int
    affine_degree: int
    punctures: int
    single_normalization_boundary: bool
    log_pure: bool
    exhaustive_pullback: bool
    target_transfer_certified: bool
    one_place_at_infinity: bool = True
    line_component_obstruction_available: bool = True

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("a normalization certificate needs a name")
        if self.generic_degree <= 0:
            raise ValueError("generic degree must be positive")
        if self.transverse_index <= 0:
            raise ValueError("transverse index must be positive")
        if self.residue_degree <= 0:
            raise ValueError("residue degree must be positive")
        if self.affine_degree <= 0:
            raise ValueError("affine degree must be positive")
        if self.punctures not in (1, 2):
            raise ValueError("the restricted gate accepts one or two punctures")


@dataclass(frozen=True)
class OneDicriticalNormalizationAudit:
    """Result of the finite-normalization residue/conductor gate."""

    name: str
    degree_identity_holds: bool
    boundary_degree: int
    affine_degree: int
    expected_generic_degree: int
    residue_immersion_certified: bool
    hurwitz_forced_affine_ramification: int | None
    minimum_conductor_collision_length: int
    conductor_length_deficit: int
    sheet_deficient: bool
    status: str
    reasons: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def audit_one_dicritical_normalization(
    certificate: OneDicriticalNormalizationCertificate,
) -> OneDicriticalNormalizationAudit:
    """Apply the theorem without inferring missing target-side data.

    The generic degree identity is

    ``d = e * f + a``,

    where ``a`` is the total contribution of affine primes over the
    nonproper-value curve.  Once residue immersion forces ``f=1``, a
    nonnormal image needs a conductor collision of length at least ``2e``.
    A normal image is excluded by the affine-line component theorem.
    """

    boundary_degree = (
        certificate.transverse_index * certificate.residue_degree
    )
    expected_degree = boundary_degree + certificate.affine_degree
    degree_identity = certificate.generic_degree == expected_degree
    conductor_length = 2 * certificate.transverse_index
    conductor_deficit = conductor_length - certificate.generic_degree
    sheet_deficient = (
        certificate.affine_degree < certificate.transverse_index
    )
    residue_immersion = (
        certificate.log_pure
        and certificate.exhaustive_pullback
        and certificate.target_transfer_certified
    )
    reasons: list[str] = []

    def finish(
        status: str,
        audit_reasons: Sequence[str],
        forced_affine_ramification: int | None,
        *,
        immersion_certified: bool = residue_immersion,
    ) -> OneDicriticalNormalizationAudit:
        return OneDicriticalNormalizationAudit(
            name=certificate.name,
            degree_identity_holds=degree_identity,
            boundary_degree=boundary_degree,
            affine_degree=certificate.affine_degree,
            expected_generic_degree=expected_degree,
            residue_immersion_certified=immersion_certified,
            hurwitz_forced_affine_ramification=forced_affine_ramification,
            minimum_conductor_collision_length=conductor_length,
            conductor_length_deficit=conductor_deficit,
            sheet_deficient=sheet_deficient,
            status=status,
            reasons=tuple(audit_reasons),
        )

    if not certificate.target_transfer_certified:
        reasons.append(
            "source data have not been transferred to the finite "
            "normalization over the original target"
        )
    if not certificate.single_normalization_boundary:
        reasons.append(
            "one source dicritical has not been shown to give one "
            "finite-normalization boundary prime"
        )
    if not certificate.exhaustive_pullback:
        reasons.append("the generic pullback over the image curve is incomplete")
    if not degree_identity:
        reasons.append(
            "generic degree does not equal boundary plus affine contributions"
        )
    if not certificate.one_place_at_infinity:
        reasons.append("the one-place-at-infinity input is not certified")
    if not certificate.log_pure:
        reasons.append("residual special-point ramification is present")

    if reasons:
        return finish("incomplete", reasons, None)

    if certificate.punctures == 2:
        budgets = two_puncture_budgets(certificate.residue_degree)
        if certificate.residue_degree == 1:
            # A finite map with two punctures cannot have degree one: both
            # punctures must have positive pole order.
            forced_affine_ramification = 1
        else:
            forced_affine_ramification = min(
                budget.forced_affine_ramification for budget in budgets
            )
        return finish(
            "excluded",
            (
                "an immersive finite two-puncture residue cover contradicts "
                "Riemann--Hurwitz",
            ),
            forced_affine_ramification,
            immersion_certified=True,
        )

    budget = one_puncture_budget(certificate.residue_degree)
    forced_affine_ramification = budget.forced_affine_ramification
    if forced_affine_ramification:
        return finish(
            "excluded",
            (
                "an immersive finite A1-to-A1 residue map must have degree one",
            ),
            forced_affine_ramification,
            immersion_certified=True,
        )

    if certificate.transverse_index == 1:
        return finish(
            "excluded",
            (
                "index one and log purity would make the finite cover etale",
            ),
            0,
            immersion_certified=True,
        )

    if not certificate.line_component_obstruction_available:
        return finish(
            "incomplete",
            ("the normal-image branch has not been excluded",),
            0,
            immersion_certified=True,
        )

    if conductor_deficit > 0:
        return finish(
            "excluded",
            (
                "a normal image is a forbidden affine-line component",
                "a nonnormal image needs a conductor fiber of length at least 2e",
            ),
            0,
            immersion_certified=True,
        )

    return finish(
        "survives_conductor_budget",
        (
            "the affine sheet contribution can pay for conductor gluing",
        ),
        0,
        immersion_certified=True,
    )


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
        "typed_normalization_examples": [
            audit_one_dicritical_normalization(certificate).as_dict()
            for certificate in (
                OneDicriticalNormalizationCertificate(
                    name="primitive minimal link",
                    generic_degree=3,
                    transverse_index=2,
                    residue_degree=1,
                    affine_degree=1,
                    punctures=1,
                    single_normalization_boundary=True,
                    log_pure=True,
                    exhaustive_pullback=True,
                    target_transfer_certified=True,
                ),
                OneDicriticalNormalizationCertificate(
                    name="72,108 Case 1 numerical preview",
                    generic_degree=29,
                    transverse_index=3,
                    residue_degree=1,
                    affine_degree=26,
                    punctures=1,
                    single_normalization_boundary=False,
                    log_pure=True,
                    exhaustive_pullback=True,
                    target_transfer_certified=False,
                ),
                OneDicriticalNormalizationCertificate(
                    name="72,108 Case 2 numerical preview",
                    generic_degree=29,
                    transverse_index=5,
                    residue_degree=1,
                    affine_degree=24,
                    punctures=1,
                    single_normalization_boundary=False,
                    log_pure=False,
                    exhaustive_pullback=True,
                    target_transfer_certified=False,
                ),
            )
        ],
        "first_free_depth_package": first_free_depth_package(),
    }
    print(json.dumps(report, indent=2, sort_keys=True))
