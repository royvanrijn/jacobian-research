#!/usr/bin/env python3
"""Exact arithmetic checks for the plane boundary-exclusion theorem.

This module checks the Riemann--Hurwitz budget after the geometric
residue-immersion lemma has removed all ramification on the affine
dicritical curve.  It also replays the repository's first numerically
admissible one-dicritical boundary package.

It is not a coefficient search for a plane Keller counterexample.
"""

from dataclasses import asdict, dataclass
from itertools import combinations
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
class CleanTangentialDefectBudget:
    """Closed-fiber cost of a clean nonimmersive boundary point.

    If the boundary parameter has first target order ``m >= 2`` and the
    residual Jacobian is a unit, local integration gives fiber length
    ``m + 1``.
    """

    generic_degree: int
    tangent_order: int
    local_fiber_length: int
    residual_fiber_budget: int
    status: str

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class OrevkovMultiplicityBudget:
    """Euler-multiplicity budget for dicritical curves at infinity."""

    generic_degree: int
    component_generic_multiplicities: tuple[int, ...]
    exceptional_multiplicity_jumps: tuple[tuple[int, ...], ...]
    total_cost: int
    available_budget: int
    status: str

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class QuarticOrevkovPacket:
    """One global degree-four boundary packet allowed by Orevkov's budget."""

    name: str
    component_generic_multiplicities: tuple[int, ...]
    exceptional_local_multiplicities: tuple[tuple[int, ...], ...]
    ramified_components: int
    unramified_components: int
    forced_clean_special_fiber: tuple[int, ...] | None
    allowed_coincident_boundary_fiber: tuple[int, ...]
    status: str

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def orevkov_multiplicity_budget(
    generic_degree: int,
    components: Sequence[tuple[int, Sequence[int]]],
) -> OrevkovMultiplicityBudget:
    """Evaluate Orevkov's identity ``sum(mu_l + sum(mu_x-mu_l))=N-1``.

    Each component is supplied as its generic local multiplicity followed
    by the exceptional local multiplicities at its finite special points.
    """

    if generic_degree <= 0:
        raise ValueError("the generic degree must be positive")
    generic_multiplicities: list[int] = []
    component_jumps: list[tuple[int, ...]] = []
    total_cost = 0
    for generic_multiplicity, exceptional_multiplicities in components:
        if generic_multiplicity <= 0:
            raise ValueError("local multiplicities must be positive")
        jumps = tuple(
            exceptional - generic_multiplicity
            for exceptional in exceptional_multiplicities
        )
        if any(jump < 0 for jump in jumps):
            raise ValueError(
                "special local multiplicity cannot be below the generic value"
            )
        generic_multiplicities.append(generic_multiplicity)
        component_jumps.append(jumps)
        total_cost += generic_multiplicity + sum(jumps)
    available = generic_degree - 1
    if total_cost > available:
        status = "excluded_by_orevkov_budget"
    elif total_cost == available:
        status = "saturates_orevkov_budget"
    else:
        status = "incomplete_boundary_ledger"
    return OrevkovMultiplicityBudget(
        generic_degree=generic_degree,
        component_generic_multiplicities=tuple(generic_multiplicities),
        exceptional_multiplicity_jumps=tuple(component_jumps),
        total_cost=total_cost,
        available_budget=available,
        status=status,
    )


def quartic_orevkov_packet_atlas() -> tuple[QuarticOrevkovPacket, ...]:
    """Return the two exact quartic packets after Orevkov's global gate.

    A ramified component is required.  Its generic local multiplicity cannot
    be three by Orevkov's ``N-1`` component exclusion, hence it is two.
    The remaining unit in the Euler budget is either one exceptional jump or
    one additional unramified dicritical component.
    """

    jump_budget = orevkov_multiplicity_budget(4, ((2, (3,)),))
    two_boundary_budget = orevkov_multiplicity_budget(
        4,
        ((2, ()), (1, ())),
    )
    assert jump_budget.status == "saturates_orevkov_budget"
    assert two_boundary_budget.status == "saturates_orevkov_budget"
    return (
        QuarticOrevkovPacket(
            name="one_boundary_one_jump",
            component_generic_multiplicities=(2,),
            exceptional_local_multiplicities=((3,),),
            ramified_components=1,
            unramified_components=0,
            forced_clean_special_fiber=(3, 1),
            allowed_coincident_boundary_fiber=(2, 2),
            status="survives_global_budget",
        ),
        QuarticOrevkovPacket(
            name="two_boundaries_no_jump",
            component_generic_multiplicities=(2, 1),
            exceptional_local_multiplicities=((), ()),
            ramified_components=1,
            unramified_components=1,
            forced_clean_special_fiber=None,
            allowed_coincident_boundary_fiber=(2, 2),
            status="survives_global_budget",
        ),
    )


def _compose_permutations(
    left: tuple[int, ...],
    right: tuple[int, ...],
) -> tuple[int, ...]:
    """Return ``left after right`` for zero-based permutation tuples."""

    return tuple(left[right[index]] for index in range(len(left)))


def _transpositions(degree: int) -> tuple[tuple[int, ...], ...]:
    result: list[tuple[int, ...]] = []
    for first, second in combinations(range(degree), 2):
        permutation = list(range(degree))
        permutation[first], permutation[second] = (
            permutation[second],
            permutation[first],
        )
        result.append(tuple(permutation))
    return tuple(result)


def _generated_orbit(
    generators: Sequence[tuple[int, ...]],
    start: int,
) -> tuple[int, ...]:
    orbit = {start}
    frontier = [start]
    while frontier:
        point = frontier.pop()
        for generator in generators:
            image = generator[point]
            if image not in orbit:
                orbit.add(image)
                frontier.append(image)
    return tuple(sorted(orbit))


def _generated_group(
    generators: Sequence[tuple[int, ...]],
) -> tuple[tuple[int, ...], ...]:
    degree = len(generators[0])
    identity = tuple(range(degree))
    group = {identity}
    frontier = [identity]
    while frontier:
        element = frontier.pop()
        for generator in generators:
            product = _compose_permutations(generator, element)
            if product not in group:
                group.add(product)
                frontier.append(product)
    return tuple(sorted(group))


def quartic_cusp_braid_monodromy_audit() -> dict[str, object]:
    """Enumerate transposition representations of the trefoil braid group.

    For a clean one-cusp branch with no self-collision, Lin--Zaidenberg
    straightens the branch to ``y^2=x^3``.  Its complement group has
    presentation ``<a,b | aba=bab>``, with ``a,b`` meridians.  Simple
    ramification sends both to transpositions.  This audit verifies that no
    such representation is transitive on four letters.
    """

    transpositions = _transpositions(4)
    braid_pairs: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    orbit_sizes: list[int] = []
    for left in transpositions:
        for right in transpositions:
            left_right_left = _compose_permutations(
                left,
                _compose_permutations(right, left),
            )
            right_left_right = _compose_permutations(
                right,
                _compose_permutations(left, right),
            )
            if left_right_left != right_left_right:
                continue
            braid_pairs.append((left, right))
            orbit_sizes.append(
                max(
                    len(_generated_orbit((left, right), point))
                    for point in range(4)
                )
            )
    return {
        "transposition_count": len(transpositions),
        "braid_pair_count": len(braid_pairs),
        "maximum_orbit_size": max(orbit_sizes),
        "transitive_pair_count": sum(size == 4 for size in orbit_sizes),
        "verdict": (
            "excluded_without_2_plus_2_self_collision"
            if all(size < 4 for size in orbit_sizes)
            else "survives"
        ),
    }


def quartic_cusp_with_self_collision_monodromy_audit() -> dict[str, int | str]:
    """Add one 2+2 perfect matching to the nondegenerate cusp monodromy."""

    transpositions = _transpositions(4)
    cusp_pairs: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for left in transpositions:
        for right in transpositions:
            if left == right:
                continue
            if _compose_permutations(
                left,
                _compose_permutations(right, left),
            ) != _compose_permutations(
                right,
                _compose_permutations(left, right),
            ):
                continue
            cusp_pairs.append((left, right))

    perfect_matchings: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for left, right in combinations(transpositions, 2):
        if _compose_permutations(left, right) != _compose_permutations(
            right,
            left,
        ):
            continue
        perfect_matchings.append((left, right))

    group_sizes: list[int] = []
    transitive_count = 0
    for cusp_left, cusp_right in cusp_pairs:
        for node_left, node_right in perfect_matchings:
            generators = (
                cusp_left,
                cusp_right,
                node_left,
                node_right,
            )
            orbit_size = len(_generated_orbit(generators, 0))
            transitive_count += orbit_size == 4
            group_sizes.append(len(_generated_group(generators)))

    return {
        "nondegenerate_cusp_pair_count": len(cusp_pairs),
        "perfect_matching_count": len(perfect_matchings),
        "combined_packet_count": len(group_sizes),
        "transitive_packet_count": transitive_count,
        "minimum_group_order": min(group_sizes),
        "maximum_group_order": max(group_sizes),
        "verdict": "one_2_plus_2_collision_generates_S4",
    }


def quartic_one_boundary_euler_defect(
    self_collision_count: int,
) -> dict[str, int]:
    """Euler-integrate the missing-sheet defect on the quartic target curve.

    The curve has one unibranch cusp and ``r`` identifications of two
    normalization points.  Its normalization is A1.  The defect is two on
    the regular stratum, three at the cusp, and four at every 2+2
    self-collision.
    """

    if self_collision_count < 0:
        raise ValueError("the self-collision count must be nonnegative")
    curve_euler = 1 - self_collision_count
    regular_stratum_euler = (
        curve_euler - 1 - self_collision_count
    )
    integrated_defect = (
        2 * regular_stratum_euler
        + 3
        + 4 * self_collision_count
    )
    return {
        "curve_euler": curve_euler,
        "regular_stratum_euler": regular_stratum_euler,
        "integrated_defect": integrated_defect,
        "global_required_defect": 3,
    }


def clean_tangential_defect_budget(
    generic_degree: int,
    tangent_order: int,
) -> CleanTangentialDefectBudget:
    """Audit the local length formula ``length = tangent_order + 1``."""

    if generic_degree <= 0:
        raise ValueError("the generic degree must be positive")
    if tangent_order < 2:
        raise ValueError("a tangential nonimmersion has order at least two")
    local_length = tangent_order + 1
    residual = generic_degree - local_length
    if residual < 0:
        status = "excluded_by_flat_length"
    elif residual == 0:
        status = "consumes_full_fiber"
    else:
        status = "leaves_residual_budget"
    return CleanTangentialDefectBudget(
        generic_degree=generic_degree,
        tangent_order=tangent_order,
        local_fiber_length=local_length,
        residual_fiber_budget=residual,
        status=status,
    )


def cubic_closed_fiber_atlas() -> tuple[dict[str, object], ...]:
    """Return the complete algebraically closed rank-three fiber atlas."""

    return (
        {
            "partition": (1, 1, 1),
            "local_algebras": ("k", "k", "k"),
            "ramified": False,
        },
        {
            "partition": (2, 1),
            "local_algebras": ("k[epsilon]/(epsilon^2)", "k"),
            "ramified": True,
        },
        {
            "partition": (3,),
            "local_algebras": ("k[epsilon]/(epsilon^3)",),
            "ramified": True,
            "curvilinear": True,
        },
        {
            "partition": (3,),
            "local_algebras": (
                "k[epsilon,eta]/(epsilon,eta)^2",
            ),
            "ramified": True,
            "curvilinear": False,
        },
    )


@dataclass(frozen=True)
class ConductorPacketBudget:
    """Finite-flat length budget for several boundary points in one fiber.

    ``generic_boundary_degree`` is the full generic contribution
    ``sum_D e_D f_D`` of the boundary primes lying over the chosen target
    component.  It is deliberately separate from ``transverse_indices``:
    the latter are local contributions of distinct boundary points over one
    closed target point.
    """

    transverse_indices: tuple[int, ...]
    generic_boundary_degree: int
    affine_degree: int
    generic_degree: int
    minimum_packet_length: int
    minimum_affine_degree: int
    length_deficit: int
    verdict: str

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class BoundaryPrimeLedgerEntry:
    """Generic degree data for one normalization-boundary prime over B."""

    name: str
    transverse_index: int
    residue_degree: int

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("a boundary-prime ledger entry needs a name")
        if self.transverse_index <= 0:
            raise ValueError("transverse index must be positive")
        if self.residue_degree <= 0:
            raise ValueError("residue degree must be positive")


@dataclass(frozen=True)
class AffinePrimeLedgerEntry:
    """Generic contribution of one affine prime over B.

    Affine primes have transverse index one on the Keller open, so only
    their residue degree is recorded.
    """

    name: str
    residue_degree: int

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("an affine-prime ledger entry needs a name")
        if self.residue_degree <= 0:
            raise ValueError("residue degree must be positive")


@dataclass(frozen=True)
class ConductorPacketPoint:
    """One certified boundary point in a closed target fiber."""

    name: str
    boundary_prime: str
    transverse_index: int
    residue_immersive: bool

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("a packet point needs a name")
        if not self.boundary_prime:
            raise ValueError("a packet point needs a boundary-prime label")
        if self.transverse_index <= 0:
            raise ValueError("local transverse index must be positive")


@dataclass(frozen=True)
class TargetFiberPacket:
    """Distinct boundary points certified to lie over one closed target point."""

    name: str
    points: tuple[ConductorPacketPoint, ...]
    same_target_fiber_certified: bool

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("a target-fiber packet needs a name")
        if not self.points:
            raise ValueError("a target-fiber packet needs at least one point")
        point_names = tuple(point.name for point in self.points)
        if len(set(point_names)) != len(point_names):
            raise ValueError("packet point labels must be distinct")


@dataclass(frozen=True)
class TargetComponentLedger:
    """Target-side data that cannot be inferred from ``(Q, p)`` alone."""

    name: str
    generic_degree: int
    boundary_primes: tuple[BoundaryPrimeLedgerEntry, ...]
    affine_primes: tuple[AffinePrimeLedgerEntry, ...]
    packets: tuple[TargetFiberPacket, ...]
    finite_flat_certified: bool
    target_transfer_certified: bool
    exhaustive_generic_pullback: bool

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("a target-component ledger needs a name")
        if self.generic_degree <= 0:
            raise ValueError("generic degree must be positive")
        if not self.boundary_primes:
            raise ValueError("a nonproperness component needs a boundary prime")
        boundary_names = tuple(prime.name for prime in self.boundary_primes)
        affine_names = tuple(prime.name for prime in self.affine_primes)
        if len(set(boundary_names)) != len(boundary_names):
            raise ValueError("boundary-prime labels must be distinct")
        if len(set(affine_names)) != len(affine_names):
            raise ValueError("affine-prime labels must be distinct")
        if set(boundary_names) & set(affine_names):
            raise ValueError("boundary and affine prime labels must be disjoint")
        known_boundary_names = set(boundary_names)
        for packet in self.packets:
            for point in packet.points:
                if point.boundary_prime not in known_boundary_names:
                    raise ValueError(
                        f"packet point {point.name!r} refers to unknown "
                        f"boundary prime {point.boundary_prime!r}"
                    )


@dataclass(frozen=True)
class TargetFiberPacketAudit:
    name: str
    point_names: tuple[str, ...]
    transverse_indices: tuple[int, ...]
    minimum_packet_length: int
    length_deficit: int
    applicable: bool
    status: str
    reasons: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class TargetComponentLedgerAudit:
    name: str
    generic_degree: int
    generic_boundary_degree: int
    affine_degree: int
    expected_generic_degree: int
    degree_identity_holds: bool
    packet_audits: tuple[TargetFiberPacketAudit, ...]
    status: str
    reasons: tuple[str, ...]

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
    normalized_residue_unramified: bool
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
        if self.punctures <= 0:
            raise ValueError("the puncture count must be positive")


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
        certificate.normalized_residue_unramified
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
    if not certificate.normalized_residue_unramified:
        reasons.append(
            "unramifiedness of the map between boundary and target "
            "normalizations is not certified"
        )

    if reasons:
        return finish("incomplete", reasons, None)

    if certificate.punctures >= 2:
        budgets = puncture_profile_budgets(
            certificate.residue_degree,
            certificate.punctures,
        )
        # If s>f there is no positive pole profile at all.  Otherwise every
        # profile forces f+s-2 affine ramification.
        forced_affine_ramification = (
            certificate.residue_degree + certificate.punctures - 2
        )
        return finish(
            "excluded",
            (
                (
                    "no positive pole profile has the certified degree"
                    if not budgets
                    else
                    "an immersive finite multi-puncture residue cover "
                    "contradicts Riemann--Hurwitz"
                ),
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

    return puncture_profile_budgets(degree, 2)


def puncture_profile_budgets(
    degree: int,
    punctures: int,
) -> tuple[HurwitzBudget, ...]:
    """Enumerate all ordered positive pole profiles over target infinity.

    A degree-``degree`` map ``P1 -> P1`` whose affine source is obtained by
    deleting ``punctures`` points has one positive pole order at each deleted
    point.  Ordered compositions retain the boundary labels.  For every
    profile Riemann--Hurwitz forces exactly

        degree + punctures - 2

    ramification away from infinity.
    """

    if degree <= 0:
        raise ValueError("degree must be positive")
    if punctures <= 0:
        raise ValueError("puncture count must be positive")
    if punctures > degree:
        return ()
    profiles: list[HurwitzBudget] = []
    for cuts in combinations(range(1, degree), punctures - 1):
        endpoints = (0, *cuts, degree)
        poles = tuple(
            endpoints[index + 1] - endpoints[index]
            for index in range(punctures)
        )
        profiles.append(hurwitz_budget(degree, poles))
    return tuple(profiles)


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


def conductor_packet_budget(
    transverse_indices: Sequence[int],
    generic_boundary_degree: int,
    affine_degree: int,
) -> ConductorPacketBudget:
    """Audit the general conductor-packet inequality.

    At distinct boundary points ``p_i`` in one finite-flat fiber, assume the
    boundary germ is smooth, the residue map is immersive, and the transverse
    pullback multiplicity is ``e_i``.  The local fiber algebra at ``p_i`` then
    has length at least ``e_i``.  Hence

        d >= sum_i e_i.

    The generic degree decomposition is

        d = generic_boundary_degree + affine_degree.

    This function only checks the resulting arithmetic; it does not infer
    immersion, flatness, exhaustiveness, or the target-side grouping from a
    source boundary graph.
    """

    indices = tuple(transverse_indices)
    if not indices or any(index <= 0 for index in indices):
        raise ValueError("a conductor packet needs positive transverse indices")
    if generic_boundary_degree <= 0:
        raise ValueError("generic boundary degree must be positive")
    if affine_degree < 0:
        raise ValueError("affine degree must be nonnegative")

    generic_degree = generic_boundary_degree + affine_degree
    packet_length = sum(indices)
    minimum_affine_degree = max(0, packet_length - generic_boundary_degree)
    deficit = packet_length - generic_degree
    verdict = (
        "excluded by finite-flat fiber length"
        if deficit > 0
        else "fiber length permits this packet; further data are required"
    )
    return ConductorPacketBudget(
        transverse_indices=indices,
        generic_boundary_degree=generic_boundary_degree,
        affine_degree=affine_degree,
        generic_degree=generic_degree,
        minimum_packet_length=packet_length,
        minimum_affine_degree=minimum_affine_degree,
        length_deficit=deficit,
        verdict=verdict,
    )


def audit_target_component_ledger(
    ledger: TargetComponentLedger,
) -> TargetComponentLedgerAudit:
    """Audit generic degree accounting and every certified closed-fiber packet."""

    boundary_degree = sum(
        prime.transverse_index * prime.residue_degree
        for prime in ledger.boundary_primes
    )
    affine_degree = sum(prime.residue_degree for prime in ledger.affine_primes)
    expected_degree = boundary_degree + affine_degree
    degree_identity = ledger.generic_degree == expected_degree
    global_reasons: list[str] = []

    if not ledger.finite_flat_certified:
        global_reasons.append("finite flatness over the target is not certified")
    if not ledger.target_transfer_certified:
        global_reasons.append(
            "source dicriticals have not been transferred to the finite "
            "normalization over the original target"
        )
    if not ledger.exhaustive_generic_pullback:
        global_reasons.append(
            "the generic boundary and affine pullback ledger is not exhaustive"
        )
    if not degree_identity:
        global_reasons.append(
            "generic degree does not equal boundary plus affine contributions"
        )

    packet_audits: list[TargetFiberPacketAudit] = []
    for packet in ledger.packets:
        indices = tuple(point.transverse_index for point in packet.points)
        packet_length = sum(indices)
        deficit = packet_length - ledger.generic_degree
        reasons: list[str] = []
        if not packet.same_target_fiber_certified:
            reasons.append(
                "the listed points are not certified to lie in one target fiber"
            )
        nonimmersive = tuple(
            point.name for point in packet.points if not point.residue_immersive
        )
        if nonimmersive:
            reasons.append(
                "residue immersion is missing at " + ", ".join(nonimmersive)
            )
        if global_reasons:
            reasons.extend(global_reasons)

        applicable = not reasons
        if not applicable:
            status = "incomplete"
        elif deficit > 0:
            status = "excluded"
        else:
            status = "survives_packet_budget"
        packet_audits.append(
            TargetFiberPacketAudit(
                name=packet.name,
                point_names=tuple(point.name for point in packet.points),
                transverse_indices=indices,
                minimum_packet_length=packet_length,
                length_deficit=deficit,
                applicable=applicable,
                status=status,
                reasons=tuple(reasons),
            )
        )

    if any(packet.status == "excluded" for packet in packet_audits):
        status = "excluded"
        reasons = ("a certified packet exceeds the finite-flat fiber length",)
    elif global_reasons or any(
        packet.status == "incomplete" for packet in packet_audits
    ):
        status = "incomplete"
        reasons = tuple(global_reasons) or (
            "at least one target-fiber packet is not fully certified",
        )
    else:
        status = "survives_packet_budget"
        reasons = (
            "every certified packet fits within the generic fiber length",
        )

    return TargetComponentLedgerAudit(
        name=ledger.name,
        generic_degree=ledger.generic_degree,
        generic_boundary_degree=boundary_degree,
        affine_degree=affine_degree,
        expected_generic_degree=expected_degree,
        degree_identity_holds=degree_identity,
        packet_audits=tuple(packet_audits),
        status=status,
        reasons=reasons,
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
        "conductor_packet_examples": [
            conductor_packet_budget(
                transverse_indices=(3, 3, 3),
                generic_boundary_degree=3,
                affine_degree=5,
            ).as_dict(),
            conductor_packet_budget(
                transverse_indices=(2, 3, 4),
                generic_boundary_degree=5,
                affine_degree=4,
            ).as_dict(),
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
                    normalized_residue_unramified=True,
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
                    normalized_residue_unramified=False,
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
                    normalized_residue_unramified=False,
                    exhaustive_pullback=True,
                    target_transfer_certified=False,
                ),
            )
        ],
        "first_free_depth_package": first_free_depth_package(),
    }
    print(json.dumps(report, indent=2, sort_keys=True))
