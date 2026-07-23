#!/usr/bin/env python3
"""Compile certified branch scales into a complete toroidal boundary.

The compiler intentionally does not guess branch scales from Newton corners.
Its input is a chart-aware certificate containing the local monomial ratios
which a normal-form theorem has already proved exhaustive.  For every scale

    [u^a : v^b]

it inserts the primitive equality ray ``(b/g,a/g)``, ``g=gcd(a,b)``, in the
regular fan at the selected boundary point.  Star subdivisions are ordinary
point blowups, so the result includes the proximity graph, every boundary
prime, its Picard class, the full intersection matrix, coordinate valuations,
and elementary tame different/conductor data.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import sympy as sp

from boundary_lattice_prefilter import (
    BoundaryConfiguration,
    LocalizationInvariants,
    boundary_intersection_matrix,
    localization_invariants,
    standard_completion,
)
from intrinsic_a2_boundary import (
    A2BoundaryAudit,
    IntrinsicA2Boundary,
    KellerPoleAudit,
    audit_a2_boundary,
    audit_keller_pole_vector,
)


Ray = tuple[int, int]


def _det(first: Ray, second: Ray) -> int:
    return first[0] * second[1] - first[1] * second[0]


def _primitive(ray: Ray) -> Ray:
    if ray[0] < 0 or ray[1] < 0 or ray == (0, 0):
        raise ValueError(f"a toroidal ray must be nonzero and nonnegative: {ray}")
    divisor = math.gcd(ray[0], ray[1])
    return ray[0] // divisor, ray[1] // divisor


@dataclass(frozen=True)
class BranchScale:
    """A local rational map ``[u^a:v^b]`` at a boundary point."""

    name: str
    u_power: int
    v_power: int

    def __post_init__(self) -> None:
        if self.u_power <= 0 or self.v_power <= 0:
            raise ValueError("branch-scale exponents must be positive")

    @property
    def branch_count(self) -> int:
        return math.gcd(self.u_power, self.v_power)

    @property
    def equality_ray(self) -> Ray:
        divisor = self.branch_count
        return self.v_power // divisor, self.u_power // divisor

    @property
    def common_base_order(self) -> int:
        return math.lcm(self.u_power, self.v_power)

    @property
    def branch_orders(self) -> Ray:
        """Orders of ``(u,v)`` on each normalized binomial branch."""

        return self.equality_ray

    @property
    def semigroup_conductor(self) -> int:
        """Conductor of ``k[[tau^p,tau^q]]`` for one primitive branch."""

        p, q = self.branch_orders
        return (p - 1) * (q - 1)

    @property
    def projection_differents(self) -> Ray:
        """Tame different exponents for the two coordinate projections."""

        p, q = self.branch_orders
        return p - 1, q - 1


@dataclass(frozen=True)
class LaurentTranslationRecord:
    """One local scale proved by a Laurent translation in a Newton case tree.

    At the SNC crossing ``x=0, w=1/y=0``, a translation
    ``y -> y + lambda*x^-q`` selects the branch ``w ~ lambda^-1*x^q``.
    In local parameters ``(u,v)=(w,x)`` this is the scale ``[w:x^q]``.

    These records deliberately do not assert that their case paths have
    already been glued into a complete global boundary certificate.
    """

    name: str
    order: int
    case_scope: str
    source: str
    common_after_branching: bool = False

    def __post_init__(self) -> None:
        if self.order <= 0:
            raise ValueError("a Laurent translation order must be positive")

    @property
    def branch_scale(self) -> BranchScale:
        return BranchScale(
            name=self.name,
            u_power=1,
            v_power=self.order,
        )


@dataclass(frozen=True)
class PullbackMonomial:
    """A target-boundary uniformizer pulled back to ``unit*u^i*v^j``."""

    name: str
    u_power: int
    v_power: int

    def __post_init__(self) -> None:
        if self.u_power < 0 or self.v_power < 0:
            raise ValueError("target-boundary pullback exponents must be nonnegative")
        if self.u_power == 0 and self.v_power == 0:
            raise ValueError("a target-boundary pullback cannot be a unit")

    def order_on(self, ray: Ray) -> int:
        return self.u_power * ray[0] + self.v_power * ray[1]


@dataclass(frozen=True)
class ToroidalCluster:
    """All relative branch scales based at one boundary point."""

    name: str
    root_boundary: str
    scales: tuple[BranchScale, ...]
    transverse_boundary: str | None = None
    pullbacks: tuple[PullbackMonomial, ...] = ()
    local_parameters: tuple[str, str] = ("u", "v")
    source: str = ""

    def __post_init__(self) -> None:
        if not self.scales:
            raise ValueError("a toroidal cluster needs at least one branch scale")
        if len(self.local_parameters) != 2:
            raise ValueError("a cluster needs boundary and transverse parameters")
        if len(set(self.local_parameters)) != 2 or not all(self.local_parameters):
            raise ValueError("local parameter names must be nonempty and distinct")
        scale_names = tuple(scale.name for scale in self.scales)
        if len(set(scale_names)) != len(scale_names):
            raise ValueError("branch-scale names must be distinct within a cluster")
        pullback_names = tuple(pullback.name for pullback in self.pullbacks)
        if len(set(pullback_names)) != len(pullback_names):
            raise ValueError("pullback names must be distinct within a cluster")


@dataclass(frozen=True)
class NewtonBoundaryCertificate:
    """The theorem-bearing input between a Newton chain and its boundary."""

    name: str
    chart: str
    clusters: tuple[ToroidalCluster, ...]
    transformations: tuple[str, ...]
    exhaustive: bool
    corners: tuple[tuple[sp.Rational, sp.Rational], ...] = ()
    theorem_source: str = ""
    missing_data: tuple[str, ...] = ()

    @property
    def frontend_complete(self) -> bool:
        return (
            self.exhaustive
            and bool(self.clusters)
            and bool(self.transformations)
            and bool(self.theorem_source)
            and not self.missing_data
        )


@dataclass(frozen=True)
class BlowupStep:
    index: int
    cluster: str
    exceptional: str
    center: tuple[str, ...]
    ray: Ray


@dataclass(frozen=True)
class BoundaryValuation:
    divisor: str
    ray: Ray
    coordinate_orders: Ray
    pullback_orders: tuple[tuple[str, int], ...]
    tame_differents: tuple[tuple[str, int], ...]


@dataclass(frozen=True)
class ScaleReport:
    name: str
    ratio: tuple[int, int]
    equality_ray: Ray
    dicritical_divisor: str
    common_base_order: int
    residue_degree: int
    branch_count: int
    branch_orders: Ray
    projection_differents: Ray
    semigroup_conductor: int


@dataclass(frozen=True)
class ClusterReport:
    name: str
    root_boundary: str
    transverse_boundary: str | None
    local_parameters: tuple[str, str]
    source: str
    rays: tuple[tuple[Ray, str | None], ...]
    valuations: tuple[BoundaryValuation, ...]
    scales: tuple[ScaleReport, ...]


@dataclass(frozen=True)
class LogBoundaryCompilation:
    source: NewtonBoundaryCertificate
    boundary: BoundaryConfiguration
    initial_intersection_form: sp.Matrix
    expected_unit_rank: int
    blowups: tuple[BlowupStep, ...]
    clusters: tuple[ClusterReport, ...]
    localization: LocalizationInvariants
    intersection_matrix: sp.Matrix

    @property
    def passes_prefilter(self) -> bool:
        return self.localization.passes(self.expected_unit_rank)

    @property
    def intrinsic_a2_audit(self) -> A2BoundaryAudit | None:
        """Run adjunction and Noether checks when the compiled chart is A2."""

        if self.source.chart != "A2":
            return None
        return audit_a2_boundary(
            IntrinsicA2Boundary(
                names=self.boundary.names,
                intersection_matrix=self.intersection_matrix,
                genera=(0,) * len(self.boundary.names),
            )
        )

    def as_dict(self) -> dict[str, object]:
        intrinsic_a2 = self.intrinsic_a2_audit
        return {
            "name": self.source.name,
            "chart": self.source.chart,
            "frontend_complete": self.source.frontend_complete,
            "exhaustive": self.source.exhaustive,
            "theorem_source": self.source.theorem_source,
            "transformations": list(self.source.transformations),
            "corners": [
                [str(first), str(second)] for first, second in self.source.corners
            ],
            "passes_prefilter": self.passes_prefilter,
            "expected_unit_rank": self.expected_unit_rank,
            "boundary_names": list(self.boundary.names),
            "boundary_matrix": [
                [int(value) for value in row]
                for row in self.boundary.class_matrix.tolist()
            ],
            "intersection_matrix": [
                [int(value) for value in row]
                for row in self.intersection_matrix.tolist()
            ],
            "intersection_graph": sorted(
                sorted(edge) for edge in self.boundary.intersections
            ),
            "smith_diagonal": list(self.localization.smith_diagonal),
            "unit_rank": self.localization.unit_rank,
            "picard_free_rank": self.localization.picard_free_rank,
            "picard_torsion": list(self.localization.picard_torsion),
            "intrinsic_a2_boundary": (
                intrinsic_a2.as_dict() if intrinsic_a2 is not None else None
            ),
            "blowups": [
                {
                    "index": step.index,
                    "cluster": step.cluster,
                    "exceptional": step.exceptional,
                    "center": list(step.center),
                    "ray": list(step.ray),
                }
                for step in self.blowups
            ],
            "clusters": [
                {
                    "name": cluster.name,
                    "root_boundary": cluster.root_boundary,
                    "transverse_boundary": cluster.transverse_boundary,
                    "local_parameters": list(cluster.local_parameters),
                    "source": cluster.source,
                    "rays": [
                        {"ray": list(ray), "divisor": divisor}
                        for ray, divisor in cluster.rays
                    ],
                    "valuations": [
                        {
                            "divisor": valuation.divisor,
                            "ray": list(valuation.ray),
                            "coordinate_orders": list(
                                valuation.coordinate_orders
                            ),
                            "pullback_orders": dict(valuation.pullback_orders),
                            "tame_differents": dict(valuation.tame_differents),
                        }
                        for valuation in cluster.valuations
                    ],
                    "scales": [
                        {
                            "name": scale.name,
                            "ratio": list(scale.ratio),
                            "equality_ray": list(scale.equality_ray),
                            "dicritical_divisor": scale.dicritical_divisor,
                            "common_base_order": scale.common_base_order,
                            "residue_degree": scale.residue_degree,
                            "branch_count": scale.branch_count,
                            "branch_orders": list(scale.branch_orders),
                            "projection_differents": list(
                                scale.projection_differents
                            ),
                            "semigroup_conductor": scale.semigroup_conductor,
                        }
                        for scale in cluster.scales
                    ],
                }
                for cluster in self.clusters
            ],
        }


@dataclass(frozen=True)
class FilledBoundaryReport:
    """Localization audit after restoring a temporary Laurent divisor."""

    filled_component: str
    boundary: BoundaryConfiguration
    localization: LocalizationInvariants
    intersection_matrix: sp.Matrix
    expected_unit_rank: int

    @property
    def passes_prefilter(self) -> bool:
        return self.localization.passes(self.expected_unit_rank)

    @property
    def intrinsic_a2_audit(self) -> A2BoundaryAudit | None:
        if self.expected_unit_rank != 0:
            return None
        return audit_a2_boundary(
            IntrinsicA2Boundary(
                names=self.boundary.names,
                intersection_matrix=self.intersection_matrix,
                genera=(0,) * len(self.boundary.names),
            )
        )

    def as_dict(self) -> dict[str, object]:
        intrinsic_a2 = self.intrinsic_a2_audit
        return {
            "filled_component": self.filled_component,
            "boundary_names": list(self.boundary.names),
            "boundary_matrix": [
                [int(value) for value in row]
                for row in self.boundary.class_matrix.tolist()
            ],
            "intersection_matrix": [
                [int(value) for value in row]
                for row in self.intersection_matrix.tolist()
            ],
            "smith_diagonal": list(self.localization.smith_diagonal),
            "unit_rank": self.localization.unit_rank,
            "picard_free_rank": self.localization.picard_free_rank,
            "picard_torsion": list(self.localization.picard_torsion),
            "expected_unit_rank": self.expected_unit_rank,
            "passes_prefilter": self.passes_prefilter,
            "intrinsic_a2_boundary": (
                intrinsic_a2.as_dict() if intrinsic_a2 is not None else None
            ),
        }


@dataclass
class _FanRay:
    ray: Ray
    divisor: str | None


def _contains(left: Ray, target: Ray, right: Ray) -> bool:
    return _det(left, target) >= 0 and _det(target, right) >= 0


def _find_cone(fan: Sequence[_FanRay], target: Ray) -> int:
    for index, (left, right) in enumerate(zip(fan, fan[1:])):
        if _contains(left.ray, target, right.ray):
            return index
    raise ArithmeticError(f"ray {target} is outside the local first quadrant")


def _insert_target(
    target: Ray,
    cluster_name: str,
    fan: list[_FanRay],
    boundary: BoundaryConfiguration,
    steps: list[BlowupStep],
) -> BoundaryConfiguration:
    """Reach one primitive ray by regular star subdivisions."""

    target = _primitive(target)
    while all(entry.ray != target for entry in fan):
        cone = _find_cone(fan, target)
        left, right = fan[cone], fan[cone + 1]
        if _det(left.ray, right.ray) != 1:
            raise ArithmeticError("the current toroidal fan is not regular")
        mediant = (
            left.ray[0] + right.ray[0],
            left.ray[1] + right.ray[1],
        )
        center = tuple(
            divisor
            for divisor in (left.divisor, right.divisor)
            if divisor is not None
        )
        if not center:
            raise ArithmeticError("a toroidal blowup center missed the boundary")
        exceptional = f"E{len(steps) + 1}"
        boundary = boundary.blow_up(center, exceptional)
        steps.append(
            BlowupStep(
                index=len(steps) + 1,
                cluster=cluster_name,
                exceptional=exceptional,
                center=center,
                ray=mediant,
            )
        )
        fan.insert(cone + 1, _FanRay(mediant, exceptional))
    return boundary


def compile_log_boundary(
    certificate: NewtonBoundaryCertificate,
) -> LogBoundaryCompilation:
    """Compile a complete certificate or refuse an underdetermined chain."""

    if not certificate.frontend_complete:
        missing = ", ".join(certificate.missing_data) or (
            "exhaustive clusters, transformations, and theorem source"
        )
        raise ValueError(
            f"log-boundary certificate for {certificate.name} is incomplete: {missing}"
        )

    boundary, initial_form, expected_unit_rank = standard_completion(
        certificate.chart
    )
    steps: list[BlowupStep] = []
    reports: list[ClusterReport] = []
    seen_cluster_names: set[str] = set()

    for cluster in certificate.clusters:
        if cluster.name in seen_cluster_names:
            raise ValueError(f"duplicate toroidal cluster name {cluster.name!r}")
        seen_cluster_names.add(cluster.name)
        if cluster.root_boundary not in boundary.names:
            raise ValueError(
                f"cluster {cluster.name}: root {cluster.root_boundary!r} is not "
                "a boundary prime available at this stage"
            )
        if (
            cluster.transverse_boundary is not None
            and cluster.transverse_boundary not in boundary.names
        ):
            raise ValueError(
                f"cluster {cluster.name}: transverse boundary "
                f"{cluster.transverse_boundary!r} is not available"
            )
        if cluster.transverse_boundary == cluster.root_boundary:
            raise ValueError("the two local boundary axes must be distinct")
        if (
            cluster.transverse_boundary is not None
            and frozenset(
                (cluster.root_boundary, cluster.transverse_boundary)
            )
            not in boundary.intersections
        ):
            raise ValueError(
                f"cluster {cluster.name}: the two boundary axes do not meet"
            )

        # The right endpoint is the strict transform of a transverse curve
        # germ, unless the certificate identifies it with a second boundary
        # component at an SNC crossing.
        fan = [
            _FanRay((1, 0), cluster.root_boundary),
            _FanRay((0, 1), cluster.transverse_boundary),
        ]
        scales = sorted(
            cluster.scales,
            key=lambda scale: (
                sp.Rational(scale.equality_ray[1], scale.equality_ray[0]),
                scale.name,
            ),
        )
        for scale in scales:
            boundary = _insert_target(
                scale.equality_ray,
                cluster.name,
                fan,
                boundary,
                steps,
            )

        divisor_for_ray = {
            entry.ray: entry.divisor for entry in fan if entry.divisor is not None
        }
        valuations = tuple(
            BoundaryValuation(
                divisor=entry.divisor,
                ray=entry.ray,
                coordinate_orders=entry.ray,
                pullback_orders=tuple(
                    (pullback.name, pullback.order_on(entry.ray))
                    for pullback in cluster.pullbacks
                ),
                tame_differents=tuple(
                    (
                        pullback.name,
                        max(pullback.order_on(entry.ray) - 1, 0),
                    )
                    for pullback in cluster.pullbacks
                    if pullback.order_on(entry.ray) > 0
                ),
            )
            for entry in fan
            if entry.divisor is not None
        )
        scale_reports = tuple(
            ScaleReport(
                name=scale.name,
                ratio=(scale.u_power, scale.v_power),
                equality_ray=scale.equality_ray,
                dicritical_divisor=divisor_for_ray[scale.equality_ray],
                common_base_order=scale.common_base_order,
                residue_degree=scale.branch_count,
                branch_count=scale.branch_count,
                branch_orders=scale.branch_orders,
                projection_differents=scale.projection_differents,
                semigroup_conductor=scale.semigroup_conductor,
            )
            for scale in cluster.scales
        )
        reports.append(
            ClusterReport(
                name=cluster.name,
                root_boundary=cluster.root_boundary,
                transverse_boundary=cluster.transverse_boundary,
                local_parameters=cluster.local_parameters,
                source=cluster.source,
                rays=tuple((entry.ray, entry.divisor) for entry in fan),
                valuations=valuations,
                scales=scale_reports,
            )
        )

    localization = localization_invariants(boundary.class_matrix.tolist())
    intersection = boundary_intersection_matrix(boundary, initial_form)
    result = LogBoundaryCompilation(
        source=certificate,
        boundary=boundary,
        initial_intersection_form=initial_form,
        expected_unit_rank=expected_unit_rank,
        blowups=tuple(steps),
        clusters=tuple(reports),
        localization=localization,
        intersection_matrix=intersection,
    )
    if not result.passes_prefilter:
        raise ArithmeticError(
            "the compiled complete boundary fails the declared chart's "
            "localization invariants"
        )
    return result


def fill_temporary_boundary(
    compilation: LogBoundaryCompilation,
    component: str,
    expected_unit_rank: int,
) -> FilledBoundaryReport:
    """Recompute the complete boundary after de-localizing one divisor."""

    boundary = compilation.boundary.fill_boundary_component(component)
    localization = localization_invariants(boundary.class_matrix.tolist())
    intersection = boundary_intersection_matrix(
        boundary,
        compilation.initial_intersection_form,
    )
    report = FilledBoundaryReport(
        filled_component=component,
        boundary=boundary,
        localization=localization,
        intersection_matrix=intersection,
        expected_unit_rank=expected_unit_rank,
    )
    if not report.passes_prefilter:
        raise ArithmeticError(
            f"filling {component!r} fails the requested localization invariants"
        )
    return report


def corner_direction(
    first: tuple[sp.Rational, sp.Rational],
    second: tuple[sp.Rational, sp.Rational],
) -> Ray:
    """Return the primitive signed normal equating two Newton weights."""

    dx = sp.Rational(second[0]) - sp.Rational(first[0])
    dy = sp.Rational(second[1]) - sp.Rational(first[1])
    rho, sigma = dy, -dx
    denominator = math.lcm(int(rho.q), int(sigma.q))
    integral = (int(rho * denominator), int(sigma * denominator))
    divisor = math.gcd(abs(integral[0]), abs(integral[1]))
    normal = (integral[0] // divisor, integral[1] // divisor)
    for coordinate in normal:
        if coordinate:
            return normal if coordinate > 0 else (-normal[0], -normal[1])
    raise ValueError("a repeated corner has no direction")


def frontier_72_108_translation_records() -> tuple[LaurentTranslationRecord, ...]:
    """Return the local scales explicitly present in the published case tree.

    The primary proof of Case (8,28) has alternative translations of orders
    two and three, followed in all three intermediate polygon cases by the
    order-four translation.  The records are a certified *local skeleton*,
    not an exhaustive global proximity graph: the proof does not label the
    original completion components or the identifications of the transformed
    SNC crossings with their earlier strict transforms.
    """

    source = (
        "GGHV, Increasing the lower bound for the two-dimensional "
        "Jacobian Conjecture, Proposition Case (8,28)"
    )
    return (
        LaurentTranslationRecord(
            name="pred_1_minus_2",
            order=2,
            case_scope="predecessor direction (1,-2)",
            source=source,
        ),
        LaurentTranslationRecord(
            name="pred_1_minus_3",
            order=3,
            case_scope="one-factor or selected two-factor (1,-3) branch",
            source=source,
        ),
        LaurentTranslationRecord(
            name="common_edge_minus_1_4",
            order=4,
            case_scope="common edge after intermediate cases a, b, c",
            source=source,
            common_after_branching=True,
        ),
    )


@dataclass(frozen=True)
class LaurentTranslationBaseIdealAudit:
    """Exact local distinction between a branch fan and the map base ideal."""

    order: int
    branch_fan_length: int
    base_ideal_length: int
    adapted_parameter: str
    original_generators: tuple[str, str]
    adapted_generators: tuple[str, str]

    @property
    def target_x_orders(self) -> tuple[int, ...]:
        """Orders of the unchanged target coordinate ``x`` on exceptionals."""

        return (1,) * self.base_ideal_length

    @property
    def target_yinf_orders(self) -> tuple[int, ...]:
        """Orders of the resolved target infinity coordinate.

        For the first ``q`` rays ``nu(w,x)=(k,1)``, the numerator and
        denominator orders of
        ``y'=(x^q+lambda*w)/(w*x^q)`` are ``k`` and ``q+k``.
        On the residual rays ``nu(x,t_q)=(1,k)`` they are ``q+k`` and
        ``2q``.  Removing their common base order gives the sequence below.
        """

        q = self.order
        return (q,) * q + tuple(range(q - 1, -1, -1))

    def verify_identity(self) -> bool:
        """Verify the adapted-generator identity symbolically."""

        x, w, t = sp.symbols("x w t")
        lam = sp.symbols("lambda", nonzero=True)
        q = self.order
        adapted_w = (t - x**q) / lam
        second = sp.expand(lam * (w * x**q).subs(w, adapted_w))
        return (
            sp.expand((x**q + lam * w).subs(w, adapted_w) - t) == 0
            and sp.expand(second - (t * x**q - x ** (2 * q))) == 0
        )


@dataclass(frozen=True)
class LaurentTranslationCompositionAudit:
    """Base-ideal audit for a sum of Laurent translations."""

    orders: tuple[int, ...]
    leading_order: int
    adapted_parameter: str
    adapted_generators: tuple[str, str]
    graph_length: int

    def verify_identity(self) -> bool:
        x, w, t = sp.symbols("x w t")
        coefficients = sp.symbols(
            " ".join(f"lambda_{order}" for order in self.orders),
            nonzero=True,
        )
        if len(self.orders) == 1:
            coefficients = (coefficients,)
        F = sp.Add(
            *(
                coefficient * x ** (self.leading_order - order)
                for order, coefficient in zip(self.orders, coefficients)
            )
        )
        adapted_w = (t - x**self.leading_order) / F
        numerator = x**self.leading_order + w * F
        denominator = w * x**self.leading_order
        transformed_denominator = sp.cancel(
            F * denominator.subs(w, adapted_w)
        )
        return (
            sp.cancel(numerator.subs(w, adapted_w) - t) == 0
            and sp.cancel(
                transformed_denominator
                - (
                    t * x**self.leading_order
                    - x ** (2 * self.leading_order)
                )
            )
            == 0
        )


@dataclass(frozen=True)
class UnselectedFactorContinuationAudit:
    """Continuation of the second order-three factor in Newton case (c)."""

    case_presence: tuple[tuple[str, bool], ...]
    residue_difference: str
    after_selected_order_three: str
    after_common_order_four_numerator: str
    order_four_center_value: str
    after_hirzebruch_transition: str
    creates_additional_basepoint: bool
    meets_filled_x_zero: bool

    def verify_identities(self) -> bool:
        """Replay both Laurent shifts and the ``F_4`` transition exactly."""

        x, y, capital_x, capital_y = sp.symbols("x y X Y")
        alpha_1, alpha_2, alpha = sp.symbols(
            "alpha_1 alpha_2 alpha",
            nonzero=True,
        )
        beta = sp.symbols("beta", nonzero=True)
        selected = sp.expand(
            (x**3 * y - alpha_2).subs(
                y,
                y + alpha_1 * x**-3,
            )
        )
        selected = selected.subs(alpha_2, alpha_1 + beta)
        after_order_four = sp.expand(
            selected.subs(y, y + alpha * x**-4)
        )
        numerator = sp.cancel(x * after_order_four)
        final_factor = sp.cancel(
            after_order_four.subs(
                {
                    x: capital_x**-1,
                    y: capital_y * capital_x**4,
                }
            )
        )
        return (
            sp.expand(selected - (x**3 * y - beta)) == 0
            and sp.expand(
                numerator
                - (x**4 * y + alpha - beta * x)
            )
            == 0
            and sp.expand(
                final_factor
                - (capital_x * (capital_y + alpha) - beta)
            )
            == 0
        )


@dataclass(frozen=True)
class FrontierExceptionalPoleAudit:
    """Target-coordinate pole orders on the common eight exceptionals."""

    divisors: tuple[str, ...]
    p_pole_orders: tuple[int, ...]
    q_pole_orders: tuple[int, ...]
    target_infinity_orders: tuple[int, ...]
    first_fan_r_orders: tuple[int, ...]
    residual_r_orders: tuple[int, ...]
    complete_global_pullback: bool

    def verify_orders(self) -> bool:
        """Recompute the valuation minima from the forced edge and vertices."""

        first_r = tuple(
            24 - 8 * k + min(4, k) for k in range(1, 5)
        )
        first_p = tuple(
            -min(2 * r_order, -1) for r_order in first_r
        )
        first_q = tuple(
            -min(3 * r_order, 2 - k, 0)
            for k, r_order in enumerate(first_r, start=1)
        )
        residual_r = tuple(-4 + k for k in range(1, 5))
        residual_p = tuple(
            -min(2 * r_order, -1) for r_order in residual_r
        )
        residual_q = tuple(
            -min(3 * r_order, -2, 0)
            for r_order in residual_r
        )
        p_orders = first_p + residual_p
        q_orders = first_q + residual_q
        return (
            first_r == self.first_fan_r_orders
            and residual_r == self.residual_r_orders
            and p_orders == self.p_pole_orders
            and q_orders == self.q_pole_orders
            and tuple(
                max(p_order, q_order)
                for p_order, q_order in zip(p_orders, q_orders)
            )
            == self.target_infinity_orders
        )


@dataclass(frozen=True)
class MinimalDicriticalExtensionAudit:
    """Enumeration of one-blowup zero-pole boundary extensions."""

    tested_centers: tuple[tuple[str, ...], ...]
    admissible_centers: tuple[tuple[str, ...], ...]
    unique_center: tuple[str, ...]
    pole_vector: tuple[int, ...]
    hyperplane_intersections: tuple[int, ...]
    canonical_coefficients: tuple[int, ...]
    geometric_degree: int
    ordinary_ramification: tuple[int, ...]
    log_ramification: tuple[int, ...]
    dicritical_candidates: tuple[str, ...]


@dataclass(frozen=True)
class TwoStepDicriticalWitness:
    """A positive-pole preparatory blowup followed by a dicritical blowup."""

    first_center: str
    first_pole: int
    second_center: str
    geometric_degree: int
    dicritical_degree: int
    dicritical_canonical_coefficient: int


@dataclass(frozen=True)
class ForcedFirstBlockClusterAudit:
    """The base cluster forced at ``E3 intersect E4`` by equation (J4)."""

    crossing: tuple[str, str]
    local_parameters: tuple[str, str]
    laurent_coordinates: tuple[str, str]
    normalized_sections: tuple[str, str, str]
    reduced_wronskian: str
    u_degree: int
    v_degree: int
    root_count: int
    first_base_order: int
    first_exceptional_pole: int
    child_base_order: int
    child_exceptional_pole: int
    child_target_degree: int
    smooth_e3_basepoints: int
    boundary_audit: A2BoundaryAudit
    pole_audit: KellerPoleAudit


@dataclass(frozen=True)
class Case2BoundaryPackageAudit:
    """The complete forced boundary package for terminal polygon Case 2."""

    first_center: tuple[str, ...]
    local_parameters: tuple[str, str]
    first_base_order: int
    first_exceptional_pole: int
    second_direction: str
    second_base_order: int
    dicritical_pole: int
    dicritical_residue: str
    dicritical_degree: int
    dicritical_ramification_intersection: int
    boundary_audit: A2BoundaryAudit
    pole_audit: KellerPoleAudit


@dataclass(frozen=True)
class Case1BoundaryPackageAudit:
    """The split-factor symmetry selecting the terminal Case-1 package."""

    selected_partition: tuple[int, ...]
    final_chart_relation: str
    reciprocal_relation: str
    first_center: tuple[str, ...]
    first_base_order: int
    first_exceptional_pole: int
    first_exceptional_coordinate: str
    residual_point: str
    adapted_parameter: str
    second_base_order: int
    dicritical_pole: int
    dicritical_residue: str
    dicritical_degree: int
    dicritical_ramification_intersection: int
    boundary_audit: A2BoundaryAudit
    pole_audit: KellerPoleAudit


@dataclass(frozen=True)
class PoissonJacobianRamificationAudit:
    """Actual ramification after the final bracket becomes ``X^2``."""

    terminal_case: str
    jacobian: str
    dicritical_chart: tuple[str, str]
    coordinate_jacobian: str
    pulled_back_jacobian: str
    generic_leading_factor: str
    boundary_names: tuple[str, ...]
    x_valuations: tuple[int, ...]
    canonical_plus_target: tuple[int, ...]
    actual_boundary_ramification: tuple[int, ...]
    affine_ramification_component: str
    affine_component_multiplicity: int
    affine_boundary_intersections: tuple[int, ...]
    dicritical_divisor: str
    dicritical_ramification_coefficient: int
    dicritical_generic_ramification_index: int
    dicritical_boundary_ramification_intersection: int
    dicritical_total_ramification_intersection: int
    constant_keller_ramification_applicable: bool


@dataclass(frozen=True)
class DicriticalResidueDegreeAudit:
    """Normalization-cover possibilities for a degree-(8,12) residue."""

    terminal_case: str
    residue: str
    coordinate_degrees: tuple[int, int]
    hyperplane_degree: int
    a_priori_normalization_cover_degrees: tuple[int, ...]
    excluded_normalization_cover_degrees: tuple[int, ...]
    possible_normalization_cover_degrees: tuple[int, ...]
    possible_image_degrees: tuple[int, ...]
    possible_normalization_differents: tuple[int, ...]
    homogeneous_initial_cover_degree: int
    forced_normalization_cover_degree: int | None
    forced_image_degree: int | None
    generic_ramification_index: int
    valuation_contributions_to_degree_29: tuple[int, ...]
    remaining_valuation_degrees: tuple[int, ...]


@dataclass(frozen=True)
class Case2ResidueStratumExclusionAudit:
    """Exact certificates excluding decomposable Case-2 residues."""

    a_priori_cover_degrees: tuple[int, ...]
    excluded_cover_degrees: tuple[int, ...]
    surviving_cover_degrees: tuple[int, ...]
    right_component_normal_forms: tuple[str, ...]
    equation_scope: str
    constraint_counts: tuple[tuple[int, int], ...]
    constraint_degrees: tuple[tuple[int, tuple[int, ...]], ...]
    singular_input_sha256: tuple[tuple[int, str], ...]
    certificate_command: str
    conclusion: str


@dataclass(frozen=True)
class Case2J1EndpointExclusionAudit:
    """Exact exclusion of the remaining degree-twelve Case-2 residue."""

    forced_endpoint: str
    endpoint_condition: str
    compatibility_equation_count: int
    compatibility_term_counts: tuple[int, ...]
    compatibility_parameter_degrees: tuple[int, ...]
    localization_equation: str
    equation_scope: str
    singular_input_sha256: str
    generic_infinity_characteristic: tuple[int, int]
    generic_infinity_exceptional_rays: tuple[tuple[int, int], ...]
    generic_infinity_self_intersections: tuple[int, ...]
    certificate_command: str
    unit_ideal: bool
    conclusion: str


@dataclass(frozen=True)
class Case2MaximalGcdExclusionAudit:
    """Exact certificate excluding ``deg(gcd(C',G'))=7`` in Case 2."""

    gcd_degree: int
    divisibility_condition: str
    total_remainder_coefficients: int
    selected_remainder_degrees: tuple[int, ...]
    j0_coefficient_degree: int
    selected_constraint_term_counts: tuple[int, ...]
    selected_constraint_parameter_degrees: tuple[int, ...]
    equation_scope: str
    singular_input_sha256: str
    excluded_exact_gcd_degrees: tuple[int, ...]
    surviving_exact_gcd_degrees: tuple[int, ...]
    certificate_command: str
    conclusion: str


@dataclass(frozen=True)
class Case2GcdSixExclusionAudit:
    """Exact pre-compatibility certificate excluding derivative gcd six."""

    gcd_degree: int
    linear_cofactor_parameter: str
    selected_factor_conditions: tuple[str, ...]
    total_g_remainder_coefficients: int
    selected_g_remainder_degrees: tuple[int, ...]
    j0_coefficient_degree: int
    selected_constraint_term_counts: tuple[int, ...]
    selected_constraint_parameter_degrees: tuple[int, ...]
    equation_scope: str
    singular_input_sha256: str
    excluded_exact_gcd_degrees: tuple[int, ...]
    surviving_precompatibility_gcd_degrees: tuple[int, ...]
    certificate_command: str
    conclusion: str


@dataclass(frozen=True)
class Case2LowerJetAudit:
    """The exact gcd reduction of the two lowest Case-2 equations."""

    equations: tuple[str, str]
    derivative_gcd_decomposition: tuple[str, ...]
    j0_conclusion: tuple[str, str]
    multiplier_degree_bound: str
    multiplier_origin_constraint: str
    reduced_j1_equation: str
    wronskian_gcd_factorization: str
    reduced_j1_wronskian_identity: str
    forced_origin_jets: tuple[str, ...]
    cover_degree_gcd_bounds: tuple[tuple[int, int], ...]
    degree_one_gcd_orders: tuple[str, ...]
    degree_one_gcd_excluded_order_pairs: tuple[tuple[int, int], ...]
    linear_series_degree: int
    degree_one_vanishing_sequences: tuple[
        tuple[str, tuple[int, int, int]], ...
    ]
    degree_one_plucker_weights: tuple[tuple[str, int], ...]
    degree_one_affine_wronskian: str
    degree_one_remaining_plucker_weight: int
    gcd_degree_residual_wronskian_degrees: tuple[
        tuple[int, int], ...
    ]
    unresolved_input: str


@dataclass(frozen=True)
class PlaneReturnEdgeAudit:
    """Poisson-square factorization on the final plane-return edge."""

    local_parameters: tuple[str, str]
    edge_sections: tuple[str, str]
    jacobian_top_layer: str
    weighted_wronskian: str
    forced_factorization: tuple[str, str]
    common_factor_degree: int
    case2_root_partition: tuple[int, ...]
    case1_root_partitions: tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class PlaneReturnRootChainAudit:
    """Regular toric resolution over one root of the quartic factor."""

    root_multiplicity: int
    equality_ray: Ray
    blowup_rays: tuple[Ray, ...]
    centers: tuple[tuple[str, ...], ...]
    exceptional_names: tuple[str, ...]
    section_orders: tuple[tuple[int, int, int], ...]
    target_infinity_poles: tuple[int, ...]
    dicritical_divisor: str
    dicritical_residue: tuple[int, int]
    hyperplane_degree: int
    normalization_cover_degree: int
    normalization_different: int
    coordinate_projection_differents: tuple[int, int]
    source_branch_conductor: int
    image_cusp_conductor: int


@dataclass(frozen=True)
class PlaneReturnPartitionAudit:
    """One homogeneous-edge boundary model indexed by a partition of four."""

    partition: tuple[int, ...]
    case_scope: tuple[str, ...]
    first_center: tuple[str, ...]
    first_base_order: int
    first_exceptional: str
    first_exceptional_pole: int
    root_chains: tuple[PlaneReturnRootChainAudit, ...]
    boundary: BoundaryConfiguration
    boundary_audit: A2BoundaryAudit
    pole_audit: KellerPoleAudit


@dataclass(frozen=True)
class HirzebruchTransitionAudit:
    """Toric audit of ``(x,y)->(x^-1,x^n*y)``."""

    degree: int
    chart: str
    valuation_matrix: tuple[tuple[int, int], tuple[int, int]]
    boundary_correspondence: tuple[tuple[str, str], ...]
    determinant: int
    involutive: bool


def laurent_translation_branch_certificate(
    record: LaurentTranslationRecord,
) -> NewtonBoundaryCertificate:
    """Build the regular fan resolving the visible branch ``w~c*x^q``.

    This is not the resolution of the Jonquieres rational map itself.  Its
    base ideal is nonmonomial at the equality residue point; see
    ``laurent_translation_base_ideal_audit``.
    """

    return NewtonBoundaryCertificate(
        name=f"Laurent translation of order {record.order}: {record.name}",
        chart="GmA1",
        clusters=(
            ToroidalCluster(
                name=record.name,
                root_boundary="Yinf",
                transverse_boundary="X0",
                local_parameters=("w=1/y", "x"),
                scales=(record.branch_scale,),
                source=record.case_scope,
            ),
        ),
        transformations=(
            f"y -> y + lambda*x^-{record.order}",
        ),
        theorem_source=record.source,
        exhaustive=True,
    )


def laurent_translation_base_ideal_audit(
    record: LaurentTranslationRecord,
) -> LaurentTranslationBaseIdealAudit:
    """Return the exact adapted monomial ideal of the rational translation.

    With ``w=1/y``, the target homogeneous ``y``-coordinates are

        [x^q + lambda*w : w*x^q].

    Put ``t=x^q+lambda*w``.  Since

        lambda*w*x^q = t*x^q - x^(2q),

    the base ideal is exactly ``(t,x^(2q))``.  It therefore has a
    ``2q``-step monomial base chain, rather than the ``q`` steps in the fan
    of the smooth equality branch ``w-c*x^q``.  The adapted curve ``t=0`` is
    not an original boundary component, so the current boundary-labelled
    compiler deliberately does not pretend to compile this chain.
    """

    q = record.order
    return LaurentTranslationBaseIdealAudit(
        order=q,
        branch_fan_length=q,
        base_ideal_length=2 * q,
        adapted_parameter=f"t=x^{q}+lambda*w",
        original_generators=(f"x^{q}+lambda*w", f"w*x^{q}"),
        adapted_generators=("t", f"x^{2 * q}"),
    )


def laurent_translation_composition_audit(
    orders: Sequence[int],
) -> LaurentTranslationCompositionAudit:
    """Show that lower Laurent orders do not change the final graph.

    A composition of translations is

        y -> y + x^-m F(x),

    where ``m=max(orders)`` and ``F(0)`` is the nonzero coefficient of the
    leading translation.  The same adapted-parameter calculation as for one
    step gives the base ideal ``(t,x^(2m))``.
    """

    normalized = tuple(sorted(set(int(order) for order in orders)))
    if not normalized or normalized[0] <= 0:
        raise ValueError("translation orders must be positive")
    m = normalized[-1]
    terms = "+".join(
        f"lambda_{order}*x^{m - order}" for order in normalized
    )
    return LaurentTranslationCompositionAudit(
        orders=normalized,
        leading_order=m,
        adapted_parameter=f"t=x^{m}+w*({terms})",
        adapted_generators=("t", f"x^{2 * m}"),
        graph_length=2 * m,
    )


def frontier_72_108_unselected_factor_audit(
) -> UnselectedFactorContinuationAudit:
    """Prove that the unselected order-three residue adds no base point.

    In the split case, translating by ``alpha_1*x^-3`` changes the second
    factor to ``x^3*y-beta``, where ``beta=alpha_2-alpha_1`` is nonzero.
    The common order-four translation changes it to

        x^-1 * (x^4*y + alpha - beta*x).

    The numerator restricts to the nonzero common coefficient ``alpha`` at
    the order-four center ``x=x^4*y=0``.  Thus the factor is a unit there.
    Under ``X=x^-1, Y=x^4*y`` it becomes
    ``X*(Y+alpha)-beta``, which also avoids the filled divisor ``X=0``.
    """

    audit = UnselectedFactorContinuationAudit(
        case_presence=(("a", False), ("b", False), ("c", True)),
        residue_difference="beta=alpha_2-alpha_1 != 0",
        after_selected_order_three="x^3*y-beta",
        after_common_order_four_numerator=(
            "x^4*y+alpha-beta*x"
        ),
        order_four_center_value="alpha != 0",
        after_hirzebruch_transition="X*(Y+alpha)-beta",
        creates_additional_basepoint=False,
        meets_filled_x_zero=False,
    )
    if not audit.verify_identities():
        raise AssertionError("unselected-factor continuation identity failed")
    return audit


def frontier_72_108_exceptional_pole_audit(
) -> FrontierExceptionalPoleAudit:
    """Compute the target-infinity pullback on ``E1,...,E8``.

    Write ``w=1/y`` after the common order-four translation.  The forced
    common edge is generated by

        R = x^24*w^-8*(x^4+alpha*w),

    with leading powers ``P~R^2`` and ``Q~R^3``.  The forced low vertices
    contribute ``x^-1`` to ``P`` and ``x^2/w`` (plus a constant) to ``Q``.
    On the first fan, ``nu_k(w,x)=(k,1)``.  On the residual fan,
    ``x^4+alpha*w=x^4*tau`` and ``nu_k(x,tau)=(1,k)``.  The competing
    valuations are distinct whenever the minimum is negative, so the pole
    orders below cannot cancel.

    The result is deliberately partial: it does not supply the two original
    completion-boundary orders or prove that no other global base cluster is
    present.
    """

    audit = FrontierExceptionalPoleAudit(
        divisors=tuple(f"E{index}" for index in range(1, 9)),
        p_pole_orders=(1, 1, 1, 8, 6, 4, 2, 1),
        q_pole_orders=(0, 0, 1, 12, 9, 6, 3, 2),
        target_infinity_orders=(1, 1, 1, 12, 9, 6, 3, 2),
        first_fan_r_orders=(17, 10, 3, -4),
        residual_r_orders=(-3, -2, -1, 0),
        complete_global_pullback=False,
    )
    if not audit.verify_orders():
        raise AssertionError("frontier exceptional pole audit failed")
    return audit


def _frontier_72_108_common_affine_completion(
) -> tuple[LogBoundaryCompilation, FilledBoundaryReport]:
    common_record = frontier_72_108_translation_records()[-1]
    common_graph = compile_log_boundary(
        laurent_translation_graph_certificate(
            common_record,
            chart="GmA1_F4",
        )
    )
    affine_fill = fill_temporary_boundary(
        common_graph,
        component="Xinf",
        expected_unit_rank=0,
    )
    return common_graph, affine_fill


def frontier_72_108_common_graph_pole_audit() -> KellerPoleAudit:
    """Audit the full pole vector on the common translation graph.

    The original pre-transition ``X0`` component has generic pole orders
    ``(1,0)`` for ``(P,Q)``; it is post-transition ``Xinf``.  Generic
    ``Yinf`` has pole orders ``(16,24)``.  Together with the exceptional
    audit this gives

        p = (1,24,1,1,1,12,9,6,3,2).

    The vector passes nefness and ramification effectivity, but every entry
    is positive.  Hence the common graph has no dicritical component and
    cannot be the complete resolution of a hypothetical nonproper Keller
    map.  Any such resolution must contain at least one additional global
    boundary cluster.
    """

    _, affine_fill = _frontier_72_108_common_affine_completion()
    exceptional = frontier_72_108_exceptional_pole_audit()
    pole_vector = (
        1,
        24,
        *exceptional.target_infinity_orders,
    )
    audit = audit_keller_pole_vector(
        affine_fill.intrinsic_a2_audit,
        pole_vector,
        require_nonproper=True,
    )
    expected_failure = (
        "no boundary prime is dicritical over the affine target "
        "(p_i=0 and (Qp)_i>0)"
    )
    if (
        audit.pole_vector != pole_vector
        or audit.geometric_degree != 427
        or audit.hyperplane_intersections
        != (0, 12, 0, 0, 11, 10, 0, 0, 2, 1)
        or audit.dicritical_candidates
        or audit.failures != (expected_failure,)
    ):
        raise AssertionError("unexpected common-graph Keller pole audit")
    return audit


def _frontier_72_108_first_block_boundary(
) -> tuple[LogBoundaryCompilation, BoundaryConfiguration, A2BoundaryAudit]:
    """Attach the source base cluster forced by the first Laurent block."""

    common_graph, affine_fill = _frontier_72_108_common_affine_completion()
    boundary = affine_fill.boundary.blow_up(
        ("E3", "E4"),
        exceptional="B0",
    )
    for index in range(1, 11):
        boundary = boundary.blow_up(
            ("B0",),
            exceptional=f"B{index}",
        )
    intersection = boundary_intersection_matrix(
        boundary,
        common_graph.initial_intersection_form,
    )
    boundary_audit = audit_a2_boundary(
        IntrinsicA2Boundary(
            names=boundary.names,
            intersection_matrix=intersection,
            genera=(0,) * len(boundary.names),
        )
    )
    return common_graph, boundary, boundary_audit


def frontier_72_108_forced_first_block_cluster_audit(
) -> ForcedFirstBlockClusterAudit:
    """Resolve the source base ideal forced at ``E3 intersect E4``.

    Put ``s=Y`` and ``t=1/(XY)`` on the final ``F4`` completion.  These are
    local equations for ``E3`` and ``E4`` and the exact Laurent variables are

        T=s/t,  z=1/s.

    Write ``A=T*U`` and ``D=T^2*V`` in the first block.  After multiplying
    ``[1:P:Q]`` by the local target-infinity equation ``s*t^12``, its lowest
    pieces are

        [s*t^12, t^11*U(s/t), t^10*V(s/t)].

    Thus the crossing has base order ten and its exceptional has pole three.
    Equation (J4) becomes

        U*V + T*(2*U*V' - 3*U'*V) = 1.

    Since ``deg(V)=10``, all ten roots of ``V`` are nonzero and simple, and
    ``U`` is nonzero at each root.  They are ten distinct free basepoints,
    each with ideal ``(epsilon,v)``.  Blowing them up gives ten pole-two
    exceptional curves of target-line degree one.

    The apparent numerical option of a free zero-pole blowup on ``E3`` is
    not a source basepoint: on the smooth torus of ``E3`` the forced low
    monomial ``X`` is a unit after removing its pole.
    """

    T = sp.symbols("T")
    U = sp.Function("U")(T)
    V = sp.Function("V")(T)
    A = T * U
    D = T**2 * V
    reduced = sp.simplify(
        (2 * A * sp.diff(D, T) - 3 * sp.diff(A, T) * D) / T**2
    )
    expected = U * V + T * (
        2 * U * sp.diff(V, T) - 3 * sp.diff(U, T) * V
    )
    if sp.simplify(reduced - expected) != 0:
        raise AssertionError("first-block Wronskian reduction failed")

    _, _, boundary_audit = _frontier_72_108_first_block_boundary()
    pole_vector = (
        frontier_72_108_common_graph_pole_audit().pole_vector
        + (3,)
        + (2,) * 10
    )
    pole_audit = audit_keller_pole_vector(
        boundary_audit,
        pole_vector,
        require_nonproper=True,
    )
    expected_failure = (
        "no boundary prime is dicritical over the affine target "
        "(p_i=0 and (Qp)_i>0)"
    )
    if (
        not boundary_audit.passes
        or pole_audit.geometric_degree != 317
        or pole_audit.hyperplane_intersections[-11:]
        != (0,) + (1,) * 10
        or pole_audit.dicritical_candidates
        or pole_audit.failures != (expected_failure,)
    ):
        raise AssertionError("unexpected forced first-block cluster")

    return ForcedFirstBlockClusterAudit(
        crossing=("E3", "E4"),
        local_parameters=("s=Y", "t=1/(X*Y)"),
        laurent_coordinates=("T=s/t", "z=1/s"),
        normalized_sections=(
            "s*t^12",
            "t^11*U(s/t)+higher",
            "t^10*V(s/t)+higher",
        ),
        reduced_wronskian=(
            "U*V+T*(2*U*V'-3*U'*V)=1"
        ),
        u_degree=7,
        v_degree=10,
        root_count=10,
        first_base_order=10,
        first_exceptional_pole=3,
        child_base_order=1,
        child_exceptional_pole=2,
        child_target_degree=1,
        smooth_e3_basepoints=0,
        boundary_audit=boundary_audit,
        pole_audit=pole_audit,
    )


def frontier_72_108_case2_boundary_package_audit(
) -> Case2BoundaryPackageAudit:
    """Compile the complete forced package for the nonvertical polygon.

    At the smooth point ``X=0`` of ``Yinf``, put ``z=1/Y``.  Multiplying
    ``[1:P:Q]`` by ``z^24`` gives common order twelve.  On the first blowup
    ``X=e*u, z=e``, every Case-2 monomial has second-stage order ``12+k``,
    where ``k=2i-j`` is its Laurent-band index.  Case 2 has only ``k>=0``,
    so the second common base order is again twelve.  The new exceptional
    has pole zero and residue map ``[1:C(r):G(r)]`` of degree twelve.

    This function includes the disjoint forced first-block cluster before
    applying the intrinsic affine-plane and boundary divisor-class gates.
    The separate transformed-Poisson audit corrects the final ``X^2``
    Jacobian.
    """

    p_polygon = ((0, 0), (1, 0), (8, 14), (8, 16))
    q_polygon = ((0, 0), (2, 1), (12, 21), (12, 24))

    def cross(
        first: tuple[int, int],
        second: tuple[int, int],
        point: tuple[int, int],
    ) -> int:
        return (
            (second[0] - first[0]) * (point[1] - first[1])
            - (second[1] - first[1]) * (point[0] - first[0])
        )

    def lattice_points(
        polygon: tuple[tuple[int, int], ...],
    ) -> tuple[tuple[int, int], ...]:
        points = []
        for i in range(max(point[0] for point in polygon) + 1):
            for j in range(max(point[1] for point in polygon) + 1):
                signs = tuple(
                    value
                    for value in (
                        cross(first, second, (i, j))
                        for first, second in zip(
                            polygon,
                            polygon[1:] + polygon[:1],
                        )
                    )
                    if value
                )
                if not signs or all(value > 0 for value in signs) or all(
                    value < 0 for value in signs
                ):
                    points.append((i, j))
        return tuple(points)

    p_support = lattice_points(p_polygon)
    q_support = lattice_points(q_polygon)
    first_orders = (24,) + tuple(
        i + 24 - j for i, j in p_support + q_support
    )
    second_orders = (12,) + tuple(
        2 * i + 12 - j for i, j in p_support + q_support
    )
    residue_p_degree = max(
        i for i, j in p_support if 2 * i + 12 - j == 12
    )
    residue_q_degree = max(
        i for i, j in q_support if 2 * i + 12 - j == 12
    )
    if (
        len(p_support) != 25
        or len(q_support) != 47
        or min(first_orders) != 12
        or min(second_orders) != 12
        or (residue_p_degree, residue_q_degree) != (8, 12)
    ):
        raise AssertionError("terminal Case-2 support audit failed")

    common_graph, boundary, _ = _frontier_72_108_first_block_boundary()
    boundary = boundary.blow_up(("Yinf",), exceptional="T1")
    boundary = boundary.blow_up(("T1",), exceptional="T2")
    intersection = boundary_intersection_matrix(
        boundary,
        common_graph.initial_intersection_form,
    )
    boundary_audit = audit_a2_boundary(
        IntrinsicA2Boundary(
            names=boundary.names,
            intersection_matrix=intersection,
            genera=(0,) * len(boundary.names),
        )
    )
    pole_vector = (
        frontier_72_108_common_graph_pole_audit().pole_vector
        + (3,)
        + (2,) * 10
        + (12, 0)
    )
    pole_audit = audit_keller_pole_vector(
        boundary_audit,
        pole_vector,
        require_nonproper=True,
    )
    dicritical_index = boundary.names.index("T2")
    ramification_intersection = int(
        (
            boundary_audit.boundary.intersection_matrix
            * sp.Matrix(pole_audit.ramification_coefficients)
        )[dicritical_index]
    )
    if (
        not boundary_audit.passes
        or not pole_audit.passes
        or pole_audit.geometric_degree != 29
        or pole_audit.dicritical_candidates != ("T2",)
        or pole_audit.hyperplane_intersections[-2:] != (0, 12)
        or int(boundary_audit.canonical_coefficients[-1]) != 0
        or ramification_intersection != 35
    ):
        raise AssertionError("unexpected terminal Case-2 boundary package")

    return Case2BoundaryPackageAudit(
        first_center=("Yinf",),
        local_parameters=("X", "z=1/Y"),
        first_base_order=12,
        first_exceptional_pole=12,
        second_direction="u=X/z=0",
        second_base_order=12,
        dicritical_pole=0,
        dicritical_residue="[1:C(r):G(r)]",
        dicritical_degree=12,
        dicritical_ramification_intersection=ramification_intersection,
        boundary_audit=boundary_audit,
        pole_audit=pole_audit,
    )


def frontier_72_108_plane_return_edge_audit() -> PlaneReturnEdgeAudit:
    """Recover the missing corner cluster from the final Newton edge.

    At ``X=0, Y=infinity`` put ``z=1/Y`` and ``s=X/z=XY``.  The upper
    parallel edge has

        P_edge=Y^8*A(s),  Q_edge=Y^12*C(s).

    Its top Jacobian layer is

        4*Y^20*(3*A'(s)*C(s)-2*A(s)*C'(s)).

    Unique factorization in characteristic zero therefore gives
    ``A=a*r^2`` and ``C=c*r^3``.  Both edge degrees force ``deg(r)=4``.
    In terminal Case 2 the absent vertical vertices force ``r=s^4``.  In
    terminal Case 1 the constant and leading coefficients of ``r`` are
    nonzero, and its root multiplicities have one of the five partitions of
    four.
    """

    s = sp.symbols("s")
    coefficients = sp.symbols("r0:5")
    r = sum(coefficient * s**index for index, coefficient in enumerate(
        coefficients
    ))
    A = r**2
    C = r**3
    wronskian = sp.expand(
        3 * sp.diff(A, s) * C - 2 * A * sp.diff(C, s)
    )
    if wronskian != 0 or sp.degree(A, s) != 8 or sp.degree(C, s) != 12:
        raise AssertionError("plane-return Poisson-square identity failed")

    partitions = (
        (4,),
        (3, 1),
        (2, 2),
        (2, 1, 1),
        (1, 1, 1, 1),
    )
    return PlaneReturnEdgeAudit(
        local_parameters=("z=1/Y", "s=X/z=X*Y"),
        edge_sections=("Y^8*A(s)", "Y^12*C(s)"),
        jacobian_top_layer="4*Y^20*(3*A'*C-2*A*C')",
        weighted_wronskian="3*A'*C-2*A*C'=0",
        forced_factorization=("A=a*r^2", "C=c*r^3"),
        common_factor_degree=4,
        case2_root_partition=(4,),
        case1_root_partitions=partitions,
    )


def _insert_plane_return_equality_ray(
    boundary: BoundaryConfiguration,
    *,
    root_index: int,
    multiplicity: int,
) -> tuple[
    BoundaryConfiguration,
    tuple[tuple[Ray, str, tuple[str, ...]], ...],
    str,
]:
    """Insert the regular fan from ``D0`` to one quartic-root direction."""

    divisor = math.gcd(multiplicity, 4)
    target = (multiplicity // divisor, 4 // divisor)
    fan: list[tuple[Ray, str | None]] = [
        ((1, 0), "D0"),
        ((0, 1), None),
    ]
    steps: list[tuple[Ray, str, tuple[str, ...]]] = []
    while all(ray != target for ray, _ in fan):
        cone = next(
            index
            for index, ((left, _), (right, _)) in enumerate(
                zip(fan, fan[1:])
            )
            if _contains(left, target, right)
        )
        left_ray, left_divisor = fan[cone]
        right_ray, right_divisor = fan[cone + 1]
        if _det(left_ray, right_ray) != 1:
            raise AssertionError("plane-return fan ceased to be regular")
        ray = (
            left_ray[0] + right_ray[0],
            left_ray[1] + right_ray[1],
        )
        center = tuple(
            name
            for name in (left_divisor, right_divisor)
            if name is not None
        )
        exceptional = f"R{root_index}_{len(steps) + 1}"
        boundary = boundary.blow_up(center, exceptional=exceptional)
        fan.insert(cone + 1, (ray, exceptional))
        steps.append((ray, exceptional, center))

    equality_divisor = next(
        name for ray, name in fan if ray == target and name is not None
    )
    return boundary, tuple(steps), equality_divisor


def frontier_72_108_plane_return_partition_audits(
) -> tuple[PlaneReturnPartitionAudit, ...]:
    """Compile the five homogeneous models allowed by the final edge alone.

    After the first blowup at the smooth ``Yinf`` point, a root of ``r`` of
    multiplicity ``m`` has the monomial homogeneous ideal

        (e^12, e^4*u^(2m), u^(3m)).

    Its three valuation orders agree on the primitive ray
    ``(m/g,4/g)``, where ``g=gcd(m,4)``.  A regular fan containing that ray
    principalizes the ideal.  The equality divisor has residue
    ``[1:t^(2g):t^(3g)]``; all other new toric divisors are contracted.

    A multiple root can acquire transverse terms, so this edge-only
    classification is not by itself a proof that all five models occur as
    resolutions of the full polynomial pair.  Each homogeneous model is
    nevertheless attached to the forced first-block cluster, and its class
    matrix, intersection matrix, canonical vector, target pole vector, and
    ordinary/log ramification vectors are audited.  The separate Case-1
    symmetry audit proves that the primary split-factor source selects the
    partition ``(4,)`` and supplies an exact adapted transverse parameter.
    """

    edge = frontier_72_108_plane_return_edge_audit()
    common_graph, first_block_boundary, _ = (
        _frontier_72_108_first_block_boundary()
    )
    common_poles = frontier_72_108_common_graph_pole_audit()
    base_poles = (
        common_poles.pole_vector
        + (3,)
        + (2,) * 10
    )
    results: list[PlaneReturnPartitionAudit] = []

    for partition in edge.case1_root_partitions:
        boundary = first_block_boundary.blow_up(
            ("Yinf",),
            exceptional="D0",
        )
        pole_vector = list(base_poles + (12,))
        root_chains: list[PlaneReturnRootChainAudit] = []

        for root_index, multiplicity in enumerate(partition, start=1):
            boundary, steps, dicritical = (
                _insert_plane_return_equality_ray(
                    boundary,
                    root_index=root_index,
                    multiplicity=multiplicity,
                )
            )
            section_orders = []
            exceptional_poles = []
            for (a, b), _, _ in steps:
                orders = (
                    12 * a,
                    4 * a + 2 * multiplicity * b,
                    3 * multiplicity * b,
                )
                section_orders.append(orders)
                pole = orders[0] - min(orders)
                exceptional_poles.append(pole)
                pole_vector.append(pole)

            divisor = math.gcd(multiplicity, 4)
            target = (
                multiplicity // divisor,
                4 // divisor,
            )
            if (
                not steps
                or target not in tuple(ray for ray, _, _ in steps)
                or exceptional_poles[-1] != 0
            ):
                raise AssertionError("plane-return equality ray was not resolved")
            root_chains.append(
                PlaneReturnRootChainAudit(
                    root_multiplicity=multiplicity,
                    equality_ray=target,
                    blowup_rays=tuple(ray for ray, _, _ in steps),
                    centers=tuple(center for _, _, center in steps),
                    exceptional_names=tuple(
                        name for _, name, _ in steps
                    ),
                    section_orders=tuple(section_orders),
                    target_infinity_poles=tuple(exceptional_poles),
                    dicritical_divisor=dicritical,
                    dicritical_residue=(2 * divisor, 3 * divisor),
                    hyperplane_degree=3 * divisor,
                    normalization_cover_degree=divisor,
                    normalization_different=divisor - 1,
                    coordinate_projection_differents=(
                        2 * divisor - 1,
                        3 * divisor - 1,
                    ),
                    source_branch_conductor=0,
                    image_cusp_conductor=2,
                )
            )

        intersection = boundary_intersection_matrix(
            boundary,
            common_graph.initial_intersection_form,
        )
        boundary_audit = audit_a2_boundary(
            IntrinsicA2Boundary(
                names=boundary.names,
                intersection_matrix=intersection,
                genera=(0,) * len(boundary.names),
            )
        )
        pole_audit = audit_keller_pole_vector(
            boundary_audit,
            tuple(pole_vector),
            require_nonproper=True,
        )
        expected_dicriticals = tuple(
            chain.dicritical_divisor for chain in root_chains
        )
        dicritical_indices = tuple(
            boundary.names.index(name) for name in expected_dicriticals
        )
        if (
            abs(int(boundary.class_matrix.det())) != 1
            or not boundary_audit.passes
            or not pole_audit.passes
            or pole_audit.geometric_degree != 29
            or pole_audit.dicritical_candidates != expected_dicriticals
            or tuple(
                pole_audit.hyperplane_intersections[index]
                for index in dicritical_indices
            )
            != tuple(chain.hyperplane_degree for chain in root_chains)
        ):
            raise AssertionError(
                f"unexpected plane-return audit for partition {partition}"
            )

        scope = ("edge-only homogeneous model",)
        if partition == edge.case2_root_partition:
            scope += (
                "terminal Case 2",
                "terminal Case 1 after split-factor symmetry",
            )
        results.append(
            PlaneReturnPartitionAudit(
                partition=partition,
                case_scope=scope,
                first_center=("Yinf",),
                first_base_order=12,
                first_exceptional="D0",
                first_exceptional_pole=12,
                root_chains=tuple(root_chains),
                boundary=boundary,
                boundary_audit=boundary_audit,
                pole_audit=pole_audit,
            )
        )

    return tuple(results)


def frontier_72_108_case1_boundary_package_audit(
) -> Case1BoundaryPackageAudit:
    """Select and resolve the actual split-factor Case-1 return package.

    The five partitions in ``frontier_72_108_plane_return_partition_audits``
    classify the homogeneous quartic models allowed by the final edge alone.
    The primary split-factor formula is stronger.  After selecting
    ``alpha_1`` it gives, up to nonzero constants,

        P_edge=(Y+alpha)^8*(X*(Y+alpha)-beta)^8,
        Q_edge=(Y+alpha)^12*(X*(Y+alpha)-beta)^12.

    Hence the common quartic is ``(s-beta)^4`` and the source selects the
    single partition ``(4,)``.

    There is also no unrecorded transverse-jet choice.  Repeat the legal
    construction with ``alpha_2`` selected.  The two final plane charts
    differ by

        Y2=Y1-beta/X+delta.

    With ``z=1/Y1``, this says

        1/Y2 = X*z/(X-beta*z+delta*X*z).

    On the first blowup ``X=e*u, z=e``, the denominator is
    ``e*(u-beta+delta*e*u)``.  Thus
    ``v=u-beta+delta*e*u`` is an exact regular adapted parameter.  Blowing
    up ``e=v=0`` gives ``Y2=r/u`` on the new exceptional, so the residue is
    the restriction of the alternate terminal polynomial pair to ``X=0``.
    Its two degrees are eight and twelve.  This proves base order twelve,
    pole zero, and target-line degree twelve without assuming that transverse
    terms vanish.
    """

    X, z, e, u, beta, delta = sp.symbols(
        "X z e u beta delta",
        nonzero=True,
    )
    reciprocal = sp.cancel(
        1 / (1 / z - beta / X + delta)
    )
    expected_reciprocal = sp.cancel(
        X * z / (X - beta * z + delta * X * z)
    )
    denominator_after_blowup = sp.expand(
        (X - beta * z + delta * X * z).subs(
            {X: e * u, z: e}
        )
        / e
    )
    if (
        sp.cancel(reciprocal - expected_reciprocal) != 0
        or denominator_after_blowup != u - beta + delta * e * u
    ):
        raise AssertionError("alternate split-factor chart identity failed")

    partition = next(
        audit
        for audit in frontier_72_108_plane_return_partition_audits()
        if audit.partition == (4,)
    )
    if (
        partition.case_scope
        != (
            "edge-only homogeneous model",
            "terminal Case 2",
            "terminal Case 1 after split-factor symmetry",
        )
        or len(partition.root_chains) != 1
        or partition.root_chains[0].target_infinity_poles != (0,)
        or partition.root_chains[0].hyperplane_degree != 12
        or partition.pole_audit.dicritical_candidates != ("R1_1",)
    ):
        raise AssertionError("unexpected source-selected Case-1 package")
    dicritical_index = partition.boundary.names.index("R1_1")
    ramification_intersection = int(
        (
            partition.boundary_audit.boundary.intersection_matrix
            * sp.Matrix(partition.pole_audit.ramification_coefficients)
        )[dicritical_index]
    )
    if ramification_intersection != 35:
        raise AssertionError("unexpected Case-1 dicritical ramification degree")

    return Case1BoundaryPackageAudit(
        selected_partition=(4,),
        final_chart_relation="Y2=Y1-beta/X+delta",
        reciprocal_relation=(
            "1/Y2=X*z/(X-beta*z+delta*X*z)"
        ),
        first_center=("Yinf",),
        first_base_order=12,
        first_exceptional_pole=12,
        first_exceptional_coordinate="u=X/z",
        residual_point="u=beta",
        adapted_parameter="v=u-beta+delta*e*u",
        second_base_order=12,
        dicritical_pole=0,
        dicritical_residue="[1:P2(0,r/u):Q2(0,r/u)]",
        dicritical_degree=12,
        dicritical_ramification_intersection=ramification_intersection,
        boundary_audit=partition.boundary_audit,
        pole_audit=partition.pole_audit,
    )


def _frontier_72_108_x_valuations(
    terminal_case: str,
) -> tuple[int, ...]:
    """Orders of the final polynomial coordinate ``X`` on the boundary."""

    common = (-1, 0) + (-1,) * 8
    first_block = (-2,) + (-2,) * 10
    if terminal_case == "Case 1":
        # The return point u=beta is away from the strict transform X=0.
        return common + first_block + (1, 1)
    if terminal_case == "Case 2":
        # The second return blowup is centered where X=0 meets D0.
        return common + first_block + (1, 2)
    raise ValueError(f"unknown terminal case {terminal_case!r}")


def _frontier_72_108_poisson_ramification_audit(
    *,
    terminal_case: str,
    package: Case1BoundaryPackageAudit | Case2BoundaryPackageAudit,
) -> PoissonJacobianRamificationAudit:
    """Correct the constant-Keller representative by ``div(X^2)``.

    The final polynomial pair in Proposition Case (8,28) satisfies
    ``[P,Q]=X^2``, not a nonzero constant.  If ``K+3H`` denotes the
    boundary-supported divisor representative used by the Keller gate, then
    the actual ramification divisor is

        R = K + 3H + div(X^2).

    Writing ``div(X)=L_X+sum nu_i D_i``, where ``L_X`` is the closure of the
    affine line ``X=0``, changes the boundary coefficients by ``2*nu`` and
    adds ``2*L_X``.  The principal-divisor identity also recovers exactly
    where ``L_X`` meets the boundary.
    """

    pole_audit = package.pole_audit
    boundary_audit = package.boundary_audit
    names = boundary_audit.boundary.names
    x_valuations = _frontier_72_108_x_valuations(terminal_case)
    if len(x_valuations) != len(names):
        raise AssertionError("X-valuation vector does not index the boundary")

    class_ramification = pole_audit.ramification_coefficients
    actual_boundary = tuple(
        coefficient + 2 * valuation
        for coefficient, valuation in zip(
            class_ramification,
            x_valuations,
        )
    )
    intersection = boundary_audit.boundary.intersection_matrix
    affine_intersections = tuple(
        -int(value)
        for value in intersection * sp.Matrix(x_valuations)
    )
    dicritical = pole_audit.dicritical_candidates[0]
    index = names.index(dicritical)
    boundary_intersection = int(
        (
            intersection * sp.Matrix(actual_boundary)
        )[index]
    )
    total_intersection = (
        boundary_intersection
        + 2 * affine_intersections[index]
    )

    e, r, beta, delta = sp.symbols("e r beta delta")
    if terminal_case == "Case 1":
        u_chart = (beta + e * r) / (1 + delta * e)
        X_chart = e * u_chart
    else:
        u_chart = None
        X_chart = e**2 * r
    Y_chart = 1 / e
    coordinate_jacobian = sp.cancel(
        sp.diff(X_chart, e) * sp.diff(Y_chart, r)
        - sp.diff(X_chart, r) * sp.diff(Y_chart, e)
    )
    pulled_back_jacobian = sp.cancel(
        X_chart**2 * coordinate_jacobian
    )
    leading_factor = sp.cancel(
        pulled_back_jacobian / e**actual_boundary[index]
    ).subs(e, 0)
    if terminal_case == "Case 1":
        expected_local = (
            1 / (1 + delta * e),
            e**2 * u_chart**2 / (1 + delta * e),
            beta**2,
        )
    else:
        expected_local = (sp.S.One, e**4 * r**2, r**2)
    if any(
        sp.cancel(actual - wanted) != 0
        for actual, wanted in zip(
            (
                coordinate_jacobian,
                pulled_back_jacobian,
                leading_factor,
            ),
            expected_local,
        )
    ):
        raise AssertionError(
            f"unexpected dicritical chart Jacobian in {terminal_case}"
        )

    expected = {
        "Case 1": {
            "dicritical": "R1_1",
            "coefficient": 2,
            "chart": ("e", "r"),
            "coordinate_jacobian": "1/(1+delta*e)",
            "pulled_back_jacobian": (
                "e^2*u(e,r)^2/(1+delta*e)"
            ),
            "generic_leading_factor": "beta^2",
            "affine_meets": (("D0", 1),),
            "boundary_intersection": 35,
        },
        "Case 2": {
            "dicritical": "T2",
            "coefficient": 4,
            "chart": ("e", "r"),
            "coordinate_jacobian": "1",
            "pulled_back_jacobian": "e^4*r^2",
            "generic_leading_factor": "r^2",
            "affine_meets": (("T2", 1),),
            "boundary_intersection": 33,
        },
    }[terminal_case]
    nonzero_affine_intersections = tuple(
        (name, value)
        for name, value in zip(names, affine_intersections)
        if value
    )
    if (
        min(actual_boundary) < 0
        or dicritical != expected["dicritical"]
        or actual_boundary[index] != expected["coefficient"]
        or nonzero_affine_intersections != expected["affine_meets"]
        or boundary_intersection != expected["boundary_intersection"]
        or total_intersection != 35
    ):
        raise AssertionError(
            f"unexpected transformed-Poisson ramification in {terminal_case}"
        )

    return PoissonJacobianRamificationAudit(
        terminal_case=terminal_case,
        jacobian="X^2",
        dicritical_chart=expected["chart"],
        coordinate_jacobian=expected["coordinate_jacobian"],
        pulled_back_jacobian=expected["pulled_back_jacobian"],
        generic_leading_factor=expected["generic_leading_factor"],
        boundary_names=names,
        x_valuations=x_valuations,
        canonical_plus_target=class_ramification,
        actual_boundary_ramification=actual_boundary,
        affine_ramification_component="closure of X=0",
        affine_component_multiplicity=2,
        affine_boundary_intersections=affine_intersections,
        dicritical_divisor=dicritical,
        dicritical_ramification_coefficient=actual_boundary[index],
        dicritical_generic_ramification_index=actual_boundary[index] + 1,
        dicritical_boundary_ramification_intersection=boundary_intersection,
        dicritical_total_ramification_intersection=total_intersection,
        constant_keller_ramification_applicable=False,
    )


def frontier_72_108_poisson_ramification_audits(
) -> tuple[PoissonJacobianRamificationAudit, ...]:
    """Return the actual Case-1 and Case-2 ramification ledgers."""

    return (
        _frontier_72_108_poisson_ramification_audit(
            terminal_case="Case 1",
            package=frontier_72_108_case1_boundary_package_audit(),
        ),
        _frontier_72_108_poisson_ramification_audit(
            terminal_case="Case 2",
            package=frontier_72_108_case2_boundary_package_audit(),
        ),
    )


def frontier_72_108_dicritical_residue_degree_audits(
) -> tuple[DicriticalResidueDegreeAudit, ...]:
    """Enumerate the surviving exact cover/image degree split.

    A polynomial parametrization of coordinate degrees ``(8,12)`` has
    normalization-cover degree dividing ``gcd(8,12)=4``.  Hence the only
    possibilities are ``delta=1,2,4``, with image degrees ``12,6,3`` and
    total normalization differents ``2*delta-2``.

    The homogeneous quartic initial form has cover degree four, but in
    terminal Case 1 the exact residue is the restriction of the alternate
    polynomial pair and may contain lower transverse terms.  In Case 2,
    exact polynomial-decomposition remainder ideals exclude ``delta=2,4``
    without using the J1 compatibility equations or J0.  Its exact residue
    is therefore birational onto a degree-twelve image.
    """

    ramification = {
        audit.terminal_case: audit
        for audit in frontier_72_108_poisson_ramification_audits()
    }
    a_priori = (1, 2, 4)
    results = []
    for terminal_case, residue, excluded in (
        (
            "Case 1",
            "[1:P2(0,r/u):Q2(0,r/u)]",
            (),
        ),
        (
            "Case 2",
            "[1:C(r):G(r)]",
            (2, 4),
        ),
    ):
        covers = tuple(
            degree for degree in a_priori if degree not in excluded
        )
        image_degrees = tuple(12 // degree for degree in covers)
        differents = tuple(2 * degree - 2 for degree in covers)
        ramification_index = ramification[
            terminal_case
        ].dicritical_generic_ramification_index
        contributions = tuple(
            ramification_index * degree for degree in covers
        )
        remaining = tuple(29 - value for value in contributions)
        if any(value <= 0 for value in remaining):
            raise AssertionError("dicritical valuation exceeds degree 29")
        results.append(
            DicriticalResidueDegreeAudit(
                terminal_case=terminal_case,
                residue=residue,
                coordinate_degrees=(8, 12),
                hyperplane_degree=12,
                a_priori_normalization_cover_degrees=a_priori,
                excluded_normalization_cover_degrees=excluded,
                possible_normalization_cover_degrees=covers,
                possible_image_degrees=image_degrees,
                possible_normalization_differents=differents,
                homogeneous_initial_cover_degree=4,
                forced_normalization_cover_degree=(
                    covers[0] if len(covers) == 1 else None
                ),
                forced_image_degree=(
                    image_degrees[0] if len(covers) == 1 else None
                ),
                generic_ramification_index=ramification_index,
                valuation_contributions_to_degree_29=contributions,
                remaining_valuation_degrees=remaining,
            )
        )
    return tuple(results)


def frontier_72_108_case2_residue_stratum_audit(
) -> Case2ResidueStratumExclusionAudit:
    """Return the pinned exact exclusion data for Case-2 covers 2 and 4."""

    return Case2ResidueStratumExclusionAudit(
        a_priori_cover_degrees=(1, 2, 4),
        excluded_cover_degrees=(2, 4),
        surviving_cover_degrees=(1,),
        right_component_normal_forms=(
            "delta=2: h=t^2+a*t, a=C_7/4",
            (
                "delta=4: h=t^4+a*t^3+b*t^2+c*t, "
                "a=C_7/2, b=(C_6-a^2)/2, "
                "c=(C_5-2*a*b)/2"
            ),
        ),
        equation_scope=(
            "J3,J2 and the twelve G coefficients determined by J1; "
            "no J1 compatibility equation and no J0 equation"
        ),
        constraint_counts=((2, 9), (4, 12)),
        constraint_degrees=(
            (2, (14, 10, 6, 25, 21, 17, 13, 9, 5)),
            (
                4,
                (14, 12, 10, 25, 23, 21, 17, 15, 13, 9, 7, 5),
            ),
        ),
        singular_input_sha256=(
            (
                2,
                (
                    "8bd56f701554cfe909f5c5dbdc78139e"
                    "68583043dc1783ce1bffafd48da0297a"
                ),
            ),
            (
                4,
                (
                    "7c4aac7c464441872edc525d8b337bf3"
                    "53f427763cfcbcdcf17c4fd5254bc93a"
                ),
            ),
        ),
        certificate_command=(
            ".venv/bin/python "
            "plane-jc/cas/audit_case2_residue_strata.py"
        ),
        conclusion=(
            "the Case-2 dicritical residue is birational onto its "
            "degree-twelve image"
        ),
    )


def frontier_72_108_case2_maximal_gcd_audit(
) -> Case2MaximalGcdExclusionAudit:
    """Return the pinned exact exclusion data for the maximal gcd row."""

    return Case2MaximalGcdExclusionAudit(
        gcd_degree=7,
        divisibility_condition="C' divides G'",
        total_remainder_coefficients=7,
        selected_remainder_degrees=(0, 1, 2),
        j0_coefficient_degree=19,
        selected_constraint_term_counts=(155, 155, 155, 9),
        selected_constraint_parameter_degrees=(13, 13, 13, 4),
        equation_scope=(
            "three coefficients of remainder(G',C') and the t^19 "
            "coefficient of J0; no residual J1 compatibility equation"
        ),
        singular_input_sha256=(
            "1ac0b4db7ddd0b50fcbef6c93d49c28f"
            "7f80cb4133a73b5f2158af6c78f3b069"
        ),
        excluded_exact_gcd_degrees=(7,),
        surviving_exact_gcd_degrees=(1, 2, 3, 4, 5, 6),
        certificate_command=(
            ".venv/bin/python "
            "plane-jc/cas/audit_case2_maximal_gcd.py"
        ),
        conclusion=(
            "the maximal derivative-gcd stratum is empty; exact gcd "
            "degrees one through six remain"
        ),
    )


def frontier_72_108_case2_gcd_six_audit(
) -> Case2GcdSixExclusionAudit:
    """Return the pinned exact exclusion data for derivative gcd six."""

    return Case2GcdSixExclusionAudit(
        gcd_degree=6,
        linear_cofactor_parameter="v in C'=H*(t+v)",
        selected_factor_conditions=("C'(0)=0", "H(0)=0"),
        total_g_remainder_coefficients=6,
        selected_g_remainder_degrees=(4, 5),
        j0_coefficient_degree=19,
        selected_constraint_term_counts=(5, 30, 656, 352, 9),
        selected_constraint_parameter_degrees=(2, 7, 16, 15, 4),
        equation_scope=(
            "C'(0), H(0), the degree-4 and degree-5 coefficients of "
            "remainder(G',H), and the t^19 coefficient of J0; no "
            "residual J1 compatibility equation"
        ),
        singular_input_sha256=(
            "08bca008eae9dfe29fd383406a1d5b80"
            "a5e75b934a453dff7560cc05cdd30033"
        ),
        excluded_exact_gcd_degrees=(6,),
        surviving_precompatibility_gcd_degrees=(1, 2, 3, 4, 5),
        certificate_command=(
            ".venv/bin/python plane-jc/cas/audit_case2_gcd6.py"
        ),
        conclusion=(
            "before the global J1 endpoint exclusion, the derivative-gcd "
            "degree-six stratum is already empty"
        ),
    )


def frontier_72_108_case2_j1_endpoint_audit(
) -> Case2J1EndpointExclusionAudit:
    """Return the pinned unit certificate on the forced ``G_12 != 0`` open."""

    return Case2J1EndpointExclusionAudit(
        forced_endpoint="G_12",
        endpoint_condition="G_12 != 0 because deg(G)=12",
        compatibility_equation_count=7,
        compatibility_term_counts=(8, 8, 8, 8, 8, 8, 8),
        compatibility_parameter_degrees=(3, 3, 3, 3, 3, 3, 3),
        localization_equation="w*G_12-1",
        equation_scope=(
            "J3,J2, the twelve G coefficients determined by J1, and all "
            "seven residual J1 compatibility equations; no J0 equation"
        ),
        singular_input_sha256=(
            "1618285117e313afb04de0ab935d4a7c"
            "f182f4baabbb0a1346f336f88aaf394e"
        ),
        generic_infinity_characteristic=(4, 13),
        generic_infinity_exceptional_rays=(
            (1, 1),
            (1, 2),
            (1, 3),
            (4, 13),
            (3, 10),
            (2, 7),
            (1, 4),
        ),
        generic_infinity_self_intersections=(
            -2,
            -2,
            -5,
            -1,
            -2,
            -2,
            -2,
        ),
        certificate_command=(
            ".venv/bin/python "
            "plane-jc/cas/case2_infinity_resolution.py"
        ),
        unit_ideal=True,
        conclusion=(
            "the degree-twelve endpoint open is empty after J1 "
            "compatibility, so terminal Case 2 is excluded before J0"
        ),
    )


def frontier_72_108_case2_lower_jet_audit() -> Case2LowerJetAudit:
    """Reduce ``(J0),(J1)`` to coprime derivative data.

    Put ``H=gcd(C',G')``, ``C'=H*c`` and ``G'=H*g`` with
    ``gcd(c,g)=1``.  The bottom equation forces

        B=K*c,  F=K*g.

    Since the Case-2 supports have ``deg(B)<=8``, ``deg(F)<=12`` and exact
    derivative degrees seven and eleven, ``deg(K)<=deg(H)+1``.  Moreover
    ``B(0)=F(0)=0``; coprimality shows that either ``K=0`` or ``t`` divides
    ``K``.  Substitution turns ``(J1)`` into one explicit equation.

    The first coefficients of ``(J4)--(J0)`` give a further triangular
    consequence:

        D_2=1, B_1=2E_2, F_1=F_2=C_1=G_1=G_2=0.

    Thus ``t`` always divides ``H``.  A normalization cover of degree
    ``delta`` also supplies a polynomial right component of degree
    ``delta``, whose derivative divides both ``C'`` and ``G'``.  The combined
    bounds are therefore ``deg(H)>=max(1,delta-1)``.
    """

    t = sp.symbols("t")
    H = sp.Function("H")(t)
    c = sp.Function("c")(t)
    g = sp.Function("g")(t)
    K = sp.Function("K")(t)
    A = sp.Function("A")(t)
    E = sp.Function("E")(t)
    Cprime = H * c
    Gprime = H * g
    B = K * c
    F = K * g
    j0 = sp.expand(B * Gprime - Cprime * F)
    j1 = sp.expand(
        2 * A * Gprime
        + B * sp.diff(F, t)
        - sp.diff(B, t) * F
        - 2 * Cprime * E
    )
    expected_j1 = sp.expand(
        2 * H * (A * g - c * E)
        + K**2 * (c * sp.diff(g, t) - sp.diff(c, t) * g)
    )
    if sp.simplify(j0) != 0 or sp.simplify(j1 - expected_j1) != 0:
        raise AssertionError("Case-2 lower-jet factorization failed")
    wronskian = sp.expand(
        Cprime * sp.diff(Gprime, t)
        - sp.diff(Cprime, t) * Gprime
    )
    relative_wronskian = (
        c * sp.diff(g, t) - sp.diff(c, t) * g
    )
    if sp.simplify(wronskian - H**2 * relative_wronskian) != 0:
        raise AssertionError("Case-2 Wronskian gcd factorization failed")
    if sp.simplify(
        H**2 * expected_j1
        - (
            2 * H**3 * (A * g - c * E)
            + K**2 * wronskian
        )
    ) != 0:
        raise AssertionError("Case-2 J1/Wronskian identity failed")

    a2, d2, b1, b2, e2, c0, c1, c2 = sp.symbols(
        "a2 d2 b1 b2 e2 c0 c1 c2"
    )
    f1, f2, g0, g1, g2 = sp.symbols("f1 f2 g0 g1 g2")
    A0 = t + a2 * t**2
    D0 = d2 * t**2
    B0 = b1 * t + b2 * t**2
    E0 = e2 * t**2
    C0 = c0 + c1 * t + c2 * t**2
    F0 = f1 * t + f2 * t**2
    G0 = g0 + g1 * t + g2 * t**2
    layers = (
        2 * A0 * sp.diff(D0, t) - 3 * sp.diff(A0, t) * D0,
        2 * (A0 * sp.diff(E0, t) - sp.diff(A0, t) * E0)
        + B0 * sp.diff(D0, t)
        - 3 * sp.diff(B0, t) * D0,
        2 * A0 * sp.diff(F0, t)
        - sp.diff(A0, t) * F0
        + B0 * sp.diff(E0, t)
        - 2 * sp.diff(B0, t) * E0
        - 3 * sp.diff(C0, t) * D0,
        2 * A0 * sp.diff(G0, t)
        + B0 * sp.diff(F0, t)
        - sp.diff(B0, t) * F0
        - 2 * sp.diff(C0, t) * E0,
        B0 * sp.diff(G0, t) - sp.diff(C0, t) * F0,
    )
    substitutions = {
        d2: 1,
        b1: 2 * e2,
        f1: 0,
        g1: 0,
        f2: c1,
        g2: 0,
    }
    triangular_checks = (
        (sp.expand(layers[0]).coeff(t, 2), d2),
        (
            sp.expand(layers[1].subs(d2, 1)).coeff(t, 2),
            2 * e2 - b1,
        ),
        (sp.expand(layers[2]).coeff(t, 1), f1),
        (sp.expand(layers[3].subs(f1, 0)).coeff(t, 1), 2 * g1),
        (
            sp.expand(
                layers[2].subs(
                    {d2: 1, b1: 2 * e2, f1: 0}
                )
            ).coeff(t, 2),
            3 * (f2 - c1),
        ),
        (
            sp.expand(layers[3].subs(substitutions)).coeff(t, 2),
            0,
        ),
        (
            sp.expand(layers[4].subs(substitutions)).coeff(t, 2),
            -c1**2,
        ),
    )
    if any(
        sp.expand(actual - expected) != 0
        for actual, expected in triangular_checks
    ):
        raise AssertionError("Case-2 origin-jet triangular audit failed")

    # If deg(H)=1, the forced factor t is all of H.  Coprimality then gives
    # ord(C')=1 and, writing q=ord(g), q>=1.  Put b=ord(K)=ord(B).
    # The leading orders in (J2) are b+q, 2b (absent when b=1), and 3.
    j2_order_pairs = []
    j1_order_pairs = []
    for b_order in range(1, 13):
        for g_cofactor_order in range(1, 11):
            j2_orders = [b_order + g_cofactor_order, 3]
            if b_order != 1:
                j2_orders.append(2 * b_order)
            if j2_orders.count(min(j2_orders)) >= 2:
                pair = (b_order, g_cofactor_order)
                j2_order_pairs.append(pair)
                # In (J1) the leading orders are q+2, 2b+q-1, b+2.
                j1_orders = (
                    g_cofactor_order + 2,
                    2 * b_order + g_cofactor_order - 1,
                    b_order + 2,
                )
                if j1_orders.count(min(j1_orders)) >= 2:
                    j1_order_pairs.append(pair)
    if j2_order_pairs != [(1, 2), (2, 1)] or j1_order_pairs != [(1, 2)]:
        raise AssertionError("Case-2 degree-one gcd order audit failed")

    # On this row the residue net <1,C,G> is a basepoint-free g^2_12.
    # At t=0 the forced orders are 0,2,4; homogenization at infinity gives
    # 0,4,12.  Pluecker's formula has total weight 3*(12-2)=30.
    vanishing_sequences = ((0, 2, 4), (0, 4, 12))
    plucker_weights = tuple(
        sum(order - index for index, order in enumerate(sequence))
        for sequence in vanishing_sequences
    )
    if plucker_weights != (3, 13):
        raise AssertionError("Case-2 Pluecker-weight audit failed")
    remaining_plucker_weight = 3 * (12 - 2) - sum(plucker_weights)
    if remaining_plucker_weight != 14:
        raise AssertionError("Case-2 residual Pluecker weight drifted")
    gcd_wronskian_degrees = tuple(
        (degree, 17 - 2 * degree) for degree in range(1, 8)
    )

    return Case2LowerJetAudit(
        equations=(
            "B*G'-C'*F=0",
            "2*A*G'+B*F'-B'*F-2*C'*E=0",
        ),
        derivative_gcd_decomposition=(
            "H=gcd(C',G')",
            "C'=H*c",
            "G'=H*g",
            "gcd(c,g)=1",
        ),
        j0_conclusion=("B=K*c", "F=K*g"),
        multiplier_degree_bound="deg(K)<=deg(H)+1",
        multiplier_origin_constraint="K=0 or t divides K",
        reduced_j1_equation=(
            "2*H*(A*g-c*E)+K^2*(c*g'-c'*g)=0"
        ),
        wronskian_gcd_factorization=(
            "C'*G''-C''*G'=H^2*(c*g'-c'*g)"
        ),
        reduced_j1_wronskian_identity=(
            "K^2*(C'*G''-C''*G')=-2*H^3*(A*g-c*E)"
        ),
        forced_origin_jets=(
            "D_2=1",
            "B_1=2*E_2",
            "F_1=F_2=0",
            "C_1=0",
            "G_1=G_2=0",
            "t divides H=gcd(C',G')",
        ),
        cover_degree_gcd_bounds=((1, 1), (2, 1), (4, 3)),
        degree_one_gcd_orders=(
            "ord(H)=ord(C')=1",
            "ord(K)=ord(B)=1",
            "ord(E)=2",
            "ord(g)=2",
            "ord(G')=ord(F)=3",
        ),
        degree_one_gcd_excluded_order_pairs=((2, 1),),
        linear_series_degree=12,
        degree_one_vanishing_sequences=(
            ("forced affine point t=0", vanishing_sequences[0]),
            ("point at infinity", vanishing_sequences[1]),
        ),
        degree_one_plucker_weights=(
            ("forced affine point t=0", plucker_weights[0]),
            ("point at infinity", plucker_weights[1]),
        ),
        degree_one_affine_wronskian=(
            "C'*G''-C''*G'=t^3*W_14(t), deg(W_14)=14"
        ),
        degree_one_remaining_plucker_weight=remaining_plucker_weight,
        gcd_degree_residual_wronskian_degrees=(
            gcd_wronskian_degrees
        ),
        unresolved_input=(
            "none globally after the J1 endpoint exclusion; exact gcd "
            "degrees one through five remain only as precompatibility "
            "compact-certificate refinements"
        ),
    )


def frontier_72_108_minimal_dicritical_extension_audit(
) -> MinimalDicriticalExtensionAudit:
    """Enumerate all one-blowup extensions with a zero-pole exceptional.

    Every smooth point of a boundary component and every existing boundary
    crossing is tested.  The new exceptional is assigned pole coefficient
    zero, then the intrinsic ``A2`` and Keller pole-vector gates are rerun.
    A smooth point of ``E3`` is the unique admissible center.

    This is only a numerical counterfactual on the unextended common graph.
    The exact first-block audit shows that the actual source base ideal is at
    ``E3 intersect E4`` and that the low monomial ``X`` excludes a smooth
    ``E3`` basepoint.  The result remains useful as a demonstration that the
    intrinsic numerical gate alone cannot recover residue geometry.
    """

    common_record = frontier_72_108_translation_records()[-1]
    common_graph = compile_log_boundary(
        laurent_translation_graph_certificate(
            common_record,
            chart="GmA1_F4",
        )
    )
    affine_fill = fill_temporary_boundary(
        common_graph,
        component="Xinf",
        expected_unit_rank=0,
    )
    common_poles = frontier_72_108_common_graph_pole_audit()
    smooth_centers = tuple((name,) for name in affine_fill.boundary.names)
    crossing_centers = tuple(
        sorted(
            (
                tuple(sorted(edge))
                for edge in affine_fill.boundary.intersections
            )
        )
    )
    centers = smooth_centers + crossing_centers
    passing: list[
        tuple[
            tuple[str, ...],
            A2BoundaryAudit,
            KellerPoleAudit,
        ]
    ] = []
    for center in centers:
        extended = affine_fill.boundary.blow_up(
            center,
            exceptional="D1",
        )
        intersection = boundary_intersection_matrix(
            extended,
            common_graph.initial_intersection_form,
        )
        boundary_audit = audit_a2_boundary(
            IntrinsicA2Boundary(
                names=extended.names,
                intersection_matrix=intersection,
                genera=(0,) * len(extended.names),
            )
        )
        pole_audit = audit_keller_pole_vector(
            boundary_audit,
            common_poles.pole_vector + (0,),
            require_nonproper=True,
        )
        if pole_audit.passes:
            passing.append((center, boundary_audit, pole_audit))

    if len(passing) != 1 or passing[0][0] != ("E3",):
        raise AssertionError("expected a unique minimal dicritical center E3")
    center, boundary_audit, pole_audit = passing[0]
    if (
        pole_audit.geometric_degree != 426
        or pole_audit.dicritical_candidates != ("D1",)
        or boundary_audit.canonical_coefficients[-1] != 0
        or pole_audit.ramification_coefficients[-1] != 0
        or pole_audit.log_ramification_coefficients[-1] != 1
        or pole_audit.hyperplane_intersections[-1] != 1
    ):
        raise AssertionError("unexpected minimal dicritical extension data")
    return MinimalDicriticalExtensionAudit(
        tested_centers=centers,
        admissible_centers=tuple(item[0] for item in passing),
        unique_center=center,
        pole_vector=pole_audit.pole_vector,
        hyperplane_intersections=pole_audit.hyperplane_intersections,
        canonical_coefficients=tuple(
            int(value) for value in boundary_audit.canonical_coefficients
        ),
        geometric_degree=pole_audit.geometric_degree,
        ordinary_ramification=pole_audit.ramification_coefficients,
        log_ramification=pole_audit.log_ramification_coefficients,
        dicritical_candidates=pole_audit.dicritical_candidates,
    )


def frontier_72_108_two_step_dicritical_witnesses(
) -> tuple[TwoStepDicriticalWitness, ...]:
    """Return exact connected witnesses showing longer extensions are not unique."""

    common_record = frontier_72_108_translation_records()[-1]
    common_graph = compile_log_boundary(
        laurent_translation_graph_certificate(
            common_record,
            chart="GmA1_F4",
        )
    )
    affine_fill = fill_temporary_boundary(
        common_graph,
        component="Xinf",
        expected_unit_rank=0,
    )
    common_poles = frontier_72_108_common_graph_pole_audit()
    specifications = (
        ("Yinf", 12),
        ("E4", 2),
        ("E7", 1),
        ("E8", 1),
    )
    witnesses = []
    for first_center, first_pole in specifications:
        first = affine_fill.boundary.blow_up(
            (first_center,),
            exceptional="D1",
        )
        second = first.blow_up(("D1",), exceptional="D2")
        intersection = boundary_intersection_matrix(
            second,
            common_graph.initial_intersection_form,
        )
        boundary_audit = audit_a2_boundary(
            IntrinsicA2Boundary(
                names=second.names,
                intersection_matrix=intersection,
                genera=(0,) * len(second.names),
            )
        )
        pole_audit = audit_keller_pole_vector(
            boundary_audit,
            common_poles.pole_vector + (first_pole, 0),
            require_nonproper=True,
        )
        if (
            not pole_audit.passes
            or pole_audit.dicritical_candidates != ("D2",)
        ):
            raise AssertionError("two-step dicritical witness failed")
        witnesses.append(
            TwoStepDicriticalWitness(
                first_center=first_center,
                first_pole=first_pole,
                second_center="D1",
                geometric_degree=pole_audit.geometric_degree,
                dicritical_degree=pole_audit.hyperplane_intersections[-1],
                dicritical_canonical_coefficient=int(
                    boundary_audit.canonical_coefficients[-1]
                ),
            )
        )
    expected = (
        ("Yinf", 12, 139, 12, 0),
        ("E4", 2, 323, 2, 0),
        ("E7", 1, 422, 1, 3),
        ("E8", 1, 425, 1, 4),
    )
    if tuple(
        (
            witness.first_center,
            witness.first_pole,
            witness.geometric_degree,
            witness.dicritical_degree,
            witness.dicritical_canonical_coefficient,
        )
        for witness in witnesses
    ) != expected:
        raise AssertionError("unexpected two-step dicritical witness data")
    return tuple(witnesses)


def laurent_translation_graph_certificate(
    record: LaurentTranslationRecord,
    chart: str = "GmA1",
) -> NewtonBoundaryCertificate:
    """Compile the full local source graph resolving one translation map.

    The first cluster follows the order-``q`` equality branch through the
    original ``X0 intersect Yinf`` crossing.  After ``q`` blowups, the
    adapted curve ``t=x^q+lambda*w=0`` separates from ``Yinf`` and meets
    ``E_q`` at a smooth point.  There the residual ideal is ``(t_q,x^q)``,
    represented by the nested scale ``[x^q:t_q]``.  The two clusters have
    lengths ``q`` and ``q`` and principalize the exact base ideal
    ``(t,x^(2q))``.

    This certificate resolves the source base ideal of one isolated
    translation.  It still does not build a compatible target completion or
    identify this graph inside the full Newton case tree.
    """

    q = record.order
    return NewtonBoundaryCertificate(
        name=f"graph of Laurent translation order {q}: {record.name}",
        chart=chart,
        clusters=(
            ToroidalCluster(
                name=f"{record.name}_branch",
                root_boundary="Yinf",
                transverse_boundary="X0",
                local_parameters=("w=1/y", "x"),
                scales=(record.branch_scale,),
                source=record.case_scope,
            ),
            ToroidalCluster(
                name=f"{record.name}_residual_base",
                root_boundary=f"E{q}",
                local_parameters=("x", "t_q"),
                scales=(
                    BranchScale(
                        name=f"{record.name}_residual_x{q}_over_t",
                        u_power=q,
                        v_power=1,
                    ),
                ),
                source="residual ideal (t_q,x^q) at the adapted-curve point",
            ),
        ),
        transformations=(
            f"y -> y + lambda*x^-{q}",
            f"adapt t=x^{q}+lambda*w; principalize (t,x^{2 * q})",
        ),
        theorem_source=record.source,
        exhaustive=True,
    )


def hirzebruch_transition_audit(degree: int) -> HirzebruchTransitionAudit:
    """Verify the final monomial map as an ``F_degree`` chart transition."""

    if degree < 0:
        raise ValueError("the Hirzebruch degree must be nonnegative")
    matrix = sp.Matrix([[-1, 0], [degree, 1]])
    return HirzebruchTransitionAudit(
        degree=degree,
        chart=f"GmA1_F{degree}",
        valuation_matrix=tuple(
            tuple(int(value) for value in row)
            for row in matrix.tolist()
        ),
        boundary_correspondence=(
            ("pre:X0", "post:Xinf"),
            ("pre:Xinf", "post:X0"),
            ("pre:Yinf", "post:Yinf"),
        ),
        determinant=int(matrix.det()),
        involutive=matrix * matrix == sp.eye(2),
    )


def frontier_72_108_incomplete_certificate() -> NewtonBoundaryCertificate:
    """Return the legacy aggregate certificate rejected by the generic IR.

    The explicit frontier audits now construct the source-selected terminal
    architecture in both cases.  ``NewtonBoundaryCertificate`` still accepts
    only named toroidal scale clusters and cannot serialize the nonmonomial
    first-block root cluster as one aggregate input.  This is an IR
    limitation, not missing boundary geometry.
    """

    return NewtonBoundaryCertificate(
        name="published (72,108) Laurent case tree",
        chart="A2",
        corners=(
            (sp.Rational(8), sp.Rational(28)),
            (sp.Rational(11, 4), sp.Rational(7)),
        ),
        clusters=(),
        transformations=(
            "coordinate flip x <-> y",
            "case-dependent y -> y + lambda*x^-2",
            "case-dependent y -> y + lambda*x^-3",
            "common y -> y + alpha*x^-4",
            "final x -> x^-1, y -> x^4*y",
        ),
        theorem_source=frontier_72_108_translation_records()[0].source,
        exhaustive=False,
        missing_data=(
            "generic NewtonBoundaryCertificate serialization of the "
            "explicit nonmonomial first-block and plane-return clusters",
        ),
    )


def frontier_72_108_residue_case_tree() -> dict[str, object]:
    """Return the symbolic factor/residue labels stated in the primary proof."""

    return {
        "pred_1_minus_2": {
            "translation_label": "lambda_1",
            "successors": ["(2,-7)", "(1,-3)"],
        },
        "pred_1_minus_3_one_factor": {
            "translation_label": "lambda",
            "resulting_case": "b",
        },
        "pred_1_minus_3_two_factors": {
            "factor_labels": ["alpha_1", "alpha_2"],
            "constraint": "alpha_1 != alpha_2",
            "selected_translation_label": "alpha_1",
            "unselected_factor": "alpha_2",
            "resulting_case": "c",
            "unselected_continuation": {
                "residue": "beta=alpha_2-alpha_1 != 0",
                "after_common_order_4_numerator": (
                    "x^4*y+alpha-beta*x"
                ),
                "order_4_center_value": "alpha != 0",
                "additional_basepoint": False,
                "after_F4_transition": "X*(Y+alpha)-beta",
                "meets_filled_X0": False,
            },
        },
        "common_order_4": {
            "translation_label": "alpha",
            "applies_after": ["a", "b", "c"],
            "constraint": "alpha != 0",
        },
    }


def frontier_72_108_local_report() -> dict[str, object]:
    """Emit every proved local graph together with the typed global gaps."""

    translations = []
    for record in frontier_72_108_translation_records():
        branch = compile_log_boundary(
            laurent_translation_branch_certificate(record)
        )
        graph = compile_log_boundary(
            laurent_translation_graph_certificate(record)
        )
        base = laurent_translation_base_ideal_audit(record)
        translations.append(
            {
                "name": record.name,
                "order": record.order,
                "case_scope": record.case_scope,
                "common_after_branching": record.common_after_branching,
                "branch_scale": {
                    "ratio": [
                        record.branch_scale.u_power,
                        record.branch_scale.v_power,
                    ],
                    "equality_ray": list(
                        record.branch_scale.equality_ray
                    ),
                    "semigroup_conductor": (
                        record.branch_scale.semigroup_conductor
                    ),
                },
                "base_ideal": {
                    "adapted_parameter": base.adapted_parameter,
                    "original_generators": list(base.original_generators),
                    "adapted_generators": list(base.adapted_generators),
                    "branch_fan_length": base.branch_fan_length,
                    "graph_length": base.base_ideal_length,
                    "target_x_orders": list(base.target_x_orders),
                    "target_yinf_orders": list(
                        base.target_yinf_orders
                    ),
                },
                "branch_fan": branch.as_dict(),
                "source_graph": graph.as_dict(),
            }
        )
    incomplete = frontier_72_108_incomplete_certificate()
    unselected = frontier_72_108_unselected_factor_audit()
    exceptional_poles = frontier_72_108_exceptional_pole_audit()
    common_poles = frontier_72_108_common_graph_pole_audit()
    first_block = frontier_72_108_forced_first_block_cluster_audit()
    case2_package = frontier_72_108_case2_boundary_package_audit()
    case1_package = frontier_72_108_case1_boundary_package_audit()
    plane_return_edge = frontier_72_108_plane_return_edge_audit()
    plane_return_partitions = (
        frontier_72_108_plane_return_partition_audits()
    )
    poisson_ramification = (
        frontier_72_108_poisson_ramification_audits()
    )
    residue_degrees = (
        frontier_72_108_dicritical_residue_degree_audits()
    )
    case2_residue_strata = (
        frontier_72_108_case2_residue_stratum_audit()
    )
    case2_j1_endpoint = frontier_72_108_case2_j1_endpoint_audit()
    case2_maximal_gcd = frontier_72_108_case2_maximal_gcd_audit()
    case2_gcd_six = frontier_72_108_case2_gcd_six_audit()
    case2_lower_jet = frontier_72_108_case2_lower_jet_audit()
    minimal_dicritical = (
        frontier_72_108_minimal_dicritical_extension_audit()
    )
    two_step_witnesses = (
        frontier_72_108_two_step_dicritical_witnesses()
    )
    common_record = frontier_72_108_translation_records()[-1]
    common_graph = compile_log_boundary(
        laurent_translation_graph_certificate(
            common_record,
            chart="GmA1_F4",
        )
    )
    final_transition = hirzebruch_transition_audit(4)
    affine_fill = fill_temporary_boundary(
        common_graph,
        component="Xinf",
        expected_unit_rank=0,
    )
    composition_cases = {
        "a": (2, 4),
        "b_direct": (3, 4),
        "b_after_q2": (2, 3, 4),
        "c_direct": (3, 4),
        "c_after_q2": (2, 3, 4),
    }
    return {
        "name": incomplete.name,
        "status": (
            "source-selected boundary architecture complete in both "
            "terminal cases; transformed-Poisson ramification corrected "
            "and residue cover degrees reduced to 1, 2, or 4"
        ),
        "translations": translations,
        "known_transformations": list(incomplete.transformations),
        "residue_case_tree": frontier_72_108_residue_case_tree(),
        "unselected_factor_continuation": {
            "case_presence": {
                case: present for case, present in unselected.case_presence
            },
            "residue_difference": unselected.residue_difference,
            "after_selected_order_three": (
                unselected.after_selected_order_three
            ),
            "after_common_order_four_numerator": (
                unselected.after_common_order_four_numerator
            ),
            "order_four_center_value": (
                unselected.order_four_center_value
            ),
            "after_hirzebruch_transition": (
                unselected.after_hirzebruch_transition
            ),
            "creates_additional_basepoint": (
                unselected.creates_additional_basepoint
            ),
            "meets_filled_x_zero": unselected.meets_filled_x_zero,
        },
        "exceptional_target_poles": {
            "divisors": list(exceptional_poles.divisors),
            "P": list(exceptional_poles.p_pole_orders),
            "Q": list(exceptional_poles.q_pole_orders),
            "target_infinity": list(
                exceptional_poles.target_infinity_orders
            ),
            "first_fan_R_orders": list(
                exceptional_poles.first_fan_r_orders
            ),
            "residual_R_orders": list(
                exceptional_poles.residual_r_orders
            ),
            "complete_global_pullback": (
                exceptional_poles.complete_global_pullback
            ),
        },
        "common_graph_keller_pole_audit": {
            **common_poles.as_dict(),
            "pretransition_names": list(
                common_poles.boundary_audit.boundary.names
            ),
            "names": [
                "Xinf",
                "Yinf",
                *(f"E{index}" for index in range(1, 9)),
            ],
            "conclusion": (
                "the common graph has no dicritical; any hypothetical "
                "nonproper Keller resolution needs an additional cluster"
            ),
        },
        "forced_first_block_cluster": {
            "crossing": list(first_block.crossing),
            "local_parameters": list(first_block.local_parameters),
            "laurent_coordinates": list(first_block.laurent_coordinates),
            "normalized_sections": list(first_block.normalized_sections),
            "reduced_wronskian": first_block.reduced_wronskian,
            "U_degree": first_block.u_degree,
            "V_degree": first_block.v_degree,
            "root_count": first_block.root_count,
            "first_base_order": first_block.first_base_order,
            "first_exceptional_pole": (
                first_block.first_exceptional_pole
            ),
            "child_base_order": first_block.child_base_order,
            "child_exceptional_pole": (
                first_block.child_exceptional_pole
            ),
            "child_target_degree": first_block.child_target_degree,
            "smooth_E3_basepoints": first_block.smooth_e3_basepoints,
            "pole_audit": first_block.pole_audit.as_dict(),
            "conclusion": (
                "the source-forced cluster has no dicritical and excludes "
                "the numerical free-E3 residue"
            ),
        },
        "terminal_case2_boundary_package": {
            "first_center": list(case2_package.first_center),
            "local_parameters": list(case2_package.local_parameters),
            "first_base_order": case2_package.first_base_order,
            "first_exceptional_pole": (
                case2_package.first_exceptional_pole
            ),
            "second_direction": case2_package.second_direction,
            "second_base_order": case2_package.second_base_order,
            "dicritical_pole": case2_package.dicritical_pole,
            "dicritical_residue": case2_package.dicritical_residue,
            "dicritical_degree": case2_package.dicritical_degree,
            "dicritical_ramification_intersection": (
                case2_package.dicritical_ramification_intersection
            ),
            "pole_audit": case2_package.pole_audit.as_dict(),
            "conclusion": (
                "all source and target clusters forced by terminal Case 2 "
                "pass the intrinsic A2 and divisor-class gates; the "
                "actual X^2 ramification is recorded separately"
            ),
        },
        "plane_return_edge": {
            "local_parameters": list(plane_return_edge.local_parameters),
            "edge_sections": list(plane_return_edge.edge_sections),
            "jacobian_top_layer": (
                plane_return_edge.jacobian_top_layer
            ),
            "weighted_wronskian": plane_return_edge.weighted_wronskian,
            "forced_factorization": list(
                plane_return_edge.forced_factorization
            ),
            "common_factor_degree": (
                plane_return_edge.common_factor_degree
            ),
            "terminal_case2_root_partition": list(
                plane_return_edge.case2_root_partition
            ),
            "edge_only_root_partitions": [
                list(partition)
                for partition in plane_return_edge.case1_root_partitions
            ],
        },
        "terminal_case1_boundary_package": {
            "selected_partition": list(case1_package.selected_partition),
            "final_chart_relation": case1_package.final_chart_relation,
            "reciprocal_relation": case1_package.reciprocal_relation,
            "first_center": list(case1_package.first_center),
            "first_base_order": case1_package.first_base_order,
            "first_exceptional_pole": (
                case1_package.first_exceptional_pole
            ),
            "first_exceptional_coordinate": (
                case1_package.first_exceptional_coordinate
            ),
            "residual_point": case1_package.residual_point,
            "adapted_parameter": case1_package.adapted_parameter,
            "second_base_order": case1_package.second_base_order,
            "dicritical_pole": case1_package.dicritical_pole,
            "dicritical_residue": case1_package.dicritical_residue,
            "dicritical_degree": case1_package.dicritical_degree,
            "dicritical_ramification_intersection": (
                case1_package.dicritical_ramification_intersection
            ),
            "pole_audit": case1_package.pole_audit.as_dict(),
            "conclusion": (
                "selecting the other q=3 factor gives an exact transverse "
                "chart, forces partition [4], and removes the last "
                "boundary-center ambiguity; lower residue coefficients "
                "still control its normalization degree"
            ),
        },
        "transformed_poisson_ramification": [
            {
                "terminal_case": audit.terminal_case,
                "jacobian": audit.jacobian,
                "dicritical_chart": list(audit.dicritical_chart),
                "coordinate_jacobian": audit.coordinate_jacobian,
                "pulled_back_jacobian": audit.pulled_back_jacobian,
                "generic_leading_factor": audit.generic_leading_factor,
                "boundary_names": list(audit.boundary_names),
                "X_valuations": list(audit.x_valuations),
                "canonical_plus_target": list(
                    audit.canonical_plus_target
                ),
                "actual_boundary_ramification": list(
                    audit.actual_boundary_ramification
                ),
                "affine_ramification_component": (
                    audit.affine_ramification_component
                ),
                "affine_component_multiplicity": (
                    audit.affine_component_multiplicity
                ),
                "affine_boundary_intersections": list(
                    audit.affine_boundary_intersections
                ),
                "dicritical_divisor": audit.dicritical_divisor,
                "dicritical_ramification_coefficient": (
                    audit.dicritical_ramification_coefficient
                ),
                "dicritical_generic_ramification_index": (
                    audit.dicritical_generic_ramification_index
                ),
                "dicritical_boundary_ramification_intersection": (
                    audit.dicritical_boundary_ramification_intersection
                ),
                "dicritical_total_ramification_intersection": (
                    audit.dicritical_total_ramification_intersection
                ),
                "constant_keller_ramification_applicable": (
                    audit.constant_keller_ramification_applicable
                ),
            }
            for audit in poisson_ramification
        ],
        "dicritical_residue_degree_splits": [
            {
                "terminal_case": audit.terminal_case,
                "residue": audit.residue,
                "coordinate_degrees": list(audit.coordinate_degrees),
                "hyperplane_degree": audit.hyperplane_degree,
                "a_priori_normalization_cover_degrees": list(
                    audit.a_priori_normalization_cover_degrees
                ),
                "excluded_normalization_cover_degrees": list(
                    audit.excluded_normalization_cover_degrees
                ),
                "possible_normalization_cover_degrees": list(
                    audit.possible_normalization_cover_degrees
                ),
                "possible_image_degrees": list(
                    audit.possible_image_degrees
                ),
                "possible_normalization_differents": list(
                    audit.possible_normalization_differents
                ),
                "homogeneous_initial_cover_degree": (
                    audit.homogeneous_initial_cover_degree
                ),
                "forced_normalization_cover_degree": (
                    audit.forced_normalization_cover_degree
                ),
                "forced_image_degree": audit.forced_image_degree,
                "generic_ramification_index": (
                    audit.generic_ramification_index
                ),
                "valuation_contributions_to_degree_29": list(
                    audit.valuation_contributions_to_degree_29
                ),
                "remaining_valuation_degrees": list(
                    audit.remaining_valuation_degrees
                ),
            }
            for audit in residue_degrees
        ],
        "terminal_case2_residue_stratum_exclusion": {
            "a_priori_cover_degrees": list(
                case2_residue_strata.a_priori_cover_degrees
            ),
            "excluded_cover_degrees": list(
                case2_residue_strata.excluded_cover_degrees
            ),
            "surviving_cover_degrees": list(
                case2_residue_strata.surviving_cover_degrees
            ),
            "right_component_normal_forms": list(
                case2_residue_strata.right_component_normal_forms
            ),
            "equation_scope": case2_residue_strata.equation_scope,
            "constraint_counts": [
                [degree, count]
                for degree, count in (
                    case2_residue_strata.constraint_counts
                )
            ],
            "constraint_degrees": [
                [degree, list(degrees)]
                for degree, degrees in (
                    case2_residue_strata.constraint_degrees
                )
            ],
            "singular_input_sha256": [
                [degree, digest]
                for degree, digest in (
                    case2_residue_strata.singular_input_sha256
                )
            ],
            "certificate_command": (
                case2_residue_strata.certificate_command
            ),
            "conclusion": case2_residue_strata.conclusion,
        },
        "terminal_case2_J1_endpoint_exclusion": {
            "forced_endpoint": case2_j1_endpoint.forced_endpoint,
            "endpoint_condition": case2_j1_endpoint.endpoint_condition,
            "compatibility_equation_count": (
                case2_j1_endpoint.compatibility_equation_count
            ),
            "compatibility_term_counts": list(
                case2_j1_endpoint.compatibility_term_counts
            ),
            "compatibility_parameter_degrees": list(
                case2_j1_endpoint.compatibility_parameter_degrees
            ),
            "localization_equation": (
                case2_j1_endpoint.localization_equation
            ),
            "equation_scope": case2_j1_endpoint.equation_scope,
            "singular_input_sha256": (
                case2_j1_endpoint.singular_input_sha256
            ),
            "generic_infinity_characteristic": list(
                case2_j1_endpoint.generic_infinity_characteristic
            ),
            "generic_infinity_exceptional_rays": [
                list(ray)
                for ray in (
                    case2_j1_endpoint.generic_infinity_exceptional_rays
                )
            ],
            "generic_infinity_self_intersections": list(
                case2_j1_endpoint.generic_infinity_self_intersections
            ),
            "certificate_command": case2_j1_endpoint.certificate_command,
            "unit_ideal": case2_j1_endpoint.unit_ideal,
            "conclusion": case2_j1_endpoint.conclusion,
        },
        "terminal_case2_maximal_gcd_exclusion": {
            "gcd_degree": case2_maximal_gcd.gcd_degree,
            "divisibility_condition": (
                case2_maximal_gcd.divisibility_condition
            ),
            "total_remainder_coefficients": (
                case2_maximal_gcd.total_remainder_coefficients
            ),
            "selected_remainder_degrees": list(
                case2_maximal_gcd.selected_remainder_degrees
            ),
            "j0_coefficient_degree": (
                case2_maximal_gcd.j0_coefficient_degree
            ),
            "selected_constraint_term_counts": list(
                case2_maximal_gcd.selected_constraint_term_counts
            ),
            "selected_constraint_parameter_degrees": list(
                case2_maximal_gcd.selected_constraint_parameter_degrees
            ),
            "equation_scope": case2_maximal_gcd.equation_scope,
            "singular_input_sha256": (
                case2_maximal_gcd.singular_input_sha256
            ),
            "excluded_exact_gcd_degrees": list(
                case2_maximal_gcd.excluded_exact_gcd_degrees
            ),
            "surviving_exact_gcd_degrees": list(
                case2_maximal_gcd.surviving_exact_gcd_degrees
            ),
            "certificate_command": case2_maximal_gcd.certificate_command,
            "conclusion": case2_maximal_gcd.conclusion,
        },
        "terminal_case2_gcd_six_exclusion": {
            "gcd_degree": case2_gcd_six.gcd_degree,
            "linear_cofactor_parameter": (
                case2_gcd_six.linear_cofactor_parameter
            ),
            "selected_factor_conditions": list(
                case2_gcd_six.selected_factor_conditions
            ),
            "total_g_remainder_coefficients": (
                case2_gcd_six.total_g_remainder_coefficients
            ),
            "selected_g_remainder_degrees": list(
                case2_gcd_six.selected_g_remainder_degrees
            ),
            "j0_coefficient_degree": case2_gcd_six.j0_coefficient_degree,
            "selected_constraint_term_counts": list(
                case2_gcd_six.selected_constraint_term_counts
            ),
            "selected_constraint_parameter_degrees": list(
                case2_gcd_six.selected_constraint_parameter_degrees
            ),
            "equation_scope": case2_gcd_six.equation_scope,
            "singular_input_sha256": (
                case2_gcd_six.singular_input_sha256
            ),
            "excluded_exact_gcd_degrees": list(
                case2_gcd_six.excluded_exact_gcd_degrees
            ),
            "surviving_precompatibility_gcd_degrees": list(
                case2_gcd_six.surviving_precompatibility_gcd_degrees
            ),
            "certificate_command": case2_gcd_six.certificate_command,
            "conclusion": case2_gcd_six.conclusion,
        },
        "terminal_case2_lower_jet_reduction": {
            "equations": list(case2_lower_jet.equations),
            "derivative_gcd_decomposition": list(
                case2_lower_jet.derivative_gcd_decomposition
            ),
            "J0_conclusion": list(case2_lower_jet.j0_conclusion),
            "multiplier_degree_bound": (
                case2_lower_jet.multiplier_degree_bound
            ),
            "multiplier_origin_constraint": (
                case2_lower_jet.multiplier_origin_constraint
            ),
            "reduced_J1_equation": (
                case2_lower_jet.reduced_j1_equation
            ),
            "wronskian_gcd_factorization": (
                case2_lower_jet.wronskian_gcd_factorization
            ),
            "reduced_J1_wronskian_identity": (
                case2_lower_jet.reduced_j1_wronskian_identity
            ),
            "forced_origin_jets": list(
                case2_lower_jet.forced_origin_jets
            ),
            "cover_degree_gcd_bounds": [
                list(pair)
                for pair in case2_lower_jet.cover_degree_gcd_bounds
            ],
            "degree_one_gcd_orders": list(
                case2_lower_jet.degree_one_gcd_orders
            ),
            "degree_one_gcd_excluded_order_pairs": [
                list(pair)
                for pair in (
                    case2_lower_jet
                    .degree_one_gcd_excluded_order_pairs
                )
            ],
            "linear_series_degree": (
                case2_lower_jet.linear_series_degree
            ),
            "degree_one_vanishing_sequences": [
                [point, list(sequence)]
                for point, sequence in (
                    case2_lower_jet.degree_one_vanishing_sequences
                )
            ],
            "degree_one_plucker_weights": [
                [point, weight]
                for point, weight in (
                    case2_lower_jet.degree_one_plucker_weights
                )
            ],
            "degree_one_affine_wronskian": (
                case2_lower_jet.degree_one_affine_wronskian
            ),
            "degree_one_remaining_plucker_weight": (
                case2_lower_jet.degree_one_remaining_plucker_weight
            ),
            "gcd_degree_residual_wronskian_degrees": [
                [degree, residual_degree]
                for degree, residual_degree in (
                    case2_lower_jet
                    .gcd_degree_residual_wronskian_degrees
                )
            ],
            "unresolved_input": case2_lower_jet.unresolved_input,
        },
        "plane_return_partition_boundaries": [
            {
                "partition": list(audit.partition),
                "case_scope": list(audit.case_scope),
                "first_center": list(audit.first_center),
                "first_base_order": audit.first_base_order,
                "first_exceptional": audit.first_exceptional,
                "first_exceptional_pole": (
                    audit.first_exceptional_pole
                ),
                "root_chains": [
                    {
                        "root_multiplicity": chain.root_multiplicity,
                        "equality_ray": list(chain.equality_ray),
                        "blowup_rays": [
                            list(ray) for ray in chain.blowup_rays
                        ],
                        "centers": [
                            list(center) for center in chain.centers
                        ],
                        "exceptional_names": list(
                            chain.exceptional_names
                        ),
                        "section_orders": [
                            list(orders)
                            for orders in chain.section_orders
                        ],
                        "target_infinity_poles": list(
                            chain.target_infinity_poles
                        ),
                        "dicritical_divisor": (
                            chain.dicritical_divisor
                        ),
                        "dicritical_residue_exponents": list(
                            chain.dicritical_residue
                        ),
                        "hyperplane_degree": chain.hyperplane_degree,
                        "normalization_cover_degree": (
                            chain.normalization_cover_degree
                        ),
                        "normalization_different": (
                            chain.normalization_different
                        ),
                        "coordinate_projection_differents": list(
                            chain.coordinate_projection_differents
                        ),
                        "source_branch_conductor": (
                            chain.source_branch_conductor
                        ),
                        "image_cusp_conductor": (
                            chain.image_cusp_conductor
                        ),
                    }
                    for chain in audit.root_chains
                ],
                "boundary_names": list(audit.boundary.names),
                "boundary_matrix": [
                    [int(value) for value in row]
                    for row in audit.boundary.class_matrix.tolist()
                ],
                "intersection_matrix": [
                    [int(value) for value in row]
                    for row in (
                        audit.boundary_audit.boundary
                        .intersection_matrix.tolist()
                    )
                ],
                "canonical_coefficients": [
                    int(value)
                    for value in audit.boundary_audit.canonical_coefficients
                ],
                "pole_audit": audit.pole_audit.as_dict(),
            }
            for audit in plane_return_partitions
        ],
        "minimal_dicritical_extension": {
            "tested_center_count": len(
                minimal_dicritical.tested_centers
            ),
            "admissible_centers": [
                list(center)
                for center in minimal_dicritical.admissible_centers
            ],
            "unique_center": list(minimal_dicritical.unique_center),
            "pole_vector": list(minimal_dicritical.pole_vector),
            "hyperplane_intersections": list(
                minimal_dicritical.hyperplane_intersections
            ),
            "canonical_coefficients": list(
                minimal_dicritical.canonical_coefficients
            ),
            "geometric_degree": minimal_dicritical.geometric_degree,
            "ordinary_ramification": list(
                minimal_dicritical.ordinary_ramification
            ),
            "log_ramification": list(
                minimal_dicritical.log_ramification
            ),
            "dicritical_candidates": list(
                minimal_dicritical.dicritical_candidates
            ),
            "scope": (
                "numerical counterfactual on the common graph; the exact "
                "source audit excludes its smooth E3 center"
            ),
        },
        "two_step_dicritical_witnesses": [
            {
                "first_center": witness.first_center,
                "first_pole": witness.first_pole,
                "second_center": witness.second_center,
                "geometric_degree": witness.geometric_degree,
                "dicritical_degree": witness.dicritical_degree,
                "dicritical_canonical_coefficient": (
                    witness.dicritical_canonical_coefficient
                ),
            }
            for witness in two_step_witnesses
        ],
        "translation_compositions": {
            name: {
                "orders": list(audit.orders),
                "leading_order": audit.leading_order,
                "adapted_generators": list(
                    audit.adapted_generators
                ),
                "graph_length": audit.graph_length,
            }
            for name, orders in composition_cases.items()
            for audit in (laurent_translation_composition_audit(orders),)
        },
        "common_translation_graph": common_graph.as_dict(),
        "final_monomial_transition": {
            "chart": final_transition.chart,
            "valuation_matrix": [
                list(row) for row in final_transition.valuation_matrix
            ],
            "boundary_correspondence": [
                list(pair)
                for pair in final_transition.boundary_correspondence
            ],
            "determinant": final_transition.determinant,
            "involutive": final_transition.involutive,
            "additional_blowups": 0,
        },
        "affine_plane_fill": {
            **affine_fill.as_dict(),
            "posttransition_filled_component": "X0",
            "posttransition_boundary_names": [
                "Xinf",
                "Yinf",
                *(f"E{index}" for index in range(1, 9)),
            ],
        },
        "boundary_architecture_missing_data": [],
        "generic_ir_limitations": list(incomplete.missing_data),
    }


def certificate_from_dict(data: dict[str, object]) -> NewtonBoundaryCertificate:
    clusters = []
    for raw_cluster in data.get("clusters", []):
        if not isinstance(raw_cluster, dict):
            raise ValueError("every cluster must be a JSON object")
        raw_scales = raw_cluster.get("scales", [])
        raw_pullbacks = raw_cluster.get("pullbacks", [])
        if not isinstance(raw_scales, list) or not all(
            isinstance(scale, dict) for scale in raw_scales
        ):
            raise ValueError("cluster scales must be JSON objects")
        if not isinstance(raw_pullbacks, list) or not all(
            isinstance(pullback, dict) for pullback in raw_pullbacks
        ):
            raise ValueError("cluster pullbacks must be JSON objects")
        scales = tuple(
            BranchScale(
                name=str(scale["name"]),
                u_power=int(scale["u_power"]),
                v_power=int(scale["v_power"]),
            )
            for scale in raw_scales
        )
        pullbacks = tuple(
            PullbackMonomial(
                name=str(pullback["name"]),
                u_power=int(pullback["u_power"]),
                v_power=int(pullback["v_power"]),
            )
            for pullback in raw_pullbacks
        )
        clusters.append(
            ToroidalCluster(
                name=str(raw_cluster["name"]),
                root_boundary=str(raw_cluster["root_boundary"]),
                scales=scales,
                transverse_boundary=(
                    str(raw_cluster["transverse_boundary"])
                    if raw_cluster.get("transverse_boundary") is not None
                    else None
                ),
                pullbacks=pullbacks,
                local_parameters=tuple(
                    str(value)
                    for value in raw_cluster.get("local_parameters", ("u", "v"))
                ),
                source=str(raw_cluster.get("source", "")),
            )
        )
    corners = tuple(
        (sp.Rational(first), sp.Rational(second))
        for first, second in data.get("corners", [])
    )
    return NewtonBoundaryCertificate(
        name=str(data.get("name", "unnamed Newton boundary")),
        chart=str(data["chart"]),
        clusters=tuple(clusters),
        transformations=tuple(str(item) for item in data.get("transformations", [])),
        exhaustive=bool(data.get("exhaustive", False)),
        corners=corners,
        theorem_source=str(data.get("theorem_source", "")),
        missing_data=tuple(str(item) for item in data.get("missing_data", [])),
    )


def compile_json(path: Path) -> LogBoundaryCompilation:
    data = json.loads(path.read_text())
    if not isinstance(data, dict):
        raise ValueError("the compiler input must be one JSON object")
    return compile_log_boundary(certificate_from_dict(data))


def branch_scale_regression_certificate() -> NewtonBoundaryCertificate:
    """The local normalized blowup of ``(u^3,v^2)`` from the branch-scale audit."""

    return NewtonBoundaryCertificate(
        name="degree-five branch-scale [u^3:v^2]",
        chart="A2",
        corners=(),
        clusters=(
            ToroidalCluster(
                name="collision_0_1",
                root_boundary="L",
                scales=(BranchScale("lambda_0_over_lambda_1", 3, 2),),
                pullbacks=(
                    PullbackMonomial("target_u", 1, 0),
                    PullbackMonomial("target_v", 0, 1),
                ),
                source="branch-scale asymptotics [x^3:y^2]",
            ),
        ),
        transformations=("normalized blowup of the monomial branch-scale ideal",),
        theorem_source="global-comparison-obstruction.tex, Lemma 2.1 and Corollary 2.3",
        exhaustive=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path)
    parser.add_argument(
        "--frontier-72-108",
        action="store_true",
        help="emit the audited isolated Laurent graphs and global gap ledger",
    )
    args = parser.parse_args()
    if args.frontier_72_108:
        if args.input is not None:
            parser.error("--frontier-72-108 does not accept an input file")
        print(
            json.dumps(
                frontier_72_108_local_report(),
                indent=2,
                sort_keys=True,
            )
        )
        return 0
    compilation = (
        compile_json(args.input)
        if args.input is not None
        else compile_log_boundary(branch_scale_regression_certificate())
    )
    print(json.dumps(compilation.as_dict(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
