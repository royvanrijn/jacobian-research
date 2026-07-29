#!/usr/bin/env python3
"""Stage-one compiler for abstract finite-normalization boundary packages.

This module checks exact necessary conditions only.  Passing every check is
reported as ``unknown``, never as existence of a polynomial Keller map.
``realized`` is reserved for a future stage-two certificate containing an
explicit root equation, reconstruction identities, polynomial-ring
isomorphism, and constant-Jacobian verification.

The implementation is dependency-free so that its group, lattice, conductor,
and adjunction certificates can be replayed with ``python3``.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
from fractions import Fraction
from typing import Sequence


Permutation = tuple[int, ...]


class PackageStatus(str, Enum):
    """The compiler's deliberately conservative three-way result."""

    REALIZED = "realized"
    OBSTRUCTED = "obstructed"
    UNKNOWN = "unknown"


class SheetColor(str, Enum):
    """Whether a height-one sheet lies in the distinguished affine open."""

    AFFINE = "affine"
    BOUNDARY = "boundary"


def identity(size: int) -> Permutation:
    return tuple(range(size))


def compose(left: Permutation, right: Permutation) -> Permutation:
    """Return ``left o right``."""

    if len(left) != len(right):
        raise ValueError("permutations must have the same degree")
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation: Permutation) -> Permutation:
    result = [0] * len(permutation)
    for source, target in enumerate(permutation):
        result[target] = source
    return tuple(result)


def permutation_from_cycles(
    degree: int, cycles: Sequence[Sequence[int]]
) -> Permutation:
    """Build a zero-based permutation from one-based disjoint cycles."""

    result = list(range(degree))
    occupied: set[int] = set()
    for cycle in cycles:
        if not cycle:
            raise ValueError("cycles must be nonempty")
        converted = [entry - 1 for entry in cycle]
        if any(entry < 0 or entry >= degree for entry in converted):
            raise ValueError("cycle entry outside the permutation degree")
        if occupied.intersection(converted):
            raise ValueError("cycles must be disjoint")
        occupied.update(converted)
        for source, target in zip(converted, converted[1:] + converted[:1]):
            result[source] = target
    return tuple(result)


def cycle_lengths(permutation: Permutation) -> tuple[int, ...]:
    seen: set[int] = set()
    lengths: list[int] = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return tuple(sorted(lengths))


def permutation_order(permutation: Permutation) -> int:
    def gcd(left: int, right: int) -> int:
        while right:
            left, right = right, left % right
        return abs(left)

    def lcm(left: int, right: int) -> int:
        return abs(left * right) // gcd(left, right)

    result = 1
    for length in cycle_lengths(permutation):
        result = lcm(result, length)
    return result


def generated_group(generators: Sequence[Permutation]) -> frozenset[Permutation]:
    if not generators:
        raise ValueError("at least one monodromy generator is required")
    degree = len(generators[0])
    if any(len(generator) != degree for generator in generators):
        raise ValueError("monodromy generators have inconsistent degrees")
    unit = identity(degree)
    group = {unit}
    frontier = [unit]
    expanded_generators = tuple(generators) + tuple(
        inverse(generator) for generator in generators
    )
    while frontier:
        element = frontier.pop()
        for generator in expanded_generators:
            candidate = compose(generator, element)
            if candidate not in group:
                group.add(candidate)
                frontier.append(candidate)
    return frozenset(group)


def matrix_rank(rows: Sequence[Sequence[int]]) -> int:
    if not rows:
        return 0
    matrix = [[Fraction(value) for value in row] for row in rows]
    column_count = len(matrix[0])
    if any(len(row) != column_count for row in matrix):
        raise ValueError("matrix rows have inconsistent lengths")
    rank = 0
    for column in range(column_count):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [value / pivot_value for value in matrix[rank]]
        for index, row in enumerate(matrix):
            if index == rank or not row[column]:
                continue
            factor = row[column]
            matrix[index] = [
                value - factor * pivot_entry
                for value, pivot_entry in zip(row, matrix[rank])
            ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def determinant(matrix: Sequence[Sequence[int]]) -> int:
    """Exact Bareiss determinant."""

    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("determinant requires a square matrix")
    if size == 0:
        return 1
    work = [list(row) for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        pivot_row = next(
            (
                index
                for index in range(pivot_index, size)
                if work[index][pivot_index] != 0
            ),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != pivot_index:
            work[pivot_index], work[pivot_row] = (
                work[pivot_row],
                work[pivot_index],
            )
            sign *= -1
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                if previous != 1 and numerator % previous:
                    raise ArithmeticError("Bareiss division was not exact")
                work[row][column] = numerator // previous
        previous = pivot
        for row in range(pivot_index + 1, size):
            work[row][pivot_index] = 0
    return sign * work[-1][-1]


@dataclass(frozen=True)
class LocalPrime:
    name: str
    ramification_index: int
    residue_degree: int
    color: SheetColor
    selected: bool = False

    def __post_init__(self) -> None:
        if self.ramification_index <= 0 or self.residue_degree <= 0:
            raise ValueError("ramification and residue degrees must be positive")

    @property
    def different_order(self) -> int:
        return self.ramification_index - 1


@dataclass(frozen=True)
class BranchDivisor:
    name: str
    monodromy: Permutation
    primes: tuple[LocalPrime, ...]

    @property
    def local_degree(self) -> int:
        return sum(
            prime.ramification_index * prime.residue_degree
            for prime in self.primes
        )

    @property
    def tame_ramification_contribution(self) -> int:
        return sum(
            prime.residue_degree * prime.different_order for prime in self.primes
        )


@dataclass(frozen=True)
class CoverDatum:
    degree: int
    group_name: str
    expected_group_order: int
    branches: tuple[BranchDivisor, ...]
    target_genus: int | None = None
    source_genus: int | None = None
    compact_product_one: bool = False
    commuting_pairs: tuple[tuple[str, str], ...] = ()

    def __post_init__(self) -> None:
        if self.degree <= 0:
            raise ValueError("cover degree must be positive")
        if not self.branches:
            raise ValueError("a cover needs branch data")


@dataclass(frozen=True)
class UnitCongruence:
    coefficients: tuple[int, ...]
    modulus: int

    def __post_init__(self) -> None:
        if self.modulus <= 1:
            raise ValueError("a finite conductor character needs modulus > 1")


@dataclass(frozen=True)
class UnitPresentation:
    """Presentation of the normalization's puncture-unit lattice.

    ``linear_relations`` encode non-torsion conductor characters and lower
    the free rank.  ``congruences`` encode torsion conductor characters;
    they preserve free rank but may replace the lattice by a finite-index
    sublattice.
    """

    generators: tuple[str, ...]
    linear_relations: tuple[tuple[int, ...], ...] = ()
    congruences: tuple[UnitCongruence, ...] = ()

    def free_rank(self) -> int:
        return len(self.generators) - matrix_rank(self.linear_relations)

    def congruence_index(self) -> int | None:
        """Index of the congruence kernel when no linear relation is present."""

        if self.linear_relations:
            return None
        if not self.congruences:
            return 1
        rank = len(self.generators)
        if any(len(item.coefficients) != rank for item in self.congruences):
            raise ValueError("unit congruence has the wrong number of coefficients")
        moduli = tuple(item.modulus for item in self.congruences)
        zero = tuple(0 for _ in moduli)
        column_images = tuple(
            tuple(
                congruence.coefficients[column] % congruence.modulus
                for congruence in self.congruences
            )
            for column in range(rank)
        )
        image = {zero}
        frontier = [zero]
        while frontier:
            value = frontier.pop()
            for generator in column_images:
                candidate = tuple(
                    (entry + step) % modulus
                    for entry, step, modulus in zip(value, generator, moduli)
                )
                if candidate not in image:
                    image.add(candidate)
                    frontier.append(candidate)
        return len(image)


@dataclass(frozen=True)
class ConductorAutomorphism:
    name: str
    point_images: tuple[tuple[str, str], ...]

    def image_map(self) -> dict[str, str]:
        return dict(self.point_images)


@dataclass(frozen=True)
class ConductorDatum:
    delta_length: int
    node_pairs: tuple[tuple[str, str], ...] = ()
    automorphisms: tuple[ConductorAutomorphism, ...] = ()

    def __post_init__(self) -> None:
        if self.delta_length < 0:
            raise ValueError("conductor delta length must be nonnegative")


@dataclass(frozen=True)
class AdjunctionDatum:
    intersection_matrix: tuple[tuple[int, ...], ...]
    canonical_class: tuple[int, ...]
    curve_class: tuple[int, ...]

    def lhs(self) -> int:
        size = len(self.intersection_matrix)
        if len(self.canonical_class) != size or len(self.curve_class) != size:
            raise ValueError("adjunction vectors do not match the lattice rank")
        if any(len(row) != size for row in self.intersection_matrix):
            raise ValueError("intersection matrix must be square")
        curve_plus_canonical = tuple(
            curve + canonical
            for curve, canonical in zip(self.curve_class, self.canonical_class)
        )
        return sum(
            self.curve_class[row]
            * self.intersection_matrix[row][column]
            * curve_plus_canonical[column]
            for row in range(size)
            for column in range(size)
        )


@dataclass(frozen=True)
class SelectedCurve:
    name: str
    normalization_genus: int
    puncture_count: int
    conductor: ConductorDatum
    units: UnitPresentation
    adjunction: AdjunctionDatum
    declared_arithmetic_genus: int

    def arithmetic_genus(self) -> int:
        return self.normalization_genus + self.conductor.delta_length


@dataclass(frozen=True)
class DeterminantLedgerRow:
    divisor: str
    source_jacobian_order: int
    controlled_pullback_order: int
    controlled_exponent: int
    target_jacobian_order: int

    @property
    def left_order(self) -> int:
        return (
            self.source_jacobian_order
            + self.controlled_exponent * self.controlled_pullback_order
        )


@dataclass(frozen=True)
class AffineSemigroup:
    """An actual finitely generated value semigroup, not its saturation."""

    generators: tuple[tuple[int, ...], ...]

    @property
    def dimension(self) -> int:
        return len(self.generators[0]) if self.generators else 0

    def validate(self) -> None:
        if not self.generators:
            raise ValueError("an affine semigroup needs generators")
        dimension = self.dimension
        if dimension == 0:
            raise ValueError("semigroup value vectors must be nonempty")
        for generator in self.generators:
            if len(generator) != dimension:
                raise ValueError("semigroup generators have inconsistent dimensions")
            if any(value < 0 for value in generator):
                raise ValueError("this stage-one semigroup must lie in N^r")
            if not any(generator):
                raise ValueError("the zero vector is not a semigroup generator")

    def contains(self, target: tuple[int, ...]) -> bool:
        """Exact nonnegative-combination membership by bounded recursion."""

        self.validate()
        if len(target) != self.dimension or any(value < 0 for value in target):
            return False
        memo: dict[tuple[int, tuple[int, ...]], bool] = {}

        def solve(index: int, remainder: tuple[int, ...]) -> bool:
            key = (index, remainder)
            if key in memo:
                return memo[key]
            if index == len(self.generators):
                result = not any(remainder)
                memo[key] = result
                return result
            generator = self.generators[index]
            bounds = [
                remainder[coordinate] // value
                for coordinate, value in enumerate(generator)
                if value > 0
            ]
            maximum = min(bounds)
            for coefficient in range(maximum + 1):
                next_remainder = tuple(
                    value - coefficient * step
                    for value, step in zip(remainder, generator)
                )
                if solve(index + 1, next_remainder):
                    memo[key] = True
                    return True
            memo[key] = False
            return False

        return solve(0, target)


@dataclass(frozen=True)
class SemigroupRequirement:
    name: str
    value: tuple[int, ...]
    required_member: bool = True


@dataclass(frozen=True)
class ReconstructionCoordinate:
    name: str
    valuation_orders: tuple[int, ...]
    pole_bounds: tuple[int, ...]


@dataclass(frozen=True)
class PolarLedger:
    """Finite pole box and compressed filtered-ring certificates."""

    divisors: tuple[str, ...]
    boundary_divisors: tuple[bool, ...]
    coordinates: tuple[ReconstructionCoordinate, ...]
    semigroup: AffineSemigroup
    semigroup_requirements: tuple[SemigroupRequirement, ...]
    polar_completeness_certificate: str
    initial_presentation_certificate: str
    rees_strictness_certificate: str

    def certificate_complete(self) -> bool:
        return all(
            value.strip()
            for value in (
                self.polar_completeness_certificate,
                self.initial_presentation_certificate,
                self.rees_strictness_certificate,
            )
        )


@dataclass(frozen=True)
class AffineSourceNecessaryData:
    """The boundary-class part of the ``U = A^n`` claim."""

    boundary_class_matrix: tuple[tuple[int, ...], ...]
    expected_unit_rank: int = 0
    expected_class_group_rank: int = 0


@dataclass(frozen=True)
class StageTwoRealizationCertificate:
    """References to the proof-bearing outputs of the symbolic stage.

    These fields name certificates rather than restating their mathematics.
    Stage one checks completeness of this interface but never promotes a
    package to ``realized``.  A future symbolic verifier must replay the
    referenced certificates before emitting that status.
    """

    root_equation: str
    local_factorization_certificate: str
    reconstruction_identities: str
    polynomial_ring_isomorphism: str
    constant_jacobian_certificate: str
    monodromy_certificate: str

    def complete(self) -> bool:
        return all(
            value.strip()
            for value in (
                self.root_equation,
                self.local_factorization_certificate,
                self.reconstruction_identities,
                self.polynomial_ring_isomorphism,
                self.constant_jacobian_certificate,
                self.monodromy_certificate,
            )
        )


@dataclass(frozen=True)
class BoundaryPackage:
    name: str
    cover: CoverDatum
    selected_curves: tuple[SelectedCurve, ...]
    determinant_ledger: tuple[DeterminantLedgerRow, ...]
    polar_ledger: PolarLedger
    affine_source: AffineSourceNecessaryData
    stage_two_realization: StageTwoRealizationCertificate | None = None


@dataclass(frozen=True)
class Diagnostic:
    code: str
    message: str
    obstruction: bool


@dataclass(frozen=True)
class Compilation:
    package: str
    status: PackageStatus
    diagnostics: tuple[Diagnostic, ...]
    invariants: tuple[tuple[str, object], ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "package": self.package,
            "status": self.status.value,
            "diagnostics": [asdict(item) for item in self.diagnostics],
            "invariants": dict(self.invariants),
        }


def _validate_cover(
    cover: CoverDatum,
    diagnostics: list[Diagnostic],
    invariants: dict[str, object],
) -> None:
    degree = cover.degree
    branch_names = {branch.name for branch in cover.branches}
    if len(branch_names) != len(cover.branches):
        diagnostics.append(
            Diagnostic("cover.duplicate_branch", "branch names are not unique", True)
        )

    for branch in cover.branches:
        if len(branch.monodromy) != degree:
            diagnostics.append(
                Diagnostic(
                    "cover.permutation_degree",
                    f"{branch.name} has the wrong permutation degree",
                    True,
                )
            )
            continue
        if sorted(branch.monodromy) != list(range(degree)):
            diagnostics.append(
                Diagnostic(
                    "cover.not_permutation",
                    f"{branch.name} is not a permutation",
                    True,
                )
            )
            continue
        if branch.local_degree != degree:
            diagnostics.append(
                Diagnostic(
                    "cover.local_degree",
                    f"{branch.name} has local degree {branch.local_degree}, not {degree}",
                    True,
                )
            )
        expanded_indices = sorted(
            prime.ramification_index
            for prime in branch.primes
            for _ in range(prime.residue_degree)
        )
        if expanded_indices != list(cycle_lengths(branch.monodromy)):
            diagnostics.append(
                Diagnostic(
                    "cover.cycle_profile",
                    f"{branch.name} cycle lengths do not match its (e,f) profile",
                    True,
                )
            )
        for prime in branch.primes:
            if prime.ramification_index > 1 and prime.color is SheetColor.AFFINE:
                diagnostics.append(
                    Diagnostic(
                        "keller.affine_ramification",
                        f"ramified prime {prime.name} is colored affine",
                        True,
                    )
                )
    invariants["selected_ramified_prime_count"] = sum(
        1
        for branch in cover.branches
        for prime in branch.primes
        if prime.selected and prime.ramification_index > 1
    )

    generators = tuple(branch.monodromy for branch in cover.branches)
    if all(
        len(generator) == degree
        and sorted(generator) == list(range(degree))
        for generator in generators
    ):
        group = generated_group(generators)
        invariants["generated_group_order"] = len(group)
        orbit = {element[0] for element in group}
        invariants["transitive"] = len(orbit) == degree
        if len(orbit) != degree:
            diagnostics.append(
                Diagnostic(
                    "cover.disconnected",
                    "the branch permutations do not act transitively",
                    True,
                )
            )
        if len(group) != cover.expected_group_order:
            diagnostics.append(
                Diagnostic(
                    "cover.group_order",
                    f"generated order {len(group)} does not match "
                    f"declared {cover.group_name} order {cover.expected_group_order}",
                    True,
                )
            )

        if cover.compact_product_one:
            branch_product = identity(degree)
            for generator in generators:
                branch_product = compose(branch_product, generator)
            invariants["branch_product_identity"] = branch_product == identity(degree)
            if branch_product != identity(degree):
                diagnostics.append(
                    Diagnostic(
                        "cover.product_relation",
                        "compact curve branch cycles do not have product one",
                        True,
                    )
                )

        by_name = {branch.name: branch.monodromy for branch in cover.branches}
        for left_name, right_name in cover.commuting_pairs:
            if left_name not in by_name or right_name not in by_name:
                diagnostics.append(
                    Diagnostic(
                        "cover.unknown_commuting_branch",
                        f"unknown commuting pair ({left_name}, {right_name})",
                        True,
                    )
                )
                continue
            left = by_name[left_name]
            right = by_name[right_name]
            if compose(left, right) != compose(right, left):
                diagnostics.append(
                    Diagnostic(
                        "cover.noncommuting_snc_inertia",
                        f"inertia at intersecting divisors {left_name} and "
                        f"{right_name} does not commute",
                        True,
                    )
                )

    if cover.target_genus is not None or cover.source_genus is not None:
        if cover.target_genus is None or cover.source_genus is None:
            diagnostics.append(
                Diagnostic(
                    "cover.partial_genus",
                    "source and target genus must be supplied together",
                    True,
                )
            )
        else:
            ramification = sum(
                branch.tame_ramification_contribution for branch in cover.branches
            )
            lhs = 2 * cover.source_genus - 2
            rhs = degree * (2 * cover.target_genus - 2) + ramification
            invariants["riemann_hurwitz"] = {
                "lhs": lhs,
                "rhs": rhs,
                "ramification": ramification,
            }
            if lhs != rhs:
                diagnostics.append(
                    Diagnostic(
                        "cover.riemann_hurwitz",
                        f"Riemann-Hurwitz gives {rhs}, not declared {lhs}",
                        True,
                    )
                )


def _validate_selected_curve(
    curve: SelectedCurve,
    diagnostics: list[Diagnostic],
    invariants: dict[str, object],
) -> None:
    actual_genus = curve.arithmetic_genus()
    curve_invariants: dict[str, object] = {
        "normalization_genus": curve.normalization_genus,
        "arithmetic_genus": actual_genus,
        "unit_rank": curve.units.free_rank(),
        "unit_congruence_index": curve.units.congruence_index(),
        "adjunction_lhs": curve.adjunction.lhs(),
    }
    invariants[f"curve:{curve.name}"] = curve_invariants

    if actual_genus != curve.declared_arithmetic_genus:
        diagnostics.append(
            Diagnostic(
                "curve.arithmetic_genus",
                f"{curve.name} has arithmetic genus {actual_genus}, not "
                f"declared {curve.declared_arithmetic_genus}",
                True,
            )
        )
    expected_adjunction = 2 * actual_genus - 2
    if curve.adjunction.lhs() != expected_adjunction:
        diagnostics.append(
            Diagnostic(
                "curve.adjunction",
                f"{curve.name} adjunction gives {curve.adjunction.lhs()}, not "
                f"{expected_adjunction}",
                True,
            )
        )

    generator_count = len(curve.units.generators)
    for row in curve.units.linear_relations:
        if len(row) != generator_count:
            diagnostics.append(
                Diagnostic(
                    "curve.unit_relation_width",
                    f"{curve.name} has a malformed unit relation",
                    True,
                )
            )
    for congruence in curve.units.congruences:
        if len(congruence.coefficients) != generator_count:
            diagnostics.append(
                Diagnostic(
                    "curve.unit_congruence_width",
                    f"{curve.name} has a malformed unit congruence",
                    True,
                )
            )

    node_pair_sets = {frozenset(pair) for pair in curve.conductor.node_pairs}
    if curve.conductor.delta_length < len(curve.conductor.node_pairs):
        diagnostics.append(
            Diagnostic(
                "curve.conductor_delta",
                f"{curve.name} has {len(curve.conductor.node_pairs)} node pairs "
                f"but conductor delta length {curve.conductor.delta_length}",
                True,
            )
        )
    for automorphism in curve.conductor.automorphisms:
        images = automorphism.image_map()
        for pair in curve.conductor.node_pairs:
            if any(point not in images for point in pair):
                diagnostics.append(
                    Diagnostic(
                        "curve.partial_conductor_action",
                        f"{automorphism.name} is not defined on conductor pair {pair}",
                        True,
                    )
                )
                continue
            image_pair = frozenset(images[point] for point in pair)
            if image_pair not in node_pair_sets:
                diagnostics.append(
                    Diagnostic(
                        "curve.conductor_not_preserved",
                        f"{automorphism.name} sends conductor pair {pair} to "
                        f"{tuple(sorted(image_pair))}",
                        True,
                    )
                )


def compile_boundary_package(package: BoundaryPackage) -> Compilation:
    diagnostics: list[Diagnostic] = []
    invariants: dict[str, object] = {
        "degree": package.cover.degree,
        "declared_group": package.cover.group_name,
    }
    _validate_cover(package.cover, diagnostics, invariants)

    for curve in package.selected_curves:
        _validate_selected_curve(curve, diagnostics, invariants)

    for row in package.determinant_ledger:
        if row.left_order != row.target_jacobian_order:
            diagnostics.append(
                Diagnostic(
                    "ledger.determinant_balance",
                    f"{row.divisor} has left order {row.left_order} and target "
                    f"order {row.target_jacobian_order}",
                    True,
                )
            )

    polar = package.polar_ledger
    if len(polar.divisors) != len(polar.boundary_divisors):
        diagnostics.append(
            Diagnostic(
                "polar.divisor_colors",
                "polar divisor names and colors have different lengths",
                True,
            )
        )
    try:
        polar.semigroup.validate()
    except ValueError as error:
        diagnostics.append(Diagnostic("semigroup.presentation", str(error), True))
    else:
        for requirement in polar.semigroup_requirements:
            actual = polar.semigroup.contains(requirement.value)
            if actual != requirement.required_member:
                diagnostics.append(
                    Diagnostic(
                        "semigroup.membership",
                        f"{requirement.name} has semigroup membership {actual}, "
                        f"not declared {requirement.required_member}",
                        True,
                    )
                )

    negative_slots = 0
    for coordinate in polar.coordinates:
        if (
            len(coordinate.valuation_orders) != len(polar.divisors)
            or len(coordinate.pole_bounds) != len(polar.divisors)
        ):
            diagnostics.append(
                Diagnostic(
                    "polar.coordinate_width",
                    f"{coordinate.name} does not have one order and bound per divisor",
                    True,
                )
            )
            continue
        for index, (order, bound) in enumerate(
            zip(coordinate.valuation_orders, coordinate.pole_bounds)
        ):
            if bound < 0:
                diagnostics.append(
                    Diagnostic(
                        "polar.negative_bound",
                        f"{coordinate.name} has a negative pole bound",
                        True,
                    )
                )
                continue
            negative_slots += bound
            if order < -bound:
                diagnostics.append(
                    Diagnostic(
                        "polar.bound_exceeded",
                        f"{coordinate.name} has order {order} below bound {-bound} "
                        f"at {polar.divisors[index]}",
                        True,
                    )
                )
            if (
                order < 0
                and index < len(polar.boundary_divisors)
                and not polar.boundary_divisors[index]
            ):
                diagnostics.append(
                    Diagnostic(
                        "polar.affine_pole",
                        f"{coordinate.name} has a pole on affine divisor "
                        f"{polar.divisors[index]}",
                        True,
                    )
                )
    invariants["negative_principal_part_slots"] = negative_slots
    if not polar.certificate_complete():
        diagnostics.append(
            Diagnostic(
                "polar.incomplete_certificate",
                "polar completeness, initial presentation, and Rees strictness "
                "must all be certified",
                False,
            )
        )

    class_matrix = package.affine_source.boundary_class_matrix
    class_determinant = determinant(class_matrix)
    invariants["boundary_class_determinant"] = class_determinant
    if abs(class_determinant) != 1:
        diagnostics.append(
            Diagnostic(
                "affine_source.boundary_class_map",
                f"boundary class determinant is {class_determinant}, not a unit",
                True,
            )
        )
    if package.affine_source.expected_unit_rank != 0:
        diagnostics.append(
            Diagnostic(
                "affine_source.units",
                "an affine-space source cannot have nonconstant units",
                True,
            )
        )
    if package.affine_source.expected_class_group_rank != 0:
        diagnostics.append(
            Diagnostic(
                "affine_source.class_group",
                "an affine-space source cannot have nonzero class-group rank",
                True,
            )
        )

    obstruction = any(item.obstruction for item in diagnostics)
    if obstruction:
        status = PackageStatus.OBSTRUCTED
    else:
        status = PackageStatus.UNKNOWN
        if package.stage_two_realization is not None:
            if package.stage_two_realization.complete():
                diagnostics.append(
                    Diagnostic(
                        "stage_two.unverified",
                        "stage-two certificate references are complete but have "
                        "not been replayed by this stage-one compiler",
                        False,
                    )
                )
            else:
                diagnostics.append(
                    Diagnostic(
                        "stage_two.incomplete",
                        "a partial stage-two realization certificate cannot "
                        "establish existence",
                        False,
                    )
                )
        diagnostics.append(
            Diagnostic(
                "stage_two.missing",
                "all implemented necessary filters pass; an explicit root "
                "equation and affine-space reconstruction are still required",
                False,
            )
        )
    return Compilation(
        package=package.name,
        status=status,
        diagnostics=tuple(diagnostics),
        invariants=tuple(sorted(invariants.items())),
    )


def _primes_from_cycle_profile(
    branch_name: str,
    indices: Sequence[int],
    *,
    ramified_boundary: bool = True,
    select_index: int | None = None,
) -> tuple[LocalPrime, ...]:
    primes = []
    for index, ramification in enumerate(indices):
        color = (
            SheetColor.BOUNDARY
            if ramification > 1 and ramified_boundary
            else SheetColor.AFFINE
        )
        primes.append(
            LocalPrime(
                name=f"{branch_name}_{index + 1}",
                ramification_index=ramification,
                residue_degree=1,
                color=color,
                selected=select_index == index,
            )
        )
    return tuple(primes)


def _three_puncture_curve() -> SelectedCurve:
    return SelectedCurve(
        name="three_punctured_rational",
        normalization_genus=0,
        puncture_count=3,
        conductor=ConductorDatum(delta_length=0),
        units=UnitPresentation(generators=("t", "t-1")),
        adjunction=AdjunctionDatum(
            intersection_matrix=((0, 1), (1, 0)),
            canonical_class=(-2, -2),
            curve_class=(1, 0),
        ),
        declared_arithmetic_genus=0,
    )


def _standard_polar_ledger(divisors: tuple[str, ...]) -> PolarLedger:
    dimension = len(divisors)
    basis = tuple(
        tuple(1 if row == column else 0 for column in range(dimension))
        for row in range(dimension)
    )
    coordinates = tuple(
        ReconstructionCoordinate(
            name=f"x_{index + 1}",
            valuation_orders=tuple(
                -1 if index == column else 0 for column in range(dimension)
            ),
            pole_bounds=tuple(
                1 if index == column else 0 for column in range(dimension)
            ),
        )
        for index in range(dimension)
    )
    return PolarLedger(
        divisors=divisors,
        boundary_divisors=tuple(True for _ in divisors),
        coordinates=coordinates,
        semigroup=AffineSemigroup(generators=basis),
        semigroup_requirements=tuple(
            SemigroupRequirement(f"basis_{index + 1}", vector)
            for index, vector in enumerate(basis)
        ),
        polar_completeness_certificate="declared denominator support",
        initial_presentation_certificate="free coordinate initial algebra",
        rees_strictness_certificate="free saturated Rees module",
    )


def a4_three_puncture_package() -> BoundaryPackage:
    """Degree-four tetrahedral seed with a double-transposition branch."""

    a = permutation_from_cycles(4, ((1, 2), (3, 4)))
    b = permutation_from_cycles(4, ((2, 3, 4),))
    c = inverse(compose(a, b))
    cover = CoverDatum(
        degree=4,
        group_name="A4",
        expected_group_order=12,
        branches=(
            BranchDivisor(
                "double_transposition",
                a,
                (
                    LocalPrime("E2_left", 2, 1, SheetColor.BOUNDARY, True),
                    LocalPrime("E2_right", 2, 1, SheetColor.BOUNDARY, True),
                ),
            ),
            BranchDivisor(
                "three_cycle_one",
                b,
                _primes_from_cycle_profile("B3a", (1, 3)),
            ),
            BranchDivisor(
                "three_cycle_two",
                c,
                _primes_from_cycle_profile("B3b", (1, 3)),
            ),
        ),
        target_genus=0,
        source_genus=0,
        compact_product_one=True,
    )
    return BoundaryPackage(
        name="a4_double_ramification_three_puncture",
        cover=cover,
        selected_curves=(_three_puncture_curve(),),
        determinant_ledger=(
            # The exact pure-target A4 lift: det Phi has W^2 K^3 L and
            # the auxiliary coordinate contributes WL, while B(Phi) has
            # W^3 K^3 L^2.
            DeterminantLedgerRow("W", 2, 1, 1, 3),
            DeterminantLedgerRow("K", 3, 0, 1, 3),
            DeterminantLedgerRow("L", 1, 1, 1, 2),
        ),
        polar_ledger=_standard_polar_ledger(("W", "K", "L")),
        affine_source=AffineSourceNecessaryData(
            boundary_class_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1))
        ),
    )


def gl3f2_triangle_package() -> BoundaryPackage:
    """Degree-seven (2,3,7) Fano action with two ramified involution sheets."""

    a = permutation_from_cycles(7, ((1, 4), (3, 6)))
    b = permutation_from_cycles(7, ((1, 2, 5), (3, 7, 4)))
    c = inverse(compose(a, b))
    cover = CoverDatum(
        degree=7,
        group_name="GL3(F2)",
        expected_group_order=168,
        branches=(
            BranchDivisor(
                "involution",
                a,
                (
                    LocalPrime("I2_left", 2, 1, SheetColor.BOUNDARY, True),
                    LocalPrime("I2_right", 2, 1, SheetColor.BOUNDARY, True),
                    LocalPrime("I1_a", 1, 1, SheetColor.AFFINE),
                    LocalPrime("I1_b", 1, 1, SheetColor.AFFINE),
                    LocalPrime("I1_c", 1, 1, SheetColor.AFFINE),
                ),
            ),
            BranchDivisor(
                "order_three",
                b,
                _primes_from_cycle_profile("G3", (1, 3, 3)),
            ),
            BranchDivisor(
                "order_seven",
                c,
                _primes_from_cycle_profile("G7", (7,)),
            ),
        ),
        target_genus=0,
        source_genus=0,
        compact_product_one=True,
    )
    return BoundaryPackage(
        name="gl3f2_fano_triangle",
        cover=cover,
        selected_curves=(),
        determinant_ledger=(
            DeterminantLedgerRow("involution_different", 2, 0, 1, 2),
            DeterminantLedgerRow("order_three_different", 4, 0, 1, 4),
            DeterminantLedgerRow("order_seven_different", 6, 0, 1, 6),
        ),
        polar_ledger=_standard_polar_ledger(
            ("involution", "order_three", "order_seven")
        ),
        affine_source=AffineSourceNecessaryData(
            boundary_class_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1))
        ),
    )


def nodal_rational_package(*, preserve_conductor: bool = True) -> BoundaryPackage:
    """A4 seed with ``G_m`` normalization and one nontrivial node pairing."""

    package = a4_three_puncture_package()
    automorphism = ConductorAutomorphism(
        name="node_involution",
        point_images=(
            (("p", "q"), ("q", "p"))
            if preserve_conductor
            else (("p", "q"), ("q", "r"))
        ),
    )
    curve = SelectedCurve(
        name="nodal_rational_gm",
        normalization_genus=0,
        puncture_count=2,
        conductor=ConductorDatum(
            delta_length=1,
            node_pairs=(("p", "q"),),
            automorphisms=(automorphism,),
        ),
        units=UnitPresentation(
            generators=("t",),
            congruences=(UnitCongruence((1,), 2),),
        ),
        # A plane cubic has arithmetic genus one, whether smooth elliptic or
        # rational with one ordinary node.
        adjunction=AdjunctionDatum(
            intersection_matrix=((1,),),
            canonical_class=(-3,),
            curve_class=(3,),
        ),
        declared_arithmetic_genus=1,
    )
    return BoundaryPackage(
        name=(
            "nodal_rational_conductor_preserved"
            if preserve_conductor
            else "nodal_rational_conductor_mismatch"
        ),
        cover=package.cover,
        selected_curves=(curve,),
        determinant_ledger=package.determinant_ledger,
        polar_ledger=package.polar_ledger,
        affine_source=package.affine_source,
    )


def elliptic_selected_boundary_package() -> BoundaryPackage:
    """Degree-three genus-one cover with six simple branch values."""

    transpositions = (
        permutation_from_cycles(3, ((1, 2),)),
        permutation_from_cycles(3, ((1, 2),)),
        permutation_from_cycles(3, ((1, 3),)),
        permutation_from_cycles(3, ((1, 3),)),
        permutation_from_cycles(3, ((2, 3),)),
        permutation_from_cycles(3, ((2, 3),)),
    )
    branches = tuple(
        BranchDivisor(
            f"simple_{index + 1}",
            permutation,
            (
                LocalPrime(
                    f"elliptic_ramified_{index + 1}",
                    2,
                    1,
                    SheetColor.BOUNDARY,
                    index == 0,
                ),
                LocalPrime(
                    f"elliptic_unramified_{index + 1}",
                    1,
                    1,
                    SheetColor.AFFINE,
                ),
            ),
        )
        for index, permutation in enumerate(transpositions)
    )
    curve = SelectedCurve(
        name="elliptic_one_puncture",
        normalization_genus=1,
        puncture_count=1,
        conductor=ConductorDatum(delta_length=0),
        # For E minus one point, Div^0 supported at that point is zero.
        units=UnitPresentation(generators=()),
        adjunction=AdjunctionDatum(
            intersection_matrix=((1,),),
            canonical_class=(-3,),
            curve_class=(3,),
        ),
        declared_arithmetic_genus=1,
    )
    return BoundaryPackage(
        name="elliptic_selected_boundary_s3",
        cover=CoverDatum(
            degree=3,
            group_name="S3",
            expected_group_order=6,
            branches=branches,
            target_genus=0,
            source_genus=1,
            compact_product_one=True,
        ),
        selected_curves=(curve,),
        determinant_ledger=tuple(
            DeterminantLedgerRow(branch.name, 1, 0, 1, 1)
            for branch in branches
        ),
        polar_ledger=_standard_polar_ledger(("elliptic_boundary",)),
        affine_source=AffineSourceNecessaryData(
            boundary_class_matrix=((1,),)
        ),
    )


def affine_ramification_obstruction_package() -> BoundaryPackage:
    package = a4_three_puncture_package()
    first_branch = package.cover.branches[0]
    bad_primes = (
        LocalPrime("bad_affine_ramification", 2, 1, SheetColor.AFFINE, True),
        first_branch.primes[1],
    )
    bad_cover = CoverDatum(
        degree=package.cover.degree,
        group_name=package.cover.group_name,
        expected_group_order=package.cover.expected_group_order,
        branches=(
            BranchDivisor(first_branch.name, first_branch.monodromy, bad_primes),
            *package.cover.branches[1:],
        ),
        target_genus=package.cover.target_genus,
        source_genus=package.cover.source_genus,
        compact_product_one=package.cover.compact_product_one,
    )
    return BoundaryPackage(
        name="affine_ramification_obstruction",
        cover=bad_cover,
        selected_curves=package.selected_curves,
        determinant_ledger=package.determinant_ledger,
        polar_ledger=package.polar_ledger,
        affine_source=package.affine_source,
    )


def a4_cone_branch_obstruction_package() -> BoundaryPackage:
    """Exact global branch profile extracted from the symbolic A4 cone.

    The target B-divisor has one prime L with (e,f)=(2,2).  Its two
    geometric inertia cycles agree with a double transposition, but L lies
    in the affine source, so the Keller coloring gate rejects the package.
    """

    package = a4_three_puncture_package()
    first_branch = package.cover.branches[0]
    cone_branch = BranchDivisor(
        first_branch.name,
        first_branch.monodromy,
        (
            LocalPrime(
                "L",
                ramification_index=2,
                residue_degree=2,
                color=SheetColor.AFFINE,
                selected=True,
            ),
        ),
    )
    cone_cover = CoverDatum(
        degree=package.cover.degree,
        group_name=package.cover.group_name,
        expected_group_order=package.cover.expected_group_order,
        branches=(cone_branch, *package.cover.branches[1:]),
        target_genus=package.cover.target_genus,
        source_genus=package.cover.source_genus,
        compact_product_one=package.cover.compact_product_one,
    )
    return BoundaryPackage(
        name="a4_symbolic_cone_branch_obstruction",
        cover=cone_cover,
        selected_curves=package.selected_curves,
        determinant_ledger=package.determinant_ledger,
        polar_ledger=package.polar_ledger,
        affine_source=package.affine_source,
    )


def determinant_ledger_obstruction_package() -> BoundaryPackage:
    package = a4_three_puncture_package()
    return BoundaryPackage(
        name="determinant_ledger_obstruction",
        cover=package.cover,
        selected_curves=package.selected_curves,
        determinant_ledger=(
            DeterminantLedgerRow("W", 2, 1, 1, 4),
            *package.determinant_ledger[1:],
        ),
        polar_ledger=package.polar_ledger,
        affine_source=package.affine_source,
    )


def semigroup_hole_obstruction_package() -> BoundaryPackage:
    """Require a saturation point absent from the certified actual semigroup."""

    package = a4_three_puncture_package()
    polar = PolarLedger(
        divisors=("D1", "D2"),
        boundary_divisors=(True, True),
        coordinates=(
            ReconstructionCoordinate("x", (-1, 0), (1, 0)),
            ReconstructionCoordinate("y", (0, -1), (0, 1)),
        ),
        semigroup=AffineSemigroup(generators=((2, 0), (1, 1), (0, 2))),
        semigroup_requirements=(
            SemigroupRequirement("saturation_hole_(1,0)", (1, 0), True),
        ),
        polar_completeness_certificate="two declared polar divisors",
        initial_presentation_certificate="quadratic semigroup initial algebra",
        rees_strictness_certificate="declared saturated relation module",
    )
    return BoundaryPackage(
        name="actual_semigroup_hole_obstruction",
        cover=package.cover,
        selected_curves=package.selected_curves,
        determinant_ledger=package.determinant_ledger,
        polar_ledger=polar,
        affine_source=package.affine_source,
    )


def all_benchmark_packages() -> tuple[BoundaryPackage, ...]:
    return (
        a4_three_puncture_package(),
        gl3f2_triangle_package(),
        nodal_rational_package(preserve_conductor=True),
        elliptic_selected_boundary_package(),
        nodal_rational_package(preserve_conductor=False),
        affine_ramification_obstruction_package(),
        a4_cone_branch_obstruction_package(),
        determinant_ledger_obstruction_package(),
        semigroup_hole_obstruction_package(),
    )


__all__ = [
    "BoundaryPackage",
    "Compilation",
    "PackageStatus",
    "StageTwoRealizationCertificate",
    "a4_cone_branch_obstruction_package",
    "a4_three_puncture_package",
    "affine_ramification_obstruction_package",
    "all_benchmark_packages",
    "compile_boundary_package",
    "determinant_ledger_obstruction_package",
    "elliptic_selected_boundary_package",
    "gl3f2_triangle_package",
    "nodal_rational_package",
    "semigroup_hole_obstruction_package",
]
