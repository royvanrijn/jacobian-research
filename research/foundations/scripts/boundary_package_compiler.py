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

from dataclasses import asdict, dataclass, replace
from enum import Enum
from fractions import Fraction
from itertools import combinations, product
from math import gcd
from typing import Sequence

from retained_root_euler_gate import (
    RetainedRootEulerDatum,
    RetainedRootEulerStatus,
    audit_retained_root_euler,
)
from conductor_jet_truncation import (
    BoundaryOutputExpressionDatum,
    ConductorBranchJetDatum,
    ConductorBranchSensitivityDatum,
    ConductorJetStatus,
    ContactExpression,
    NormalJetInputDatum,
    audit_conductor_sensitivity_ledger,
    audit_conductor_jet_truncation,
)


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


class TropicalFeasibilityStatus(str, Enum):
    """Outcome of the optional toroidal necessary-condition front end."""

    NOT_DECLARED = "not_declared"
    UNCERTIFIED = "uncertified"
    OBSTRUCTED = "obstructed"
    FEASIBLE = "feasible"


@dataclass(frozen=True)
class TropicalRay:
    """A primitive ray of a source or target valuation fan."""

    name: str
    vector: tuple[int, ...]


@dataclass(frozen=True)
class ToroidalFanDatum:
    """Finite fan skeleton used by the boundary compiler.

    The incidence certificate is theorem-bearing input.  The compiler checks
    ray primitivity, cone ranks, and smoothness, but deliberately does not
    infer a global toroidal embedding from a list of cones.
    """

    lattice_basis: tuple[str, ...]
    rays: tuple[TropicalRay, ...]
    maximal_cones: tuple[tuple[str, ...], ...]
    incidence_certificate: str
    require_smooth_cones: bool = True


@dataclass(frozen=True)
class BoundaryColor:
    """A height-one valuation carried by a ray of the fan.

    Several colors may lie over one ray.  This is how a toroidal extraction
    records distinct residue branches with the same monomial scale.
    """

    name: str
    color: SheetColor
    carrier_ray: str
    residual_label: str = ""


@dataclass(frozen=True)
class ValuationIdentity:
    """A coefficientwise divisor identity in the valuation matrix."""

    name: str
    function_coefficients: tuple[int, ...]
    expected_orders: tuple[int, ...]


@dataclass(frozen=True)
class TropicalRegularityRequirement:
    """Require a named function to have no unlisted boundary poles."""

    function: str
    allowed_pole_colors: tuple[str, ...] = ()


@dataclass(frozen=True)
class UnimodularBlockRequirement:
    """Require a named valuation submatrix to be square unimodular."""

    name: str
    boundary_colors: tuple[str, ...]
    functions: tuple[str, ...]


@dataclass(frozen=True)
class ToroidalAffineScreen:
    """Unit/class necessary conditions for one proposed affine completion."""

    boundary_colors: tuple[str, ...]
    unit_functions: tuple[str, ...]
    normal_core_certificate: str
    require_constant_units: bool = True
    require_trivial_class_group: bool = True
    affine_toric_sufficiency_certificate: str = ""


@dataclass(frozen=True)
class IntegralVariable:
    name: str
    lower_bound: int
    upper_bound: int

    def __post_init__(self) -> None:
        if self.upper_bound < self.lower_bound:
            raise ValueError(
                f"integral variable {self.name!r} has an empty interval"
            )


@dataclass(frozen=True)
class LinearConstraint:
    name: str
    coefficients: tuple[int, ...]
    relation: str
    right_hand_side: int

    def __post_init__(self) -> None:
        if self.relation not in {"==", "<=", ">="}:
            raise ValueError(
                f"linear constraint {self.name!r} has invalid relation "
                f"{self.relation!r}"
            )


@dataclass(frozen=True)
class ValuationFeasibilityProblem:
    """A finite integral problem whose equations include all boundary rows.

    If ``x_i`` is an integral variable and ``c_i`` its row of function
    coefficients, the valuation equation is

        V * (fixed + sum_i x_i*c_i) = target.

    Extra scalar inequalities can encode support-function positivity,
    deletion choices, or bounded ansatz restrictions.  Infeasibility is an
    obstruction only when the input explicitly marks the finite search as
    exhaustive and supplies its scope certificate.
    """

    name: str
    variables: tuple[IntegralVariable, ...]
    fixed_function_coefficients: tuple[int, ...]
    variable_function_coefficients: tuple[tuple[int, ...], ...]
    target_orders: tuple[int, ...]
    constraints: tuple[LinearConstraint, ...] = ()
    infeasibility_is_obstruction: bool = False
    exhaustive_scope_certificate: str = ""


@dataclass(frozen=True)
class ColoredDivisorSpanProblem:
    """Unbounded integral-span gate for a proposed divisor architecture.

    The named generator functions define an integer column lattice in the
    group of colored boundary divisors.  A target outside that lattice cannot
    be realized by any Laurent monomial in those generators, independently
    of exponent bounds or signs.  A conclusive obstruction therefore needs
    an explicit certificate that the named generators exhaust the proposed
    architecture.
    """

    name: str
    generator_functions: tuple[str, ...]
    target_orders: tuple[int, ...]
    infeasibility_is_obstruction: bool = False
    exhaustive_scope_certificate: str = ""


@dataclass(frozen=True)
class ToroidalBoundaryDatum:
    """Colored fan, valuation matrix, and integral affine-completion screens."""

    fan: ToroidalFanDatum
    boundary_colors: tuple[BoundaryColor, ...]
    valuation_functions: tuple[str, ...]
    valuation_matrix: tuple[tuple[int, ...], ...]
    valuation_certificate: str
    identities: tuple[ValuationIdentity, ...] = ()
    regularity_requirements: tuple[TropicalRegularityRequirement, ...] = ()
    unimodular_blocks: tuple[UnimodularBlockRequirement, ...] = ()
    affine_screen: ToroidalAffineScreen | None = None
    divisor_span_problems: tuple[ColoredDivisorSpanProblem, ...] = ()
    feasibility_problems: tuple[ValuationFeasibilityProblem, ...] = ()
    nonlinear_residue: tuple[str, ...] = ()


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


def smith_invariant_factors(
    matrix: Sequence[Sequence[int]],
) -> tuple[int, ...]:
    """Return nonzero Smith factors via gcds of minors.

    The stage-one compiler is dependency-free and its lattices are small.
    For an integer matrix, the gcd of the ``k``-minors is the product of the
    first ``k`` Smith factors, which gives an exact implementation without a
    computer-algebra dependency.
    """

    if not matrix:
        raise ValueError("an integer matrix must have at least one row")
    rows = tuple(tuple(int(value) for value in row) for row in matrix)
    column_count = len(rows[0])
    if any(len(row) != column_count for row in rows):
        raise ValueError("matrix rows have inconsistent lengths")
    if column_count == 0:
        return ()
    rank = matrix_rank(rows)
    determinantal_divisors = [1]
    for size in range(1, rank + 1):
        divisor = 0
        for row_indices in combinations(range(len(rows)), size):
            for column_indices in combinations(range(column_count), size):
                minor = tuple(
                    tuple(rows[row][column] for column in column_indices)
                    for row in row_indices
                )
                divisor = gcd(divisor, abs(determinant(minor)))
        if divisor == 0:
            raise ArithmeticError("nonzero determinantal divisor was not found")
        determinantal_divisors.append(divisor)
    return tuple(
        determinantal_divisors[index] // determinantal_divisors[index - 1]
        for index in range(1, len(determinantal_divisors))
    )


def factorial_core_invariants(
    valuation_matrix: Sequence[Sequence[int]],
) -> dict[str, object]:
    """Compute the exact unit and Weil-class profile of a factorial core."""

    if not valuation_matrix or not valuation_matrix[0]:
        raise ValueError("a factorial-core matrix must have at least one entry")
    rows = tuple(tuple(int(value) for value in row) for row in valuation_matrix)
    character_rank = len(rows[0])
    if any(len(row) != character_rank for row in rows):
        raise ValueError("factorial-core valuation rows have inconsistent lengths")
    rank = matrix_rank(rows)
    smith = smith_invariant_factors(rows)
    return {
        "boundary_count": len(rows),
        "character_rank": character_rank,
        "matrix_rank": rank,
        "unit_rank": character_rank - rank,
        "class_group_free_rank": len(rows) - rank,
        "class_group_torsion": tuple(value for value in smith if value > 1),
        "smith_diagonal": smith,
    }


torus_core_invariants = factorial_core_invariants


def core_class_order(
    valuation_matrix: Sequence[Sequence[int]],
    class_vector: Sequence[int],
) -> int | None:
    """Return a represented class's exact order, or ``None`` if infinite."""

    rows = tuple(tuple(int(value) for value in row) for row in valuation_matrix)
    if not rows:
        raise ValueError("a core presentation matrix must have at least one row")
    if len(class_vector) != len(rows):
        raise ValueError("a class vector needs one entry per boundary prime")
    vector = tuple(int(value) for value in class_vector)
    augmented = tuple(row + (vector[index],) for index, row in enumerate(rows))
    if matrix_rank(augmented) > matrix_rank(rows):
        return None
    base_divisor = 1
    for value in smith_invariant_factors(rows):
        base_divisor *= value
    enlarged_divisor = 1
    for value in smith_invariant_factors(augmented):
        enlarged_divisor *= value
    if enlarged_divisor == 0 or base_divisor % enlarged_divisor:
        raise ArithmeticError("determinantal-divisor quotient is not integral")
    return base_divisor // enlarged_divisor


def presented_core_invariants(
    unit_valuation_matrix: Sequence[Sequence[int]],
    core_relation_matrix: Sequence[Sequence[int]],
    relation_boundary_corrections: Sequence[Sequence[int]],
) -> dict[str, object]:
    """Compile ``Cl(U)`` from a presented core class group and its lifts.

    If ``R`` presents ``Cl(W)`` and ``A`` records the boundary part of chosen
    rational witnesses for the lifted relations, then ``Cl(U)`` is presented
    by the block matrix ``[[V,A],[0,R]]``.  The unit group of ``U`` remains
    the kernel of ``V``.
    """

    valuation = tuple(
        tuple(int(value) for value in row) for row in unit_valuation_matrix
    )
    relations = tuple(
        tuple(int(value) for value in row) for row in core_relation_matrix
    )
    corrections = tuple(
        tuple(int(value) for value in row)
        for row in relation_boundary_corrections
    )
    if not valuation:
        raise ValueError("a presented core needs at least one boundary row")
    if not relations:
        raise ValueError("a presented core needs at least one class-generator row")
    boundary_count = len(valuation)
    unit_count = len(valuation[0])
    core_generator_count = len(relations)
    relation_count = len(relations[0])
    if any(len(row) != unit_count for row in valuation):
        raise ValueError("presented-core valuation rows have inconsistent widths")
    if any(len(row) != relation_count for row in relations):
        raise ValueError("presented-core relation rows have inconsistent widths")
    if len(corrections) != boundary_count or any(
        len(row) != relation_count for row in corrections
    ):
        raise ValueError(
            "relation-boundary corrections must have shape boundary x relations"
        )

    block = tuple(
        valuation[index] + corrections[index]
        for index in range(boundary_count)
    ) + tuple(
        (0,) * unit_count + relations[index]
        for index in range(core_generator_count)
    )
    class_rank = matrix_rank(block)
    smith = smith_invariant_factors(block)
    return {
        "boundary_count": boundary_count,
        "unit_generator_count": unit_count,
        "core_class_generator_count": core_generator_count,
        "core_relation_count": relation_count,
        "unit_rank": unit_count - matrix_rank(valuation),
        "class_group_free_rank": (
            boundary_count + core_generator_count - class_rank
        ),
        "class_group_torsion": tuple(value for value in smith if value > 1),
        "smith_diagonal": smith,
        "presentation_matrix": block,
    }


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
    factorial_core: "FactorialCoreDatum | None" = None
    presented_core: "PresentedCoreDatum | None" = None


@dataclass(frozen=True)
class CoreClassQuery:
    """A boundary vector representing a rank-one reflexive class."""

    name: str
    vector: tuple[int, ...]
    required_trivial: bool = False


@dataclass(frozen=True)
class FactorialCoreDatum:
    """Certified localization data for a dense affine UFD open in ``U``."""

    boundary_primes: tuple[str, ...]
    unit_generators: tuple[str, ...]
    valuation_matrix: tuple[tuple[int, ...], ...]
    normality_certificate: str
    factoriality_certificate: str
    complete_boundary_certificate: str
    class_queries: tuple[CoreClassQuery, ...] = ()

    def scope_complete(self) -> bool:
        return bool(
            self.normality_certificate.strip()
            and self.factoriality_certificate.strip()
            and self.complete_boundary_certificate.strip()
        )


@dataclass(frozen=True)
class PresentedCoreDatum:
    """A core with finitely presented ``Cl(W)`` and lifted relations."""

    boundary_primes: tuple[str, ...]
    unit_generators: tuple[str, ...]
    unit_valuation_matrix: tuple[tuple[int, ...], ...]
    core_class_generators: tuple[str, ...]
    core_relations: tuple[str, ...]
    core_relation_matrix: tuple[tuple[int, ...], ...]
    relation_boundary_corrections: tuple[tuple[int, ...], ...]
    normality_certificate: str
    core_class_presentation_certificate: str
    relation_lift_certificate: str
    complete_boundary_certificate: str
    class_queries: tuple[CoreClassQuery, ...] = ()

    def scope_complete(self) -> bool:
        return all(
            value.strip()
            for value in (
                self.normality_certificate,
                self.core_class_presentation_certificate,
                self.relation_lift_certificate,
                self.complete_boundary_certificate,
            )
        )


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
    retained_root_euler: RetainedRootEulerDatum | None = None
    conductor_jet_branches: tuple[ConductorBranchJetDatum, ...] | None = None
    conductor_jet_sensitivity: (
        tuple[ConductorBranchSensitivityDatum, ...] | None
    ) = None
    toroidal_boundary: ToroidalBoundaryDatum | None = None
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


@dataclass(frozen=True)
class ToroidalAudit:
    status: TropicalFeasibilityStatus
    diagnostics: tuple[Diagnostic, ...]
    invariants: tuple[tuple[str, object], ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "status": self.status.value,
            "diagnostics": [asdict(item) for item in self.diagnostics],
            "invariants": dict(self.invariants),
        }


def _submatrix(
    matrix: Sequence[Sequence[int]],
    row_indices: Sequence[int],
    column_indices: Sequence[int],
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(int(matrix[row][column]) for column in column_indices)
        for row in row_indices
    )


def _primitive_row_with_scale(
    row: Sequence[int],
) -> tuple[tuple[int, ...], int] | None:
    """Return ``(primitive_row, signed_scale)`` for a nonzero integer row."""

    values = tuple(int(value) for value in row)
    divisor = gcd(*(abs(value) for value in values))
    if divisor == 0:
        return None
    primitive = tuple(value // divisor for value in values)
    first = next(value for value in primitive if value)
    if first < 0:
        primitive = tuple(-value for value in primitive)
        divisor = -divisor
    return primitive, divisor


def colored_proportionality_witnesses(
    generator_matrix: Sequence[Sequence[int]],
    target_orders: Sequence[int],
    row_names: Sequence[str] | None = None,
) -> tuple[dict[str, object], ...]:
    """Return exact row relations violated by a colored target divisor.

    If rows ``A_i=a*p`` and ``A_j=b*p`` are proportional, every divisor in
    the column span of ``A`` satisfies ``b*t_i=a*t_j``.  A nonzero mismatch
    is therefore an obstruction over the rationals, hence also over the
    integers and every nonnegative semigroup.  One witness is returned per
    proportionality class; zero generator rows with a nonzero target are
    recorded separately.
    """

    rows = tuple(tuple(int(value) for value in row) for row in generator_matrix)
    targets = tuple(int(value) for value in target_orders)
    if len(rows) != len(targets):
        raise ValueError("a colored target needs one order per generator row")
    if not rows:
        raise ValueError("a colored divisor-span matrix needs at least one row")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("colored divisor-span rows have inconsistent widths")
    names = (
        tuple(row_names)
        if row_names is not None
        else tuple(f"row_{index}" for index in range(len(rows)))
    )
    if len(names) != len(rows):
        raise ValueError("colored divisor-span row names have the wrong length")

    decompositions = tuple(_primitive_row_with_scale(row) for row in rows)
    witnesses: list[dict[str, object]] = []
    for index, decomposition in enumerate(decompositions):
        if decomposition is None and targets[index] != 0:
            witnesses.append(
                {
                    "kind": "zero_generator_row",
                    "colors": (names[index],),
                    "generator_rows": (rows[index],),
                    "target_orders": (targets[index],),
                    "mismatch": targets[index],
                }
            )
    proportional_classes: dict[tuple[int, ...], list[int]] = {}
    for index, decomposition in enumerate(decompositions):
        if decomposition is None:
            continue
        proportional_classes.setdefault(decomposition[0], []).append(index)
    for indices in proportional_classes.values():
        witness_found = False
        for left_offset, left in enumerate(indices):
            _left_primitive, left_scale = decompositions[left]  # type: ignore[misc]
            for right in indices[left_offset + 1 :]:
                _right_primitive, right_scale = decompositions[right]  # type: ignore[misc]
                mismatch = (
                    right_scale * targets[left]
                    - left_scale * targets[right]
                )
                if mismatch == 0:
                    continue
                witnesses.append(
                    {
                        "kind": "proportional_generator_rows",
                        "colors": (names[left], names[right]),
                        "generator_rows": (rows[left], rows[right]),
                        "row_scales": (left_scale, right_scale),
                        "target_orders": (targets[left], targets[right]),
                        "mismatch": mismatch,
                    }
                )
                witness_found = True
                break
            if witness_found:
                break
    return tuple(witnesses)


def _constraint_holds(
    constraint: LinearConstraint, values: tuple[int, ...]
) -> bool:
    left = sum(
        coefficient * value
        for coefficient, value in zip(constraint.coefficients, values)
    )
    if constraint.relation == "==":
        return left == constraint.right_hand_side
    if constraint.relation == "<=":
        return left <= constraint.right_hand_side
    return left >= constraint.right_hand_side


def _pareto_minimal_models(
    models: Sequence[tuple[int, ...]],
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        model
        for model in models
        if not any(
            other != model
            and all(left <= right for left, right in zip(other, model))
            for other in models
        )
    )


def audit_toroidal_boundary(
    datum: ToroidalBoundaryDatum | None,
) -> ToroidalAudit:
    """Compile the colored fan and valuation matrix before nonlinear algebra.

    A feasible audit is only a tropical necessary-condition result.  It does
    not certify a global toroidal embedding, affineness outside the separately
    certified affine-toric case, polynomial adjugate division, or a Keller
    realization.
    """

    if datum is None:
        return ToroidalAudit(
            status=TropicalFeasibilityStatus.NOT_DECLARED,
            diagnostics=(),
            invariants=(),
        )

    diagnostics: list[Diagnostic] = []
    invariants: dict[str, object] = {}
    uncertified = False
    fan = datum.fan
    lattice_rank = len(fan.lattice_basis)
    ray_names = tuple(ray.name for ray in fan.rays)
    ray_index = {name: index for index, name in enumerate(ray_names)}
    color_names = tuple(color.name for color in datum.boundary_colors)
    color_index = {name: index for index, name in enumerate(color_names)}
    function_names = datum.valuation_functions
    function_index = {
        name: index for index, name in enumerate(function_names)
    }

    malformed = False
    if lattice_rank == 0 or len(set(fan.lattice_basis)) != lattice_rank:
        diagnostics.append(
            Diagnostic(
                "toroidal.fan_lattice",
                "the fan needs a nonempty named lattice basis",
                True,
            )
        )
        malformed = True
    if not fan.rays or len(ray_index) != len(ray_names):
        diagnostics.append(
            Diagnostic(
                "toroidal.fan_ray_names",
                "fan ray names must be nonempty and distinct",
                True,
            )
        )
        malformed = True

    ray_profiles: dict[str, object] = {}
    for ray in fan.rays:
        primitive = bool(ray.vector) and gcd(
            *(abs(value) for value in ray.vector)
        ) == 1
        well_formed = (
            len(ray.vector) == lattice_rank
            and any(ray.vector)
            and primitive
        )
        ray_profiles[ray.name] = {
            "vector": ray.vector,
            "primitive": primitive,
        }
        if not well_formed:
            diagnostics.append(
                Diagnostic(
                    "toroidal.fan_ray",
                    f"ray {ray.name} is not a primitive nonzero vector of "
                    f"rank {lattice_rank}",
                    True,
                )
            )
            malformed = True

    cone_profiles: list[dict[str, object]] = []
    used_rays: set[str] = set()
    for cone in fan.maximal_cones:
        if not cone or len(set(cone)) != len(cone):
            diagnostics.append(
                Diagnostic(
                    "toroidal.fan_cone",
                    f"fan cone {cone} is empty or repeats a ray",
                    True,
                )
            )
            malformed = True
            continue
        unknown = tuple(name for name in cone if name not in ray_index)
        if unknown:
            diagnostics.append(
                Diagnostic(
                    "toroidal.fan_cone_ray",
                    f"fan cone {cone} contains unknown rays {unknown}",
                    True,
                )
            )
            malformed = True
            continue
        used_rays.update(cone)
        vectors = tuple(fan.rays[ray_index[name]].vector for name in cone)
        rank = matrix_rank(vectors)
        smith = smith_invariant_factors(vectors) if rank else ()
        smooth = rank == len(vectors) and all(value == 1 for value in smith)
        cone_profiles.append(
            {
                "rays": cone,
                "rank": rank,
                "smith_diagonal": smith,
                "smooth": smooth,
            }
        )
        if rank != len(vectors):
            diagnostics.append(
                Diagnostic(
                    "toroidal.fan_cone_rank",
                    f"fan cone {cone} is not simplicial of its declared size",
                    True,
                )
            )
        elif fan.require_smooth_cones and not smooth:
            diagnostics.append(
                Diagnostic(
                    "toroidal.fan_cone_singular",
                    f"fan cone {cone} has Smith diagonal {smith}",
                    True,
                )
            )
    unused_rays = tuple(name for name in ray_names if name not in used_rays)
    if unused_rays:
        diagnostics.append(
            Diagnostic(
                "toroidal.fan_unused_ray",
                f"fan rays {unused_rays} occur in no maximal cone",
                True,
            )
        )
    if not fan.incidence_certificate.strip():
        diagnostics.append(
            Diagnostic(
                "toroidal.fan_uncertified",
                "the cone list needs a certificate that it is the relevant "
                "toroidal fan or cone complex",
                False,
            )
        )
        uncertified = True
    invariants["fan"] = {
        "lattice_basis": fan.lattice_basis,
        "rays": ray_profiles,
        "maximal_cones": tuple(cone_profiles),
        "incidence_certified": bool(fan.incidence_certificate.strip()),
    }

    if not datum.boundary_colors or len(color_index) != len(color_names):
        diagnostics.append(
            Diagnostic(
                "toroidal.boundary_colors",
                "boundary color names must be nonempty and distinct",
                True,
            )
        )
        malformed = True
    for color in datum.boundary_colors:
        if color.carrier_ray not in ray_index:
            diagnostics.append(
                Diagnostic(
                    "toroidal.boundary_carrier",
                    f"boundary color {color.name} has unknown carrier ray "
                    f"{color.carrier_ray}",
                    True,
                )
            )
            malformed = True
    invariants["boundary_colors"] = tuple(
        {
            "name": color.name,
            "color": color.color.value,
            "carrier_ray": color.carrier_ray,
            "residual_label": color.residual_label,
        }
        for color in datum.boundary_colors
    )

    if not function_names or len(function_index) != len(function_names):
        diagnostics.append(
            Diagnostic(
                "toroidal.valuation_functions",
                "valuation function names must be nonempty and distinct",
                True,
            )
        )
        malformed = True
    if len(datum.valuation_matrix) != len(color_names) or any(
        len(row) != len(function_names) for row in datum.valuation_matrix
    ):
        diagnostics.append(
            Diagnostic(
                "toroidal.valuation_shape",
                "the valuation matrix must have boundary-color rows and "
                "named-function columns",
                True,
            )
        )
        malformed = True
    if not datum.valuation_certificate.strip():
        diagnostics.append(
            Diagnostic(
                "toroidal.valuation_uncertified",
                "the valuation matrix needs exact divisor-order certificates",
                False,
            )
        )
        uncertified = True
    invariants["valuation_matrix"] = {
        "rows": color_names,
        "columns": function_names,
        "matrix": datum.valuation_matrix,
        "certified": bool(datum.valuation_certificate.strip()),
    }

    identity_profiles: dict[str, object] = {}
    if not malformed:
        for identity_datum in datum.identities:
            if (
                len(identity_datum.function_coefficients)
                != len(function_names)
                or len(identity_datum.expected_orders) != len(color_names)
            ):
                diagnostics.append(
                    Diagnostic(
                        "toroidal.valuation_identity_shape",
                        f"identity {identity_datum.name} has incompatible width",
                        True,
                    )
                )
                continue
            actual = tuple(
                sum(
                    order * coefficient
                    for order, coefficient in zip(
                        row, identity_datum.function_coefficients
                    )
                )
                for row in datum.valuation_matrix
            )
            balanced = actual == identity_datum.expected_orders
            identity_profiles[identity_datum.name] = {
                "actual_orders": actual,
                "expected_orders": identity_datum.expected_orders,
                "balanced": balanced,
            }
            if not balanced:
                diagnostics.append(
                    Diagnostic(
                        "toroidal.valuation_identity",
                        f"identity {identity_datum.name} has orders {actual}, "
                        f"not {identity_datum.expected_orders}",
                        True,
                    )
                )
    invariants["valuation_identities"] = identity_profiles

    if not malformed:
        for requirement in datum.regularity_requirements:
            if requirement.function not in function_index:
                diagnostics.append(
                    Diagnostic(
                        "toroidal.regularity_function",
                        f"regularity requirement names unknown function "
                        f"{requirement.function}",
                        True,
                    )
                )
                continue
            unknown = tuple(
                name
                for name in requirement.allowed_pole_colors
                if name not in color_index
            )
            if unknown:
                diagnostics.append(
                    Diagnostic(
                        "toroidal.regularity_color",
                        f"regularity requirement for {requirement.function} "
                        f"names unknown colors {unknown}",
                        True,
                    )
                )
                continue
            column = function_index[requirement.function]
            allowed = set(requirement.allowed_pole_colors)
            for row, boundary in enumerate(datum.boundary_colors):
                order = datum.valuation_matrix[row][column]
                if order >= 0:
                    continue
                if boundary.color is SheetColor.AFFINE:
                    diagnostics.append(
                        Diagnostic(
                            "toroidal.affine_pole",
                            f"{requirement.function} has order {order} on "
                            f"affine color {boundary.name}",
                            True,
                        )
                    )
                elif boundary.name not in allowed:
                    diagnostics.append(
                        Diagnostic(
                            "toroidal.unlisted_pole",
                            f"{requirement.function} has unlisted pole order "
                            f"{order} on {boundary.name}",
                            True,
                        )
                    )

    block_profiles: dict[str, object] = {}
    if not malformed:
        for block in datum.unimodular_blocks:
            unknown_rows = tuple(
                name for name in block.boundary_colors if name not in color_index
            )
            unknown_columns = tuple(
                name for name in block.functions if name not in function_index
            )
            if unknown_rows or unknown_columns:
                diagnostics.append(
                    Diagnostic(
                        "toroidal.unimodular_block_names",
                        f"block {block.name} has unknown rows {unknown_rows} "
                        f"or columns {unknown_columns}",
                        True,
                    )
                )
                continue
            rows = tuple(color_index[name] for name in block.boundary_colors)
            columns = tuple(function_index[name] for name in block.functions)
            matrix = _submatrix(datum.valuation_matrix, rows, columns)
            square = len(rows) == len(columns)
            block_determinant = determinant(matrix) if square else None
            unimodular = square and abs(block_determinant) == 1
            block_profiles[block.name] = {
                "matrix": matrix,
                "determinant": block_determinant,
                "unimodular": unimodular,
            }
            if not unimodular:
                diagnostics.append(
                    Diagnostic(
                        "toroidal.unimodular_block",
                        f"block {block.name} has determinant "
                        f"{block_determinant}, not a unit",
                        True,
                    )
                )
    invariants["unimodular_blocks"] = block_profiles

    if datum.affine_screen is not None and not malformed:
        screen = datum.affine_screen
        unknown_rows = tuple(
            name for name in screen.boundary_colors if name not in color_index
        )
        unknown_columns = tuple(
            name for name in screen.unit_functions if name not in function_index
        )
        if not screen.boundary_colors or not screen.unit_functions:
            diagnostics.append(
                Diagnostic(
                    "toroidal.affine_screen_shape",
                    "the affine screen needs at least one boundary row and "
                    "one unit-function column",
                    True,
                )
            )
        elif unknown_rows or unknown_columns:
            diagnostics.append(
                Diagnostic(
                    "toroidal.affine_screen_names",
                    f"affine screen has unknown rows {unknown_rows} or "
                    f"columns {unknown_columns}",
                    True,
                )
            )
        elif not screen.normal_core_certificate.strip():
            diagnostics.append(
                Diagnostic(
                    "toroidal.affine_screen_uncertified",
                    "the affine unit/class screen needs normal-core and "
                    "complete-boundary certificates",
                    False,
                )
            )
            uncertified = True
        else:
            rows = tuple(color_index[name] for name in screen.boundary_colors)
            columns = tuple(
                function_index[name] for name in screen.unit_functions
            )
            matrix = _submatrix(datum.valuation_matrix, rows, columns)
            profile = factorial_core_invariants(matrix)
            affine_space_certified = bool(
                screen.affine_toric_sufficiency_certificate.strip()
                and profile["unit_rank"] == 0
                and profile["class_group_free_rank"] == 0
                and not profile["class_group_torsion"]
            )
            invariants["affine_screen"] = {
                **profile,
                "affine_space_certified": affine_space_certified,
            }
            if screen.require_constant_units and profile["unit_rank"] != 0:
                diagnostics.append(
                    Diagnostic(
                        "toroidal.affine_units",
                        "the toroidal affine screen leaves nonconstant units",
                        True,
                    )
                )
            if screen.require_trivial_class_group and (
                profile["class_group_free_rank"] != 0
                or profile["class_group_torsion"]
            ):
                diagnostics.append(
                    Diagnostic(
                        "toroidal.affine_class_group",
                        "the toroidal affine screen has nontrivial class group",
                        True,
                    )
                )

    span_profiles: dict[str, object] = {}
    if not malformed:
        for problem in datum.divisor_span_problems:
            unknown_columns = tuple(
                name
                for name in problem.generator_functions
                if name not in function_index
            )
            bad_shape = (
                not problem.generator_functions
                or len(set(problem.generator_functions))
                != len(problem.generator_functions)
                or bool(unknown_columns)
                or len(problem.target_orders) != len(color_names)
            )
            if bad_shape:
                diagnostics.append(
                    Diagnostic(
                        "toroidal.divisor_span_shape",
                        f"divisor-span problem {problem.name} has unknown "
                        f"generators {unknown_columns} or incompatible widths",
                        True,
                    )
                )
                continue
            columns = tuple(
                function_index[name] for name in problem.generator_functions
            )
            matrix = _submatrix(
                datum.valuation_matrix,
                tuple(range(len(color_names))),
                columns,
            )
            matrix_rank_value = matrix_rank(matrix)
            augmented = tuple(
                row + (problem.target_orders[index],)
                for index, row in enumerate(matrix)
            )
            augmented_rank = matrix_rank(augmented)
            target_class_order = core_class_order(
                matrix, problem.target_orders
            )
            integral_span = target_class_order == 1
            witnesses = colored_proportionality_witnesses(
                matrix,
                problem.target_orders,
                color_names,
            )
            span_profiles[problem.name] = {
                "generator_functions": problem.generator_functions,
                "matrix_rank": matrix_rank_value,
                "augmented_rank": augmented_rank,
                "target_class_order": target_class_order,
                "in_integral_span": integral_span,
                "proportionality_witnesses": witnesses,
            }
            if integral_span:
                continue
            if (
                problem.infeasibility_is_obstruction
                and problem.exhaustive_scope_certificate.strip()
            ):
                diagnostics.append(
                    Diagnostic(
                        "toroidal.divisor_span_obstruction",
                        f"target of exact divisor-span problem {problem.name} "
                        "is outside the generator lattice",
                        True,
                    )
                )
            else:
                diagnostics.append(
                    Diagnostic(
                        "toroidal.divisor_span_uncertified",
                        f"target of divisor-span problem {problem.name} is "
                        "outside the named generator lattice, but the input "
                        "does not certify that the generators are exhaustive",
                        False,
                    )
                )
                uncertified = True
    invariants["divisor_span_problems"] = span_profiles

    problem_profiles: dict[str, object] = {}
    if not malformed:
        for problem in datum.feasibility_problems:
            variable_count = len(problem.variables)
            variable_names = tuple(variable.name for variable in problem.variables)
            bad_shape = (
                len(set(variable_names)) != variable_count
                or len(problem.fixed_function_coefficients)
                != len(function_names)
                or len(problem.variable_function_coefficients)
                != variable_count
                or any(
                    len(row) != len(function_names)
                    for row in problem.variable_function_coefficients
                )
                or len(problem.target_orders) != len(color_names)
                or any(
                    len(constraint.coefficients) != variable_count
                    for constraint in problem.constraints
                )
            )
            if bad_shape:
                diagnostics.append(
                    Diagnostic(
                        "toroidal.feasibility_shape",
                        f"feasibility problem {problem.name} has incompatible "
                        "variable, function, or boundary widths",
                        True,
                    )
                )
                continue
            search_size = 1
            for variable in problem.variables:
                search_size *= variable.upper_bound - variable.lower_bound + 1
            if search_size > 1_000_000:
                diagnostics.append(
                    Diagnostic(
                        "toroidal.feasibility_too_large",
                        f"problem {problem.name} has {search_size} bounded "
                        "assignments; split or tighten it before exact replay",
                        False,
                    )
                )
                uncertified = True
                continue

            models: list[tuple[int, ...]] = []
            domains = tuple(
                range(variable.lower_bound, variable.upper_bound + 1)
                for variable in problem.variables
            )
            for values in product(*domains):
                if not all(
                    _constraint_holds(constraint, values)
                    for constraint in problem.constraints
                ):
                    continue
                coefficients = tuple(
                    fixed
                    + sum(
                        value * variable_row[column]
                        for value, variable_row in zip(
                            values, problem.variable_function_coefficients
                        )
                    )
                    for column, fixed in enumerate(
                        problem.fixed_function_coefficients
                    )
                )
                orders = tuple(
                    sum(order * coefficient for order, coefficient in zip(row, coefficients))
                    for row in datum.valuation_matrix
                )
                if orders == problem.target_orders:
                    models.append(values)

            minimal = _pareto_minimal_models(models)
            problem_profiles[problem.name] = {
                "variables": variable_names,
                "search_size": search_size,
                "model_count": len(models),
                "minimal_models": tuple(
                    dict(zip(variable_names, model)) for model in minimal
                ),
            }
            if not models:
                if (
                    problem.infeasibility_is_obstruction
                    and problem.exhaustive_scope_certificate.strip()
                ):
                    diagnostics.append(
                        Diagnostic(
                            "toroidal.feasibility_infeasible",
                            f"exact exhaustive problem {problem.name} has no "
                            "integral model",
                            True,
                        )
                    )
                else:
                    diagnostics.append(
                        Diagnostic(
                            "toroidal.feasibility_bounded_empty",
                            f"bounded problem {problem.name} has no model, but "
                            "its input does not certify a global obstruction",
                            False,
                        )
                    )
                    uncertified = True
    invariants["feasibility_problems"] = problem_profiles
    invariants["nonlinear_residue"] = datum.nonlinear_residue

    if any(diagnostic.obstruction for diagnostic in diagnostics):
        status = TropicalFeasibilityStatus.OBSTRUCTED
    elif uncertified:
        status = TropicalFeasibilityStatus.UNCERTIFIED
    else:
        status = TropicalFeasibilityStatus.FEASIBLE
        if datum.nonlinear_residue:
            diagnostics.append(
                Diagnostic(
                    "toroidal.nonlinear_residue",
                    "all declared tropical gates pass; the listed nonlinear "
                    "residue remains",
                    False,
                )
            )
    return ToroidalAudit(
        status=status,
        diagnostics=tuple(diagnostics),
        invariants=tuple(sorted(invariants.items())),
    )


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

    toroidal = audit_toroidal_boundary(package.toroidal_boundary)
    invariants["toroidal_boundary"] = toroidal.to_dict()
    diagnostics.extend(toroidal.diagnostics)

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

    retained_root_euler = audit_retained_root_euler(
        package.retained_root_euler
    )
    invariants["retained_root_euler_gate"] = retained_root_euler.to_dict()
    if retained_root_euler.status is RetainedRootEulerStatus.OBSTRUCTED:
        diagnostics.append(
            Diagnostic(
                "affine_source.retained_root_euler",
                retained_root_euler.reason,
                True,
            )
        )
    elif retained_root_euler.status is RetainedRootEulerStatus.UNCERTIFIED:
        diagnostics.append(
            Diagnostic(
                "affine_source.retained_root_euler_uncertified",
                retained_root_euler.reason,
                False,
            )
        )

    conductor_jet = audit_conductor_jet_truncation(
        package.conductor_jet_branches
    )
    conductor_sensitivity = audit_conductor_sensitivity_ledger(
        package.conductor_jet_sensitivity
    )
    if (
        package.conductor_jet_branches is not None
        and package.conductor_jet_sensitivity is not None
    ):
        raise ValueError(
            "declare either the scalar conductor jet ledger or the "
            "dependency-sensitive ledger, not both"
        )
    invariants["conductor_jet_truncation"] = conductor_jet.to_dict()
    invariants["conductor_jet_sensitivity"] = conductor_sensitivity.to_dict()
    if conductor_jet.status in {
        ConductorJetStatus.UNCERTIFIED,
        ConductorJetStatus.INSUFFICIENT,
    }:
        diagnostics.append(
            Diagnostic(
                "boundary_module.conductor_jet_truncation",
                conductor_jet.reason,
                False,
            )
        )
    if conductor_sensitivity.status in {
        ConductorJetStatus.UNCERTIFIED,
        ConductorJetStatus.INSUFFICIENT,
    }:
        diagnostics.append(
            Diagnostic(
                "boundary_module.conductor_jet_sensitivity",
                conductor_sensitivity.reason,
                False,
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

    factorial_core = package.affine_source.factorial_core
    presented_core = package.affine_source.presented_core
    if factorial_core is not None and presented_core is not None:
        diagnostics.append(
            Diagnostic(
                "affine_source.multiple_cores",
                "supply either a factorial core or a presented-class core",
                True,
            )
        )
    if factorial_core is not None:
        row_count = len(factorial_core.valuation_matrix)
        row_widths = {len(row) for row in factorial_core.valuation_matrix}
        malformed = False
        if row_count != len(factorial_core.boundary_primes):
            diagnostics.append(
                Diagnostic(
                    "affine_source.core_boundary_width",
                    "factorial-core rows must index all listed boundary primes",
                    True,
                )
            )
            malformed = True
        if row_widths != {len(factorial_core.unit_generators)}:
            diagnostics.append(
                Diagnostic(
                    "affine_source.core_unit_width",
                    "factorial-core columns must index the unit-lattice basis",
                    True,
                )
            )
            malformed = True
        if len(set(factorial_core.boundary_primes)) != len(
            factorial_core.boundary_primes
        ):
            diagnostics.append(
                Diagnostic(
                    "affine_source.core_boundary_names",
                    "factorial-core boundary-prime names must be distinct",
                    True,
                )
            )
            malformed = True
        if len(set(factorial_core.unit_generators)) != len(
            factorial_core.unit_generators
        ):
            diagnostics.append(
                Diagnostic(
                    "affine_source.core_unit_names",
                    "factorial-core unit-generator names must be distinct",
                    True,
                )
            )
            malformed = True
        for query in factorial_core.class_queries:
            if len(query.vector) != row_count:
                diagnostics.append(
                    Diagnostic(
                        "affine_source.core_class_width",
                        f"class {query.name} needs one entry per boundary prime",
                        True,
                    )
                )
                malformed = True
        if not factorial_core.scope_complete():
            diagnostics.append(
                Diagnostic(
                    "affine_source.core_scope_unproved",
                    "normality, affine UFD structure of the core, and the complete "
                    "codimension-one complement must be certified before the "
                    "valuation matrix is conclusive",
                    False,
                )
            )
        elif not malformed:
            profile = factorial_core_invariants(factorial_core.valuation_matrix)
            invariants["factorial_core"] = profile
            if profile["unit_rank"] != 0:
                diagnostics.append(
                    Diagnostic(
                        "affine_source.core_units",
                        "the factorial-core valuation map leaves nonconstant units",
                        True,
                    )
                )
            if (
                profile["class_group_free_rank"] != 0
                or profile["class_group_torsion"]
            ):
                diagnostics.append(
                    Diagnostic(
                        "affine_source.core_class_group",
                        "the factorial-core valuation map has nontrivial Weil class "
                        "group",
                        True,
                    )
                )
            class_orders = {
                query.name: core_class_order(
                    factorial_core.valuation_matrix, query.vector
                )
                for query in factorial_core.class_queries
            }
            if class_orders:
                invariants["core_class_orders"] = class_orders
            for query in factorial_core.class_queries:
                if query.required_trivial and class_orders[query.name] != 1:
                    order = class_orders[query.name]
                    description = "infinite" if order is None else str(order)
                    diagnostics.append(
                        Diagnostic(
                            "affine_source.core_required_class",
                            f"required-trivial class {query.name} has order "
                            f"{description}",
                            True,
                        )
                    )

    if presented_core is not None:
        boundary_count = len(presented_core.boundary_primes)
        core_generator_count = len(presented_core.core_class_generators)
        malformed = False
        named_lists = (
            ("boundary", presented_core.boundary_primes),
            ("unit", presented_core.unit_generators),
            ("core class", presented_core.core_class_generators),
            ("core relation", presented_core.core_relations),
        )
        for label, names in named_lists:
            if len(set(names)) != len(names):
                diagnostics.append(
                    Diagnostic(
                        "affine_source.presented_core_names",
                        f"presented-core {label} names must be distinct",
                        True,
                    )
                )
                malformed = True
        if len(presented_core.unit_valuation_matrix) != boundary_count:
            malformed = True
        if len(presented_core.core_relation_matrix) != core_generator_count:
            malformed = True
        if len(presented_core.relation_boundary_corrections) != boundary_count:
            malformed = True
        unit_count = len(presented_core.unit_generators)
        relation_count = len(presented_core.core_relations)
        if any(
            len(row) != unit_count
            for row in presented_core.unit_valuation_matrix
        ):
            malformed = True
        if any(
            len(row) != relation_count
            for row in presented_core.core_relation_matrix
        ):
            malformed = True
        if any(
            len(row) != relation_count
            for row in presented_core.relation_boundary_corrections
        ):
            malformed = True
        if malformed:
            diagnostics.append(
                Diagnostic(
                    "affine_source.presented_core_shape",
                    "presented-core matrices do not match their named rows",
                    True,
                )
            )
        expected_query_width = boundary_count + core_generator_count
        for query in presented_core.class_queries:
            if len(query.vector) != expected_query_width:
                diagnostics.append(
                    Diagnostic(
                        "affine_source.presented_core_class_width",
                        f"class {query.name} needs boundary plus core coordinates",
                        True,
                    )
                )
                malformed = True
        if not presented_core.scope_complete():
            diagnostics.append(
                Diagnostic(
                    "affine_source.presented_core_scope_unproved",
                    "normality, the core class presentation, every lifted-relation "
                    "witness, and the complete boundary must be certified",
                    False,
                )
            )
        elif not malformed:
            try:
                profile = presented_core_invariants(
                    presented_core.unit_valuation_matrix,
                    presented_core.core_relation_matrix,
                    presented_core.relation_boundary_corrections,
                )
            except ValueError as error:
                diagnostics.append(
                    Diagnostic(
                        "affine_source.presented_core_shape", str(error), True
                    )
                )
            else:
                invariants["presented_core"] = profile
                if profile["unit_rank"] != 0:
                    diagnostics.append(
                        Diagnostic(
                            "affine_source.presented_core_units",
                            "the presented-core valuation map leaves "
                            "nonconstant units",
                            True,
                        )
                    )
                if (
                    profile["class_group_free_rank"] != 0
                    or profile["class_group_torsion"]
                ):
                    diagnostics.append(
                        Diagnostic(
                            "affine_source.presented_core_class_group",
                            "the lifted block presentation has nontrivial "
                            "Weil class group",
                            True,
                        )
                    )
                class_orders = {
                    query.name: core_class_order(
                        profile["presentation_matrix"], query.vector
                    )
                    for query in presented_core.class_queries
                }
                if class_orders:
                    invariants["presented_core_class_orders"] = class_orders
                for query in presented_core.class_queries:
                    if query.required_trivial and class_orders[query.name] != 1:
                        order = class_orders[query.name]
                        description = "infinite" if order is None else str(order)
                        diagnostics.append(
                            Diagnostic(
                                "affine_source.presented_core_required_class",
                                f"required-trivial class {query.name} has order "
                                f"{description}",
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


def torus_class_group_obstruction_package(
    degree: int = 5,
) -> BoundaryPackage:
    """Attach the balanced wild valuation block to a factorial torus core."""

    if degree <= 2:
        raise ValueError("the obstruction fixture requires N>2")
    package = a4_three_puncture_package()
    return BoundaryPackage(
        name=f"balanced_wild_torus_class_obstruction_N{degree}",
        cover=package.cover,
        selected_curves=package.selected_curves,
        determinant_ledger=package.determinant_ledger,
        polar_ledger=package.polar_ledger,
        affine_source=AffineSourceNecessaryData(
            boundary_class_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
            factorial_core=FactorialCoreDatum(
                boundary_primes=("A0", "L1", "Z0"),
                unit_generators=("x", "u", "z"),
                valuation_matrix=(
                    (1, 0, 0),
                    (degree - 2, degree - 1, 0),
                    (0, 0, 1),
                ),
                normality_certificate="balanced normalization theorem",
                factoriality_certificate="Laurent polynomial core",
                complete_boundary_certificate=(
                    "complete codimension-one complement of the Laurent core"
                ),
                class_queries=(
                    CoreClassQuery(
                        "L1", (0, 1, 0), required_trivial=True
                    ),
                    CoreClassQuery("div(x)", (1, degree - 2, 0)),
                ),
            ),
        ),
    )


def presented_core_extension_obstruction_package() -> BoundaryPackage:
    """A lifted core relation producing a nonsplit ``Z/4`` class group."""

    package = a4_three_puncture_package()
    return BoundaryPackage(
        name="presented_core_nonsplit_z4_obstruction",
        cover=package.cover,
        selected_curves=package.selected_curves,
        determinant_ledger=package.determinant_ledger,
        polar_ledger=package.polar_ledger,
        affine_source=AffineSourceNecessaryData(
            boundary_class_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
            presented_core=PresentedCoreDatum(
                boundary_primes=("D",),
                unit_generators=("t",),
                unit_valuation_matrix=((2,),),
                core_class_generators=("G",),
                core_relations=("2G",),
                core_relation_matrix=((2,),),
                relation_boundary_corrections=((1,),),
                normality_certificate="normal reconstruction theorem",
                core_class_presentation_certificate="Cl(W)=<G|2G=0>",
                relation_lift_certificate="div(f)=D+2G",
                complete_boundary_certificate="D is the full complement",
                class_queries=(
                    CoreClassQuery("G", (0, 1), required_trivial=True),
                ),
            ),
        ),
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


def retained_root_euler_obstruction_package(
    retained_degree: int = 4,
) -> BoundaryPackage:
    """Attach the certified balanced retained-root gate to a passing package."""

    package = a4_three_puncture_package()
    return BoundaryPackage(
        name=f"retained_root_euler_obstruction_degree_{retained_degree}",
        cover=package.cover,
        selected_curves=package.selected_curves,
        determinant_ledger=package.determinant_ledger,
        polar_ledger=package.polar_ledger,
        affine_source=package.affine_source,
        retained_root_euler=RetainedRootEulerDatum(
            retained_degree=retained_degree,
            squarefree=True,
            nonzero_constant_term=True,
            omitted_fierce_boundary_count=1,
            balanced_chart_certificate=(
                "P=x*u,T=x^2*u,Q=A(x^2*u)*(u-x^(N-1)) over D(A)"
            ),
            different_support_certificate=(
                "relative different supported on the omitted fierce boundary"
            ),
            root_fibre_certificate=(
                "each simple nonzero retained root replaces Gm by two A1s"
            ),
            omitted_boundary_certificate=(
                "E is the single A1 parametrized by "
                "P=s^N,T=s^(N+1),Q=0"
            ),
            theorem_source=(
                "extended-geometry/PLANE_WILD_BOUNDARY_ATLAS.md, "
                "Theorem 7.12"
            ),
        ),
    )


def a4_toroidal_ledger_datum() -> ToroidalBoundaryDatum:
    """Colored-orthant encoding of the exact pure-target A4 ledger."""

    return ToroidalBoundaryDatum(
        fan=ToroidalFanDatum(
            lattice_basis=("ord_W", "ord_K", "ord_L"),
            rays=(
                TropicalRay("W_ray", (1, 0, 0)),
                TropicalRay("K_ray", (0, 1, 0)),
                TropicalRay("L_ray", (0, 0, 1)),
            ),
            maximal_cones=(("W_ray", "K_ray", "L_ray"),),
            incidence_certificate=(
                "A4_PURE_TARGET_LEDGER_LIFT: divisorial valuation orthant"
            ),
        ),
        boundary_colors=(
            BoundaryColor("W", SheetColor.BOUNDARY, "W_ray"),
            BoundaryColor("K", SheetColor.BOUNDARY, "K_ray"),
            BoundaryColor("L", SheetColor.BOUNDARY, "L_ray"),
        ),
        valuation_functions=(
            "det_DPhi",
            "auxiliary_mask",
            "pullback_target_boundary",
        ),
        valuation_matrix=(
            (2, 1, 3),
            (3, 0, 3),
            (1, 1, 2),
        ),
        valuation_certificate=(
            "A4_PURE_TARGET_LEDGER_LIFT, valuation table (W,K,L)"
        ),
        identities=(
            ValuationIdentity(
                "pure_target_log_keller_balance",
                (1, 1, -1),
                (0, 0, 0),
            ),
        ),
        nonlinear_residue=(
            "realize two source-dependent masks with entrywise polynomial inverse",
            "certify smooth affine-space source and target rather than chart opens",
            "preserve the oriented quartic field and complete regular fibres",
        ),
    )


def d5_toroidal_ledger_datum() -> ToroidalBoundaryDatum:
    """Newton-fan and colored-branch encoding of the D5 valuation gate."""

    return ToroidalBoundaryDatum(
        fan=ToroidalFanDatum(
            lattice_basis=("ord_a", "ord_u"),
            rays=(
                TropicalRay("a_axis", (1, 0)),
                TropicalRay("first_blowup", (1, 1)),
                TropicalRay("parabolic_scale", (1, 2)),
                TropicalRay("u_axis", (0, 1)),
            ),
            maximal_cones=(
                ("a_axis", "first_blowup"),
                ("first_blowup", "parabolic_scale"),
                ("parabolic_scale", "u_axis"),
            ),
            incidence_certificate=(
                "D5-N3: two toric blowups separate C,R_plus,R_minus at "
                "their distinct parabolic residues"
            ),
        ),
        boundary_colors=(
            BoundaryColor(
                "C", SheetColor.BOUNDARY, "parabolic_scale", "u/a^2=1/4"
            ),
            BoundaryColor(
                "R_plus",
                SheetColor.BOUNDARY,
                "parabolic_scale",
                "u/a^2=2/(3+sqrt(5))",
            ),
            BoundaryColor(
                "R_minus",
                SheetColor.BOUNDARY,
                "parabolic_scale",
                "u/a^2=2/(3-sqrt(5))",
            ),
            BoundaryColor(
                "E_11", SheetColor.BOUNDARY, "first_blowup"
            ),
            BoundaryColor(
                "E_12", SheetColor.BOUNDARY, "parabolic_scale"
            ),
        ),
        valuation_functions=(
            "derivative_Q",
            "pullback_Delta",
            "source_C",
            "source_R_plus",
            "source_R_minus",
        ),
        valuation_matrix=(
            (0, 1, 1, 0, 0),
            (1, 2, 0, 1, 0),
            (1, 2, 0, 0, 1),
            (2, 5, 1, 1, 1),
            (4, 10, 2, 2, 2),
        ),
        valuation_certificate=(
            "D5-N3 and Delta=C*R_plus^2*R_minus^2; exceptional rows use "
            "weights (1,1) and (1,2)"
        ),
        feasibility_problems=(
            ValuationFeasibilityProblem(
                name="branch_supported_log_balance_m_le_4",
                variables=(
                    IntegralVariable("m", 0, 4),
                    IntegralVariable("s_C", 0, 8),
                    IntegralVariable("s_plus", 0, 8),
                    IntegralVariable("s_minus", 0, 8),
                ),
                fixed_function_coefficients=(1, 0, 0, 0, 0),
                variable_function_coefficients=(
                    (0, -1, 0, 0, 0),
                    (0, 0, 1, 0, 0),
                    (0, 0, 0, 1, 0),
                    (0, 0, 0, 0, 1),
                ),
                target_orders=(0, 0, 0, 0, 0),
                exhaustive_scope_certificate=(
                    "D5-N3 proves s=(m,2m-1,2m-1); this finite replay "
                    "contains its primitive model and first three translates"
                ),
            ),
        ),
        nonlinear_residue=(
            "move both ramification colors and the old incidence graph",
            "clear every inverse-adjugate entry, not only its determinant",
            "recognize smooth factorial affine spaces with constant units",
            "prove finite flat rank five and the natural D5 normal closure",
        ),
    )


def davenport_toroidal_ledger_datum() -> ToroidalBoundaryDatum:
    """One-ray colored fan for the Davenport three-column Cox ledger."""

    return ToroidalBoundaryDatum(
        fan=ToroidalFanDatum(
            lattice_basis=("ord_Delta",),
            rays=(TropicalRay("Delta_ray", (1,)),),
            maximal_cones=(("Delta_ray",),),
            incidence_certificate=(
                "DAVENPORT_COX_BOUNDARY_OBSTRUCTION: generic Delta-normal fan"
            ),
        ),
        boundary_colors=(
            BoundaryColor(
                "E3", SheetColor.BOUNDARY, "Delta_ray", "unramified degree 3"
            ),
            BoundaryColor(
                "E6", SheetColor.BOUNDARY, "Delta_ray", "unramified degree 6"
            ),
            BoundaryColor(
                "J", SheetColor.BOUNDARY, "Delta_ray", "ramification color"
            ),
        ),
        valuation_functions=(
            "pullback_Delta",
            "jacobian_J",
            "primitive_E3_character",
        ),
        valuation_matrix=(
            (1, 0, 1),
            (1, 0, 0),
            (2, 1, 0),
        ),
        valuation_certificate=(
            "DAVENPORT_COX_BOUNDARY_OBSTRUCTION, equations (2.5)-(2.8)"
        ),
        unimodular_blocks=(
            UnimodularBlockRequirement(
                "cox_lattice_completion",
                ("E3", "E6", "J"),
                (
                    "pullback_Delta",
                    "jacobian_J",
                    "primitive_E3_character",
                ),
            ),
        ),
        nonlinear_residue=(
            "realize the primitive E3 character without the new divisor L(Y)",
            "fill Delta while controlling all three source colors",
            "preserve the GL3(F2) Gassmann closure on the common open",
            "prove that the modified source and target are affine spaces",
        ),
    )


def f20_toroidal_ledger_datum() -> ToroidalBoundaryDatum:
    """Finite toroidal colors of the corrected Lecacheux root cover."""

    # The three q-r tangencies form one cubic Galois orbit.  Their root-cover
    # calculation is identical over the residue field, so retain one template
    # instead of three independent family ledgers.
    qr_tangent_templates = (
        (
            "E1_A_ramified_2",
            "first exceptional, index-two part of the triple cluster",
            (0, 2, 2, 2),
        ),
        (
            "E1_A_unramified",
            "first exceptional, unramified part of the triple cluster",
            (0, 1, 1, 1),
        ),
        (
            "E1_B_ramified_2",
            "first exceptional, index-two double cluster",
            (0, 2, 2, 1),
        ),
        (
            "E2_A_sheet_1",
            "second exceptional, first triple-cluster sheet",
            (0, 2, 2, 2),
        ),
        (
            "E2_A_sheet_2",
            "second exceptional, second triple-cluster sheet",
            (0, 2, 2, 2),
        ),
        (
            "E2_A_sheet_3",
            "second exceptional, third triple-cluster sheet",
            (0, 2, 2, 2),
        ),
        (
            "E2_B_sheet_1",
            "second exceptional, first double-cluster sheet",
            (0, 2, 2, 1),
        ),
        (
            "E2_B_sheet_2",
            "second exceptional, second double-cluster sheet",
            (0, 2, 2, 1),
        ),
    )
    qr_tangent_colors = tuple(
        BoundaryColor(
            f"qr_tangent_{center}_{suffix}",
            SheetColor.BOUNDARY,
            "qr_tangent_ray",
            f"cubic-orbit center {center}: {label}",
        )
        for center in range(1, 4)
        for suffix, label, _row in qr_tangent_templates
    )
    qr_tangent_rows = tuple(
        row
        for _center in range(1, 4)
        for _suffix, _label, row in qr_tangent_templates
    )

    base = ToroidalBoundaryDatum(
        fan=ToroidalFanDatum(
            lattice_basis=("ord_d", "ord_q", "ord_r"),
            rays=(
                TropicalRay("d_ray", (1, 0, 0)),
                TropicalRay("q_ray", (0, 1, 0)),
                TropicalRay("r_ray", (0, 0, 1)),
                TropicalRay("qr_tangent_ray", (0, 1, 1)),
                TropicalRay("triple_E1_ray", (1, 1, 1)),
                TropicalRay("triple_E2_ray", (2, 2, 1)),
            ),
            maximal_cones=(
                ("triple_E1_ray", "q_ray", "qr_tangent_ray"),
                ("triple_E1_ray", "qr_tangent_ray", "r_ray"),
                ("d_ray", "triple_E1_ray", "r_ray"),
                ("triple_E2_ray", "q_ray", "triple_E1_ray"),
                ("d_ray", "triple_E2_ray", "triple_E1_ray"),
                ("d_ray", "q_ray", "triple_E2_ray"),
            ),
            incidence_certificate=(
                "corrected Lecacheux F20 discriminant: generic "
                "divisorial valuation orthant with the two smooth star "
                "subdivisions forced by the conjugate triple tangencies and "
                "the cubic orbit of q-r tangencies"
            ),
        ),
        boundary_colors=(
            BoundaryColor(
                "d_unramified",
                SheetColor.BOUNDARY,
                "d_ray",
                "X=1/4, e=1",
            ),
            BoundaryColor(
                "d_ramified_4",
                SheetColor.BOUNDARY,
                "d_ray",
                "X=1+s/2, e=4",
            ),
            BoundaryColor(
                "q_collision_plus",
                SheetColor.BOUNDARY,
                "q_ray",
                "first normalization branch through X=a_q, e=1",
            ),
            BoundaryColor(
                "q_collision_minus",
                SheetColor.BOUNDARY,
                "q_ray",
                "second normalization branch through X=a_q, e=1",
            ),
            BoundaryColor(
                "q_residual_1",
                SheetColor.BOUNDARY,
                "q_ray",
                "first geometric residual-cubic sheet, e=1",
            ),
            BoundaryColor(
                "q_residual_2",
                SheetColor.BOUNDARY,
                "q_ray",
                "second geometric residual-cubic sheet, e=1",
            ),
            BoundaryColor(
                "q_residual_3",
                SheetColor.BOUNDARY,
                "q_ray",
                "third geometric residual-cubic sheet, e=1",
            ),
            BoundaryColor(
                "q_node_slope_1",
                SheetColor.BOUNDARY,
                "q_ray",
                "first geometric slope over the q-node blowup, e=1",
            ),
            BoundaryColor(
                "q_node_slope_2",
                SheetColor.BOUNDARY,
                "q_ray",
                "second geometric slope over the q-node blowup, e=1",
            ),
            BoundaryColor(
                "q_node_slope_3",
                SheetColor.BOUNDARY,
                "q_ray",
                "third geometric slope over the q-node blowup, e=1",
            ),
            BoundaryColor(
                "q_node_slope_4",
                SheetColor.BOUNDARY,
                "q_ray",
                "fourth geometric slope over the q-node blowup, e=1",
            ),
            BoundaryColor(
                "q_node_simple",
                SheetColor.BOUNDARY,
                "q_ray",
                "simple X=-1 sheet over the q-node blowup, e=1",
            ),
            BoundaryColor(
                "r_ramified_2_plus",
                SheetColor.BOUNDARY,
                "r_ray",
                "first geometric root of the repeated quadratic, e=2",
            ),
            BoundaryColor(
                "r_ramified_2_minus",
                SheetColor.BOUNDARY,
                "r_ray",
                "second geometric root of the repeated quadratic, e=2",
            ),
            BoundaryColor(
                "r_unramified",
                SheetColor.BOUNDARY,
                "r_ray",
                "simple residual root, e=1",
            ),
            BoundaryColor(
                "r_cusp_E1_total_5",
                SheetColor.BOUNDARY,
                "r_ray",
                "first cusp exceptional, total index five",
            ),
            BoundaryColor(
                "r_cusp_E2_total_5",
                SheetColor.BOUNDARY,
                "r_ray",
                "second cusp exceptional, total index five",
            ),
            BoundaryColor(
                "r_cusp_E3_unramified",
                SheetColor.BOUNDARY,
                "r_ray",
                "third cusp exceptional, one unramified sheet",
            ),
            BoundaryColor(
                "r_cusp_E3_ramified_plus",
                SheetColor.BOUNDARY,
                "r_ray",
                "third cusp exceptional, first index-two color",
            ),
            BoundaryColor(
                "r_cusp_E3_ramified_minus",
                SheetColor.BOUNDARY,
                "r_ray",
                "third cusp exceptional, second index-two color",
            ),
            BoundaryColor(
                "r_cusp_E4_sheet_1",
                SheetColor.BOUNDARY,
                "r_ray",
                "fourth cusp exceptional, first unramified sheet",
            ),
            BoundaryColor(
                "r_cusp_E4_sheet_2",
                SheetColor.BOUNDARY,
                "r_ray",
                "fourth cusp exceptional, second unramified sheet",
            ),
            BoundaryColor(
                "r_cusp_E4_sheet_3",
                SheetColor.BOUNDARY,
                "r_ray",
                "fourth cusp exceptional, third unramified sheet",
            ),
            BoundaryColor(
                "r_cusp_E4_sheet_4",
                SheetColor.BOUNDARY,
                "r_ray",
                "fourth cusp exceptional, fourth unramified sheet",
            ),
            BoundaryColor(
                "r_cusp_E4_sheet_5",
                SheetColor.BOUNDARY,
                "r_ray",
                "fourth cusp exceptional, fifth unramified sheet",
            ),
            BoundaryColor(
                "triple_plus_E1_ramified_4",
                SheetColor.BOUNDARY,
                "triple_E1_ray",
                "positive conjugate center, index-four cluster",
            ),
            BoundaryColor(
                "triple_plus_E1_simple",
                SheetColor.BOUNDARY,
                "triple_E1_ray",
                "positive conjugate center, simple sheet",
            ),
            BoundaryColor(
                "triple_minus_E1_ramified_4",
                SheetColor.BOUNDARY,
                "triple_E1_ray",
                "negative conjugate center, index-four cluster",
            ),
            BoundaryColor(
                "triple_minus_E1_simple",
                SheetColor.BOUNDARY,
                "triple_E1_ray",
                "negative conjugate center, simple sheet",
            ),
            BoundaryColor(
                "triple_plus_E2_cluster_1",
                SheetColor.BOUNDARY,
                "triple_E2_ray",
                "positive conjugate center, first cluster sheet",
            ),
            BoundaryColor(
                "triple_plus_E2_cluster_2",
                SheetColor.BOUNDARY,
                "triple_E2_ray",
                "positive conjugate center, second cluster sheet",
            ),
            BoundaryColor(
                "triple_plus_E2_cluster_3",
                SheetColor.BOUNDARY,
                "triple_E2_ray",
                "positive conjugate center, third cluster sheet",
            ),
            BoundaryColor(
                "triple_plus_E2_cluster_4",
                SheetColor.BOUNDARY,
                "triple_E2_ray",
                "positive conjugate center, fourth cluster sheet",
            ),
            BoundaryColor(
                "triple_plus_E2_simple",
                SheetColor.BOUNDARY,
                "triple_E2_ray",
                "positive conjugate center, simple sheet",
            ),
            BoundaryColor(
                "triple_minus_E2_cluster_1",
                SheetColor.BOUNDARY,
                "triple_E2_ray",
                "negative conjugate center, first cluster sheet",
            ),
            BoundaryColor(
                "triple_minus_E2_cluster_2",
                SheetColor.BOUNDARY,
                "triple_E2_ray",
                "negative conjugate center, second cluster sheet",
            ),
            BoundaryColor(
                "triple_minus_E2_cluster_3",
                SheetColor.BOUNDARY,
                "triple_E2_ray",
                "negative conjugate center, third cluster sheet",
            ),
            BoundaryColor(
                "triple_minus_E2_cluster_4",
                SheetColor.BOUNDARY,
                "triple_E2_ray",
                "negative conjugate center, fourth cluster sheet",
            ),
            BoundaryColor(
                "triple_minus_E2_simple",
                SheetColor.BOUNDARY,
                "triple_E2_ray",
                "negative conjugate center, simple sheet",
            ),
        )
        + qr_tangent_colors,
        valuation_functions=(
            "pullback_d",
            "pullback_q",
            "pullback_r",
            "derivative_P_X",
        ),
        valuation_matrix=(
            (1, 0, 0, 0),
            (4, 0, 0, 3),
            (0, 1, 0, 1),
            (0, 1, 0, 1),
            (0, 1, 0, 0),
            (0, 1, 0, 0),
            (0, 1, 0, 0),
            (0, 2, 0, 1),
            (0, 2, 0, 1),
            (0, 2, 0, 1),
            (0, 2, 0, 1),
            (0, 2, 0, 0),
            (0, 0, 2, 1),
            (0, 0, 2, 1),
            (0, 0, 1, 0),
            (0, 0, 10, 4),
            (0, 0, 20, 8),
            (0, 0, 5, 2),
            (0, 0, 10, 4),
            (0, 0, 10, 4),
            (0, 0, 10, 4),
            (0, 0, 10, 4),
            (0, 0, 10, 4),
            (0, 0, 10, 4),
            (0, 0, 10, 4),
            (4, 4, 4, 7),
            (1, 1, 1, 0),
            (4, 4, 4, 7),
            (1, 1, 1, 0),
            (2, 2, 1, 3),
            (2, 2, 1, 3),
            (2, 2, 1, 3),
            (2, 2, 1, 3),
            (2, 2, 1, 0),
            (2, 2, 1, 3),
            (2, 2, 1, 3),
            (2, 2, 1, 3),
            (2, 2, 1, 3),
            (2, 2, 1, 0),
        )
        + qr_tangent_rows,
        valuation_certificate=(
            "Disc_X(P)=d^3*q^2*r^2/256; modulo d the root profile is "
            "(4,1), q is a transverse unramified double-root collision, "
            "the q-node blowup has four simple slope colors plus X=-1, "
            "modulo r the generic profile is (2,2,1), and the four cusp "
            "exceptionals have profiles 5, 5, (1,2,2), and (1,1,1,1,1); "
            "the conjugate triple centers have E1 profile (4,1) and E2 "
            "profile (1,1,1,1,1); at each cubic-orbit q-r tangency the "
            "center profile is (3,2), E1 has colors (2,1,2), and E2 has "
            "five unramified colors"
        ),
        nonlinear_residue=(),
    )

    # These four proposed mask columns are deliberately separated from the
    # proved divisor columns.  The first three are the base-factor monomials;
    # the last is the conductor selector w-1.  Its divisor is supported over
    # the boundary at infinity, so its orders on the complete finite atlas
    # are zero.  The q-normalization calculation in the exact verifier proves
    # that this selector completes the conductor unit lattice unimodularly.
    return replace(
        base,
        valuation_functions=base.valuation_functions
        + (
            "mask_d",
            "mask_q",
            "mask_r",
            "q_selector_w_minus_1",
        ),
        valuation_matrix=tuple(
            row + row[:3] + (0,) for row in base.valuation_matrix
        ),
        valuation_certificate=(
            base.valuation_certificate
            + "; the candidate mask columns are d, q, r, and the q-conductor "
            "selector w-1, whose finite orders vanish"
        ),
        nonlinear_residue=(
            "realize a non-base-factor Cox or mask divisor after the minimal "
            "base-factor-plus-selector ansatz fails",
            "clear every inverse-adjugate entry for each surviving mask model",
            "recognize affine-space source and target while preserving F20",
        ),
    )


def f20_base_factor_mask_datum() -> ToroidalBoundaryDatum:
    """Minimal F20 mask ansatz using d, q, r and one q selector.

    The q-conductor unit computation reduces selector powers modulo pullback
    units to exponent zero or one.  Nonnegative base-factor exponents above
    eight cannot meet a derivative target whose largest order is eight, so
    the finite box is exhaustive for this Laurent-monomial architecture.
    """

    base = f20_toroidal_ledger_datum()
    function_count = len(base.valuation_functions)

    def coefficient(column: int) -> tuple[int, ...]:
        return tuple(
            1 if index == column else 0
            for index in range(function_count)
        )

    scope_certificate = (
        "all Laurent-monomial divisors generated by d, q, r, and w-1; "
        "Norm(w-1)=8/(y+3) reduces selector powers modulo pullback units"
    )
    span_problem = ColoredDivisorSpanProblem(
        name="uncolored_base_factor_plus_q_selector_span",
        generator_functions=(
            "mask_d",
            "mask_q",
            "mask_r",
            "q_selector_w_minus_1",
        ),
        target_orders=tuple(row[3] for row in base.valuation_matrix),
        infeasibility_is_obstruction=True,
        exhaustive_scope_certificate=scope_certificate,
    )
    problem = ValuationFeasibilityProblem(
        name="base_factor_plus_q_selector_derivative_principalization",
        variables=(
            IntegralVariable("exponent_d", 0, 8),
            IntegralVariable("exponent_q", 0, 8),
            IntegralVariable("exponent_r", 0, 8),
            IntegralVariable("selector_parity", 0, 1),
        ),
        fixed_function_coefficients=(0,) * function_count,
        variable_function_coefficients=(
            coefficient(4),
            coefficient(5),
            coefficient(6),
            coefficient(7),
        ),
        target_orders=tuple(row[3] for row in base.valuation_matrix),
        infeasibility_is_obstruction=True,
        exhaustive_scope_certificate=(
            "nonnegative part of the Laurent-monomial architecture: "
            "d^a*q^b*r^c*(w-1)^e "
            "with a,b,c nonnegative; Norm(w-1)=8/(y+3) reduces e modulo "
            "pullback units, and any exponent above eight overshoots a "
            "positive factor row of the derivative target"
        ),
    )
    return replace(
        base,
        divisor_span_problems=base.divisor_span_problems + (span_problem,),
        feasibility_problems=base.feasibility_problems + (problem,),
        nonlinear_residue=(
            "no model survives to inverse-adjugate or affine-space tests "
            "inside the certified base-factor-plus-selector architecture",
            "a broader construction must add a genuinely colored Cox divisor",
        ),
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
        torus_class_group_obstruction_package(),
        presented_core_extension_obstruction_package(),
        semigroup_hole_obstruction_package(),
        retained_root_euler_obstruction_package(),
    )


__all__ = [
    "BoundaryColor",
    "ColoredDivisorSpanProblem",
    "BoundaryOutputExpressionDatum",
    "BoundaryPackage",
    "Compilation",
    "ConductorBranchJetDatum",
    "ConductorBranchSensitivityDatum",
    "ContactExpression",
    "NormalJetInputDatum",
    "PackageStatus",
    "IntegralVariable",
    "LinearConstraint",
    "StageTwoRealizationCertificate",
    "ToroidalAffineScreen",
    "ToroidalBoundaryDatum",
    "ToroidalFanDatum",
    "TropicalFeasibilityStatus",
    "TropicalRay",
    "TropicalRegularityRequirement",
    "UnimodularBlockRequirement",
    "ValuationFeasibilityProblem",
    "ValuationIdentity",
    "CoreClassQuery",
    "FactorialCoreDatum",
    "PresentedCoreDatum",
    "RetainedRootEulerDatum",
    "a4_cone_branch_obstruction_package",
    "a4_toroidal_ledger_datum",
    "a4_three_puncture_package",
    "affine_ramification_obstruction_package",
    "all_benchmark_packages",
    "compile_boundary_package",
    "colored_proportionality_witnesses",
    "audit_toroidal_boundary",
    "d5_toroidal_ledger_datum",
    "davenport_toroidal_ledger_datum",
    "determinant_ledger_obstruction_package",
    "elliptic_selected_boundary_package",
    "f20_base_factor_mask_datum",
    "f20_toroidal_ledger_datum",
    "gl3f2_triangle_package",
    "nodal_rational_package",
    "presented_core_extension_obstruction_package",
    "presented_core_invariants",
    "retained_root_euler_obstruction_package",
    "smith_invariant_factors",
    "semigroup_hole_obstruction_package",
    "torus_class_group_obstruction_package",
    "core_class_order",
    "factorial_core_invariants",
    "torus_core_invariants",
]
