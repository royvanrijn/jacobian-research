#!/usr/bin/env python3
"""Exact groundwork for structural elliptic-curve rank searches.

The repository already has deep parameter and point searches.  This module
supplies small, auditable primitives for the complementary structural lanes:

* residual Selmer bookkeeping and the 2-division cubic;
* quotient-aware point promotion (implemented separately);
* bounded Neron--Severi lattice enumeration for K3 multisection/fibration
  candidates;
* V4-cover quotient and genus bookkeeping;
* quadratic-twist/base-change planning; and
* projective p-adic balls, including neighborhoods of infinity.

These routines deliberately stop before the hard external mathematics.
Enumerating a lattice vector does not prove effectiveness, a Selmer class may
belong to Sha, and a cover decomposition does not determine rational points.
The corresponding task manifests must retain those boundaries.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import gcd, isqrt, lcm
from typing import Iterable, Mapping, Sequence


Q = Fraction


def _is_prime(value: int) -> bool:
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, isqrt(value) + 1))


def _primitive_integer_coefficients(
    coefficients_ascending: Sequence[Fraction],
) -> tuple[int, ...]:
    if not coefficients_ascending:
        raise ValueError("a polynomial needs at least one coefficient")
    denominator = 1
    for coefficient in coefficients_ascending:
        denominator = lcm(denominator, Q(coefficient).denominator)
    integers = [int(Q(coefficient) * denominator) for coefficient in coefficients_ascending]
    content = 0
    for value in integers:
        content = gcd(content, abs(value))
    if content:
        integers = [value // content for value in integers]
    while len(integers) > 1 and integers[-1] == 0:
        integers.pop()
    if integers[-1] < 0:
        integers = [-value for value in integers]
    return tuple(integers)


def two_division_cubic(ainvs: Sequence[Fraction]) -> tuple[int, ...]:
    """Return the primitive integral 2-division cubic, ascending.

    For ``[a1,a2,a3,a4,a6]`` the x-coordinates of nonzero 2-torsion satisfy

    ``4*x^3 + b2*x^2 + 2*b4*x + b6 = 0``.

    The returned polynomial is scaled to primitive integer coefficients.  It
    defines the same cubic algebra but is not itself a class-group or Selmer
    computation.
    """

    if len(ainvs) != 5:
        raise ValueError("five generalized Weierstrass coefficients are required")
    a1, a2, a3, a4, a6 = (Q(value) for value in ainvs)
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    return _primitive_integer_coefficients((b6, 2 * b4, b2, Q(4)))


@dataclass(frozen=True)
class ResidualSelmerBudget:
    """Dimension bookkeeping after a certified subgroup is known.

    ``selmer_dimension`` means the dimension of ``Sel_ell(E)``.
    ``known_free_rank`` counts certified independent non-torsion points, and
    ``rational_ell_torsion_dimension`` is ``dim E(Q)[ell]``.  The unexplained
    quotient may contain additional Mordell--Weil directions, Sha classes, or
    both.
    """

    ell: int
    selmer_dimension: int
    known_free_rank: int
    rational_ell_torsion_dimension: int = 0

    def __post_init__(self) -> None:
        if not _is_prime(self.ell):
            raise ValueError("ell must be prime")
        if min(
            self.selmer_dimension,
            self.known_free_rank,
            self.rational_ell_torsion_dimension,
        ) < 0:
            raise ValueError("Selmer dimensions and ranks must be nonnegative")
        known_kummer_dimension = (
            self.known_free_rank + self.rational_ell_torsion_dimension
        )
        if known_kummer_dimension > self.selmer_dimension:
            raise ValueError("the known Kummer image cannot exceed the Selmer group")

    @property
    def rank_upper_bound(self) -> int:
        return self.selmer_dimension - self.rational_ell_torsion_dimension

    @property
    def unexplained_selmer_dimension(self) -> int:
        return (
            self.selmer_dimension
            - self.known_free_rank
            - self.rational_ell_torsion_dimension
        )

    def target_not_excluded(self, target_rank: int) -> bool:
        if target_rank < 0:
            raise ValueError("target rank must be nonnegative")
        return target_rank <= self.rank_upper_bound

    def to_record(self) -> dict[str, object]:
        return {
            "ell": self.ell,
            "selmer_dimension": self.selmer_dimension,
            "known_free_rank": self.known_free_rank,
            "rational_ell_torsion_dimension": self.rational_ell_torsion_dimension,
            "rank_upper_bound": self.rank_upper_bound,
            "unexplained_selmer_dimension": self.unexplained_selmer_dimension,
            "claim_boundary": (
                "unexplained Selmer classes can be Mordell--Weil or Sha; "
                "a positive residual dimension is not a rank lower bound"
            ),
        }


@dataclass(frozen=True)
class IntegralLattice:
    gram: tuple[tuple[int, ...], ...]
    labels: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        dimension = len(self.gram)
        if dimension == 0:
            raise ValueError("the lattice must have positive dimension")
        if any(len(row) != dimension for row in self.gram):
            raise ValueError("the Gram matrix must be square")
        if any(
            self.gram[left][right] != self.gram[right][left]
            for left in range(dimension)
            for right in range(dimension)
        ):
            raise ValueError("the Gram matrix must be symmetric")
        if self.labels and len(self.labels) != dimension:
            raise ValueError("lattice labels must match the Gram dimension")

    @property
    def dimension(self) -> int:
        return len(self.gram)

    def pair(self, left: Sequence[int], right: Sequence[int]) -> int:
        if len(left) != self.dimension or len(right) != self.dimension:
            raise ValueError("lattice vectors have the wrong dimension")
        return sum(
            int(left[row]) * self.gram[row][column] * int(right[column])
            for row in range(self.dimension)
            for column in range(self.dimension)
        )

    def norm(self, vector: Sequence[int]) -> int:
        return self.pair(vector, vector)


@dataclass(frozen=True)
class K3DivisorCandidate:
    vector: tuple[int, ...]
    fiber_degree: int
    self_intersection: int
    arithmetic_genus: int

    def to_record(self) -> dict[str, object]:
        return {
            "vector": list(self.vector),
            "fiber_degree": self.fiber_degree,
            "self_intersection": self.self_intersection,
            "arithmetic_genus": self.arithmetic_genus,
            "claim_boundary": (
                "a lattice class is only a candidate until effectiveness, "
                "nefness and geometric realization are checked"
            ),
        }


def _primitive_vector(vector: Sequence[int]) -> bool:
    content = 0
    for value in vector:
        content = gcd(content, abs(int(value)))
    return content == 1


def _canonical_sign(vector: tuple[int, ...]) -> tuple[int, ...]:
    first = next((value for value in vector if value), 0)
    return tuple(-value for value in vector) if first < 0 else vector


def enumerate_k3_divisor_candidates(
    lattice: IntegralLattice,
    *,
    fiber_vector: Sequence[int],
    coefficient_bound: int,
    fiber_degrees: Iterable[int],
    self_intersections: Iterable[int],
    quotient_by_sign: bool = False,
) -> tuple[K3DivisorCandidate, ...]:
    """Enumerate primitive lattice classes matching K3 adjunction data.

    On a K3 surface, ``p_a(D)=1+D^2/2``.  This finite coefficient box is a
    search manifest, not a complete effective-cone computation.
    """

    if coefficient_bound < 1:
        raise ValueError("coefficient_bound must be positive")
    fiber = tuple(int(value) for value in fiber_vector)
    if len(fiber) != lattice.dimension:
        raise ValueError("fiber_vector has the wrong dimension")
    degrees = {int(value) for value in fiber_degrees}
    intersections = {int(value) for value in self_intersections}
    if any(value <= 0 for value in degrees):
        raise ValueError("multisection degrees must be positive")
    if any(value % 2 for value in intersections):
        raise ValueError("K3 self-intersections are even")

    seen: set[tuple[int, ...]] = set()
    candidates: list[K3DivisorCandidate] = []
    values = range(-coefficient_bound, coefficient_bound + 1)
    for raw in product(values, repeat=lattice.dimension):
        vector = tuple(raw)
        if not any(vector) or not _primitive_vector(vector):
            continue
        key = _canonical_sign(vector) if quotient_by_sign else vector
        if key in seen:
            continue
        degree = lattice.pair(vector, fiber)
        self_intersection = lattice.norm(vector)
        if degree not in degrees or self_intersection not in intersections:
            continue
        genus_numerator = 2 + self_intersection
        if genus_numerator < 0 or genus_numerator % 2:
            continue
        seen.add(key)
        candidates.append(
            K3DivisorCandidate(
                vector=vector,
                fiber_degree=degree,
                self_intersection=self_intersection,
                arithmetic_genus=genus_numerator // 2,
            )
        )
    candidates.sort(
        key=lambda item: (
            item.fiber_degree,
            item.arithmetic_genus,
            item.self_intersection,
            item.vector,
        )
    )
    return tuple(candidates)


def enumerate_isotropic_fibration_candidates(
    lattice: IntegralLattice,
    *,
    ample_vector: Sequence[int],
    coefficient_bound: int,
) -> tuple[tuple[int, ...], ...]:
    """Enumerate primitive positive isotropic rays as fibration candidates."""

    if coefficient_bound < 1:
        raise ValueError("coefficient_bound must be positive")
    ample = tuple(int(value) for value in ample_vector)
    if len(ample) != lattice.dimension:
        raise ValueError("ample_vector has the wrong dimension")
    rays: set[tuple[int, ...]] = set()
    values = range(-coefficient_bound, coefficient_bound + 1)
    for raw in product(values, repeat=lattice.dimension):
        vector = tuple(raw)
        if not any(vector) or not _primitive_vector(vector):
            continue
        if lattice.norm(vector) or lattice.pair(vector, ample) <= 0:
            continue
        rays.add(vector)
    return tuple(sorted(rays))


@dataclass(frozen=True)
class BranchDivisor:
    """An even branch divisor on P1, represented by named geometric places."""

    places: tuple[str, ...]

    def __post_init__(self) -> None:
        if any(not place for place in self.places):
            raise ValueError("branch-place labels must be nonempty")
        if len(set(self.places)) != len(self.places):
            raise ValueError("branch places must be distinct")
        if len(self.places) % 2:
            raise ValueError("a quadratic cover of P1 has an even branch count")

    @classmethod
    def squarefree_polynomial(
        cls, degree: int, *, prefix: str
    ) -> "BranchDivisor":
        if degree < 1 or not prefix:
            raise ValueError("a positive degree and prefix are required")
        places = [f"{prefix}:root:{index}" for index in range(1, degree + 1)]
        if degree % 2:
            places.append("infinity")
        return cls(tuple(places))

    @property
    def genus(self) -> int:
        if not self.places:
            raise ValueError("an unramified split squareclass is not a connected cover")
        return (len(self.places) - 2) // 2

    def product(self, other: "BranchDivisor") -> "BranchDivisor":
        return BranchDivisor(
            tuple(sorted(set(self.places).symmetric_difference(other.places)))
        )

    def to_record(self) -> dict[str, object]:
        return {
            "places": list(self.places),
            "branch_count": len(self.places),
            "genus": self.genus,
        }


@dataclass(frozen=True)
class V4CoverDecomposition:
    first: BranchDivisor
    second: BranchDivisor

    def __post_init__(self) -> None:
        if self.first.places == self.second.places:
            raise ValueError("the two quadratic characters must be independent")
        if not self.third.places:
            raise ValueError("the product character must be nontrivial")

    @property
    def third(self) -> BranchDivisor:
        return self.first.product(self.second)

    @property
    def quotient_genera(self) -> tuple[int, int, int]:
        return self.first.genus, self.second.genus, self.third.genus

    @property
    def cover_genus(self) -> int:
        # Kani--Rosen dimension identity for a connected V4 cover.
        return sum(self.quotient_genera)

    def to_record(self) -> dict[str, object]:
        return {
            "first_character": self.first.to_record(),
            "second_character": self.second.to_record(),
            "product_character": self.third.to_record(),
            "quotient_genera": list(self.quotient_genera),
            "cover_genus": self.cover_genus,
            "jacobian_dimension_identity": (
                "dim J(C)=dim J(C1)+dim J(C2)+dim J(C3)"
            ),
            "claim_boundary": (
                "the decomposition organizes downstream Selmer/Chabauty work; "
                "it does not determine the rational points"
            ),
        }


@dataclass(frozen=True)
class TwistCharacter:
    name: str
    branch_divisor: BranchDivisor

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("a twist character needs a name")

    @property
    def base_change_genus(self) -> int:
        return self.branch_divisor.genus

    def rank_lower_bound_after_base_change(
        self, *, base_rank_lower_bound: int, twist_rank_lower_bound: int
    ) -> int:
        if min(base_rank_lower_bound, twist_rank_lower_bound) < 0:
            raise ValueError("rank lower bounds must be nonnegative")
        return base_rank_lower_bound + twist_rank_lower_bound

    def to_record(self) -> dict[str, object]:
        return {
            "name": self.name,
            "branch_divisor": self.branch_divisor.to_record(),
            "base_change_genus": self.base_change_genus,
            "search_principle": (
                "seek a new non-torsion section on the quadratic twist; "
                "division of an old section does not increase rank"
            ),
        }


@dataclass(frozen=True)
class ProjectiveRational:
    numerator: int
    denominator: int

    def __post_init__(self) -> None:
        if self.denominator == 0 and self.numerator == 0:
            raise ValueError("(0:0) is not a projective point")
        content = gcd(abs(self.numerator), abs(self.denominator))
        if content != 1:
            raise ValueError("projective coordinates must be primitive")
        if self.denominator < 0 or (
            self.denominator == 0 and self.numerator < 0
        ):
            raise ValueError("use the canonical projective sign")

    @classmethod
    def normalized(cls, numerator: int, denominator: int) -> "ProjectiveRational":
        numerator = int(numerator)
        denominator = int(denominator)
        if numerator == 0 and denominator == 0:
            raise ValueError("(0:0) is not a projective point")
        content = gcd(abs(numerator), abs(denominator))
        numerator //= content
        denominator //= content
        if denominator < 0 or (denominator == 0 and numerator < 0):
            numerator = -numerator
            denominator = -denominator
        return cls(numerator, denominator)

    def to_record(self) -> dict[str, int]:
        return {"numerator": self.numerator, "denominator": self.denominator}


@dataclass(frozen=True)
class ProjectivePadicBall:
    """A residue ball in either affine chart of P1(Q_p).

    In the affine chart the local coordinate is ``T=a/b`` and matching means
    ``a-residue*b = 0 mod p^k`` with ``b`` a unit.  In the infinity chart the
    coordinate is ``S=b/a=1/T`` and matching means
    ``b-residue*a = 0 mod p^k`` with ``a`` a unit.  Residue zero in the latter
    is a genuine p-adic neighborhood of infinity.
    """

    prime: int
    exponent: int
    chart: str
    residue: int = 0

    def __post_init__(self) -> None:
        if not _is_prime(self.prime):
            raise ValueError("p-adic balls require a prime")
        if self.exponent < 1:
            raise ValueError("the p-adic exponent must be positive")
        if self.chart not in {"affine", "infinity"}:
            raise ValueError("chart must be 'affine' or 'infinity'")
        object.__setattr__(self, "residue", self.residue % self.modulus)

    @property
    def modulus(self) -> int:
        return self.prime**self.exponent

    def matches(self, point: ProjectiveRational) -> bool:
        a_value, b_value = point.numerator, point.denominator
        if self.chart == "affine":
            if b_value % self.prime == 0:
                return False
            return (a_value - self.residue * b_value) % self.modulus == 0
        if a_value % self.prime == 0:
            return False
        return (b_value - self.residue * a_value) % self.modulus == 0

    @property
    def linear_form(self) -> tuple[int, int]:
        """Return ``(A,B)`` for ``A*a+B*b == 0 mod p^k``."""

        if self.chart == "affine":
            return 1, -self.residue
        return -self.residue, 1

    def to_record(self) -> dict[str, object]:
        return {
            "prime": self.prime,
            "exponent": self.exponent,
            "chart": self.chart,
            "residue": self.residue,
            "modulus": self.modulus,
            "linear_form_A_a_plus_B_b": list(self.linear_form),
            "interpretation": (
                "a projective local condition; exact local reduction must "
                "still be recomputed after specialization"
            ),
        }


def _crt_many(residues: Sequence[int], moduli: Sequence[int]) -> tuple[int, int]:
    if len(residues) != len(moduli) or not residues:
        raise ValueError("CRT needs equally many nonempty residues and moduli")
    residue = 0
    modulus = 1
    for next_residue, next_modulus in zip(residues, moduli, strict=True):
        if next_modulus <= 0 or gcd(modulus, next_modulus) != 1:
            raise ValueError("projective CRT moduli must be positive and coprime")
        step = (
            (int(next_residue) - residue)
            * pow(modulus, -1, next_modulus)
            % next_modulus
        )
        residue += modulus * step
        modulus *= next_modulus
        residue %= modulus
    return residue, modulus


@dataclass(frozen=True)
class ProjectiveCongruenceLattice:
    """The rank-two lattice cut out by projective local linear forms."""

    basis: tuple[tuple[int, int], tuple[int, int]]
    modulus: int
    index: int
    combined_linear_form: tuple[int, int]
    balls: tuple[ProjectivePadicBall, ...]

    def __post_init__(self) -> None:
        determinant = abs(
            self.basis[0][0] * self.basis[1][1]
            - self.basis[0][1] * self.basis[1][0]
        )
        if determinant != self.index or not determinant:
            raise ValueError("the projective congruence basis has the wrong index")
        for vector in self.basis:
            for ball in self.balls:
                coefficient_a, coefficient_b = ball.linear_form
                if (coefficient_a * vector[0] + coefficient_b * vector[1]) % ball.modulus:
                    raise ValueError("a lattice basis vector misses a local congruence")

    def contains(self, vector: Sequence[int]) -> bool:
        if len(vector) != 2:
            raise ValueError("projective lattice vectors have dimension two")
        return all(
            (ball.linear_form[0] * int(vector[0]) + ball.linear_form[1] * int(vector[1]))
            % ball.modulus
            == 0
            for ball in self.balls
        )

    def point_matches_charts(self, point: ProjectiveRational) -> bool:
        """Also enforce the unit condition selecting each affine chart."""

        return all(ball.matches(point) for ball in self.balls)

    def to_record(self) -> dict[str, object]:
        return {
            "basis": [list(vector) for vector in self.basis],
            "modulus": self.modulus,
            "index": self.index,
            "combined_linear_form": list(self.combined_linear_form),
            "balls": [ball.to_record() for ball in self.balls],
            "claim_boundary": (
                "lattice membership encodes the homogeneous congruences; "
                "primitive-vector and chart-unit checks remain separate"
            ),
        }


def projective_congruence_lattice(
    balls: Sequence[ProjectivePadicBall],
) -> ProjectiveCongruenceLattice:
    """Intersect one projective p-adic linear condition per prime exactly.

    Pairwise coprime moduli are combined coefficientwise by CRT into one
    congruence ``A*a+B*b=0 mod M``.  A Hermite-style basis is then constructed
    without floating-point lattice reduction.
    """

    if not balls:
        raise ValueError("at least one projective p-adic ball is required")
    if len({ball.prime for ball in balls}) != len(balls):
        raise ValueError("combine at most one already-compressed ball per prime")
    moduli = tuple(ball.modulus for ball in balls)
    coefficient_as = tuple(ball.linear_form[0] for ball in balls)
    coefficient_bs = tuple(ball.linear_form[1] for ball in balls)
    combined_a, modulus_a = _crt_many(coefficient_as, moduli)
    combined_b, modulus_b = _crt_many(coefficient_bs, moduli)
    if modulus_a != modulus_b:
        raise AssertionError("coefficientwise CRT returned inconsistent moduli")
    modulus = modulus_a

    divisor_a = gcd(combined_a, modulus)
    common = gcd(divisor_a, combined_b)
    first = (modulus // divisor_a, 0)
    second_denominator = divisor_a // common
    reduced_modulus = modulus // divisor_a
    if reduced_modulus == 1:
        second_numerator = 0
    else:
        second_numerator = (
            -(combined_b // common)
            * pow(combined_a // divisor_a, -1, reduced_modulus)
            % reduced_modulus
        )
    second = (second_numerator, second_denominator)
    return ProjectiveCongruenceLattice(
        basis=(first, second),
        modulus=modulus,
        index=modulus // common,
        combined_linear_form=(combined_a, combined_b),
        balls=tuple(balls),
    )


@dataclass(frozen=True)
class ExternalComputationTask:
    task_id: str
    lane: str
    exact_inputs: Mapping[str, object]
    engine_options: tuple[str, ...]
    required_outputs: tuple[str, ...]
    success_gate: str
    claim_boundary: str

    def __post_init__(self) -> None:
        if not self.task_id or not self.lane:
            raise ValueError("external tasks need identifiers and lanes")
        if not self.engine_options or not self.required_outputs:
            raise ValueError("external tasks need engines and required outputs")

    def to_record(self) -> dict[str, object]:
        return {
            "task_id": self.task_id,
            "lane": self.lane,
            "exact_inputs": dict(self.exact_inputs),
            "engine_options": list(self.engine_options),
            "required_outputs": list(self.required_outputs),
            "success_gate": self.success_gate,
            "claim_boundary": self.claim_boundary,
        }
