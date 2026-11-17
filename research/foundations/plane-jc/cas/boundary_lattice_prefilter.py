#!/usr/bin/env python3
"""Exact boundary-lattice prefilter for plane Newton candidates.

For a smooth projective compactification X with free Picard group and boundary
primes D_i, the columns of ``boundary_matrix`` are the classes [D_i].  Smith
normal form reads the divisor-localization sequence as

    units/k^* -> Z^{boundary primes} -> Pic(X) -> Pic(U) -> 0.

The checker is deliberately independent of Newton coefficient equations.  It
only audits a *complete* proposed boundary list.  Corners, dicritical primes,
or a selected subgraph are not complete inputs.

There is a dual input for a normal reconstruction open containing a dense
affine UFD core.  The rows of ``factorial_core_valuation_matrix`` are all codimension-one
primes outside the core and the columns are a basis of core units modulo
constants.  Its kernel is the unit lattice and its cokernel is the Weil class
group.  ``torus_core_valuation_matrix`` remains a supported special-case key.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


@dataclass(frozen=True)
class LocalizationInvariants:
    """Free and torsion invariants read from a boundary-class matrix."""

    picard_rank: int
    boundary_count: int
    matrix_rank: int
    unit_rank: int
    picard_free_rank: int
    picard_torsion: tuple[int, ...]
    smith_diagonal: tuple[int, ...]

    def passes(self, expected_unit_rank: int) -> bool:
        """Test Pic(U)=0 and the requested rank of O(U)^*/k^*."""

        return (
            self.unit_rank == expected_unit_rank
            and self.picard_free_rank == 0
            and not self.picard_torsion
        )


@dataclass(frozen=True)
class TorusCoreInvariants:
    """Units and Weil-class invariants of a complete class-trivial core."""

    boundary_count: int
    character_rank: int
    matrix_rank: int
    unit_rank: int
    class_free_rank: int
    class_torsion: tuple[int, ...]
    smith_diagonal: tuple[int, ...]

    def passes_affine_space_gate(self) -> bool:
        """Test the unit/class-group conditions forced by affine space."""

        return (
            self.unit_rank == 0
            and self.class_free_rank == 0
            and not self.class_torsion
        )


@dataclass(frozen=True)
class MultipleFiberInvariants:
    """Vertical class group forced by a certified multiple-fiber packet."""

    multiplicities: tuple[int, ...]
    fiber_count: int
    class_torsion: tuple[int, ...]
    smith_diagonal: tuple[int, ...]
    relation_matrix: tuple[tuple[int, ...], ...]
    full_class_group_certified: bool

    @property
    def multiplicity(self) -> int | None:
        """Return the common multiplicity, or ``None`` for a mixed packet."""

        if len(set(self.multiplicities)) == 1:
            return self.multiplicities[0]
        return None

    def passes_affine_space_gate(self) -> bool:
        """A nontrivial forced vertical subgroup rejects affine space."""

        return not self.class_torsion


@dataclass(frozen=True)
class BoundaryConfiguration:
    """A complete boundary in a fixed total-transform Picard basis.

    ``class_matrix`` has boundary-prime classes as columns.  ``intersections``
    records the unordered pairs of boundary primes which currently meet.
    A one-parent blowup is understood to take place at a smooth point of that
    component away from all recorded boundary crossings.
    """

    class_matrix: sp.Matrix
    names: tuple[str, ...]
    intersections: frozenset[frozenset[str]]

    def __post_init__(self) -> None:
        if self.class_matrix.cols != len(self.names):
            raise ValueError("boundary names must index the class-matrix columns")
        if len(set(self.names)) != len(self.names):
            raise ValueError("boundary-prime names must be distinct")
        known = set(self.names)
        for edge in self.intersections:
            if len(edge) != 2 or not set(edge) <= known:
                raise ValueError(f"invalid boundary intersection {set(edge)}")

    def blow_up(
        self,
        center: Sequence[str],
        exceptional: str | None = None,
    ) -> "BoundaryConfiguration":
        """Blow up a smooth boundary point with one or two boundary parents."""

        parents = tuple(center)
        if not 1 <= len(parents) <= 2:
            raise ValueError("a boundary center meets one or two primes")
        if len(set(parents)) != len(parents):
            raise ValueError("a boundary center cannot repeat a prime")
        missing = [name for name in parents if name not in self.names]
        if missing:
            raise ValueError(f"unknown boundary components {missing}")
        crossing = frozenset(parents)
        if len(parents) == 2 and crossing not in self.intersections:
            raise ValueError(f"components {parents} do not meet in the current boundary")

        name = exceptional or f"E{self.class_matrix.rows}"
        if name in self.names:
            raise ValueError(f"exceptional name {name!r} is already in use")

        old = self.class_matrix
        matrix = sp.zeros(old.rows + 1, old.cols + 1)
        matrix[: old.rows, : old.cols] = old
        for parent in parents:
            matrix[old.rows, self.names.index(parent)] = -1
        matrix[old.rows, old.cols] = 1

        intersections = set(self.intersections)
        if len(parents) == 2:
            intersections.remove(crossing)
        for parent in parents:
            intersections.add(frozenset((parent, name)))
        return BoundaryConfiguration(
            class_matrix=matrix,
            names=self.names + (name,),
            intersections=frozenset(intersections),
        )

    def fill_boundary_component(self, name: str) -> "BoundaryConfiguration":
        """Put one temporary boundary divisor back into the open surface."""

        if name not in self.names:
            raise ValueError(f"unknown boundary component {name!r}")
        index = self.names.index(name)
        matrix = self.class_matrix.copy()
        matrix.col_del(index)
        return BoundaryConfiguration(
            class_matrix=matrix,
            names=tuple(
                boundary_name
                for boundary_name in self.names
                if boundary_name != name
            ),
            intersections=frozenset(
                edge for edge in self.intersections if name not in edge
            ),
        )


def standard_completion(chart: str) -> tuple[BoundaryConfiguration, sp.Matrix, int]:
    """Return a standard completion, its intersection form, and unit rank.

    Supported charts are ``A2`` compactified by ``P^2``, ``GmA1``
    compactified by ``P^1 x P^1``, and ``GmA1_Fn`` compactified by the
    Hirzebruch surface whose retained infinity section has self-intersection
    ``n``.
    """

    if chart == "A2":
        return (
            BoundaryConfiguration(sp.Matrix([[1]]), ("L",), frozenset()),
            sp.Matrix([[1]]),
            0,
        )
    if chart == "GmA1":
        names = ("X0", "Xinf", "Yinf")
        intersections = frozenset(
            (frozenset(("X0", "Yinf")), frozenset(("Xinf", "Yinf")))
        )
        return (
            BoundaryConfiguration(
                sp.Matrix([[1, 1, 0], [0, 0, 1]]),
                names,
                intersections,
            ),
            sp.Matrix([[0, 1], [1, 0]]),
            1,
        )
    hirzebruch = re.fullmatch(r"GmA1_F([0-9]+)", chart)
    if hirzebruch:
        degree = int(hirzebruch.group(1))
        names = ("X0", "Xinf", "Yinf")
        intersections = frozenset(
            (frozenset(("X0", "Yinf")), frozenset(("Xinf", "Yinf")))
        )
        return (
            BoundaryConfiguration(
                sp.Matrix([[1, 1, 0], [0, 0, 1]]),
                names,
                intersections,
            ),
            sp.Matrix([[0, 1], [1, degree]]),
            1,
        )
    raise ValueError(
        "chart must be 'A2', 'GmA1', or 'GmA1_Fn' for n>=0"
    )


def boundary_intersection_matrix(
    configuration: BoundaryConfiguration,
    initial_intersection_form: sp.Matrix,
) -> sp.Matrix:
    """Compute the full boundary intersection matrix after point blowups."""

    initial_rank = initial_intersection_form.rows
    if initial_intersection_form.cols != initial_rank:
        raise ValueError("the initial intersection form must be square")
    exceptional_count = configuration.class_matrix.rows - initial_rank
    if exceptional_count < 0:
        raise ValueError("the class matrix is smaller than its initial completion")
    form = sp.diag(initial_intersection_form, *([-1] * exceptional_count))
    return configuration.class_matrix.T * form * configuration.class_matrix


def localization_invariants(
    boundary_matrix: Sequence[Sequence[int]],
) -> LocalizationInvariants:
    """Compute exact localization invariants from integral class columns."""

    matrix = sp.Matrix(boundary_matrix)
    if matrix.rows == 0 or matrix.cols == 0:
        raise ValueError("the boundary matrix must have at least one row and column")
    if any(value.is_Integer is not True for value in matrix):
        raise ValueError("boundary classes must be integral")

    smith = smith_normal_form(matrix, domain=ZZ)
    full_diagonal = tuple(
        abs(int(smith[index, index]))
        for index in range(min(smith.rows, smith.cols))
    )
    nonzero_diagonal = tuple(value for value in full_diagonal if value != 0)
    rank = len(nonzero_diagonal)
    return LocalizationInvariants(
        picard_rank=matrix.rows,
        boundary_count=matrix.cols,
        matrix_rank=rank,
        unit_rank=matrix.cols - rank,
        picard_free_rank=matrix.rows - rank,
        picard_torsion=tuple(value for value in nonzero_diagonal if value > 1),
        smith_diagonal=full_diagonal,
    )


def torus_core_invariants(
    valuation_matrix: Sequence[Sequence[int]],
) -> TorusCoreInvariants:
    """Read units and ``Cl(U)`` from a complete torus-core valuation map.

    The matrix has one row for every codimension-one prime in ``U - T`` and
    one column for every character in a basis of ``T = G_m^n``.  This
    function checks only the integer lattice.  Normality, the torus-chart
    identification, and completeness of the codimension-one boundary are
    proof-bearing input conditions which no integer matrix can certify.
    """

    profile = localization_invariants(valuation_matrix)
    return TorusCoreInvariants(
        boundary_count=profile.picard_rank,
        character_rank=profile.boundary_count,
        matrix_rank=profile.matrix_rank,
        unit_rank=profile.unit_rank,
        class_free_rank=profile.picard_free_rank,
        class_torsion=profile.picard_torsion,
        smith_diagonal=profile.smith_diagonal,
    )


def factorial_core_invariants(
    valuation_matrix: Sequence[Sequence[int]],
) -> TorusCoreInvariants:
    """Alias for the affine-UFD core case, which has trivial class group."""

    return torus_core_invariants(valuation_matrix)


def multiple_fiber_invariants(
    multiplicity: int | Sequence[int],
    fiber_count: int | None = None,
    *,
    reduced_sum_principal: bool,
    generic_unit_rank: int,
    generic_class_group_trivial: bool = False,
    other_vertical_classes_trivial: bool = False,
) -> MultipleFiberInvariants:
    """Compile the exact vertical Smith block for a multiple-fiber packet.

    Geometry, not the integer matrix, must certify that the reduced fiber sum
    is principal and that the generic-fiber unit rank is one.  Under those
    hypotheses the relation columns are ``diag(m_1,...,m_r)`` and the
    all-ones column, giving ``(direct_sum Z/m_i)/<(1,...,1)>``.  Passing an
    integer and ``fiber_count`` recovers the common-multiplicity shortcut
    ``(Z/m)^(r-1)``.  Trivial generic class group and triviality of every
    other vertical prime class promote this subgroup to the full class group.
    Irreducibility alone does not imply the latter over a general base curve.
    """

    if isinstance(multiplicity, int):
        if fiber_count is None:
            raise ValueError("common multiplicity input requires a fiber count")
        multiplicities = (multiplicity,) * fiber_count
    else:
        multiplicities = tuple(int(value) for value in multiplicity)
        if fiber_count is not None and fiber_count != len(multiplicities):
            raise ValueError("fiber count does not match the multiplicity list")
        fiber_count = len(multiplicities)
    if fiber_count < 1:
        raise ValueError("a multiple-fiber packet needs at least one fiber")
    if any(value < 2 for value in multiplicities):
        raise ValueError("multiple-fiber multiplicities must be at least two")
    if not reduced_sum_principal:
        raise ValueError("the reduced fiber sum must be certified principal")
    if generic_unit_rank != 1:
        raise ValueError("the generic-fiber unit rank must be certified as one")

    relations = sp.Matrix.hstack(
        sp.diag(*multiplicities),
        sp.ones(fiber_count, 1),
    )
    profile = localization_invariants(relations.tolist())
    expected_order = math.prod(multiplicities) // math.lcm(*multiplicities)
    actual_order = math.prod(profile.picard_torsion)
    if profile.picard_free_rank != 0 or actual_order != expected_order:
        raise ArithmeticError("unexpected multiple-fiber Smith profile")
    primes = sorted(
        {
            int(prime)
            for value in multiplicities
            for prime in sp.factorint(value)
        }
    )
    for prime in primes:
        expected_exponents = sorted(
            int(sp.factorint(value).get(prime, 0)) for value in multiplicities
        )[:-1]
        actual_exponents = sorted(
            exponent
            for value in profile.picard_torsion
            if (exponent := int(sp.factorint(value).get(prime, 0)))
        )
        if actual_exponents != [value for value in expected_exponents if value]:
            raise ArithmeticError("unexpected multiple-fiber primary profile")
    return MultipleFiberInvariants(
        multiplicities=multiplicities,
        fiber_count=fiber_count,
        class_torsion=profile.picard_torsion,
        smith_diagonal=profile.smith_diagonal,
        relation_matrix=tuple(
            tuple(int(relations[row, column]) for column in range(relations.cols))
            for row in range(relations.rows)
        ),
        full_class_group_certified=(
            generic_class_group_trivial and other_vertical_classes_trivial
        ),
    )


def core_class_order(
    valuation_matrix: Sequence[Sequence[int]],
    class_vector: Sequence[int],
) -> int | None:
    """Return the order of a boundary class in the core-localization cokernel.

    The return value is a positive integer, with ``1`` meaning trivial.
    ``None`` means that the class has infinite order.  If ``L`` is the column
    lattice of ``V`` and ``L' = L + Z*delta`` has the same rank, the order is
    ``[L':L]``.  Determinantal divisors compute this as the quotient of the
    products of the nonzero Smith factors of ``V`` and ``[V|delta]``.
    """

    matrix = sp.Matrix(valuation_matrix)
    vector = sp.Matrix(tuple(class_vector))
    if vector.cols != 1 or vector.rows != matrix.rows:
        raise ValueError("a class vector needs one entry per boundary prime")
    if any(value.is_Integer is not True for value in vector):
        raise ValueError("boundary class coordinates must be integral")
    if matrix.cols == 0:
        return 1 if all(value == 0 for value in vector) else None
    augmented = matrix.row_join(vector)
    base = torus_core_invariants(matrix.tolist())
    enlarged = torus_core_invariants(augmented.tolist())
    if enlarged.matrix_rank > base.matrix_rank:
        return None
    base_divisor = math.prod(
        value for value in base.smith_diagonal if value != 0
    )
    enlarged_divisor = math.prod(
        value for value in enlarged.smith_diagonal if value != 0
    )
    if enlarged_divisor == 0 or base_divisor % enlarged_divisor:
        raise ArithmeticError("determinantal-divisor quotient is not integral")
    return base_divisor // enlarged_divisor


def presented_core_invariants(
    unit_valuation_matrix: Sequence[Sequence[int]],
    core_relation_matrix: Sequence[Sequence[int]],
    relation_boundary_corrections: Sequence[Sequence[int]],
) -> dict[str, object]:
    """Compile a nontrivial core class presentation by a lifted block map."""

    valuation = sp.Matrix(unit_valuation_matrix)
    relations = sp.Matrix(core_relation_matrix)
    corrections = sp.Matrix(relation_boundary_corrections)
    if valuation.rows == 0:
        raise ValueError("a presented core needs at least one boundary row")
    if relations.rows == 0:
        raise ValueError("a presented core needs at least one class-generator row")
    if corrections.shape != (valuation.rows, relations.cols):
        raise ValueError(
            "relation-boundary corrections must have shape boundary x relations"
        )
    block = valuation.row_join(corrections).col_join(
        sp.zeros(relations.rows, valuation.cols).row_join(relations)
    )
    if block.cols == 0:
        class_free_rank = block.rows
        class_torsion: tuple[int, ...] = ()
        smith_diagonal: tuple[int, ...] = ()
    else:
        profile = localization_invariants(block.tolist())
        class_free_rank = profile.picard_free_rank
        class_torsion = profile.picard_torsion
        smith_diagonal = profile.smith_diagonal
    return {
        "unit_rank": valuation.cols - valuation.rank(),
        "class_free_rank": class_free_rank,
        "class_torsion": class_torsion,
        "smith_diagonal": smith_diagonal,
        "presentation_matrix": tuple(
            tuple(int(block[row, column]) for column in range(block.cols))
            for row in range(block.rows)
        ),
    }


def boundary_blowup_matrix(
    centers: Sequence[Sequence[str]],
) -> tuple[sp.Matrix, tuple[str, ...]]:
    """Construct boundary classes for blowups of P^2 along its boundary.

    The initial boundary is ``L`` with class H.  Each center is the list of
    one or two current boundary components containing the blowup point.  The
    returned rows use the total-transform basis H,E1,...,Es; columns are the
    final strict transforms followed by exceptional primes in creation order.
    """

    configuration, _, _ = standard_completion("A2")
    for index, raw_center in enumerate(centers, start=1):
        try:
            configuration = configuration.blow_up(raw_center, f"E{index}")
        except ValueError as error:
            raise ValueError(f"blowup {index}: {error}") from error
    return configuration.class_matrix, configuration.names


def _format_report(name: str, matrix: sp.Matrix, expected_unit_rank: int) -> str:
    invariants = localization_invariants(matrix.tolist())
    status = "PASS" if invariants.passes(expected_unit_rank) else "REJECT"
    return (
        f"{status} {name}: shape={matrix.rows}x{matrix.cols}, "
        f"SNF={invariants.smith_diagonal}, units_rank={invariants.unit_rank}, "
        f"Pic_free_rank={invariants.picard_free_rank}, "
        f"Pic_torsion={invariants.picard_torsion}"
    )


def _format_factorial_core_report(
    name: str,
    matrix: sp.Matrix,
    class_vectors: dict[str, Sequence[int]] | None = None,
) -> str:
    invariants = factorial_core_invariants(matrix.tolist())
    status = "PASS" if invariants.passes_affine_space_gate() else "REJECT"
    report = (
        f"{status} {name}: shape={matrix.rows}x{matrix.cols}, "
        f"SNF={invariants.smith_diagonal}, units_rank={invariants.unit_rank}, "
        f"Cl_free_rank={invariants.class_free_rank}, "
        f"Cl_torsion={invariants.class_torsion}"
    )
    if class_vectors:
        orders = {
            class_name: core_class_order(matrix.tolist(), vector)
            for class_name, vector in class_vectors.items()
        }
        report += f", class_orders={orders}"
    return report


def _format_torus_core_report(name: str, matrix: sp.Matrix) -> str:
    """Backward-compatible formatter for the torus special case."""

    return _format_factorial_core_report(name, matrix)


def _format_multiple_fiber_report(
    name: str,
    profile: MultipleFiberInvariants,
) -> str:
    status = "PASS" if profile.passes_affine_space_gate() else "REJECT"
    scope = "full_Cl" if profile.full_class_group_certified else "vertical_subgroup"
    multiplicities = (
        str(profile.multiplicity)
        if profile.multiplicity is not None
        else str(profile.multiplicities)
    )
    return (
        f"{status} {name}: multiple_fibers={profile.fiber_count}, "
        f"multiplicities={multiplicities}, "
        f"SNF={profile.smith_diagonal}, "
        f"Cl_torsion={profile.class_torsion}, scope={scope}"
    )


def run_regressions() -> None:
    """Check affine-plane, Laurent-chart, and obstruction fixtures."""

    base, _ = boundary_blowup_matrix([])
    once, _ = boundary_blowup_matrix([("L",)])
    crossing, names = boundary_blowup_matrix([("L",), ("L", "E1")])
    assert names == ("L", "E1", "E2")
    assert base == sp.Matrix([[1]])
    assert once == sp.Matrix([[1, 0], [-1, 1]])
    assert crossing == sp.Matrix([[1, 0, 0], [-1, 1, 0], [-1, -1, 1]])
    try:
        boundary_blowup_matrix([("L",), ("L", "E1"), ("L", "E1")])
    except ValueError as error:
        assert "do not meet" in str(error)
    else:
        raise AssertionError("a stale boundary crossing must be rejected")

    marked_n2 = sp.Matrix([[1, 1], [1, 1]])
    marked_n3 = sp.Matrix([[2, 1], [1, 1]])
    wrong_multiplicity = sp.Matrix([[2, 0], [-2, 1]])
    laurent_chart = sp.Matrix([[1, 1, 0], [0, 0, 1]])
    hirzebruch, hirzebruch_form, hirzebruch_units = standard_completion(
        "GmA1_F4"
    )

    for matrix in (base, once, crossing, marked_n3):
        assert localization_invariants(matrix.tolist()).passes(expected_unit_rank=0)
    assert not localization_invariants(marked_n2.tolist()).passes(0)
    assert not localization_invariants(wrong_multiplicity.tolist()).passes(0)
    assert localization_invariants(laurent_chart.tolist()).passes(expected_unit_rank=1)
    assert not localization_invariants(laurent_chart.tolist()).passes(0)
    assert hirzebruch.class_matrix == laurent_chart
    assert hirzebruch_form == sp.Matrix([[0, 1], [1, 4]])
    assert hirzebruch_units == 1
    assert boundary_intersection_matrix(
        hirzebruch, hirzebruch_form
    ) == sp.Matrix(
        [
            [0, 0, 1],
            [0, 0, 1],
            [1, 1, 4],
        ]
    )
    filled_hirzebruch = hirzebruch.fill_boundary_component("X0")
    assert filled_hirzebruch.names == ("Xinf", "Yinf")
    assert localization_invariants(
        filled_hirzebruch.class_matrix.tolist()
    ).passes(expected_unit_rank=0)

    torus_affine_plane = sp.eye(2)
    torus_free_unit = sp.Matrix([[1, 0]])
    assert torus_core_invariants(
        torus_affine_plane.tolist()
    ).passes_affine_space_gate()
    free_unit_profile = torus_core_invariants(torus_free_unit.tolist())
    assert free_unit_profile.unit_rank == 1
    assert not free_unit_profile.passes_affine_space_gate()
    assert core_class_order(torus_affine_plane.tolist(), (1, 0)) == 1
    assert core_class_order([[1], [0]], (0, 1)) is None
    split_core = presented_core_invariants([[2]], [[2]], [[0]])
    nonsplit_core = presented_core_invariants([[2]], [[2]], [[1]])
    free_core = presented_core_invariants([[]], [[]], [[]])
    assert split_core["class_torsion"] == (2, 2)
    assert nonsplit_core["class_torsion"] == (4,)
    assert free_core["unit_rank"] == 0
    assert free_core["class_free_rank"] == 2
    assert free_core["smith_diagonal"] == ()
    assert core_class_order(free_core["presentation_matrix"], (0, 0)) == 1
    assert core_class_order(free_core["presentation_matrix"], (1, 0)) is None
    assert nonsplit_core["presentation_matrix"] == ((2, 1), (0, 2))
    assert core_class_order(
        nonsplit_core["presentation_matrix"], (0, 1)
    ) == 4
    for degree in (3, 5, 6, 7):
        balanced = sp.Matrix(
            [[1, 0], [degree - 2, degree - 1]]
        )
        balanced_profile = torus_core_invariants(balanced.tolist())
        assert balanced_profile.unit_rank == 0
        assert balanced_profile.class_free_rank == 0
        assert balanced_profile.class_torsion == (degree - 1,)
        assert core_class_order(balanced.tolist(), (0, 1)) == degree - 1
        assert core_class_order(balanced.tolist(), (1, degree - 2)) == 1
        assert not balanced_profile.passes_affine_space_gate()
    double_fibers = multiple_fiber_invariants(
        2,
        3,
        reduced_sum_principal=True,
        generic_unit_rank=1,
        generic_class_group_trivial=True,
        other_vertical_classes_trivial=True,
    )
    triple_fibers = multiple_fiber_invariants(
        3,
        4,
        reduced_sum_principal=True,
        generic_unit_rank=1,
    )
    assert double_fibers.class_torsion == (2, 2)
    assert double_fibers.full_class_group_certified
    assert not double_fibers.passes_affine_space_gate()
    assert triple_fibers.class_torsion == (3, 3, 3)
    assert not triple_fibers.full_class_group_certified
    mixed_fibers = multiple_fiber_invariants(
        (4, 6, 9),
        reduced_sum_principal=True,
        generic_unit_rank=1,
    )
    assert mixed_fibers.class_torsion == (6,)
    assert math.prod(mixed_fibers.class_torsion) == (
        math.prod(mixed_fibers.multiplicities)
        // math.lcm(*mixed_fibers.multiplicities)
    )
    coprime_fibers = multiple_fiber_invariants(
        (4, 9, 25),
        reduced_sum_principal=True,
        generic_unit_rank=1,
    )
    assert not coprime_fibers.class_torsion
    assert coprime_fibers.passes_affine_space_gate()
    try:
        multiple_fiber_invariants(
            2,
            3,
            reduced_sum_principal=False,
            generic_unit_rank=1,
        )
    except ValueError as error:
        assert "certified principal" in str(error)
    else:
        raise AssertionError("uncertified reduced-sum relations must be rejected")

    print(_format_report("A2 after two boundary blowups", crossing, 0))
    print(_format_report("marked-point n=2 matrix", marked_n2, 0))
    print(_format_report("marked-point n=3 matrix", marked_n3, 0))
    print(_format_report("incorrect doubled boundary class", wrong_multiplicity, 0))
    print(_format_report("Gm x A1 standard completion", laurent_chart, 1))
    print(_format_report("Gm x A1 Hirzebruch F4 completion", hirzebruch.class_matrix, 1))
    print(_format_torus_core_report("unimodular Gm^2 core", torus_affine_plane))
    print(_format_torus_core_report("one-boundary Gm^2 core", torus_free_unit))
    print(
        "REJECT presented-core nonsplit extension: "
        f"SNF={nonsplit_core['smith_diagonal']}, "
        "Cl_torsion=(4,), class_orders={'G': 4}"
    )
    print(_format_multiple_fiber_report("three double fibers", double_fibers))
    print(_format_multiple_fiber_report("four triple fibers", triple_fibers))
    print(_format_multiple_fiber_report("mixed fibers 4,6,9", mixed_fibers))
    print(_format_multiple_fiber_report("coprime fibers 4,9,25", coprime_fibers))
    for degree in (3, 5, 6, 7):
        print(
            _format_factorial_core_report(
                f"balanced wild core N={degree}",
                sp.Matrix(
                    [[1, 0], [degree - 2, degree - 1]]
                ),
                {"L1": (0, 1), "div(x)": (1, degree - 2)},
            )
        )
    print("PASS boundary-lattice regressions: chart-aware localization invariants agree")


def check_json(path: Path) -> int:
    """Check a user-supplied matrix or boundary-blowup description."""

    data = json.loads(path.read_text())
    name = str(data.get("name", path.name))
    if "multiple_fiber_profile" in data:
        packet = data["multiple_fiber_profile"]
        multiplicity_input = (
            tuple(int(value) for value in packet["multiplicities"])
            if "multiplicities" in packet
            else int(packet["multiplicity"])
        )
        fiber_count = (
            None
            if "multiplicities" in packet
            else int(packet["fiber_count"])
        )
        profile = multiple_fiber_invariants(
            multiplicity_input,
            fiber_count,
            reduced_sum_principal=bool(packet["reduced_sum_principal"]),
            generic_unit_rank=int(packet["generic_unit_rank"]),
            generic_class_group_trivial=bool(
                packet.get("generic_class_group_trivial", False)
            ),
            other_vertical_classes_trivial=bool(
                packet.get("other_vertical_classes_trivial", False)
            ),
        )
        report = _format_multiple_fiber_report(name, profile)
    elif "factorial_core_valuation_matrix" in data:
        matrix = sp.Matrix(data["factorial_core_valuation_matrix"])
        report = _format_factorial_core_report(
            name, matrix, data.get("core_class_vectors")
        )
    elif "torus_core_valuation_matrix" in data:
        matrix = sp.Matrix(data["torus_core_valuation_matrix"])
        report = _format_factorial_core_report(
            name, matrix, data.get("core_class_vectors")
        )
    elif "boundary_matrix" in data:
        expected_unit_rank = int(data["expected_unit_rank"])
        matrix = sp.Matrix(data["boundary_matrix"])
        report = _format_report(name, matrix, expected_unit_rank)
    elif "boundary_blowups" in data:
        expected_unit_rank = int(data["expected_unit_rank"])
        matrix, _ = boundary_blowup_matrix(data["boundary_blowups"])
        report = _format_report(name, matrix, expected_unit_rank)
    else:
        raise ValueError(
            "JSON needs boundary_matrix, boundary_blowups, or "
            "factorial_core_valuation_matrix, or multiple_fiber_profile"
        )
    print(report)
    return 0 if report.startswith("PASS") else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        help="optional JSON candidate; without it, run exact regressions",
    )
    args = parser.parse_args()
    if args.input is None:
        run_regressions()
        return 0
    return check_json(args.input)


if __name__ == "__main__":
    raise SystemExit(main())
