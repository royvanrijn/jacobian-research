#!/usr/bin/env python3
"""Exact boundary-lattice prefilter for plane Newton candidates.

For a smooth projective compactification X with free Picard group and boundary
primes D_i, the columns of ``boundary_matrix`` are the classes [D_i].  Smith
normal form reads the divisor-localization sequence as

    units/k^* -> Z^{boundary primes} -> Pic(X) -> Pic(U) -> 0.

The checker is deliberately independent of Newton coefficient equations.  It
only audits a *complete* proposed boundary list.  Corners, dicritical primes,
or a selected subgraph are not complete inputs.
"""

from __future__ import annotations

import argparse
import json
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

    print(_format_report("A2 after two boundary blowups", crossing, 0))
    print(_format_report("marked-point n=2 matrix", marked_n2, 0))
    print(_format_report("marked-point n=3 matrix", marked_n3, 0))
    print(_format_report("incorrect doubled boundary class", wrong_multiplicity, 0))
    print(_format_report("Gm x A1 standard completion", laurent_chart, 1))
    print(_format_report("Gm x A1 Hirzebruch F4 completion", hirzebruch.class_matrix, 1))
    print("PASS boundary-lattice regressions: chart-aware localization invariants agree")


def check_json(path: Path) -> int:
    """Check a user-supplied matrix or boundary-blowup description."""

    data = json.loads(path.read_text())
    name = str(data.get("name", path.name))
    expected_unit_rank = int(data["expected_unit_rank"])
    if "boundary_matrix" in data:
        matrix = sp.Matrix(data["boundary_matrix"])
    elif "boundary_blowups" in data:
        matrix, _ = boundary_blowup_matrix(data["boundary_blowups"])
    else:
        raise ValueError("JSON needs boundary_matrix or boundary_blowups")
    report = _format_report(name, matrix, expected_unit_rank)
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
