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

    components: dict[str, list[int]] = {"L": [1]}
    order = ["L"]
    intersections: set[frozenset[str]] = set()
    for index, raw_center in enumerate(centers, start=1):
        center = tuple(raw_center)
        if not 1 <= len(center) <= 2:
            raise ValueError(f"blowup {index}: a boundary center meets one or two primes")
        if len(set(center)) != len(center):
            raise ValueError(f"blowup {index}: repeated component in center")
        missing = [name for name in center if name not in components]
        if missing:
            raise ValueError(f"blowup {index}: unknown components {missing}")
        if len(center) == 2 and frozenset(center) not in intersections:
            raise ValueError(
                f"blowup {index}: components {center} do not meet in the current boundary"
            )

        for vector in components.values():
            vector.append(0)
        for name in center:
            components[name][-1] -= 1

        exceptional = f"E{index}"
        components[exceptional] = [0] * index + [1]
        order.append(exceptional)
        if len(center) == 2:
            intersections.remove(frozenset(center))
        for name in center:
            intersections.add(frozenset((name, exceptional)))

    matrix = sp.Matrix.hstack(*(sp.Matrix(components[name]) for name in order))
    return matrix, tuple(order)


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

    for matrix in (base, once, crossing, marked_n3):
        assert localization_invariants(matrix.tolist()).passes(expected_unit_rank=0)
    assert not localization_invariants(marked_n2.tolist()).passes(0)
    assert not localization_invariants(wrong_multiplicity.tolist()).passes(0)
    assert localization_invariants(laurent_chart.tolist()).passes(expected_unit_rank=1)
    assert not localization_invariants(laurent_chart.tolist()).passes(0)

    print(_format_report("A2 after two boundary blowups", crossing, 0))
    print(_format_report("marked-point n=2 matrix", marked_n2, 0))
    print(_format_report("marked-point n=3 matrix", marked_n3, 0))
    print(_format_report("incorrect doubled boundary class", wrong_multiplicity, 0))
    print(_format_report("Gm x A1 standard completion", laurent_chart, 1))
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
