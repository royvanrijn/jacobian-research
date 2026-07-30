#!/usr/bin/env python3
"""Exact filtered-syzygy frontier for the universal cubic cotangent module.

The universal and central unit-pruned presentations agree through collision
order two, but this does not imply that their syzygies agree.  This checker
computes the central minimal resolution and then applies the *unchanged*
central syzygies of the 6-by-25 presentation to the universal matrix.
Their exact remainders modulo the central image are retained as the first
finite obstruction to a coefficient-independent strictness proof.

The calculation does not decide universal cotangent saturation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

import verify_cubic_symbol_double_saturation as cubic  # noqa: E402
from research_universal_cubic_quartic_kernel_saturation import (  # noqa: E402
    unit_pruned_differential_relations,
    universal_tensor,
)


DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "universal_cubic_filtered_syzygy_frontier.json"
)


def singular_polynomial(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> str:
    """Serialize an expanded rational polynomial without ambiguous division."""

    terms: list[str] = []
    for monomial, coefficient in sp.Poly(
        sp.expand(expression), *variables
    ).terms():
        rational = sp.Rational(coefficient)
        scalar = (
            f"({rational.p}/{rational.q})"
            if rational.q != 1
            else str(rational.p)
        )
        monomial_string = "*".join(
            (
                f"{variable}^{degree}"
                if degree != 1
                else str(variable)
            )
            for variable, degree in zip(variables, monomial)
            if degree
        )
        terms.append(
            scalar + (f"*{monomial_string}" if monomial_string else "")
        )
    return ("+".join(terms).replace("+-", "-")) or "0"


def singular_module(
    columns: list[list[sp.Expr]],
    variables: tuple[sp.Symbol, ...],
) -> str:
    return ",".join(
        "["
        + ",".join(
            singular_polynomial(entry, variables) for entry in column
        )
        + "]"
        for column in columns
    )


def collision_degree(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> int | None:
    if expression == 0:
        return None
    return min(
        sum(monomial[:3])
        for monomial, _coefficient in sp.Poly(
            sp.expand(expression), *variables
        ).terms()
    )


def entrywise_filtration_audit(
    universal_columns: list[list[sp.Expr]],
    central_columns: list[list[sp.Expr]],
    variables: tuple[sp.Symbol, ...],
) -> dict[str, Any]:
    counts: dict[tuple[int | None, int | None], int] = {}
    violations: list[dict[str, int]] = []
    for column_index, (universal_column, central_column) in enumerate(
        zip(universal_columns, central_columns)
    ):
        for row_index, (universal_entry, central_entry) in enumerate(
            zip(universal_column, central_column)
        ):
            difference = sp.expand(universal_entry - central_entry)
            central_degree = collision_degree(
                central_entry, variables
            )
            difference_degree = collision_degree(
                difference, variables
            )
            counts[(central_degree, difference_degree)] = (
                counts.get((central_degree, difference_degree), 0) + 1
            )
            if (
                central_degree is not None
                and difference_degree is not None
                and difference_degree <= central_degree
            ):
                violations.append(
                    {
                        "row": row_index,
                        "column": column_index,
                        "central_collision_degree": central_degree,
                        "difference_collision_degree": difference_degree,
                    }
                )
    return {
        "strict_raise_on_entries_with_nonzero_central_part": not violations,
        "violations": violations,
        "degree_pair_counts": [
            {
                "central_collision_degree": central_degree,
                "difference_collision_degree": difference_degree,
                "entry_count": count,
            }
            for (central_degree, difference_degree), count in sorted(
                counts.items(),
                key=lambda item: (
                    -1 if item[0][0] is None else item[0][0],
                    -1 if item[0][1] is None else item[0][1],
                ),
            )
        ],
    }


def marker(stdout: str, name: str) -> int:
    prefix = f"@@{name}="
    values = [
        line[len(prefix) :].strip()
        for line in stdout.splitlines()
        if line.startswith(prefix)
    ]
    assert len(values) == 1, (name, values, stdout)
    return int(values[0])


def section(stdout: str, name: str) -> list[str]:
    begin = f"@@BEGIN_{name}"
    end = f"@@END_{name}"
    lines = stdout.splitlines()
    start = lines.index(begin) + 1
    stop = lines.index(end)
    return [line.rstrip() for line in lines[start:stop] if line.strip()]


def remainder_columns(lines: list[str], rank: int) -> list[dict[str, Any]]:
    columns: list[dict[str, Any]] = []
    index = 0
    while index < len(lines):
        prefix = "@@SOURCE_SYZYGY_COLUMN="
        assert lines[index].startswith(prefix), lines[index]
        source_column = int(lines[index][len(prefix) :])
        entries = lines[index + 1 : index + 1 + rank]
        assert len(entries) == rank
        columns.append(
            {
                "source_syzygy_column": source_column,
                "entries": entries,
            }
        )
        index += rank + 1
    return columns


def singular_version(singular: str) -> str:
    completed = subprocess.run(
        [singular, "-v"],
        text=True,
        capture_output=True,
        check=True,
    )
    return completed.stdout.splitlines()[0]


def run_singular(
    universal_columns: list[list[sp.Expr]],
    central_columns: list[list[sp.Expr]],
    variables: tuple[sp.Symbol, ...],
    timeout: int,
) -> tuple[dict[str, Any], str]:
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    program = "\n".join(
        [
            (
                "ring filtered_cubic=0,"
                f"({','.join(map(str, variables))}),(dp(3),dp(24));"
            ),
            (
                "module CentralInput="
                + singular_module(central_columns, variables)
                + ";"
            ),
            (
                "module UniversalInput="
                + singular_module(universal_columns, variables)
                + ";"
            ),
            "module CentralBasis=std(CentralInput);",
            "resolution CentralResolution=mres(CentralBasis,0);",
            "module CentralSyzygies=syz(CentralInput);",
            (
                "matrix ObstructionMatrix="
                "matrix(UniversalInput)*matrix(CentralSyzygies);"
            ),
            "module Obstruction=module(ObstructionMatrix);",
            (
                "module ObstructionRemainder="
                "reduce(Obstruction,CentralBasis);"
            ),
            "poly BoundaryForm=x+y+z;",
            (
                "module CentralColon=std("
                "quotient(CentralBasis,ideal(BoundaryForm)));"
            ),
            'print("@@CENTRAL_STANDARD_BASIS_SIZE="+'
            "string(size(CentralBasis)));",
            'print("@@MINIMAL_D1_COLUMNS="+'
            "string(ncols(CentralResolution[1])));",
            'print("@@MINIMAL_D2_COLUMNS="+'
            "string(ncols(CentralResolution[2])));",
            'print("@@MINIMAL_D3_ZERO="+'
            "string(CentralResolution[3][1]==0));",
            'print("@@CENTRAL_SYZYGY_COUNT="+'
            "string(size(CentralSyzygies)));",
            'print("@@OBSTRUCTION_COLUMN_COUNT="+'
            "string(size(Obstruction)));",
            'print("@@NONZERO_REMAINDER_COUNT="+'
            "string(size(ObstructionRemainder)));",
            'print("@@OBSTRUCTION_TOTAL_DEGREE="+'
            "string(deg(Obstruction)));",
            'print("@@REMAINDER_TOTAL_DEGREE="+'
            "string(deg(ObstructionRemainder)));",
            (
                'print("@@CENTRAL_BOUNDARY_FORM_REGULAR="+string('
                "size(simplify(reduce("
                "CentralColon,std(CentralBasis)),2))==0"
                " and size(simplify(reduce("
                "CentralBasis,std(CentralColon)),2))==0));"
            ),
            'print("@@BEGIN_OBSTRUCTION_REMAINDER");',
            "int remainder_column;",
            "int remainder_row;",
            (
                "for (remainder_column=1;"
                "remainder_column<=ncols(ObstructionRemainder);"
                "remainder_column++)"
            ),
            "{",
            "  if (ObstructionRemainder[remainder_column]<>0)",
            "  {",
            (
                '    print("@@SOURCE_SYZYGY_COLUMN="+'
                "string(remainder_column));"
            ),
            (
                "    for (remainder_row=1;"
                "remainder_row<=nrows(ObstructionRemainder);"
                "remainder_row++)"
            ),
            "    {",
            (
                "      print(string(ObstructionRemainder["
                "remainder_row,remainder_column]));"
            ),
            "    }",
            "  }",
            "}",
            'print("@@END_OBSTRUCTION_REMAINDER");',
            'print("@@COMPLETE=1");',
            "quit;",
        ]
    )
    completed = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=True,
    )
    assert "@@COMPLETE=1" in completed.stdout
    assert "? error occurred" not in completed.stdout
    remainder = remainder_columns(
        section(completed.stdout, "OBSTRUCTION_REMAINDER"),
        rank=6,
    )
    remainder_serialization = json.dumps(
        remainder,
        separators=(",", ":"),
        sort_keys=True,
    )
    assert len(remainder) == marker(
        completed.stdout, "NONZERO_REMAINDER_COUNT"
    )
    data = {
        "central_standard_basis_size": marker(
            completed.stdout, "CENTRAL_STANDARD_BASIS_SIZE"
        ),
        "central_minimal_cokernel_resolution_ranks": [
            6,
            marker(completed.stdout, "MINIMAL_D1_COLUMNS"),
            marker(completed.stdout, "MINIMAL_D2_COLUMNS"),
        ],
        "central_resolution_stops_after_d2": bool(
            marker(completed.stdout, "MINIMAL_D3_ZERO")
        ),
        "central_syzygy_count_for_pruned_input": marker(
            completed.stdout, "CENTRAL_SYZYGY_COUNT"
        ),
        "unmodified_syzygy_obstruction_column_count": marker(
            completed.stdout, "OBSTRUCTION_COLUMN_COUNT"
        ),
        "nonzero_remainder_count_modulo_central_image": marker(
            completed.stdout, "NONZERO_REMAINDER_COUNT"
        ),
        "obstruction_total_degree": marker(
            completed.stdout, "OBSTRUCTION_TOTAL_DEGREE"
        ),
        "remainder_total_degree": marker(
            completed.stdout, "REMAINDER_TOTAL_DEGREE"
        ),
        "central_boundary_form": "x+y+z",
        "central_boundary_form_regular": bool(
            marker(completed.stdout, "CENTRAL_BOUNDARY_FORM_REGULAR")
        ),
        "obstruction_remainder_columns": remainder,
        "obstruction_remainder_sha256": hashlib.sha256(
            remainder_serialization.encode()
        ).hexdigest(),
    }
    return data, singular_version(singular)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout", type=int, default=180)
    args = parser.parse_args()

    parameters, tensor = universal_tensor()
    universal_columns, universal_pivots = (
        unit_pruned_differential_relations(
            cubic.differential_relations(
                cubic.CUBIC_STRATA["smooth"], tensor
            )
        )
    )
    central_columns, central_pivots = unit_pruned_differential_relations(
        cubic.differential_relations(cubic.CUBIC_STRATA["smooth"])
    )
    assert universal_pivots == central_pivots
    assert len(universal_columns) == len(central_columns) == 25
    assert len(universal_columns[0]) == len(central_columns[0]) == 6

    # Put collision variables in the first Singular block while retaining
    # the parameter variables as genuine polynomial variables.
    variables = tuple(cubic.BASE_VARIABLES) + tuple(parameters)
    filtration = entrywise_filtration_audit(
        universal_columns,
        central_columns,
        variables,
    )
    assert filtration[
        "strict_raise_on_entries_with_nonzero_central_part"
    ]

    singular_data, singular_release = run_singular(
        universal_columns,
        central_columns,
        variables,
        args.timeout,
    )
    assert singular_data[
        "central_minimal_cokernel_resolution_ranks"
    ] == [6, 13, 7]
    assert singular_data["central_resolution_stops_after_d2"]
    assert singular_data["central_boundary_form_regular"]
    assert singular_data[
        "nonzero_remainder_count_modulo_central_image"
    ] == 12

    artifact = {
        "schema": "universal-cubic-filtered-syzygy-frontier.v1",
        "mathematical_status": "exact finite syzygy obstruction",
        "universal_cotangent_saturation": "open",
        "ring": "Q[u1,...,u24,x,y,z]",
        "boundary_ideal": ["x", "y", "z"],
        "presentation": {
            "central_rows": 6,
            "central_columns": 25,
            "universal_rows": 6,
            "universal_columns": 25,
            "unit_pivots": [
                {"row": row, "column": column, "value": value}
                for row, column, value in universal_pivots
            ],
        },
        "entrywise_filtration": filtration,
        "central_resolution_and_syzygy_test": singular_data,
        "proved": [
            (
                "the central smooth cotangent module has a free cokernel "
                "resolution with ranks 7 -> 13 -> 6 and projective "
                "dimension at most two"
            ),
            (
                "x+y+z is regular on the central cotangent module"
            ),
            (
                "although every universal perturbation entry with a "
                "nonzero central part raises collision degree strictly, "
                "the unchanged central syzygies have twelve nonzero exact "
                "remainders modulo the central image"
            ),
            (
                "the submodule generated by those remainder classes "
                "inherits regularity of x+y+z from the central cotangent "
                "module, so it has no (x,y,z)-supported torsion"
            ),
        ],
        "interpretation": (
            "entrywise low-jet agreement does not by itself provide Rees "
            "strictness or a coefficient-independent lift of the central "
            "resolution; the recorded remainder module is a horizontal "
            "presentation-gauge mismatch, so the next finite problem is "
            "to construct corrected syzygies rather than saturate these "
            "classes away"
        ),
        "not_proved": [
            "Rees strictness of the universal cotangent presentation",
            "vanishing of the syzygy defects after allowed corrections",
            "regularity of a boundary form on the universal module",
            "universal cotangent saturation",
        ],
        "software": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
            "singular": singular_release,
        },
        "reproduce": (
            ".venv/bin/python "
            "scripts/verify_universal_cubic_filtered_syzygy_frontier.py"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        "PASS: central cotangent resolution ranks are 7 -> 13 -> 6"
    )
    print("PASS: central boundary form x+y+z is regular")
    print(
        "PASS: twelve exact unmodified-syzygy remainders are nonzero"
    )
    print(f"PASS: wrote {args.output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
