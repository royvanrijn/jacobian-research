#!/usr/bin/env python3
"""Close the complete eight-support bidegree-(3,3) coefficient census.

The full through-12 census has fourteen nonunit coefficient-torus systems
and one timeout.  This checker does three exact things:

* validates the stored characteristic-zero unit certificate for the timed
  out odd-parity system, rerunnable with ``--rerun-parity``;
* requests full complex rational parametrizations for the fourteen nonunit
  systems and checks their scheme and squarefree-eliminant degrees; and
* verifies explicit one-sided normal forms giving as many distinct complex
  points as those eliminants permit.

Consequently every complex point on the fourteen nonunit systems is in the
pair-linear one-sided nullcone.  Together with the separate size-at-most-
seven certificate, an actual SIC counterexample in the standard monomial
basis must have at least nine nonzero coefficients.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

import sympy as sp

from research_two_pair_sic_bidegree33_sparse_six_counterexample import (
    screen_support,
)


ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "artifacts" / "generated-results"
CENSUS = (
    GENERATED
    / "two_pair_sic_bidegree33_sparse_support8_screen.json"
)
PARITY = (
    GENERATED
    / "two_pair_sic_bidegree33_sparse_support8_parity_msolve14_char0.json"
)
OUTPUT = (
    GENERATED
    / "two_pair_sic_bidegree33_sparse_support8_closure.json"
)
SPECIAL_CUBE_ROOT = (
    (0, 0),
    (0, 3),
    (1, 1),
    (1, 2),
    (2, 1),
    (2, 2),
    (3, 0),
    (3, 3),
)
SPECIAL_TWO_FACTOR = (
    (0, 1),
    (0, 2),
    (1, 0),
    (1, 3),
    (2, 0),
    (2, 3),
    (3, 1),
    (3, 2),
)
ODD_PARITY = (
    (0, 1),
    (0, 3),
    (1, 0),
    (1, 2),
    (2, 1),
    (2, 3),
    (3, 0),
    (3, 2),
)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rerun-parity",
        action="store_true",
        help=(
            "rerun the approximately eight-minute exact msolve unit "
            "calculation and refresh its pinned raw artifact"
        ),
    )
    return parser.parse_args()


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def support_key(support: object) -> tuple[tuple[int, int], ...]:
    return tuple(tuple(map(int, position)) for position in support)


def row_support(rows: tuple[int, int]) -> tuple[tuple[int, int], ...]:
    return tuple((i, j) for i in rows for j in range(4))


def column_support(columns: tuple[int, int]) -> tuple[tuple[int, int], ...]:
    return tuple((i, j) for i in range(4) for j in columns)


def validate_census() -> dict[str, object]:
    payload = json.loads(CENSUS.read_text(encoding="utf-8"))
    if not payload["complete_mixed_enumeration"]:
        raise AssertionError("the size-eight census is not complete")
    if payload["through"] != 12 or payload["support_size"] != 8:
        raise AssertionError("unexpected size-eight census scope")
    if payload["counts"] != {
        "excluded_on_coefficient_torus": 12765,
        "msolve_nonempty": 14,
        "timeout": 1,
    }:
        raise AssertionError("unexpected size-eight census counts")

    rectangles = {
        *(row_support(rows) for rows in combinations(range(4), 2)),
        *(column_support(columns) for columns in combinations(range(4), 2)),
    }
    expected_nonunit = rectangles | {
        SPECIAL_CUBE_ROOT,
        SPECIAL_TWO_FACTOR,
    }
    actual_nonunit = {
        support_key(record["support"])
        for record in payload["records"]
        if record["status"] == "msolve_nonempty"
    }
    actual_timeout = {
        support_key(record["support"])
        for record in payload["records"]
        if record["status"] == "timeout"
    }
    if actual_nonunit != expected_nonunit:
        raise AssertionError("unexpected nonunit supports in size-eight census")
    if actual_timeout != {ODD_PARITY}:
        raise AssertionError("unexpected timeout support in size-eight census")
    return {
        "mixed_supports": payload["mixed_support_count"],
        "unit_systems_through_12": 12765,
        "nonunit_systems": 14,
        "timed_out_systems": 1,
        "sha256": digest(CENSUS),
    }


def rerun_parity() -> None:
    with tempfile.TemporaryDirectory(
        prefix="sic33-support8-parity-"
    ) as directory:
        path = Path(directory) / "parity.json"
        subprocess.run(
            [
                sys.executable,
                "scripts/research_two_pair_sic_bidegree33_sparse_six_counterexample.py",
                "--support-size",
                "8",
                "--through",
                "14",
                "--start",
                "8384",
                "--limit",
                "1",
                "--timeout",
                "600",
                "--threads",
                "4",
                "--output",
                str(path),
            ],
            cwd=ROOT,
            check=True,
        )
        payload = json.loads(path.read_text(encoding="utf-8"))
    payload["records"][0]["seconds"] = round(
        float(payload["records"][0]["seconds"]),
        6,
    )
    PARITY.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def validate_parity() -> dict[str, object]:
    payload = json.loads(PARITY.read_text(encoding="utf-8"))
    if (
        payload["support_size"] != 8
        or payload["through"] != 14
        or payload["start"] != 8384
        or payload["stop"] != 8385
    ):
        raise AssertionError("unexpected parity artifact scope")
    record = payload["records"][0]
    if support_key(record["support"]) != ODD_PARITY:
        raise AssertionError("unexpected parity support")
    if (
        record["status"] != "excluded_on_coefficient_torus"
        or record["msolve"]["returncode"] != 0
        or record["msolve"]["result"] not in ("[-1]", "[-1]:")
    ):
        raise AssertionError("the parity system is not an exact unit")
    return {
        "mixed_support_index": 8384,
        "support": [list(position) for position in ODD_PARITY],
        "moments": [1, 14],
        "result": "unit ideal over QQ",
        "sha256": digest(PARITY),
    }


def rur_degrees(
    support: tuple[tuple[int, int], ...],
    msolve: str,
) -> tuple[int, int, str]:
    record = screen_support(
        support,
        through=12,
        timeout=60,
        msolve=msolve,
        threads=1,
        rational_parametrization=True,
    )
    if record["status"] != "msolve_nonempty":
        raise AssertionError(f"unexpected RUR status on {support}")
    text = record["msolve"]["result_head"]
    header = re.match(
        r"\[0,\s*\[0,\s*(\d+),\s*(\d+),",
        text,
    )
    eliminant = re.search(
        r"\[1,\s*\[\s*\[\s*(\d+),",
        text,
    )
    if header is None or eliminant is None:
        raise AssertionError(f"could not parse RUR on {support}")
    variable_count, scheme_degree = map(int, header.groups())
    if variable_count != 8:
        raise AssertionError(f"unexpected RUR variable count on {support}")
    return scheme_degree, int(eliminant.group(1)), hashlib.sha256(
        " ".join(text.split()).encode("utf-8")
    ).hexdigest()


def coefficient(
    polynomial: sp.Poly,
    i: int,
    j: int,
    W: sp.Symbol,
    V: sp.Symbol,
    Z: sp.Symbol,
    Y: sp.Symbol,
) -> sp.Expr:
    return polynomial.coeff_monomial(
        W ** (3 - i) * V**i * Z ** (3 - j) * Y**j
    )


def verify_explicit_nullcone_forms() -> dict[str, object]:
    W, V, Z, Y, q, r = sp.symbols("W V Z Y q r")

    row_records: list[dict[str, object]] = []
    for a, b in combinations(range(4), 2):
        gap = b - a
        form = (
            W ** (3 - a) * V**a
            - 3**gap * W ** (3 - b) * V**b
        ) * (Z + Y / 3) ** 3
        polynomial = sp.Poly(sp.expand(form), W, V, Z, Y)
        support = {
            (i, j)
            for i in range(4)
            for j in range(4)
            if coefficient(polynomial, i, j, W, V, Z, Y) != 0
        }
        if support != set(row_support((a, b))):
            raise AssertionError("row normal form has unexpected support")
        if coefficient(polynomial, a, 0, W, V, Z, Y) != 1:
            raise AssertionError("row normal form misses first normalization")
        if coefficient(polynomial, a, 1, W, V, Z, Y) != 1:
            raise AssertionError("row normal form misses second normalization")
        row_records.append(
            {
                "rows": [a, b],
                "gap": gap,
                "factorization": (
                    f"(W^{3-a}*V^{a}-3^{gap}*W^{3-b}*V^{b})"
                    "*(Z+Y/3)^3"
                ),
                "distinct_complex_points": 1,
            }
        )

    column_records: list[dict[str, object]] = []
    for c, d in combinations(range(4), 2):
        gap = d - c
        relation = r**gap - (-1) ** (gap + 1)
        form = (W - r * V) ** 3 * (
            Z ** (3 - c) * Y**c + Z ** (3 - d) * Y**d
        )
        polynomial = sp.Poly(sp.expand(form), W, V, Z, Y)
        for i in range(4):
            if coefficient(polynomial, i, c, W, V, Z, Y) != coefficient(
                polynomial, i, d, W, V, Z, Y
            ):
                raise AssertionError("column normal form columns differ")
        if coefficient(polynomial, 0, c, W, V, Z, Y) != 1:
            raise AssertionError("column normal form misses normalization")
        if sp.Poly(relation, r).degree() != gap:
            raise AssertionError("unexpected column root count")
        column_records.append(
            {
                "columns": [c, d],
                "gap": gap,
                "factorization": (
                    f"(W-r*V)^3*(Z^{3-c}*Y^{c}+"
                    f"Z^{3-d}*Y^{d})"
                ),
                "parameter_equation": (
                    f"r^{gap}={(-1) ** (gap + 1)}"
                ),
                "distinct_complex_points": gap,
            }
        )

    positive_positions = tuple(
        (i, j) for i in range(4) for j in range(4) if i > j
    )
    parameters = (-3 * q**2, -3 * q, 0, -1, 3 * q, -3 * q**2)
    cube_root_form = sp.expand(
        sum(
            parameter
            * W ** (3 - i)
            * (V - q * W) ** i
            * (Z + q * Y) ** (3 - j)
            * Y**j
            for parameter, (i, j) in zip(
                parameters,
                positive_positions,
                strict=True,
            )
        )
    )
    cube_polynomial = sp.Poly(cube_root_form, W, V, Z, Y)
    cube_matrix = [
        [
            sp.rem(
                sp.Poly(
                    coefficient(
                        cube_polynomial,
                        i,
                        j,
                        W,
                        V,
                        Z,
                        Y,
                    ),
                    q,
                ),
                sp.Poly(q**3 - 1, q),
            ).as_expr()
            for j in range(4)
        ]
        for i in range(4)
    ]
    expected_cube_matrix = [
        [1, 0, 0, 1],
        [0, 9, 9 * q, 0],
        [0, -9 * q**2, -9, 0],
        [-1, 0, 0, -1],
    ]
    if cube_matrix != expected_cube_matrix:
        raise AssertionError("cube-root one-sided normal form mismatch")

    A, B, P, Q = W - V, W + V, Z + Y, Z - Y
    two_factor_form = sp.expand(
        -A * P * (A * Q - B * P) * (A * Q + B * P) / 4
    )
    two_factor_polynomial = sp.Poly(two_factor_form, W, V, Z, Y)
    expected_two_factor = sp.zeros(4)
    for value, position in zip(
        (1, 1, 1, 1, -1, -1, -1, -1),
        SPECIAL_TWO_FACTOR,
        strict=True,
    ):
        expected_two_factor[position] = value
    actual_two_factor = sp.Matrix(
        4,
        4,
        lambda i, j: coefficient(
            two_factor_polynomial,
            i,
            j,
            W,
            V,
            Z,
            Y,
        ),
    )
    if actual_two_factor != expected_two_factor:
        raise AssertionError("two-factor one-sided factorization mismatch")
    if sp.expand(
        (W - V) * (Z - Y) / 2
        + (W + V) * (Z + Y) / 2
        - (W * Z + V * Y)
    ) != 0:
        raise AssertionError("the displayed pair change is not symplectic")

    return {
        "row_rectangles": row_records,
        "column_rectangles": column_records,
        "cube_root_exception": {
            "support": [list(position) for position in SPECIAL_CUBE_ROOT],
            "parameter_equation": "q^3=1",
            "distinct_complex_points": 3,
            "one_sided_parameters_i_gt_j": [
                "-3*q^2",
                "-3*q",
                "0",
                "-1",
                "3*q",
                "-3*q^2",
            ],
        },
        "two_factor_exception": {
            "support": [list(position) for position in SPECIAL_TWO_FACTOR],
            "distinct_complex_points": 1,
            "factorization": "-A*P*(A*Q-B*P)*(A*Q+B*P)/4",
            "pair_change": (
                "A=W-V, B=W+V, P=Z+Y, Q=Z-Y; "
                "W'=A,V'=B,Z'=Q/2,Y'=P/2"
            ),
            "transformed_support": ["(0,1)", "(2,3)"],
        },
    }


def verify_rurs(msolve: str) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for rows in combinations(range(4), 2):
        support = row_support(rows)
        scheme, eliminant, checksum = rur_degrees(support, msolve)
        if (scheme, eliminant) != (5, 1):
            raise AssertionError(f"unexpected row RUR degrees on {rows}")
        records.append(
            {
                "type": "row rectangle",
                "indices": list(rows),
                "scheme_degree": scheme,
                "squarefree_eliminant_degree": eliminant,
                "normalized_rur_sha256": checksum,
            }
        )
    for columns in combinations(range(4), 2):
        support = column_support(columns)
        gap = columns[1] - columns[0]
        scheme, eliminant, checksum = rur_degrees(support, msolve)
        if (scheme, eliminant) != (5 * gap, gap):
            raise AssertionError(
                f"unexpected column RUR degrees on {columns}"
            )
        records.append(
            {
                "type": "column rectangle",
                "indices": list(columns),
                "scheme_degree": scheme,
                "squarefree_eliminant_degree": eliminant,
                "normalized_rur_sha256": checksum,
            }
        )
    for label, support, expected in (
        ("cube-root exception", SPECIAL_CUBE_ROOT, (24, 3)),
        ("two-factor exception", SPECIAL_TWO_FACTOR, (28, 1)),
    ):
        scheme, eliminant, checksum = rur_degrees(support, msolve)
        if (scheme, eliminant) != expected:
            raise AssertionError(f"unexpected RUR degrees on {label}")
        records.append(
            {
                "type": label,
                "scheme_degree": scheme,
                "squarefree_eliminant_degree": eliminant,
                "normalized_rur_sha256": checksum,
            }
        )
    return records


def main() -> None:
    options = arguments()
    if options.rerun_parity:
        rerun_parity()
    msolve = shutil.which("msolve")
    if msolve is None:
        raise SystemExit("msolve is required")

    census = validate_census()
    parity = validate_parity()
    normal_forms = verify_explicit_nullcone_forms()
    rurs = verify_rurs(msolve)
    total_points = sum(
        int(record["squarefree_eliminant_degree"]) for record in rurs
    )
    if total_points != 20:
        raise AssertionError("unexpected total number of normalized points")

    payload = {
        "calculation": "two_pair_sic_bidegree33_sparse_support8_closure",
        "status": "proved",
        "field": "characteristic zero",
        "scope": (
            "complete standard-monomial coefficient-torus classification "
            "for supports of size eight"
        ),
        "census": census,
        "parity_unit": parity,
        "rational_univariate_representations": rurs,
        "explicit_one_sided_normal_forms": normal_forms,
        "distinct_normalized_complex_points_on_nonunit_systems": 20,
        "conclusion": (
            "the sole timeout is a unit ideal through order 14, and every "
            "complex point on each of the fourteen nonunit through-12 "
            "systems is pair-linearly one-sided; with the separate "
            "size-at-most-seven result, every actual SIC counterexample "
            "has standard monomial support at least nine"
        ),
    }
    OUTPUT.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("PASS complete size-eight coefficient-torus census")
    print("PASS odd-parity timeout resolved by a QQ unit ideal")
    print("PASS all 20 normalized complex points are one-sided")
    print("PASS any standard-basis SIC counterexample has support at least 9")


if __name__ == "__main__":
    main()
