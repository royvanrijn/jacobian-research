#!/usr/bin/env python3
"""Exact quartic completion obstruction for five K12 graph families."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import sympy as sp

from audit_k12_coordinate_pair_frontier import build_k12
from audit_k12_parameterized_completion import quadratic_graph_family
from audit_k12_z8_cubic_completion import (
    PRIME,
    encode,
    fixed_minor_certificate,
    high_degree_group,
    multiply,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "k12_single_defect_quartic_completion_frontier.json"
)


def audit_family(
    k12: list[sp.Expr],
    z: tuple[sp.Symbol, ...],
    pivot: int,
    bad_component: int,
) -> dict[str, object]:
    parameters, coefficients, source_variables, restricted = (
        quadratic_graph_family(k12, z, pivot)
    )
    encoded = [
        encode(component, source_variables + parameters)
        for component in restricted
    ]
    available = [
        component
        for index, component in enumerate(encoded)
        if index != bad_component
    ]
    cache = {
        (index,): component
        for index, component in enumerate(available)
    }
    columns = []
    full_basis_count = 0
    for degree in (1, 2, 3, 4):
        for indices in itertools.combinations_with_replacement(
            range(10), degree
        ):
            full_basis_count += 1
            if degree == 1:
                product = cache[indices]
            else:
                product = multiply(
                    cache[indices[:-1]], available[indices[-1]]
                )
                if degree < 4:
                    cache[indices] = product
            column = high_degree_group(product)
            if column:
                columns.append(column)
    target = high_degree_group(encoded[bad_component])
    certificate, determinant = fixed_minor_certificate(
        columns,
        target,
        parameters,
        (sp.Integer(0),) * len(parameters),
    )
    assert full_basis_count == 1_000
    assert len(columns) == 990
    assert not determinant.free_symbols
    all_rows = set(target)
    for column in columns:
        all_rows.update(column)
    return {
        "source_pivot": pivot + 1,
        "parameter_count": len(parameters),
        "linear_target_coefficients": [
            str(value) for value in coefficients
        ],
        "selected_obstructed_component": (
            bad_component + 1
            if bad_component < pivot
            else bad_component + 2
        ),
        "maximum_target_degree": 4,
        "full_target_basis_count": full_basis_count,
        "nonzero_high_degree_columns": len(columns),
        "high_degree_row_count": len(all_rows),
        "certificate": certificate,
        "constant_column_determinant": str(determinant),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    z, k12 = build_k12()
    # (zero-based source pivot, zero-based retained bad component).
    specifications = [(6, 2), (8, 1), (9, 2), (10, 0), (11, 0)]
    records = [
        audit_family(k12, z, pivot, bad)
        for pivot, bad in specifications
    ]
    assert [record["source_pivot"] for record in records] == [
        7,
        9,
        10,
        11,
        12,
    ]
    artifact = {
        "format": "k12-single-defect-quartic-completion-frontier-v1",
        "status": "exact bounded obstruction over Q",
        "selection_prime": PRIME,
        "families": records,
        "proof_logic": (
            "For each parameterized quadratic graph family, sparse modular "
            "elimination at the literal point selects 990 column rows and "
            "one additional augmented row. Reconstruction over the full "
            "rational parameter ring gives a nonzero constant 990-by-990 "
            "column determinant and a nonzero constant multiple as the "
            "991-by-991 augmented determinant. Therefore the target defect "
            "is outside the entire degree-at-most-four completion span at "
            "every parameter value."
        ),
        "conclusion": (
            "None of the five single-defect quadratic graph-coordinate "
            "families z7,z9,z10,z11,z12 admits a one-stage target "
            "completion of degree at most four restoring degree at most "
            "three."
        ),
        "scope": (
            "The multi-defect z8 family at target degree four is not "
            "covered. Nor are cubic graph corrections, nonlinear target "
            "coordinates, or ordered multi-stage automorphisms."
        ),
    }
    serialized = json.dumps(artifact, indent=2) + "\n"
    if args.write:
        OUTPUT.write_text(serialized)
    else:
        assert OUTPUT.exists(), f"missing {OUTPUT.relative_to(ROOT)}"
        assert OUTPUT.read_text() == serialized, (
            f"{OUTPUT.relative_to(ROOT)} is stale; regenerate with --write"
        )
    digest = hashlib.sha256(serialized.encode()).hexdigest()
    for record in records:
        print(
            "PASS quartic completion: "
            f"pivot z{record['source_pivot']}, 990 constant-minor columns"
        )
    print(f"PASS checked {OUTPUT.relative_to(ROOT)}")
    print(f"SHA256 {digest}")


if __name__ == "__main__":
    main()
