#!/usr/bin/env sage-python
"""Compute Hhat^(-1)(C2,L) from a full integral involution lattice.

The input is deliberately strict: it must describe the *full* free lattice,
not visible invariant and anti-invariant sublattices.  With columns as lattice
coordinates, ``involution_matrix_rows`` is the matrix S acting by ``v |-> S*v``.

Input schema::

    {
      "schema": "elkies-k3.full-involution-lattice.v1",
      "basis_labels": ["..."],
      "involution_matrix_rows": [[...]],
      "gram_matrix_rows": [[...]],                 # optional
      "candidate_vectors_in_basis": [[...]]        # optional
    }

For A-=ker_Z(1+S), choose a saturated column basis B-.  The script computes
the integral matrix D defined by

    1-S = B- D.

Then Hhat^(-1)=A-/(1-S)A=coker(D).  Its Smith invariants must all be 1 or 2.
The returned quotient-functional rows give canonical-enough F2 signatures:
a vector c in anti-basis coordinates maps to (ell*c)_ell.

This script cannot compute the seventeen product quotients from the currently
stored visible R17 lattice; it becomes applicable only after a full base-change
MW lattice (or an equivalent integral presentation) is supplied.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from sage.all import GF, QQ, ZZ, identity_matrix, matrix, vector


SCHEMA = "elkies-k3.full-involution-lattice.v1"
OUTPUT_SCHEMA = "elkies-k3.involution-tate-hminus1.v1"


def integral_matrix(value, name):
    try:
        result = matrix(ZZ, value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} is not an integral matrix") from error
    return result


def integral_solution(left, right, name):
    solution = left.change_ring(QQ).solve_right(right.change_ring(QQ))
    if any(entry.denominator() != 1 for entry in solution.list()):
        raise ArithmeticError(f"{name} is not integral")
    return matrix(ZZ, solution.nrows(), solution.ncols(), [ZZ(x) for x in solution.list()])


def compute(record):
    if record.get("schema") != SCHEMA:
        raise ValueError(f"expected input schema {SCHEMA}")
    labels = record.get("basis_labels")
    if not isinstance(labels, list) or not labels:
        raise ValueError("basis_labels must be a nonempty list")
    rank = len(labels)
    if len(set(map(str, labels))) != rank:
        raise ValueError("basis_labels must be distinct")

    involution = integral_matrix(record.get("involution_matrix_rows"), "involution")
    if involution.dimensions() != (rank, rank):
        raise ValueError("involution matrix has the wrong dimensions")
    identity = identity_matrix(ZZ, rank)
    if involution * involution != identity:
        raise ArithmeticError("the supplied matrix is not an involution")

    gram = None
    if "gram_matrix_rows" in record:
        gram = integral_matrix(record["gram_matrix_rows"], "Gram")
        if gram.dimensions() != (rank, rank) or gram != gram.transpose():
            raise ValueError("Gram matrix must be symmetric of full rank dimensions")
        if gram.det() == 0:
            raise ArithmeticError("Gram matrix is singular")
        if involution.transpose() * gram * involution != gram:
            raise ArithmeticError("involution does not preserve the Gram matrix")

    norm = identity + involution
    difference = identity - involution
    anti_basis = norm.right_kernel().basis_matrix().transpose()
    anti_rank = anti_basis.ncols()
    if anti_rank == 0:
        difference_coordinates = matrix(ZZ, 0, rank)
    else:
        difference_coordinates = integral_solution(
            anti_basis, difference, "coordinates of (1-S)A in A-"
        )
        if anti_basis * difference_coordinates != difference:
            raise ArithmeticError("anti-basis coordinate reconstruction failed")

    if anti_rank:
        smith_diagonal, _left_change, _right_change = difference_coordinates.smith_form()
        smith_invariants = [abs(int(smith_diagonal[index, index])) for index in range(anti_rank)]
    else:
        smith_invariants = []
    if any(value not in (1, 2) for value in smith_invariants):
        raise ArithmeticError(
            "Hhat^(-1) must be killed by two, but the Smith form has another invariant"
        )

    reduced_difference = matrix(GF(2), difference_coordinates)
    glue_dimension = int(reduced_difference.rank())
    cohomology_dimension = anti_rank - glue_dimension
    functionals = reduced_difference.transpose().right_kernel().basis_matrix()
    if functionals.nrows() != cohomology_dimension:
        raise ArithmeticError("quotient-functional dimension mismatch")

    candidates = []
    for index, raw in enumerate(record.get("candidate_vectors_in_basis", [])):
        point = vector(ZZ, raw)
        if len(point) != rank:
            raise ValueError(f"candidate {index} has the wrong rank")
        if norm * point:
            raise ArithmeticError(f"candidate {index} is not anti-invariant")
        if anti_rank:
            coordinates_matrix = integral_solution(
                anti_basis,
                matrix(ZZ, rank, 1, list(point)),
                f"candidate {index} anti-basis coordinates",
            )
            coordinates = vector(ZZ, coordinates_matrix.column(0))
        else:
            coordinates = vector(ZZ, [])
        signature = [
            int(sum(functionals[row, column] * GF(2)(coordinates[column])
                    for column in range(anti_rank)))
            for row in range(cohomology_dimension)
        ]
        item = {
            "index": index,
            "anti_basis_coordinates": [int(value) for value in coordinates],
            "hminus1_signature": signature,
            "is_coboundary": not any(signature),
        }
        if gram is not None:
            item["height"] = int(point * gram * point)
        candidates.append(item)

    return {
        "schema": OUTPUT_SCHEMA,
        "status": "PASS_EXACT_TATE_HMINUS1_FROM_FULL_INTEGRAL_INVOLUTION",
        "rank": rank,
        "anti_rank": anti_rank,
        "character_glue_dimension": glue_dimension,
        "hminus1_dimension": cohomology_dimension,
        "hminus1_smith_invariants": [
            value for value in smith_invariants if value != 1
        ],
        "anti_basis_columns": [
            [int(anti_basis[row, column]) for row in range(rank)]
            for column in range(anti_rank)
        ],
        "difference_coordinates_rows": [
            [int(value) for value in row] for row in difference_coordinates.rows()
        ],
        "quotient_functional_rows_mod2": [
            [int(value) for value in row] for row in functionals.rows()
        ],
        "candidates": candidates,
        "claim_boundary": (
            "Exact only for the supplied full free lattice. A visible sublattice gives "
            "neither the full character-glue group nor Hhat^(-1)."
        ),
    }


def self_test():
    split = compute(
        {
            "schema": SCHEMA,
            "basis_labels": ["plus", "minus"],
            "involution_matrix_rows": [[1, 0], [0, -1]],
            "gram_matrix_rows": [[2, 0], [0, 2]],
            "candidate_vectors_in_basis": [[0, 1]],
        }
    )
    if (
        split["anti_rank"] != 1
        or split["character_glue_dimension"] != 0
        or split["hminus1_dimension"] != 1
        or split["candidates"][0]["is_coboundary"]
    ):
        raise ArithmeticError("split-eigensum control failed")

    swapped = compute(
        {
            "schema": SCHEMA,
            "basis_labels": ["left", "right"],
            "involution_matrix_rows": [[0, 1], [1, 0]],
            "gram_matrix_rows": [[2, 0], [0, 2]],
            "candidate_vectors_in_basis": [[1, -1]],
        }
    )
    if (
        swapped["anti_rank"] != 1
        or swapped["character_glue_dimension"] != 1
        or swapped["hminus1_dimension"] != 0
        or not swapped["candidates"][0]["is_coboundary"]
    ):
        raise ArithmeticError("swapped-pair graph-glue control failed")
    return {
        "status": "PASS_SELF_TESTS",
        "controls": {
            "split_eigensum_hminus1_dimension": split["hminus1_dimension"],
            "swapped_pair_hminus1_dimension": swapped["hminus1_dimension"],
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    arguments = parser.parse_args()
    if arguments.self_test:
        print(json.dumps(self_test(), indent=2, sort_keys=True))
        return
    if arguments.input is None:
        parser.error("--input is required unless --self-test is used")
    record = json.loads(arguments.input.read_text())
    result = compute(record)
    result["input"] = {
        "path": str(arguments.input),
        "sha256": sha256(arguments.input.read_bytes()).hexdigest(),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if arguments.output is None or not arguments.output.exists():
            raise ArithmeticError("--check requires an existing --output")
        if arguments.output.read_text() != serialized:
            raise ArithmeticError("stored output differs from exact replay")
    elif arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(serialized)
    else:
        print(serialized, end="")


if __name__ == "__main__":
    main()
