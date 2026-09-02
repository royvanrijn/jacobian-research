#!/usr/bin/env python3
"""Extract one rank-seven catalogue surface for direct low-MW source search.

The sequential prescribed-root engine needs one ordered auxiliary Gram and one
same-surface target-frame Gram.  This adapter selects both from the exact
surface-first catalogue without changing their lattice data.  The
``support_design`` field is a compatibility alias used by the current ordered
embedding engine; it does not assert that a general auxiliary is a Golay
design.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CATALOGUE = (
    ROOT / "artifacts/generated-results/elkies-k3-rank7-auxiliary-catalogue-v1.json"
)
DEFAULT_T_ARITHMETIC = (
    ROOT / "artifacts/generated-results/elkies-k3-rank7-t-arithmetic-v1.json"
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_key(value: list[list[int]]) -> tuple:
    diagonal = tuple(row[index] for index, row in enumerate(value))
    return max(diagonal), sum(diagonal), tuple(tuple(row) for row in value)


def determinant_bareiss(value: list[list[int]]) -> int:
    work = [list(map(int, row)) for row in value]
    size = len(work)
    sign = 1
    denominator = 1
    for column in range(size - 1):
        pivot = next(
            (row for row in range(column, size) if work[row][column]), None
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        pivot_value = work[column][column]
        for row in range(column + 1, size):
            for index in range(column + 1, size):
                numerator = (
                    work[row][index] * pivot_value
                    - work[row][column] * work[column][index]
                )
                if numerator % denominator:
                    raise ArithmeticError("Bareiss division was not exact")
                work[row][index] = numerator // denominator
        denominator = pivot_value
    return sign * work[-1][-1]


def validate_gram(value: list[list[int]], rank: int, determinant: int) -> None:
    if len(value) != rank or any(len(row) != rank for row in value):
        raise ArithmeticError(f"selected Gram is not rank {rank}")
    if any(value[left][right] != value[right][left] for left in range(rank) for right in range(rank)):
        raise ArithmeticError("selected Gram is not symmetric")
    if any(value[index][index] % 2 for index in range(rank)):
        raise ArithmeticError("selected Gram is not even")
    if determinant_bareiss(value) != determinant:
        raise ArithmeticError("selected Gram determinant mismatch")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue", type=Path, default=DEFAULT_CATALOGUE)
    parser.add_argument("--t-arithmetic", type=Path, default=DEFAULT_T_ARITHMETIC)
    parser.add_argument("--surface-id", required=True)
    parser.add_argument("--partner-index", type=int)
    parser.add_argument("--frame-id")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    catalogue_path = arguments.catalogue.resolve()
    arithmetic_path = arguments.t_arithmetic.resolve()
    catalogue = json.loads(catalogue_path.read_text())
    if catalogue.get("schema") != "elkies-k3.rank7-auxiliary-catalogue.v1":
        raise ValueError("unexpected rank-seven catalogue schema")
    surface = next(
        (row for row in catalogue["surfaces"] if row["surface_id"] == arguments.surface_id),
        None,
    )
    if surface is None:
        raise ValueError(f"unknown surface id: {arguments.surface_id}")
    arithmetic_ledger = json.loads(arithmetic_path.read_text())
    if arithmetic_ledger.get("schema") != "elkies-k3.rank7-t-arithmetic.v1":
        raise ValueError("unexpected rank-seven T-arithmetic schema")
    if arithmetic_ledger["input"]["catalogue_sha256"] != digest(catalogue_path):
        raise ArithmeticError("T-arithmetic ledger does not match catalogue hash")
    arithmetic = next(
        (
            row
            for row in arithmetic_ledger["surfaces"]
            if row["surface_id"] == arguments.surface_id
        ),
        None,
    )
    if arithmetic is None:
        raise ArithmeticError("selected surface has no T-arithmetic row")
    if not arithmetic["pre_solver_gate"]["arithmetic_attempt_recorded"]:
        raise ArithmeticError("selected surface has not passed the T-arithmetic gate")
    if not arithmetic["pre_solver_gate"]["equation_solver_may_launch"]:
        raise ArithmeticError(
            "selected surface has typed-open T-arithmetic curve identification; "
            "equation-target extraction remains blocked"
        )

    partners = surface["partner_auxiliaries"]
    if arguments.partner_index is None:
        partner_index, partner = min(
            enumerate(partners, start=1),
            key=lambda item: matrix_key(item[1]["reduced_gram"]),
        )
    else:
        partner_index = arguments.partner_index
        if not 1 <= partner_index <= len(partners):
            raise ValueError("--partner-index is outside the surface partner list")
        partner = partners[partner_index - 1]

    frames = surface["frames"]
    if arguments.frame_id is None:
        frame = min(
            frames,
            key=lambda row: (-int(row["mw_rank_for_rho_19"]), row["frame_id"]),
        )
    else:
        frame = next((row for row in frames if row["frame_id"] == arguments.frame_id), None)
        if frame is None:
            raise ValueError("--frame-id is not attached to the selected surface")

    auxiliary = partner["reduced_gram"]
    if int(partner["determinant"]) != int(surface["determinant"]):
        raise ArithmeticError("surface and auxiliary determinants differ")
    if int(frame["determinant"]) != int(surface["determinant"]):
        raise ArithmeticError("surface and target-frame determinants differ")
    determinant = int(surface["determinant"])
    validate_gram(auxiliary, 7, determinant)
    validate_gram(frame["gram"], 17, determinant)

    output = {
        "schema": "elkies-k3.rank7-catalogue-source-search-target.v1",
        "status": "PASS_EXACT_CATALOGUE_SOURCE_SEARCH_TARGET_EXTRACTION",
        "input": {
            "catalogue": relative(catalogue_path),
            "catalogue_sha256": digest(catalogue_path),
            "t_arithmetic": relative(arithmetic_path),
            "t_arithmetic_sha256": digest(arithmetic_path),
        },
        "surface_id": surface["surface_id"],
        "legacy_ns_ids": surface["legacy_ns_ids"],
        "determinant": int(surface["determinant"]),
        "auxiliary": {
            "partner_index_one_based": partner_index,
            "gram_sha256": partner["gram_sha256"],
            "catalogue_gram": partner["gram"],
            "ordered_reduced_gram": auxiliary,
            "ordered_basis_diagonal": [row[index] for index, row in enumerate(auxiliary)],
            "provenance": partner["provenance"],
        },
        "frame": {
            "frame_id": frame["frame_id"],
            "gram": frame["gram"],
            "gram_sha256": frame["gram_sha256"],
            "root_type": frame["root_type"],
            "root_rank": int(frame["root_rank"]),
            "mw_rank_for_rho_19": int(frame["mw_rank_for_rho_19"]),
            "determinant": int(frame["determinant"]),
        },
        "support_design": {
            "raw_octad_intersection_gram": auxiliary,
            "compatibility_status": (
                "ORDERED_AUXILIARY_GRAM_ALIAS_FOR_CURRENT_SEQUENTIAL_ENGINE; "
                "NO_GOLAY_DESIGN_CLAIM"
            ),
        },
        "t_arithmetic_pre_solver_gate": arithmetic,
        "proof_boundary": {
            "proved": (
                "The ordered auxiliary and target frame are exact records on the "
                "same catalogue (T,NS) surface and have the same determinant."
                " A hash-matched T-arithmetic attempt is attached before this target "
                "can be emitted."
            ),
            "not_proved": (
                "This adapter performs no new embedding search, rational marking, "
                "equation construction, or neighbour certification. Typed open "
                "arithmetic fields in the attached row are not promoted."
            ),
        },
        "reproduce": (
            "python3 elkies-k3/scripts/extract_rank7_catalogue_source_search_target.py "
            f"--surface-id {surface['surface_id']} --partner-index {partner_index} "
            f"--frame-id {frame['frame_id']} --output {relative(arguments.output.resolve())}"
        ),
    }
    serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if output_path.read_text() != serialized:
            raise SystemExit("rank-seven source-search target adapter is stale")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "RANK7SOURCETARGET|"
        f"surface={surface['surface_id']}|partner={partner_index}|"
        f"target_mw={frame['mw_rank_for_rho_19']}|status=PASS"
    )


if __name__ == "__main__":
    main()
