#!/usr/bin/env python3
"""Combine a global dimension bound with explicit norm/local obstruction witnesses.

This is the non-enumerative companion to the quotient-native matrix builder.
Suppose a certified argument gives ``dim V <= D`` for a global squareclass
space containing the 2-Selmer group and the known MW image.  Exact global
elements need not span V: if their images in the direct sum of norm/local
condition cokernels have rank r, then the full condition map has rank at least
r and

    dim Sel_2(E) / image(MW_known) <= D - r - dim(MW_known).

Thus partial relation collection can monotonically improve a rigorous upper
bound.  No anonymous class-group direction is silently discarded.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys
from typing import Any, Mapping


CAS = Path(__file__).resolve().parent
sys.path.insert(0, str(CAS))
from residual_selmer_quotient import F2Error, f2_rank_rows  # noqa: E402


INPUT_SCHEMA = "elliptic-curves.mw29-relative-selmer-witness-input.v1"
OUTPUT_SCHEMA = "elliptic-curves.mw29-relative-selmer-witness-bound.v1"
PROTOCOL = "MW29REL2WIT-v1"


def _binary_row(value: object, width: int, name: str) -> list[int]:
    if not isinstance(value, list) or len(value) != width:
        raise F2Error(f"{name} must be a binary row of width {width}")
    row = [int(bit) for bit in value]
    if any(bit not in (0, 1) for bit in row):
        raise F2Error(f"{name} contains a non-binary entry")
    return row


def _rank_with_blocks(
    rows: list[list[int]], blocks: list[dict[str, Any]], chosen: tuple[int, ...]
) -> int:
    coordinates = [
        coordinate
        for block_index in chosen
        for coordinate in blocks[block_index]["coordinate_indices"]
    ]
    return f2_rank_rows(
        ([row[coordinate] for coordinate in coordinates] for row in rows),
        len(coordinates),
    )


def audit(
    manifest: Mapping[str, object], *, maximum_cut_size: int = 8,
    maximum_cut_subsets: int = 1_000_000,
) -> Mapping[str, object]:
    if manifest.get("schema") != INPUT_SCHEMA:
        raise F2Error(f"input schema must be {INPUT_SCHEMA}")
    global_upper = int(manifest.get("global_ambient_dimension_upper_bound", -1))
    known_dimension = int(manifest.get("known_mw_dimension", 29))
    if global_upper < 0 or known_dimension < 0 or global_upper < known_dimension:
        raise F2Error("global and known dimensions are inconsistent")

    block_input = manifest.get("condition_blocks")
    if not isinstance(block_input, list):
        raise F2Error("condition_blocks must be a list")
    blocks: list[dict[str, Any]] = []
    cursor = 0
    for index, block in enumerate(block_input):
        if not isinstance(block, Mapping):
            raise F2Error(f"condition block {index} is not an object")
        width = int(block.get("width", -1))
        if width < 0:
            raise F2Error(f"condition block {index} has negative width")
        blocks.append(
            {
                "place": str(block.get("place", f"block-{index}")),
                "width": width,
                "coordinate_indices": list(range(cursor, cursor + width)),
            }
        )
        cursor += width
    if len({block["place"] for block in blocks}) != len(blocks):
        raise F2Error("condition block labels are not unique")

    witness_input = manifest.get("witnesses")
    if not isinstance(witness_input, list):
        raise F2Error("witnesses must be a list")
    labels = []
    rows = []
    generators = []
    for index, witness in enumerate(witness_input):
        if not isinstance(witness, Mapping):
            raise F2Error(f"witness {index} is not an object")
        labels.append(str(witness.get("label", f"witness-{index}")))
        generators.append(witness.get("generator"))
        rows.append(_binary_row(witness.get("condition_syndrome"), cursor, f"witness {index}"))
    if len(set(labels)) != len(labels):
        raise F2Error("witness labels are not unique")

    full_rank = f2_rank_rows(rows, cursor)
    raw_residual_upper = global_upper - full_rank - known_dimension
    if raw_residual_upper < 0:
        raise F2Error(
            "condition-image rank contradicts the supplied global dimension upper bound"
        )
    all_blocks = tuple(range(len(blocks)))

    independent_labels = []
    accepted_rows: list[list[int]] = []
    current = 0
    for label, row in zip(labels, rows):
        rank = f2_rank_rows([*accepted_rows, row], cursor)
        if rank > current:
            independent_labels.append(label)
            accepted_rows.append(row)
            current = rank

    greedy = []
    chosen: list[int] = []
    remaining = list(all_blocks)
    current_rank = 0
    while remaining:
        scored = []
        for block in remaining:
            rank = _rank_with_blocks(rows, blocks, tuple([*chosen, block]))
            scored.append((rank - current_rank, -block, rank, block))
        gain, _tie, rank, selected = max(scored)
        chosen.append(selected)
        remaining.remove(selected)
        current_rank = rank
        greedy.append(
            {
                "step": len(chosen),
                "place": blocks[selected]["place"],
                "rank_gain": gain,
                "cumulative_certified_condition_rank": rank,
                "raw_relative_selmer_dimension_upper_bound": (
                    global_upper - known_dimension - rank
                ),
            }
        )

    delete_one = []
    for omitted in all_blocks:
        rank = _rank_with_blocks(
            rows, blocks, tuple(index for index in all_blocks if index != omitted)
        )
        delete_one.append(
            {
                "deleted_place": blocks[omitted]["place"],
                "certified_condition_rank": rank,
                "rank_drop": full_rank - rank,
                "raw_relative_selmer_dimension_upper_bound": (
                    global_upper - known_dimension - rank
                ),
            }
        )

    certification = manifest.get("certification")
    if not isinstance(certification, Mapping):
        raise F2Error("certification must be an object")
    required = {
        "global_dimension_upper_bound_certified": certification.get(
            "global_dimension_upper_bound_certified"
        ) is True,
        "witnesses_are_exact_global_squareclasses": certification.get(
            "witnesses_are_exact_global_squareclasses"
        ) is True,
        "witnesses_lie_in_global_ambient_certified": certification.get(
            "witnesses_lie_in_global_ambient_certified"
        ) is True,
        "condition_syndromes_certified": certification.get(
            "condition_syndromes_certified"
        ) is True,
        "condition_blocks_are_necessary_selmer_conditions": certification.get(
            "condition_blocks_are_necessary_selmer_conditions"
        ) is True,
        "known_mw_in_condition_kernel_certified": certification.get(
            "known_mw_in_condition_kernel_certified"
        ) is True,
    }
    method = str(certification.get("method", "")).strip()
    hypothesis_value = certification.get("hypothesis")
    hypothesis = None if hypothesis_value is None else str(hypothesis_value).strip()
    certified = all(required.values()) and bool(method)
    parity_certified = certification.get("residual_dimension_parity_certified") is True
    parity_value = manifest.get("residual_selmer_dimension_parity")
    if parity_certified:
        if parity_value is None or int(parity_value) not in (0, 1):
            raise F2Error("a certified residual parity must be 0 or 1")
        residual_parity = int(parity_value)
        residual_upper = raw_residual_upper - (
            (raw_residual_upper - residual_parity) & 1
        )
        if residual_upper < 0:
            raise F2Error("residual parity contradicts the known Selmer lower bound")
    else:
        residual_parity = None
        residual_upper = raw_residual_upper

    for record in [*greedy, *delete_one]:
        raw = int(record["raw_relative_selmer_dimension_upper_bound"])
        if certified:
            record["parity_sharpened_relative_selmer_dimension_upper_bound"] = (
                raw
                if residual_parity is None
                else raw - ((raw - residual_parity) & 1)
            )
        else:
            record["raw_relative_selmer_dimension_upper_bound"] = None
            record["parity_sharpened_relative_selmer_dimension_upper_bound"] = None

    def closes(rank: int) -> bool:
        raw = global_upper - known_dimension - rank
        if raw < 0:
            raise F2Error("condition rank exceeds the global residual bound")
        if residual_parity is None:
            return raw == 0
        return raw - ((raw - residual_parity) & 1) == 0

    minimum_closing_cut = None
    subsets_examined = 0
    complete_through = -1
    truncated = False
    if certified and residual_upper == 0:
        limit = min(maximum_cut_size, len(blocks))
        for size in range(limit + 1):
            completed_size = True
            for subset in combinations(all_blocks, size):
                if subsets_examined >= maximum_cut_subsets:
                    completed_size = False
                    truncated = True
                    break
                subsets_examined += 1
                if closes(_rank_with_blocks(rows, blocks, subset)):
                    minimum_closing_cut = {
                        "size": size,
                        "places": [blocks[index]["place"] for index in subset],
                        "minimality_proved": True,
                    }
                    break
            if minimum_closing_cut is not None:
                complete_through = size
                break
            if not completed_size:
                break
            complete_through = size

    if certified and residual_upper == 0:
        status = (
            "CERTIFIED_RELATIVE_2SELMER_QUOTIENT_ZERO_UNDER_HYPOTHESIS"
            if hypothesis
            else "CERTIFIED_RELATIVE_2SELMER_QUOTIENT_ZERO"
        )
    elif certified:
        status = (
            "CERTIFIED_RELATIVE_2SELMER_DIMENSION_UPPER_BOUND_UNDER_HYPOTHESIS"
            if hypothesis
            else "CERTIFIED_RELATIVE_2SELMER_DIMENSION_UPPER_BOUND"
        )
    else:
        status = "INCOMPLETE_WITNESS_BOUND"

    return {
        "schema": OUTPUT_SCHEMA,
        "protocol": PROTOCOL,
        "case_id": manifest.get("case_id"),
        "status": status,
        "global_ambient_dimension_upper_bound": global_upper,
        "known_mw_dimension": known_dimension,
        "witness_count": len(rows),
        "condition_target_dimension": cursor,
        "certified_condition_image_rank": full_rank if certified else None,
        "independent_condition_witness_labels": independent_labels,
        "witnesses": [
            {"label": label, "generator": generator, "condition_syndrome": row}
            for label, generator, row in zip(labels, generators, rows)
        ],
        "relative_selmer_dimension_upper_bound_raw": (
            raw_residual_upper if certified else None
        ),
        "certified_residual_dimension_parity": (
            residual_parity if parity_certified else None
        ),
        "relative_selmer_dimension_upper_bound": residual_upper if certified else None,
        "greedy_place_order": greedy,
        "delete_one_place_ranks": delete_one,
        "minimum_closing_place_cut": minimum_closing_cut,
        "minimum_cut_search_complete_through_size": complete_through,
        "minimum_cut_subsets_examined": subsets_examined,
        "minimum_cut_search_truncated_by_subset_budget": truncated,
        "certification": {
            **required,
            "residual_dimension_parity_certified": parity_certified,
            "method": method or None,
            "hypothesis": hypothesis or None,
            "all_requirements_met": certified,
        },
        "claim_boundary": (
            "The bound uses only a global dimension upper bound and a lower "
            "bound on the rank of the norm/local obstruction map. It neither "
            "enumerates the global ambient space nor treats missing witnesses "
            "as absent squareclasses. Auxiliary fingerprints that are not "
            "necessary Selmer conditions are ineligible."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--maximum-cut-size", type=int, default=8)
    parser.add_argument("--maximum-cut-subsets", type=int, default=1_000_000)
    args = parser.parse_args()
    manifest = json.loads(args.input.read_text())
    if not isinstance(manifest, Mapping):
        raise F2Error("top-level JSON value must be an object")
    result = dict(
        audit(
            manifest,
            maximum_cut_size=args.maximum_cut_size,
            maximum_cut_subsets=args.maximum_cut_subsets,
        )
    )
    result["input"] = str(args.input)
    result["input_sha256"] = sha256(args.input.read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|status={result['status']}"
        f"|condition_rank={result['certified_condition_image_rank']}"
        f"|residual_upper={result['relative_selmer_dimension_upper_bound']}"
        f"|output={args.output}"
    )


if __name__ == "__main__":
    main()
