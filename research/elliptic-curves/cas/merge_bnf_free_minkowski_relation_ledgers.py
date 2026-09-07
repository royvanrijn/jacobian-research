#!/usr/bin/env python3
"""Merge exact sparse-hypergraph edges from independent Minkowski runs.

A single special-q run can end with a forest even though its edges close
cycles with edges from another run.  The collector therefore stores every
exact partial relation, not only dependencies found locally.  This script
replays those edges in one GF(2) sparse eliminator, retains exact generator
provenance for new cross-run cycles, and emits another relation ledger.

The output is still a relation-collection checkpoint.  It is not a proof that
the factor base generates the class group and is not a 2-Selmer computation.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path

from run_fermigier_rank20_fixedfb_quadratic_specialq import (
    SparseLargePrimeEliminator,
)


PROTOCOL = "BNFFREEMERGE"
SCHEMA = "elliptic-curves.bnf-free-principal-relation-ledger.v1"
CROSS_RUN_SOURCE = "merged_minkowski_cross_run"


def insert_row(pivots: dict[int, int], row: int) -> bool:
    while row:
        pivot = row.bit_length() - 1
        previous = pivots.get(pivot)
        if previous is None:
            pivots[pivot] = row
            return True
        row ^= previous
    return False


def packed_rank(rows: list[int]) -> int:
    pivots: dict[int, int] = {}
    for row in rows:
        insert_row(pivots, row)
    return len(pivots)


def vertex_key(value: object) -> tuple[int, str]:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError("large-prime vertex must be [rational_prime, ideal_hnf]")
    return int(value[0]), str(value[1])


def projective_key(generator: dict) -> tuple[str, ...]:
    values = generator.get("primitive_projective_power_basis")
    if not isinstance(values, list) or len(values) != 3:
        raise ValueError(
            "every partial-relation generator must store three primitive "
            "projective power-basis coordinates"
        )
    return tuple(str(value) for value in values)


def compatible_value(ledger: dict, name: str) -> str:
    """Use canonical JSON so nested records compare without hashability tricks."""

    return json.dumps(ledger.get(name), sort_keys=True, separators=(",", ":"))


def validate_compatible(ledgers: list[dict]) -> None:
    if len(ledgers) < 2:
        raise ValueError("at least two relation ledgers are required")
    compatibility_fields = (
        "schema",
        "curve_preset",
        "factor_hint_certificate",
        "defining_polynomial_ascending",
        "field_discriminant",
        "generator_coordinate_order",
        "factor_base_bound",
        "factor_base_completion",
        "selmer_rational_primes",
        "factor_base",
        "S_columns",
    )
    reference = ledgers[0]
    if reference.get("schema") != SCHEMA:
        raise ValueError("unexpected relation-ledger schema")
    for ledger_index, ledger in enumerate(ledgers):
        if ledger.get("schema") != SCHEMA:
            raise ValueError(f"ledger {ledger_index} has an unexpected schema")
        if ledger.get("large_prime_merge_mode") != "sparse-hypergraph":
            raise ValueError(
                f"ledger {ledger_index} was not collected with sparse-hypergraph merging"
            )
        if not isinstance(ledger.get("generators"), list):
            raise ValueError(f"ledger {ledger_index} has no generator list")
        if not isinstance(ledger.get("partial_relations"), list):
            raise ValueError(
                f"ledger {ledger_index} predates exact partial-relation retention"
            )
        if not isinstance(ledger.get("closed_relations"), list):
            raise ValueError(f"ledger {ledger_index} has no closed-relation list")
        for field in compatibility_fields:
            if compatible_value(ledger, field) != compatible_value(reference, field):
                raise ValueError(
                    f"ledger {ledger_index} is incompatible in field {field}"
                )


def merge_loaded_ledgers(ledgers: list[dict]) -> dict:
    """Return one exact ledger with all cross-run dependencies replayed."""

    validate_compatible(ledgers)
    merged = deepcopy(ledgers[0])
    merged_generators: list[dict] = []
    merged_closed: list[dict] = []
    merged_partial: list[dict] = []
    merged_unresolved: list[dict] = []
    offsets: list[int] = []
    generator_source: dict[int, int] = {}

    for ledger_index, ledger in enumerate(ledgers):
        offset = len(merged_generators)
        offsets.append(offset)
        for local_index, generator in enumerate(ledger["generators"]):
            global_index = offset + local_index
            generator_source[global_index] = ledger_index
            merged_generators.append(deepcopy(generator))
        for relation in ledger["closed_relations"]:
            shifted = deepcopy(relation)
            indices = [int(value) for value in shifted.get("generator_indices", ())]
            if len(indices) != len(set(indices)) or any(
                value < 0 or value >= len(ledger["generators"]) for value in indices
            ):
                raise ValueError(
                    f"ledger {ledger_index} has an invalid closed-relation witness"
                )
            shifted["generator_indices"] = [offset + value for value in indices]
            shifted["source_ledger_index"] = ledger_index
            merged_closed.append(shifted)
        for record in ledger.get("unresolved_cofactors", ()):
            shifted = deepcopy(record)
            local_index = int(shifted["generator_index"])
            if local_index < 0 or local_index >= len(ledger["generators"]):
                raise ValueError(
                    f"ledger {ledger_index} has an invalid unresolved-cofactor witness"
                )
            shifted["generator_index"] = offset + local_index
            shifted["source_ledger_index"] = ledger_index
            merged_unresolved.append(shifted)

    eliminator = SparseLargePrimeEliminator()
    seen_projective: dict[tuple[str, ...], int] = {}
    per_ledger_seen: list[set[tuple[str, ...]]] = [set() for _ in ledgers]
    cross_run_relations: list[dict] = []
    skipped_cross_run_duplicates = 0

    for ledger_index, ledger in enumerate(ledgers):
        offset = offsets[ledger_index]
        for relation_index, relation in enumerate(ledger["partial_relations"]):
            if not isinstance(relation, dict):
                raise ValueError(
                    f"ledger {ledger_index} partial relation {relation_index} is not an object"
                )
            local_generator = int(relation.get("generator_index", -1))
            if local_generator < 0 or local_generator >= len(ledger["generators"]):
                raise ValueError(
                    f"ledger {ledger_index} partial relation {relation_index} "
                    "has an invalid generator"
                )
            key = projective_key(ledger["generators"][local_generator])
            if key in per_ledger_seen[ledger_index]:
                raise ValueError(
                    f"ledger {ledger_index} repeats a projective generator; "
                    "refuse a possible rational-multiple fake cycle"
                )
            per_ledger_seen[ledger_index].add(key)
            if key in seen_projective:
                skipped_cross_run_duplicates += 1
                continue

            global_generator = offset + local_generator
            seen_projective[key] = global_generator
            vertices = [
                vertex_key(value)
                for value in relation.get("large_prime_vertices", ())
            ]
            row = int(relation.get("fb_parity_mask_hex", "0"), 16)
            if row >> len(merged["factor_base"]):
                raise ValueError(
                    f"ledger {ledger_index} partial relation {relation_index} "
                    "exceeds the factor-base dimension"
                )
            shifted = deepcopy(relation)
            shifted["generator_index"] = global_generator
            shifted["source_ledger_index"] = ledger_index
            merged_partial.append(shifted)
            cycle, provenance = eliminator.add(vertices, row, global_generator)
            if cycle is None:
                continue
            source_ledgers = sorted({generator_source[value] for value in provenance})
            if len(source_ledgers) < 2:
                # The source ledger already retained this local dependency.
                continue
            closed = {
                "fb_parity_mask_hex": hex(cycle),
                "generator_indices": sorted(provenance),
                "kind": "unit_dependency" if cycle == 0 else "minkowski_lp_cycle",
                "source": CROSS_RUN_SOURCE,
                "source_ledger_indices": source_ledgers,
            }
            merged_closed.append(closed)
            cross_run_relations.append(closed)

    original_rows = [
        int(relation["fb_parity_mask_hex"], 16)
        for relation in merged_closed
        if relation.get("source") != CROSS_RUN_SOURCE
    ]
    cross_rows = [
        int(relation["fb_parity_mask_hex"], 16)
        for relation in cross_run_relations
    ]
    s_rows = [1 << int(column) for column in merged["S_columns"]]
    rank_before = packed_rank(original_rows + s_rows)
    rank_after = packed_rank(original_rows + cross_rows + s_rows)

    merged["status"] = "exact_merged_minkowski_relations_not_class_group_completion"
    merged["special_ideal_mode"] = "merged-independent-runs"
    merged["large_prime_merge_mode"] = "sparse-hypergraph"
    merged["generators"] = merged_generators
    merged["partial_relations"] = merged_partial
    merged["closed_relations"] = merged_closed
    merged["unresolved_cofactors"] = merged_unresolved
    merged["large_prime_elimination"] = {
        "vertex_count": len(eliminator.vertex_columns),
        "edge_count": eliminator.edge_count,
        "rank": len(eliminator.pivots),
        "dependency_count": eliminator.dependency_count,
        "nullity": eliminator.edge_count - len(eliminator.pivots),
    }
    merged["merged_relation_collection"] = {
        "input_ledger_count": len(ledgers),
        "accepted_projectively_distinct_partial_relation_count": len(merged_partial),
        "skipped_cross_run_projective_duplicate_count": skipped_cross_run_duplicates,
        "new_cross_run_closed_relation_count": len(cross_run_relations),
        "rank_modulo_existing_closed_relations_and_S_before": rank_before,
        "rank_modulo_existing_closed_relations_and_S_after": rank_after,
        "cross_run_quotient_rank_gain": rank_after - rank_before,
        "claim_boundary": (
            "Exact principal-relation replay only; no factor-base generation, "
            "class-group completion, norm-kernel, local-solubility, or Selmer claim."
        ),
    }
    merged.pop("batch_gcd_resolution", None)
    return merged


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--relation-ledger", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    input_paths = [path.resolve() for path in args.relation_ledger]
    if len(set(input_paths)) != len(input_paths):
        raise ValueError("input relation-ledger paths must be distinct")
    if args.output.resolve() in input_paths:
        raise ValueError("write the merged ledger to a distinct path")

    ledgers = []
    sources = []
    for path in args.relation_ledger:
        payload = path.read_bytes()
        ledger = json.loads(payload)
        if not isinstance(ledger, dict):
            raise ValueError(f"relation ledger {path} is not a JSON object")
        ledgers.append(ledger)
        sources.append(
            {
                "path": str(path.resolve()),
                "sha256": sha256(payload).hexdigest(),
            }
        )

    merged = merge_loaded_ledgers(ledgers)
    merged["merged_relation_collection"]["input_ledgers"] = sources
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(merged, indent=2, sort_keys=True) + "\n")
    summary = merged["merged_relation_collection"]
    elimination = merged["large_prime_elimination"]
    print(
        f"{PROTOCOL}|stage=complete|inputs={summary['input_ledger_count']}"
        f"|partial={elimination['edge_count']}|vertices={elimination['vertex_count']}"
        f"|nullity={elimination['nullity']}"
        f"|cross_cycles={summary['new_cross_run_closed_relation_count']}"
        f"|quotient_rank_gain={summary['cross_run_quotient_rank_gain']}"
        f"|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
