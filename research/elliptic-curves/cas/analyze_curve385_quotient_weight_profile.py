#!/usr/bin/env python3
"""Certify the quotient-weight profile of the frozen curve-385 blind ledger.

This is a post-search structural analysis.  It uses only the exact integral
coordinates and chart provenance already stored in the blind ledger; it does
not rerun a point search or load the public rank-29 fixture.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
from typing import Any, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
BLIND = ART / "curve385_iterated_half_lattice_blind_v1.json"
OUTPUT = ART / "curve385_quotient_weight_profile_v1.json"

EXPECTED_BLIND_SHA256 = "356001898f738f607d984e081663a015825e11de0c606d35055af156eb2d7502"
GENERIC_RANK = 17
INITIAL_RANK = 20
FINAL_RANK = 29
OLD_CLASS_COUNT = 43


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def canonical_point(record: dict[str, str]) -> tuple[Fraction, Fraction, int]:
    x = Fraction(record["x"])
    y = Fraction(record["y"])
    if (x, -y) < (x, y):
        return x, -y, -1
    return x, y, 1


def point_key(record: dict[str, str]) -> str:
    x, y, unused_sign = canonical_point(record)
    return f"{x}|{y}"


def rational_rank(rows: Iterable[Sequence[int]]) -> int:
    matrix = [list(map(Fraction, row)) for row in rows if any(row)]
    if not matrix:
        return 0
    column_count = len(matrix[0])
    if any(len(row) != column_count for row in matrix):
        raise ArithmeticError("rank matrix is ragged")
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [entry / scale for entry in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or not matrix[row][column]:
                continue
            scale = matrix[row][column]
            matrix[row] = [
                entry - scale * pivot_entry
                for entry, pivot_entry in zip(matrix[row], matrix[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return pivot_row


def gf2_independent(words: Sequence[int]) -> bool:
    pivots: list[int] = []
    for word in words:
        reduced = int(word)
        for pivot in pivots:
            reduced = min(reduced, reduced ^ pivot)
        if not reduced:
            return False
        pivots.append(reduced)
        pivots.sort(reverse=True)
    return True


def source_priority(source: str) -> int | None:
    prefix = "iteration:1:priority:"
    if not source.startswith(prefix):
        return None
    return int(source[len(prefix) :])


def build_payload() -> dict[str, Any]:
    if digest(BLIND) != EXPECTED_BLIND_SHA256:
        raise ArithmeticError("the frozen curve-385 blind ledger changed")
    blind = json.loads(BLIND.read_text())
    if blind.get("status") != "STOPPED_AT_DECLARED_LIFT_LIMIT":
        raise ArithmeticError("the curve-385 blind ledger is not frozen at its limit")
    if blind["blindness_boundary"]["public_rank29_fixture_loaded"] is not False:
        raise ArithmeticError("the blind ledger crossed its public-fixture boundary")
    if len(blind["iterations"]) != 1:
        raise ArithmeticError("the frozen ledger stopped having one lift iteration")

    initial = blind["initial_transition"]
    iteration = blind["iterations"][0]
    if (initial["rank_before"], initial["rank_after"]) != (GENERIC_RANK, INITIAL_RANK):
        raise ArithmeticError("the primitive M17-to-M20 transition changed")
    if initial["discovered_group_saturation"]["events"]:
        raise ArithmeticError("the initial transition acquired a saturation event")
    if (iteration["basis_rank_before"], iteration["basis_rank_after"]) != (
        INITIAL_RANK,
        FINAL_RANK,
    ):
        raise ArithmeticError("the M20-to-M29 transition changed")
    if iteration["finite_index_saturation_event_count"]:
        raise ArithmeticError("the lift round acquired a finite-index event")
    if iteration["basis_before"] != iteration["basis_after"][:INITIAL_RANK]:
        raise ArithmeticError("M20 is no longer the first direct summand of the M29 basis")
    if blind["curve"]["generic_points"] != iteration["basis_before"][:GENERIC_RANK]:
        raise ArithmeticError("M17 is no longer the first direct summand of the M20 basis")

    coordinate_rows: dict[str, list[int]] = {}
    for index, point in enumerate(iteration["basis_after"]):
        key = point_key(point)
        unused_x, unused_y, sign = canonical_point(point)
        row = [0] * FINAL_RANK
        row[index] = sign
        coordinate_rows[key] = row
    saturation = iteration["discovered_group_saturation"]
    if saturation["status"] != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
        raise ArithmeticError("the discovered-group classifier is not exact")
    for relation in saturation["exact_integral_relations"]:
        coordinates = list(map(int, relation["coordinates"]))
        if len(coordinates) != FINAL_RANK:
            raise ArithmeticError("an exact coordinate row has the wrong length")
        coordinate_rows[point_key(relation["point"])] = coordinates

    cover_by_priority = {
        int(row["priority"]): row for row in iteration["cover_records"]
    }
    if len(cover_by_priority) != 301:
        raise ArithmeticError("the frozen lift census stopped having 301 charts")
    if set(cover_by_priority) != set(range(1, 302)):
        raise ArithmeticError("the frozen chart priorities are not contiguous")

    discovery_rows: list[dict[str, Any]] = []
    for discovery in blind["discoveries"]:
        priorities = sorted(
            priority
            for source in discovery["sources"]
            if (priority := source_priority(source)) is not None
        )
        if not priorities:
            continue
        key = point_key(discovery["point"])
        if key not in coordinate_rows:
            raise ArithmeticError("an iteration discovery lacks exact final-basis coordinates")
        words = sorted(
            {int(cover_by_priority[priority]["quotient_word"]) for priority in priorities}
        )
        discovery_rows.append(
            {
                "key": key,
                "priorities": priorities,
                "words": words,
                "coordinates_mod_M20": coordinate_rows[key][INITIAL_RANK:],
            }
        )

    def quotient_rank_for_words(words: Iterable[int]) -> int:
        accepted = set(map(int, words))
        return rational_rank(
            row["coordinates_mod_M20"]
            for row in discovery_rows
            if accepted.intersection(row["words"])
        )

    exact_weight_rows = []
    cumulative_weight_rows = []
    for weight in range(1, 4):
        exact_words = {word for word in range(1, 8) if word.bit_count() == weight}
        cumulative_words = {word for word in range(1, 8) if word.bit_count() <= weight}
        exact_covers = [
            row
            for row in iteration["cover_records"]
            if int(row["quotient_word"]) in exact_words
        ]
        cumulative_covers = [
            row
            for row in iteration["cover_records"]
            if int(row["quotient_word"]) in cumulative_words
        ]
        exact_discoveries = [
            row for row in discovery_rows if exact_words.intersection(row["words"])
        ]
        cumulative_discoveries = [
            row for row in discovery_rows if cumulative_words.intersection(row["words"])
        ]
        exact_weight_rows.append(
            {
                "weight": weight,
                "quotient_word_count": len(exact_words),
                "chart_count": len(exact_covers),
                "charts_with_finite_points": sum(
                    bool(row["search"]["finite_curve_points"]) for row in exact_covers
                ),
                "finite_point_occurrences": sum(
                    len(row["search"]["finite_curve_points"]) for row in exact_covers
                ),
                "distinct_iteration_discoveries": len(exact_discoveries),
                "quotient_rank_over_M20": quotient_rank_for_words(exact_words),
            }
        )
        cumulative_weight_rows.append(
            {
                "maximum_weight": weight,
                "quotient_word_count": len(cumulative_words),
                "chart_count": len(cumulative_covers),
                "charts_with_finite_points": sum(
                    bool(row["search"]["finite_curve_points"])
                    for row in cumulative_covers
                ),
                "finite_point_occurrences": sum(
                    len(row["search"]["finite_curve_points"])
                    for row in cumulative_covers
                ),
                "distinct_iteration_discoveries": len(cumulative_discoveries),
                "quotient_rank_over_M20": quotient_rank_for_words(cumulative_words),
            }
        )

    word_profile = []
    for word in range(1, 8):
        covers = [
            row for row in iteration["cover_records"] if row["quotient_word"] == word
        ]
        discoveries = [row for row in discovery_rows if word in row["words"]]
        word_profile.append(
            {
                "word": word,
                "binary": f"{word:03b}",
                "weight": word.bit_count(),
                "chart_count": len(covers),
                "charts_with_finite_points": sum(
                    bool(row["search"]["finite_curve_points"]) for row in covers
                ),
                "finite_point_occurrences": sum(
                    len(row["search"]["finite_curve_points"]) for row in covers
                ),
                "distinct_iteration_discoveries": len(discoveries),
                "quotient_rank_over_M20": quotient_rank_for_words({word}),
            }
        )

    event_profile = []
    for event in saturation["events"]:
        if event["type"] != "NEW_Q_INDEPENDENT_DIRECTION":
            raise ArithmeticError("the frozen lift round has a non-rank event")
        priorities = sorted(
            priority
            for source in event["sources"]
            if (priority := source_priority(source)) is not None
        )
        words = sorted({cover_by_priority[priority]["quotient_word"] for priority in priorities})
        event_profile.append(
            {
                "basis_rank_after": int(event["basis_rank_after"]),
                "source_priorities": priorities,
                "quotient_words": words,
                "minimum_quotient_weight": min(word.bit_count() for word in words),
            }
        )
    event_profile.sort(key=lambda row: min(row["source_priorities"]))
    if len(event_profile) != FINAL_RANK - INITIAL_RANK:
        raise ArithmeticError("the lift round stopped having nine basis extensions")

    quotient_bases = [
        basis for basis in combinations(range(1, 8), 3) if gf2_independent(basis)
    ]
    if len(quotient_bases) != 28:
        raise ArithmeticError("the unordered GL(3,2) basis census changed")
    weight_one_histogram = Counter(
        quotient_rank_for_words(basis) for basis in quotient_bases
    )
    omitted_rows = []
    successful_weight_two_bases = 0
    for omitted_word in range(1, 8):
        rank = quotient_rank_for_words(set(range(1, 8)) - {omitted_word})
        basis_count = sum(
            basis[0] ^ basis[1] ^ basis[2] == omitted_word
            for basis in quotient_bases
        )
        if rank == FINAL_RANK - INITIAL_RANK:
            successful_weight_two_bases += basis_count
        omitted_rows.append(
            {
                "omitted_physical_word": omitted_word,
                "omitted_physical_word_binary": f"{omitted_word:03b}",
                "basis_count": basis_count,
                "weight_at_most_two_quotient_rank_over_M20": rank,
            }
        )

    twelve_bit_schedule = []
    cumulative_word_count = 0
    for maximum_weight in range(1, 5):
        cumulative_word_count += comb(12, maximum_weight)
        twelve_bit_schedule.append(
            {
                "maximum_weight": maximum_weight,
                "quotient_word_count": cumulative_word_count,
                "chart_count": OLD_CLASS_COUNT * cumulative_word_count,
            }
        )
    twelve_bit_schedule.append(
        {
            "maximum_weight": "all_nonzero",
            "quotient_word_count": (1 << 12) - 1,
            "chart_count": OLD_CLASS_COUNT * ((1 << 12) - 1),
        }
    )

    if cumulative_weight_rows[0]["quotient_rank_over_M20"] != 7:
        raise ArithmeticError("weight one stopped recovering quotient rank seven")
    if cumulative_weight_rows[1]["quotient_rank_over_M20"] != 9:
        raise ArithmeticError("weight at most two stopped recovering all nine directions")
    if cumulative_weight_rows[2]["quotient_rank_over_M20"] != 9:
        raise ArithmeticError("the full cylinder stopped recovering quotient rank nine")
    if any(row["minimum_quotient_weight"] > 2 for row in event_profile):
        raise ArithmeticError("a basis extension now requires quotient weight three")

    return {
        "schema": "elliptic-curves.curve385-quotient-weight-profile.v1",
        "status": "PASS_EXACT_POSTHOC_WEIGHT_PROFILE",
        "input": {
            "path": relative(BLIND),
            "sha256": EXPECTED_BLIND_SHA256,
            "public_rank29_fixture_loaded": False,
        },
        "primitive_split": {
            "generic_rank": GENERIC_RANK,
            "initial_discovered_rank": INITIAL_RANK,
            "new_complement_rank": INITIAL_RANK - GENERIC_RANK,
            "initial_basis_is_generic_basis_followed_by_complement": True,
            "initial_finite_index_saturation_event_count": 0,
            "interpretation": "M20/2M20 = M17/2M17 direct-sum F2^3 for the discovered subgroup",
        },
        "cylinder": {
            "distinguished_generic_class_count": OLD_CLASS_COUNT,
            "all_word_chart_count_including_previously_searched_zero_slice": 344,
            "fresh_nonzero_word_chart_count": 301,
            "exact_weight_profile": exact_weight_rows,
            "cumulative_weight_profile": cumulative_weight_rows,
            "word_profile": word_profile,
            "basis_extension_events_in_search_priority_order": event_profile,
            "all_nine_basis_extensions_have_quotient_weight_at_most_two": True,
        },
        "basis_dependence_audit": {
            "unordered_F2_basis_count": len(quotient_bases),
            "weight_one_rank_histogram": {
                str(rank): count for rank, count in sorted(weight_one_histogram.items())
            },
            "natural_basis_weight_one_rank": quotient_rank_for_words({1, 2, 4}),
            "best_weight_one_rank": max(weight_one_histogram),
            "best_weight_one_unordered_bases": [
                [f"{word:03b}" for word in basis]
                for basis in quotient_bases
                if quotient_rank_for_words(basis) == max(weight_one_histogram)
            ],
            "weight_at_most_two_omitted_word_profile": omitted_rows,
            "weight_at_most_two_full_rank_basis_count": successful_weight_two_bases,
            "weight_at_most_two_deficient_basis_count": len(quotient_bases)
            - successful_weight_two_bases,
        },
        "rank32_scale_projection": {
            "quotient_bit_count": 12,
            "old_class_count": OLD_CLASS_COUNT,
            "staged_chart_counts": twelve_bit_schedule,
        },
        "claim_boundary": [
            "This is an exact coordinate-rank analysis of a completed blind ledger, not a new point search.",
            "The weight profile was inspected after the 301-chart outcome and is therefore structural posthoc evidence, not a prospective success probability.",
            "Hamming weight depends on the chosen quotient basis; the complete GL(3,2) audit records that dependence.",
            "No bounded miss is a rank upper bound, and no saturation in the unknown full E(Q) is claimed.",
        ],
        "input_hashes": {
            relative(BLIND): EXPECTED_BLIND_SHA256,
            relative(Path(__file__)): digest(Path(__file__)),
        },
        "generator": relative(Path(__file__)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    encoded = canonical_bytes(build_payload())
    if args.check:
        if not args.output.exists() or args.output.read_bytes() != encoded:
            raise SystemExit(f"stale or missing artifact: {args.output}")
        print(f"C385WEIGHT|status=PASS|sha256={sha256(encoded).hexdigest()}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(encoded)
    print(f"C385WEIGHT|status=WROTE|sha256={sha256(encoded).hexdigest()}")


if __name__ == "__main__":
    main()
