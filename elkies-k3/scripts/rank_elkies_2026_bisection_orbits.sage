#!/usr/bin/env sage-python
"""Rank the pinned R17 bisection orbits by exact equation-input cost.

The published compact R17 model gives a particularly useful Mordell--Weil
basis: every basis point is supplied either by full coordinates or by a short
quadratic chord from an earlier point.  This script enumerates the complete
norm-ten shell in the repository's reduced short basis, removes the cosets
having norm-six representatives, transports every remaining representative
to the published basis, and keeps the cheapest representative in each of the
39,120 section-translation orbits.

The primary score is a deterministic upper bound for the number of elliptic
group additions needed to form the height-ten trace.  Exact support,
dependency-closure, serialized coordinate-bit, and lattice-coordinate costs
are retained as tie breakers rather than collapsed into a floating score.
The first ``--pool-size`` representatives are then annotated by the exact
norm-four-mask disjointness graph.

This is an equation-priority computation, not a bisection equation or a
quadratic-cover collision certificate.  Use
``construct_elkies_2026_bisections.sage`` for the equation stage.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from sage.all import QuadraticForm, ZZ, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
PINNED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
SHORT_COORDS = ROOT / "elkies-k3/data/lattice/short_vector_basis_coords.txt"
SHORT_GRAM = ROOT / "elkies-k3/data/lattice/short_vector_basis_gram.txt"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
TARGET = ROOT / "artifacts/generated-results/elkies-2026-published-r17-target.json"
FIRST_COVER = ROOT / "elkies-k3/data/fibrations/elkies_2026_rank18_first_cover.json"
PAIRED_COVER = ROOT / "elkies-k3/data/fibrations/elkies_2026_rank19_paired_cover.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-2026-bisection-equation-priority.json"
DEFAULT_TABLE = ROOT / "artifacts/generated-results/elkies-2026-bisection-equation-priority.tsv"
DEFAULT_PAIRS = ROOT / "artifacts/generated-results/elkies-2026-bisection-equation-priority-disjoint-pairs.tsv"


def load_matrix(path: Path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def parity_mask(value) -> int:
    return sum((int(entry) % 2) << index for index, entry in enumerate(value))


def entries(value) -> str:
    return " ".join(str(int(entry)) for entry in value)


def rational_bit_cost(value: str) -> int:
    text = str(value)
    numerator, separator, denominator = text.partition("/")
    cost = abs(int(numerator)).bit_length()
    if separator:
        cost += abs(int(denominator)).bit_length()
    return cost


def section_input_costs(section_data: dict) -> tuple[list[int], list[set[int]]]:
    """Return exact serialized bit costs and recursive chord dependencies."""

    records = section_data["sections"]
    direct_costs = []
    closures: list[set[int]] = []
    for index, record in enumerate(records):
        fields = [record["x_coefficients_low_to_high"]]
        closure = {index}
        if index == 0:
            fields.append(record["y_coefficients_low_to_high"])
        else:
            chord = record["chord"]
            reference = int(chord["reference_basis_index"])
            assert 0 <= reference < index
            closure.update(closures[reference])
            fields.append(chord["slope_coefficients_low_to_high"])
        direct_costs.append(sum(rational_bit_cost(value) for field in fields for value in field))
        closures.append(closure)
    return direct_costs, closures


def binary_scalar_additions(coefficient: int) -> int:
    """Standard double-and-add upper bound; negation has zero addition cost."""

    value = abs(int(coefficient))
    if value <= 1:
        return 0
    return value.bit_length() - 1 + value.bit_count() - 1


def equation_score(
    published_vector: tuple[int, ...], direct_costs: list[int], closures: list[set[int]]
) -> tuple[tuple, dict]:
    support = [index for index, coefficient in enumerate(published_vector) if coefficient]
    dependency_closure: set[int] = set()
    for index in support:
        dependency_closure.update(closures[index])
    scalar_additions = sum(binary_scalar_additions(published_vector[index]) for index in support)
    combination_additions = max(0, len(support) - 1)
    group_additions = scalar_additions + combination_additions
    coordinate_input_bits = sum(direct_costs[index] for index in dependency_closure)
    maximum = max(abs(value) for value in published_vector)
    l1 = sum(abs(value) for value in published_vector)
    # The tuple is the score.  Its components are deliberately kept separate
    # so that changing a proxy weight cannot silently change the mathematical
    # meaning of the ranking.
    key = (
        group_additions,
        len(support),
        len(dependency_closure),
        coordinate_input_bits,
        maximum,
        l1,
        published_vector,
    )
    return key, {
        "group_addition_upper_bound": group_additions,
        "scalar_additions": scalar_additions,
        "combination_additions": combination_additions,
        "support_count": len(support),
        "support_one_based": [index + 1 for index in support],
        "dependency_count": len(dependency_closure),
        "dependency_closure_one_based": [index + 1 for index in sorted(dependency_closure)],
        "coordinate_input_bits": coordinate_input_bits,
        "maximum_absolute_coefficient": maximum,
        "coefficient_l1": l1,
    }


def canonical_scored_orientation(value, short_to_published, direct_costs, closures):
    published = tuple(int(entry) for entry in value * short_to_published)
    choices = []
    for sign in (1, -1):
        candidate = tuple(sign * entry for entry in published)
        key, score = equation_score(candidate, direct_costs, closures)
        choices.append((key, score, candidate, sign))
    return min(choices, key=lambda item: item[0])


def known_cover_audit(ranked_by_mask: dict[int, dict], published_to_short) -> list[dict]:
    audits = []
    sources = [
        (FIRST_COVER, "trace_recovery", "published_basis_vector", "published-first-cover"),
        (PAIRED_COVER, "second_cover", "trace_published_basis_vector", "published-second-cover"),
    ]
    for path, container, key, label in sources:
        if not path.is_file():
            continue
        payload = json.loads(path.read_text())
        published = vector(ZZ, payload[container][key])
        short = published * published_to_short
        assert all(entry in ZZ for entry in short)
        orbit = parity_mask(short)
        row = ranked_by_mask[orbit]
        audits.append({
            "label": label,
            "source": relative(path),
            "orbit_mask": orbit,
            "orbit_hex": f"0x{orbit:05x}",
            "equation_rank": row["equation_rank"],
            "published_basis_vector": [int(entry) for entry in published],
            "best_published_basis_vector": list(row["published_basis_w"]),
            "same_minimal_representative_up_to_sign": (
                tuple(int(entry) for entry in published) == row["published_basis_w"]
                or tuple(-int(entry) for entry in published) == row["published_basis_w"]
            ),
        })
    return audits


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pool-size", type=int, default=1000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--table-output", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--pairs-output", type=Path, default=DEFAULT_PAIRS)
    parser.add_argument("--pari-stack-gb", type=int, default=4)
    arguments = parser.parse_args()
    if arguments.pool_size <= 0:
        parser.error("--pool-size must be positive")

    pinned = load_matrix(PINNED)
    short_coordinates = load_matrix(SHORT_COORDS)
    short_gram = load_matrix(SHORT_GRAM)
    assert short_gram == short_coordinates * pinned * short_coordinates.transpose()
    target = json.loads(TARGET.read_text())
    assert target["status"] == "PASS_EXACT_PUBLISHED_R17_IS_PINNED_R17"
    basis_change = matrix(ZZ, target["pinned_identification"]["basis_change_matrix"])
    orientation = target["pinned_identification"]["gram_identity_orientation"]
    assert orientation == "M^T*Gpub*M=Gpinned"
    short_to_published = short_coordinates * basis_change.transpose()
    published_to_short = short_to_published.inverse()
    assert all(entry.denominator() == 1 for entry in published_to_short.list())
    published_to_short = matrix(ZZ, published_to_short)

    section_data = json.loads(SECTIONS.read_text())
    direct_costs, closures = section_input_costs(section_data)
    assert len(direct_costs) == len(closures) == 17

    coefficients = []
    for row in range(17):
        for column in range(row, 17):
            coefficients.append(
                short_gram[row, row] // 2 if row == column else short_gram[row, column]
            )
    pari.allocatemem(arguments.pari_stack_gb * 1024**3)
    shells = QuadraticForm(ZZ, 17, coefficients).short_vector_list_up_to_length(6, True)
    norm_four_masks = {parity_mask(value) for value in shells[2]}
    excluded_masks = {parity_mask(value) for value in shells[3]}
    assert len(norm_four_masks) == 1311
    assert len(excluded_masks) == 26672

    best: dict[int, tuple] = {}
    multiplicities = Counter()
    for raw_value in shells[5]:
        value = vector(ZZ, raw_value)
        orbit = parity_mask(value)
        if orbit in excluded_masks:
            continue
        multiplicities[orbit] += 1
        key, score, published, sign = canonical_scored_orientation(
            value, short_to_published, direct_costs, closures
        )
        signed_short = tuple(sign * int(entry) for entry in value)
        prior = best.get(orbit)
        candidate = (key, score, published, signed_short)
        if prior is None or candidate[0] < prior[0]:
            best[orbit] = candidate
    assert len(best) == 39120
    assert sum(multiplicities.values()) == 806238

    ranked = []
    for equation_rank, (orbit, (key, score, published, short)) in enumerate(
        sorted(best.items(), key=lambda item: (item[1][0], item[0])), start=1
    ):
        short_vector = vector(ZZ, short)
        pinned_vector = short_vector * short_coordinates
        assert pinned_vector * pinned * pinned_vector == 10
        ranked.append({
            "equation_rank": equation_rank,
            "orbit_mask": orbit,
            "orbit_hex": f"0x{orbit:05x}",
            "published_basis_w": published,
            "pinned_rank17_w": tuple(int(entry) for entry in pinned_vector),
            "short_basis_w": short,
            "minimal_unoriented_count": multiplicities[orbit],
            **score,
            "_sort_key": key,
        })

    pool = ranked[: min(arguments.pool_size, len(ranked))]
    pool_masks = {row["orbit_mask"] for row in pool}
    degrees = Counter()
    pair_rows = []
    by_mask = {row["orbit_mask"]: row for row in pool}
    for left in sorted(pool_masks):
        for delta in norm_four_masks:
            right = left ^ delta
            if left < right and right in pool_masks:
                degrees[left] += 1
                degrees[right] += 1
                pair_rows.append((left, right, delta))
    for row in pool:
        row["disjoint_degree_in_pool"] = degrees[row["orbit_mask"]]
    priority = sorted(
        pool,
        key=lambda row: (
            row["group_addition_upper_bound"],
            -row["disjoint_degree_in_pool"],
            row["_sort_key"][1:],
            row["orbit_mask"],
        ),
    )
    for priority_rank, row in enumerate(priority, start=1):
        row["priority_rank"] = priority_rank
    pair_rows.sort(
        key=lambda item: (
            max(by_mask[item[0]]["equation_rank"], by_mask[item[1]]["equation_rank"]),
            by_mask[item[0]]["equation_rank"] + by_mask[item[1]]["equation_rank"],
            item,
        )
    )

    arguments.table_output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "priority_rank", "equation_rank", "orbit_mask", "orbit_hex",
        "group_addition_upper_bound", "support_count", "dependency_count",
        "coordinate_input_bits", "maximum_absolute_coefficient", "coefficient_l1",
        "disjoint_degree_in_pool", "minimal_unoriented_count", "published_basis_w",
        "pinned_rank17_w", "short_basis_w",
    ]
    with arguments.table_output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in sorted(pool, key=lambda item: item["priority_rank"]):
            serialized = dict(row)
            for key in ("published_basis_w", "pinned_rank17_w", "short_basis_w"):
                serialized[key] = entries(serialized[key])
            writer.writerow({field: serialized[field] for field in fields})

    arguments.pairs_output.parent.mkdir(parents=True, exist_ok=True)
    with arguments.pairs_output.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow([
            "left_orbit_mask", "left_hex", "left_priority_rank", "left_equation_rank",
            "right_orbit_mask", "right_hex", "right_priority_rank", "right_equation_rank",
            "norm_four_difference_mask", "difference_hex",
        ])
        for left, right, delta in pair_rows:
            left_row, right_row = by_mask[left], by_mask[right]
            writer.writerow([
                left, f"0x{left:05x}", left_row["priority_rank"], left_row["equation_rank"],
                right, f"0x{right:05x}", right_row["priority_rank"], right_row["equation_rank"],
                delta, f"0x{delta:05x}",
            ])

    ranked_by_mask = {row["orbit_mask"]: row for row in ranked}
    payload = {
        "schema": "elkies-k3.elkies-2026-bisection-equation-priority.v1",
        "status": "PASS_EXACT_R17_BISECTION_EQUATION_PRIORITY",
        "scope": (
            "Complete exact ranking of the 39,120 surviving lattice orbits by a "
            "published-coordinate equation-input score; no bisection equation or "
            "quadratic-extension collision is asserted here."
        ),
        "inputs": {
            relative(path): digest(path)
            for path in (PINNED, SHORT_COORDS, SHORT_GRAM, SECTIONS, TARGET)
        },
        "score": {
            "sort_key": [
                "group_addition_upper_bound", "support_count", "dependency_count",
                "coordinate_input_bits", "maximum_absolute_coefficient", "coefficient_l1",
                "lexicographic_published_basis_vector",
            ],
            "scalar_rule": (
                "For nonzero n, floor(log2(abs(n)))+popcount(abs(n))-1; then add "
                "support_count-1 to combine the nonzero scalar multiples."
            ),
            "section_representation_rule": (
                "Recursive closure of the exact published full-coordinate/chord records; "
                "the bit cost is the sum of numerator and denominator bit lengths."
            ),
            "orientation": "minimize the same score over w and -w",
        },
        "complete_enumeration": {
            "surviving_translation_orbits": len(ranked),
            "surviving_unoriented_norm_ten_representatives": sum(multiplicities.values()),
            "excluded_norm_six_masks": len(excluded_masks),
            "norm_four_difference_masks": len(norm_four_masks),
            "group_addition_histogram": {
                str(cost): count
                for cost, count in sorted(Counter(row["group_addition_upper_bound"] for row in ranked).items())
            },
        },
        "priority_pool": {
            "size": len(pool),
            "equation_rank_boundary": len(pool),
            "disjoint_pairs_inside_pool": len(pair_rows),
            "table": relative(arguments.table_output),
            "table_sha256": digest(arguments.table_output),
            "disjoint_pairs_table": relative(arguments.pairs_output),
            "disjoint_pairs_table_sha256": digest(arguments.pairs_output),
            "priority_rule": (
                "Within each group-addition cost, prefer larger exact disjoint degree "
                "inside the equation-cheapest pool, then use the remaining equation score."
            ),
        },
        "published_cover_audit": known_cover_audit(ranked_by_mask, published_to_short),
        "proof_boundary": (
            "The shell and mask calculations are complete. The cost is an exact "
            "deterministic arithmetic proxy, not a theorem that coefficient height or "
            "wall-clock reconstruction time is monotone in this order."
        ),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "ELKIES2026BISECTIONPRIORITY|orbits={}|pool={}|disjoint_pairs={}|"
        "status={}|output={}".format(
            len(ranked), len(pool), len(pair_rows), payload["status"], arguments.output
        )
    )


if __name__ == "__main__":
    main()
