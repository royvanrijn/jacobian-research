#!/usr/bin/env sage-python
"""Rank the complete 39,120 hidden-103b2 bisection classes by equation cost."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
import runpy

from sage.all import ZZ, matrix, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
SOURCE_ORBITS = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-orbits-full.tsv"
ENUMERATOR = ROOT / "elkies-k3/scripts/enumerate_rootless_bisection_orbits.sage"
DIRECT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit103b2-direct-fibration-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-103b2-bisection-priority-v1.json"
TABLE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-103b2-bisection-priority-v1.tsv"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def entries(value) -> str:
    return " ".join(str(int(entry)) for entry in value)


def binary_scalar_additions(coefficient: int) -> int:
    value = abs(int(coefficient))
    if value <= 1:
        return 0
    return value.bit_length() - 1 + value.bit_count() - 1


def score(coefficients) -> tuple:
    coefficients = tuple(map(int, coefficients))
    support = tuple(index for index, value in enumerate(coefficients) if value)
    additions = sum(binary_scalar_additions(coefficients[index]) for index in support)
    additions += max(0, len(support) - 1)
    return (
        additions,
        len(support),
        max(map(abs, coefficients)),
        sum(map(abs, coefficients)),
        coefficients,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--table-output", type=Path, default=TABLE)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    direct = json.loads(DIRECT.read_text())
    if direct["sections"]["status"] != "PASS_EXACT_SATURATED_RANK17_BASIS":
        raise ArithmeticError("hidden equation section basis is not saturated")

    section_coordinates = matrix(
        ZZ, direct["sections"]["coordinate_matrix_in_compiled_frame"]
    )
    direct_frame = matrix(ZZ, direct["frame_certificate"]["frame_gram"])
    section_gram = matrix(ZZ, direct["sections"]["height_gram"])
    if abs(section_coordinates.det()) != 1:
        raise ArithmeticError("hidden section coordinate matrix is not unimodular")

    short_change = direct_frame.LLL_gram().transpose()
    short_gram = short_change * direct_frame * short_change.transpose()
    short_to_section = short_change * section_coordinates.inverse()
    if any(value.denominator() != 1 for value in short_to_section):
        raise ArithmeticError("hidden short-to-section transport is not integral")
    short_to_section = matrix(ZZ, short_to_section)
    enumerator = runpy.run_path(str(ENUMERATOR))

    def representative_key(short_tuple):
        return score(vector(ZZ, short_tuple) * short_to_section)

    streaming = enumerator["streaming_short_vectors"](
        short_gram, bound=10, representative_key=representative_key
    )
    excluded = streaming["masks_by_norm"][2] | streaming["masks_by_norm"][6]
    candidates = set(streaming["representatives"]) - excluded
    if len(candidates) != 39120:
        raise ArithmeticError(
            f"expected 39120 complete hidden translation classes, obtained {len(candidates)}"
        )

    rows = []
    cost_histogram = Counter()
    for orbit in candidates:
        short_vector = vector(ZZ, streaming["representatives"][orbit])
        direct_vector = short_vector * short_change
        section_vector = short_vector * short_to_section
        if direct_vector * direct_frame * direct_vector != 10:
            raise ArithmeticError("transported hidden bisection has wrong norm")
        if section_vector * section_gram * section_vector != 10:
            raise ArithmeticError("hidden section word has wrong height")
        row_score = score(section_vector)
        cost_histogram[row_score[0]] += 1
        rows.append({
            "orbit_mask": orbit,
            "orbit_hex": f"0x{orbit:05x}",
            "short_basis_w": tuple(map(int, short_vector)),
            "direct_hidden_w": tuple(map(int, direct_vector)),
            "section_basis_w": tuple(map(int, section_vector)),
            "minimal_unoriented_count": streaming["unoriented_multiplicities"][orbit],
            "group_addition_upper_bound": row_score[0],
            "support_count": row_score[1],
            "maximum_absolute_coefficient": row_score[2],
            "coefficient_l1": row_score[3],
            "_score": row_score,
        })
    rows.sort(key=lambda row: (row["_score"], row["orbit_mask"]))
    for rank, row in enumerate(rows, start=1):
        row["priority_rank"] = rank

    fields = [
        "priority_rank", "orbit_mask", "orbit_hex",
        "group_addition_upper_bound", "support_count",
        "maximum_absolute_coefficient", "coefficient_l1",
        "minimal_unoriented_count", "section_basis_w", "direct_hidden_w", "short_basis_w",
    ]
    lines = ["\t".join(fields)]
    for row in rows:
        serialized = dict(row)
        for key in ("section_basis_w", "direct_hidden_w", "short_basis_w"):
            serialized[key] = entries(serialized[key])
        lines.append("\t".join(str(serialized[field]) for field in fields))
    table_text = "\n".join(lines) + "\n"

    result = {
        "schema": "elkies-k3.r17-norm12-103b2-bisection-priority.v1",
        "status": "PASS_EXACT_COMPLETE_103B2_BISECTION_EQUATION_PRIORITY",
        "class_count": len(rows),
        "minimum_norm": 10,
        "rootless_frame": [list(map(int, row)) for row in direct_frame.rows()],
        "representative_selection": {
            "complete_norm_ten_shell": True,
            "score": [
                "group_addition_upper_bound", "support_count",
                "maximum_absolute_coefficient", "coefficient_l1",
                "lexicographic_section_basis_word",
            ],
            "cost_histogram": {
                str(cost): count for cost, count in sorted(cost_histogram.items())
            },
        },
        "coordinate_chain": (
            "LLL-reduced hidden frame -> hidden direct frame -> saturated hidden equation basis"
        ),
        "priority_table": relative(arguments.table_output),
        "priority_table_sha256": hashlib.sha256(table_text.encode()).hexdigest(),
        "inputs": {
            relative(path): digest(path) for path in (ENUMERATOR, DIRECT)
        },
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "required_features": ["exact lattice arithmetic", "complete Fincke-Pohst traversal"],
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/rank_r17_norm12_103b2_bisection_orbits.sage"
        ),
        "proof_boundary": (
            "This is a complete exact equation-cost ranking of all 39120 hidden-marking "
            "translation classes, choosing the cheapest hidden-section word over the complete norm-ten shell. "
            "It does not construct their hidden equations or hash their covers."
        ),
    }
    serialized_result = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not arguments.table_output.exists() or arguments.table_output.read_text() != table_text:
            raise ArithmeticError("stored hidden priority table differs from replay")
        if not arguments.output.exists() or arguments.output.read_text() != serialized_result:
            raise ArithmeticError("stored hidden priority artifact differs from replay")
    else:
        arguments.table_output.parent.mkdir(parents=True, exist_ok=True)
        arguments.table_output.write_text(table_text)
        arguments.output.write_text(serialized_result)
    print(
        "HIDDEN103B2PRIORITY|classes=39120|min_cost={}|max_cost={}|output={}".format(
            rows[0]["group_addition_upper_bound"],
            rows[-1]["group_addition_upper_bound"],
            relative(arguments.output),
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
