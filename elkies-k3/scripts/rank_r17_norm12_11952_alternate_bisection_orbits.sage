#!/usr/bin/env sage-python
"""Rank all 39,147 alternate-Q80 bisection classes by exact group-law cost."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
import runpy

from sage.all import QQ, ZZ, matrix, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
ENUMERATOR = ROOT / "elkies-k3/scripts/enumerate_rootless_bisection_orbits.sage"
ORBIT_CERTIFICATE = ROOT / "artifacts/generated-results/elkies-k3-q80-alternate-rootless-bisection-orbits.json"
HISTORICAL_FRAME = ROOT / "artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json"
DIRECT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
DIRECT_08F72 = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit08f72-direct-fibration-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisection-priority-v1.json"
TABLE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisection-priority-v1.tsv"
OUTPUT_08F72 = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-08f72-alternate-bisection-priority-v1.json"
TABLE_08F72 = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-08f72-alternate-bisection-priority-v1.tsv"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def parity_mask(value) -> int:
    return sum((int(entry) % 2) << index for index, entry in enumerate(value))


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
    parser.add_argument(
        "--source-label",
        choices=("norm12-orbit-11952", "norm12-orbit-08f72"),
        default="norm12-orbit-11952",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--table-output", type=Path)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    is_primary = arguments.source_label == "norm12-orbit-11952"
    direct_path = DIRECT if is_primary else DIRECT_08F72
    output = arguments.output or (OUTPUT if is_primary else OUTPUT_08F72)
    table_output = arguments.table_output or (TABLE if is_primary else TABLE_08F72)

    orbit_certificate = json.loads(ORBIT_CERTIFICATE.read_text())
    historical = json.loads(HISTORICAL_FRAME.read_text())
    direct = json.loads(direct_path.read_text())
    if orbit_certificate["status"] != "PASS_ALTERNATE_ROOTLESS_LATTICE_BISECTION_ORBITS":
        raise ArithmeticError("alternate bisection lattice certificate is not complete")
    if direct["sections"]["status"] != "PASS_EXACT_SATURATED_RANK17_BASIS":
        raise ArithmeticError("alternate equation section basis is not saturated")

    historical_frame = matrix(ZZ, historical["rootless_frame"])
    direct_frame = matrix(ZZ, direct["frame_certificate"]["frame_gram"])
    historical_to_direct_isometry = matrix(
        ZZ, direct["frame_certificate"]["integral_isometry_to_alternate_Q80"]
    )
    section_coordinates = matrix(
        ZZ, direct["sections"]["coordinate_matrix_in_compiled_frame"]
    )
    if historical_to_direct_isometry * historical_frame * historical_to_direct_isometry.transpose() != direct_frame:
        raise ArithmeticError("historical-to-direct alternate isometry changed")
    if abs(section_coordinates.det()) != 1:
        raise ArithmeticError("alternate section coordinate matrix is not unimodular")

    short_change = matrix(
        ZZ, orbit_certificate["input"]["row_short_basis_change"]
    )
    short_gram = short_change * historical_frame * short_change.transpose()
    short_to_section = (
        short_change * historical_to_direct_isometry.inverse() * section_coordinates.inverse()
    )
    if any(value.denominator() != 1 for value in short_to_section):
        raise ArithmeticError("short-to-section transport is not integral")
    short_to_section = matrix(ZZ, short_to_section)

    enumerator = runpy.run_path(str(ENUMERATOR))
    streaming_short_vectors = enumerator["streaming_short_vectors"]

    def representative_key(short_tuple):
        section_vector = vector(ZZ, short_tuple) * short_to_section
        return score(section_vector)

    streaming = streaming_short_vectors(
        short_gram, bound=10, representative_key=representative_key
    )
    masks_by_norm = streaming["masks_by_norm"]
    excluded = masks_by_norm[2] | masks_by_norm[6]
    candidates = set(streaming["representatives"]) - excluded
    if len(candidates) != 39147:
        raise ArithmeticError(f"expected 39147 alternate classes, obtained {len(candidates)}")

    rows = []
    cost_histogram = Counter()
    for orbit in candidates:
        short_vector = vector(ZZ, streaming["representatives"][orbit])
        historical_vector = short_vector * short_change
        direct_vector_qq = historical_vector * historical_to_direct_isometry.inverse()
        section_vector_qq = direct_vector_qq * section_coordinates.inverse()
        if any(value.denominator() != 1 for value in direct_vector_qq) or any(
            value.denominator() != 1 for value in section_vector_qq
        ):
            raise ArithmeticError("alternate bisection transport is not integral")
        direct_vector = vector(ZZ, direct_vector_qq)
        section_vector = vector(ZZ, section_vector_qq)
        if historical_vector * historical_frame * historical_vector != 10:
            raise ArithmeticError("selected alternate representative has wrong norm")
        if section_vector * matrix(ZZ, direct["sections"]["height_gram"]) * section_vector != 10:
            raise ArithmeticError("selected section word has wrong height")
        row_score = score(section_vector)
        cost_histogram[row_score[0]] += 1
        rows.append(
            {
                "orbit_mask": orbit,
                "orbit_hex": f"0x{orbit:05x}",
                "short_basis_w": tuple(map(int, short_vector)),
                "historical_alternate_w": tuple(map(int, historical_vector)),
                "direct_alternate_w": tuple(map(int, direct_vector)),
                "section_basis_w": tuple(map(int, section_vector)),
                "minimal_unoriented_count": streaming["unoriented_multiplicities"][orbit],
                "group_addition_upper_bound": row_score[0],
                "support_count": row_score[1],
                "maximum_absolute_coefficient": row_score[2],
                "coefficient_l1": row_score[3],
                "_score": row_score,
            }
        )
    rows.sort(key=lambda row: (row["_score"], row["orbit_mask"]))
    for rank, row in enumerate(rows, start=1):
        row["priority_rank"] = rank

    fields = [
        "priority_rank",
        "orbit_mask",
        "orbit_hex",
        "group_addition_upper_bound",
        "support_count",
        "maximum_absolute_coefficient",
        "coefficient_l1",
        "minimal_unoriented_count",
        "section_basis_w",
        "direct_alternate_w",
        "historical_alternate_w",
        "short_basis_w",
    ]
    lines = ["\t".join(fields)]
    for row in rows:
        serialized = dict(row)
        for key in (
            "section_basis_w",
            "direct_alternate_w",
            "historical_alternate_w",
            "short_basis_w",
        ):
            serialized[key] = entries(serialized[key])
        lines.append("\t".join(str(serialized[field]) for field in fields))
    table_text = "\n".join(lines) + "\n"

    payload = {
        "schema": (
            "elkies-k3.r17-norm12-11952-alternate-bisection-priority.v1"
            if is_primary
            else "elkies-k3.r17-norm12-08f72-alternate-bisection-priority.v1"
        ),
        "status": "PASS_EXACT_COMPLETE_ALTERNATE_BISECTION_EQUATION_PRIORITY",
        "class_count": len(rows),
        "minimum_norm": 10,
        "representative_selection": {
            "complete_norm_ten_shell": True,
            "score": [
                "group_addition_upper_bound",
                "support_count",
                "maximum_absolute_coefficient",
                "coefficient_l1",
                "lexicographic_section_basis_word",
            ],
            "cost_histogram": {
                str(cost): count for cost, count in sorted(cost_histogram.items())
            },
        },
        "coordinate_chain": (
            "short -> historical alternate -> direct compiled frame -> saturated equation section basis"
        ),
        "priority_table": relative(table_output),
        "priority_table_sha256": hashlib.sha256(table_text.encode()).hexdigest(),
        "inputs": {
            relative(path): digest(path)
            for path in (ENUMERATOR, ORBIT_CERTIFICATE, HISTORICAL_FRAME, direct_path)
        },
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "required_features": ["exact lattice arithmetic", "complete Fincke-Pohst traversal"],
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "rank_r17_norm12_11952_alternate_bisection_orbits.sage"
            + ("" if is_primary else " --source-label norm12-orbit-08f72")
        ),
        "proof_boundary": (
            "This is a complete exact equation-cost ranking of all 39147 alternate "
            "translation classes. It does not construct their equations or hash their "
            "quadratic extensions."
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not table_output.exists() or table_output.read_text() != table_text:
            raise ArithmeticError("stored alternate priority table differs from replay")
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored alternate priority artifact differs from replay")
    else:
        table_output.parent.mkdir(parents=True, exist_ok=True)
        table_output.write_text(table_text)
        output.write_text(serialized)
    print(
        "ALTBISECTPRIORITY|classes=39147|min_cost={}|max_cost={}|output={}".format(
            rows[0]["group_addition_upper_bound"],
            rows[-1]["group_addition_upper_bound"],
            relative(output),
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
