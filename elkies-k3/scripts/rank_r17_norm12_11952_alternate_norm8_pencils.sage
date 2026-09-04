#!/usr/bin/env sage-python
"""Enumerate and equation-rank the alternate-Q80 norm-eight pencil classes.

For the rootless alternate fibration write ``NS=U+(-M)`` with fibre
``F=(1,0,0)`` and zero ``O=(-1,1,0)``.  A degree-two isotropic class that is
disjoint from the zero has the form

    D_w=(2,2,w),                    (w,w)=8.

Translation by the section ``x`` sends ``w`` to ``w+2x``.  The class is
nonnegative on every section precisely when the coset ``w+2M`` has minimum
norm eight.  This script enumerates those cosets completely and chooses the
cheapest norm-eight trace in the saturated equation section basis.

This is the finite lattice input for product-character bisection inversion.
It does not construct a pencil equation or compare a branch squareclass.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
import runpy

from sage.all import ZZ, matrix, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
ENUMERATOR = ROOT / "elkies-k3/scripts/enumerate_rootless_bisection_orbits.sage"
HISTORICAL_FRAME = ROOT / "artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json"
DIRECT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
DIRECT_08F72 = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit08f72-direct-fibration-v1.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.json"
DEFAULT_TABLE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv"
OUTPUT_08F72 = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-08f72-alternate-norm8-pencil-priority-v1.json"
TABLE_08F72 = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-08f72-alternate-norm8-pencil-priority-v1.tsv"


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
    output = arguments.output or (DEFAULT_OUTPUT if is_primary else OUTPUT_08F72)
    table_output = arguments.table_output or (DEFAULT_TABLE if is_primary else TABLE_08F72)

    historical = json.loads(HISTORICAL_FRAME.read_text())
    direct = json.loads(direct_path.read_text())
    historical_frame = matrix(ZZ, historical["rootless_frame"])
    direct_frame = matrix(ZZ, direct["frame_certificate"]["frame_gram"])
    historical_to_direct = matrix(
        ZZ, direct["frame_certificate"]["integral_isometry_to_alternate_Q80"]
    )
    section_coordinates = matrix(
        ZZ, direct["sections"]["coordinate_matrix_in_compiled_frame"]
    )
    if historical_to_direct * historical_frame * historical_to_direct.transpose() != direct_frame:
        raise ArithmeticError("historical/direct alternate-frame isometry changed")
    if abs(section_coordinates.det()) != 1:
        raise ArithmeticError("equation section coordinates are not unimodular")

    # Work in the same deterministic LLL row basis as the complete alternate
    # rational-bisection census.
    short_change = historical_frame.LLL_gram().transpose()
    short_gram = short_change * historical_frame * short_change.transpose()
    short_to_section_qq = (
        short_change * historical_to_direct.inverse() * section_coordinates.inverse()
    )
    if any(value.denominator() != 1 for value in short_to_section_qq):
        raise ArithmeticError("short-to-section coordinate map is not integral")
    short_to_section = matrix(ZZ, short_to_section_qq)

    enumerator = runpy.run_path(str(ENUMERATOR))
    streaming_short_vectors = enumerator["streaming_short_vectors"]

    def representative_key(short_tuple):
        return score(vector(ZZ, short_tuple) * short_to_section)

    streaming = streaming_short_vectors(
        short_gram, bound=8, representative_key=representative_key
    )
    masks_by_norm = streaming["masks_by_norm"]
    excluded = masks_by_norm[2] | masks_by_norm[4] | masks_by_norm[6]
    candidates = set(streaming["representatives"]) - excluded
    if len(candidates) != 63917:
        raise ArithmeticError(
            f"expected 63917 alternate minimum-norm-eight cosets, got {len(candidates)}"
        )

    rows = []
    cost_histogram = Counter()
    for orbit in candidates:
        short_vector = vector(ZZ, streaming["representatives"][orbit])
        historical_vector = short_vector * short_change
        direct_vector_qq = historical_vector * historical_to_direct.inverse()
        section_vector_qq = direct_vector_qq * section_coordinates.inverse()
        if any(value.denominator() != 1 for value in direct_vector_qq) or any(
            value.denominator() != 1 for value in section_vector_qq
        ):
            raise ArithmeticError("norm-eight coordinate transport is not integral")
        direct_vector = vector(ZZ, direct_vector_qq)
        section_vector = vector(ZZ, section_vector_qq)
        if historical_vector * historical_frame * historical_vector != 8:
            raise ArithmeticError("selected representative has wrong historical norm")
        if section_vector * matrix(ZZ, direct["sections"]["height_gram"]) * section_vector != 8:
            raise ArithmeticError("selected representative has wrong equation height")
        row_score = score(section_vector)
        cost_histogram[row_score[0]] += 1
        rows.append(
            {
                "orbit_mask": int(orbit),
                "orbit_hex": f"0x{orbit:05x}",
                "short_basis_w": tuple(map(int, short_vector)),
                "historical_alternate_w": tuple(map(int, historical_vector)),
                "direct_alternate_w": tuple(map(int, direct_vector)),
                "section_basis_w": tuple(map(int, section_vector)),
                "minimal_unoriented_count": int(
                    streaming["unoriented_multiplicities"][orbit]
                ),
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
            "elkies-k3.r17-norm12-11952-alternate-norm8-pencil-priority.v1"
            if is_primary
            else "elkies-k3.r17-norm12-08f72-alternate-norm8-pencil-priority.v1"
        ),
        "status": "PASS_EXACT_COMPLETE_ALTERNATE_NORM8_PENCIL_PRIORITY",
        "class_count": len(rows),
        "lattice_dictionary": {
            "class": "D_w=(2,2,w)",
            "conditions": ["D_w^2=0", "D_w.F=2", "D_w.O=0", "w.M.w=8"],
            "section_translation": "w -> w+2x",
            "section_intersection": "D_w.S_x=(w-2x).M.(w-2x)/4-2",
            "section_nonnegative_iff": "the coset w+2M has minimum norm at least 8",
            "enumerated_layer": "minimum norm exactly 8",
        },
        "complete_enumeration": {
            "method": "LLL-reduced Fincke-Pohst traversal with exact leaf norms",
            "bound": 8,
            "pari_exact_signed_count_through_bound": int(
                streaming["pari_signed_count"]
            ),
            "signed_shell_counts": {
                str(value): int(streaming["signed_counts"].get(value, 0))
                for value in range(2, 9, 2)
            },
            "parity_cosets_hit_by_shell": {
                str(value): len(masks_by_norm[value]) for value in range(2, 9, 2)
            },
            "excluded_lower_norm_cosets": len(excluded),
            "surviving_minimum_norm_eight_cosets": len(rows),
        },
        "representative_selection": {
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
            "short -> historical alternate -> direct compiled frame -> "
            "saturated equation section basis"
        ),
        "priority_table": relative(table_output),
        "priority_table_sha256": hashlib.sha256(table_text.encode()).hexdigest(),
        "inputs": {
            relative(path): digest(path)
            for path in (ENUMERATOR, HISTORICAL_FRAME, direct_path)
        },
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "required_features": [
                "exact lattice arithmetic",
                "complete Fincke-Pohst traversal",
            ],
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "rank_r17_norm12_11952_alternate_norm8_pencils.sage"
            + ("" if is_primary else " --source-label norm12-orbit-08f72")
        ),
        "proof_boundary": (
            "This is the complete section-nonnegative norm-eight/pole-order-zero "
            "translation layer. It constructs neither chord pencils nor branch matches."
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not table_output.exists() or table_output.read_text() != table_text:
            raise ArithmeticError("stored norm-eight priority table differs from replay")
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored norm-eight priority certificate differs from replay")
    else:
        table_output.parent.mkdir(parents=True, exist_ok=True)
        table_output.write_text(table_text)
        output.write_text(serialized)
    print(
        "ALTNORM8|classes={}|min_cost={}|max_cost={}|output={}".format(
            len(rows),
            rows[0]["group_addition_upper_bound"],
            rows[-1]["group_addition_upper_bound"],
            relative(output),
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
