#!/usr/bin/env sage-python
"""Enumerate and equation-rank the hidden-103b2 norm-eight pencil classes.

For a rootless rank-17 Mordell--Weil lattice ``M``, the degree-two isotropic
classes ``D_w=(2,2,w)`` in the singular genus-one-pencil layer are indexed by
cosets ``w+2M`` whose minimum norm is eight.  This script enumerates that
finite layer exactly in the saturated hidden-103b2 equation basis and records
the number of minimum representatives up to sign in every coset.

Those multiplicities are the independent lattice input used to prove that
the known split members exhaust the even-multiplicity part of each pencil
discriminant.  No pencil equation or branch character is constructed here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
import runpy

from sage.all import ZZ, matrix, pari, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
ENUMERATOR = ROOT / "elkies-k3/scripts/enumerate_rootless_bisection_orbits.sage"
DIRECT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-orbit103b2-direct-fibration-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-103b2-norm8-pencil-priority-v1.json"
)
DEFAULT_TABLE = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-103b2-norm8-pencil-priority-v1.tsv"
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def parity_mask(value) -> int:
    return sum((int(entry) & 1) << index for index, entry in enumerate(value))


def entries(value) -> str:
    return " ".join(str(int(entry)) for entry in value)


def score(coefficients) -> tuple:
    oriented = min(tuple(map(int, coefficients)), tuple(-int(x) for x in coefficients))
    return (
        sum(abs(entry) for entry in oriented),
        sum(bool(entry) for entry in oriented),
        max(abs(entry) for entry in oriented),
        oriented,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=DIRECT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--table-output", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    model_path = arguments.model.resolve()
    direct = json.loads(model_path.read_text())
    if direct.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
        raise ArithmeticError("hidden direct equation certificate is not exact")
    if direct["sections"].get("status") != "PASS_EXACT_SATURATED_RANK17_BASIS":
        raise ArithmeticError("hidden equation section basis is not saturated")
    gram = matrix(ZZ, direct["sections"]["height_gram"])
    if gram.nrows() != 17 or gram.ncols() != 17 or gram.det() != 948:
        raise ArithmeticError("unexpected hidden Mordell--Weil lattice")
    if pari(gram).qfminim(2)[0] != 0:
        raise ArithmeticError("hidden Mordell--Weil lattice is not rootless")

    short_change = gram.LLL_gram().transpose()
    if abs(short_change.det()) != 1:
        raise ArithmeticError("LLL coordinate change is not unimodular")
    short_gram = short_change * gram * short_change.transpose()
    enumerator = runpy.run_path(str(ENUMERATOR))

    def representative_key(short_tuple):
        return score(vector(ZZ, short_tuple) * short_change)

    streaming = enumerator["streaming_short_vectors"](
        short_gram, bound=8, representative_key=representative_key
    )
    masks_by_norm = streaming["masks_by_norm"]
    excluded = masks_by_norm[2] | masks_by_norm[4] | masks_by_norm[6]
    candidates = set(streaming["representatives"]) - excluded
    if len(candidates) != 63925:
        raise ArithmeticError(
            f"expected 63925 hidden minimum-norm-eight cosets, got {len(candidates)}"
        )

    rows = []
    l1_histogram = Counter()
    multiplicity_histogram = Counter()
    seen_section_masks = set()
    for short_mask in candidates:
        short_vector = vector(ZZ, streaming["representatives"][short_mask])
        section_vector = short_vector * short_change
        if section_vector * gram * section_vector != 8:
            raise ArithmeticError("selected representative has wrong equation height")
        section_mask = parity_mask(section_vector)
        if section_mask in seen_section_masks:
            raise ArithmeticError("duplicate section-basis parity mask")
        seen_section_masks.add(section_mask)
        row_score = score(section_vector)
        minimum_count = int(streaming["unoriented_multiplicities"][short_mask])
        l1_histogram[row_score[0]] += 1
        multiplicity_histogram[minimum_count] += 1
        rows.append(
            {
                "orbit_mask": section_mask,
                "orbit_hex": f"0x{section_mask:05x}",
                "short_orbit_mask": int(short_mask),
                "coefficient_l1": row_score[0],
                "support_count": row_score[1],
                "maximum_absolute_coefficient": row_score[2],
                "minimal_unoriented_count": minimum_count,
                "section_basis_w": row_score[3],
                "short_basis_w": tuple(map(int, short_vector)),
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
        "short_orbit_mask",
        "coefficient_l1",
        "support_count",
        "maximum_absolute_coefficient",
        "minimal_unoriented_count",
        "section_basis_w",
        "short_basis_w",
    ]
    lines = ["\t".join(fields)]
    for row in rows:
        serialized = dict(row)
        for key in ("section_basis_w", "short_basis_w"):
            serialized[key] = entries(serialized[key])
        lines.append("\t".join(str(serialized[field]) for field in fields))
    table_text = "\n".join(lines) + "\n"

    payload = {
        "schema": "elkies-k3.r17-norm12-103b2-norm8-pencil-priority.v1",
        "status": "PASS_EXACT_COMPLETE_103B2_NORM8_PENCIL_PRIORITY",
        "class_count": len(rows),
        "lattice_dictionary": {
            "class": "D_w=(2,2,w)",
            "conditions": ["D_w^2=0", "D_w.F=2", "D_w.O=0", "w.M.w=8"],
            "section_translation": "w -> w+2x",
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
                "coefficient_l1",
                "support_count",
                "maximum_absolute_coefficient",
                "lexicographic_oriented_section_basis_word",
            ],
            "l1_histogram": {
                str(key): value for key, value in sorted(l1_histogram.items())
            },
        },
        "minimal_unoriented_count_histogram": {
            str(key): value for key, value in sorted(multiplicity_histogram.items())
        },
        "coordinate_chain": "LLL-reduced section basis -> saturated hidden equation basis",
        "priority_table": relative(arguments.table_output),
        "priority_table_sha256": hashlib.sha256(table_text.encode()).hexdigest(),
        "inputs": {
            relative(path): digest(path) for path in (ENUMERATOR, model_path)
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
            "rank_r17_norm12_103b2_norm8_pencils.sage "
            f"--model {relative(model_path)} --output {relative(arguments.output)} "
            f"--table-output {relative(arguments.table_output)}"
        ),
        "proof_boundary": (
            "This is the complete hidden-103b2 minimum-norm-eight translation "
            "layer and its exact minimum-vector multiplicities. It constructs "
            "neither chord pencils nor branch characters."
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not arguments.table_output.exists() or arguments.table_output.read_text() != table_text:
            raise ArithmeticError("stored hidden norm-eight table differs from replay")
        if not arguments.output.exists() or arguments.output.read_text() != serialized:
            raise ArithmeticError("stored hidden norm-eight certificate differs from replay")
    else:
        arguments.table_output.parent.mkdir(parents=True, exist_ok=True)
        arguments.table_output.write_text(table_text)
        arguments.output.write_text(serialized)
    print(
        "HIDDEN103B2NORM8|classes={}|min_l1={}|max_l1={}|output={}".format(
            len(rows), rows[0]["coefficient_l1"], rows[-1]["coefficient_l1"],
            relative(arguments.output),
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
