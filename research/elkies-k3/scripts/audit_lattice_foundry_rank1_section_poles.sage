#!/usr/bin/env sage-python
"""Audit minimum section pole orders for prescribed-root MW-rank-one sources.

For a primitive rank-16 root lattice ``R`` in a positive-definite rank-17
frame ``M``, a root-adapted basis writes the Gram matrix as

    [ R  b ]
    [ b' c ].

Vectors with last coordinate ``n`` represent the nonzero Mordell--Weil class
``nP`` modulo the root lattice.  Completing the square reduces the exact
minimum frame norm in that class to an affine closest-vector problem in
``R``.  The associated section has pole order ``(norm-4)/2``.  Once
``(n+1)^2 h`` exceeds the best norm found, where ``h`` is the rank-one height,
the positive Schur complement proves that no later multiple can improve it.

The CVP branch decisions are independently repeated with double-double and
256-bit MPFR GSO arithmetic.  Every returned norm and the height/Schur
identity are recomputed exactly over ``QQ``.  Nonprimitive root rows are
retained with a typed open status because their torsion glue must first be
resolved; no height or pole is inferred for them.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from fpylll import Enumeration, FPLLL, GSO, IntegerMatrix
from sage.all import QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUTS = [
    ROOT
    / (
        "artifacts/generated-results/"
        "elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-"
        f"group-{group}-v1.json"
    )
    for group in "abcd"
]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-rank1-section-poles-v1.json"
)
ROOT_RANK = 16
FRAME_RANK = 17


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_norm(values: list[int], gram) -> int:
    value = vector(ZZ, values)
    return int(value * gram * value)


def make_gso(root_gram, float_type: str, precision: int):
    if float_type == "mpfr":
        FPLLL.set_precision(precision)
    integer_matrix = IntegerMatrix.from_matrix(
        [[int(entry) for entry in row] for row in root_gram.rows()]
    )
    gso = GSO.Mat(
        integer_matrix, gram=True, float_type=float_type, update=True
    )
    mu = [
        [gso.get_mu(i, j) if i > j else 0.0 for j in range(ROOT_RANK)]
        for i in range(ROOT_RANK)
    ]
    return gso, mu


def affine_cvp(
    frame_gram,
    root_gram,
    inverse_root_times_b,
    multiple: int,
    gso,
    mu,
) -> dict:
    target_original = -multiple * inverse_root_times_b
    target_gso = [
        float(target_original[i])
        + sum(
            float(target_original[j]) * mu[j][i]
            for j in range(i + 1, ROOT_RANK)
        )
        for i in range(ROOT_RANK)
    ]
    zero_distance = target_original * root_gram * target_original
    solutions = Enumeration(gso).enumerate(
        0,
        ROOT_RANK,
        float(zero_distance) + 1.0,
        0,
        target=target_gso,
    )
    if not solutions:
        raise RuntimeError("affine root-lattice CVP returned no solution")
    reported_distance, coordinates = solutions[0]
    root_coordinates = [int(round(value)) for value in coordinates]
    if any(
        abs(value - integer) > 1e-7
        for value, integer in zip(coordinates, root_coordinates)
    ):
        raise RuntimeError("affine CVP coordinates are not integral")
    full_coordinates = root_coordinates + [multiple]
    norm = exact_norm(full_coordinates, frame_gram)
    exact_distance = (
        (vector(QQ, root_coordinates) - target_original)
        * root_gram
        * (vector(QQ, root_coordinates) - target_original)
    )
    return {
        "norm": norm,
        "root_coordinates": root_coordinates,
        "reported_distance": float(reported_distance),
        "exact_distance": exact_distance,
        "distance_error": abs(float(reported_distance) - float(exact_distance)),
    }


def audit_source(entry: dict) -> dict:
    source = entry["source"]
    base = {
        "ns_id": entry["ns_id"],
        "source_id": entry["source_id"],
        "source_gram_sha256": source["gram_sha256"],
        "root_type": source["root_type"],
        "root_lattice_primitive": bool(source["root_lattice_primitive"]),
        "torsion": source["torsion"],
    }
    if not source["root_lattice_primitive"]:
        return {
            **base,
            "status": "OPEN_NONPRIMITIVE_ROOT_GLUE_ANALYSIS_REQUIRED",
            "minimum_section_pole_order": None,
        }
    if source["mw_rank_for_rho_19"] != 1 or source["root_rank"] != ROOT_RANK:
        raise ValueError("section-pole audit expects root rank 16 and MW rank 1")
    if source["torsion"] != 1 or not source["mw_height_gram"]:
        raise ValueError("primitive rank-one source lacks its exact height")

    frame_gram = matrix(ZZ, source["root_adapted_gram"])
    if frame_gram.nrows() != FRAME_RANK:
        raise ValueError("unexpected root-adapted frame rank")
    root_gram = matrix(QQ, frame_gram[:ROOT_RANK, :ROOT_RANK])
    b = vector(QQ, [frame_gram[i, ROOT_RANK] for i in range(ROOT_RANK)])
    c = QQ(frame_gram[ROOT_RANK, ROOT_RANK])
    inverse_root_times_b = root_gram.inverse() * b
    height = QQ(source["mw_height_gram"][0][0])
    schur_height = c - b * inverse_root_times_b
    if schur_height != height:
        raise ValueError("rank-one height does not equal the exact Schur complement")

    primary_gso, primary_mu = make_gso(root_gram, "dd", 0)
    audit_gso, audit_mu = make_gso(root_gram, "mpfr", 256)
    tested = []
    best = None
    multiple = 1
    maximum_error = 0.0
    while True:
        primary = affine_cvp(
            frame_gram,
            root_gram,
            inverse_root_times_b,
            multiple,
            primary_gso,
            primary_mu,
        )
        audit = affine_cvp(
            frame_gram,
            root_gram,
            inverse_root_times_b,
            multiple,
            audit_gso,
            audit_mu,
        )
        if primary["norm"] != audit["norm"]:
            raise ValueError("cross-precision affine CVP norm mismatch")
        norm = primary["norm"]
        if norm < 4 or norm % 2:
            raise ValueError(f"invalid nonzero-section frame norm {norm}")
        maximum_error = max(
            maximum_error,
            primary["distance_error"],
            audit["distance_error"],
        )
        tested.append(
            {
                "multiple": multiple,
                "minimum_frame_norm": norm,
                "section_pole_order": (norm - 4) // 2,
                "root_coordinates": primary["root_coordinates"],
            }
        )
        candidate = (norm, multiple, primary["root_coordinates"])
        if best is None or candidate[:2] < best[:2]:
            best = candidate
        next_multiple = multiple + 1
        if next_multiple * next_multiple * height > best[0]:
            break
        multiple = next_multiple

    assert best is not None
    return {
        **base,
        "status": "PASS_EXACT_NORMS_CROSS_PRECISION_AFFINE_CVP",
        "mw_height": str(height),
        "minimum_section_pole_order": (best[0] - 4) // 2,
        "minimum_section_frame_norm": best[0],
        "minimizing_multiple": best[1],
        "root_coordinates": best[2],
        "multiples_tested": tested,
        "next_multiple_height_lower_bound": str(
            (multiple + 1) * (multiple + 1) * height
        ),
        "maximum_reported_distance_error": maximum_error,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, action="append", default=[])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    input_paths = [
        path.resolve() for path in (arguments.input or DEFAULT_INPUTS)
    ]
    rows = []
    inputs = {}
    for path in input_paths:
        payload = json.loads(path.read_text())
        if payload.get("schema") != (
            "elkies-k3.lattice-foundry-prescribed-root-sources.v1"
        ):
            raise ValueError(f"unexpected source schema in {path}")
        inputs[relative(path)] = digest(path)
        for entry in payload["sources"]:
            row = audit_source(entry)
            row["source_artifact"] = relative(path)
            rows.append(row)
    rows.sort(key=lambda row: (row["ns_id"], row["source_artifact"], row["source_id"]))

    passed = [row for row in rows if row["status"].startswith("PASS_")]
    open_rows = [row for row in rows if row["status"].startswith("OPEN_")]
    pole_histogram = {}
    for row in passed:
        key = str(row["minimum_section_pole_order"])
        pole_histogram[key] = pole_histogram.get(key, 0) + 1
    output = {
        "schema": "elkies-k3.lattice-foundry-rank1-section-poles.v1",
        "status": "PASS_PRIMITIVE_ROWS_WITH_TYPED_NONPRIMITIVE_OPEN_ROWS",
        "inputs": inputs,
        "method": {
            "primary_cvp": "fplll double-double affine enumeration",
            "audit_cvp": "independent fplll MPFR-256 affine enumeration",
            "exact_checks": (
                "integral frame norm, rational Schur height, and future-multiple "
                "height lower bound"
            ),
        },
        "accounting": {
            "source_rows": len(rows),
            "primitive_rows_passed": len(passed),
            "nonprimitive_rows_open": len(open_rows),
            "minimum_section_pole_histogram_primitive_rows": dict(
                sorted(pole_histogram.items(), key=lambda item: int(item[0]))
            ),
        },
        "proof_boundary": {
            "proved": (
                "For every primitive-root row, all MW multiples capable of "
                "improving the best exact frame norm are checked."
            ),
            "numerical": (
                "CVP branch decisions are cross-precision audited, not formally "
                "verified; every returned norm and stopping bound is exact."
            ),
            "open": (
                "Nonprimitive root rows require torsion glue analysis before a "
                "Mordell--Weil generator height or pole is assigned."
            ),
        },
        "sources": rows,
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/audit_lattice_foundry_rank1_section_poles.sage"
        ),
    }
    serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if output_path.read_text() != serialized:
            raise SystemExit("rank-one section-pole artifact is stale")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "FOUNDRYRANK1POLES|"
        f"sources={len(rows)}|passed={len(passed)}|open={len(open_rows)}|"
        f"minimum={min(row['minimum_section_pole_order'] for row in passed)}|"
        "status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
