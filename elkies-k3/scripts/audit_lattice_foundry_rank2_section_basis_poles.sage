#!/usr/bin/env sage-python
"""Audit the cheapest complete section basis for primitive MW-rank-two sources.

For a primitive rank-15 root lattice ``R`` in a positive-definite rank-17
frame ``M``, a root-adapted basis writes the Gram matrix as

    [ R  C ]
    [ C' D ].

A tail vector ``z in ZZ^2`` specifies a Mordell--Weil class modulo the root
lattice.  Completing the square reduces the exact minimum frame norm in that
class to an affine closest-vector problem in ``R``; its section pole order is
``(norm-4)/2``.  The Schur complement ``H = D-C'R^-1C`` is the exact MW height
Gram and gives the lower bound ``z'Hz <= norm``.

The displayed tail basis gives a finite initial upper bound.  We enumerate
every tail class satisfying the exact height bound, solve its affine CVP with
both double-double and MPFR-256 GSO arithmetic, and then test every pair with
tail determinant one.  This proves the lexicographically cheapest complete
basis (minimum maximum pole, then minimum other pole) within the numerical
CVP audit boundary.  Every returned norm and the height/Schur identity are
recomputed exactly over ``QQ``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from math import isqrt
from pathlib import Path

from fpylll import Enumeration, FPLLL, GSO, IntegerMatrix
from sage.all import QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-prescribed-root-sources-all-ns-3e8-all-a-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-rank2-section-basis-poles-v1.json"
)
ROOT_RANK = 15
TAIL_RANK = 2
FRAME_RANK = ROOT_RANK + TAIL_RANK


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
    inverse_root_times_c,
    tail: tuple[int, int],
    gso,
    mu,
) -> dict:
    tail_vector = vector(QQ, tail)
    target_original = -(inverse_root_times_c * tail_vector)
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
    norm = exact_norm(root_coordinates + list(tail), frame_gram)
    displacement = vector(QQ, root_coordinates) - target_original
    exact_distance = displacement * root_gram * displacement
    return {
        "norm": norm,
        "root_coordinates": root_coordinates,
        "reported_distance": float(reported_distance),
        "exact_distance": exact_distance,
        "distance_error": abs(float(reported_distance) - float(exact_distance)),
    }


def floor_sqrt_rational(value) -> int:
    value = QQ(value)
    if value < 0:
        raise ValueError("cannot take a real square root of a negative bound")
    return isqrt(int(value.numerator() // value.denominator()))


def canonical_tail(x: int, y: int) -> bool:
    """Retain one representative of each sign pair {z,-z}."""

    return x > 0 or (x == 0 and y > 0)


def determinant(left: tuple[int, int], right: tuple[int, int]) -> int:
    return left[0] * right[1] - left[1] * right[0]


def audit_source(entry: dict) -> dict:
    source = entry["source"]
    if source["mw_rank_for_rho_19"] != TAIL_RANK:
        raise ValueError("rank-two audit received a source of another MW rank")
    if source["root_rank"] != ROOT_RANK:
        raise ValueError("rank-two audit expects root rank 15")
    if not source["root_lattice_primitive"] or source["torsion"] != 1:
        raise ValueError("rank-two audit requires primitive torsion-free root quotient")

    frame_gram = matrix(ZZ, source["root_adapted_gram"])
    if frame_gram.nrows() != FRAME_RANK:
        raise ValueError("unexpected root-adapted frame rank")
    root_gram = matrix(QQ, frame_gram[:ROOT_RANK, :ROOT_RANK])
    c = matrix(
        QQ,
        ROOT_RANK,
        TAIL_RANK,
        [
            frame_gram[i, ROOT_RANK + j]
            for i in range(ROOT_RANK)
            for j in range(TAIL_RANK)
        ],
    )
    d = matrix(QQ, frame_gram[ROOT_RANK:, ROOT_RANK:])
    inverse_root_times_c = root_gram.inverse() * c
    schur_height = d - c.transpose() * inverse_root_times_c
    stored_height = matrix(QQ, source["mw_height_gram"])
    if schur_height != stored_height:
        raise ValueError("rank-two height does not equal the exact Schur complement")
    if not schur_height.is_positive_definite():
        raise ValueError("rank-two MW height is not positive definite")

    primary_gso, primary_mu = make_gso(root_gram, "dd", 0)
    audit_gso, audit_mu = make_gso(root_gram, "mpfr", 256)
    cache = {}
    maximum_error = 0.0

    def solve_tail(tail: tuple[int, int]) -> dict:
        nonlocal maximum_error
        canonical = tail
        sign = 1
        if not canonical_tail(*canonical):
            canonical = (-canonical[0], -canonical[1])
            sign = -1
        if canonical in cache:
            result = cache[canonical]
            if sign == 1:
                return result
            return {
                **result,
                "tail": list(tail),
                "root_coordinates": [-value for value in result["root_coordinates"]],
            }
        primary = affine_cvp(
            frame_gram,
            root_gram,
            inverse_root_times_c,
            canonical,
            primary_gso,
            primary_mu,
        )
        audit = affine_cvp(
            frame_gram,
            root_gram,
            inverse_root_times_c,
            canonical,
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
        tail_vector = vector(QQ, canonical)
        height = tail_vector * schur_height * tail_vector
        result = {
            "tail": list(canonical),
            "height": str(height),
            "minimum_frame_norm": norm,
            "section_pole_order": (norm - 4) // 2,
            "root_coordinates": primary["root_coordinates"],
        }
        cache[canonical] = result
        return result

    displayed = [solve_tail((1, 0)), solve_tail((0, 1))]
    initial_norm_bound = max(row["minimum_frame_norm"] for row in displayed)

    # Completing the square a second time gives exact independent bounds on
    # each tail coordinate among z with z'Hz <= initial_norm_bound.
    h00 = schur_height[0, 0]
    h01 = schur_height[0, 1]
    h11 = schur_height[1, 1]
    x_lower_coefficient = h00 - h01 * h01 / h11
    y_lower_coefficient = h11 - h01 * h01 / h00
    x_bound = floor_sqrt_rational(initial_norm_bound / x_lower_coefficient)
    y_bound = floor_sqrt_rational(initial_norm_bound / y_lower_coefficient)

    height_eligible = []
    for x in range(-x_bound, x_bound + 1):
        for y in range(-y_bound, y_bound + 1):
            if not canonical_tail(x, y):
                continue
            tail_vector = vector(QQ, [x, y])
            height = tail_vector * schur_height * tail_vector
            if height <= initial_norm_bound:
                height_eligible.append((x, y))
    height_eligible.sort(key=lambda tail: (abs(tail[0]) + abs(tail[1]), tail))

    norm_eligible = []
    for tail in height_eligible:
        row = solve_tail(tail)
        if row["minimum_frame_norm"] <= initial_norm_bound:
            norm_eligible.append(row)

    best = None
    best_key = None
    for left_index, left in enumerate(norm_eligible):
        left_tail = tuple(left["tail"])
        for right in norm_eligible[left_index + 1 :]:
            right_tail = tuple(right["tail"])
            if abs(determinant(left_tail, right_tail)) != 1:
                continue
            ordered = sorted(
                [left, right],
                key=lambda row: (
                    row["minimum_frame_norm"],
                    row["tail"],
                ),
            )
            key = (
                ordered[1]["minimum_frame_norm"],
                ordered[0]["minimum_frame_norm"],
                tuple(ordered[0]["tail"]),
                tuple(ordered[1]["tail"]),
            )
            if best_key is None or key < best_key:
                best_key = key
                best = ordered
    if best is None:
        raise RuntimeError("displayed unimodular tail basis was not recovered")

    basis_maximum_norm = max(row["minimum_frame_norm"] for row in best)
    # Any omitted tail class has height, hence frame norm, above the initial
    # bound.  Since the displayed basis lies within that bound, the exhaustive
    # determinant-one comparison proves the optimum maximum norm.
    if basis_maximum_norm > initial_norm_bound:
        raise AssertionError("optimized basis exceeded its initial upper bound")

    minimum = min(
        norm_eligible,
        key=lambda row: (row["minimum_frame_norm"], row["tail"]),
    )
    return {
        "ns_id": entry["ns_id"],
        "source_id": entry["source_id"],
        "source_gram_sha256": source["gram_sha256"],
        "root_type": source["root_type"],
        "status": "PASS_EXACT_NORMS_CROSS_PRECISION_COMPLETE_TAIL_ENUMERATION",
        "mw_height_gram": [[str(value) for value in row] for row in schur_height.rows()],
        "minimum_nonzero_section": minimum,
        "minimum_basis_maximum_frame_norm": basis_maximum_norm,
        "minimum_basis_maximum_pole_order": (basis_maximum_norm - 4) // 2,
        "minimum_basis_sorted_pole_profile": [
            row["section_pole_order"] for row in best
        ],
        "minimum_basis": best,
        "displayed_tail_basis": displayed,
        "initial_frame_norm_bound": initial_norm_bound,
        "height_coordinate_box": {
            "x_absolute_bound": x_bound,
            "y_absolute_bound": y_bound,
            "x_lower_coefficient": str(x_lower_coefficient),
            "y_lower_coefficient": str(y_lower_coefficient),
        },
        "height_eligible_tail_classes_modulo_sign": len(height_eligible),
        "norm_eligible_tail_classes_modulo_sign": len(norm_eligible),
        "eligible_tail_classes": sorted(
            norm_eligible,
            key=lambda row: (
                row["minimum_frame_norm"],
                row["tail"],
            ),
        ),
        "maximum_reported_distance_error": maximum_error,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    input_path = arguments.input.resolve()
    payload = json.loads(input_path.read_text())
    if payload.get("schema") != (
        "elkies-k3.lattice-foundry-prescribed-root-sources.v1"
    ):
        raise ValueError("unexpected prescribed-root source schema")
    entries = [
        entry
        for entry in payload["sources"]
        if entry["source"]["mw_rank_for_rho_19"] == TAIL_RANK
        and entry["source"]["root_lattice_primitive"]
        and entry["source"]["torsion"] == 1
    ]
    entries.sort(key=lambda entry: (entry["ns_id"], entry["source_id"]))
    if arguments.limit is not None:
        if arguments.check:
            raise SystemExit("--limit cannot be combined with --check")
        entries = entries[: arguments.limit]

    rows = [audit_source(entry) for entry in entries]
    histogram = {}
    for row in rows:
        key = str(row["minimum_basis_maximum_pole_order"])
        histogram[key] = histogram.get(key, 0) + 1
    output = {
        "schema": "elkies-k3.lattice-foundry-rank2-section-basis-poles.v1",
        "status": "PASS_PRIMITIVE_MW2_ROWS_EXHAUSTIVE_TAIL_CLASSES",
        "inputs": {relative(input_path): digest(input_path)},
        "method": {
            "primary_cvp": "fplll double-double affine enumeration",
            "audit_cvp": "independent fplll MPFR-256 affine enumeration",
            "tail_enumeration": (
                "all z in ZZ^2 modulo sign with exact z^T H z no larger "
                "than the displayed-basis frame-norm upper bound"
            ),
            "basis_test": "all eligible tail pairs with absolute determinant one",
            "basis_objective": (
                "minimum maximum section pole, then minimum other section pole"
            ),
        },
        "accounting": {
            "source_rows": len(rows),
            "minimum_basis_maximum_pole_histogram": dict(
                sorted(histogram.items(), key=lambda item: int(item[0]))
            ),
        },
        "proof_boundary": {
            "proved": (
                "The exact height bound makes the tail search finite; every tail "
                "class capable of improving the displayed basis is included, and "
                "every determinant-one pair is compared using exact frame norms."
            ),
            "numerical": (
                "CVP branch decisions are cross-precision audited, not formally "
                "verified; every returned norm and Schur height is exact."
            ),
            "scope": (
                "Primitive torsion-free MW-rank-two rows in the declared all-A "
                "prescribed-root inventory only."
            ),
        },
        "sources": rows,
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/audit_lattice_foundry_rank2_section_basis_poles.sage"
        ),
    }
    serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if output_path.read_text() != serialized:
            raise SystemExit("rank-two section-basis-pole artifact is stale")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    minimum = min(row["minimum_basis_maximum_pole_order"] for row in rows)
    print(
        "FOUNDRYRANK2BASISPOLES|"
        f"sources={len(rows)}|minimum={minimum}|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
