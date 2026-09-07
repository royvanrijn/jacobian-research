#!/usr/bin/env python3
"""Build the quotient-first fingerprint of Nagao's rank-20 section-7 fibre."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from math import comb
from pathlib import Path
import sys
from typing import Sequence


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
CAS = ELLIPTIC / "cas"
sys.path[:0] = [str(ELLIPTIC), str(CAS)]

from certify_nagao_rank20_t5081 import CONSTRUCTION, PARAMETER_T  # noqa: E402
from nagao_1994 import (  # noqa: E402
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    short_jacobian_coefficients,
)
from nagao_1994_section7 import SECTION7_LINEAR_COMPANION_SECTIONS  # noqa: E402
from latent_lattice import EllipticCurve, height_gram  # noqa: E402
from latent_lattice.pari import (  # noqa: E402
    recover_exact_embedding,
    row_embedding_smith_invariant_factors,
)


R17_BUILDER = CAS / "build_elkies_2026_rank_jump_fingerprints.py"
SPEC = importlib.util.spec_from_file_location("r17_rank_jump_fingerprints", R17_BUILDER)
R17 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(R17)

CERTIFICATE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elliptic_nagao_rank20_t5081_rank20_certificate.json"
)
GENERIC_RANK = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elliptic_nagao_section7_picard_bound.json"
)
LINEAR_SECTIONS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elliptic_nagao_section7_linear_sections.json"
)
DIRECTION = (
    ROOT
    / "archive/elliptic-curves/artifacts/generated-results"
    / "elliptic_nagao_rank20_t5081_direction.json"
)
COVER_SKEW = (
    ROOT
    / "archive/elliptic-curves/artifacts/generated-results"
    / "elliptic_nagao_rank20_t5081_cover_skew.json"
)
ELL2COVER = (
    ROOT
    / "archive/elliptic-curves/artifacts/generated-results"
    / "elliptic_nagao_rank20_t5081_ell2cover.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "nagao_section7_rank_jump_fingerprint_v1.json"
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def display(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def rational_rank(rows: Sequence[Sequence[int]]) -> int:
    if not rows:
        return 0
    matrix = [[Fraction(value) for value in row] for row in rows]
    rank = 0
    for column in range(len(matrix[0])):
        pivot = next(
            (row for row in range(rank, len(matrix)) if matrix[row][column]), None
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        scale = matrix[rank][column]
        matrix[rank] = [value / scale for value in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank or not matrix[row][column]:
                continue
            scale = matrix[row][column]
            matrix[row] = [
                left - scale * right
                for left, right in zip(matrix[row], matrix[rank])
            ]
        rank += 1
    return rank


def binary_rref(rows: Sequence[Sequence[int]]) -> tuple[list[list[int]], list[int]]:
    reduced: list[list[int]] = []
    pivots: list[int] = []
    for source in rows:
        row = [int(value) % 2 for value in source]
        for pivot, basis in zip(pivots, reduced):
            if row[pivot]:
                row = [left ^ right for left, right in zip(row, basis)]
        if not any(row):
            continue
        pivot = row.index(1)
        for index, basis in enumerate(reduced):
            if basis[pivot]:
                reduced[index] = [
                    left ^ right for left, right in zip(basis, row)
                ]
        insertion = sum(existing < pivot for existing in pivots)
        pivots.insert(insertion, pivot)
        reduced.insert(insertion, row)
    return reduced, pivots


def binary_coset_reducer(rows: Sequence[Sequence[int]]):
    reduced, pivots = binary_rref(rows)

    def reduce(vector: Sequence[int]) -> tuple[int, ...]:
        answer = [int(value) % 2 for value in vector]
        for pivot, basis in zip(pivots, reduced):
            if answer[pivot]:
                answer = [left ^ right for left, right in zip(answer, basis)]
        return tuple(answer)

    return reduce, len(reduced)


def degree_visibility(
    generic_rows: Sequence[Sequence[int]], quotient_tensor_dimension: int
) -> list[dict[str, object]]:
    direction = json.loads(DIRECTION.read_text())
    skew = json.loads(COVER_SKEW.read_text())
    ell2 = json.loads(ELL2COVER.read_text())
    direct_relations = tuple(
        record["basis_relation"] for record in direction["results"]["candidate_points"]
    )
    skew_relations = tuple(
        record["basis_relation"] for record in skew["results"]["candidate_points"]
    )
    relations = direct_relations + skew_relations
    if len(relations) != len(set(map(tuple, relations))) or len(relations) != 224:
        raise AssertionError("the bounded degree-two relation census changed")
    generic_rank = rational_rank(generic_rows)
    augmented_rank = rational_rank(tuple(generic_rows) + relations)
    reduce, generic_mod2_rank = binary_coset_reducer(generic_rows)
    direct_classes = tuple(reduce(relation) for relation in direct_relations)
    skew_classes = tuple(reduce(relation) for relation in skew_relations)
    classes = Counter(direct_classes + skew_classes)
    zero = (0,) * len(generic_rows[0])
    distinct_classes = tuple(classes)
    tensor_span = R17.rank_mod_prime(distinct_classes, 2)
    if tensor_span != quotient_tensor_dimension:
        raise AssertionError("the bounded degree-two classes stopped spanning mod 2")
    return [
        {
            "cover_degree": 2,
            "atlas_status": "bounded_archived_cover_search_exact_relation_replay",
            "ambient_candidate_point_count": len(relations),
            "distinct_exact_ambient_relations": len(set(map(tuple, relations))),
            "visible_free_quotient_span_dimension": augmented_rank - generic_rank,
            "visible_tensor_quotient_span_dimension_over_f2": tensor_span,
            "generic_image_dimension_in_ambient_mod2": generic_mod2_rank,
            "distinct_quotient_classes_over_f2": len(classes),
            "zero_quotient_class_point_count": classes.get(zero, 0),
            "nonzero_quotient_class_count": len(classes) - (zero in classes),
            "class_multiplicity_histogram": {
                str(multiplicity): count
                for multiplicity, count in sorted(Counter(classes.values()).items())
            },
            "same_class_candidate_pair_count": sum(
                comb(multiplicity, 2) for multiplicity in classes.values()
            ),
            "class_support_weight_histogram_in_canonical_ambient_representatives": {
                str(weight): count
                for weight, count in sorted(
                    Counter(sum(vector) for vector in classes).items()
                )
            },
            "cover_packets": [
                {
                    "id": "direct-and-alternate-quartic-covers",
                    "candidate_count": len(direct_relations),
                    "distinct_quotient_classes_over_f2": len(set(direct_classes)),
                    "selected_cover_count": direction["full_mod2_class_scan"][
                        "selected_cover_count_after_weight_diversification"
                    ],
                    "ambient_mod2_classes_scored": direction["full_mod2_class_scan"][
                        "nonzero_classes_scored"
                    ],
                },
                {
                    "id": "skew-quartic-cover",
                    "candidate_count": len(skew_relations),
                    "distinct_quotient_classes_over_f2": len(set(skew_classes)),
                },
            ],
            "packet_class_intersection_count": len(
                set(direct_classes) & set(skew_classes)
            ),
            "exact_relation_replay": True,
            "scope": (
                "complete for the archived bounded searches and their returned "
                "points, not a complete atlas of every degree-two cover"
            ),
        },
        {
            "cover_degree": 3,
            "atlas_status": "not_available",
            "visible_free_quotient_span_dimension": None,
            "reason": "no complete descended degree-three cover atlas; missing is not zero",
        },
        {
            "cover_degree": 4,
            "atlas_status": "strict_timeout_no_cover_output",
            "visible_free_quotient_span_dimension": None,
            "source_status": ell2["status"],
            "reason": "the archived ell2cover attempt timed out; missing is not zero",
        },
    ]


def build_payload(
    *, digits: int, maximum_vectors: int, shortest_vector_count: int, timeout: float
) -> dict[str, object]:
    certificate = json.loads(CERTIFICATE.read_text())
    ambient_basis = tuple(
        (Fraction(record["jacobian_x"]), Fraction(record["jacobian_y"]))
        for record in certificate["exact_rank_certificate"]["saturated_basis"]
    )
    visible = tuple(
        quartic_point_to_short_jacobian(CONSTRUCTION, PARAMETER_T, point)
        for point in primitive_visible_points(CONSTRUCTION, PARAMETER_T)
    )
    generic_basis = visible[:11] + (
        SECTION7_LINEAR_COMPANION_SECTIONS[0].jacobian_point(PARAMETER_T),
    )
    model = short_jacobian_coefficients(CONSTRUCTION, PARAMETER_T)
    curve = EllipticCurve(model)
    embedding = recover_exact_embedding(
        curve,
        ambient_basis,
        generic_basis,
        digits=max(digits, 100),
        timeout=timeout,
    )
    factors = tuple(sorted(row_embedding_smith_invariant_factors(embedding.columns)))
    quotient = R17.smith_quotient_record(
        ambient_rank=20, subgroup_rank=12, smith_factors=factors
    )
    gram = height_gram(curve, ambient_basis, digits=digits, timeout=timeout)
    geometry = R17.quotient_height_geometry(
        gram,
        embedding.columns,
        digits=digits,
        maximum_vectors=maximum_vectors,
        shortest_vector_count=shortest_vector_count,
        timeout=timeout,
    )
    visibility = degree_visibility(
        embedding.columns,
        quotient["tensor_dimensions_over_f_ell"]["2"],
    )
    fingerprint = {
        "label": "nagao-section7-rank20-t5081",
        "constructor_parameter_T": "5081/47",
        "paper_parameter_t": "5081/94",
        "global_or_canonical_model": [str(value) for value in model],
        "certified_rank_lower_bound": 20,
        "generic_rank": 12,
        "generic_basis": "visible sections 0,...,10 followed by plus-7/27",
        "exact_generic_embedding": {
            "columns": [list(column) for column in embedding.columns],
            "smith_invariant_factors": list(factors),
            "maximum_absolute_coordinate": embedding.max_abs_coordinate,
            "nonzero_coordinate_count": embedding.nonzero_coordinates,
            "height_dual_numerical_residual_max": embedding.numerical_residual_max,
            "exact_group_law_replay": True,
        },
        "quotient_structure": quotient,
        "quotient_height_geometry": geometry,
        "degree_visibility": visibility,
        "response_variables": {
            "free_quotient_rank_lower_bound": quotient[
                "free_quotient_rank_lower_bound"
            ],
            "tensor_dimensions_over_f_ell": quotient[
                "tensor_dimensions_over_f_ell"
            ],
            "degree_two_bounded_visible_free_span": visibility[0][
                "visible_free_quotient_span_dimension"
            ],
            "degree_two_bounded_visible_f2_tensor_span": visibility[0][
                "visible_tensor_quotient_span_dimension_over_f2"
            ],
        },
    }
    return {
        "schema": "elliptic-curves.nagao-section7-rank-jump-fingerprint.v1",
        "status": "PASS_CERTIFIED_SUBGROUP_QUOTIENT_FINGERPRINT",
        "target": "escape from the generic rank-12 Mordell-Weil lattice",
        "fingerprints": [fingerprint],
        "height_precision_decimal_digits": digits,
        "inputs": {
            display(CERTIFICATE): digest(CERTIFICATE),
            display(GENERIC_RANK): digest(GENERIC_RANK),
            display(LINEAR_SECTIONS): digest(LINEAR_SECTIONS),
            display(DIRECTION): digest(DIRECTION),
            display(COVER_SKEW): digest(COVER_SKEW),
            display(ELL2COVER): digest(ELL2COVER),
            display(R17_BUILDER): digest(R17_BUILDER),
        },
        "proof_boundary": (
            "generic rank 12 and specialized rank at least 20 are exact; the "
            "embedding, Smith structure, and returned cover-point relations replay "
            "exactly. Canonical heights are high-precision numerical values. The "
            "ambient rank-20 subgroup is certified but is not asserted to be all E(Q)."
        ),
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves:elliptic-curves/cas python3 "
            "elliptic-curves/cas/build_nagao_section7_rank_jump_fingerprint.py"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--digits", type=int, default=80)
    parser.add_argument("--maximum-vectors", type=int, default=2_000_000)
    parser.add_argument("--shortest-vector-count", type=int, default=16)
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.digits < 50:
        raise SystemExit("--digits must be at least 50")
    payload = build_payload(
        digits=args.digits,
        maximum_vectors=args.maximum_vectors,
        shortest_vector_count=args.shortest_vector_count,
        timeout=args.timeout,
    )
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != rendered:
            raise SystemExit(f"stale or missing fingerprint artifact: {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    row = payload["fingerprints"][0]
    print(
        "NAGAOSECTION7RANKJUMPFINGERPRINT|free="
        f"{row['quotient_structure']['free_quotient_rank_lower_bound']}|mod2="
        f"{row['quotient_structure']['tensor_dimensions_over_f_ell']['2']}|degree2="
        f"{row['degree_visibility'][0]['visible_free_quotient_span_dimension']}|status="
        f"{payload['status']}"
    )


if __name__ == "__main__":
    main()
