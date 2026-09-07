#!/usr/bin/env python3
"""Build quotient-first fingerprints for the two certified Fermigier anchors."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from math import comb, factorial
from pathlib import Path
import sys
from typing import Sequence


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
CAS = ELLIPTIC / "cas"
sys.path[:0] = [str(ELLIPTIC), str(CAS)]

from ecsearch.fermigier_rank import specialize_fermigier_rank_sections  # noqa: E402
from elliptic_candidate_record import (  # noqa: E402
    WeierstrassChange,
    change_weierstrass_model,
    source_point_to_target,
)
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

E22_POINTS = ELLIPTIC / "data/fermigier_e22_points.json"
E22_CERTIFICATE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elliptic_fermigier_rank22_points.json"
)
RANK20_CERTIFICATE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json"
)
DIRECTION_BALL = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elliptic_fermigier_exceptional_quotient_ball.json"
)
GENERIC_RANK = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elliptic_fermigier_generic_rank_exact.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "fermigier_rank_jump_fingerprints_v1.json"
)
E22_MINIMAL_CHANGE = WeierstrassChange.from_values(
    (
        "7/1521",
        "-54787899485230240/771147",
        "-757/1521",
        "124998592673793420851/3518743761",
    )
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def display(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def point(record: dict[str, str]) -> tuple[Fraction, Fraction]:
    return Fraction(record["x"]), Fraction(record["y"])


def binary_rank(rows: Sequence[Sequence[int]]) -> int:
    return R17.rank_mod_prime(rows, 2)


def cycle_matroid_intersections(
    *, dimension: int, coefficient_vectors: Sequence[Sequence[int]]
) -> dict[str, object]:
    """Exact class/circuit census for signed weight-one/two directions.

    Modulo sign these are the edges of ``K_(dimension+1)``: basis vectors are
    edges to a distinguished root and pair sums are the remaining edges.
    Minimal dependencies are therefore precisely simple cycles.
    """

    classes = Counter(tuple(abs(int(value)) % 2 for value in row) for row in coefficient_vectors)
    if any(sum(vector) not in (1, 2) for vector in classes):
        raise AssertionError("the direction ball is no longer weight one/two")
    if binary_rank(tuple(classes)) != dimension:
        raise AssertionError("the bounded direction ball stopped spanning")
    expected_classes = dimension + comb(dimension, 2)
    if len(classes) != expected_classes:
        raise AssertionError("the direction-class census changed")
    circuits = {}
    for length in range(3, dimension + 2):
        class_cycles_with_root = comb(dimension, length - 1) * factorial(length - 1) // 2
        class_cycles_without_root = (
            comb(dimension, length) * factorial(length - 1) // 2
            if length <= dimension
            else 0
        )
        labelled_lifts_with_root = class_cycles_with_root * 4 ** (length - 1)
        labelled_lifts_without_root = class_cycles_without_root * 4**length
        circuits[str(length)] = {
            "distinct_class_circuits": class_cycles_with_root + class_cycles_without_root,
            "signed_label_circuits": labelled_lifts_with_root + labelled_lifts_without_root,
        }
    return {
        "matroid_identification": f"cycle matroid of complete graph K_{dimension + 1}",
        "visible_free_quotient_span_dimension": dimension,
        "signed_direction_count": len(coefficient_vectors),
        "distinct_unoriented_mod2_classes": len(classes),
        "class_multiplicity_histogram": {
            str(multiplicity): count
            for multiplicity, count in sorted(Counter(classes.values()).items())
        },
        "minimal_dependency_circuit_census_by_size": circuits,
        "intersection_interpretation": (
            "two direction packets overlap in quotient support according to the "
            "graphic matroid; simple cycles are exactly the minimal intersections"
        ),
    }


def visibility_record(anchor: dict[str, object], dimension: int) -> list[dict[str, object]]:
    records = anchor["records"]
    coefficients = tuple(record["coefficient_vector"] for record in records)
    if anchor["exceptional_basis_dimension"] != dimension:
        raise AssertionError("direction dimension and quotient rank differ")
    return [
        {
            "cover_degree": 2,
            "atlas_status": "complete_bounded_signed_weight_at_most_two_direction_ball",
            "basis_labels": anchor["exceptional_basis_labels"],
            **cycle_matroid_intersections(
                dimension=dimension, coefficient_vectors=coefficients
            ),
            "minimum_quartic_x_projective_height": anchor[
                "minimum_quartic_x_projective_height"
            ],
            "maximum_quartic_x_projective_height": anchor[
                "maximum_quartic_x_projective_height"
            ],
            "scope": (
                "complete for the declared coefficient alphabet and support bound; "
                "not a complete census of every degree-two cover"
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
            "atlas_status": "not_available",
            "visible_free_quotient_span_dimension": None,
            "reason": "no complete descended degree-four cover atlas; missing is not zero",
        },
    ]


def embedding_record(embedding) -> dict[str, object]:
    factors = tuple(sorted(row_embedding_smith_invariant_factors(embedding.columns)))
    return {
        "columns": [list(column) for column in embedding.columns],
        "smith_invariant_factors": list(factors),
        "maximum_absolute_coordinate": embedding.max_abs_coordinate,
        "nonzero_coordinate_count": embedding.nonzero_coordinates,
        "height_dual_numerical_residual_max": embedding.numerical_residual_max,
        "exact_group_law_replay": True,
    }


def build_one(
    *,
    label: str,
    parameter_u: Fraction,
    model: Sequence[Fraction],
    ambient_basis,
    generic_basis,
    certified_rank: int,
    direction_anchor: dict[str, object],
    digits: int,
    maximum_vectors: int,
    shortest_vector_count: int,
    timeout: float,
) -> dict[str, object]:
    curve = EllipticCurve(tuple(model))
    embedding = recover_exact_embedding(
        curve, ambient_basis, generic_basis, digits=max(digits, 100), timeout=timeout
    )
    exact_embedding = embedding_record(embedding)
    factors = exact_embedding["smith_invariant_factors"]
    quotient = R17.smith_quotient_record(
        ambient_rank=certified_rank,
        subgroup_rank=12,
        smith_factors=factors,
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
    visibility = visibility_record(
        direction_anchor, quotient["free_quotient_rank_lower_bound"]
    )
    return {
        "label": label,
        "canonical_parameter_u": str(parameter_u),
        "global_or_canonical_model": [str(value) for value in model],
        "certified_rank_lower_bound": certified_rank,
        "generic_rank": 12,
        "exact_generic_embedding": exact_embedding,
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
            "degree_two_bounded_ball_visible_free_span": visibility[0][
                "visible_free_quotient_span_dimension"
            ],
        },
    }


def build_payload(
    *, digits: int, maximum_vectors: int, shortest_vector_count: int, timeout: float
) -> dict[str, object]:
    e22_data = json.loads(E22_POINTS.read_text())
    e22_model = tuple(map(Fraction, e22_data["weierstrass_coefficients"]))
    e22_basis = tuple(
        (Fraction(x), Fraction(y)) for x, y in e22_data["points"]
    )
    e22_specialization = specialize_fermigier_rank_sections(Fraction(19754, 39))
    if change_weierstrass_model(
        e22_specialization.canonical_model, E22_MINIMAL_CHANGE
    ) != e22_model:
        raise AssertionError("the exact E22 minimal-model change changed")
    e22_generic = tuple(
        source_point_to_target(point_value, E22_MINIMAL_CHANGE)
        for point_value in e22_specialization.section_differences
    )

    rank20_data = json.loads(RANK20_CERTIFICATE.read_text())
    rank20_model = tuple(
        map(
            Fraction,
            rank20_data["models"]["canonical_generalized"][
                "coefficients_a1_a2_a3_a4_a6"
            ],
        )
    )
    rank20_basis = tuple(
        point(record["points"]["canonical_generalized"])
        for record in rank20_data["imported_selected_twenty_basis"]["basis"]
    )
    rank20_specialization = specialize_fermigier_rank_sections(Fraction(28917, 20))
    if rank20_specialization.canonical_model != rank20_model:
        raise AssertionError("the rank-20 canonical model changed")
    directions = json.loads(DIRECTION_BALL.read_text())["direction_balls"]
    fingerprints = [
        build_one(
            label="fermigier-E22",
            parameter_u=Fraction(19754, 39),
            model=e22_model,
            ambient_basis=e22_basis,
            generic_basis=e22_generic,
            certified_rank=22,
            direction_anchor=directions["E22"],
            digits=digits,
            maximum_vectors=maximum_vectors,
            shortest_vector_count=shortest_vector_count,
            timeout=timeout,
        ),
        build_one(
            label="fermigier-rank20-near-miss",
            parameter_u=Fraction(28917, 20),
            model=rank20_model,
            ambient_basis=rank20_basis,
            generic_basis=rank20_specialization.section_differences,
            certified_rank=20,
            direction_anchor=directions["rank20"],
            digits=digits,
            maximum_vectors=maximum_vectors,
            shortest_vector_count=shortest_vector_count,
            timeout=timeout,
        ),
    ]
    return {
        "schema": "elliptic-curves.fermigier-rank-jump-fingerprints.v1",
        "status": "PASS_CERTIFIED_SUBGROUP_QUOTIENT_FINGERPRINTS",
        "target": "escape from the generic rank-12 Mordell-Weil lattice",
        "fingerprints": fingerprints,
        "height_precision_decimal_digits": digits,
        "inputs": {
            display(E22_POINTS): digest(E22_POINTS),
            display(E22_CERTIFICATE): digest(E22_CERTIFICATE),
            display(RANK20_CERTIFICATE): digest(RANK20_CERTIFICATE),
            display(DIRECTION_BALL): digest(DIRECTION_BALL),
            display(GENERIC_RANK): digest(GENERIC_RANK),
            display(R17_BUILDER): digest(R17_BUILDER),
        },
        "proof_boundary": (
            "rank and generic rank are exact lower/theorem data; embeddings, Smith "
            "structure, and direction relations replay exactly. Canonical heights "
            "are high-precision numerical values. Quotients are relative to the "
            "displayed certified ambient subgroups, not asserted full E(Q)."
        ),
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves:elliptic-curves/cas python3 "
            "elliptic-curves/cas/build_fermigier_rank_jump_fingerprints.py"
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
    summary = ",".join(
        f"{row['label']}:{row['quotient_structure']['free_quotient_rank_lower_bound']}"
        for row in payload["fingerprints"]
    )
    print(f"FERMIGIERRANKJUMPFINGERPRINTS|quotients={summary}|status={payload['status']}")


if __name__ == "__main__":
    main()
