#!/usr/bin/env python3
"""Reverse-engineer construction fingerprints of the rank-28 R17 quotient.

This is deliberately a positive-control analysis.  It starts from the eleven
public complement points already certified independent modulo the specialized
generic R17 subgroup.  It does not search for a new point or infer a rank upper
bound.

The exact layers are the point coordinates, Kummer representatives, division
polynomials, their factorization/Galois data, Frobenius factorization patterns,
and the rational isogeny graph.  Canonical heights and the closest-vector shell
calculation use the existing high-precision PARI height Gram.  The latter is
replayed by complete Fincke--Pohst enumeration on three rounded Gram matrices;
agreement is a stability audit, not interval certification of canonical heights.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, getcontext
from fractions import Fraction
from hashlib import sha256
import importlib
import json
from math import gcd, isqrt, log, sqrt
from pathlib import Path
import sys
from typing import Any

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sage.all import (  # type: ignore
    EllipticCurve,
    GF,
    PolynomialRing,
    QQ,
    ZZ,
    matrix,
    pari,
    prime_range,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path[:0] = [str(ELLIPTIC), str(ELLIPTIC / "cas")]

RELATIONS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_exceptional_specialization_relations_v1.json"
)
TRUTH = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "latent_lattice_calibration_truth_v1.json"
)
SELMER = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_rank28_public11_selmer_candidates_v1.json"
)
LOCAL_SIGNATURES = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_rank28_generic17_local_signature_v1.json"
)
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
PUBLIC = ELLIPTIC / "cas/elkies_rank28.py"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_rank28_construction_fingerprints_v1.json"
)

PARAMETER = "-9529/5471"
ROUNDING_SCALES = (10_000, 100_000, 1_000_000)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def display(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def dtext(value: Decimal, digits: int = 34) -> str:
    if not value:
        return "0"
    return format(value, f".{digits}g")


def bit_height(value: Fraction | int) -> int:
    value = Fraction(str(value))
    return max(abs(value.numerator).bit_length(), value.denominator.bit_length())


def rational_is_square(value: Fraction) -> bool:
    if value < 0:
        return False
    numerator_root = isqrt(value.numerator)
    denominator_root = isqrt(value.denominator)
    return (
        numerator_root * numerator_root == value.numerator
        and denominator_root * denominator_root == value.denominator
    )


def polynomial_key(polynomial) -> list[str]:
    return [str(coefficient) for coefficient in polynomial.list()]


def polynomial_digest(polynomial) -> str:
    payload = json.dumps(polynomial_key(polynomial), separators=(",", ":"))
    return sha256(payload.encode()).hexdigest()


def pari_columns(gen_matrix, dimension: int) -> list[list[int]]:
    return [
        [int(gen_matrix[row, column]) for row in range(dimension)]
        for column in range(int(gen_matrix.ncols()))
    ]


def decimal_quadratic(gram: list[list[Decimal]], coordinates: list[int]) -> Decimal:
    return sum(
        Decimal(coordinates[left])
        * gram[left][right]
        * Decimal(coordinates[right])
        for left in range(len(coordinates))
        for right in range(len(coordinates))
    )


def rounded_matrix(gram: list[list[Decimal]], scale: int):
    return matrix(
        ZZ,
        [
            [int((value * Decimal(scale)).to_integral_value()) for value in row]
            for row in gram
        ],
    )


def first_shell_stability(
    generic_gram: list[list[Decimal]],
    cross: list[Decimal],
    raw: Decimal,
) -> dict[str, Any]:
    """Find the shortest representative of Q+M on rounded height forms."""

    dimension = len(generic_gram)
    corrections = []
    rounded_values = []
    generic_shell_vectors = []
    generic_shell_values = []
    for scale in ROUNDING_SCALES:
        generic_integer = rounded_matrix(generic_gram, scale)
        generic_minimum = pari(generic_integer).qfminim()
        generic_vectors = pari_columns(generic_minimum[2], dimension)
        generic_shell_vectors.append(generic_vectors[0])
        generic_shell_values.append(Decimal(str(generic_minimum[1])) / Decimal(scale))

        full_gram = [
            list(generic_gram[row]) + [cross[row]] for row in range(dimension)
        ] + [list(cross) + [raw]]
        full_integer = rounded_matrix(full_gram, scale)
        # Q itself is an a priori bound, and these rank-18 enumerations retain
        # only a few hundred signed vectors at the present control.
        bound = int((raw * Decimal(scale)).to_integral_value()) + 2
        enumeration = pari(full_integer).qfminim(bound, None, 2)
        candidates = []
        for candidate in pari_columns(enumeration[2], dimension + 1):
            if abs(candidate[-1]) != 1:
                continue
            if candidate[-1] < 0:
                candidate = [-value for value in candidate]
            value = int(vector(ZZ, candidate) * full_integer * vector(ZZ, candidate))
            candidates.append((value, candidate[:-1]))
        if not candidates:
            raise ArithmeticError("rounded coset enumeration found no Q+M representative")
        value, correction = min(candidates)
        corrections.append(correction)
        rounded_values.append(Decimal(value) / Decimal(scale))

    stable = len({tuple(row) for row in corrections}) == 1
    generic_stable = len({tuple(row) for row in generic_shell_vectors}) == 1
    if not stable or not generic_stable:
        raise ArithmeticError("closest-vector correction changed across rounding scales")
    correction = corrections[0]
    full_gram = [
        list(generic_gram[row]) + [cross[row]] for row in range(dimension)
    ] + [list(cross) + [raw]]
    representative = correction + [1]
    height = decimal_quadratic(full_gram, representative)
    generic_height = decimal_quadratic(generic_gram, generic_shell_vectors[0])
    return {
        "correction_in_specialized_generic_R17_basis": correction,
        "stable_across_integer_gram_scales": list(ROUNDING_SCALES),
        "rounded_minimum_heights": [dtext(value) for value in rounded_values],
        "canonical_height_recomputed_from_90_digit_gram": dtext(height),
        "specialized_generic_first_shell_height": dtext(generic_height),
        "specialized_generic_first_shell_vector": generic_shell_vectors[0],
        "in_or_below_specialized_generic_first_shell": bool(height <= generic_height),
    }


def division_record(polynomial) -> dict[str, Any]:
    factorization = polynomial.factor()
    group = polynomial.galois_group(pari_group=True)
    discriminant = polynomial.discriminant()
    return {
        "degree": int(polynomial.degree()),
        "defining_polynomial_coefficients_low_to_high": polynomial_key(polynomial),
        "defining_polynomial_sha256": polynomial_digest(polynomial),
        "factor_degrees_over_Q": [
            int(factor.degree()) for factor, exponent in factorization for _ in range(int(exponent))
        ],
        "irreducible_over_Q": bool(polynomial.is_irreducible()),
        "galois_group": {
            "order": int(group.order()),
            "pari_sign": int(group.signature()),
            "transitive_group_id": int(group.transitive_number()),
            "name": str(group.label()),
        },
        "maximum_coefficient_projective_height_bits": max(
            bit_height(Fraction(str(value))) for value in polynomial
        ),
        "polynomial_discriminant": str(discriminant),
        "polynomial_discriminant_bits": bit_height(discriminant),
    }


def common_frobenius_primes(records: list[dict[str, Any]], count: int = 32) -> list[int]:
    answer = []
    for prime in prime_range(5, 2000):
        if all(
            Fraction(record[f"preimage_{degree}"]["polynomial_discriminant"]).numerator
            % int(prime)
            for record in records
            for degree in (2, 3)
        ):
            answer.append(int(prime))
            if len(answer) == count:
                return answer
    raise ArithmeticError("not enough common unramified Frobenius primes")


def factor_pattern_mod_prime(coefficients: list[str], prime: int) -> str:
    ring = PolynomialRing(GF(prime), "z")
    polynomial = ring([GF(prime)(QQ(value)) for value in coefficients])
    degrees = sorted(
        int(factor.degree())
        for factor, exponent in polynomial.factor()
        for _ in range(int(exponent))
    )
    return ".".join(map(str, degrees))


def normalize_distance(distance: np.ndarray) -> np.ndarray:
    values = distance[np.triu_indices_from(distance, 1)]
    positive = values[values > 0]
    scale = float(np.median(positive)) if len(positive) else 1.0
    return distance / scale


def silhouette(distance: np.ndarray, labels: np.ndarray) -> float:
    scores = []
    for index, label in enumerate(labels):
        same = np.where(labels == label)[0]
        same = same[same != index]
        a = float(np.mean(distance[index, same])) if len(same) else 0.0
        alternatives = [
            float(np.mean(distance[index, np.where(labels == other)[0]]))
            for other in sorted(set(labels))
            if other != label
        ]
        b = min(alternatives)
        scores.append((b - a) / max(a, b) if max(a, b) else 0.0)
    return float(np.mean(scores))


def cluster_records(records: list[dict[str, Any]], quotient_gram: np.ndarray) -> dict[str, Any]:
    size = len(records)
    numeric = np.array(
        [
            [
                record["coordinate_size"]["x_projective_height_bits"],
                record["coordinate_size"]["y_projective_height_bits"],
                float(record["relative_height"]["raw_canonical_height"]),
                float(record["relative_height"]["quotient_height_defect"]),
                float(record["relative_height"]["first_shell_representative"]["canonical_height_recomputed_from_90_digit_gram"]),
                record["kummer_two_cover"]["generator_maximum_height_bits"],
                record["parent_transport"]["parent_crossratio_base_projective_height_bits"],
                float(record["parent_transport"]["parent_point_canonical_height_128bit"]),
                record["preimage_2"]["maximum_coefficient_projective_height_bits"],
                record["preimage_3"]["maximum_coefficient_projective_height_bits"],
            ]
            for record in records
        ],
        dtype=float,
    )
    deviations = numeric.std(axis=0)
    deviations[deviations == 0] = 1
    numeric = (numeric - numeric.mean(axis=0)) / deviations
    numeric_distance = np.linalg.norm(numeric[:, None, :] - numeric[None, :, :], axis=2)

    local_bits = np.array(
        [
            [int(bit) for bit in f"{int(record['local_kummer_code'], 16):053b}"]
            for record in records
        ],
        dtype=float,
    )
    local_distance = np.mean(local_bits[:, None, :] != local_bits[None, :, :], axis=2)

    division_codes = [
        [
            record[f"preimage_{degree}"]["frobenius_factorization_patterns"][prime]
            for degree in (2, 3)
            for prime in record[f"preimage_{degree}"]["frobenius_factorization_patterns"]
        ]
        for record in records
    ]
    division_distance = np.zeros((size, size), dtype=float)
    for left in range(size):
        for right in range(size):
            division_distance[left, right] = np.mean(
                [x != y for x, y in zip(division_codes[left], division_codes[right])]
            )

    diagonal = np.sqrt(np.diag(quotient_gram))
    correlations = quotient_gram / np.outer(diagonal, diagonal)
    height_distance = np.sqrt(np.maximum(0.0, 2.0 - 2.0 * correlations))

    views = {
        "arithmetic_size": normalize_distance(numeric_distance),
        "bad_place_kummer": normalize_distance(local_distance),
        "division_field_frobenius": normalize_distance(division_distance),
        "quotient_height_angle": normalize_distance(height_distance),
    }
    consensus = sum(views.values()) / len(views)
    consensus = (consensus + consensus.T) / 2.0
    np.fill_diagonal(consensus, 0.0)
    hierarchy = linkage(squareform(consensus, checks=True), method="average")
    candidates = []
    for cluster_count in range(2, 6):
        labels = fcluster(hierarchy, cluster_count, criterion="maxclust")
        actual = len(set(labels))
        if actual != cluster_count:
            continue
        candidates.append((silhouette(consensus, labels), cluster_count, labels))
    score, cluster_count, labels = max(candidates, key=lambda row: (row[0], -row[1]))
    clusters = []
    for label in sorted(set(labels)):
        members = [records[index]["label"] for index in range(size) if labels[index] == label]
        clusters.append({"cluster": int(label), "members": members, "size": len(members)})
    return {
        "status": "HEURISTIC_CONSENSUS_CLUSTERING_NOT_A_THEOREM",
        "views_equal_weight": list(views),
        "linkage": "average",
        "selected_cluster_count": int(cluster_count),
        "selected_mean_silhouette": f"{score:.17g}",
        "candidate_silhouettes": {
            str(count): f"{value:.17g}" for value, count, _labels in candidates
        },
        "clusters": clusters,
        "consensus_distance_matrix": [
            [f"{value:.12g}" for value in row] for row in consensus
        ],
        "boundary": (
            "The eleven points are a selected positive-control basis. Feature scaling, equal view "
            "weights, linkage, and silhouette selection are exploratory choices; cluster membership "
            "does not certify a shared construction."
        ),
    }


def generic_preimage_templates(model: dict[str, Any]) -> dict[str, Any]:
    coefficient_ring = PolynomialRing(QQ, ["A", "B", "X"])
    A, B, X = coefficient_ring.gens()
    field = coefficient_ring.fraction_field()
    curve = EllipticCurve(field, [0, 0, 0, A, B])

    template_records = {}
    for degree in (2, 3):
        x_map = curve.multiplication_by_m(degree)[0]
        template = x_map.numerator() - X * x_map.denominator()
        template_records[str(degree)] = {
            "formula": str(template),
            "degree_in_preimage_x": int(template.degree(template.parent().gens()[0])),
            "factor_count_over_Q_A_B_X": len(list(template.factor())),
            "irreducible_over_Q_A_B_X": len(list(template.factor())) == 1,
        }

    # Rebuild over Q(t,X) and factor after substituting the actual published
    # generic R17 coefficients.  Store hashes rather than enormous expansions.
    base = PolynomialRing(QQ, ["t", "X"])
    t, generic_x = base.gens()
    function_field = base.fraction_field()
    A_t = sum(
        function_field(QQ(value)) * function_field(t) ** index
        for index, value in enumerate(model["A_coefficients_low_to_high"])
    )
    B_t = sum(
        function_field(QQ(value)) * function_field(t) ** index
        for index, value in enumerate(model["B_coefficients_low_to_high"])
    )
    generic_curve = EllipticCurve(function_field, [0, 0, 0, A_t, B_t])
    for degree in (2, 3):
        x_map = generic_curve.multiplication_by_m(degree)[0]
        polynomial = x_map.numerator() - function_field(generic_x) * x_map.denominator()
        factors = list(polynomial.factor())
        serialized = str(polynomial)
        template_records[str(degree)].update(
            {
                "published_R17_factor_count_over_Q_t_X": len(factors),
                "published_R17_irreducible_over_Q_t_X": len(factors) == 1,
                "published_R17_expanded_polynomial_sha256": sha256(serialized.encode()).hexdigest(),
                "published_R17_expanded_serialized_bytes": len(serialized.encode()),
            }
        )
    return {
        "model": display(MODEL),
        "meaning": "F_n(z;X)=numerator(x([n]R))-X*denominator(x([n]R))",
        "templates": template_records,
        "interpretation": (
            "These irreducible generic multiplication covers reproduce the degree-4 and degree-9 "
            "preimage objects symbolically on R17. They are universal in X and therefore do not "
            "produce an exceptional specialization section by themselves."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    getcontext().prec = 100

    relation_payload = json.loads(RELATIONS.read_text())
    truth_payload = json.loads(TRUTH.read_text())
    selmer_payload = json.loads(SELMER.read_text())
    local_payload = json.loads(LOCAL_SIGNATURES.read_text())
    model_payload = json.loads(MODEL.read_text())
    module = importlib.import_module("elkies_rank28")

    relation = next(
        row for row in relation_payload["fibres"] if row["parameter"] == PARAMETER
    )
    truth = next(
        row for row in truth_payload["positive_controls"] if row["parameter"] == PARAMETER
    )
    selmer_by_label = {
        row["label"].removeprefix("public-complement-"): row
        for row in selmer_payload["candidates"]
    }
    local_by_label = {
        row["label"]: row
        for row in local_payload["public_positive_control_complement_images"]
    }

    generic_gram = [
        [Decimal(str(value)) for value in row] for row in truth["canonical_height_gram"]
    ]
    quotient_gram_decimal = [
        [Decimal(str(value)) for value in row] for row in relation["quotient_height_gram"]
    ]
    quotient_gram = np.array(
        [[float(value) for value in row] for row in quotient_gram_decimal], dtype=float
    )

    public_curve = EllipticCurve(QQ, list(module.GENERAL_WEIERSTRASS_COEFFICIENTS))
    multiplication_maps = {
        degree: public_curve.multiplication_by_m(degree) for degree in (2, 3)
    }
    records = []
    division_discriminants: dict[int, list[Fraction]] = {2: [], 3: []}
    for index, relation_record in enumerate(relation["exceptional_basis"]):
        label = relation_record["label"]
        point = tuple(Fraction(value) for value in relation_record["public_point"])
        selmer = selmer_by_label[label]
        local = local_by_label[label]
        cross = [Decimal(value) for value in relation_record["specialized_generic_height_pairings"]]
        raw = Decimal(relation_record["raw_canonical_height"])
        defect = Decimal(relation_record["quotient_height_specialization_defect"])
        shell = first_shell_stability(generic_gram, cross, raw)

        division = {}
        for degree in (2, 3):
            x_map, y_map = multiplication_maps[degree]
            polynomial = (
                x_map.numerator() - QQ(point[0]) * x_map.denominator()
            ).univariate_polynomial().monic()
            division[f"preimage_{degree}"] = division_record(polynomial)
            coordinate_ring = y_map.numerator().parent()
            _source_x, source_y = coordinate_ring.gens()
            y_coefficient = y_map.numerator().coefficient({source_y: 1})
            if y_map.numerator().degree(source_y) != 1 or y_map.denominator().degree(source_y):
                raise ArithmeticError("multiplication y-map is not affine-linear in source y")
            y_coefficient_univariate = y_coefficient.univariate_polynomial()
            division[f"preimage_{degree}"][
                "selected_preimage_y_rational_in_x_field"
            ] = bool(y_coefficient_univariate.gcd(polynomial).degree() == 0)
            division_discriminants[degree].append(
                Fraction(str(polynomial.discriminant()))
            )

        diagonal = sqrt(float(quotient_gram_decimal[index][index]))
        correlations = [
            float(quotient_gram_decimal[index][other])
            / (diagonal * sqrt(float(quotient_gram_decimal[other][other])))
            for other in range(len(quotient_gram_decimal))
        ]
        other_indices = [other for other in range(len(correlations)) if other != index]
        strongest_positive = max(other_indices, key=lambda other: correlations[other])
        strongest_negative = min(other_indices, key=lambda other: correlations[other])

        generator_coefficients = [Fraction(value) for value in selmer["generator_coefficients"]]
        record = {
            "label": label,
            "source_public_point_index_one_based": relation_record[
                "source_public_point_index_one_based"
            ],
            "public_point": [str(value) for value in point],
            "coordinate_size": {
                "x_algebraic_degree_over_Q": 1,
                "x_projective_height_bits": bit_height(point[0]),
                "y_projective_height_bits": bit_height(point[1]),
                "logarithmic_x_height_natural": f"{log(max(abs(point[0].numerator), point[0].denominator)):.17g}",
            },
            "relative_height": {
                "raw_canonical_height": relation_record["raw_canonical_height"],
                "quotient_height_defect": relation_record[
                    "quotient_height_specialization_defect"
                ],
                "orthogonal_projection_coefficients_in_generic_basis": relation_record[
                    "orthogonal_projection_coefficients_in_generic_basis"
                ],
                "first_shell_representative": shell,
                "smallest_positive_multiplier_entering_first_shell": (
                    1 if shell["in_or_below_specialized_generic_first_shell"] else None
                ),
                "all_multipliers_n_at_least_2_excluded_from_first_shell_by_continuous_quotient_lower_bound": bool(
                    Decimal(4) * defect
                    > Decimal(shell["specialized_generic_first_shell_height"])
                ),
                "intersection_analogue": {
                    "status": "HEIGHT_PAIRING_PROFILE_NOT_K3_INTERSECTION",
                    "quotient_height_correlations": [f"{value:.17g}" for value in correlations],
                    "strongest_positive": {
                        "label": relation["exceptional_basis"][strongest_positive]["label"],
                        "correlation": f"{correlations[strongest_positive]:.17g}",
                    },
                    "strongest_negative": {
                        "label": relation["exceptional_basis"][strongest_negative]["label"],
                        "correlation": f"{correlations[strongest_negative]:.17g}",
                    },
                },
            },
            "divisor_class": relation_record["divisor_class_status"],
            "kummer_two_cover": {
                "normalization": selmer["normalization"],
                "generator_coefficients": [str(value) for value in generator_coefficients],
                "generator_maximum_height_bits": max(
                    bit_height(value) for value in generator_coefficients
                ),
                "norm": selmer["norm"],
                "norm_square_root": selmer["norm_square_root"],
                "verified_rational_cover_witness": selmer["rational_cover_witness"],
                "exposes_direction": True,
            },
            "local_kummer_code": local["local"],
            "parent_transport": relation_record["parent_transport"],
            "rational_bisection_relations_involving_direction": relation_record[
                "rational_bisection_relations_involving_direction"
            ],
            **division,
        }
        records.append(record)

    frobenius_primes = common_frobenius_primes(records)
    for record in records:
        for degree in (2, 3):
            division = record[f"preimage_{degree}"]
            division["frobenius_factorization_patterns"] = {
                str(prime): factor_pattern_mod_prime(
                    division["defining_polynomial_coefficients_low_to_high"], prime
                )
                for prime in frobenius_primes
            }

    squareclass_audit = {}
    for degree in (2, 3):
        first = division_discriminants[degree][0]
        squareclass_audit[str(degree)] = {
            "all_polynomial_discriminants_same_squareclass": all(
                rational_is_square(value / first)
                for value in division_discriminants[degree]
            ),
            "common_galois_group": records[0][f"preimage_{degree}"]["galois_group"],
            "all_irreducible": all(
                record[f"preimage_{degree}"]["irreducible_over_Q"] for record in records
            ),
        }

    isogeny_matrix = pari(public_curve).ellisomat(0, 1)[1]
    isogeny_class_size = int(isogeny_matrix.nrows())
    clustering = cluster_records(records, quotient_gram)
    largest_cluster = max(cluster["size"] for cluster in clustering["clusters"])

    payload = {
        "schema": "elliptic-curves.elkies-2026-rank28-construction-fingerprints.v1",
        "status": "PASS_EXACT_KUMMER_DIVISION_FIELDS_AND_STABLE_NUMERICAL_RELATIVE_HEIGHT_FINGERPRINTS",
        "parameter": PARAMETER,
        "generic_subgroup": "specialized saturated R17 subgroup M_t of rank 17",
        "exceptional_quotient": "displayed free quotient L_t/M_t of rank 11",
        "inputs": {
            display(path): digest(path)
            for path in (RELATIONS, TRUTH, SELMER, LOCAL_SIGNATURES, MODEL, PUBLIC, Path(__file__))
        },
        "records": records,
        "division_field_common_audit": {
            "common_unramified_frobenius_primes": frobenius_primes,
            "levels": squareclass_audit,
            "full_preimage_point_field_note": (
                "The multiplication y-map is linear in the source y-coordinate, so after fixing "
                "[n]R=Q the selected y(R) is rational in x(R) and Q; the displayed degree-4/9 "
                "x-polynomial therefore defines a full preimage point field."
            ),
        },
        "rational_isogeny_graph": {
            "engine": "PARI ellisomat(E,0,1)",
            "isogeny_class_size": isogeny_class_size,
            "degree_matrix": str(isogeny_matrix),
            "nontrivial_Q_rational_isogeny_available": isogeny_class_size > 1,
        },
        "generic_R17_symbolic_reproduction": generic_preimage_templates(model_payload),
        "clustering": clustering,
        "mechanism_assessment": {
            "operation_exposing_all_eleven": (
                "The cubic 2-descent Kummer map Q |-> x(Q)-theta exposes 11/11 directions and "
                "materializes a rational two-cover for each."
            ),
            "why_this_is_not_yet_a_rank32_search_mechanism": (
                "The Kummer map represents any known rational point; it does not construct the "
                "point before it is known. The generic R17 multiplication covers are irreducible "
                "universal covers in X, not new sections."
            ),
            "small_division_field_result": (
                "All eleven directions have the same maximal observed degree/Galois type at [2] "
                "and [3], and the curve has no nontrivial Q-rational isogeny. These operations do "
                "not isolate a 7--10 direction construction class."
            ),
            "largest_exploratory_consensus_cluster": largest_cluster,
            "single_operation_explains_7_to_10_as_a_construction": False,
            "best_next_geometric_test": (
                "Search higher-genus degree-two curves or targeted higher-degree multisections whose "
                "specializations realize the exact eleven Kummer classes, using the Kummer class as "
                "the target constraint rather than enumerating parameters by Nagao score."
            ),
        },
        "claim_boundary": [
            "The public complement and generic subgroup independence are inherited exact certificates.",
            "Kummer, defining-polynomial, factorization, Galois-group, Frobenius, and isogeny data are exact.",
            "Canonical heights are PARI decimals, not interval-certified reals.",
            "Closest-vector results are complete for each rounded integral Gram and stable at three scales; this is a numerical stability certificate, not an exact canonical-height CVP theorem.",
            "A rational point on one fibre has no intrinsic divisor class or K3 intersection profile without a chosen multisection through it.",
            "The clustering is exploratory and the selected eleven-point basis is not a population sample.",
        ],
        "reproducing_command": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elliptic-curves/cas/analyze_elkies_2026_rank28_construction_fingerprints.py"
        ),
    }

    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit(f"stale or missing output: {display(args.output)}")
        print(
            "ELKIESR17R28CONSTRUCTION|status=PASS_CHECK|"
            f"output={display(args.output)}"
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        "ELKIESR17R28CONSTRUCTION|"
        f"directions={len(records)}|clusters={clustering['selected_cluster_count']}|"
        f"largest_cluster={largest_cluster}|status={payload['status']}|"
        f"output={display(args.output)}"
    )


if __name__ == "__main__":
    main()
