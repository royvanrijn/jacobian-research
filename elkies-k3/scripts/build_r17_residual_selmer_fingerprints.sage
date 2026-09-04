#!/usr/bin/env sage
"""Build exact residual-Kummer fingerprints for the pinned R17 controls.

This is deliberately richer than a carrier score or a single Selmer
dimension.  It retains the localization maps of every known residual point
class, their simultaneous kernels, delete-one-place ranks, local reduction
and component data, and every available pairing record.  It also records
finite-control p-adic residue strata and CRT prototype cylinders.  The latter
are search hypotheses, not proofs that a whole cylinder has constant Selmer
conditions.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
import json
from pathlib import Path
from typing import Any

from sage.all import ZZ, crt
from sage.version import version as sage_version
from sage.all import pari


ROOT = Path(__file__).resolve().parents[2]
SOURCE_SCRIPT = ROOT / "elkies-k3/scripts/certify_r17_074d9_local_kummer_meet.sage"
OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-residual-selmer-fingerprints-v1.json"
)
SCHEMA = "elkies-k3.r17-residual-selmer-fingerprints.v1"
PROTOCOL = "R17RESIDUALSELMERFINGERPRINT"
TARGET_IDS = (351, 356, 376, 377, 385, 12)
HIGH_GAIN_IDS = (356, 385)
LOW_GAIN_ID = 376
RESIDUAL_GAINS = {351: 8, 356: 12, 376: 5, 377: 6, 385: 12, 12: 12}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def fingerprint_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def load_source_module():
    loader = SourceFileLoader(
        "r17_074d9_local_kummer_source", str(SOURCE_SCRIPT)
    )
    spec = spec_from_loader(loader.name, loader)
    if spec is None:
        raise ImportError(f"cannot load exact Kummer source {SOURCE_SCRIPT}")
    module = module_from_spec(spec)
    loader.exec_module(module)
    return module


def stacked_localization_rank(source, places, dimension, omitted_index=None):
    rows = [[] for _ in range(dimension)]
    for index, place in enumerate(places):
        if index == omitted_index:
            continue
        local_rows = place["known_residual_localization"][
            "canonical_quotient_basis_image_rows"
        ]
        if len(local_rows) != dimension:
            raise ArithmeticError("a localization map has the wrong source dimension")
        for target, local in zip(rows, local_rows):
            target.extend(int(bit) for bit in local)
    return source.f2_rank(rows)


def localization_intersections(source, places, dimension):
    cumulative = []
    for stop in range(1, len(places) + 1):
        rank = stacked_localization_rank(source, places[:stop], dimension)
        cumulative.append(
            {
                "through_place": places[stop - 1]["rational_prime"],
                "place_kind": places[stop - 1]["place_kind"],
                "stacked_localization_rank": rank,
                "simultaneous_localization_kernel_dimension": dimension - rank,
            }
        )
    full_rank = stacked_localization_rank(source, places, dimension)
    deleted = []
    for index, place in enumerate(places):
        rank = stacked_localization_rank(source, places, dimension, index)
        deleted.append(
            {
                "deleted_place": place["rational_prime"],
                "deleted_place_kind": place["place_kind"],
                "stacked_localization_rank_after_deletion": rank,
                "rank_drop": full_rank - rank,
                "simultaneous_localization_kernel_dimension_after_deletion": (
                    dimension - rank
                ),
            }
        )
    return {
        "matrix_convention": (
            "Rows are the known residual generators and columns are concatenated "
            "coordinates in each exact local Kummer image. This is a localization "
            "matrix on the known subgroup, not the global Selmer-condition matrix."
        ),
        "source_dimension": dimension,
        "full_stacked_localization_rank": full_rank,
        "full_simultaneous_localization_kernel_dimension": dimension - full_rank,
        "cumulative_by_increasing_rational_prime": cumulative,
        "delete_one_place": deleted,
    }


def block_support_spectrum(
    places, dimension, localization_key, source_description
):
    """Enumerate the exact place-block support code on the known quotient.

    Delete-one and delete-two ranks can both be identically uninformative when
    every nonzero class is visible at several places.  Since the displayed
    residual dimensions are at most twelve, enumerating all ``2^dimension-1``
    nonzero combinations is cheap and gives the precise deletion threshold.
    """

    local_rows = [
        [
            tuple(int(bit) for bit in row)
            for row in place[localization_key][
                "canonical_quotient_basis_image_rows"
            ]
        ]
        for place in places
    ]
    if any(len(rows) != dimension for rows in local_rows):
        raise ArithmeticError("a localization map has the wrong source dimension")

    weight_counts = {}
    minimum_weight = None
    minimum_words = []
    for coefficient_mask in range(1, 1 << dimension):
        support = []
        for place, rows in zip(places, local_rows):
            width = len(rows[0]) if rows else 0
            image = [0] * width
            for source_index, row in enumerate(rows):
                if coefficient_mask & (1 << source_index):
                    image = [left ^ right for left, right in zip(image, row)]
            if any(image):
                support.append(
                    {
                        "rational_prime": place["rational_prime"],
                        "place_kind": place["place_kind"],
                    }
                )
        weight = len(support)
        weight_counts[weight] = weight_counts.get(weight, 0) + 1
        record = {
            "coefficient_mask_hex": format(coefficient_mask, "x"),
            "coefficient_vector": [
                (coefficient_mask >> index) & 1 for index in range(dimension)
            ],
            "place_block_support": support,
        }
        if minimum_weight is None or weight < minimum_weight:
            minimum_weight = weight
            minimum_words = [record]
        elif weight == minimum_weight:
            minimum_words.append(record)

    if minimum_weight is None:
        minimum_weight = 0
    return {
        "code_convention": (
            f"The image of each nonzero {source_description} class is treated as a "
            "block word with one block per audited rational place. Block weight "
            "counts places with nonzero localization, not local coordinates."
        ),
        "source_dimension": dimension,
        "nonzero_class_count": (1 << dimension) - 1,
        "audited_place_count": len(places),
        "block_weight_enumerator": {
            str(weight): count for weight, count in sorted(weight_counts.items())
        },
        "minimum_nonzero_place_block_weight": minimum_weight,
        "every_deletion_of_at_most_this_many_places_preserves_injectivity": max(
            0, minimum_weight - 1
        ),
        "minimum_support_word_count": len(minimum_words),
        "minimum_support_words": minimum_words,
        "claim_boundary": (
            f"This is the exact localization support code of the certified "
            f"{source_description} point subgroup, not of the complete residual "
            "Selmer group."
        ),
    }


def place_block_dual_components(
    places, dimension, localization_key, source_description
):
    """Compute connected components of the block-localization dual matroid."""

    def independent_basis(vectors):
        pivots = {}
        basis = []
        for vector in vectors:
            reduced = vector
            while reduced and reduced.bit_length() - 1 in pivots:
                reduced ^= pivots[reduced.bit_length() - 1]
            if reduced:
                pivots[reduced.bit_length() - 1] = reduced
                basis.append(vector)
        return basis

    elements = []
    zero_places = []
    for place_index, place in enumerate(places):
        rows = place[localization_key]["canonical_quotient_basis_image_rows"]
        if len(rows) != dimension:
            raise ArithmeticError("a localization map has the wrong source dimension")
        width = len(rows[0]) if rows else 0
        dual_vectors = [
            sum((int(rows[index][column]) & 1) << index for index in range(dimension))
            for column in range(width)
        ]
        local_basis = independent_basis(vector for vector in dual_vectors if vector)
        if not local_basis:
            zero_places.append(place["rational_prime"])
        for vector in local_basis:
            elements.append((vector, place_index))

    parents = list(range(len(elements)))

    def find(index):
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left, right):
        left = find(left)
        right = find(right)
        if left != right:
            parents[right] = left

    elements_by_place = {}
    for element_index, (_vector, place_index) in enumerate(elements):
        elements_by_place.setdefault(place_index, []).append(element_index)
    for local_elements in elements_by_place.values():
        for element_index in local_elements[1:]:
            union(local_elements[0], element_index)

    pivots = {}
    global_basis_elements = []
    for element_index, (vector, _place_index) in enumerate(elements):
        reduced = vector
        combination = 0
        while reduced and reduced.bit_length() - 1 in pivots:
            pivot_vector, pivot_combination = pivots[reduced.bit_length() - 1]
            reduced ^= pivot_vector
            combination ^= pivot_combination
        if reduced:
            basis_index = len(global_basis_elements)
            global_basis_elements.append(element_index)
            pivots[reduced.bit_length() - 1] = (
                reduced,
                combination ^ (1 << basis_index),
            )
        else:
            for basis_index, basis_element in enumerate(global_basis_elements):
                if combination & (1 << basis_index):
                    union(element_index, basis_element)

    grouped = {}
    for element_index in range(len(elements)):
        grouped.setdefault(find(element_index), []).append(element_index)
    components = []
    component_dimension_sum = 0
    for component_elements in grouped.values():
        vectors = [elements[index][0] for index in component_elements]
        component_dimension = len(independent_basis(vectors))
        component_dimension_sum += component_dimension
        place_indices = sorted({elements[index][1] for index in component_elements})
        components.append(
            {
                "dual_dimension": component_dimension,
                "local_dual_basis_vector_count": len(component_elements),
                "places": [
                    {
                        "rational_prime": places[index]["rational_prime"],
                        "place_kind": places[index]["place_kind"],
                    }
                    for index in place_indices
                ],
            }
        )
    components.sort(
        key=lambda component: (
            -component["dual_dimension"],
            canonical_text(component["places"]),
        )
    )
    global_dual_rank = len(global_basis_elements)
    if component_dimension_sum != global_dual_rank:
        raise ArithmeticError("dual-matroid component dimensions do not add up")
    return {
        "construction": (
            "Take a basis of image(localization_v^*) for each place, force all "
            "basis vectors from one place into one block, and compute the connected "
            "components of their binary vector matroid via fundamental circuits."
        ),
        "source_description": source_description,
        "source_dimension": dimension,
        "global_localization_dual_rank": global_dual_rank,
        "local_dual_basis_vector_count": len(elements),
        "zero_localization_places": zero_places,
        "component_count": len(components),
        "components": components,
        "place_block_matroid_indecomposable": (
            global_dual_rank == dimension and len(components) == 1
        ),
        "claim_boundary": (
            "Indecomposable here excludes a direct-sum partition separated by "
            "audited place blocks on the certified point subspace. It is not a "
            "Cassels-pairing or complete residual-Selmer indecomposability theorem."
        ),
    }


def two_part(value: int) -> int:
    value = abs(int(value))
    result = 1
    while value and value % 2 == 0:
        result *= 2
        value //= 2
    return result


def place_feature(place):
    local = place["known_residual_localization"]
    selected = place["selected_block_localization"]
    directions = place["directions"]
    return {
        "place_kind": place["place_kind"],
        "reduction_kind": place["reduction_kind"],
        "kodaira_symbol": place["kodaira_symbol"],
        "conductor_exponent": int(place["conductor_exponent"]),
        "minimal_discriminant_valuation": int(place["minimal_discriminant_valuation"]),
        "tamagawa_two_part": two_part(place["tamagawa_number"]),
        "ambient_local_kummer_dimension": int(place["ambient_local_kummer_dimension"]),
        "known_residual_image_dimension": int(local["quotient_image_dimension"]),
        "known_residual_localization_kernel_dimension": len(
            local["kernel_in_canonical_quotient_coordinates"]
        ),
        "selected_block_image_dimension": int(selected["quotient_image_dimension"]),
        "component_image_order_multiset": sorted(
            int(row["component_group_image_order"]) for row in directions
        ),
        "nontrivial_component_hilbert_pair_count": sum(
            any(
                item["value"]["hilbert_symbol"] == -1
                for item in row["component_hilbert_symbol_multiset"]
            )
            for row in place["pairwise_hilbert_symbols"]
        ),
    }


def compact_place(place):
    """Retain subspaces and invariant summaries; hash the bulky pair table."""

    feature = place_feature(place)
    return {
        "rational_prime": place["rational_prime"],
        **feature,
        "local_factor_descriptors": place["local_factor_descriptors"],
        "known_residual_localization": place["known_residual_localization"],
        "selected_block_localization": place["selected_block_localization"],
        "direction_invariants": [
            {
                "label": row["label"],
                "valuation_parity_support": row["valuation_parity_support"],
                "component_group_image_order": row["component_group_image_order"],
                "odd_local_squareclass_multiset": row.get(
                    "odd_local_squareclass_multiset"
                ),
                "self_component_hilbert_symbol_multiset": row[
                    "self_component_hilbert_symbol_multiset"
                ],
            }
            for row in place["directions"]
        ],
        "pairing_summary": {
            "all_corestricted_local_tate_symbols_trivial": place[
                "all_pairwise_corestricted_local_tate_symbols_trivial"
            ],
            "unordered_pair_count": len(place["pairwise_hilbert_symbols"]),
            "nontrivial_component_hilbert_pair_count": feature[
                "nontrivial_component_hilbert_pair_count"
            ],
            "full_pair_table_sha256": fingerprint_hash(
                place["pairwise_hilbert_symbols"]
            ),
            "full_pair_table_source": (
                "artifacts/generated-results/"
                "elkies-k3-r17-074d9-local-kummer-meet-v1.json"
            ),
        },
    }


def p_adic_residue(parameter: str, prime: int, exponent: int):
    value = Fraction(parameter)
    numerator = value.numerator
    denominator = value.denominator
    modulus = prime**exponent
    if denominator % prime:
        residue = numerator * pow(denominator, -1, modulus) % modulus
        return {"chart": "affine_t", "residue": residue, "modulus": modulus}
    if numerator % prime:
        residue = denominator * pow(numerator, -1, modulus) % modulus
        return {"chart": "infinity_s_equals_1_over_t", "residue": residue, "modulus": modulus}
    raise ArithmeticError("a reduced projective parameter has no unit p-adic chart")


def parameter_strata(curves, common_primes):
    records = []
    for prime in common_primes:
        for exponent in (1, 2, 3):
            groups = {}
            samples = []
            for curve_id in TARGET_IDS:
                curve = curves[curve_id]
                residue = p_adic_residue(curve["family_parameter"], prime, exponent)
                place = next(
                    row
                    for row in curve["local_places"]
                    if int(row["rational_prime"]) == prime
                )
                feature = place_feature(place)
                key = (residue["chart"], residue["residue"])
                groups.setdefault(key, []).append(curve_id)
                samples.append(
                    {
                        "curve_id": curve_id,
                        "residual_gain": RESIDUAL_GAINS[curve_id],
                        **residue,
                        "local_feature_hash": fingerprint_hash(feature),
                    }
                )
            records.append(
                {
                    "prime": prime,
                    "exponent": exponent,
                    "samples": samples,
                    "shared_residue_groups": [
                        {"chart": key[0], "residue": key[1], "curve_ids": ids}
                        for key, ids in sorted(groups.items())
                        if len(ids) > 1
                    ],
                }
            )
    return records


def high12_vs_low5(curves):
    indices = {
        curve_id: {
            int(place["rational_prime"]): place
            for place in curves[curve_id]["local_places"]
        }
        for curve_id in (*HIGH_GAIN_IDS, LOW_GAIN_ID)
    }
    common = sorted(set.intersection(*(set(row) for row in indices.values())))
    scalar_fields = (
        "reduction_kind",
        "kodaira_symbol",
        "conductor_exponent",
        "minimal_discriminant_valuation",
        "tamagawa_two_part",
        "ambient_local_kummer_dimension",
        "known_residual_image_dimension",
        "known_residual_localization_kernel_dimension",
        "selected_block_image_dimension",
        "nontrivial_component_hilbert_pair_count",
    )
    comparisons = []
    for prime in common:
        features = {
            curve_id: place_feature(indices[curve_id][prime])
            for curve_id in (*HIGH_GAIN_IDS, LOW_GAIN_ID)
        }
        separating = []
        for field in scalar_fields:
            high_values = {canonical_text(features[curve_id][field]) for curve_id in HIGH_GAIN_IDS}
            low_value = canonical_text(features[LOW_GAIN_ID][field])
            if low_value not in high_values:
                separating.append(
                    {
                        "feature": field,
                        "high_gain_values": {
                            str(curve_id): features[curve_id][field]
                            for curve_id in HIGH_GAIN_IDS
                        },
                        "low_gain_value": features[LOW_GAIN_ID][field],
                    }
                )
        comparisons.append(
            {
                "prime": prime,
                "features": {str(key): value for key, value in features.items()},
                "scalar_features_separating_both_plus12_samples_from_plus5": separating,
            }
        )
    return {
        "high_gain_controls": list(HIGH_GAIN_IDS),
        "high_gain": 12,
        "low_gain_control": LOW_GAIN_ID,
        "low_gain": 5,
        "common_exactly_audited_places": common,
        "placewise_comparisons": comparisons,
        "boundary": (
            "Two high-gain samples and one low-gain sample identify discriminating "
            "features, not a statistical law or a sufficient Selmer criterion."
        ),
    }


def crt_prototypes(curves, comparison):
    discriminating_primes = [
        row["prime"]
        for row in comparison["placewise_comparisons"]
        if row["scalar_features_separating_both_plus12_samples_from_plus5"]
    ]
    prototypes = []
    for curve_id in HIGH_GAIN_IDS:
        parameter = curves[curve_id]["family_parameter"]
        chosen = []
        residues = []
        moduli = []
        for prime in discriminating_primes:
            residue = p_adic_residue(parameter, prime, 3)
            if residue["chart"] != "affine_t":
                continue
            chosen.append(prime)
            residues.append(ZZ(residue["residue"]))
            moduli.append(ZZ(residue["modulus"]))
            if len(chosen) == 5:
                break
        if not chosen:
            continue
        combined = ZZ(crt(residues, moduli))
        modulus = ZZ(1)
        for value in moduli:
            modulus *= value
        if combined > modulus // 2:
            combined -= modulus
        prototypes.append(
            {
                "prototype_curve_id": curve_id,
                "prime_power_conditions": [
                    {"prime": prime, "exponent": 3, "residue": int(residue), "modulus": int(modulus_i)}
                    for prime, residue, modulus_i in zip(chosen, residues, moduli)
                ],
                "combined_integer_parameter_class": {
                    "residue": str(combined),
                    "modulus": str(modulus),
                    "parameter_family": f"t = {combined} + {modulus}*n",
                },
                "status": "EXACT_CRT_CLASS_HEURISTIC_FINGERPRINT_PRESERVATION",
            }
        )
    return {
        "prototypes": prototypes,
        "claim_boundary": (
            "CRT exactly preserves the displayed p-adic parameter residues. It does "
            "not yet prove that the residual Kummer or Selmer fingerprints are "
            "constant on those cylinders; candidates require exact re-audit."
        ),
    }


def build():
    source = load_source_module()
    full = source.build()
    curves = {int(row["curve_id"]): row for row in full["curves"]}
    if set(curves) != set(TARGET_IDS):
        raise ArithmeticError("the exact local-control set changed")
    fingerprints = []
    for curve_id in TARGET_IDS:
        curve = curves[curve_id]
        dimension = len(curve["certified_displayed_exceptional_quotient_basis"])
        if dimension != RESIDUAL_GAINS[curve_id]:
            raise ArithmeticError("a certified displayed residual gain changed")
        fingerprints.append(
            {
                "curve_id": curve_id,
                "role": curve["role"],
                "native_chart": curve["native_chart"],
                "family_parameter": curve["family_parameter"],
                "certified_known_residual_dimension": dimension,
                "known_residual_basis": curve[
                    "certified_displayed_exceptional_quotient_basis"
                ],
                "selected_comparison_block": curve["selected_comparison_block"],
                "two_division_etale_algebra": curve["two_division_etale_algebra"],
                "kummer_images": curve["kummer_images"],
                "local_places": [compact_place(place) for place in curve["local_places"]],
                "localization_intersections": localization_intersections(
                    source, curve["local_places"], dimension
                ),
                "known_residual_place_block_support_code": block_support_spectrum(
                    curve["local_places"],
                    dimension,
                    "known_residual_localization",
                    "known residual",
                ),
                "selected_comparison_block_place_support_code": block_support_spectrum(
                    curve["local_places"],
                    len(curve["selected_comparison_block"]),
                    "selected_block_localization",
                    "selected comparison block",
                ),
                "known_residual_place_block_dual_components": place_block_dual_components(
                    curve["local_places"],
                    dimension,
                    "known_residual_localization",
                    "known residual",
                ),
                "selected_comparison_block_dual_components": place_block_dual_components(
                    curve["local_places"],
                    len(curve["selected_comparison_block"]),
                    "selected_block_localization",
                    "selected comparison block",
                ),
                "complete_two_selmer_status": "UNKNOWN_BNF_OR_COMPLETE_DESCENT_REQUIRED",
                "cassels_pairing_status": (
                    {
                        "status": full["record_quotient_arithmetic_blocks"][
                            str(curve_id)
                        ]["status"],
                        "pairing_descent_obstruction_span": full[
                            "record_quotient_arithmetic_blocks"
                        ][str(curve_id)]["pairing_descent_obstruction_span"],
                        "pairing_descent_obstruction_primes": sorted(
                            {
                                row["rational_prime"]
                                for row in full["record_quotient_arithmetic_blocks"][
                                    str(curve_id)
                                ]["pairing_descent_obstructions"]
                            },
                            key=ZZ,
                        ),
                        "corestricted_local_tate_pairing_control": full[
                            "record_quotient_arithmetic_blocks"
                        ][str(curve_id)]["corestricted_local_tate_pairing_control"],
                        "indecomposable_components": full[
                            "record_quotient_arithmetic_blocks"
                        ][str(curve_id)]["indecomposable_components"],
                        "full_pairing_source": (
                            "artifacts/generated-results/"
                            "elkies-k3-r17-074d9-quotient-arithmetic-blocks-v1.json"
                        ),
                    }
                    if curve_id in HIGH_GAIN_IDS
                    else {
                        "status": "NOT_COMPUTED_ON_COMPLETE_RESIDUAL_SELMER_GROUP",
                        "known_point_local_tate_pairings": "all corestricted local pairings are zero",
                    }
                ),
            }
        )

    common_primes = sorted(
        set.intersection(
            *(
                {int(place["rational_prime"]) for place in curves[curve_id]["local_places"]}
                for curve_id in TARGET_IDS
            )
        )
    )
    comparison = high12_vs_low5(curves)
    output = {
        "schema": SCHEMA,
        "status": "PASS_EXACT_KNOWN_RESIDUAL_KUMMER_FINGERPRINTS_FULL_SELMER_UNKNOWN",
        "summary": {
            "control_curve_ids": list(TARGET_IDS),
            "known_residual_dimensions": {
                str(curve_id): RESIDUAL_GAINS[curve_id] for curve_id in TARGET_IDS
            },
            "complete_two_selmer_groups_computed": False,
            "full_selmer_leave_one_place_out_matrices_computed": False,
            "known_residual_localization_matrices_computed": True,
            "high12_vs_low5_common_place_count": len(
                comparison["common_exactly_audited_places"]
            ),
        },
        "fingerprints": fingerprints,
        "plus12_vs_plus5": comparison,
        "control_parameter_residue_stratification": {
            "common_exactly_audited_primes": common_primes,
            "exponents": [1, 2, 3],
            "records": parameter_strata(curves, common_primes),
            "scope": (
                "Exact residues and exact fingerprints at the finitely many controls; "
                "no local constancy on an entire residue cylinder is inferred."
            ),
        },
        "crt_search_prototypes": crt_prototypes(curves, comparison),
        "monotone_sieve_interface": {
            "finite_upper_bound_policy": (
                "Only a proved global residual upper bound may reject a candidate."
            ),
            "missing_bnf_policy": (
                "No finite upper bound yet; bounded search may proceed but no Selmer "
                "or rank claim follows."
            ),
            "gate_module": "elliptic-curves/cas/elkies_residual_selmer_gate.py",
        },
        "claim_boundary": [
            "These are exact Kummer localizations of certified known point classes, not complete local conditions on the unknown global 2-Selmer ambient space.",
            "The simultaneous kernels and delete-one-place ranks are for the known residual subgroup; the complete Selmer analogues remain contingent on the checkpointed BNF/descent worker.",
            "The componentwise Hilbert forms for curves 356 and 385 fail to descend through the rigid plane, so no intrinsic indecomposable ten-block claim is made.",
            "Sampled p-adic residue correlations and CRT classes guide bounded search only and must be re-certified at every generated parameter.",
        ],
        "inputs": {
            str(SOURCE_SCRIPT.relative_to(ROOT)): digest(SOURCE_SCRIPT),
            **full["inputs"],
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "build_r17_residual_selmer_fingerprints.sage --check"
        ),
        "software_assumptions": {
            "sage": str(sage_version),
            "pari": ".".join(str(part) for part in pari.version()),
        },
    }
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    document = build()
    serialized = json.dumps(document, indent=2, sort_keys=True) + "\n"
    output = args.output.resolve()
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored residual-Selmer fingerprints differ from replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        f"{PROTOCOL}|controls={len(TARGET_IDS)}|status={document['status']}|"
        f"output={output.relative_to(ROOT)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
