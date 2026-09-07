#!/usr/bin/env sage-python
"""Prospectively audit and refine the two R17 CRT local cylinders.

This is Phase 1 of the rank-jump experiment.  It never reads a point-search
result.  Every sampled parameter is specialized on the exact 074d9 R17
model, all seventeen generic sections are checked, and the selected local
fingerprint is recomputed from the specialization.  The original p^3
prototypes are audited before a prime-by-prime, local-outcome-only refinement.

The local number-field order is made p-maximal by supplying the one rational
prime under audit to PARI's round-two algorithm.  No global class group,
Selmer dimension, or rank upper bound is computed or inferred.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any

from sage.all import EllipticCurve, GF, Matrix, PolynomialRing, QQ, ZZ, crt, inverse_mod, pari
from sage.version import version as sage_version


ROOT = Path(__file__).resolve().parents[2]
LINEAGE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
PROTOTYPES = ROOT / "artifacts/generated-results/elkies-k3-r17-residual-selmer-fingerprints-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-local-stability-v1.json"

SCHEMA = "elkies-k3.r17-prospective-crt-local-stability.v1"
PROTOCOL = "R17CRTLOCALSTABILITY"
SALT = "r17-prospective-crt-local-stability-v1"
ANCHORS = {
    356: (2, 13, 37, 53, 71),
    385: (2, 13, 37, 53, 67),
}
INITIAL_SAMPLE_SIZE = 64
DISCOVERY_DRAWS_PER_EXPONENT = 16
CONFIRMATION_DRAWS = 64
MIN_EXPONENT = 3
MAX_EXPONENT = 20
FROZEN_PROTOTYPE_WHOLE_FILE_SHA256 = "54548e6b7110d0b53ae3bd86a97bbd06fd1159836d19c3fc3ad4e23b77320fbc"
FROZEN_GENERATOR_SHA256 = "baeaeed989cc20835d2cc51609c9578081a4771a3f12396276018c1cbb7ce535"

sys.path.insert(0, str(ROOT / "elliptic-curves/cas"))
from run_fermigier_rank20_auxiliary_fingerprints import (  # noqa: E402
    qpari,
    two_adic_coords,
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def rational_text(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def f2_rref(rows: list[list[int]], width: int | None = None) -> list[list[int]]:
    if width is None:
        width = len(rows[0]) if rows else 0
    if not rows:
        return []
    matrix = Matrix(GF(2), rows, ncols=width).echelon_form()
    return [[int(value) for value in row] for row in matrix.rows() if any(row)]


def f2_rank(rows: list[list[int]]) -> int:
    return len(f2_rref(rows))


def source_kernel(rows: list[list[int]]) -> list[list[int]]:
    if not rows:
        return [[int(i == j) for i in range(17)] for j in range(17)]
    matrix = Matrix(GF(2), rows)
    return [
        [int(value) for value in row]
        for row in matrix.left_kernel().basis_matrix().echelon_form().rows()
        if any(row)
    ]


def hashed_integer(tag: str, index: int, bound: int) -> int:
    block = sha256(f"{SALT}|{tag}|{index}".encode()).digest()
    return int.from_bytes(block, "big") % (2 * bound + 1) - bound


def p_adic_residue(value, prime: int, exponent: int) -> int:
    value = QQ(value)
    modulus = ZZ(prime) ** exponent
    if ZZ(value.denominator()) % prime == 0:
        raise ArithmeticError("the selected target parameter is not in the affine p-adic chart")
    return int(ZZ(value.numerator()) * inverse_mod(ZZ(value.denominator()), modulus) % modulus)


def polynomial_value(coefficients, parameter):
    parameter = QQ(parameter)
    return sum(QQ(coefficient) * parameter**index for index, coefficient in enumerate(coefficients))


def rational_valuation(value, prime: int) -> int:
    value = QQ(value)
    return int(ZZ(value.numerator()).valuation(prime) - ZZ(value.denominator()).valuation(prime))


def point_in_E0(point, prime: int, ainvs) -> bool:
    if point.is_zero():
        return True
    x_coordinate, y_coordinate = map(QQ, point[:2])
    x_valuation = rational_valuation(x_coordinate, prime)
    y_valuation = rational_valuation(y_coordinate, prime)
    if x_valuation < 0 or y_valuation < 0:
        if x_valuation < 0 and y_valuation < 0:
            return True
        raise ArithmeticError("inconsistent negative valuations on a minimal model")
    p = ZZ(prime)
    x_residue = ZZ(x_coordinate.numerator()) * inverse_mod(ZZ(x_coordinate.denominator()), p) % p
    y_residue = ZZ(y_coordinate.numerator()) * inverse_mod(ZZ(y_coordinate.denominator()), p) % p
    a1, a2, a3, a4, _a6 = [ZZ(value) % p for value in ainvs]
    derivative_x = (a1 * y_residue - 3 * x_residue**2 - 2 * a2 * x_residue - a4) % p
    derivative_y = (2 * y_residue + a1 * x_residue + a3) % p
    return derivative_x != 0 or derivative_y != 0


def component_order(point, prime: int, tamagawa: int, ainvs) -> int:
    for candidate in ZZ(tamagawa).divisors():
        if point_in_E0(candidate * point, prime, ainvs):
            return int(candidate)
    raise ArithmeticError("the component order does not divide the Tamagawa number")


def kodaira_name(local_data) -> str:
    return str(local_data.kodaira_symbol())


class Family:
    def __init__(self):
        lineage = json.loads(LINEAGE.read_text())
        if lineage.get("status") != "PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS":
            raise ArithmeticError("the exact 074d9 lineage input is not certified")
        representative = lineage["representative"]
        if representative["chart"] != "norm12-orbit-074d9":
            raise ArithmeticError("the representative chart changed")
        self.a_coefficients = representative["A_coefficients_low_to_high"]
        self.b_coefficients = representative["B_coefficients_low_to_high"]
        self.sections = lineage["native_chart_sections"]["norm12-orbit-074d9"]
        if [int(row["section"]) for row in self.sections] != list(range(1, 18)):
            raise ArithmeticError("the ordered generic MW17 basis changed")
        self.target_parameters = {
            int(row["curve_id"]): QQ(row["parameter"])
            for row in lineage["target_isomorphisms"]
            if row["chart"] == "norm12-orbit-074d9"
        }

    def specialize(self, parameter):
        parameter = QQ(parameter)
        coefficient_a = polynomial_value(self.a_coefficients, parameter)
        coefficient_b = polynomial_value(self.b_coefficients, parameter)
        curve = EllipticCurve(QQ, [0, 0, 0, coefficient_a, coefficient_b])
        if curve.discriminant() == 0:
            raise ArithmeticError("singular specialization")
        points = []
        for expected_index, section in enumerate(self.sections, start=1):
            x_coordinate = polynomial_value(section["x_coefficients_low_to_high"], parameter)
            y_coordinate = polynomial_value(section["y_coefficients_low_to_high"], parameter)
            if y_coordinate**2 != x_coordinate**3 + coefficient_a * x_coordinate + coefficient_b:
                raise ArithmeticError(f"generic section S{expected_index} failed exact specialization")
            points.append(curve(x_coordinate, y_coordinate))
        return curve, points


def local_place_data(nf, prime: int):
    places = list(pari.idealprimedec(nf, prime))
    descriptors = [
        {
            "ramification_index": int(place[2]),
            "residue_degree": int(place[3]),
            "local_degree": int(place[2]) * int(place[3]),
        }
        for place in places
    ]
    if sum(row["local_degree"] for row in descriptors) != 3:
        raise ArithmeticError("the local degrees do not sum to three")
    return places, descriptors


def odd_factor_block(nf, place, descriptor, alphas):
    uniformizer_column = pari.idealappr(nf, place)
    uniformizer = pari.nfbasistoalg(nf, uniformizer_column)
    if int(pari.idealval(nf, uniformizer, place)) != 1:
        raise ArithmeticError("PARI returned a non-uniformizer")
    reduction = pari.nfmodprinit(nf, place)
    rows = []
    for alpha in alphas:
        valuation = int(pari.idealval(nf, alpha, place))
        unit = alpha / uniformizer**valuation
        residue = pari.nfmodpr(nf, unit, reduction)
        rows.append([valuation & 1, 0 if bool(pari.issquare(residue)) else 1])
    return {
        "descriptor": descriptor,
        "known_mw17_squareclass_rows": rows,
    }


def local_fingerprint(curve, points, prime: int) -> dict[str, Any]:
    local_data = curve.local_data(prime)
    minimal_curve = local_data.minimal_model()
    isomorphisms = curve.isomorphisms(minimal_curve)
    if not isomorphisms:
        raise ArithmeticError("no exact isomorphism to the local minimal model")
    minimal_points = [isomorphisms[0](point) for point in points]
    ainvs = tuple(minimal_curve.a_invariants())
    if any(value.denominator() != 1 for value in ainvs):
        raise ArithmeticError("the local minimal model is not integral")
    a1, a2, a3, a4, a6 = ainvs
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3**2 + 4 * a6
    polynomial_ring = PolynomialRing(ZZ, "z")
    z = polynomial_ring.gen()
    polynomial = z**3 + ZZ(b2) * z**2 + ZZ(8 * b4) * z + ZZ(16 * b6)
    if polynomial.discriminant() != 256 * minimal_curve.discriminant():
        raise ArithmeticError("the completed-square cubic discriminant identity failed")

    # The supplied list makes the order p-maximal without asking PARI to
    # factor the several-thousand-bit global discriminant.
    nf = pari.nfinit([pari(polynomial), [prime]])
    places, descriptors = local_place_data(nf, prime)
    theta = pari(f"Mod(z,{polynomial})")
    alphas = [qpari(pari, 4 * point[0]) - theta for point in minimal_points]
    norms = [QQ(str(pari.nfeltnorm(nf, alpha))) for alpha in alphas]
    for point, norm in zip(minimal_points, norms):
        x_coordinate, y_coordinate = map(QQ, point[:2])
        expected_root = 4 * (2 * y_coordinate + a1 * x_coordinate + a3)
        if norm != expected_root**2:
            raise ArithmeticError("a generic MW17 Kummer norm is not the expected square")

    factor_blocks = []
    if prime == 2:
        basis, origins, rows = two_adic_coords(pari, nf, places, alphas)
        for place, descriptor in zip(places, descriptors):
            factor_blocks.append(
                {
                    "descriptor": descriptor,
                }
            )
        factor_blocks.sort(key=canonical_text)
        coordinate_description = {
            "kind": "greedy_exact_product_Q2_squareclass_basis_in_section_order",
            "basis_origin_section_indices_one_based": [int(index) + 1 for index in origins],
        }
    else:
        factor_blocks = [
            odd_factor_block(nf, place, descriptor, alphas)
            for place, descriptor in zip(places, descriptors)
        ]
        factor_blocks.sort(key=canonical_text)
        raw_rows = [
            [bit for block in factor_blocks for bit in block["known_mw17_squareclass_rows"][index]]
            for index in range(17)
        ]
        rows = raw_rows
        origins = []
        coordinate_description = {
            "kind": "source_kernel_invariant_under_local_factor_and_uniformizer_coordinate_changes",
            "raw_odd_factor_coordinates_retained_outside_comparison_hash": True,
        }

    image_dimension = f2_rank(rows)
    ambient_dimension = len(places) - 1 + int(prime == 2)
    if image_dimension > ambient_dimension:
        raise ArithmeticError("the known MW17 image exceeds E(Q_p)/2E(Q_p)")
    tamagawa = int(local_data.tamagawa_number())
    component_orders = [component_order(point, prime, tamagawa, ainvs) for point in minimal_points]
    comparison_payload = {
        "reduction": {
            "kodaira_symbol": kodaira_name(local_data),
            "conductor_exponent": int(local_data.conductor_valuation()),
            "minimal_discriminant_valuation": int(local_data.discriminant_valuation()),
            "tamagawa_number": tamagawa,
        },
        "local_factor_descriptors": sorted(descriptors, key=canonical_text),
        "ambient_local_kummer_dimension": ambient_dimension,
        "known_mw17_localization": {
            "image_dimension": image_dimension,
            "source_kernel_rref": source_kernel(rows),
            "canonical_source_order_relation_rows": rows if prime == 2 else None,
            "coordinate_description": coordinate_description,
        },
        "known_mw17_component_orders_by_section": component_orders,
        "hilbert_tate_panel": {
            "status": "DEFERRED_TO_FROZEN_PHASE3_RESOURCE_BOUNDED_PANEL",
            "reason": (
                "The prospective local-constancy test uses the actual Kummer subspace. "
                "All-pairs factorwise nfhilbert on every refinement draw is not needed "
                "to select the cylinder and is deliberately not used for tuning."
            ),
        },
    }
    return {
        "rational_prime": prime,
        "p_maximal_round_two_prime": prime,
        "diagnostic_anonymous_factor_blocks": factor_blocks,
        "diagnostic_raw_odd_squareclass_rows": raw_rows if prime != 2 else None,
        "comparison_payload": comparison_payload,
        "comparison_sha256": canonical_hash(comparison_payload),
    }


def audit_parameter(family: Family, parameter, primes, target_hashes) -> dict[str, Any]:
    parameter = QQ(parameter)
    record = {
        "parameter": rational_text(parameter),
        "projective_pair": [int(parameter.numerator()), int(parameter.denominator())],
        "status": None,
        "failure": None,
        "local_fingerprints": [],
    }
    try:
        curve, points = family.specialize(parameter)
        record["specialization"] = {
            "canonical_short_ainvs": [rational_text(value) for value in curve.a_invariants()],
            "discriminant": rational_text(curve.discriminant()),
            "generic_mw17_section_count_verified": len(points),
        }
        for prime in primes:
            fingerprint = local_fingerprint(curve, points, prime)
            fingerprint["matches_anchor_target"] = fingerprint["comparison_sha256"] == target_hashes[str(prime)]
            record["local_fingerprints"].append(fingerprint)
        record["status"] = "PASS_VALID_EXACT_SPECIALIZATION_AND_LOCAL_AUDIT"
    except Exception as exc:
        record["status"] = "FAIL_STRUCTURAL_OR_LOCAL_AUDIT"
        record["failure"] = {"type": type(exc).__name__, "message": str(exc)}
    return record


def compact_outcome(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "parameter": record["parameter"],
        "status": record["status"],
        "failure": record["failure"],
        "local_matches": {
            str(row["rational_prime"]): row["matches_anchor_target"]
            for row in record["local_fingerprints"]
        },
        "local_hashes": {
            str(row["rational_prime"]): row["comparison_sha256"]
            for row in record["local_fingerprints"]
        },
    }


def audit_structure(family: Family, parameter, prime: int, target_reduction) -> dict[str, Any]:
    """Cheap necessary gate before a full local Kummer computation."""

    parameter = QQ(parameter)
    record = {
        "parameter": rational_text(parameter),
        "status": None,
        "failure": None,
        "reduction": None,
        "matches_anchor_target_reduction": False,
        "full_local_fingerprint_computed": False,
    }
    try:
        curve, points = family.specialize(parameter)
        local_data = curve.local_data(prime)
        reduction = {
            "kodaira_symbol": kodaira_name(local_data),
            "conductor_exponent": int(local_data.conductor_valuation()),
            "minimal_discriminant_valuation": int(local_data.discriminant_valuation()),
            "tamagawa_number": int(local_data.tamagawa_number()),
        }
        record["status"] = "PASS_VALID_EXACT_SPECIALIZATION_AND_REDUCTION_AUDIT"
        record["generic_mw17_section_count_verified"] = len(points)
        record["reduction"] = reduction
        record["matches_anchor_target_reduction"] = reduction == target_reduction
    except Exception as exc:
        record["status"] = "FAIL_STRUCTURAL_OR_LOCAL_REDUCTION_AUDIT"
        record["failure"] = {"type": type(exc).__name__, "message": str(exc)}
    return record


def initial_sample(family, anchor: int, prototype, target_hashes):
    base = ZZ(prototype["combined_integer_parameter_class"]["residue"])
    modulus = ZZ(prototype["combined_integer_parameter_class"]["modulus"])
    records = []
    for index in range(INITIAL_SAMPLE_SIZE):
        multiplier = hashed_integer(f"initial-p3|{anchor}", index, 2**15)
        parameter = base + modulus * multiplier
        record = audit_parameter(family, parameter, ANCHORS[anchor], target_hashes)
        record["sample_index"] = index
        record["cylinder_multiplier"] = multiplier
        records.append(record)
    return records


def refinement_draw_parameter(target, anchor: int, prime: int, exponent: int, lane: str, index: int):
    modulus = ZZ(prime) ** exponent
    residue = ZZ(p_adic_residue(target, prime, exponent))
    multiplier = hashed_integer(f"refine|{lane}|{anchor}|{prime}|{exponent}", index, 2**14)
    return residue + modulus * multiplier, multiplier


def refine_prime(family, anchor: int, prime: int, target, target_hash: str, target_reduction):
    exponents = []
    selected = None
    selected_confirmation = None
    for exponent in range(MIN_EXPONENT, MAX_EXPONENT + 1):
        structural_discovery = []
        for index in range(DISCOVERY_DRAWS_PER_EXPONENT):
            parameter, multiplier = refinement_draw_parameter(
                target, anchor, prime, exponent, "discovery", index
            )
            record = audit_structure(family, parameter, prime, target_reduction)
            record["sample_index"] = index
            record["cylinder_multiplier"] = multiplier
            structural_discovery.append(record)
        structural_passes = all(
            row["status"] == "PASS_VALID_EXACT_SPECIALIZATION_AND_REDUCTION_AUDIT"
            and row["matches_anchor_target_reduction"]
            for row in structural_discovery
        )
        exponent_record = {
            "exponent": exponent,
            "modulus": str(ZZ(prime) ** exponent),
            "target_residue": p_adic_residue(target, prime, exponent),
            "structural_discovery": structural_discovery,
            "all_structural_discovery_draws_match": structural_passes,
        }
        exponents.append(exponent_record)
        if not structural_passes:
            continue
        discovery = []
        for index, structural in enumerate(structural_discovery):
            parameter = QQ(structural["parameter"])
            record = audit_parameter(family, parameter, (prime,), {str(prime): target_hash})
            compact = compact_outcome(record)
            compact["sample_index"] = index
            compact["cylinder_multiplier"] = structural["cylinder_multiplier"]
            compact["full_local_fingerprint_computed"] = True
            discovery.append(compact)
        discovery_passes = all(
            row["status"] == "PASS_VALID_EXACT_SPECIALIZATION_AND_LOCAL_AUDIT"
            and row["local_matches"].get(str(prime), False)
            for row in discovery
        )
        exponent_record["full_fingerprint_discovery"] = discovery
        exponent_record["all_full_fingerprint_discovery_draws_match"] = discovery_passes
        if not discovery_passes:
            continue
        confirmation = []
        for index in range(CONFIRMATION_DRAWS):
            parameter, multiplier = refinement_draw_parameter(
                target, anchor, prime, exponent, "confirmation", index
            )
            record = audit_parameter(family, parameter, (prime,), {str(prime): target_hash})
            record["sample_index"] = index
            record["cylinder_multiplier"] = multiplier
            confirmation.append(record)
        confirmation_passes = all(
            row["status"] == "PASS_VALID_EXACT_SPECIALIZATION_AND_LOCAL_AUDIT"
            and row["local_fingerprints"][0]["matches_anchor_target"]
            for row in confirmation
        )
        exponent_record["confirmation_count"] = len(confirmation)
        exponent_record["all_confirmation_draws_match"] = confirmation_passes
        exponent_record["confirmation_sha256"] = canonical_hash(
            [compact_outcome(row) for row in confirmation]
        )
        if confirmation_passes:
            selected = exponent
            selected_confirmation = confirmation
            break
    if selected is None:
        return {
            "prime": prime,
            "status": "NO_EMPIRICALLY_STABLE_EXPONENT_THROUGH_BOUND",
            "tested_exponents": exponents,
            "selected_exponent": None,
            "selected_confirmation_records": [],
        }
    return {
        "prime": prime,
        "status": "EMPIRICALLY_STABLE_ON_DISCOVERY_AND_CONFIRMATION_SAMPLES",
        "tested_exponents": exponents,
        "selected_exponent": selected,
        "selected_residue": p_adic_residue(target, prime, selected),
        "selected_modulus": str(ZZ(prime) ** selected),
        "selected_confirmation_records": [
            compact_outcome(row) for row in selected_confirmation
        ],
        "claim_boundary": "Sample stability is not a proof of local constancy on the entire p-adic cylinder.",
    }


def build():
    family = Family()
    prototype_artifact = json.loads(PROTOTYPES.read_text())
    prototypes = {
        int(row["prototype_curve_id"]): row
        for row in prototype_artifact["crt_search_prototypes"]["prototypes"]
    }
    if set(prototypes) != set(ANCHORS):
        raise ArithmeticError("the exact CRT prototype pair changed")

    target_records = {}
    target_hashes = {}
    for anchor, primes in ANCHORS.items():
        target = family.target_parameters[anchor]
        record = audit_parameter(family, target, primes, {str(prime): "TARGET" for prime in primes})
        if record["status"] != "PASS_VALID_EXACT_SPECIALIZATION_AND_LOCAL_AUDIT":
            raise ArithmeticError(f"anchor {anchor} failed the independent generic-MW17 audit")
        hashes = {
            str(row["rational_prime"]): row["comparison_sha256"]
            for row in record["local_fingerprints"]
        }
        for row in record["local_fingerprints"]:
            row["matches_anchor_target"] = True
        target_records[str(anchor)] = record
        target_hashes[str(anchor)] = hashes

    initial = {}
    refinements = {}
    cylinder_rows = []
    for anchor, primes in ANCHORS.items():
        print(f"{PROTOCOL}|anchor={anchor}|stage=initial-p3|samples={INITIAL_SAMPLE_SIZE}", flush=True)
        initial_records = initial_sample(
            family, anchor, prototypes[anchor], target_hashes[str(anchor)]
        )
        initial[str(anchor)] = initial_records
        prime_refinements = []
        for prime in primes:
            print(f"{PROTOCOL}|anchor={anchor}|prime={prime}|stage=refinement", flush=True)
            result = refine_prime(
                family,
                anchor,
                prime,
                family.target_parameters[anchor],
                target_hashes[str(anchor)][str(prime)],
                target_records[str(anchor)]["local_fingerprints"][
                    list(primes).index(prime)
                ]["comparison_payload"]["reduction"],
            )
            prime_refinements.append(result)
        refinements[str(anchor)] = prime_refinements
        if any(row["selected_exponent"] is None for row in prime_refinements):
            continue
        residues = [ZZ(row["selected_residue"]) for row in prime_refinements]
        moduli = [ZZ(row["selected_modulus"]) for row in prime_refinements]
        combined = ZZ(crt(residues, moduli))
        modulus = ZZ.prod(moduli)
        if combined > modulus // 2:
            combined -= modulus
        cylinder_rows.append(
            {
                "anchor_curve_id": anchor,
                "chart": "norm12-orbit-074d9",
                "prime_power_conditions": [
                    {
                        "prime": row["prime"],
                        "exponent": row["selected_exponent"],
                        "residue": row["selected_residue"],
                        "modulus": row["selected_modulus"],
                        "target_local_fingerprint_sha256": target_hashes[str(anchor)][str(row["prime"])],
                    }
                    for row in prime_refinements
                ],
                "combined_integer_parameter_class": {
                    "residue": str(combined),
                    "modulus": str(modulus),
                    "parameter_family": f"t = {combined} + {modulus}*n",
                },
                "empirical_only": True,
            }
        )

    frozen_definition = {
        "schema": "elkies-k3.r17-prospective-crt-frozen-cylinder-definition.v1",
        "selection_algorithm": (
            "For each anchor/place, scan k=3..20. Select the first k for which all 16 salted "
            "local-only discovery draws and all 64 disjoint salted confirmation draws exactly "
            "match the anchor's independently recomputed generic-MW17 local fingerprint."
        ),
        "selection_uses_point_search_outcomes": False,
        "salt": SALT,
        "cylinders": cylinder_rows,
        "claim_boundary": (
            "The congruences are exact and the stored samples match. Constancy of the Kummer/local "
            "fingerprint on either entire residue cylinder is not proved."
        ),
    }
    frozen_definition_hash = canonical_hash(frozen_definition)

    initial_summary = {}
    for anchor, rows in initial.items():
        counts = Counter()
        for row in rows:
            if row["status"] != "PASS_VALID_EXACT_SPECIALIZATION_AND_LOCAL_AUDIT":
                counts["structural_or_local_failure"] += 1
                continue
            matched = sum(local["matches_anchor_target"] for local in row["local_fingerprints"])
            counts[f"matched_{matched}_of_{len(row['local_fingerprints'])}_places"] += 1
        initial_summary[anchor] = dict(sorted(counts.items()))

    return {
        "schema": SCHEMA,
        "status": (
            "FROZEN_EMPIRICALLY_REFINED_LOCAL_CYLINDERS"
            if len(cylinder_rows) == len(ANCHORS)
            else "LOCAL_REFINEMENT_INCOMPLETE"
        ),
        "phase_boundary": {
            "phase": 1,
            "point_search_outcomes_read": False,
            "nagao_scores_read": False,
            "public_rank_used_for_candidate_selection": False,
            "historical_exceptional_points_used_in_prospective_fingerprint": False,
        },
        "predeclared_structural_conditions": [
            "exact nonzero specialized discriminant",
            "all seventeen saturated generic R17 sections satisfy the specialized equation",
            "exact isomorphism to an integral p-minimal model at each audited place",
            "completed-square 2-division cubic has discriminant 256 times the curve discriminant",
            "every generic MW17 Kummer norm is the expected exact square",
        ],
        "prospective_fingerprint_definition": {
            "uses_only_generic_MW17_and_local_curve_data": True,
            "features": [
                "Kodaira symbol, conductor exponent, minimal-discriminant valuation, Tamagawa number",
                "local factor degrees of the completed-square 2-division algebra",
                "ambient E(Q_p)/2E(Q_p) dimension",
                "actual generic-MW17 local squareclass image rows and source kernel",
                "generic-section component orders",
                "explicit NOT_COMPUTED marker for the resource-bounded Phase-3 Hilbert/Tate panel",
            ],
            "does_not_compute": [
                "complete 2-Selmer group",
                "Selmer upper bound",
                "Mordell-Weil rank upper bound",
                "exceptional point directions on prospective fibres",
            ],
        },
        "anchors": target_records,
        "original_p3_prototypes": {
            "source": relative(PROTOTYPES),
            "initial_sample_size_per_anchor": INITIAL_SAMPLE_SIZE,
            "records": initial,
            "summary": initial_summary,
        },
        "local_outcome_only_refinement": {
            "minimum_exponent": MIN_EXPONENT,
            "maximum_exponent": MAX_EXPONENT,
            "discovery_draws_per_exponent": DISCOVERY_DRAWS_PER_EXPONENT,
            "confirmation_draws_at_first_discovery_pass": CONFIRMATION_DRAWS,
            "results": refinements,
        },
        "frozen_cylinder_definition": frozen_definition,
        "frozen_cylinder_definition_sha256": frozen_definition_hash,
        "claim_boundary": [
            "A fingerprint match is an exact equality of the stored local data and subspaces, not a Selmer or rank statement.",
            "A sampled mismatch disproves constancy at that sampled point; sampled agreement does not prove cylinder-wide constancy.",
            "The refinement reads only local fingerprints and is frozen before any Mordell-Weil point-search outcome.",
            "No bounded miss is generated or interpreted in Phase 1.",
        ],
        "inputs": {
            relative(LINEAGE): digest(LINEAGE),
            relative(PROTOTYPES): digest(PROTOTYPES),
        },
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "command": "sage -python elkies-k3/scripts/audit_r17_prospective_crt_local_stability.sage",
        },
        "software_assumptions": {
            "sage": str(sage_version),
            "pari": ".".join(str(part) for part in pari.version()),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    document = build()
    serialized = json.dumps(document, indent=2, sort_keys=True) + "\n"
    output = args.output.resolve()
    if args.check:
        if not output.exists():
            raise ArithmeticError("stored prospective CRT local-stability audit is missing")
        stored = json.loads(output.read_text())
        if stored != document:
            # The upstream fingerprint artifact was extended after this
            # experiment was frozen.  Preserve its original whole-file hash
            # rather than silently refreezing after outcomes, while requiring
            # every mathematical field recomputed from the current extension
            # to agree byte-for-byte.  This check also accounts for this
            # check-only compatibility branch changing the generator hash.
            prototype_key = relative(PROTOTYPES)
            if (
                stored["inputs"].get(prototype_key) != FROZEN_PROTOTYPE_WHOLE_FILE_SHA256
                or stored["generation"].get("script_sha256") != FROZEN_GENERATOR_SHA256
            ):
                raise ArithmeticError("the frozen Phase-1 provenance is not the reviewed provenance")
            comparable = json.loads(serialized)
            comparable["inputs"][prototype_key] = stored["inputs"][prototype_key]
            comparable["generation"]["script_sha256"] = stored["generation"]["script_sha256"]
            if comparable != stored:
                raise ArithmeticError("stored prospective CRT mathematical payload differs from replay")
            print(
                f"{PROTOCOL}|warning=FROZEN_UPSTREAM_WHOLE_FILE_HASH_DRIFT|"
                "mathematical_payload=IDENTICAL",
                flush=True,
            )
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        f"{PROTOCOL}|status={document['status']}|"
        f"cylinder_hash={document['frozen_cylinder_definition_sha256']}|"
        f"output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
