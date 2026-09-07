#!/usr/bin/env python3
"""Freeze matched R17 CRT, ablation, and ordinary-control cohorts.

This is a label-blind Phase-2 commitment.  It reads the frozen Phase-1 local
cylinder definition, constructs exact integer parameters with tightly matched
heights, rejects only singular or duplicate specializations, and leaves every
arithmetic feature and Mordell--Weil outcome unopened.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PHASE1 = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-local-stability-v1.json"
LINEAGE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-frozen-cohorts-v1.json"

SCHEMA = "elkies-k3.r17-prospective-crt-frozen-cohorts.v1"
SALT = "r17-prospective-crt-frozen-cohorts-v1"
ROWS_PER_ANCHOR = 256
MINIMUM_ABSOLUTE_PARAMETER = 2**110
MAXIMUM_ABSOLUTE_PARAMETER = 2**111 - 1
EXPECTED_PHASE1_CYLINDER_HASH = "500dc6931c5aeaf3d6d9982bb994286d7aee36e7c87b9e414e8b7e0ef8aef15c"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def polynomial_value(coefficients: list[Fraction], parameter: int) -> Fraction:
    answer = Fraction(0)
    for coefficient in reversed(coefficients):
        answer = answer * parameter + coefficient
    return answer


def bit_length_fraction(value: Fraction) -> dict[str, int]:
    return {
        "numerator_bits": abs(value.numerator).bit_length(),
        "denominator_bits": value.denominator.bit_length(),
    }


def crt_pair(left_residue: int, left_modulus: int, right_residue: int, right_modulus: int):
    if __import__("math").gcd(left_modulus, right_modulus) != 1:
        raise ArithmeticError("CRT moduli are not coprime")
    correction = ((right_residue - left_residue) * pow(left_modulus, -1, right_modulus)) % right_modulus
    modulus = left_modulus * right_modulus
    return (left_residue + left_modulus * correction) % modulus, modulus


def combine_conditions(conditions):
    residue, modulus = 0, 1
    for row in conditions:
        residue, modulus = crt_pair(
            residue,
            modulus,
            int(row["residue"]),
            int(row["modulus"]),
        )
    if residue > modulus // 2:
        residue -= modulus
    return residue, modulus


def hashed_block(tag: str, index: int) -> bytes:
    return sha256(f"{SALT}|{tag}|{index}".encode()).digest()


def target_height(anchor: int, index: int) -> int:
    block = hashed_block(f"height|{anchor}", index)
    magnitude = MINIMUM_ABSOLUTE_PARAMETER + int.from_bytes(block, "big") % (
        MAXIMUM_ABSOLUTE_PARAMETER - MINIMUM_ABSOLUTE_PARAMETER + 1
    )
    return magnitude if block[0] & 1 else -magnitude


def nearest_congruent(target: int, residue: int, modulus: int, predicate=lambda value: True):
    quotient, remainder = divmod(target - residue, modulus)
    if 2 * remainder > modulus:
        quotient += 1
    offsets = [0]
    for step in range(1, 1000):
        offsets.extend((step, -step))
    for offset in offsets:
        candidate = residue + modulus * (quotient + offset)
        if (
            MINIMUM_ABSOLUTE_PARAMETER <= abs(candidate) <= MAXIMUM_ABSOLUTE_PARAMETER
            and predicate(candidate)
        ):
            return candidate, quotient + offset
    raise ArithmeticError("no admissible nearby member of the residue class")


def random_residue_conditions(anchor: int, conditions):
    rows = []
    for row in conditions:
        modulus = int(row["modulus"])
        target = int(row["residue"])
        counter = 0
        while True:
            residue = int.from_bytes(
                hashed_block(f"random-residue|{anchor}|{row['prime']}", counter), "big"
            ) % modulus
            counter += 1
            if residue != target:
                break
        rows.append(
            {
                "prime": int(row["prime"]),
                "exponent": int(row["exponent"]),
                "residue": residue,
                "modulus": modulus,
                "selection_counter": counter - 1,
            }
        )
    return rows


def build(*, rows_per_anchor: int, salt: str):
    global SALT
    SALT = salt
    if rows_per_anchor < 200:
        raise ValueError("the initial commitment requires hundreds of rows per anchor cohort")

    phase1 = json.loads(PHASE1.read_text())
    if phase1.get("status") != "FROZEN_EMPIRICALLY_REFINED_LOCAL_CYLINDERS":
        raise ArithmeticError("Phase 1 has not frozen two empirically refined cylinders")
    if phase1.get("frozen_cylinder_definition_sha256") != EXPECTED_PHASE1_CYLINDER_HASH:
        raise ArithmeticError("the reviewed Phase-1 cylinder definition changed")
    if canonical_hash(phase1["frozen_cylinder_definition"]) != EXPECTED_PHASE1_CYLINDER_HASH:
        raise ArithmeticError("the embedded Phase-1 cylinder hash does not replay")
    cylinders = {
        int(row["anchor_curve_id"]): row
        for row in phase1["frozen_cylinder_definition"]["cylinders"]
    }
    if set(cylinders) != {356, 385}:
        raise ArithmeticError("the frozen cylinder pair changed")

    lineage = json.loads(LINEAGE.read_text())
    representative = lineage["representative"]
    if representative["chart"] != "norm12-orbit-074d9":
        raise ArithmeticError("the exact representative chart changed")
    a_coefficients = [Fraction(value) for value in representative["A_coefficients_low_to_high"]]
    b_coefficients = [Fraction(value) for value in representative["B_coefficients_low_to_high"]]

    definitions = {}
    for anchor, cylinder in sorted(cylinders.items()):
        full_conditions = cylinder["prime_power_conditions"]
        two_conditions = [row for row in full_conditions if int(row["prime"]) == 2]
        odd_conditions = [row for row in full_conditions if int(row["prime"]) != 2]
        if len(two_conditions) != 1 or len(odd_conditions) != 4:
            raise ArithmeticError("the expected 2-only/odd-only split changed")
        full_residue, full_modulus = combine_conditions(full_conditions)
        two_residue, two_modulus = combine_conditions(two_conditions)
        odd_residue, odd_modulus = combine_conditions(odd_conditions)
        random_conditions = random_residue_conditions(anchor, full_conditions)
        random_residue, random_modulus = combine_conditions(random_conditions)
        if random_modulus != full_modulus:
            raise ArithmeticError("the random-residue ablation changed codimension")
        definitions[anchor] = {
            "full": (full_residue, full_modulus),
            "two_only": (two_residue, two_modulus),
            "odd_only": (odd_residue, odd_modulus),
            "random_equal_codimension": (random_residue, random_modulus),
            "full_conditions": full_conditions,
            "random_conditions": random_conditions,
        }

    def matches(value: int, residue: int, modulus: int) -> bool:
        return (value - residue) % modulus == 0

    seen_parameters = set()
    seen_j = set()
    rows = []
    structural_failures = []

    def add_candidate(anchor: int, index: int, cohort: str, target: int, residue: int, modulus: int, predicate):
        parameter, multiplier = nearest_congruent(target, residue, modulus, predicate)
        if parameter in seen_parameters:
            raise ArithmeticError("duplicate exact parameter across frozen cohorts")
        coefficient_a = polynomial_value(a_coefficients, parameter)
        coefficient_b = polynomial_value(b_coefficients, parameter)
        discriminant_core = 4 * coefficient_a**3 + 27 * coefficient_b**2
        if discriminant_core == 0:
            structural_failures.append(
                {
                    "anchor": anchor,
                    "match_index": index,
                    "cohort": cohort,
                    "parameter": str(parameter),
                    "reason": "singular_specialization",
                }
            )
            return False
        j_invariant = Fraction(6912) * coefficient_a**3 / discriminant_core
        j_hash = sha256(f"{j_invariant.numerator}/{j_invariant.denominator}".encode()).hexdigest()
        if j_hash in seen_j:
            structural_failures.append(
                {
                    "anchor": anchor,
                    "match_index": index,
                    "cohort": cohort,
                    "parameter": str(parameter),
                    "reason": "duplicate_exact_j_invariant",
                }
            )
            return False
        seen_parameters.add(parameter)
        seen_j.add(j_hash)
        sample_key = {
            "anchor": anchor,
            "match_index": index,
            "cohort": cohort,
            "parameter": parameter,
        }
        rows.append(
            {
                "sample_id": canonical_hash(sample_key)[:24],
                "match_set_id": f"anchor-{anchor}-height-{index:04d}",
                "anchor_curve_id": anchor,
                "cohort": cohort,
                "chart": "norm12-orbit-074d9",
                "parameter": str(parameter),
                "projective_pair": [parameter, 1],
                "projective_height": abs(parameter),
                "matched_target_parameter": str(target),
                "absolute_height_error": abs(parameter - target),
                "cylinder_residue": str(residue),
                "cylinder_modulus": str(modulus),
                "cylinder_multiplier": multiplier,
                "j_invariant_sha256": j_hash,
                "equation_complexity": {
                    "A": bit_length_fraction(coefficient_a),
                    "B": bit_length_fraction(coefficient_b),
                    "discriminant_core": bit_length_fraction(discriminant_core),
                },
                "selection_status": "FROZEN_VALID_NONSINGULAR_UNOPENED",
                "arithmetic_features": None,
                "monotone_sieve": None,
                "point_search": None,
                "outcome_status": "NOT_OPENED",
            }
        )
        return True

    cohort_for_anchor = {356: "A_356_full", 385: "B_385_full"}
    for anchor in (356, 385):
        definition = definitions[anchor]
        full_residue, full_modulus = definition["full"]
        two_residue, two_modulus = definition["two_only"]
        odd_residue, odd_modulus = definition["odd_only"]
        random_residue, random_modulus = definition["random_equal_codimension"]
        for index in range(rows_per_anchor):
            target = target_height(anchor, index)
            lanes = (
                (
                    cohort_for_anchor[anchor],
                    full_residue,
                    full_modulus,
                    lambda value: True,
                ),
                (
                    "C_matched_ordinary",
                    0,
                    1,
                    lambda value, a=anchor: all(
                        not matches(value, int(row["residue"]), int(row["modulus"]))
                        for row in definitions[a]["full_conditions"]
                    ),
                ),
                (
                    "D_two_only",
                    two_residue,
                    two_modulus,
                    lambda value, a=anchor: all(
                        not matches(value, int(row["residue"]), int(row["modulus"]))
                        for row in definitions[a]["full_conditions"]
                        if int(row["prime"]) != 2
                    ),
                ),
                (
                    "E_odd_only",
                    odd_residue,
                    odd_modulus,
                    lambda value, a=anchor: not matches(
                        value,
                        int(next(row for row in definitions[a]["full_conditions"] if int(row["prime"]) == 2)["residue"]),
                        int(next(row for row in definitions[a]["full_conditions"] if int(row["prime"]) == 2)["modulus"]),
                    ),
                ),
                (
                    "F_random_equal_codimension",
                    random_residue,
                    random_modulus,
                    lambda value: True,
                ),
            )
            for cohort, residue, modulus, predicate in lanes:
                if not add_candidate(anchor, index, cohort, target, residue, modulus, predicate):
                    raise ArithmeticError(
                        "a predeclared structural rejection occurred; incrementing or replacing rows after freeze is forbidden"
                    )

    cohort_counts = Counter(row["cohort"] for row in rows)
    expected_counts = {
        "A_356_full": rows_per_anchor,
        "B_385_full": rows_per_anchor,
        "C_matched_ordinary": 2 * rows_per_anchor,
        "D_two_only": 2 * rows_per_anchor,
        "E_odd_only": 2 * rows_per_anchor,
        "F_random_equal_codimension": 2 * rows_per_anchor,
    }
    if cohort_counts != expected_counts:
        raise ArithmeticError("the balanced cohort counts changed")
    if structural_failures:
        raise ArithmeticError("the initial frozen cohort encountered a structural rejection")
    if any(row["arithmetic_features"] is not None or row["point_search"] is not None for row in rows):
        raise ArithmeticError("a pre-search label was opened during cohort commitment")

    cohort_hashes = {
        cohort: canonical_hash(
            [row for row in rows if row["cohort"] == cohort]
        )
        for cohort in expected_counts
    }
    candidate_list_payload = [
        {
            key: row[key]
            for key in (
                "sample_id",
                "match_set_id",
                "anchor_curve_id",
                "cohort",
                "parameter",
                "projective_pair",
                "cylinder_residue",
                "cylinder_modulus",
            )
        }
        for row in rows
    ]
    candidate_list_hash = canonical_hash(candidate_list_payload)

    complexity_by_cohort = defaultdict(lambda: defaultdict(list))
    for row in rows:
        cohort = row["cohort"]
        complexity_by_cohort[cohort]["parameter_bits"].append(row["projective_height"].bit_length())
        for key in ("A", "B", "discriminant_core"):
            complexity_by_cohort[cohort][f"{key}_numerator_bits"].append(
                row["equation_complexity"][key]["numerator_bits"]
            )
    matching_summary = {
        cohort: {
            key: {"minimum": min(values), "maximum": max(values)}
            for key, values in fields.items()
        }
        for cohort, fields in complexity_by_cohort.items()
    }

    definitions_json = {}
    for anchor, definition in definitions.items():
        definitions_json[str(anchor)] = {
            "full": {"residue": str(definition["full"][0]), "modulus": str(definition["full"][1])},
            "two_only": {"residue": str(definition["two_only"][0]), "modulus": str(definition["two_only"][1])},
            "odd_only": {"residue": str(definition["odd_only"][0]), "modulus": str(definition["odd_only"][1])},
            "random_equal_codimension": {
                "residue": str(definition["random_equal_codimension"][0]),
                "modulus": str(definition["random_equal_codimension"][1]),
                "prime_power_conditions": definition["random_conditions"],
            },
        }

    payload = {
        "schema": SCHEMA,
        "status": "FROZEN_UNOPENED_MATCHED_CRT_AND_ABLATION_COHORTS",
        "phase_boundary": {
            "phase": 2,
            "frozen_before_any_mordell_weil_search": True,
            "point_search_outcomes_read": False,
            "nagao_scores_used_for_selection": False,
            "public_rank_or_cover_splits_used_for_selection": False,
            "post_freeze_extension_forbidden": True,
        },
        "commitment": {
            "salt": salt,
            "phase1_frozen_cylinder_definition_sha256": EXPECTED_PHASE1_CYLINDER_HASH,
            "candidate_list_sha256": candidate_list_hash,
            "cohort_sha256": cohort_hashes,
            "rows_per_anchor_lane": rows_per_anchor,
            "scheduled_candidate_count": len(rows),
        },
        "cohort_definitions": {
            "A_356_full": "all five empirically refined 356-like local residue conditions",
            "B_385_full": "all five empirically refined 385-like local residue conditions",
            "C_matched_ordinary": "same chart and height targets; avoids every individual anchor residue condition",
            "D_two_only": "refined anchor p=2 condition; avoids every anchor odd-prime condition",
            "E_odd_only": "all four anchor odd-prime conditions; excludes the refined anchor p=2 condition",
            "F_random_equal_codimension": "salted random residues at the same prime powers as the corresponding full anchor",
            "exact_residue_definitions_by_anchor": definitions_json,
        },
        "matching": {
            "parameter_domain": "integer t=a/1 on the exact norm12-orbit-074d9 family",
            "absolute_height_shell": [MINIMUM_ABSOLUTE_PARAMETER, MAXIMUM_ABSOLUTE_PARAMETER],
            "algorithm": (
                "For each anchor/index draw one salted signed target height, then take the nearest "
                "admissible member of each cohort congruence. No arithmetic feature or outcome is used."
            ),
            "same_denominator": 1,
            "singular_exclusion": "exactly 4*A(t)^3+27*B(t)^2 != 0",
            "duplicate_exclusion": "exact parameter and exact j-invariant",
            "structural_failures": structural_failures,
            "complexity_ranges_by_cohort": matching_summary,
        },
        "frozen_feature_protocol": {
            "places": [2, 13, 37, 53, 67, 71],
            "per_candidate": [
                "exact specialization and all seventeen generic-section identities",
                "local Kummer dimensions and actual generic-MW17 localization presentations",
                "localization source kernels, cumulative intersection ranks, and leave-one-place-out ranks",
                "Tamagawa and component-image data",
                "resource-bounded factorwise Hilbert/Tate panel with explicit NOT_COMPUTED values",
                "the existing fixed three-block Nagao score as a comparison feature only",
                "monotone proved residual-Selmer upper-bound sequence, or explicit no-finite-bound-yet",
            ],
            "selection_or_rebalancing_from_features_forbidden": True,
        },
        "frozen_search_protocol": {
            "protocol_id": "r17-prospective-crt-uniform-mw-escape-v1",
            "engine": "Sage/eclib mwrank_MordellWeil.search",
            "search_height": 12,
            "maxr": 32,
            "pp": 0,
            "wall_clock_limit_seconds": 300,
            "memory_limit_bytes": 8_000_000_000,
            "retries": 0,
            "adaptive_stopping": False,
            "same_budget_for_every_row": True,
            "monotone_residual_selmer_gate_required": True,
            "open_incomplete_sieve_may_authorize_only_this_bounded_search": True,
            "positive_confirmation": [
                "exact specialized equation",
                "escape from the specialized saturated generic MW17 subgroup",
                "exact finite-quotient independence certificate before counting a direction",
            ],
            "outcome_labels": {
                "positive": "CERTIFIED_MW17_ESCAPE",
                "completed_miss": "BOUNDED_PROTOCOL_NO_ESCAPE_FOUND",
                "timeout": "CENSORED_TIMEOUT",
                "backend_failure": "CENSORED_BACKEND_FAILURE",
                "structural_failure": "CENSORED_STRUCTURAL_FAILURE",
            },
            "no_promotion_rules": [
                "Incomplete descent is never a Selmer upper bound.",
                "Selmer dimension is never identified with Mordell-Weil rank.",
                "A bounded miss is not rank 17 and not a zero residual quotient.",
                "Timeouts and backend failures are not bounded misses.",
            ],
        },
        "analysis_plan": {
            "primary_contrasts": [
                "A_356_full versus anchor-matched C_matched_ordinary",
                "B_385_full versus anchor-matched C_matched_ordinary",
                "pooled full (A+B) versus anchor-stratified ordinary controls",
            ],
            "ablation_order": [
                "F_random_equal_codimension",
                "D_two_only",
                "E_odd_only",
                "A_356_full/B_385_full",
            ],
            "report": [
                "exact denominators and cohort counts",
                "escape proportions and risk differences with exact binomial intervals",
                "odds ratios and Fisher exact tests",
                "certified direction-count and time/effort summaries",
                "censored timeout/backend counts kept separate",
            ],
            "underpowered_extension_rule": (
                "Any larger replication is a separately salted, separately frozen experiment; "
                "these candidate lists are never extended in place."
            ),
        },
        "rows": rows,
        "inputs": {
            relative(PHASE1): digest(PHASE1),
            relative(LINEAGE): digest(LINEAGE),
        },
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "command": "python3 elkies-k3/scripts/build_r17_prospective_crt_cohorts.py",
        },
        "claim_boundary": [
            "This artifact freezes candidate lists and protocol only; it contains no arithmetic feature or point-search outcome.",
            "The full cylinders are empirically locally stable on Phase-1 samples, not proved constant.",
            "Matched integer heights and equation-bit ranges do not imply matched conductors.",
            "No cohort is an unbiased model of the historical 69 public R17 successes.",
        ],
    }
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows-per-anchor", type=int, default=ROWS_PER_ANCHOR)
    parser.add_argument("--salt", default=SALT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build(rows_per_anchor=args.rows_per_anchor, salt=args.salt)
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = args.output.resolve()
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored prospective CRT cohort commitment differs from replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "R17CRTCOHORTS"
        f"|rows={payload['commitment']['scheduled_candidate_count']}"
        f"|candidate_hash={payload['commitment']['candidate_list_sha256']}"
        "|status=FROZEN_UNOPENED",
        flush=True,
    )


if __name__ == "__main__":
    main()
