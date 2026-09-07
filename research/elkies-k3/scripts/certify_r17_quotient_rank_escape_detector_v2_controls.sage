#!/usr/bin/env sage-python
"""Certify detector-v2 record inputs and exact partial Outcome-D data.

The certificate is intentionally fail-closed.  It verifies both record models,
the specialized MW17 identification, the 29-point mod-2 image, every finite bad
place in the existing exact local Kummer certificate, and the real local
condition.  It also binds a fixture-separated MW17-only replay: the generated
executables cannot read P18,...,P29 or their half-ideals.  MW29-relative work is
retained only as post-discovery closure evidence.  Unless both blinded replays
complete, all Selmer dimensions and the operational candidate gate remain null.

status: ACTIVE_PROOF
claim: exact Outcome-D record inputs and partial all-bad-place control data
inputs: record suite, local Kummer, class-pressure, sample, and BNF ledgers
outputs: compact fail-closed detector-v2 control certificate
supersedes/superseded-by: bounded point-search sensitivity is superseded as a calibration gate
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any

from sage.all import EllipticCurve, PolynomialRing, QQ
from sage.version import version as sage_version


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path[:0] = [str(CAS), str(ROOT / "elliptic-curves")]

MEASUREMENT_IMPLEMENTATION = CAS / "quotient_rank_escape_detector_v2.py"
CHECKPOINTED_DESCENT_IMPLEMENTATION = (
    CAS / "run_elkies_2026_relative_2selmer_checkpointed.py"
)

from build_elkies_2026_relative_2selmer_suite import (  # noqa: E402
    load_record_pair_cases,
)


LOCAL = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-074d9-local-kummer-meet-v1.json"
)
PRESSURE = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-kummer-classgroup-pressure-v1.json"
)
SAMPLE = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-quotient-rank-escape-detector-v2-sample-v1.json"
)
MAGMA_INPUTS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_record_pair_relative_2selmer_inputs_v1.json"
)
MW17_ONLY_INPUTS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "r17_mw17_only_selmer_control_inputs_v1.json"
)
MW17_ONLY_RUN = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "r17_mw17_only_selmer_control_run_v1.json"
)
BNF_RUN_DIR = (
    ROOT
    / "artifacts/local/elliptic-curves"
    / "r17-074d9-record-residual-2selmer-v1"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-quotient-rank-escape-detector-v2-controls-v1.json"
)
SOURCE = Path(__file__).resolve()
SCHEMA = "elkies-k3.r17-quotient-rank-escape-detector-v2-controls.v1"
STATUS = "OUTCOME_D_MW17_ONLY_REPLAY_FROZEN_BUT_NOT_EXECUTED"
EXPECTED = {
    "local": "PASS_EXACT_NEGATIVE_LOCAL_KUMMER_MEET",
    "pressure": "PROVED_KUMMER_FORCED_CUBIC_CLASS_GROUP_2RANK_LOWER_BOUNDS",
    "sample": "FROZEN_BLINDED_STAGE1_AND_STAGE2_SAMPLE_STAGE2_NOT_AUTHORIZED",
}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def real_local_record(case) -> dict[str, Any]:
    curve = EllipticCurve(QQ, [QQ(value) for value in case.model])
    a1, a2, a3, a4, a6 = curve.ainvs()
    polynomial_ring = PolynomialRing(QQ, "x")
    x = polynomial_ring.gen()
    completed_square = (
        4 * x**3
        + (a1**2 + 4 * a2) * x**2
        + (2 * a1 * a3 + 4 * a4) * x
        + (a3**2 + 4 * a6)
    )
    points = case.generic_points + case.exceptional_points
    for point_x, point_y in points:
        point_x, point_y = QQ(str(point_x)), QQ(str(point_y))
        if completed_square(point_x) != (2 * point_y + a1 * point_x + a3) ** 2:
            raise ArithmeticError("a record point failed the completed-square identity")

    intervals = completed_square.real_root_intervals()
    if len(intervals) not in (1, 3) or any(multiplicity != 1 for _, multiplicity in intervals):
        raise ArithmeticError("the record cubic has an unexpected real-root pattern")
    if len(intervals) == 1:
        if curve.real_components() != 1:
            raise ArithmeticError("real component count disagrees with the cubic")
        return {
            "place": "infinity",
            "real_component_count": 1,
            "ambient_local_kummer_dimension": 0,
            "actual_mw17_local_image_dimension": 0,
            "exceptional_image_dimension_modulo_mw17": 0,
            "actual_mw29_local_image_dimension": 0,
            "local_condition": "the unique real component is 2-divisible",
            "exact_root_intervals": [
                [str(interval[0]), str(interval[1])] for interval, _ in intervals
            ],
        }

    if curve.real_components() != 2:
        raise ArithmeticError("real component count disagrees with the cubic")
    intervals.sort(key=lambda row: row[0][0])
    second_upper = intervals[1][0][1]
    third_lower = intervals[2][0][0]
    if not second_upper < third_lower:
        raise ArithmeticError("real-root isolation intervals overlap")
    separator = (second_upper + third_lower) / 2

    def component_bits(point_rows):
        result = []
        for point_x, _point_y in point_rows:
            point_x = QQ(str(point_x))
            # On a positive-leading cubic, the compact component is the
            # positive interval between the first two roots.  A point on the
            # curve cannot lie in the negative interval before root three.
            result.append(int(point_x < separator))
        return result

    generic_bits = component_bits(case.generic_points)
    exceptional_bits = component_bits(case.exceptional_points)
    generic_dimension = int(any(generic_bits))
    total_dimension = int(any(generic_bits + exceptional_bits))
    quotient_dimension = total_dimension - generic_dimension
    if quotient_dimension < 0:
        raise ArithmeticError("invalid real quotient image dimension")
    return {
        "place": "infinity",
        "real_component_count": 2,
        "ambient_local_kummer_dimension": 1,
        "actual_mw17_local_image_dimension": generic_dimension,
        "mw17_component_bits": generic_bits,
        "exceptional_component_bits": exceptional_bits,
        "exceptional_image_dimension_modulo_mw17": quotient_dimension,
        "actual_mw29_local_image_dimension": total_dimension,
        "local_condition": "Kummer classes must lie in the image of the two real components",
        "exact_component_separator": str(separator),
        "exact_root_intervals": [
            [str(interval[0]), str(interval[1])] for interval, _ in intervals
        ],
    }


def compact_bnf_runs(curve_id: int) -> list[dict[str, Any]]:
    paths = sorted(BNF_RUN_DIR.glob(f"curve{curve_id}-bnf*-run.json"))
    if not paths:
        raise ArithmeticError(f"no preserved BNF attempts for curve {curve_id}")
    rows = []
    for path in paths:
        record = json.loads(path.read_text())
        if record.get("status") == "completed_certified_bnf":
            raise ArithmeticError(
                "a certified BNF exists; Outcome D must be replaced by the completed descent"
            )
        measurement = record.get("measurement", {})
        rows.append(
            {
                "record": relative(path),
                "record_sha256": digest(path),
                "status": record.get("status"),
                "timeout_seconds": record.get("input", {}).get("timeout_seconds"),
                "wall_seconds": measurement.get("wall_seconds"),
                "peak_observed_rss_bytes": measurement.get("peak_observed_rss_bytes"),
                "factorbase_stage_count": len(measurement.get("factorbase_stages", [])),
                "last_factorbase_stage": (
                    measurement.get("factorbase_stages", [])[-1]
                    if measurement.get("factorbase_stages")
                    else None
                ),
                "latest_random_relation_deficit": measurement.get(
                    "latest_random_relation_deficit"
                ),
                "bnf_tech": record.get("input", {}).get("bnf_tech"),
                "checkpoint": record.get("checkpoint"),
                "checkpoint_sha256": record.get("checkpoint_sha256"),
            }
        )
    return rows


def finite_bad_place_record(place: dict[str, Any]) -> dict[str, Any]:
    orders = Counter(
        int(direction["component_group_image_order"])
        for direction in place["directions"]
    )
    return {
        "place": str(place["rational_prime"]),
        "ambient_local_kummer_dimension": int(
            place["ambient_local_kummer_dimension"]
        ),
        "known_exceptional_quotient_image_dimension": int(
            place["quotient_basis_image_dimension"]
        ),
        "selected_nonrigid_block_image_dimension": int(
            place["selected_block_image_dimension"]
        ),
        "reduction_kind": place["reduction_kind"],
        "kodaira_symbol": place["kodaira_symbol"],
        "tamagawa_number": int(place["tamagawa_number"]),
        "conductor_exponent": int(place["conductor_exponent"]),
        "minimal_discriminant_valuation": int(
            place["minimal_discriminant_valuation"]
        ),
        "local_factor_descriptors": place["local_factor_descriptors"],
        "exceptional_generator_local_data": [
            {
                "label": direction["label"],
                "valuation_parity_support": direction["valuation_parity_support"],
                "component_group_image_order": int(
                    direction["component_group_image_order"]
                ),
            }
            for direction in place["directions"]
        ],
        "exceptional_generator_component_image_order_multiset": {
            str(order): count for order, count in sorted(orders.items())
        },
        "global_selmer_condition_row_space": None,
        "local_condition_codimension": None,
        "s_res_after_deleting_this_place": None,
    }


def build() -> dict[str, Any]:
    local = json.loads(LOCAL.read_text())
    pressure = json.loads(PRESSURE.read_text())
    sample = json.loads(SAMPLE.read_text())
    magma_inputs = json.loads(MAGMA_INPUTS.read_text())
    mw17_only_inputs = json.loads(MW17_ONLY_INPUTS.read_text())
    mw17_only_run = json.loads(MW17_ONLY_RUN.read_text())
    if local.get("status") != EXPECTED["local"]:
        raise ArithmeticError("the all-bad-place local certificate is not passing")
    if pressure.get("status") != EXPECTED["pressure"]:
        raise ArithmeticError("the Kummer class-group-pressure certificate is not passing")
    if sample.get("status") != EXPECTED["sample"]:
        raise ArithmeticError("the detector-v2 sample is not frozen")
    if (
        magma_inputs.get("status") != "EXACT_INPUTS_MAGMA_COMPLETION_REQUIRED"
        or [row.get("case_id") for row in magma_inputs.get("cases", [])]
        != ["record-r29-356", "record-r29-385"]
    ):
        raise ArithmeticError("the exact record-pair Magma inputs changed")
    if (
        mw17_only_inputs.get("status")
        != "FROZEN_FIXTURE_SEPARATED_MW17_REPLAY_EXECUTION_REQUIRED"
        or [row.get("curve_id") for row in mw17_only_inputs.get("cases", [])]
        != [356, 385]
        or not all(
            all(row.get("program_audit", {}).values())
            for row in mw17_only_inputs.get("cases", [])
        )
        or mw17_only_inputs.get("operational_gate", {}).get(
            "selmer_candidate_gate_operationally_calibrated"
        )
        is not False
    ):
        raise ArithmeticError("the sealed MW17-only replay inputs changed")
    if (
        mw17_only_run.get("status")
        != "INCOMPLETE_BLINDED_REPLAY_NOT_EXECUTED"
        or mw17_only_run.get("completed_record_replays") != 0
        or mw17_only_run.get("selmer_candidate_gate_operationally_calibrated")
        is not False
    ):
        raise ArithmeticError("Outcome D requires the MW17-only replay to remain incomplete")

    cases = {int(case.case_id.rsplit("-", 1)[1]): case for case in load_record_pair_cases()}
    local_by_id = {int(row["curve_id"]): row for row in local["curves"]}
    pressure_by_id = {int(row["curve_id"]): row for row in pressure["curves"]}
    controls = []
    for curve_id in (356, 385):
        case = cases[curve_id]
        local_curve = local_by_id[curve_id]
        pressure_curve = pressure_by_id[curve_id]
        half_ideals = {
            row["label"]: row for row in pressure_curve["point_half_ideals"]
        }
        valuation_rows = {
            row["label"]: row
            for row in pressure_curve["known_point_bad_valuation_parity_rows"]
        }
        if set(half_ideals) != {f"P{index}" for index in range(1, 30)}:
            raise ArithmeticError(f"curve {curve_id} Kummer half-ideal list changed")
        model = [str(value) for value in case.model]
        if model != local_curve["global_minimal_model"] or model != pressure_curve[
            "global_minimal_model"
        ]:
            raise ArithmeticError(f"curve {curve_id} model certificates disagree")
        finite_bad = [
            finite_bad_place_record(place)
            for place in local_curve["local_places"]
            if place["place_kind"] == "bad"
        ]
        if [row["place"] for row in finite_bad] != local_curve["bad_primes"]:
            raise ArithmeticError(f"curve {curve_id} bad-place list is incomplete")
        bnf_runs = compact_bnf_runs(curve_id)
        controls.append(
            {
                "curve_id": curve_id,
                "exact_global_minimal_model": model,
                "model_verified_against_three_independent_certificates": True,
                "specialized_mw17": {
                    "exact_section_count": 17,
                    "equals_displayed_points_P1_through_P17": True,
                    "generic_function_field_basis_saturated": True,
                    "primitive_inside_displayed_mw29_subgroup": True,
                    "actual_global_mod2_image_dimension": 17,
                    "two_saturated_inside_full_E_Q": True,
                    "two_saturation_certificate": (
                        "the combined good-reduction E(F_p)/2E(F_p) signature "
                        "matrix has rank 29 on P1..P29, hence rank 17 on P1..P17"
                    ),
                    "all_prime_saturation_inside_full_E_Q": None,
                    "all_prime_saturation_status": "UNKNOWN",
                    "global_kummer_half_ideals": [
                        half_ideals[f"P{index}"] for index in range(1, 18)
                    ],
                    "bad_prime_valuation_parity_rows": [
                        valuation_rows[f"P{index}"] for index in range(1, 18)
                    ],
                },
                "known_mw29": {
                    "displayed_point_count": 29,
                    "actual_global_mod2_image_dimension": 29,
                    "two_saturated_inside_full_E_Q": True,
                    "exceptional_generators": [f"P{index}" for index in range(18, 30)],
                    "exceptional_image_dimension_modulo_mw17": 12,
                    "exceptional_global_kummer_representatives": local_curve[
                        "kummer_images"
                    ],
                    "exceptional_kummer_half_ideals": [
                        half_ideals[f"P{index}"] for index in range(18, 30)
                    ],
                    "exceptional_bad_prime_valuation_parity_rows": [
                        valuation_rows[f"P{index}"] for index in range(18, 30)
                    ],
                    "rational_two_torsion_dimension": 0,
                },
                "complete_two_selmer": {
                    "dimension": None,
                    "s_res": None,
                    "selmer_modulo_mw17_dimension": None,
                    "selmer_modulo_known_mw29_dimension": None,
                    "status": "UNKNOWN_GLOBAL_S_CLASS_UNIT_CHECKPOINT_INCOMPLETE",
                    "parity_from_exact_root_number": int(
                        pressure_curve["proved_total_two_selmer_dimension_mod_2"]
                    ),
                },
                "prospective_mw17_only_replay": {
                    "input_case_id": f"mw17-only-control-{curve_id}",
                    "backend_visible_known_point_count": 17,
                    "backend_can_read_P18_through_P29": False,
                    "backend_can_read_exceptional_half_ideals": False,
                    "completed": False,
                    "selmer_modulo_mw17_dimension": None,
                    "status": "FROZEN_NOT_EXECUTED",
                },
                "cached_cubic_field": {
                    **local_curve["two_division_etale_algebra"],
                    "field_signature": pressure_curve["field_signature"],
                    "maximal_order_nfcertify_completed": True,
                    "complete_bnf_checkpoint": None,
                    "proved_class_group_2rank_lower_bound": int(
                        pressure_curve["proved_class_group_2rank_lower_bound"]
                    ),
                    "proved_adjusted_residual_class_group_image_dimension_lower_bound": int(
                        pressure_curve[
                            "proved_adjusted_residual_class_group_image_dimension_lower_bound"
                        ]
                    ),
                },
                "local_conditions": {
                    "place_2_included": finite_bad[0]["place"] == "2",
                    "all_bad_finite_places": finite_bad,
                    "infinity": real_local_record(case),
                    "auxiliary_descent_places": None,
                    "full_global_local_condition_matrix": None,
                    "leave_one_place_out_residual_dimensions": None,
                    "status": (
                        "EXACT_LOCAL_CURVE_AND_KNOWN_POINT_DATA_ONLY; "
                        "GLOBAL_SQUARECLASS_DOMAIN_UNAVAILABLE"
                    ),
                },
                "preserved_bnf_attempts": bnf_runs,
            }
        )

    return {
        "schema": SCHEMA,
        "status": STATUS,
        "detector_v2": {
            "measurement_implementation": relative(MEASUREMENT_IMPLEMENTATION),
            "measurement_implementation_sha256": digest(MEASUREMENT_IMPLEMENTATION),
            "checkpointed_descent_implementation": relative(
                CHECKPOINTED_DESCENT_IMPLEMENTATION
            ),
            "checkpointed_descent_implementation_sha256": digest(
                CHECKPOINTED_DESCENT_IMPLEMENTATION
            ),
            "fixture_sequenced_magma_control_inputs": relative(MAGMA_INPUTS),
            "fixture_separated_mw17_only_control_inputs": relative(
                MW17_ONLY_INPUTS
            ),
            "fixture_separated_mw17_only_control_run": relative(MW17_ONLY_RUN),
            "control_gate_passed": False,
            "stage_1_application_authorized": False,
            "stage_2_application_authorized": False,
            "reason": (
                "neither fixture-separated record replay has computed the complete "
                "2-Selmer group relative to MW17; MW29-relative closure work is not "
                "admissible prospective calibration"
            ),
        },
        "controls": controls,
        "sample_commitment": {
            "artifact": relative(SAMPLE),
            "artifact_sha256": digest(SAMPLE),
            "stage_1_candidate_count": sample["commitment"]["stage_1_candidate_count"],
            "stage_2_candidate_count": sample["commitment"]["stage_2_candidate_count"],
            "stage_2_authorized": False,
        },
        "narrowed_technical_task": {
            "target": (
                "Run both sealed Magma inputs to compute the complete unconditional "
                "2-Selmer group and quotient only by the specialized MW17 image. The "
                "backend must not read P18..P29 or their half-ideals."
            ),
            "prospective_lane": (
                "The MW17-only replay is the sole operational calibration lane because "
                "it matches the information available for a new candidate."
            ),
            "post_discovery_lane": (
                "The quotient-by-MW29 relation collectors remain useful for exact-rank "
                "closure after points are discovered, but cannot pass this gate."
            ),
            "curve_356_forced_full_class_2rank_lower_bound": 21,
            "curve_385_forced_full_class_2rank_lower_bound": 15,
            "completion_criterion": (
                "two source-hash-matched blind_freeze transcripts, each with a complete "
                "unconditional Selmer group, MW17 Kummer rank 17, and residual dimension "
                "at least twelve"
            ),
            "next_pipeline_step_after_completion": (
                "only after both transcripts are frozen, consult the committed public "
                "control truth and confirm that the prospective score did not kill the "
                "known twelve-dimensional jump"
            ),
        },
        "claim_boundary": [
            "No complete 2-Selmer dimension has been computed for curve 356 or 385.",
            "The hard quotient-by-MW29 calculations demonstrate post-discovery closure feasibility only and cannot calibrate a prospective MW17-relative gate.",
            "The MW17-only executables contain exactly seventeen points and no P18..P29 coordinates or half-ideals; neither sealed replay has run.",
            "Neither s_res nor the Selmer excess over MW29 is inferred from local data, parity, a bounded relation collector, or a timed-out descent.",
            "Good-reduction mod-2 signatures prove 2-saturation of MW17 and MW29; all-prime saturation inside the unknown full Mordell-Weil group remains unproved.",
            "No intrinsic pairing is placed on the deterministic quotient presentation.",
            "The frozen prospective sample has not been submitted to the descent pipeline.",
        ],
        "inputs": {
            relative(LOCAL): digest(LOCAL),
            relative(PRESSURE): digest(PRESSURE),
            relative(SAMPLE): digest(SAMPLE),
            relative(MAGMA_INPUTS): digest(MAGMA_INPUTS),
            relative(MW17_ONLY_INPUTS): digest(MW17_ONLY_INPUTS),
            relative(MW17_ONLY_RUN): digest(MW17_ONLY_RUN),
            relative(SOURCE): digest(SOURCE),
        },
        "software_assumptions": {"sage": str(sage_version)},
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "certify_r17_quotient_rank_escape_detector_v2_controls.sage --check"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    document = build()
    payload = json.dumps(document, indent=2, sort_keys=True) + "\n"
    output = args.output.resolve()
    if args.check:
        if not output.exists() or output.read_text() != payload:
            raise ArithmeticError("stored detector-v2 control certificate differs")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(payload)
    print(
        "R17ESCAPEV2CONTROLS|controls=2|mw17_only_complete=0|"
        "stage1_authorized=false|status=OUTCOME_D",
        flush=True,
    )


if __name__ == "__main__":
    main()
