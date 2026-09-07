#!/usr/bin/env python3
"""Build the transcendental-first queue for the rank-19 K3 foundry.

This planner deliberately does not inspect rootless frames while ordering the
arithmetic phase.  It begins with the rank-three T ledger, records exact full
marking decisions when available, and otherwise ranks only the next
discriminant-kernel / rational-point calculation.  A row reaches the separate
NS-complement and rootless-frame handoff only after an exact non-CM rational
point on the full marked curve has been certified.

Coarse Clifford-curve genus is retained as a heuristic when the stable curve
is unknown; it is never promoted to a full-marking genus or existence claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
T_ARITHMETIC = GENERATED / "elkies-k3-rank7-t-arithmetic-v1.json"
CLASSIFIER = GENERATED / "elkies-k3-rank19-arithmetic-marking-classifier-v1.json"
GLOBAL_DECISIONS = (
    ROOT
    / "elkies-k3/data/arithmetic/arithmetic-first-marking-decisions-v1.json"
)
OUTPUT = GENERATED / "elkies-k3-arithmetic-first-marked-t-foundry-v1.json"
H3_SURFACE_ID = "K3-8188cdcda8c57b2d"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resolve_path(payload: dict, dotted_path: str):
    value = payload
    for component in dotted_path.split("."):
        if not isinstance(value, dict) or component not in value:
            raise KeyError(f"missing assertion path {dotted_path}")
        value = value[component]
    return value


def validate_certificate(certificate: dict) -> dict:
    path = ROOT / certificate["path"]
    if not path.is_file():
        raise FileNotFoundError(certificate["path"])
    payload = json.loads(path.read_text())
    for assertion in certificate["assertions"]:
        actual = resolve_path(payload, assertion["path"])
        if actual != assertion["equals"]:
            raise AssertionError(
                f"{certificate['path']}:{assertion['path']} changed: "
                f"{actual!r} != {assertion['equals']!r}"
            )
    return {
        "path": certificate["path"],
        "sha256": digest(path),
        "assertions_replayed": len(certificate["assertions"]),
    }


def coarse_label(row: dict) -> str | None:
    base = row["arithmetic_source"]["base_curve"]
    return base.get("label") or base.get("group")


def full_curve_summary(classifier_row: dict | None) -> dict:
    if classifier_row is None:
        return {
            "status": "UNKNOWN_NO_ROOTLESS_SUBCATALOGUE_DECISION_RECORD",
            "label": None,
            "genus": None,
            "rational_non_CM_point": False,
            "stable_kernel_index_over_coarse": None,
        }
    curve = classifier_row["full_discriminant_marking_curve"]
    tests = classifier_row["arithmetic_tests"]
    has_point_witness = any(
        test["status"]
        in {
            "PASS_EXACT_NONCM_QQ_POINT",
            "PASS_EXACT_FULL_RATIONAL_MARKING_WITNESS",
        }
        for test in tests
    )
    classification = classifier_row["classification"]
    has_point = (
        True
        if has_point_witness
        else False
        if classification == "ARITHMETICALLY_EXCLUDED"
        else None
        if classifier_row.get("phase_2_certificate_status")
        == "UNRESOLVED_FOR_EXPLICIT_REASON"
        else False
    )
    kernel_indices = [
        quotient["degree"]
        for quotient in classifier_row["easy_quotient_maps"]
        if isinstance(quotient.get("degree"), int)
        and "full" in quotient.get("source", "").lower()
    ]
    direct_kernel_index = curve.get("degree_from_norm_one_curve")
    return {
        "status": curve["status"],
        "label": curve.get("label") or curve.get("explicit_model"),
        "genus": curve.get("genus"),
        "rational_non_CM_point": has_point,
        "stable_kernel_index_over_coarse": (
            int(direct_kernel_index)
            if isinstance(direct_kernel_index, int)
            else (min(kernel_indices) if kernel_indices else None)
        ),
    }


def curve_priority(row: dict, classifier_row: dict | None) -> tuple:
    full = full_curve_summary(classifier_row)
    coarse = row["arithmetic_source"]["base_curve"]
    coarse_genus = coarse.get("genus")
    normalization = row["similarity_normalization"]
    similarity_gap = int(normalization["literal_content"])
    order = row["clifford"]["integral_even_clifford_order"]
    local_index = int(order["local_level_index"])
    full_genus = full["genus"]
    kernel_index = full["stable_kernel_index_over_coarse"]
    coarse_exact = str(coarse.get("status", "")).startswith("PASS_EXACT")

    if full_genus in (0, 1):
        curve_tier = 0
    elif full_genus == 2:
        curve_tier = 1
    elif coarse_exact and coarse_genus in (0, 1):
        curve_tier = 2
    elif coarse_exact and coarse_genus == 2:
        curve_tier = 3
    elif coarse_genus in (0, 1, 2):
        curve_tier = 4
    else:
        curve_tier = 5

    if full_genus is not None:
        genus_key = int(full_genus)
    elif coarse_genus is not None:
        genus_key = int(coarse_genus)
    else:
        genus_key = 10**9

    split_penalty = 0 if row["clifford"]["quaternion_discriminant"] == 1 else 1
    return (
        curve_tier,
        0 if full["rational_non_CM_point"] else 1,
        genus_key,
        0 if kernel_index is not None else 1,
        int(kernel_index) if kernel_index is not None else 10**9,
        similarity_gap,
        local_index,
        split_penalty,
        int(row["determinant"]),
        row["surface_id"],
    )


def build(
    t_arithmetic: dict,
    classifier: dict,
    global_decisions: dict,
    t_arithmetic_path: Path = T_ARITHMETIC,
    classifier_path: Path = CLASSIFIER,
    global_decisions_path: Path = GLOBAL_DECISIONS,
) -> dict:
    if t_arithmetic["schema"] != "elkies-k3.rank7-t-arithmetic.v1":
        raise ValueError("unexpected T-arithmetic schema")
    if classifier["schema"] != "elkies-k3.rank19-arithmetic-marking-classifier.v1":
        raise ValueError("unexpected arithmetic-classifier schema")
    if global_decisions["schema"] != "elkies-k3.arithmetic-first-marking-decisions.v1":
        raise ValueError("unexpected arithmetic-first decision schema")

    classified = {row["surface_id"]: row for row in classifier["candidates"]}
    classifier_decisions = {
        surface_id: row["classification"] for surface_id, row in classified.items()
    }
    global_by_id = {
        row["surface_id"]: row for row in global_decisions["records"]
    }
    if len(global_by_id) != len(global_decisions["records"]):
        raise ValueError("duplicate arithmetic-first decision surface")
    if set(global_by_id) & set(classified):
        raise ValueError("arithmetic-first and rootless-classifier decisions overlap")
    t_surface_ids = {row["surface_id"] for row in t_arithmetic["surfaces"]}
    if not set(global_by_id) <= t_surface_ids:
        raise ValueError("arithmetic-first decision contains an unknown T row")
    global_replay = {
        surface_id: [
            validate_certificate(certificate)
            for certificate in decision["certificates"]
        ]
        for surface_id, decision in global_by_id.items()
    }
    rows_out = []
    for row in t_arithmetic["surfaces"]:
        surface_id = row["surface_id"]
        classifier_row = classified.get(surface_id)
        global_decision = global_by_id.get(surface_id)
        classification = (
            global_decision["classification"]
            if global_decision is not None
            else classifier_decisions.get(
                surface_id, "UNSCREENED_NO_FULL_MARKING_DECISION"
            )
        )
        full = (
            global_decision["full_marking_curve"]
            if global_decision is not None
            else full_curve_summary(classifier_row)
        )
        coarse = (
            classifier_row["full_discriminant_marking_curve"].get(
                "coarse_norm_one_curve",
                row["arithmetic_source"]["base_curve"],
            )
            if classifier_row is not None
            else row["arithmetic_source"]["base_curve"]
        )
        normalization = row["similarity_normalization"]
        order = row["clifford"]["integral_even_clifford_order"]
        priority = curve_priority(row, classifier_row)
        phase_2_status = (
            classifier_row.get("phase_2_certificate_status")
            if classifier_row is not None
            else global_decision.get("phase_2_certificate_status")
            if global_decision is not None
            else None
        )
        row_out = {
            "surface_id": surface_id,
            "determinant": int(row["determinant"]),
            "transcendental_gram": row["literal_transcendental_gram"],
            "rational_isotropy": row["rational_isotropy"]["isotropic"],
            "quaternion_discriminant": int(
                row["clifford"]["quaternion_discriminant"]
            ),
            "integral_order": {
                "reduced_discriminant": int(order["reduced_discriminant"]),
                "local_level_index": int(order["local_level_index"]),
                "source_status": row["arithmetic_source"]["status"],
            },
            "similarity_marking_gap": {
                "literal_content": int(normalization["literal_content"]),
                "quadratic_integrality_scale": int(
                    normalization["quadratic_integrality_scale"]
                ),
                "full_stable_kernel_still_required": (
                    full["status"]
                    not in {
                        "PASS_EXACT_STABLE_DISCRIMINANT_KERNEL_MODULAR_CURVE",
                        "PASS_EXACT_PROJECTIVE_STABLE_DISCRIMINANT_KERNEL_CURVE",
                        "PARTIAL_ABSTRACT_STABLE_ORTHOGONAL_CURVE_WITH_EXPLICIT_QQ_POINT",
                    }
                ),
            },
            "coarse_curve_diagnostic": {
                "status": coarse["status"],
                "label": coarse.get("label") or coarse.get("group"),
                "genus": coarse.get("genus"),
                "warning": (
                    "Coarse genus is a prioritization diagnostic only until the "
                    "stable discriminant-kernel curve is computed."
                ),
            },
            "full_marking_curve": full,
            "arithmetic_classification": classification,
            "phase_one_priority_key": list(priority[:-1]),
            "next_gate": (
                None
                if classification
                in {"ARITHMETICALLY_EXCLUDED", "ARITHMETICALLY_POSSIBLE"}
                else (
                    classifier_row.get("next_arithmetic_gate")
                    if phase_2_status is not None
                    and classifier_row is not None
                    and classifier_row.get("next_arithmetic_gate")
                    else "Compute the literal-lattice stable discriminant kernel, identify "
                    "a genus-0/1 (occasionally genus-2) quotient, and determine its "
                    "rational noncuspidal non-CM lifts."
                )
            ),
            "rootless_frame_data_used_in_priority": False,
        }
        if phase_2_status is not None:
            row_out["phase_2_certificate_status"] = phase_2_status
        if (
            classification in {"ARITHMETICALLY_EXCLUDED", "ARITHMETICALLY_POSSIBLE"}
            or classifier_row is not None
            and classifier_row.get("certificate_replay")
        ):
            row_out.update(
                {
                    "decision": (
                        global_decision["decision"]
                        if global_decision is not None
                        else classifier_row["classification_decision"]
                    ),
                    "certificate_replay": (
                        global_replay[surface_id]
                        if global_decision is not None
                        else classifier_row["certificate_replay"]
                    ),
                    "theorem_inputs": (
                        global_decision["theorem_inputs"]
                        if global_decision is not None
                        else classifier_row["theorem_inputs"]
                    ),
                }
            )
        rows_out.append(row_out)

    by_id = {row["surface_id"]: row for row in rows_out}
    if len(by_id) != len(rows_out):
        raise ValueError("duplicate T surface id")

    excluded = sorted(
        (row for row in rows_out if row["arithmetic_classification"] == "ARITHMETICALLY_EXCLUDED"),
        key=lambda row: row["surface_id"],
    )
    exact_positive = sorted(
        (row for row in rows_out if row["arithmetic_classification"] == "ARITHMETICALLY_POSSIBLE"),
        key=lambda row: row["surface_id"],
    )
    research_queue = sorted(
        (
            row
            for row in rows_out
            if row["arithmetic_classification"]
            not in {"ARITHMETICALLY_EXCLUDED", "ARITHMETICALLY_POSSIBLE"}
        ),
        key=lambda row: (*row["phase_one_priority_key"], row["surface_id"]),
    )
    low_genus_diagnostic = [
        row
        for row in research_queue
        if row["coarse_curve_diagnostic"]["genus"] in (0, 1, 2)
    ]

    # The only exact positive row is the already-realized H3 control.  A new
    # row must first acquire the same positive classification before its
    # orthogonal NS and rootless frames are inspected.
    assert [row["surface_id"] for row in exact_positive] == [H3_SURFACE_ID]
    new_positive_handoff = [
        row for row in exact_positive if row["surface_id"] != H3_SURFACE_ID
    ]
    assert not new_positive_handoff

    classification_counts = Counter(
        row["arithmetic_classification"] for row in rows_out
    )
    return {
        "schema": "elkies-k3.arithmetic-first-marked-t-foundry.v1",
        "status": "PASS_ARITHMETIC_FIRST_EMPTY_NEW_NS_ROOTLESS_HANDOFF",
        "policy": {
            "order": [
                "low-genus full marked Shimura/modular curve",
                "rational noncuspidal non-CM point",
                "NS = T^perp with a saturated rational marking",
                "rootless-frame test",
                "equation compilation",
            ],
            "priority_features": [
                "full marked genus 0 or 1, occasionally 2",
                "explicit rational non-CM point",
                "small stable-kernel index",
                "small literal-order similarity/marking gap",
            ],
            "priority_key_order": [
                "full/coarse curve tier",
                "rational non-CM point present",
                "best available genus",
                "stable-kernel index known",
                "stable-kernel index over coarse curve",
                "literal content similarity gap",
                "integral-order local level index",
                "split before division when otherwise tied",
                "determinant",
            ],
            "fail_closed": (
                "A low-genus coarse norm-one curve is not a full marking curve. "
                "No NS/rootless or equation handoff is emitted before an exact "
                "positive rational marking certificate."
            ),
        },
        "inputs": {
            relative(t_arithmetic_path): digest(t_arithmetic_path),
            relative(classifier_path): digest(classifier_path),
            relative(global_decisions_path): digest(global_decisions_path),
        },
        "accounting": {
            "transcendental_rows": len(rows_out),
            "classifications": dict(sorted(classification_counts.items())),
            "arithmetically_excluded": len(excluded),
            "exact_positive_controls": len(exact_positive),
            "arithmetic_research_queue": len(research_queue),
            "coarse_genus_at_most_two_diagnostics": len(low_genus_diagnostic),
            "new_positive_NS_rootless_handoff": len(new_positive_handoff),
        },
        "excluded_before_NS_or_equation_work": [
            {
                "surface_id": row["surface_id"],
                "determinant": row["determinant"],
                "full_marking_curve": row["full_marking_curve"],
                "decision": row["decision"],
                "certificate_replay": row["certificate_replay"],
                "theorem_inputs": row["theorem_inputs"],
            }
            for row in excluded
        ],
        "positive_controls": [
            {
                "surface_id": row["surface_id"],
                "determinant": row["determinant"],
                "already_realized": row["surface_id"] == H3_SURFACE_ID,
            }
            for row in exact_positive
        ],
        "curve_identification_queue": research_queue,
        "coarse_low_genus_diagnostic_shortlist": low_genus_diagnostic,
        "new_positive_NS_rootless_handoff": new_positive_handoff,
        "proof_boundary": {
            "proved": (
                "Every T-arithmetic row is ordered without using rootless-frame data, "
                "all exact marking decisions are propagated, and the downstream "
                "NS/rootless handoff is empty unless a new exact positive marking exists."
            ),
            "not_proved": (
                "Coarse low genus does not determine stable marked genus or rational "
                "points. Unscreened and UNKNOWN rows are neither constructions nor "
                "obstructions. This planner does not enumerate new ternary genera beyond "
                "the current rank-seven catalogue."
            ),
        },
        "reproduce": "python3 elkies-k3/scripts/build_arithmetic_first_marked_t_foundry.py",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t-arithmetic", type=Path, default=T_ARITHMETIC)
    parser.add_argument("--classifier", type=Path, default=CLASSIFIER)
    parser.add_argument("--decisions", type=Path, default=GLOBAL_DECISIONS)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    t_path = args.t_arithmetic.resolve()
    classifier_path = args.classifier.resolve()
    decisions_path = args.decisions.resolve()
    output_path = args.output.resolve()
    payload = build(
        json.loads(t_path.read_text()),
        json.loads(classifier_path.read_text()),
        json.loads(decisions_path.read_text()),
        t_path,
        classifier_path,
        decisions_path,
    )
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not output_path.exists() or output_path.read_text() != encoded:
            raise SystemExit(f"stale artifact: {output_path}")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(encoded)
    accounting = payload["accounting"]
    print(
        "TFIRST|T={}|excluded={}|research={}|coarse_g_le_2={}|"
        "new_NS_handoff={}|status=PASS".format(
            accounting["transcendental_rows"],
            accounting["arithmetically_excluded"],
            accounting["arithmetic_research_queue"],
            accounting["coarse_genus_at_most_two_diagnostics"],
            accounting["new_positive_NS_rootless_handoff"],
        )
    )


if __name__ == "__main__":
    main()
