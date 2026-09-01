#!/usr/bin/env sage
"""Build the cross-Niemeier component-permutation/mod-2 priority ledger.

This is a deterministic experiment scheduler, not an embedding census.  It
uses repeated root-component types to identify backends whose component
permutation envelope can contain non-scalar involutions or order-four
elements.  For A7^2 D5^2 it additionally consumes the exact Dih_4 umbral
section and the existing six-frame negative control: the section contains
2B, 2C, and 4A classes, while all observed stabilizer images contain only
1A/2A and act trivially on M/2M.

Future embeddings are accepted into the high-priority experiment only after
their *full ambient* stabilizer is computed, its induced complement action is
reduced modulo two, and rank(g-I)>0 is certified.  Repeated components alone
are therefore a scheduling signal, never a mathematical claim about the
actual stabilizer.

status: EXACT_PRIORITY_LEDGER_HEURISTIC_BACKEND_ORDER
"""

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_NIEMEIER = (
    ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
)
DEFAULT_UMBRAL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-umbral-orbits-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-cross-niemeier-mod2-priority-v1.json"
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def component_type(component):
    return f"{component['family']}{component['rank']}"


def build(niemeier, umbral):
    assert niemeier["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    assert umbral["schema"] == "elkies-k3.lattice-foundry-umbral-orbits.v1"
    assert umbral["status"] == (
        "PASS_EXACT_AMBIENT_STABILIZERS_D2_ORBITS_AND_SAMPLED_D3_ORBITS"
    )
    exact_section_classes = Counter(row["class"] for row in umbral["group_section"])
    assert exact_section_classes == Counter(
        {"1A": 1, "2A": 1, "2B": 2, "2C": 2, "4A": 2}
    )
    observed_classes = sorted(
        {
            label
            for target in umbral["targets"]
            for label in target["full_ambient_stabilizer"][
                "umbral_image_classes"
            ]
        }
    )
    assert observed_classes == ["1A", "2A"]
    for target in umbral["targets"]:
        for action in target["full_ambient_stabilizer"]["induced_actions"]:
            if "2A" in action["compatible_umbral_classes"]:
                assert action["fixed_degree_two_rational_cosets"] == target[
                    "degree_two_rational_cosets"
                ]["count"]

    rows = []
    for ambient in niemeier["rooted_niemeier_lattices"]:
        multiplicities = Counter(
            component_type(component) for component in ambient["root_components"]
        )
        repeated = {
            key: value for key, value in sorted(multiplicities.items()) if value >= 2
        }
        swap_types = [key for key, value in repeated.items() if value >= 2]
        four_cycle_types = [key for key, value in repeated.items() if value >= 4]
        envelope_order = math.prod(math.factorial(value) for value in multiplicities.values())
        if ambient["label"] == "2A7_2D5":
            tier = 0
            reason = (
                "Exact Dih_4 section contains 2B, 2C, and 4A, but the six "
                "existing full stabilizers realize only the mod-2-trivial 1A/2A control."
            )
            requested = ["2B", "2C", "4A"]
            evidence = "EXACT_GROUP_SECTION_AND_NEGATIVE_CONTROL"
        elif four_cycle_types:
            tier = 1
            reason = (
                "At least four identical root components give the permutation "
                "envelope order-four and non-scalar involution candidates."
            )
            requested = [
                "order-4 component permutation",
                "non-scalar order-2 component permutation",
            ]
            evidence = "COMPONENT_ENVELOPE_HEURISTIC_ACTUAL_GX_AND_STABILIZER_OPEN"
        elif swap_types:
            tier = 2
            reason = (
                "Repeated identical root components give non-scalar swap candidates."
            )
            requested = ["non-scalar order-2 component permutation"]
            evidence = "COMPONENT_ENVELOPE_HEURISTIC_ACTUAL_GX_AND_STABILIZER_OPEN"
        else:
            tier = 3
            reason = (
                "No repeated root-component type; component permutations are not "
                "the first mod-2 symmetry route."
            )
            requested = []
            evidence = "LOW_COMPONENT_PERMUTATION_PRIORITY"
        rows.append(
            {
                "backend_id": f"ROOTED-{ambient['label']}",
                "ambient_label": ambient["label"],
                "priority_tier": tier,
                "component_type_multiplicities": dict(sorted(multiplicities.items())),
                "repeated_component_types": repeated,
                "full_symmetric_component_permutation_envelope_order": envelope_order,
                "swap_component_types": swap_types,
                "four_cycle_component_types": four_cycle_types,
                "requested_stabilizer_action_types": requested,
                "reason": reason,
                "evidence_status": evidence,
                "acceptance_gate": {
                    "full_ambient_stabilizer_required": True,
                    "chamber_section_intersection_alone_is_insufficient": True,
                    "induced_complement_action_required": True,
                    "mod2_moved_dimension_test": "rank_GF2(g_M - I) > 0",
                    "then_compute": (
                        "fixed-point and orbit distribution on the rational subset of M/2M"
                    ),
                },
            }
        )
    rows.sort(key=lambda row: (row["priority_tier"], row["ambient_label"]))
    tier_distribution = Counter(row["priority_tier"] for row in rows)
    assert len(rows) == 23
    return {
        "schema": "elkies-k3.cross-niemeier-mod2-priority.v1",
        "status": "PASS_EXACT_PRIORITY_LEDGER_HEURISTIC_BACKEND_ORDER",
        "proof_scope": {
            "proved": (
                "The root-component multiplicities of all 23 rooted Niemeier "
                "backends are classified into a deterministic experiment order. "
                "For A7^2 D5^2, the exact Dih_4 section and six-frame negative "
                "control are imported and checked."
            ),
            "not_proved": (
                "A repeated-component envelope does not prove that its permutations "
                "preserve Niemeier glue, stabilize a future auxiliary, or act "
                "nontrivially on M/2M. Every retained experiment must pass the "
                "declared full-stabilizer and induced-action gates."
            ),
        },
        "selection_policy": {
            "primary_signal": (
                "full stabilizer contains a non-scalar component permutation"
            ),
            "required_exact_signal": "rank_GF2(g_M-I)>0 on the complement",
            "seed_classes_for_A7_2_D5_2": ["2B", "2C", "4A"],
            "negative_control_classes": ["1A", "2A"],
            "use_after_gate": (
                "prioritize orbit-resolved rational bisection cosets and source searches"
            ),
        },
        "exact_A7_2_D5_2_control": {
            "section_class_distribution": dict(sorted(exact_section_classes.items())),
            "observed_full_stabilizer_image_classes": observed_classes,
            "frames_checked": len(umbral["targets"]),
            "mod2_result": (
                "all observed complement images are generated by +/-I and act trivially on M/2M"
            ),
        },
        "accounting": {
            "rooted_backends": len(rows),
            "priority_tier_distribution": {
                str(key): value for key, value in sorted(tier_distribution.items())
            },
        },
        "backends": rows,
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--niemeier", type=Path, default=DEFAULT_NIEMEIER)
parser.add_argument("--umbral", type=Path, default=DEFAULT_UMBRAL)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
result = build(
    json.loads(arguments.niemeier.read_text()),
    json.loads(arguments.umbral.read_text()),
)
result["inputs"] = {
    str(arguments.niemeier.relative_to(ROOT)): digest(arguments.niemeier),
    str(arguments.umbral.relative_to(ROOT)): digest(arguments.umbral),
}
encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not arguments.output.exists() or arguments.output.read_text() != encoded:
        raise SystemExit("cross-Niemeier mod-2 priority artifact is stale")
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(encoded)
print(
    "MOD2PRIORITY|backends={}|tiers={}|seed=2B,2C,4A|status=PASS".format(
        result["accounting"]["rooted_backends"],
        result["accounting"]["priority_tier_distribution"],
    )
)
