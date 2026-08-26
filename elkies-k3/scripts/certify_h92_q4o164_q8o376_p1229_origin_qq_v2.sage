#!/usr/bin/env sage -python
"""Attach the preferred P1229 origin to the exact q8/orbit376 quartic.

Version 2 consumes the actual ``degree_one_candidates`` field emitted by the
known-section pointing artifact.  The equation curve is first identified as
P1229 independently by exact MW9 height/pairing data, then its degree-one
restriction is used as the quartic origin.
"""

import hashlib
import json
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
RR = LOCAL / "q4o164-q8o376-rr-p2-qq.json"
POINTING = LOCAL / "q4o164-q8o376-known-section-pointing-qq.json"
IDENTIFICATION = LOCAL / "q4o164-transported-degree1-marked-classes-qq.json"
HANDOFF = GENERATED / "elkies-k3-h3-q4o1584-route-optimization-handoff.json"
OUTPUT = LOCAL / "q4o164-q8o376-p1229-origin-qq.json"
INPUTS = (RR, POINTING, IDENTIFICATION, HANDOFF)
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def emit(status, **fields):
    payload = {
        "schema": "elkies-k3.q4o164-q8o376-p1229-origin-qq.v2",
        "status": status,
        **fields,
        "inputs": {
            "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
            "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
        },
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"Q8O376P1229|status={status}|output={OUTPUT}", flush=True)
    return payload


for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")
rr = json.loads(RR.read_text())
pointing = json.loads(POINTING.read_text())
identification = json.loads(IDENTIFICATION.read_text())
handoff = json.loads(HANDOFF.read_text())

assert rr["status"] == "PASS_EXACT_QQ_Q4O164_Q8O376_UNPOINTED_RR_AND_4A1_JACOBIAN"
assert pointing["status"] in {
    "PASS_EXACT_QQ_Q4O164_Q8O376_KNOWN_SECTION_POINTING",
    "NO_KNOWN_SECTION_Q4O164_Q8O376_DEGREE_ONE_POINTING",
}
assert identification["status"] in {
    "PASS_EXACT_QQ_Q4O164_TRANSPORTED_DEGREE1_P1229_IDENTIFIED",
    "PASS_EXACT_QQ_Q4O164_TRANSPORTED_DEGREE1_MARKED_CLASSES",
}
assert handoff["status"] == "PASS_EXACT_Q4O1584_Q4O164_Q8O376_Q12_ROOTLESS_OPTIONS_HANDOFF_PROMOTED"

q8_frames = []
for value in handoff.values():
    if not isinstance(value, dict):
        continue
    for record in value.get("compiler_frames", []):
        if record.get("edge") == "q8/o376":
            q8_frames.append(record)
if not q8_frames:
    raise ArithmeticError("handoff has no q8/o376 compiler frame")
if any(record.get("selected_child_zero") != "P1229" for record in q8_frames):
    raise ArithmeticError("certified q8 route does not select P1229 as child zero")

matches = [record for record in identification["sections"] if record["matches_P1229"]]
negative_matches = [record for record in identification["sections"] if record["negative_matches_P1229"]]
if len(matches) + len(negative_matches) != 1:
    emit(
        "OPEN_Q4O164_Q8O376_P1229_EQUATION_ORIGIN",
        P1229_direct_matches=len(matches),
        P1229_negative_matches=len(negative_matches),
        reason="transported exact sections do not yet identify a unique P1229 sign",
    )
    raise SystemExit(2)

marked = matches[0] if matches else negative_matches[0]
transported_index = int(marked["transported_index"])
required_sign = 1 if matches else -1
name = f"transported_one_node_{transported_index}"
point_candidates = [
    record for record in pointing.get("degree_one_candidates", [])
    if record.get("name") == name and int(record.get("sign", 0)) == required_sign
]
if len(point_candidates) != 1:
    emit(
        "OPEN_Q4O164_Q8O376_P1229_RESTRICTION_NOT_DEGREE_ONE",
        transported_index=transported_index,
        required_sign=required_sign,
        restriction_degrees=pointing.get("restriction_degrees", []),
    )
    raise SystemExit(3)

point = point_candidates[0]
if int(point["new_base_degree"]) != 1:
    raise ArithmeticError("P1229 restriction is not degree one")
if not point.get("exact_surface_to_quartic_identity"):
    raise ArithmeticError("P1229 restriction lacks exact surface-to-quartic identity")
if not point.get("exact_pointed_child_invariant_identity"):
    raise ArithmeticError("P1229 pointing lacks exact child invariant identity")

payload = emit(
    "PASS_EXACT_QQ_Q4O164_Q8O376_P1229_POINTED_ORIGIN",
    child={"ADE": "4A1", "MW_rank_if_rho19": 13, "selected_zero": "P1229"},
    P1229_identification={
        "transported_index": transported_index,
        "equation_sign": required_sign,
        "marked_MW9_tail": marked["marked_MW9_tail"],
        "canonical_height_on_q4o164": marked["canonical_height"],
        "exact_height_pairing_identification": True,
    },
    quartic_point={
        "old_base_as_function_of_new_base": point["old_base_as_function_of_new_base"],
        "quartic_ordinate_as_function_of_new_base": point["quartic_ordinate_as_function_of_new_base"],
        "pointed_generalized_coefficients": point["pointed_generalized_coefficients"],
        "exact_surface_to_quartic_identity": True,
        "exact_pointed_child_invariant_identity": True,
    },
    route_lock={
        "q8_edge": "q8/orbit376",
        "selected_child_zero": "P1229",
        "matches_promoted_route": True,
    },
    method={
        "large_Groebner_required": False,
        "resolved_I4_chart_required_for_zero": False,
        "runtime_seconds": time.monotonic() - started,
    },
    next_required=(
        "Prefer direct q12/orbit5867 construction from this P1229-pointed 4A1 model. "
        "Only attach child I2 component labels if the q12 RR compiler actually needs them."
    ),
    proof_boundary=(
        "The exact equation curve used as quartic origin is independently identified with "
        "the inherited marked class P1229 and has degree one over the exact q8 pencil. "
        "Thus the invariant 4A1 Jacobian is pointed by the route-selected zero."
    ),
)
