#!/usr/bin/env python3
"""Audit p-adic obstructions for the final norm-four collision seeds."""

from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
DEFAULT_SCREEN = (
    GENERATED
    / "elkies-k3-r17-norm12-11952-genus3-normalization-full-p17-v1.json"
)
DEFAULT_PAIRS = (
    GENERATED
    / "elkies-k3-r17-norm12-11952-genus3-cover-collision-intersection-v1.json"
)
DEFAULT_MIXED = (
    GENERATED
    / "elkies-k3-r17-norm12-11952-genus3-normalization-mixed-collision-970-969-p17-v1.json"
)
DEFAULT_OUTPUT = (
    GENERATED
    / "elkies-k3-r17-norm12-11952-genus3-cover-collision-p17-hensel-audit-v1.json"
)


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def cover_key(coefficients, prime):
    values = [int(value) % prime for value in coefficients]
    if len(values) != 3 or values[2] == 0:
        return None
    inverse = pow(values[2], -1, prime)
    return (
        tuple(value * inverse % prime for value in values),
        pow(values[2], (prime - 1) // 2, prime),
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--screen", type=Path, default=DEFAULT_SCREEN)
parser.add_argument("--pairs", type=Path, default=DEFAULT_PAIRS)
parser.add_argument("--mixed", type=Path, default=DEFAULT_MIXED)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

screen_path = args.screen.resolve()
pairs_path = args.pairs.resolve()
mixed_path = args.mixed.resolve()
screen = json.loads(screen_path.read_text())
pairs_payload = json.loads(pairs_path.read_text())
mixed = json.loads(mixed_path.read_text())
if screen.get("schema") != (
    "elkies-k3.r17-norm12-direct-genus3-normalization-modp-search.v1"
) or int(screen["prime"]) != 17:
    raise ValueError("expected the full p=17 norm-four screen")
if pairs_payload.get("schema") != (
    "elkies-k3.r17-norm12-direct-genus3-cover-collision-intersection.v1"
):
    raise ValueError("unexpected final-pair schema")
if mixed.get("schema") != (
    "elkies-k3.r17-norm12-direct-genus3-mixed-collision-hensel.v1"
):
    raise ValueError("unexpected mixed-collision schema")

pairs = [
    tuple(map(int, record["trace_masks"]))
    for record in pairs_payload["surviving_distinct_norm_four_trace_pairs"]
]
by_mask = defaultdict(list)
for survivor in screen["survivors"]:
    by_mask[int(survivor["translation_orbit_mask"])].append(survivor)

degree_eight_indices = {}
for mask, survivors in by_mask.items():
    filtered = [
        survivor
        for survivor in survivors
        if len(survivor["branch_coefficients_low_to_high"]) == 9
        and len(survivor["removed_square_factor_coefficients_low_to_high"]) == 4
    ]
    for index, survivor in enumerate(filtered):
        degree_eight_indices[id(survivor)] = index

audit_cache = {}
audit_paths = []


def degree_eight_audit(survivor):
    if id(survivor) not in degree_eight_indices:
        return None
    trace_index = int(survivor["trace_index"])
    seed_index = degree_eight_indices[id(survivor)]
    key = trace_index, seed_index
    if key not in audit_cache:
        path = (
            GENERATED
            / f"elkies-k3-r17-norm12-11952-genus3-normalization-trace{trace_index}"
            f"-seed{seed_index}-p17-resolved32-v1.json"
        )
        payload = json.loads(path.read_text())
        if payload.get("schema") != (
            "elkies-k3.r17-norm12-direct-genus3-normalization-hensel.v1"
        ):
            raise ValueError(f"unexpected Hensel schema: {path}")
        if len(payload["lifts"]) != 1:
            raise ValueError(f"expected one audited seed: {path}")
        lift = payload["lifts"][0]
        if int(lift["trace_index"]) != trace_index or int(
            lift["survivor_index_within_trace"]
        ) != seed_index:
            raise ValueError(f"Hensel seed mismatch: {path}")
        resolved = lift["resolved_singular_hensel_branches"]
        audit_cache[key] = {
            "path": path,
            "obstructed": (
                not resolved["truncated"]
                and int(resolved["final_extendable_state_count"]) == 0
            ),
            "jacobian_rank": int(lift["jacobian_rank"]),
            "levels": resolved["levels"],
        }
        audit_paths.append(path)
    return audit_cache[key]


collision_seed_records = []
unresolved = []
for left_mask, right_mask in pairs:
    for left_index, left in enumerate(by_mask[left_mask]):
        left_key = cover_key(left["reduced_quadratic_coefficients_low_to_high"], 17)
        for right_index, right in enumerate(by_mask[right_mask]):
            if left_key != cover_key(
                right["reduced_quadratic_coefficients_low_to_high"], 17
            ):
                continue
            left_audit = degree_eight_audit(left)
            right_audit = degree_eight_audit(right)
            eliminating_side = None
            eliminating_audit = None
            if left_audit is not None and left_audit["obstructed"]:
                eliminating_side, eliminating_audit = "left", left_audit
            elif right_audit is not None and right_audit["obstructed"]:
                eliminating_side, eliminating_audit = "right", right_audit

            record = {
                "trace_masks": [left_mask, right_mask],
                "trace_indices": [
                    int(left["trace_index"]), int(right["trace_index"])
                ],
                "survivor_indices_within_all_trace_survivors": [
                    left_index, right_index
                ],
                "raw_branch_degrees": [
                    len(left["branch_coefficients_low_to_high"]) - 1,
                    len(right["branch_coefficients_low_to_high"]) - 1,
                ],
                "slopes_m0_m1_m2": [left["m0_m1_m2"], right["m0_m1_m2"]],
            }
            if eliminating_audit is not None:
                record.update(
                    {
                        "status": "PASS_ONE_SIDE_OBSTRUCTED_BEFORE_P3",
                        "eliminating_side": eliminating_side,
                        "eliminating_hensel_artifact": relative(
                            eliminating_audit["path"]
                        ),
                    }
                )
            else:
                expected_mixed = (
                    int(left["trace_index"]) == int(mixed["degree_six_trace_index"])
                    and left_index == int(mixed["degree_six_survivor_index"])
                    and int(right["trace_index"])
                    == int(mixed["degree_eight_trace_index"])
                    and right_index == int(mixed["degree_eight_survivor_index"])
                )
                if expected_mixed and mixed["status"] == (
                    "PASS_SECOND_ORDER_OBSTRUCTION_NO_PRIME_CUBED_LIFT"
                ):
                    record.update(
                        {
                            "status": "PASS_COUPLED_MIXED_SYSTEM_OBSTRUCTED_BEFORE_P3",
                            "eliminating_hensel_artifact": relative(mixed_path),
                        }
                    )
                else:
                    record["status"] = "UNRESOLVED_P17_COLLISION_SEED"
                    unresolved.append(record)
            collision_seed_records.append(record)

if not collision_seed_records:
    raise ValueError("no p=17 collision seeds found")

output_path = args.output if args.output.is_absolute() else ROOT / args.output
output_path.parent.mkdir(parents=True, exist_ok=True)
input_paths = [screen_path, pairs_path, mixed_path] + sorted(set(audit_paths))
result = {
    "schema": "elkies-k3.r17-norm12-direct-genus3-collision-hensel-audit.v1",
    "status": (
        "PASS_ALL_DISPLAYED_P17_AFFINE_COLLISION_SEEDS_OBSTRUCTED_BEFORE_P3"
        if not unresolved
        else "UNRESOLVED_P17_AFFINE_COLLISION_SEEDS"
    ),
    "source_label": "norm12-orbit-11952",
    "prime": 17,
    "trace_pair_count": len(pairs),
    "collision_seed_pair_count": len(collision_seed_records),
    "degree_eight_seed_artifact_count": len(set(audit_paths)),
    "collision_seeds": collision_seed_records,
    "unresolved_collision_seeds": unresolved,
    "proof_boundary": (
        "For every affine p=17 cover-collision seed on the four trace pairs "
        "surviving the five-prime trace-pair intersection, at least one individual "
        "normalization system or the sole mixed coupled system has no lift to "
        "17^3. This exactly excludes those p-adic residue classes. It does not "
        "exclude rational collisions nonintegral at 17, boundary-chart points, or "
        "trace pairs removed only by the earlier modular intersections."
    ),
    "inputs": {relative(path): digest(path) for path in input_paths},
    "reproducing_command": (
        "python3 "
        "elkies-k3/scripts/analyze_r17_norm12_direct_genus3_hensel_collisions.py"
    ),
}
output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(
    "R17GENUS3HENSELAUDIT"
    f"|trace_pairs={len(pairs)}|collision_seeds={len(collision_seed_records)}"
    f"|unresolved={len(unresolved)}|output={relative(output_path)}"
    f"|status={result['status']}",
    flush=True,
)
