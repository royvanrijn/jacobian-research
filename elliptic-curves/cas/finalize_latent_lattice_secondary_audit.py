#!/usr/bin/env python3
"""Consolidate the four isolated runs of the frozen E29/398--400 audit."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results/elliptic-curves"
LOCAL = ROOT / "artifacts/local/elliptic-curves"
FREEZE = GENERATED / "latent_lattice_secondary_method_freeze_v1.json"
E29_PARTIAL = GENERATED / "latent_lattice_secondary_frozen_dimension_v1.json.partial"
CURVE399 = LOCAL / "latent_lattice_secondary_curve399_v1.json"
CURVE400 = LOCAL / "latent_lattice_secondary_curve400_v1.json"
OUTPUT = GENERATED / "latent_lattice_secondary_frozen_dimension_v1.json"

TAG = "LATENT-LATTICE-E29-398-400-FROZEN-2026-09-01-v1"
FREEZE_SHA256 = "8795cdd203ba1c698e0f0534c14a45a91c37a6f4d405e795a9f2f295f86bfcba"
INPUT_SHA256 = {
    E29_PARTIAL: "1be34294f5576cbaeca612f2dc0d340c8025338366a8280517ea3447e1aa4a94",
    CURVE399: "f0a45d9cd5279ea9fc5f40fd611cbd15f919f2679c0e10c3a78b77fd71652d92",
    CURVE400: "fc76849ee2d0b323f75435a9a7b424201a5129ed82c54a55ebb4a15928e4fcd7",
}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def checked_json(path: Path, expected_hash: str) -> dict[str, object]:
    observed = digest(path)
    if observed != expected_hash:
        raise SystemExit(f"input changed: {path}: {observed}")
    return json.loads(path.read_text())


def build_payload() -> dict[str, object]:
    if digest(FREEZE) != FREEZE_SHA256:
        raise SystemExit("secondary freeze changed")
    partial = checked_json(E29_PARTIAL, INPUT_SHA256[E29_PARTIAL])
    curve399 = checked_json(CURVE399, INPUT_SHA256[CURVE399])
    curve400 = checked_json(CURVE400, INPUT_SHA256[CURVE400])
    if partial.get("algorithm_tag") != TAG or partial.get("freeze_sha256") != FREEZE_SHA256:
        raise SystemExit("E29 partial has the wrong frozen-method identity")
    e29_fibres = partial.get("completed_fibres", [])
    if [record.get("curve_id") for record in e29_fibres] != [12]:
        raise SystemExit("E29 partial does not contain exactly curve 12")
    isolated = {record["fibres"][0]["curve_id"]: record for record in (curve399, curve400)}
    if set(isolated) != {399, 400}:
        raise SystemExit("isolated audits do not contain curves 399 and 400")
    freeze = json.loads(FREEZE.read_text())
    fibres = [
        e29_fibres[0],
        {
            "curve_id": 398,
            "ambient_rank": 30,
            "rank_lower_bound": 30,
            "status": "FAIL_FROZEN_RESOURCE_BOUND",
            "failed_stage": "adaptive_short_vector_enumeration",
            "failed_height_bound": 20,
            "pari_timeout_seconds": 600,
            "exception_type": "subprocess.TimeoutExpired",
            "selected_dimension": None,
            "interpretation": (
                "The first frozen height-bound enumeration did not finish inside the "
                "unchanged resource limit. No larger timeout or alternative bound was tried."
            ),
        },
        isolated[399]["fibres"][0],
        isolated[400]["fibres"][0],
    ]
    return {
        "schema": "elliptic-curves.latent-lattice-secondary-frozen-audit.v1",
        "algorithm_tag": TAG,
        "freeze_sha256": FREEZE_SHA256,
        "status": "FAIL_FROZEN_SECONDARY_GATE_RESOURCE_AND_DIMENSION",
        "target_source_sha256": freeze["target_source_sha256"],
        "selected_dimensions": [None, None, 12, 16],
        "dimensions_recurring_in_all_four_fibres": [],
        "component_stage": "NOT_RUN_FROZEN_DIMENSION_GATE_FAILED",
        "fibres": fibres,
        "isolated_execution": {
            "reason": "A resource failure on curve 398 was isolated so it could not suppress curves 399 and 400.",
            "algorithm_fields_changed": [],
            "resource_limit_changed": False,
        },
        "research_disposition": "PARKED_AFTER_SECOND_INDEPENDENT_FROZEN_TARGET_FAILURE",
        "proof_boundary": (
            "This is a complete fail-closed audit of the frozen selector on the four named records. "
            "It proves neither that their Mordell--Weil groups have no common primitive sublattice "
            "of rank 10 through 20 nor that curves 399 and 400 have generic ranks 12 and 16. The "
            "reported dimensions are statistical selector outputs, and curve 398 is a resource-bound "
            "failure. The precommitted all-four recurrence gate failed, so no relation-component, "
            "Hermite-shape, finite-index, height-lattice, or equation stage was authorized."
        ),
        "inputs": {
            str(FREEZE.relative_to(ROOT)): FREEZE_SHA256,
            **{str(path.relative_to(ROOT)): value for path, value in INPUT_SHA256.items()},
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = json.dumps(build_payload(), indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit("secondary frozen audit artifact is stale")
        print(f"LATENTSECONDARYAUDIT|check=PASS|sha256={sha256(rendered.encode()).hexdigest()}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        "LATENTSECONDARYAUDIT|status=FAIL_FROZEN_SECONDARY_GATE_RESOURCE_AND_DIMENSION|"
        f"output={args.output}|sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
