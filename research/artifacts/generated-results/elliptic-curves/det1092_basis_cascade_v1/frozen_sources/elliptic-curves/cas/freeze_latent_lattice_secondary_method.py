#!/usr/bin/env python3
"""Freeze the unchanged latent-lattice method for E29 and ICARM 398--400."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
PRIMARY = ARTIFACTS / "latent_lattice_target_method_freeze_v1.json"
TARGET_METADATA = ROOT / "elliptic-curves/data/elkies_2026_r17_j_recognition_targets.json"
OUTPUT = ARTIFACTS / "latent_lattice_secondary_method_freeze_v1.json"
PRIMARY_SHA256 = "ef6f8b7be7a14095efa7529fb795d237e06465ba1cec023dcb4845287609c9f4"
TAG = "LATENT-LATTICE-E29-398-400-FROZEN-2026-09-01-v1"
TARGET_HASHES = {
    12: "206ef6992f433155d349618d55c289a00cf9014eb222da059c57e0db76131c0e",
    398: "5e09b5ed49cde24d20fcf300794e58a47f7f75ac7bba98c92b68ff3654df49f4",
    399: "92125a3aafd44ff45028ade826ef62338667164a402b3264507816e3c2009ead",
    400: "704f2292b4395923e2887a4ddd7a35f03d46baed42f5c75aac2dfe8519dd4275",
}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if digest(PRIMARY) != PRIMARY_SHA256:
        raise SystemExit("primary frozen method changed")
    primary = json.loads(PRIMARY.read_text())
    metadata = json.loads(TARGET_METADATA.read_text())
    observed = {
        int(record["icarm_id"]): record["source_sha256"]
        for record in metadata["targets"]
        if int(record["icarm_id"]) in TARGET_HASHES
    }
    if observed != TARGET_HASHES:
        raise SystemExit("secondary target metadata hashes changed")
    payload = {
        "schema": "elliptic-curves.latent-lattice-secondary-method-freeze.v1",
        "algorithm_tag": TAG,
        "status": "FROZEN_UNCHANGED_SECONDARY_TARGET_METHOD_NO_TUNING",
        "created_date": "2026-09-01",
        "inherits_algorithm_tag": primary["algorithm_tag"],
        "inherits_manifest_sha256": PRIMARY_SHA256,
        "algorithm_fields_changed": [],
        "target_records": [12, 398, 399, 400],
        "target_source_sha256": {str(key): value for key, value in TARGET_HASHES.items()},
        "dimension_acceptance": "one frozen dimension must recur in all four fibres",
        "hold_out_protocol": "all four choices of three training fibres and one untouched fibre",
        "target_metadata_was_known_after_primary_freeze": True,
        "no_tuning_rule": primary["no_tuning_rule"],
        "algorithm": {
            key: primary[key]
            for key in (
                "cloud_protocol",
                "dimension_protocol",
                "finite_protocol",
                "component_protocol",
                "forbidden_under_tag",
            )
        },
        "cloud_protocol": primary["cloud_protocol"],
        "dimension_protocol": primary["dimension_protocol"],
        "inputs": {
            str(PRIMARY.relative_to(ROOT)): digest(PRIMARY),
            str(TARGET_METADATA.relative_to(ROOT)): digest(TARGET_METADATA),
        },
        "proof_boundary": (
            "This manifest changes only the target identities and the logically necessary four-of-four "
            "recurrence/leave-one-out wording. Every computational parameter and score is inherited "
            "unchanged from the pre-wgxli primary freeze."
        ),
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit("secondary latent-lattice freeze is stale")
        print(f"LATENTSECONDARYFREEZE|check=PASS|tag={TAG}|sha256={sha256(rendered.encode()).hexdigest()}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(f"LATENTSECONDARYFREEZE|status=FROZEN|tag={TAG}|output={args.output}|sha256={sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
