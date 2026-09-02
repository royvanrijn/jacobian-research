#!/usr/bin/env python3
"""Expose the certified Golay determinant-720 target to foundry spectrum tools.

The Golay-octad target predates the consolidated lattice-foundry database.
This deliberately small adapter preserves its exact Gram matrix and the
rootless invariants needed by the generic degree-2/3/4 spectrum programs.  It
does not claim that the target occurred in the database's original search.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TARGET = (
    ROOT / "artifacts/generated-results/elkies-k3-golay-octad-rank17-det720.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-foundry-adapter-v1.json"
)
FRAME_ID = "G720-F001"
NS_ID = "G720"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def gram_digest(gram: list[list[int]]) -> str:
    text = "\n".join(" ".join(map(str, row)) for row in gram) + "\n"
    return hashlib.sha256(text.encode()).hexdigest()


def build(target_path: Path) -> dict:
    target = json.loads(target_path.read_text())
    if target.get("schema") != "elkies-k3.golay-octad-rank17-design.v1":
        raise ValueError("unexpected Golay target schema")
    frame = target["frame"]
    if (
        int(frame["rank"]) != 17
        or int(frame["determinant"]) != 720
        or int(frame["minimum_squared_norm"]) != 4
    ):
        raise ValueError("unexpected Golay target invariants")
    gram = frame["gram"]
    target_gram_sha256 = gram_digest(gram)

    adapted_frame = {
        "frame_id": FRAME_ID,
        "determinant": 720,
        "gram": gram,
        "gram_sha256": target_gram_sha256,
        "root_rank": 0,
        "root_type": "rootless",
        "signed_root_count": 0,
        "mw_rank_for_rho_19": 17,
        "rootless_intrinsics": {
            "minimum_squared_norm": 4,
            "norm_four_vectors": int(frame["norm_four_vectors"]),
            "norm_four_unoriented_pairs": int(frame["norm_four_unoriented_pairs"]),
            "automorphism_group_order": 32,
        },
    }
    return {
        "schema": "elkies-k3.golay-det720-foundry-adapter.v1",
        "status": "PASS_EXACT_TARGET_ADAPTER",
        "proof_boundary": (
            "This is a schema adapter for generic spectrum programs. It copies the "
            "certified target Gram exactly and records the exact PARI qfauto order 32; "
            "it adds no new K3, equation, curve, or arithmetic claim."
        ),
        "ns_classes": [
            {
                "ns_id": NS_ID,
                "determinant": 720,
                "frames": [adapted_frame],
            }
        ],
        "rootless_targets": [
            {
                "frame_id": FRAME_ID,
                "ns_id": NS_ID,
                "determinant": 720,
                "is_existing_H3_control": False,
            }
        ],
        "inputs": {relative(target_path): digest(target_path)},
        "reproduce": (
            f"python3 {relative(Path(__file__))} --target {relative(target_path)} "
            f"--output {relative(DEFAULT_OUTPUT)}"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    target_path = arguments.target.resolve()
    output_path = arguments.output.resolve()
    payload = build(target_path)
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if output_path.read_text() != serialized:
            raise SystemExit("Golay determinant-720 foundry adapter is stale")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "GOLAY720ADAPTER|frame=G720-F001|determinant=720|"
        "mw_rank=17|automorphism_order=32|status=PASS"
    )


if __name__ == "__main__":
    main()
