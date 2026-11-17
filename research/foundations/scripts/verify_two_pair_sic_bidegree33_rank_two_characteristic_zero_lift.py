#!/usr/bin/env python3
"""Verify the finite characteristic-zero lift and its leading descent.

This wrapper combines the exact operator checks with an independent Singular
replay of the first coefficient of the rational telescoping identity.  Its
scope is deliberately finite: it does not claim the remaining 57 descent
levels, the terminal syzygy, or the endpoint identities over characteristic
zero.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts" / "generated-results"
OPERATOR = (
    ARTIFACTS
    / "two_pair_sic_bidegree33_rank_two_ore_characteristic_zero_lift.json"
)
OPERATOR_VERIFICATION = (
    ARTIFACTS
    / (
        "two_pair_sic_bidegree33_rank_two_ore_"
        "characteristic_zero_lift_verification.json"
    )
)
IMAGE_CACHE = (
    ARTIFACTS
    / "two_pair_sic_bidegree33_rank_two_ore_reconstruct_images.json"
)
CERTIFICATE = (
    ARTIFACTS
    / (
        "two_pair_sic_bidegree33_rank_two_char0_"
        "interior_divergence_certificate_m58_m57.sing"
    )
)
CHECKPOINT = (
    ARTIFACTS
    / (
        "two_pair_sic_bidegree33_rank_two_char0_"
        "interior_divergence_checkpoint_m57.json"
    )
)
RESIDUAL = CHECKPOINT.with_suffix(".poly")
RESEARCH_RESULT = (
    ARTIFACTS
    / (
        "two_pair_sic_bidegree33_rank_two_char0_"
        "interior_divergence_research_m58_m57.json"
    )
)
OUTPUT = (
    ARTIFACTS
    / "two_pair_sic_bidegree33_rank_two_characteristic_zero_lift.json"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command: list[str]) -> str:
    completed = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    output = completed.stdout + completed.stderr
    print(output, end="")
    return output


def main() -> None:
    operator_output = run(
        [
            sys.executable,
            str(
                ROOT
                / "scripts"
                / (
                    "verify_two_pair_sic_bidegree33_rank_two_ore_"
                    "characteristic_zero_lift.py"
                )
            ),
        ]
    )
    if "PASS 27 exact rational moment identities" not in operator_output:
        raise RuntimeError("exact operator verification did not complete")

    descent_output = run(
        [
            sys.executable,
            str(
                ROOT
                / "scripts"
                / (
                    "verify_two_pair_sic_bidegree33_rank_two_"
                    "interior_divergence_chunk.py"
                )
            ),
            "--operator",
            str(OPERATOR),
            "--certificate",
            str(CERTIFICATE),
            "--output-checkpoint",
            str(CHECKPOINT),
            "--timeout",
            "900",
        ]
    )
    if "PASS independently verified m^58 through m^58 over QQ" not in (
        descent_output
    ):
        raise RuntimeError("leading characteristic-zero descent did not replay")

    lift = json.loads(OPERATOR.read_text())
    verification = json.loads(OPERATOR_VERIFICATION.read_text())
    checkpoint = json.loads(CHECKPOINT.read_text())
    research = json.loads(RESEARCH_RESULT.read_text())
    if lift["operator"] != {
        "order": 14,
        "m_degree": 58,
        "coefficient_count": 885,
    }:
        raise RuntimeError("unexpected lifted-operator shape")
    if verification["all_cached_images_replayed"] != 205:
        raise RuntimeError("not all cached prime images were replayed")
    if verification["exact_rational_moment_rows"] != 27:
        raise RuntimeError("unexpected exact-moment verification range")
    if checkpoint["completed_m_degrees"] != [58, 58]:
        raise RuntimeError("unexpected characteristic-zero descent range")
    if checkpoint["next_m_degree"] != 57:
        raise RuntimeError("leading descent checkpoint is not at m^57")
    if research["steps"][0] != {
        "m_degree": 58,
        "certificate_m_degree": 57,
        "residual_terms": 3151,
        "X_terms": 5722,
        "Y_terms": 5769,
        "relation_terms": 0,
        "divergence_terms": 5631,
    }:
        raise RuntimeError("leading descent production statistics changed")

    files = [
        OPERATOR,
        OPERATOR_VERIFICATION,
        IMAGE_CACHE,
        CERTIFICATE,
        CHECKPOINT,
        RESIDUAL,
        RESEARCH_RESULT,
    ]
    result = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-"
            "characteristic-zero-lift-v1"
        ),
        "status": (
            "exact finite characteristic-zero operator lift and leading "
            "descent; all-order telescoping remains open"
        ),
        "point": 0,
        "operator": {"order": 14, "m_degree": 58},
        "operator_reconstruction": {
            "build_primes": 200,
            "fresh_holdout_primes": 5,
            "primitive_coefficient_maximum_bits": verification[
                "maximum_primitive_coefficient_bits"
            ],
            "exact_rational_moment_rows": 27,
        },
        "leading_descent": {
            "verified_m_degrees": [58, 58],
            "next_m_degree": 57,
            "certificate_X_terms": 5722,
            "certificate_Y_terms": 5769,
            "checkpoint_residual_terms": checkpoint["residual_terms"],
            "coefficient_ring": "QQ",
            "independent_singular_replay": True,
        },
        "remaining_gap": (
            "construct and independently replay levels m^57 through m^1, "
            "the terminal m^0 syzygy, and the two endpoint identities over QQ"
        ),
        "files_sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in files
        },
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS exact finite characteristic-zero operator lift")
    print("PASS exact leading m^58 descent over QQ")
    print("PASS all-order characteristic-zero certificate remains open")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
