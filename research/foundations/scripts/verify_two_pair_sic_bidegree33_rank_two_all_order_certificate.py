#!/usr/bin/env python3
"""Replay the complete fixed-fiber modular order-14 certificate."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts" / "generated-results"
OPERATOR = (
    ARTIFACTS
    / "two_pair_sic_bidegree33_rank_two_ore_reconstruct_research.json"
)
CHECKPOINT_38 = (
    ARTIFACTS
    / "two_pair_sic_bidegree33_rank_two_interior_divergence_checkpoint_m38.json"
)
CHECKPOINT_18 = (
    ARTIFACTS
    / "two_pair_sic_bidegree33_rank_two_interior_divergence_checkpoint_m18.json"
)
CHECKPOINT_0 = (
    ARTIFACTS
    / "two_pair_sic_bidegree33_rank_two_interior_divergence_checkpoint_m0.json"
)
CERTIFICATE_58_38 = (
    ARTIFACTS
    / "two_pair_sic_bidegree33_rank_two_interior_divergence_certificate_m58_m38.sing"
)
CERTIFICATE_38_18 = (
    ARTIFACTS
    / "two_pair_sic_bidegree33_rank_two_interior_divergence_certificate_m38_m18.sing"
)
CERTIFICATE_18_0 = (
    ARTIFACTS
    / "two_pair_sic_bidegree33_rank_two_interior_divergence_certificate_m18_m0.sing"
)
TERMINAL_R = (
    ARTIFACTS
    / "two_pair_sic_bidegree33_rank_two_terminal_syzygy_R.sing"
)
TERMINAL_RESULT = (
    ARTIFACTS
    / "two_pair_sic_bidegree33_rank_two_terminal_syzygy_block_research.json"
)
ENDPOINT_RESULT = (
    ARTIFACTS
    / "two_pair_sic_bidegree33_rank_two_endpoint_trace_research.json"
)
OUTPUT = (
    ARTIFACTS
    / "two_pair_sic_bidegree33_rank_two_all_order_certificate.json"
)


def run(command: list[str], timeout: int) -> None:
    completed = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=True,
    )
    if completed.stdout:
        print(completed.stdout, end="")
    if completed.stderr:
        print(completed.stderr, end="", file=sys.stderr)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        while chunk := source.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    python = sys.executable
    chunk_verifier = str(
        ROOT
        / "scripts"
        / "verify_two_pair_sic_bidegree33_rank_two_interior_divergence_chunk.py"
    )
    chunks = (
        (CERTIFICATE_58_38, None, CHECKPOINT_38),
        (CERTIFICATE_38_18, CHECKPOINT_38, CHECKPOINT_18),
        (CERTIFICATE_18_0, CHECKPOINT_18, CHECKPOINT_0),
    )
    for certificate, input_checkpoint, output_checkpoint in chunks:
        command = [
            python,
            chunk_verifier,
            "--certificate",
            str(certificate),
            "--output-checkpoint",
            str(output_checkpoint),
            "--timeout",
            str(arguments.timeout),
        ]
        if input_checkpoint is not None:
            command.extend(
                ["--input-checkpoint", str(input_checkpoint)]
            )
        run(command, arguments.timeout + 30)

    run(
        [
            python,
            str(
                ROOT
                / "scripts"
                / "verify_two_pair_sic_bidegree33_rank_two_terminal_syzygy.py"
            ),
            "--timeout",
            str(arguments.timeout),
        ],
        arguments.timeout + 30,
    )
    endpoint_command = [
        python,
        str(
            ROOT
            / "scripts"
            / "verify_two_pair_sic_bidegree33_rank_two_endpoint_trace.py"
        ),
    ]
    for certificate in (
        CERTIFICATE_58_38,
        CERTIFICATE_38_18,
        CERTIFICATE_18_0,
    ):
        endpoint_command.extend(["--certificate", str(certificate)])
    endpoint_command.extend(
        [
            "--terminal-syzygy",
            str(TERMINAL_R),
            "--timeout",
            str(arguments.timeout),
        ]
    )
    run(endpoint_command, arguments.timeout + 30)

    operator = json.loads(OPERATOR.read_text())
    terminal = json.loads(TERMINAL_RESULT.read_text())
    endpoint = json.loads(ENDPOINT_RESULT.read_text())
    checkpoint0 = json.loads(CHECKPOINT_0.read_text())
    if terminal["status"] != "exact modular final-syzygy certificate":
        raise RuntimeError("terminal certificate status is not exact")
    if endpoint["status"] != "exact modular zero endpoint trace":
        raise RuntimeError("endpoint certificate status is not exact zero")
    if checkpoint0["next_m_degree"] != 0:
        raise RuntimeError("descending certificate did not reach m^0")

    certificate_files = (
        OPERATOR,
        CHECKPOINT_38,
        CHECKPOINT_38.with_suffix(".poly"),
        CHECKPOINT_18,
        CHECKPOINT_18.with_suffix(".poly"),
        CHECKPOINT_0,
        CHECKPOINT_0.with_suffix(".poly"),
        CERTIFICATE_58_38,
        CERTIFICATE_38_18,
        CERTIFICATE_18_0,
        TERMINAL_R,
        TERMINAL_RESULT,
        ENDPOINT_RESULT,
    )
    result = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-"
            "all-order-certificate-v1"
        ),
        "status": (
            "exact all-order modular fixed-fiber telescoping certificate"
        ),
        "prime": int(operator["modular_operator"]["prime"]),
        "point": operator["point"],
        "operator": {"order": 14, "m_degree": 58},
        "coefficient_identities_verified": 58,
        "terminal_residual_before_syzygy_terms": checkpoint0[
            "residual_terms"
        ],
        "terminal_syzygy_R_terms": terminal["R_terms"],
        "terminal_residual_after_syzygy_terms": 0,
        "endpoint_trace": endpoint["status"],
        "all_order_conclusion": (
            "for every integer m>=0, the stored order-14 operator "
            "annihilates the fixed-fiber normalized moment sequence over "
            "F_1000003"
        ),
        "scope_limit": (
            "this does not reconstruct characteristic zero, prove a "
            "generic rank-two parameter identity, or classify the "
            "exceptional locus"
        ),
        "files_sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in certificate_files
        },
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS all 58 modular coefficient identities")
    print("PASS exact terminal Koszul correction")
    print("PASS all-order zero endpoint trace")
    print(
        "PASS order-14 recurrence for every m>=0 at the fixed fiber "
        "modulo 1000003"
    )
    print(f"PASS wrote {arguments.output}")


if __name__ == "__main__":
    main()
