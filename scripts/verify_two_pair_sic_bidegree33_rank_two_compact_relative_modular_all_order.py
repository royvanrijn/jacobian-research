#!/usr/bin/env python3
"""Replay the compact order-64 modular all-order certificate."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "artifacts" / "generated-results"
LOCAL = ROOT / "artifacts" / "local"

OPERATOR = (
    GENERATED
    / "two_pair_sic_bidegree33_rank_two_compact_relative_pf_research.json"
)
CHECKPOINT_7 = (
    LOCAL
    / (
        "two_pair_sic_bidegree33_rank_two_compact_relative_shift_"
        "mod1000003_checkpoint_m7.json"
    )
)
CHECKPOINT_0 = (
    LOCAL
    / (
        "two_pair_sic_bidegree33_rank_two_compact_relative_shift_"
        "mod1000003_checkpoint_m0.json"
    )
)
CERTIFICATE_8_7 = (
    LOCAL
    / (
        "two_pair_sic_bidegree33_rank_two_compact_relative_shift_"
        "mod1000003_certificate_m8_m7.sing"
    )
)
CERTIFICATE_7_0 = (
    LOCAL
    / (
        "two_pair_sic_bidegree33_rank_two_compact_relative_shift_"
        "mod1000003_certificate_m7_m0.sing"
    )
)
TERMINAL_R = (
    LOCAL
    / (
        "two_pair_sic_bidegree33_rank_two_compact_relative_shift_"
        "mod1000003_terminal_syzygy_R.sing"
    )
)
TERMINAL_RESULT = (
    LOCAL
    / (
        "two_pair_sic_bidegree33_rank_two_compact_relative_shift_"
        "mod1000003_terminal_syzygy_research.json"
    )
)
ENDPOINT_RESULT = (
    LOCAL
    / (
        "two_pair_sic_bidegree33_rank_two_compact_relative_shift_"
        "mod1000003_endpoint_trace.json"
    )
)
OUTPUT = (
    GENERATED
    / (
        "two_pair_sic_bidegree33_rank_two_compact_relative_"
        "modular_all_order_certificate.json"
    )
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
        (CERTIFICATE_8_7, None, CHECKPOINT_7),
        (CERTIFICATE_7_0, CHECKPOINT_7, CHECKPOINT_0),
    )
    for certificate, input_checkpoint, output_checkpoint in chunks:
        command = [
            python,
            chunk_verifier,
            "--operator",
            str(OPERATOR),
            "--certificate",
            str(certificate),
            "--output-checkpoint",
            str(output_checkpoint),
            "--timeout",
            str(arguments.timeout),
        ]
        if input_checkpoint is not None:
            command.extend(["--input-checkpoint", str(input_checkpoint)])
        run(command, arguments.timeout + 30)

    run(
        [
            python,
            str(
                ROOT
                / "scripts"
                / "verify_two_pair_sic_bidegree33_rank_two_terminal_syzygy.py"
            ),
            "--operator",
            str(OPERATOR),
            "--terminal",
            str(CHECKPOINT_0.with_suffix(".poly")),
            "--certificate",
            str(TERMINAL_R),
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
        "--operator",
        str(OPERATOR),
    ]
    for certificate in (CERTIFICATE_8_7, CERTIFICATE_7_0):
        endpoint_command.extend(["--certificate", str(certificate)])
    endpoint_command.extend(
        [
            "--terminal-syzygy",
            str(TERMINAL_R),
            "--timeout",
            str(arguments.timeout),
            "--output",
            str(ENDPOINT_RESULT),
        ]
    )
    run(endpoint_command, arguments.timeout + 30)

    operator = json.loads(OPERATOR.read_text())
    terminal = json.loads(TERMINAL_RESULT.read_text())
    endpoint = json.loads(ENDPOINT_RESULT.read_text())
    checkpoint0 = json.loads(CHECKPOINT_0.read_text())
    if terminal["status"] != "exact modular final-syzygy certificate":
        raise RuntimeError("terminal certificate is not exact")
    if endpoint["status"] != "exact modular zero endpoint trace":
        raise RuntimeError("endpoint trace is not exact zero")
    if checkpoint0["next_m_degree"] != 0:
        raise RuntimeError("descending certificate did not reach m^0")

    certificate_files = (
        OPERATOR,
        CHECKPOINT_7,
        CHECKPOINT_7.with_suffix(".poly"),
        CHECKPOINT_0,
        CHECKPOINT_0.with_suffix(".poly"),
        CERTIFICATE_8_7,
        CERTIFICATE_7_0,
        TERMINAL_R,
        TERMINAL_RESULT,
        ENDPOINT_RESULT,
    )
    result = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-compact-relative-"
            "modular-all-order-certificate-v1"
        ),
        "status": (
            "exact all-order modular compact relative telescoping certificate"
        ),
        "prime": int(operator["modular_operator"]["prime"]),
        "point": operator["point"],
        "operator": {"order": 64, "m_degree": 8},
        "source_differential_operator": {"order": 8, "z_degree": 72},
        "coefficient_identities_verified": 8,
        "terminal_residual_before_syzygy_terms": checkpoint0[
            "residual_terms"
        ],
        "terminal_syzygy_R_terms": terminal["R_terms"],
        "terminal_residual_after_syzygy_terms": 0,
        "endpoint_trace": endpoint["status"],
        "all_order_conclusion": (
            "for every integer m>=0, the compact order-64, m-degree-8 "
            "shift operator annihilates the fixed-fiber normalized "
            "moment sequence over F_1000003"
        ),
        "characteristic_zero_bridge": (
            "the reconstructed characteristic-zero compact operator is "
            "an exact left multiple Q_50 of the lifted G_14; lifting this "
            "divergence and its endpoints to characteristic zero remains open"
        ),
        "files_sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in certificate_files
        },
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS all eight modular coefficient identities")
    print("PASS exact terminal Koszul correction")
    print("PASS all-order zero endpoint trace")
    print("PASS compact order-64 recurrence for every m>=0 modulo 1000003")
    print(f"PASS wrote {arguments.output}")


if __name__ == "__main__":
    main()
