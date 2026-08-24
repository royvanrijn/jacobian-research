#!/usr/bin/env python3
"""Pin an exact rank-at-least-16 certificate for direct Nagao-root T=62/35.

The source rational scan selected sixteen exact points numerically at both
72- and 120-digit height precision.  Numerical rank is used only to choose
the input set.  This script independently checks exact curve membership,
runs one foreground-capped PARI small-prime saturation, and proves that the
returned sixteen exact points are independent by their combined images in
``E(F_p)/3E(F_p)``.  Rational 3-torsion is excluded at a separate good prime.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import json
import os
from pathlib import Path
import platform
import shlex
import sys
from typing import Any

from extend_nagao_u42_frontier import saturate_exact_basis
from search_mestre_02557104116148_direct_rational import (
    ROOTS,
    family_coefficients,
)
from search_mestre_root_tuple_scale import (
    point_digest,
    point_on_short_curve,
    sha256_file,
)
from search_mestre_root_tuple_scale_max100 import stable_json_digest
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

PARAMETER = Q(62, 35)
SOURCE_STAGE = "H1000000"
SOURCE_ARTIFACT = "elliptic_mestre_02557104116148_direct_rational.json"
EXPECTED_SOURCE_ARTIFACT_SHA256 = (
    "4874478c553c81ed69fffb49738b5975900a26a17d96f4dca9203a8244e75db6"
)
EXPECTED_SOURCE_RESULT_SHA256 = (
    "c8d231506669c58fbb74b7e9d19b742a15564881f72b12c02d42fc9b1dadb687"
)
EXPECTED_INPUT_POINT_SHA256 = (
    "565cd0d79b4c9ca44f27455bc6c56c1547b3a7476fff6aecfa97f21c5424a886"
)
EXPECTED_SATURATED_BASIS_SHA256 = (
    "7d2766fd064e9b6dafa604e8ced9cb1b4ed13b654eda16f0097b831f32da788d"
)
EXPECTED_CERTIFICATE_PRIMES = (
    17, 43, 47, 59, 61, 67, 73, 79,
    97, 101, 103, 131, 137, 149, 151,
)
EXPECTED_CONDUCTOR = (
    "2685930484353663888032391590962133596851202601813686353960195"
)
EXPECTED_LOG_CONDUCTOR = (
    "139.143132796676806403103988033287761319526891323675969649953"
)
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/"
    "elliptic_mestre_02557104116148_t62_35_rank16_certificate.json"
)


def load_source(path: Path) -> tuple[dict[str, Any], dict[str, Any], tuple]:
    if sha256_file(path) != EXPECTED_SOURCE_ARTIFACT_SHA256:
        raise AssertionError("the frozen direct-T source artifact changed")
    artifact = json.loads(path.read_text())
    if artifact["result_sha256"] != EXPECTED_SOURCE_RESULT_SHA256:
        raise AssertionError("the frozen direct-T result digest changed")
    records = [
        record
        for record in artifact["selected_records"]
        if Q(record["parameter"]) == PARAMETER
    ]
    if len(records) != 1:
        raise AssertionError("the exact T=62/35 source record changed")
    record = records[0]
    stage = record["point_stages"][SOURCE_STAGE]
    points = tuple(
        (Q(point["x"]), Q(point["y"]))
        for point in stage["numerical_subset"]
    )
    coefficients = family_coefficients(PARAMETER)
    if (
        stage["status"] != "completed"
        or int(stage["stable_numerical_rank"]) != 16
        or len(points) != 16
        or point_digest(points) != EXPECTED_INPUT_POINT_SHA256
        or any(not point_on_short_curve(coefficients, point) for point in points)
    ):
        raise AssertionError("the exact H1m input basis changed")
    conductor = record["conductor_phase"]
    if (
        conductor["conductor"] != EXPECTED_CONDUCTOR
        or conductor["log_conductor"] != EXPECTED_LOG_CONDUCTOR
        or conductor["root_number"] != 1
    ):
        raise AssertionError("the exact conductor record changed")
    return record, stage, points


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=root / "artifacts/generated-results" / SOURCE_ARTIFACT,
    )
    parser.add_argument("--saturation-prime-bound", type=int, default=50)
    parser.add_argument("--saturation-timeout", type=float, default=60.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=1_000)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--output", type=Path, default=root / DEFAULT_OUTPUT)
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.saturation_prime_bound != 50:
        raise SystemExit("the saturation prime bound is pinned at 50")
    if not 0 < args.saturation_timeout <= 60:
        raise SystemExit("the one saturation cap must lie in (0,60]")
    if args.certificate_prime_bound != 1_000:
        raise SystemExit("the exact finite-reduction prime bound is pinned at 1000")
    if args.stack_bytes != 512_000_000:
        raise SystemExit("the PARI stack is pinned at 512000000 bytes")


def exclusive_write(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")


def main() -> None:
    args = build_parser().parse_args()
    validate_args(args)
    if args.output.exists():
        raise SystemExit("refusing to overwrite the exact rank-16 certificate")
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]
    record, stage, input_points = load_source(args.source)
    coefficients = family_coefficients(PARAMETER)

    saturated, saturation = saturate_exact_basis(
        coefficients,
        input_points,
        prime_bound=args.saturation_prime_bound,
        timeout=args.saturation_timeout,
        stack_bytes=args.stack_bytes,
    )
    if (
        len(saturated) != 16
        or saturation["saturated_basis_sha256"]
        != EXPECTED_SATURATED_BASIS_SHA256
        or saturation["height_determinant_ratio"]
        != "1073741824.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
        or any(not point_on_short_curve(coefficients, point) for point in saturated)
    ):
        raise AssertionError("the exact small-prime saturation replay changed")

    certificate = mod3_independence_certificate(
        coefficients,
        saturated,
        prime_bound=args.certificate_prime_bound,
    )
    if (
        certificate["status"] != "certified exact algebraic rank lower bound"
        or certificate["certified_algebraic_rank_lower_bound"] != 16
        or certificate["combined_exact_rank_over_F3"] != 16
        or tuple(certificate["certificate_primes"])
        != EXPECTED_CERTIFICATE_PRIMES
        or certificate["point_sha256"] != EXPECTED_SATURATED_BASIS_SHA256
        or certificate["rational_3_torsion_exclusion"] != {
            "prime": 23,
            "group_order": 32,
            "reason": "rational prime-to-p torsion injects at good reduction",
        }
    ):
        raise AssertionError("the exact mod-3 certificate replay changed")

    conductor = record["conductor_phase"]
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "certified exact algebraic rank lower bound 16",
        "theorem": {
            "statement": (
                "The displayed elliptic curve over Q has Mordell-Weil rank "
                "at least 16."
            ),
            "proof_method": (
                "sixteen exact rational points have full combined rank in "
                "E(F_p)/3E(F_p); rational 3-torsion is excluded at good "
                "reduction p=23"
            ),
            "certified_algebraic_rank_lower_bound": 16,
        },
        "curve": {
            "family_roots": list(ROOTS),
            "parameter_T": str(PARAMETER),
            "short_weierstrass_coefficients": [str(value) for value in coefficients],
            "minimal_model": conductor["minimal_model"],
            "conductor": conductor["conductor"],
            "log_conductor": conductor["log_conductor"],
            "below_strict_log_conductor_182_72": True,
            "root_number": conductor["root_number"],
            "minimal_discriminant": conductor["minimal_discriminant"],
        },
        "input_selection": {
            "source_stage": SOURCE_STAGE,
            "source_stable_numerical_rank": stage["stable_numerical_rank"],
            "source_height_matrix_runs": stage["height_matrix_runs"],
            "input_point_count": len(input_points),
            "input_point_sha256": point_digest(input_points),
            "numerical_height_rank_used_only_to_select_exact_points": True,
        },
        "small_prime_saturation": saturation,
        "exact_finite_reduction_certificate": certificate,
        "target_assessment": {
            "rank21_log_conductor_target_hit": False,
            "rank30_target_hit": False,
            "reason": "the exact lower bound 16 is below both target ranks",
        },
        "provenance": {
            "script_path": str(script_path.relative_to(root)),
            "script_sha256": sha256_file(script_path),
            "source_artifact": str(args.source.relative_to(root)),
            "source_artifact_sha256": EXPECTED_SOURCE_ARTIFACT_SHA256,
            "source_result_sha256": EXPECTED_SOURCE_RESULT_SHA256,
            "reproducing_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "foreground_capped_saturation_calls": 1,
            "same_stage_retries": 0,
            "owned_processes_remaining": 0,
        },
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    artifact["result_sha256"] = stable_json_digest(
        {
            "theorem": artifact["theorem"],
            "curve": artifact["curve"],
            "input": artifact["input_selection"],
            "saturation": artifact["small_prime_saturation"],
            "certificate": artifact["exact_finite_reduction_certificate"],
            "target": artifact["target_assessment"],
        }
    )
    exclusive_write(args.output, artifact)
    print(
        "certified rank>=16 "
        f"T={PARAMETER} lnN={conductor['log_conductor']} "
        f"basis={EXPECTED_SATURATED_BASIS_SHA256} output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
