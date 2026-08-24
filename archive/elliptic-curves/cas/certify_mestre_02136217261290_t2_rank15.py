#!/usr/bin/env python3
"""Pin the exact rank-at-least-15 certificate found in the max-root-300 census.

The source screen selected fifteen exact points numerically at both 72- and
120-digit height precision.  Numerical rank is used only to select that input
set.  This standalone replay reconstructs the curve directly from its six
roots, checks exact curve membership, and proves independence using the
combined images of the points in ``E(F_p)/3E(F_p)``.  A separate good
reduction excludes rational 3-torsion.  No new parameter or point search is
performed here.
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

from search_mestre_root_tuple_scale import (
    point_digest,
    point_on_short_curve,
    sha256_file,
)
from search_mestre_root_tuple_scale_max100 import stable_json_digest
from search_mestre_root_tuple_scale_max200 import (
    mod3_independence_certificate,
    visible_points_and_coefficients,
)


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

ROOTS = (0, 2, 136, 217, 261, 290)
PARAMETER = Q(2)
IDENTIFIER = "r0_2_136_217_261_290_t2"
SOURCE_ARTIFACT = "elliptic_mestre_root_tuple_scale_max300.json"
EXPECTED_SOURCE_ARTIFACT_SHA256 = (
    "e4a7be774ac0cae3c636c70bde7490e7ced7e313971dfba3c9017e48d730fca7"
)
EXPECTED_SOURCE_RESULT_SHA256 = (
    "a567485ccc11cf56a11ddf83a32c809b5fd2670a820bdcf3f13a754d75b669ce"
)
EXPECTED_SOURCE_SCRIPT_SHA256 = (
    "922cd33621e882dbb5483b041c547f568d3f4fdfea2bffc8cab0d3741a3445b4"
)
EXPECTED_POINT_SHA256 = (
    "43063eff53a2764cda3b950ee257e87cac9f498be3025be582b742ca7303583c"
)
EXPECTED_CERTIFICATE_PRIMES = (
    17,
    23,
    29,
    37,
    41,
    43,
    83,
    101,
    103,
    107,
    131,
    149,
    173,
    193,
    199,
)
EXPECTED_CONDUCTOR = "27535464408096664363840114552671696329422686810"
EXPECTED_LOG_CONDUCTOR = (
    "106.931803973405473954774824140429339998614531372730705726131"
)
EXPECTED_MINIMAL_MODEL = (
    1,
    -1,
    0,
    -26132638263172325109,
    51410336506082459645986038865,
)
EXPECTED_MINIMAL_DISCRIMINANT = (
    "381552631565512125977759238930167772386274815686562212500"
)
CERTIFICATE_PRIME_BOUND = 499
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/"
    "elliptic_mestre_02136217261290_t2_rank15_certificate.json"
)


def load_source(path: Path) -> tuple[dict[str, Any], dict[str, Any], tuple]:
    if sha256_file(path) != EXPECTED_SOURCE_ARTIFACT_SHA256:
        raise AssertionError("the frozen max-root-300 screen artifact changed")
    artifact = json.loads(path.read_text())
    if artifact["result_sha256"] != EXPECTED_SOURCE_RESULT_SHA256:
        raise AssertionError("the frozen max-root-300 result digest changed")
    records = [
        record
        for record in artifact["leader_followup"]["records"]
        if record["identifier"] == IDENTIFIER
    ]
    if len(records) != 1:
        raise AssertionError("the exact T=2 source record changed")
    record = records[0]
    triage = record["point_triage"]
    points = tuple(
        (Q(point["x"]), Q(point["y"]))
        for point in triage["numerical_subset"]
    )
    _, coefficients, _ = visible_points_and_coefficients(ROOTS, PARAMETER)
    if (
        tuple(record["roots"]) != ROOTS
        or Q(record["parameter"]) != PARAMETER
        or triage["status"]
        != "completed exact H5000 point checks and numerical height triage"
        or triage["height_bound"] != 5_000
        or triage["stable_numerical_rank"] != 15
        or len(points) != 15
        or point_digest(points) != EXPECTED_POINT_SHA256
        or any(not point_on_short_curve(coefficients, point) for point in points)
    ):
        raise AssertionError("the exact H5000 rank-15 input set changed")
    conductor = record["conductor_phase"]
    if (
        conductor["conductor"] != EXPECTED_CONDUCTOR
        or conductor["log_conductor"] != EXPECTED_LOG_CONDUCTOR
        or tuple(conductor["minimal_model"]) != EXPECTED_MINIMAL_MODEL
        or conductor["minimal_discriminant"] != EXPECTED_MINIMAL_DISCRIMINANT
        or conductor["root_number"] != -1
    ):
        raise AssertionError("the exact conductor record changed")
    return record, triage, points


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=root / "artifacts/generated-results" / SOURCE_ARTIFACT,
    )
    parser.add_argument("--output", type=Path, default=root / DEFAULT_OUTPUT)
    return parser


def exclusive_write(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")


def main() -> None:
    args = build_parser().parse_args()
    if args.output.exists():
        raise SystemExit("refusing to overwrite the exact rank-15 certificate")
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]
    source_script = script_path.with_name("search_mestre_root_tuple_scale_max300.py")
    if sha256_file(source_script) != EXPECTED_SOURCE_SCRIPT_SHA256:
        raise AssertionError("the frozen max-root-300 screen script changed")

    record, triage, points = load_source(args.source)
    _, coefficients, _ = visible_points_and_coefficients(ROOTS, PARAMETER)
    certificate = mod3_independence_certificate(
        coefficients, points, prime_bound=CERTIFICATE_PRIME_BOUND
    )
    if (
        certificate["status"] != "certified exact algebraic rank lower bound"
        or certificate["certified_algebraic_rank_lower_bound"] != 15
        or certificate["combined_exact_rank_over_F3"] != 15
        or tuple(certificate["certificate_primes"])
        != EXPECTED_CERTIFICATE_PRIMES
        or certificate["point_sha256"] != EXPECTED_POINT_SHA256
        or certificate["rational_3_torsion_exclusion"]
        != {
            "prime": 13,
            "group_order": 20,
            "reason": "rational prime-to-p torsion injects at good reduction",
        }
    ):
        raise AssertionError("the exact mod-3 certificate replay changed")

    conductor = record["conductor_phase"]
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "certified exact algebraic rank lower bound 15",
        "theorem": {
            "statement": (
                "The displayed elliptic curve over Q has Mordell-Weil rank "
                "at least 15."
            ),
            "proof_method": (
                "fifteen exact rational points have full combined rank in "
                "E(F_p)/3E(F_p); rational 3-torsion is excluded at good "
                "reduction p=13"
            ),
            "certified_algebraic_rank_lower_bound": 15,
        },
        "curve": {
            "family_roots": list(ROOTS),
            "parameter_T": "2",
            "sign_equivalent_parameter": "-2",
            "short_weierstrass_coefficients": [str(value) for value in coefficients],
            "minimal_model": conductor["minimal_model"],
            "conductor": conductor["conductor"],
            "log_conductor": conductor["log_conductor"],
            "below_strict_log_conductor_182_72": True,
            "root_number": conductor["root_number"],
            "minimal_discriminant": conductor["minimal_discriminant"],
        },
        "input_selection": {
            "source_height_bound": triage["height_bound"],
            "source_stable_numerical_rank": triage["stable_numerical_rank"],
            "source_height_matrix_runs": triage["height_matrix_runs"],
            "point_count": len(points),
            "point_sha256": point_digest(points),
            "exact_curve_membership_replayed": True,
            "numerical_height_rank_used_only_to_select_exact_points": True,
            "no_new_parameter_or_point_search": True,
        },
        "points": [
            {"x": str(x_value), "y": str(y_value)}
            for x_value, y_value in points
        ],
        "exact_finite_reduction_certificate": certificate,
        "target_assessment": {
            "rank21_log_conductor_target_hit": False,
            "rank30_target_hit": False,
            "reason": "the exact lower bound 15 is below both target ranks",
        },
        "provenance": {
            "script_path": str(script_path.relative_to(root)),
            "script_sha256": sha256_file(script_path),
            "source_artifact": str(args.source.relative_to(root)),
            "source_artifact_sha256": EXPECTED_SOURCE_ARTIFACT_SHA256,
            "source_result_sha256": EXPECTED_SOURCE_RESULT_SHA256,
            "source_script": str(source_script.relative_to(root)),
            "source_script_sha256": EXPECTED_SOURCE_SCRIPT_SHA256,
            "certificate_prime_bound": CERTIFICATE_PRIME_BOUND,
            "reproducing_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "external_process_calls": 0,
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
            "points": artifact["points"],
            "certificate": artifact["exact_finite_reduction_certificate"],
            "target": artifact["target_assessment"],
        }
    )
    exclusive_write(args.output, artifact)
    print(
        "certified rank>=15 "
        f"roots={ROOTS} T=2 lnN={conductor['log_conductor']} "
        f"points={EXPECTED_POINT_SHA256} output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
