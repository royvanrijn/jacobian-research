#!/usr/bin/env python3
"""Delta=11/5 explicit-formula diagnostic for the rank-18 T=6629/174 fiber.

The prime-power sum contains every term in Bober's sinc-squared test function
through ``floor(exp(22*pi/5))``.  The common archimedean term uses the closed
dilogarithm identity already independently audited by the section-7
diagnostic.  A deliberately conservative 0.001 allowance is added.

Under GRH the resulting value bounds analytic rank.  Root number +1 then
forces an even analytic order, while identifying the exact algebraic rank
also requires BSD.  The rank-at-least-18 finite-reduction certificate itself
is unconditional.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shlex
import subprocess
import sys
from typing import Any

import mpmath

from ek_k3 import rational_to_string
from explicit_formula_rank20_t5081_delta22 import (
    EXPECTED_PRIME_LIMIT,
    archimedean_closed_form,
    explicit_formula_upper,
    gp_program,
    parse_prime_sum,
    prime_limit,
)
from explicit_formula_rank_diagnostic import BOBER_SOURCE
from pari_bridge import pari_version


Q = Fraction
REPOSITORY = Path(__file__).resolve().parents[2]
PARAMETER_T = Q(6629, 174)
DELTA = Q(11, 5)
NUMERICAL_ALLOWANCE = mpmath.mpf("0.001")
INPUT_RANK = 18
INPUT_ARTIFACT = (
    REPOSITORY
    / "artifacts/generated-results/elliptic_nagao_rank21_historical_finalists.json"
)
INPUT_ARTIFACT_SHA256 = "90fc658cdb7c39c96317ee888be1364b8c9f368859230e25161dc45cd6a3cec7"
FORMULA_HELPER = REPOSITORY / "elliptic-curves/cas/explicit_formula_rank20_t5081_delta22.py"
FORMULA_HELPER_SHA256 = "c33da285159adb34bd14dcd2caf5ebdd77bddf5454d0eb226ca44da8543a80b7"
EXPECTED_CONDUCTOR = 16852483164717580988499792436861320283505875689469771069679157006330
EXPECTED_LOG_CONDUCTOR = (
    "154.795114152373636353692290456048113196833306511053393251152"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/explicit_formula_rank18_t6629_delta22.py"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_checkpoint() -> tuple[dict[str, Any], tuple[int, ...]]:
    if sha256_file(INPUT_ARTIFACT) != INPUT_ARTIFACT_SHA256:
        raise AssertionError("the pinned rank-18 input artifact changed")
    if sha256_file(FORMULA_HELPER) != FORMULA_HELPER_SHA256:
        raise AssertionError("the audited explicit-formula helper changed")
    data = json.loads(INPUT_ARTIFACT.read_text(encoding="utf-8"))
    matches = [
        record
        for record in data["exact_checkpoints_stable_numerical_rank_at_least_18"]
        if Q(record["constructor_parameter"]) == PARAMETER_T
    ]
    if len(matches) != 1:
        raise AssertionError("the T=6629/174 checkpoint was not unique")
    checkpoint = matches[0]
    certificate = checkpoint["exact_rank_certificate"]
    conductor = checkpoint["conductor"]
    if (
        certificate["status"] != "certified"
        or certificate["certified_algebraic_rank_lower_bound"] != INPUT_RANK
        or certificate["combined_exact_rank_over_F2"] != INPUT_RANK
    ):
        raise AssertionError("the exact algebraic lower bound changed")
    if (
        conductor["conductor"] != EXPECTED_CONDUCTOR
        or conductor["log_conductor"] != EXPECTED_LOG_CONDUCTOR
        or conductor["root_number"] != 1
        or not conductor["below_strict_log_conductor_target"]
    ):
        raise AssertionError("the conductor/root checkpoint changed")
    model = tuple(int(value) for value in conductor["minimal_model"])
    if len(model) != 5:
        raise AssertionError("the minimal model changed dimension")
    return checkpoint, model


def run_prime_sum(
    model: tuple[int, ...], *, timeout: float, stack_bytes: int
) -> dict[str, str | int]:
    if not 0 < timeout <= 120 or stack_bytes < 64_000_000:
        raise ValueError("invalid PARI resource bounds")
    support = prime_limit(DELTA)
    if support != EXPECTED_PRIME_LIMIT:
        raise AssertionError("the Delta=11/5 support cutoff changed")
    result = subprocess.run(
        ["gp", "-q", "-s", str(stack_bytes)],
        input=gp_program(model, delta=DELTA, support_limit=support),
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if result.returncode != 0 or "***" in result.stderr:
        raise RuntimeError(f"PARI prime sum failed: {result.stderr.strip()}")
    record = parse_prime_sum(result.stdout)
    if int(record["conductor"]) != EXPECTED_CONDUCTOR:
        raise AssertionError("PARI conductor replay changed")
    return record


def build_artifact(args: argparse.Namespace) -> dict[str, Any]:
    checkpoint, model = load_checkpoint()
    prime_sum = run_prime_sum(
        model,
        timeout=args.timeout,
        stack_bytes=args.stack_bytes,
    )
    with mpmath.workdps(90):
        closed_integral, archimedean_contribution = archimedean_closed_form(DELTA)
        conservative_upper = explicit_formula_upper(
            log_conductor=str(prime_sum["log_conductor"]),
            prime_sum=str(prime_sum["prime_sum"]),
            delta=DELTA,
            archimedean_contribution=archimedean_contribution,
            numerical_allowance=NUMERICAL_ALLOWANCE,
        )
        unpadded_value = conservative_upper - NUMERICAL_ALLOWANCE
        if not conservative_upper < 20:
            raise AssertionError("the conservative explicit-formula value is no longer below 20")
        value_record = {
            "closed_form_full_archimedean_integral": mpmath.nstr(closed_integral, 80),
            "closed_form_archimedean_contribution": mpmath.nstr(archimedean_contribution, 80),
            "unpadded_explicit_formula_value": mpmath.nstr(unpadded_value, 80),
            "numerical_allowance": mpmath.nstr(NUMERICAL_ALLOWANCE, 20),
            "conservative_upper_value": mpmath.nstr(conservative_upper, 80),
            "conservative_upper_value_less_than_20": True,
        }
    return {
        "schema_version": 1,
        "status": "conditional explicit-formula diagnostic complete",
        "candidate": {
            "constructor_parameter": rational_to_string(PARAMETER_T),
            "minimal_model": list(model),
            "conductor": EXPECTED_CONDUCTOR,
            "log_conductor": checkpoint["conductor"]["log_conductor"],
            "root_number": 1,
            "unconditional_certified_algebraic_rank_lower_bound": INPUT_RANK,
        },
        "method": {
            "test_function": "sinc-squared",
            "delta": rational_to_string(DELTA),
            "prime_limit": EXPECTED_PRIME_LIMIT,
            "prime_limit_definition": "floor(exp(2*pi*Delta))",
            "source": BOBER_SOURCE,
            "closed_archimedean_identity": (
                "A_D=(pi*D)^-1[-EulerGamma+"
                "(pi^2/6-Li_2(exp(-2*pi*D)))/(2*pi*D)]"
            ),
            "formula_helper_path": str(FORMULA_HELPER),
            "formula_helper_sha256": FORMULA_HELPER_SHA256,
        },
        "prime_power_sum": prime_sum,
        "explicit_formula": value_record,
        "input": {
            "rank18_artifact_path": str(INPUT_ARTIFACT),
            "rank18_artifact_sha256": INPUT_ARTIFACT_SHA256,
            "saturated_basis_sha256": checkpoint["exact_rank_certificate"]["saturated_point_sha256"],
        },
        "interpretation": {
            "under_grh": (
                "the explicit-formula value is below 20, hence analytic rank is "
                "at most 19; root number +1 forces even analytic order, hence at most 18"
            ),
            "under_bsd_and_grh": (
                "the exact algebraic lower bound 18 and analytic upper bound 18 "
                "identify algebraic and analytic rank exactly 18"
            ),
            "unconditional": (
                "the finite-reduction certificate proves algebraic rank at least "
                "18; this computation supplies no unconditional rank upper bound"
            ),
        },
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "mpmath": mpmath.__version__,
        },
        "declared_budget": {
            "timeout_seconds": args.timeout,
            "pari_stack_bytes": args.stack_bytes,
            "all_prime_powers_in_support_included": True,
        },
        "reproducing_command": REPRODUCING_COMMAND,
        "actual_command": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "script_sha256": sha256_file(Path(__file__).resolve()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            REPOSITORY
            / "artifacts/generated-results/elliptic_nagao_rank18_t6629_explicit_formula_delta22.json"
        ),
    )
    args = parser.parse_args()
    artifact = build_artifact(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"wrote {args.output}: upper="
        f"{artifact['explicit_formula']['conservative_upper_value']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
