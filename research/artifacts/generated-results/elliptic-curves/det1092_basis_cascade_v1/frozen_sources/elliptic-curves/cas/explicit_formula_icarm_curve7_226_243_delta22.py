#!/usr/bin/env python3
"""Conditional Delta=11/5 diagnostics for the remaining low-conductor ICARM rank-20 curves.

The three public curves below complete the list of ICARM rank-at-least-20
curves with log conductor below 182.72 that was not already audited in the
curve-245 and curve-262/275 scripts.  The computation is conditional exactly
as there: under GRH an explicit-formula upper bound below 22, together with
root number +1 and the unconditional rank lower bound 20, forces analytic
rank 20.  BSD is additionally needed to identify the algebraic rank.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import subprocess

import mpmath

from explicit_formula_rank20_t5081_delta22 import (
    EXPECTED_PRIME_LIMIT,
    archimedean_closed_form,
    explicit_formula_upper,
    gp_program,
    parse_prime_sum,
)
from pari_bridge import pari_version


Q = Fraction
DELTA = Q(11, 5)
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/elliptic-curves"
    / "icarm_curve7_226_243_explicit_formula_delta22_v1.json"
)
CURVES = (
    {
        "id": 243,
        "ainvs": (
            1,
            0,
            1,
            -791198747812844165197303241658,
            291678735985274857428612896086571996361540568,
        ),
        "conductor": 123989445321322909065995690886185865732040186216447789402341303782551770,
    },
    {
        "id": 226,
        "ainvs": (
            1,
            0,
            0,
            -53324817965388276805370879748910,
            152600569230942786227963554343099291511940660100,
        ),
        "conductor": 282182911611488014589765703157855796043317619561544915217487973922307930,
    },
    {
        "id": 7,
        "ainvs": (
            1,
            0,
            0,
            -431092980766333677958362095891166,
            5156283555366643659035652799871176909391533088196,
        ),
        "conductor": 73813242020125452593092920037241715783942979322145629022774016712292126930,
    },
)


def run_gp(program: str, *, timeout: float, stack_bytes: int) -> str:
    process = subprocess.run(
        ["gp", "-q", "-s", str(stack_bytes)],
        input=program,
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if process.returncode != 0 or "***" in process.stderr:
        raise RuntimeError(f"PARI/GP failed: {process.stderr.strip()}")
    return process.stdout


def curve_record(curve: dict, *, timeout: float, stack_bytes: int) -> dict:
    model = tuple(map(int, curve["ainvs"]))
    prime_output = run_gp(
        gp_program(model, delta=DELTA, support_limit=EXPECTED_PRIME_LIMIT),
        timeout=timeout,
        stack_bytes=stack_bytes,
    )
    prime_record = parse_prime_sum(prime_output)
    prime_record.pop("pari_milliseconds_after_curve_setup", None)
    if int(prime_record["conductor"]) != curve["conductor"]:
        raise AssertionError(f"ICARM #{curve['id']} conductor changed")
    vector = ",".join(map(str, model))
    root_output = run_gp(
        f'E=ellinit([{vector}]);print("ROOT|",ellrootno(E));quit\n',
        timeout=min(timeout, 30),
        stack_bytes=stack_bytes,
    )
    if [row for row in root_output.splitlines() if row.startswith("ROOT|")] != [
        "ROOT|1"
    ]:
        raise AssertionError(f"ICARM #{curve['id']} root number changed")
    with mpmath.workdps(60):
        archimedean = archimedean_closed_form(DELTA)[1]
        upper = explicit_formula_upper(
            log_conductor=str(prime_record["log_conductor"]),
            prime_sum=str(prime_record["prime_sum"]),
            delta=DELTA,
            archimedean_contribution=archimedean,
        )
        if not mpmath.mpf(20) < upper < mpmath.mpf(22):
            raise AssertionError(f"ICARM #{curve['id']} enclosure changed: {upper}")
        upper_text = mpmath.nstr(upper, 60)
    return {
        "icarm_id": curve["id"],
        "source": f"https://elliptic-rank.icarm.cloud/curve/{curve['id']}",
        "ainvs": list(map(str, model)),
        "rank_lower_bound": 20,
        "conductor": str(curve["conductor"]),
        "exact_root_number": 1,
        "prime_sum": prime_record,
        "conditional_explicit_formula_upper": upper_text,
        "strictly_below_22": True,
    }


def build_payload(*, timeout: float, stack_bytes: int) -> dict:
    return {
        "schema_version": 1,
        "artifact_kind": "conditional_explicit_formula_rank_diagnostics",
        "status": "conditional_fixed_fibre_closure_for_three_rank20_curves",
        "test_function": "Bober sinc-squared",
        "delta": "11/5",
        "curves": [
            curve_record(curve, timeout=timeout, stack_bytes=stack_bytes)
            for curve in CURVES
        ],
        "conclusion": {
            "unconditional": "each displayed curve has rank at least 20 only",
            "under_grh": "each analytic rank is 20",
            "under_grh_and_bsd": "each algebraic rank is 20",
            "search_consequence": "all seven public sub-cutoff rank-20 fibres now have even root number; prioritize deformation rather than fixed-fibre point search",
        },
        "software": {"pari_gp": list(pari_version())},
        "reproduction": (
            "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
            "elliptic-curves/cas/explicit_formula_icarm_curve7_226_243_delta22.py --check"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    payload = build_payload(
        timeout=arguments.timeout, stack_bytes=arguments.stack_bytes
    )
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.write:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
    elif not arguments.output.exists() or arguments.output.read_text() != rendered:
        raise SystemExit("stale or missing ICARM #7/#226/#243 explicit-formula artifact")
    print(
        "ICARM7226243EF|"
        + "|".join(
            f"curve{row['icarm_id']}={row['conditional_explicit_formula_upper']}"
            for row in payload["curves"]
        )
        + f"|mode={'write' if arguments.write else 'check'}|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
