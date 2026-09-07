#!/usr/bin/env python3
"""Conditional Delta=11/5 fixed-fibre diagnostics for ICARM #262 and #275.

Both public curves have an unconditional rank lower bound 20.  This script
replays their exact conductors and root numbers with PARI and evaluates the
same Bober sinc-squared explicit formula used for ICARM #245.  Under GRH, an
upper value below 22 together with root number +1 forces analytic rank 20.
BSD is additionally required to identify the algebraic rank.
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
    ROOT/"artifacts/generated-results/elliptic-curves"/
    "icarm_curve262_275_explicit_formula_delta22_v1.json"
)
CURVES = (
    {
        "id": 262,
        "ainvs": (
            1, -1, 1, -21065871080015087031831279377,
            1192838816489664881520195774886398643920001,
        ),
        "conductor": 13935395240740481313503432825112888552210319085746508942961924055710,
    },
    {
        "id": 275,
        "ainvs": (
            1, 0, 1, -2034488389107661074627844285,
            35847670110541831966937994064437784692732,
        ),
        "conductor": 42943483208607336815574462222443765285847682460232958112909965535488718,
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
            raise AssertionError(f"ICARM #{curve['id']} enclosure changed")
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
        "status": "conditional_fixed_fibre_closure_for_two_rank20_curves",
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
            "search_consequence": "move next-point work from these fixed fibres to deformations",
        },
        "software": {"pari_gp": list(pari_version())},
        "reproduction": (
            "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
            "elliptic-curves/cas/explicit_formula_icarm_curve262_275_delta22.py --check"
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
    rendered = json.dumps(payload, indent=2, sort_keys=True)+"\n"
    if arguments.write:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
    elif not arguments.output.exists() or arguments.output.read_text() != rendered:
        raise SystemExit("stale or missing ICARM #262/#275 explicit-formula artifact")
    print(
        "ICARM262275EF|"
        + "|".join(
            f"curve{row['icarm_id']}={row['conditional_explicit_formula_upper']}"
            for row in payload["curves"]
        )
        + f"|mode={'write' if arguments.write else 'check'}|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
