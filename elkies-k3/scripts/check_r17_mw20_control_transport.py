#!/usr/bin/env python3
"""Exact rational splitting preflight for R17 quadratic-base controls.

Given ``q(z)=a*z^2+b*z+c`` in one of the eight exact 074d9-lineage
coordinates, evaluate the five certified record controls.  A control has a
rational point above it precisely when ``q(z0)`` is a rational square (zero
means a rational ramification point).  This only tests whether the fibre can
be transported over QQ; testing the surviving quotient requires the three
new generic sections and their exact specializations.
"""

# status: ACTIVE_COMPILER
# claim: exact rational-square splitting preflight for the five native controls
# inputs: artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json
# outputs: caller-selected compact JSON preflight certificate

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from math import isqrt
import json
from pathlib import Path
import shlex
import sys


ROOT = Path(__file__).resolve().parents[2]
LINEAGE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
)
TARGET_IDS = (351, 356, 376, 377, 385)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def rational_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def rational_square_root(value: Fraction):
    if value < 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator * numerator != value.numerator or denominator * denominator != value.denominator:
        return None
    return Fraction(numerator, denominator)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chart")
    parser.add_argument(
        "--q",
        help="comma-separated coefficients c,b,a of q(z)=a*z^2+b*z+c",
    )
    parser.add_argument(
        "--candidate-artifact",
        type=Path,
        help="read the chart and q coefficients from a stored Nagao finalist",
    )
    parser.add_argument("--finalist-rank", type=int, default=1)
    parser.add_argument("--generic-rank", type=int, default=20)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.generic_rank < 17:
        parser.error("--generic-rank must be at least 17")
    extra_inputs = {}
    if args.candidate_artifact is not None:
        if args.chart is not None or args.q is not None:
            parser.error("do not combine --candidate-artifact with --chart or --q")
        candidate_path = args.candidate_artifact
        candidate_payload = json.loads(candidate_path.read_text())
        finalists = {
            int(record["rank"]): record for record in candidate_payload["finalists"]
        }
        if args.finalist_rank not in finalists:
            parser.error("--finalist-rank is absent from --candidate-artifact")
        finalist = finalists[args.finalist_rank]
        chart = finalist["lineage_coordinate"]
        coefficients = tuple(
            Fraction(value)
            for value in finalist["q_in_lineage_coordinate_coefficients_low_to_high"]
        )
        try:
            candidate_key = str(candidate_path.resolve().relative_to(ROOT))
        except ValueError:
            candidate_key = str(candidate_path.resolve())
        extra_inputs[candidate_key] = digest(candidate_path)
    else:
        if args.chart is None or args.q is None:
            parser.error("provide --candidate-artifact or both --chart and --q")
        chart = args.chart
        coefficients = tuple(Fraction(value) for value in args.q.split(","))
    if len(coefficients) != 3:
        parser.error("--q must contain exactly c,b,a")
    c, b, a = coefficients
    if b * b == 4 * a * c:
        parser.error("the quadratic character is singular")

    lineage = json.loads(LINEAGE.read_text())
    if lineage.get("status") != "PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS":
        raise ValueError("the exact lineage certificate is unavailable")
    if chart not in lineage["chart_transports"]:
        parser.error("--chart is not in the exact published-R17 lineage")
    quotients = {
        int(record["curve_id"]): record for record in lineage["exceptional_quotients"]
    }
    parameters = {
        int(record["curve_id"]): Fraction(record["parameter"])
        for record in lineage["target_isomorphisms"]
        if record["chart"] == chart and int(record["curve_id"]) in TARGET_IDS
    }
    if tuple(sorted(parameters)) != TARGET_IDS:
        raise ValueError("the selected chart does not contain all five exact controls")

    rows = []
    for curve_id in TARGET_IDS:
        parameter = parameters[curve_id]
        value = a * parameter * parameter + b * parameter + c
        root = rational_square_root(value)
        old_tail = int(quotients[curve_id]["free_rank"])
        fibre_rank_lower_bound = 17 + old_tail
        if root is None:
            splitting = "NO_RATIONAL_PREIMAGE"
            preimages = []
        elif root == 0:
            splitting = "RATIONAL_RAMIFICATION_POINT"
            preimages = ["0/1"]
        else:
            splitting = "TWO_RATIONAL_PREIMAGES"
            preimages = [rational_text(root), rational_text(-root)]
        rows.append(
            {
                "curve_id": curve_id,
                "parameter": rational_text(parameter),
                "q_value": rational_text(value),
                "splitting": splitting,
                "cover_coordinate_values": preimages,
                "certified_fibre_rank_lower_bound": fibre_rank_lower_bound,
                "old_tail_beyond_rank_17": old_tail,
                "reference_rank_budget_beyond_assumed_generic_rank": max(
                    0, fibre_rank_lower_bound - args.generic_rank
                ),
                "rank_budget_on_this_rational_lift": (
                    None
                    if root is None
                    else max(0, fibre_rank_lower_bound - args.generic_rank)
                ),
                "exact_new_tail_status": (
                    "PENDING_SPECIALIZATION_OF_NEW_GENERIC_SECTIONS"
                    if root is not None
                    else "NOT_AVAILABLE_OVER_Q_ON_THIS_BASE_CHANGE"
                ),
            }
        )

    payload = {
        "schema": "elkies-k3.r17-mw20-control-transport-preflight.v1",
        "status": "PASS_EXACT_CONTROL_SPLITTING_PREFLIGHT",
        "character": {
            "chart": chart,
            "coefficients_c_b_a": [rational_text(value) for value in coefficients],
            "formula": "q(z)=a*z^2+b*z+c",
        },
        "assumed_generic_rank": args.generic_rank,
        "controls": rows,
        "alternate_q80_controls": {
            "status": "NOT_APPLICABLE_ON_THE_DISPLAYED_BASE_CHANGE",
            "reason": (
                "The alternate-Q80 controls lie in different rational-PGL2 j-map "
                "classes, not at base parameters of the published-R17/074d9 "
                "fibration. They require a separate explicit common-K3 birational "
                "transport before a pullback test is defined."
            ),
        },
        "inputs": {
            str(LINEAGE.relative_to(ROOT)): digest(LINEAGE),
            **extra_inputs,
        },
        "reproducing_command": shlex.join(
            argument for argument in sys.argv if argument != "--check"
        ),
        "proof_boundary": (
            "The square tests exactly determine rational points of the quadratic "
            "cover above the five controls. The rank budgets use only certified "
            "fibre lower bounds. They do not certify that any of the old quotient "
            "basis remains independent modulo the new generic sections; that final "
            "tail test requires explicit specialization coordinates for all three "
            "new twist sections."
        ),
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text() != encoded:
            raise SystemExit("stored artifact differs from replay")
        print(f"PASS check {args.output}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded)
    split_count = sum(
        row["splitting"] != "NO_RATIONAL_PREIMAGE" for row in rows
    )
    print(
        f"PASS exact control splitting preflight chart={chart} "
        f"split={split_count}/5 output={args.output}"
    )


if __name__ == "__main__":
    main()
