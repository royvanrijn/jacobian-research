#!/usr/bin/env python3
"""Reproduce Fermigier's rank-22 benchmark model, conductor, and score table.

This verifier checks exact curve arithmetic and a historical PARI score
calculation.  It does not independently recheck the 22-point regulator or
prove the rank claim from Fermigier's paper.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import shutil
import subprocess
from typing import Any

from fermigier_mestre import (
    FermigierMestreFamily,
    NORMALIZED_RECORD_PARAMETER,
    PUBLISHED_PARAMETER,
    ROOTS,
)
from pari_bridge import minimal_curve_data, pari_version


PUBLISHED_MODEL = (
    1,
    0,
    1,
    -940299517776391362903023121165864,
    10707363070719743033425295515449274534651125011362,
)
PUBLISHED_CONDUCTOR = int(
    "22720638514787473197194583889675055980109503436060704437972911338086049759883790"
)
EXPECTED_SCORE_ROUNDING = {
    50: "29.49",
    100: "44.12",
    200: "57.54",
    400: "81.51",
    1000: "105.17",
    2000: "122.76",
    4000: "143.84",
    10000: "166.47",
}
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/verify_fermigier_benchmark.py"
)


def fermigier_score_table(timeout: float) -> dict[int, dict[str, Any]]:
    """Return Fermigier's S(E,M), where M indexes primes and p=2 is omitted."""

    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    cutoffs = tuple(EXPECTED_SCORE_ROUNDING)
    conditions = "||".join(f"n=={cutoff}" for cutoff in cutoffs)
    model = ",".join(str(value) for value in PUBLISHED_MODEL)
    program = "\n".join(
        (
            "default(realprecision,60);",
            f"E=ellinit([{model}]);",
            "s=0.;p=2;",
            (
                "for(n=2,10000,p=nextprime(p+1);ap=ellap(E,p);"
                "s+=(2-ap)/(p+1-ap)*log(p);"
                f'if({conditions},print("ROW ",n," ",p," ",s)));'
            ),
            "quit",
        )
    )
    result = subprocess.run(
        [executable, "-q"],
        input=program + "\n",
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if result.returncode != 0 or "***" in result.stderr:
        raise RuntimeError(f"PARI/GP score calculation failed: {result.stderr.strip()}")
    table: dict[int, dict[str, Any]] = {}
    for line in result.stdout.splitlines():
        if not line.startswith("ROW "):
            continue
        _, cutoff, last_prime, value = line.split()
        table[int(cutoff)] = {
            "last_prime": int(last_prime),
            "value": value,
            "rounded_to_two_decimals": f"{float(value):.2f}",
        }
    if tuple(table) != cutoffs:
        raise AssertionError("PARI did not emit the complete pinned score table")
    return table


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic-curves"
        / "elliptic_fermigier_benchmark.json",
    )
    parser.add_argument("--timeout", type=float, default=60.0)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    points = FermigierMestreFamily.known_quartic_points(NORMALIZED_RECORD_PARAMETER)
    if len(points) != 13:
        raise AssertionError("the normalized quartic did not expose all 13 sections")
    jacobian_points = FermigierMestreFamily.known_jacobian_points(
        NORMALIZED_RECORD_PARAMETER
    )
    curve = minimal_curve_data(
        FermigierMestreFamily.coefficients(NORMALIZED_RECORD_PARAMETER),
        timeout=args.timeout,
        known_points=jacobian_points[1:],
    )
    if tuple(curve["minimal_model"]) != PUBLISHED_MODEL:
        raise AssertionError("the normalized family did not reproduce the published model")
    if curve["conductor"] != PUBLISHED_CONDUCTOR:
        raise AssertionError("the computed conductor changed")

    score_table = fermigier_score_table(args.timeout)
    observed_rounding = {
        cutoff: row["rounded_to_two_decimals"] for cutoff, row in score_table.items()
    }
    if observed_rounding != EXPECTED_SCORE_ROUNDING:
        raise AssertionError("the historical score-table rounding did not reproduce")

    artifact = {
        "schema_version": 1,
        "status": (
            "verified computation: exact family/model/conductor and PARI score replay; "
            "the published rank lower bound is cited, not independently certified here"
        ),
        "family": {
            "root_tuple": list(ROOTS),
            "published_parameter": str(PUBLISHED_PARAMETER),
            "internal_normalized_parameter": str(NORMALIZED_RECORD_PARAMETER),
            "visible_quartic_points_checked": len(points),
            "jacobian_images_checked": len(jacobian_points),
            "numerical_rank_seed": (
                "Jacobian images 2--13; the recorded height determinant is "
                "numerical evidence, not an independence certificate"
            ),
        },
        "published_rank_lower_bound": 22,
        "rank_certificate_in_this_artifact": None,
        "curve": curve,
        "strict_target": {"rank_at_least": 21, "log_conductor_less_than": "182.72"},
        "strict_target_met_by_benchmark": False,
        "score_definition": (
            "sum over 2 < p <= p_M of ((2-a_p)/(p+1-a_p))*log(p); "
            "M is the ordinal of p_M, not a numerical prime cutoff"
        ),
        "score_table": {str(key): value for key, value in score_table.items()},
        "software": {"python": platform.python_version(), "pari_gp": pari_version()},
        "sources": [
            "https://matwbn.icm.edu.pl/ksiazki/aa/aa82/aa8243.pdf",
            "https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html",
        ],
        "reproducing_command": REPRODUCING_COMMAND,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    print(f"minimal model reproduced: {tuple(curve['minimal_model']) == PUBLISHED_MODEL}")
    print(f"log(N)={curve['log_conductor']} (strict target is < 182.72)")
    print("historical score table reproduced through M=10000")


if __name__ == "__main__":
    main()
