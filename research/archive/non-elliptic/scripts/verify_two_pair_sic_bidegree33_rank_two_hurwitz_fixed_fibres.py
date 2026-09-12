#!/usr/bin/env python3
"""Exclude three additional exact fibres of the cubic Hurwitz chart.

On

    F=(1+x)B(y)+x^2*(lambda+x)D(y),

normalize B(0)=1, eliminate d3 with mu_1, and localize the quadratic
discriminant and channel minor M_01.  Reconstruct the exact integer
moments of 3F from the provably large Mersenne prime used by the generic
Hurwitz explorer.  Exact QQ msolve calculations prove that
mu_2,...,mu_8 generate the unit ideal on lambda=-1, 1, and 2.

These are complete fibres on the declared localization.  They do not
classify the generic lambda-line, its localization boundaries, other
channel charts, or exceptional binary-cubic pencils.
"""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from research_two_pair_sic_bidegree33_rank_two_hurwitz import (
    EXACT_LIFT_PRIME,
    exact_moment_polynomials,
    exact_msolve_source,
    exact_polynomial_string,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_hurwitz_fixed_fibres.json"
)
ORDERS = tuple(range(2, 9))
LAMBDA_VALUES = (-1, 1, 2)


def solve_exact(source: str, timeout: int) -> dict[str, object]:
    executable = shutil.which("msolve")
    if executable is None:
        raise RuntimeError("msolve is required")
    with tempfile.TemporaryDirectory(prefix="sic33-hurwitz-fixed-") as directory:
        input_path = Path(directory) / "system.ms"
        output_path = Path(directory) / "result.ms"
        input_path.write_text(source, encoding="utf-8")
        completed = subprocess.run(
            [
                executable,
                "-f",
                str(input_path),
                "-o",
                str(output_path),
                "-t",
                "4",
                "-l",
                "2",
                "-v",
                "0",
            ],
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
        result = output_path.read_text(encoding="utf-8") if output_path.exists() else ""
    compact = " ".join(result.split())
    if completed.returncode != 0 or compact not in {"[-1]", "[-1]:"}:
        raise AssertionError(
            "exact Hurwitz fibre was not a unit ideal: "
            + json.dumps(
                {
                    "returncode": completed.returncode,
                    "result_head": result[:2000],
                    "stderr_tail": completed.stderr[-2000:],
                },
                sort_keys=True,
            )
        )
    return {
        "returncode": completed.returncode,
        "status": "unit_ideal",
        "result_sha256": sha256(compact.encode()).hexdigest(),
        "result_bytes": len(result.encode()),
    }


def main() -> None:
    moments = exact_moment_polynomials(ORDERS)
    profiles = {
        str(order): {
            "terms": len(moments[order]),
            "sha256": sha256(
                exact_polynomial_string(moments[order]).encode()
            ).hexdigest(),
        }
        for order in ORDERS
    }
    records = []
    for lambda_value in LAMBDA_VALUES:
        source = exact_msolve_source(moments, "01", lambda_value)
        records.append(
            {
                "lambda": lambda_value,
                "msolve_input_sha256": sha256(source.encode()).hexdigest(),
                "certificate": solve_exact(source, 300),
            }
        )
        print(f"certified lambda={lambda_value}", flush=True)

    artifact = {
        "format": "two-pair-sic-bidegree33-rank-two-hurwitz-fixed-fibres-v1",
        "field": "characteristic zero",
        "chart": "F=(1+x)B(y)+x^2*(lambda+x)D(y)",
        "normalization": "B(0)=1 and d3 eliminated by mu_1",
        "localization": "quadratic discriminant Delta and channel minor M_01",
        "orders": list(ORDERS),
        "lambda_values": list(LAMBDA_VALUES),
        "scaled_form": "3F",
        "exact_reconstruction": {
            "prime": str(EXACT_LIFT_PRIME),
            "coefficient_bound": "(3m)!*52^m",
            "prime_exceeds_twice_bound": True,
        },
        "moment_profiles": profiles,
        "fibres": records,
        "conclusion": (
            "the complete declared lambda=-1,1,2 fibres have no common "
            "zero of mu_2 through mu_8"
        ),
        "scope": (
            "three fixed fibres on one localized generic channel chart; "
            "the generic lambda-line, localization boundaries, other channel "
            "charts, and exceptional cubic pencils remain open"
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")

    print("PASS exact QQ moments reconstructed for the cubic Hurwitz chart")
    print("PASS lambda=-1,1,2 localized fibres are unit ideals through mu_8")


if __name__ == "__main__":
    main()
