#!/usr/bin/env python3
"""Replay corrected later moments on sampled SIC(2) (3,3) trace candidates.

The input is produced by
``explore_two_pair_sic_bidegree33_boundary_coefficients.py --trace-norm
--trace-samples``.  This checker extracts the base points where multiplication
by (mu_6,mu_7) drops rank in the rank-twelve (mu_3,mu_4,mu_5) quotient.  It
then reconstructs those finite fibers directly and evaluates the corrected
moments mu_8,...,mu_12,mu_14.

Everything in this file is finite-field computation.  A successful replay is
not a characteristic-zero nullcone certificate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from explore_two_pair_sic_bidegree33_boundary_coefficients import substitute  # noqa: E402
from explore_two_pair_sic_bidegree33_full_anchor import (  # noqa: E402
    chart_expression,
    moment_terms,
    prepare_s0_branch_for_msolve,
)


CORRECTED_ORDERS = tuple(range(2, 13)) + (14,)
LATER_ORDERS = (8, 9, 10, 11, 12, 14)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=47)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--input", type=Path)
    arguments = parser.parse_args()
    prime = arguments.prime
    assert 3 * max(CORRECTED_ORDERS) < prime
    singular = shutil.which("Singular")
    assert singular is not None

    input_path = arguments.input or (
        ROOT
        / "artifacts"
        / "generated-results"
        / f"two_pair_sic_bidegree33_boundary_trace_norm_samples_mod{prime}.json"
    )
    sample_data = json.loads(input_path.read_text())
    assert sample_data["prime"] == prime
    point_order = ("s1", "t0", "Q", "t2")
    points: list[tuple[int, ...]] = []
    for sample in sample_data["samples"]:
        if sample["joint_multiplication_rank"] == sample["quotient_length"]:
            continue
        point = tuple(sample["point"][variable] for variable in point_order)
        if point not in points:
            points.append(point)
    if not points:
        raise AssertionError("the sample artifact contains no Fitting-rank drop")

    expressions = [
        chart_expression(moment_terms(order, prime), 0, prime)
        for order in CORRECTED_ORDERS
    ]
    variables, polynomials = prepare_s0_branch_for_msolve(
        singular,
        expressions,
        prime,
        "s0-boundary",
        arguments.timeout,
    )
    assert variables == ("s1", "s2", "s3", "s5", "t0", "t1", "t2", "t4")
    available = dict(zip(CORRECTED_ORDERS[1:], polynomials))
    adapted = {
        order: substitute(
            polynomial,
            (
                ("t1", "(s1*t0-L)"),
                ("s2", "(s1^2-(13/3)*t0^2-Q)"),
                ("L", "(1)"),
            ),
        )
        for order, polynomial in available.items()
    }

    program = [f"ring base={prime},(s1,t0,Q,t2,s5,t4,s3),dp;"]
    program.extend(
        f"poly a{order}={adapted[order]};"
        for order in CORRECTED_ORDERS[1:]
    )
    for point_index, point in enumerate(points, 1):
        s1, t0, quadratic, t2 = point
        for order in CORRECTED_ORDERS[1:]:
            program.append(
                f"poly q{order}x{point_index}="
                f"subst(subst(subst(subst(a{order},s1,{s1}),t0,{t0}),"
                f"Q,{quadratic}),t2,{t2});"
            )
    program.append(f"ring fiber={prime},(s5,t4,s3),dp;")
    for point_index, _point in enumerate(points, 1):
        program.extend(
            f"poly p{order}x{point_index}=imap(base,q{order}x{point_index});"
            for order in CORRECTED_ORDERS[1:]
        )
        program.append(
            f"ideal Ix{point_index}=p3x{point_index},p4x{point_index},"
            f"p5x{point_index},p6x{point_index},p7x{point_index};"
        )
        program.append(f"ideal Gx{point_index}=std(Ix{point_index});")
        program.append(
            f'print("POINT {point_index} "+string(vdim(Gx{point_index}))'
            f'+" "+string(size(Gx{point_index})));'
        )
        program.append(
            f"int ix{point_index}; "
            f"for(ix{point_index}=1;"
            f"ix{point_index}<=size(Gx{point_index});ix{point_index}++)"
            "{"
            f'print("GB {point_index} "+string(Gx{point_index}[ix{point_index}]));'
            "}"
        )
        for order in LATER_ORDERS:
            program.append(
                f'print("LATER {point_index} {order} "+'
                f"string(reduce(p{order}x{point_index},Gx{point_index})));"
            )

    completed = subprocess.run(
        [singular, "-q"],
        input="\n".join(program),
        text=True,
        capture_output=True,
        check=True,
        timeout=arguments.timeout,
    )
    if "?" in completed.stdout or "overflow" in completed.stdout:
        raise AssertionError(completed.stdout[-8000:])

    results: list[dict[str, object]] = []
    for point_index, point in enumerate(points, 1):
        quotient = re.search(
            rf"(?m)^POINT {point_index} (\d+) (\d+)$",
            completed.stdout,
        )
        if quotient is None:
            raise AssertionError(completed.stdout[-8000:])
        groebner_basis = re.findall(
            rf"(?m)^GB {point_index} (.*)$",
            completed.stdout,
        )
        later = {
            int(order): value
            for order, value in re.findall(
                rf"(?m)^LATER {point_index} (\d+) (.*)$",
                completed.stdout,
            )
        }
        if set(later) != set(LATER_ORDERS):
            raise AssertionError(completed.stdout[-8000:])
        results.append(
            {
                "base_point": dict(zip(point_order, point)),
                "common_quotient_length": int(quotient.group(1)),
                "groebner_basis": groebner_basis,
                "later_moment_normal_forms": {
                    str(order): later[order] for order in LATER_ORDERS
                },
                "killed_by_mu8": later[8] != "0",
            }
        )

    output = (
        ROOT
        / "artifacts"
        / "generated-results"
        / f"two_pair_sic_bidegree33_boundary_trace_candidates_mod{prime}.json"
    )
    payload = {
        "prime": prime,
        "normalization": "s0=1, A=B=mu2=0, L=1",
        "source_samples": str(input_path.relative_to(ROOT)),
        "corrected_moment_set": list(range(1, 13)) + [14],
        "results": results,
        "all_sampled_trace_candidates_killed_by_mu8": all(
            result["killed_by_mu8"] for result in results
        ),
        "scope": (
            "exact finite-field replay at sampled Fitting-rank-drop points; "
            "not a characteristic-zero nullcone certificate"
        ),
        "reproduction_command": (
            ".venv/bin/python "
            "scripts/verify_two_pair_sic_bidegree33_boundary_trace_candidates.py "
            f"--prime {prime} --timeout {arguments.timeout}"
        ),
    }
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"TRACE_CANDIDATES points={len(results)} "
        f"killed_by_mu8={sum(result['killed_by_mu8'] for result in results)}"
    )
    for result in results:
        print(
            f"TRACE_CANDIDATE point={result['base_point']} "
            f"fiber={result['groebner_basis']} "
            f"mu8={result['later_moment_normal_forms']['8']}"
        )
    print(f"TRACE_CANDIDATES_WROTE {output.relative_to(ROOT)}")
    print("PASS: every sampled trace candidate is killed by corrected mu8")


if __name__ == "__main__":
    main()
