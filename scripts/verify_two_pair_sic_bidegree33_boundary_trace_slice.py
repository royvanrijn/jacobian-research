#!/usr/bin/env python3
"""Exhaust a fixed two-parameter slice of the SIC(2) (3,3) boundary.

After the exact pivots and the weighted normalization L=1, this checker fixes
``s1`` and ``t0``.  It computes the complete zero-dimensional quotient of
mu_3,...,mu_7 in the five remaining variables and then adjoins corrected
mu_8.  Results at several good primes are modular evidence, not a
characteristic-zero certificate.
"""

from __future__ import annotations

import argparse
import hashlib
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


ORDERS = tuple(range(2, 9))


def compute_slice(
    singular: str,
    prime: int,
    s1_value: int,
    t0_value: int,
    timeout: int,
) -> dict[str, object]:
    expressions = [
        chart_expression(moment_terms(order, prime), 0, prime)
        for order in ORDERS
    ]
    variables, polynomials = prepare_s0_branch_for_msolve(
        singular,
        expressions,
        prime,
        "s0-boundary",
        timeout,
    )
    assert variables == ("s1", "s2", "s3", "s5", "t0", "t1", "t2", "t4")
    available = dict(zip(ORDERS[1:], polynomials))
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
    declarations = "\n".join(
        (
            f"poly a{order}={adapted[order]};\n"
            f"poly q{order}=subst(subst("
            f"a{order},s1,{s1_value}),t0,{t0_value});"
        )
        for order in ORDERS[1:]
    )
    fetch = "\n".join(
        f"poly p{order}=imap(base,q{order});" for order in ORDERS[1:]
    )
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring base={prime},(s1,t0,Q,t2,s5,t4,s3),dp;
{declarations}
ring slice={prime},(s5,t4,s3,Q,t2),dp;
option(redSB);
{fetch}
ideal I7=p3,p4,p5,p6,p7;
ideal G7=slimgb(I7);
print(
  "SLICE7 "+string(dim(G7))+" "+string(size(G7))+" "
  +string(vdim(G7))+" "+string(G7[1]==1)
);
int basisIndex;
for (basisIndex=1;basisIndex<=size(G7);basisIndex++)
{{
  print("SLICE7_LEADING "+string(leadexp(G7[basisIndex])));
}}
ideal I8=G7,p8;
ideal G8=slimgb(I8);
print(
  "SLICE8 "+string(dim(G8))+" "+string(size(G8))+" "
  +string(vdim(G8))+" "+string(G8[1]==1)
);
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    if "?" in completed.stdout or "overflow" in completed.stdout:
        raise AssertionError(completed.stdout[-8000:])
    seventh = re.search(
        r"(?m)^SLICE7 (-?\d+) (\d+) (-?\d+) ([01])$",
        completed.stdout,
    )
    eighth = re.search(
        r"(?m)^SLICE8 (-?\d+) (\d+) (-?\d+) ([01])$",
        completed.stdout,
    )
    if seventh is None or eighth is None:
        raise AssertionError(completed.stdout[-8000:])
    leading = re.findall(
        r"(?m)^SLICE7_LEADING ([0-9,]+)$",
        completed.stdout,
    )
    if len(leading) != int(seventh.group(2)):
        raise AssertionError("incomplete leading-ideal export")
    leading_text = "\n".join(leading)
    return {
        "prime": prime,
        "through_mu7": {
            "dimension": int(seventh.group(1)),
            "groebner_basis_size": int(seventh.group(2)),
            "quotient_length": int(seventh.group(3)),
            "unit_ideal": bool(int(seventh.group(4))),
            "leading_exponents_s5_t4_s3_Q_t2": [
                [int(value) for value in exponent.split(",")]
                for exponent in leading
            ],
            "leading_exponents_sha256": hashlib.sha256(
                leading_text.encode()
            ).hexdigest(),
        },
        "through_mu8": {
            "dimension": int(eighth.group(1)),
            "groebner_basis_size": int(eighth.group(2)),
            "quotient_length": int(eighth.group(3)),
            "unit_ideal": bool(int(eighth.group(4))),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", default="47,101")
    parser.add_argument("--s1", type=int, default=1)
    parser.add_argument("--t0", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=600)
    arguments = parser.parse_args()
    primes = tuple(int(value) for value in arguments.primes.split(",") if value)
    assert primes and all(3 * max(ORDERS) < prime for prime in primes)
    singular = shutil.which("Singular")
    assert singular is not None

    results = [
        compute_slice(
            singular,
            prime,
            arguments.s1 % prime,
            arguments.t0 % prime,
            arguments.timeout,
        )
        for prime in primes
    ]
    assert all(result["through_mu7"]["quotient_length"] == 1128 for result in results)
    assert all(result["through_mu8"]["unit_ideal"] for result in results)

    output = (
        ROOT
        / "artifacts"
        / "generated-results"
        / (
            "two_pair_sic_bidegree33_boundary_trace_slice_"
            f"s1_{arguments.s1}_t0_{arguments.t0}.json"
        )
    )
    payload = {
        "normalization": "s0=1, A=B=mu2=0, L=1",
        "fixed_slice": {"s1": arguments.s1, "t0": arguments.t0},
        "corrected_moments_used": list(range(3, 9)),
        "results": results,
        "scope": (
            "complete finite-field slice computation; not a global or "
            "characteristic-zero nullcone certificate"
        ),
        "reproduction_command": (
            ".venv/bin/python "
            "scripts/verify_two_pair_sic_bidegree33_boundary_trace_slice.py "
            f"--primes {','.join(str(prime) for prime in primes)} "
            f"--s1 {arguments.s1} --t0 {arguments.t0} "
            f"--timeout {arguments.timeout}"
        ),
    }
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    for result in results:
        print(
            f"TRACE_SLICE prime={result['prime']} "
            f"length7={result['through_mu7']['quotient_length']} "
            f"basis7={result['through_mu7']['groebner_basis_size']} "
            f"unit8={int(result['through_mu8']['unit_ideal'])}"
        )
    print(f"TRACE_SLICE_WROTE {output.relative_to(ROOT)}")
    print("PASS: corrected mu8 kills the complete sampled slice at every prime")


if __name__ == "__main__":
    main()
