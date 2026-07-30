#!/usr/bin/env python3
"""Exhaust a fixed two-parameter slice of the SIC(2) (3,3) boundary.

After the exact pivots and the weighted normalization L=1, this checker fixes
``s1`` and ``t0``.  It computes the complete zero-dimensional quotient of
mu_3,...,mu_7 in the five remaining variables and then adjoins corrected
mu_8.  It compares several good primes and uses verified modular
reconstruction to certify the corresponding characteristic-zero unit ideal.
The result is an exact slice certificate, not a global nullcone certificate.
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
import tempfile


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from explore_two_pair_sic_bidegree33_boundary_coefficients import substitute  # noqa: E402
from explore_two_pair_sic_bidegree33_full_anchor import (  # noqa: E402
    chart_expression,
    moment_terms,
    prepare_s0_branch_for_msolve,
)
from verify_two_pair_sic_bidegree33_boundary_generic_quotient import (  # noqa: E402
    exact_chart_expression,
    exact_moment_terms,
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


def compute_exact_slice(
    singular: str,
    s1_value: int,
    t0_value: int,
    timeout: int,
) -> dict[str, object]:
    """Prove the corrected mu_3,...,mu_8 slice ideal is the unit ideal."""

    expressions = [
        exact_chart_expression(exact_moment_terms(order))
        for order in ORDERS
    ]
    variables, polynomials = prepare_s0_branch_for_msolve(
        singular,
        expressions,
        0,
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
LIB "modstd.lib";
ring base=0,(s1,t0,Q,t2,s5,t4,s3),dp;
{declarations}
ring slice=0,(s5,t4,s3,Q,t2),dp;
option(redSB);
{fetch}
ideal I8=p3,p4,p5,p6,p7,p8;
ideal G8=modStd(I8,1);
print(
  "EXACT_SLICE8 "+string(dim(G8))+" "+string(size(G8))+" "
  +string(vdim(G8))+" "+string(G8[1]==1)
);
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    if "?" in completed.stdout:
        raise AssertionError(completed.stdout[-8000:])
    marker = re.search(
        r"(?m)^EXACT_SLICE8 (-?\d+) (\d+) (-?\d+) ([01])$",
        completed.stdout,
    )
    if marker is None:
        raise AssertionError(completed.stdout[-8000:])
    return {
        "characteristic": 0,
        "algorithm": "Singular modStd with exactness=1",
        "through_mu8": {
            "dimension": int(marker.group(1)),
            "groebner_basis_size": int(marker.group(2)),
            "quotient_length": int(marker.group(3)),
            "unit_ideal": bool(int(marker.group(4))),
        },
    }


def compute_modular_hyperslice(
    singular: str,
    msolve: str,
    prime: int,
    fixed_variable: str,
    fixed_value: int,
    timeout: int,
) -> dict[str, object]:
    """Test one coordinate-fixed hyperslice through corrected mu_8."""

    assert fixed_variable in ("s1", "t0")
    free_variable = "t0" if fixed_variable == "s1" else "s1"
    expressions = [
        chart_expression(moment_terms(order, prime), 0, prime)
        for order in ORDERS
    ]
    _variables, polynomials = prepare_s0_branch_for_msolve(
        singular,
        expressions,
        prime,
        "s0-boundary",
        timeout,
    )
    replacements = (
        ("t1", "(s1*t0-L)"),
        ("s2", "(s1^2-(13/3)*t0^2-Q)"),
        ("L", "(1)"),
        (fixed_variable, f"({fixed_value})"),
    )
    generators = [
        substitute(polynomial, replacements) for polynomial in polynomials
    ]
    with tempfile.TemporaryDirectory(prefix="sic33-trace-hyperslice-") as directory:
        input_path = Path(directory) / "hyperslice.ms"
        output_path = Path(directory) / "hyperslice.out"
        input_path.write_text(
            f"s5,t4,s3,Q,t2,{free_variable}\n"
            f"{prime}\n"
            + ",\n".join(generators)
            + "\n"
        )
        completed = subprocess.run(
            [
                msolve,
                "-f",
                str(input_path),
                "-o",
                str(output_path),
                "-t",
                "4",
                "-v",
                "1",
                "-l",
                "44",
                "-g",
                "1",
            ],
            text=True,
            capture_output=True,
            check=False,
            timeout=timeout,
        )
        if completed.returncode != 0:
            raise AssertionError(
                f"msolve exit {completed.returncode}\n"
                f"{completed.stdout[-2000:]}\n{completed.stderr[-2000:]}"
            )
        result = output_path.read_text()
    return {
        "prime": prime,
        "fixed_coordinate": {fixed_variable: fixed_value},
        "free_variables": [free_variable, "Q", "t2", "s3", "s5", "t4"],
        "through_mu8_unit_ideal": result.rstrip().endswith("[1]:"),
        "basis_size": 1 if "#length of basis:      1 element" in result else None,
        "backend": "msolve linear algebra 44",
    }


def compute_exact_hyperslice(
    singular: str,
    msolve: str,
    fixed_variable: str,
    fixed_value: int,
    timeout: int,
) -> dict[str, object]:
    """Certify one coordinate-fixed hyperslice over the rationals."""

    assert fixed_variable in ("s1", "t0")
    free_variable = "t0" if fixed_variable == "s1" else "s1"
    expressions = [
        exact_chart_expression(exact_moment_terms(order))
        for order in ORDERS
    ]
    _variables, polynomials = prepare_s0_branch_for_msolve(
        singular,
        expressions,
        0,
        "s0-boundary",
        timeout,
    )
    replacements = (
        ("t1", "(s1*t0-L)"),
        ("s2", "(s1^2-(13/3)*t0^2-Q)"),
        ("L", "(1)"),
        (fixed_variable, f"({fixed_value})"),
    )
    generators = [
        substitute(polynomial, replacements) for polynomial in polynomials
    ]
    with tempfile.TemporaryDirectory(
        prefix="sic33-trace-hyperslice-exact-"
    ) as directory:
        input_path = Path(directory) / "hyperslice.ms"
        output_path = Path(directory) / "hyperslice.out"
        input_path.write_text(
            f"s5,t4,s3,Q,t2,{free_variable}\n"
            "0\n"
            + ",\n".join(generators)
            + "\n"
        )
        completed = subprocess.run(
            [
                msolve,
                "-f",
                str(input_path),
                "-o",
                str(output_path),
                "-t",
                "4",
                "-v",
                "1",
                "-l",
                "2",
                "-g",
                "1",
            ],
            text=True,
            capture_output=True,
            check=False,
            timeout=timeout,
        )
        if completed.returncode != 0:
            raise AssertionError(
                f"msolve exit {completed.returncode}\n"
                f"{completed.stdout[-2000:]}\n{completed.stderr[-2000:]}"
            )
        result = output_path.read_text()
    return {
        "characteristic": 0,
        "fixed_coordinate": {fixed_variable: fixed_value},
        "free_variables": [free_variable, "Q", "t2", "s3", "s5", "t4"],
        "through_mu8_unit_ideal": result.rstrip().endswith("[1]:"),
        "basis_size": 1 if "#length of basis:      1 element" in result else None,
        "backend": "msolve over Q, exact sparse linear algebra 2",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", default="47,101")
    parser.add_argument("--s1", type=int, default=1)
    parser.add_argument("--t0", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--skip-characteristic-zero", action="store_true")
    arguments = parser.parse_args()
    primes = tuple(int(value) for value in arguments.primes.split(",") if value)
    assert primes and all(3 * max(ORDERS) < prime for prime in primes)
    singular = shutil.which("Singular")
    assert singular is not None
    msolve = shutil.which("msolve")
    assert msolve is not None

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
    exact_result = None
    if not arguments.skip_characteristic_zero:
        exact_result = compute_exact_slice(
            singular,
            arguments.s1,
            arguments.t0,
            arguments.timeout,
        )
        assert exact_result["through_mu8"]["unit_ideal"]
    hyperslices: list[dict[str, object]] = []
    for fixed_variable, fixed_value in (
        ("t0", arguments.t0),
        ("s1", arguments.s1),
    ):
        modular_results = [
            compute_modular_hyperslice(
                singular,
                msolve,
                prime,
                fixed_variable,
                fixed_value % prime,
                arguments.timeout,
            )
            for prime in primes
        ]
        assert all(
            result["through_mu8_unit_ideal"] for result in modular_results
        )
        exact_hyperslice_result = None
        if not arguments.skip_characteristic_zero:
            exact_hyperslice_result = compute_exact_hyperslice(
                singular,
                msolve,
                fixed_variable,
                fixed_value,
                arguments.timeout,
            )
            assert exact_hyperslice_result["through_mu8_unit_ideal"]
        hyperslices.append(
            {
                "fixed_variable": fixed_variable,
                "fixed_value": fixed_value,
                "results": modular_results,
                "characteristic_zero_result": exact_hyperslice_result,
            }
        )

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
        "characteristic_zero_result": exact_result,
        "scope": (
            "complete finite-field slice computation plus an exact "
            "characteristic-zero unit-ideal certificate when present; "
            "not a global nullcone certificate"
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
    hyperslice_outputs: list[Path] = []
    for hyperslice in hyperslices:
        fixed_variable = hyperslice["fixed_variable"]
        fixed_value = hyperslice["fixed_value"]
        hyperslice_output = (
            ROOT
            / "artifacts"
            / "generated-results"
            / (
                "two_pair_sic_bidegree33_boundary_trace_hyperslice_"
                f"{fixed_variable}_{fixed_value}.json"
            )
        )
        hyperslice_payload = {
            "normalization": "s0=1, A=B=mu2=0, L=1",
            "fixed_coordinate": {fixed_variable: fixed_value},
            "corrected_moments_used": list(range(3, 9)),
            "results": hyperslice["results"],
            "characteristic_zero_result": hyperslice[
                "characteristic_zero_result"
            ],
            "scope": (
                "complete finite-field hyperslice unit computations plus an "
                "exact characteristic-zero msolve certificate when present; "
                "not a global boundary certificate"
            ),
            "reproduction_command": payload["reproduction_command"],
        }
        hyperslice_output.write_text(
            json.dumps(hyperslice_payload, indent=2, sort_keys=True) + "\n"
        )
        hyperslice_outputs.append(hyperslice_output)
    for result in results:
        print(
            f"TRACE_SLICE prime={result['prime']} "
            f"length7={result['through_mu7']['quotient_length']} "
            f"basis7={result['through_mu7']['groebner_basis_size']} "
            f"unit8={int(result['through_mu8']['unit_ideal'])}"
        )
    if exact_result is not None:
        print(
            "TRACE_SLICE characteristic=0 "
            f"unit8={int(exact_result['through_mu8']['unit_ideal'])} "
            f"algorithm={exact_result['algorithm']}"
        )
    for hyperslice in hyperslices:
        fixed_variable = hyperslice["fixed_variable"]
        fixed_value = hyperslice["fixed_value"]
        for result in hyperslice["results"]:
            print(
                f"TRACE_HYPERSLICE prime={result['prime']} "
                f"{fixed_variable}={fixed_value} "
                f"unit8={int(result['through_mu8_unit_ideal'])}"
            )
        exact_hyperslice_result = hyperslice["characteristic_zero_result"]
        if exact_hyperslice_result is not None:
            print(
                "TRACE_HYPERSLICE characteristic=0 "
                f"{fixed_variable}={fixed_value} "
                f"unit8={int(exact_hyperslice_result['through_mu8_unit_ideal'])} "
                f"backend={exact_hyperslice_result['backend']}"
            )
    print(f"TRACE_SLICE_WROTE {output.relative_to(ROOT)}")
    for hyperslice_output in hyperslice_outputs:
        print(f"TRACE_HYPERSLICE_WROTE {hyperslice_output.relative_to(ROOT)}")
    print("PASS: corrected mu8 kills the complete slice in every checked field")


if __name__ == "__main__":
    main()
