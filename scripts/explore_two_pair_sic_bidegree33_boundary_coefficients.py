#!/usr/bin/env python3
"""Modular coefficient-ideal test on the reduced SIC(2) (3,3) boundary.

On the s0=1, A=B=0 chart, the exact mu_2, A, and B pivots leave six
base variables and the two fiber variables (s5,t4).  The certified generic
quotient by (mu_4,mu_5) has rank six.  This script reduces the later moments
in that quotient, extracts their fiber coefficients, and computes the
resulting base ideal on the L*Q*J principal open.

The output is modular evidence only.  In particular, a unit ideal must be
reconstructed over characteristic zero before it can be promoted to a
nullcone certificate.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from explore_two_pair_sic_bidegree33_full_anchor import (  # noqa: E402
    chart_expression,
    moment_terms,
    prepare_s0_branch_for_msolve,
)
from verify_two_pair_sic_bidegree33_boundary_generic_quotient import (  # noqa: E402
    T2_NUMERATOR_ADAPTED,
    substitute,
)


BASE_VARIABLES = ("s1", "s2", "s3", "t0", "t1", "t2")
LINEAR = "s1*t0-t1"
QUADRATIC = "s1^2-s2-(13/3)*t0^2"
J_DIVISOR = (
    "9801*s1^4-19602*s1^2*s2-23832*s1^2*t0^2"
    "-60840*s1*t0*t1+9801*s2^2+54252*s2*t0^2"
    "+75076*t0^4+30420*t1^2"
)


def quotient_coefficients(
    singular: str,
    prime: int,
    orders: tuple[int, ...],
    polynomials: list[str],
    timeout: int,
) -> dict[int, list[str]]:
    names = [f"p{order}" for order in orders]
    declarations = "\n".join(
        f"poly {name}={polynomial};"
        for name, polynomial in zip(names, polynomials)
    )
    reductions = []
    for order in orders:
        if order in (3, 4, 5):
            continue
        reductions.append(
            f"""
poly r{order}=reduce(p{order},G);
poly z{order}=r{order};
number c{order};
while (z{order}!=0)
{{
  c{order}=leadcoef(z{order});
  print(
    "COEFFICIENT {order} "
    +string(leadexp(z{order}))
    +" "+string(numerator(c{order}))
  );
  z{order}=z{order}-c{order}*leadmonom(z{order});
}}
"""
        )
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring fiber=({prime},{",".join(BASE_VARIABLES)}),(s5,t4),dp;
{declarations}
ideal G=std(p4,p5);
print("QUOTIENT "+string(vdim(G))+" "+string(size(G)));
int basisIndex;
for (basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print("LEADING "+string(leadexp(G[basisIndex])));
}}
{"".join(reductions)}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    quotient = re.search(r"(?m)^QUOTIENT (\d+) (\d+)$", completed.stdout)
    if quotient is None:
        raise AssertionError(completed.stdout[-4000:])
    print(
        f"QUOTIENT length={quotient.group(1)} basis={quotient.group(2)}",
        flush=True,
    )
    leading = re.findall(r"(?m)^LEADING ([0-9,]+)$", completed.stdout)
    print(f"FIBER_INITIAL {','.join(leading)}", flush=True)
    coefficients: dict[int, list[str]] = {order: [] for order in orders}
    for order, _exponent, coefficient in re.findall(
        r"(?m)^COEFFICIENT (\d+) ([0-9,]+) (.*)$",
        completed.stdout,
    ):
        coefficients[int(order)].append(coefficient)
    for order in orders:
        if order not in (3, 4, 5):
            print(
                f"REMAINDER order={order} "
                f"coefficients={len(coefficients[order])}",
                flush=True,
            )
    return coefficients


def base_ideal(
    singular: str,
    prime: int,
    orders: tuple[int, ...],
    base_polynomial: str,
    coefficients: dict[int, list[str]],
    timeout: int,
) -> None:
    generators = [base_polynomial]
    for order in orders:
        if order in (3, 4, 5):
            continue
        generators.extend(coefficients[order])
    generator_text = ",".join(generators)
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring base={prime},({",".join(BASE_VARIABLES)},u),dp;
poly L={LINEAR};
poly Q={QUADRATIC};
poly J={J_DIVISOR};
ideal I={generator_text},u*L*Q*J-1;
ideal G=std(I);
print(
  "BASE "+string(dim(G))+" "+string(size(G))+" "+string(G[1]==1)
);
int basisIndex;
for (basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print("BASE_LEADING "+string(leadexp(G[basisIndex])));
}}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    marker = re.search(r"(?m)^BASE (-?\d+) (\d+) ([01])$", completed.stdout)
    if marker is None:
        raise AssertionError(completed.stdout[-4000:])
    dimension, size, unit = marker.groups()
    print(
        f"BASE dimension={dimension} basis={size} unit={unit}",
        flush=True,
    )
    leading = re.findall(r"(?m)^BASE_LEADING ([0-9,]+)$", completed.stdout)
    print(f"BASE_INITIAL {','.join(leading)}", flush=True)


def adapted_coefficients(
    singular: str,
    prime: int,
    orders: tuple[int, ...],
    polynomials: list[str],
    timeout: int,
) -> tuple[str, dict[int, list[str]]]:
    """Use the known t2 pivot and return coefficients in adapted variables."""

    replacements = (
        ("t1", "(s1*t0-L)"),
        ("s2", "(s1^2-(13/3)*t0^2-Q)"),
        ("t2", "tt"),
    )
    adapted = [
        substitute(polynomial, replacements)
        for polynomial in polynomials
    ]
    names = [f"p{order}" for order in orders]
    declarations = "\n".join(
        f"poly {name}={polynomial};"
        for name, polynomial in zip(names, adapted)
    )
    reductions = []
    for order in orders:
        if order in (3, 4, 5):
            continue
        reductions.append(
            f"""
poly r{order}=reduce(p{order},G);
poly z{order}=r{order};
number c{order};
while (z{order}!=0)
{{
  c{order}=leadcoef(z{order});
  print(
    "ADAPTED_COEFFICIENT {order} "
    +string(leadexp(z{order}))
    +" "+string(numerator(c{order}))
  );
  z{order}=z{order}-c{order}*leadmonom(z{order});
}}
"""
        )
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring fiber=({prime},s1,s3,t0,L,Q),(s5,t4),dp;
number tt=({T2_NUMERATOR_ADAPTED})/(93366*L*Q);
{declarations}
ideal G=std(p4,p5);
print("ADAPTED_P3 "+string(numerator(leadcoef(p3))));
print("ADAPTED_QUOTIENT "+string(vdim(G))+" "+string(size(G)));
{"".join(reductions)}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    p3_marker = re.search(r"(?m)^ADAPTED_P3 \((.*)\)$", completed.stdout)
    quotient = re.search(
        r"(?m)^ADAPTED_QUOTIENT (\d+) (\d+)$",
        completed.stdout,
    )
    if p3_marker is None or quotient is None:
        raise AssertionError(completed.stdout[-4000:])
    print(
        f"ADAPTED_QUOTIENT length={quotient.group(1)} "
        f"basis={quotient.group(2)}",
        flush=True,
    )
    coefficients: dict[int, list[str]] = {order: [] for order in orders}
    for order, _exponent, coefficient in re.findall(
        r"(?m)^ADAPTED_COEFFICIENT (\d+) ([0-9,]+) (.*)$",
        completed.stdout,
    ):
        coefficients[int(order)].append(coefficient)
    for order in orders:
        if order not in (3, 4, 5):
            print(
                f"ADAPTED_REMAINDER order={order} "
                f"coefficients={len(coefficients[order])}",
                flush=True,
            )
    metadata_declarations = [f"poly c3={p3_marker.group(1)};"]
    metadata_prints = [
        'print("ADAPTED_COEFFICIENT_META 3 0 "+string(deg(c3))+" "+string(size(c3)));'
    ]
    for order in orders:
        for index, coefficient in enumerate(coefficients[order]):
            name = f"c{order}_{index}"
            metadata_declarations.append(f"poly {name}={coefficient};")
            metadata_prints.append(
                f'print("ADAPTED_COEFFICIENT_META {order} {index} "'
                f'+string(deg({name}))+" "+string(size({name})));'
            )
    metadata = subprocess.run(
        [singular, "-q"],
        input=f"""
ring metadata={prime},(s1,s3,t0,L,Q),dp;
{"".join(metadata_declarations)}
{"".join(metadata_prints)}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    for order, index, degree, terms in re.findall(
        r"(?m)^ADAPTED_COEFFICIENT_META (\d+) (\d+) (\d+) (\d+)$",
        metadata.stdout,
    ):
        print(
            f"ADAPTED_COEFFICIENT_META order={order} index={index} "
            f"degree={degree} terms={terms}",
            flush=True,
        )
    return p3_marker.group(1), coefficients


def adapted_base_ideal(
    singular: str,
    prime: int,
    orders: tuple[int, ...],
    base_polynomial: str,
    coefficients: dict[int, list[str]],
    timeout: int,
) -> None:
    generators = [base_polynomial]
    for order in orders:
        if order in (3, 4, 5):
            continue
        generators.extend(coefficients[order])
    adapted_j = "(99*Q+155*t0^2)^2+30420*L^2"
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring base={prime},(s1,s3,t0,L,Q,u),dp;
ideal I={",".join(generators)},u*L*Q*({adapted_j})-1;
ideal G=std(I);
print(
  "ADAPTED_BASE "+string(dim(G))+" "+string(size(G))+" "
  +string(G[1]==1)
);
int basisIndex;
for (basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print("ADAPTED_BASE_LEADING "+string(leadexp(G[basisIndex])));
}}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    marker = re.search(
        r"(?m)^ADAPTED_BASE (-?\d+) (\d+) ([01])$",
        completed.stdout,
    )
    if marker is None:
        raise AssertionError(completed.stdout[-4000:])
    dimension, size, unit = marker.groups()
    print(
        f"ADAPTED_BASE dimension={dimension} basis={size} unit={unit}",
        flush=True,
    )
    leading = re.findall(
        r"(?m)^ADAPTED_BASE_LEADING ([0-9,]+)$",
        completed.stdout,
    )
    print(f"ADAPTED_BASE_INITIAL {','.join(leading)}", flush=True)


def adapted_base_msolve(
    msolve: str,
    prime: int,
    orders: tuple[int, ...],
    base_polynomial: str,
    coefficients: dict[int, list[str]],
    timeout: int,
    threads: int,
    linear_algebra: int,
    coefficient_counts: dict[int, int],
    normalize_linear: bool,
) -> None:
    generators = [base_polynomial]
    for order in orders:
        if order in (3, 4, 5):
            continue
        generators.extend(
            coefficients[order][: coefficient_counts.get(order, len(coefficients[order]))]
        )
    saturation = "L*Q*((99*Q+155*t0^2)^2+30420*L^2)"
    variables = "s1,s3,t0,L,Q"
    if normalize_linear:
        replacements = (("L", "(1)"),)
        generators = [
            substitute(generator, replacements)
            for generator in generators
        ]
        saturation = substitute(saturation, replacements)
        variables = "s1,s3,t0,Q"
    generators.append(saturation)
    with tempfile.TemporaryDirectory(prefix="sic33-boundary-base-") as directory:
        input_path = Path(directory) / "base.ms"
        output_path = Path(directory) / "base.out"
        input_path.write_text(
            variables
            + "\n"
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
                str(threads),
                "-v",
                "1",
                "-l",
                str(linear_algebra),
                "-S",
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
                f"msolve exit {completed.returncode}\n{completed.stdout[-2000:]}\n"
                f"{completed.stderr[-2000:]}"
            )
        result = output_path.read_text().strip()
    if result in ("[-1]", "[-1]:"):
        status = "unit"
    elif result.startswith("[1,") and ",-1,[]" in result.replace(" ", ""):
        status = "positive-dimensional"
    else:
        status = "finite-nonempty"
    leading_lines = [
        line
        for line in completed.stdout.splitlines()
        if "leading" in line.lower() or line.startswith("[")
    ]
    print(f"ADAPTED_MSOLVE status={status}", flush=True)
    if status != "unit":
        print(
            "ADAPTED_MSOLVE_RESULT "
            + result.replace("\n", " ")[:1000],
            flush=True,
        )
    if leading_lines:
        print(
            "ADAPTED_MSOLVE_INITIAL "
            + " | ".join(leading_lines[-10:]),
            flush=True,
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=47)
    parser.add_argument("--orders", default="2,3,4,5,6")
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--adapted", action="store_true")
    parser.add_argument(
        "--backend",
        choices=("singular", "msolve"),
        default="singular",
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--linear-algebra",
        choices=(1, 2, 42, 44),
        type=int,
        default=2,
    )
    parser.add_argument(
        "--coefficient-counts",
        default="",
        help=(
            "optional comma-separated order:count limits for the adapted "
            "msolve base ideal, for example 6:3,7:1"
        ),
    )
    parser.add_argument(
        "--normalize-l",
        action="store_true",
        help="use weighted homogeneity and L!=0 to impose L=1",
    )
    arguments = parser.parse_args()
    orders = tuple(int(value) for value in arguments.orders.split(",") if value)
    coefficient_counts = {
        int(order): int(count)
        for item in arguments.coefficient_counts.split(",")
        if item
        for order, count in (item.split(":", 1),)
    }
    assert orders[:4] == (2, 3, 4, 5)
    assert 3 * max(orders) < arguments.prime
    singular = shutil.which("Singular")
    assert singular is not None
    msolve = shutil.which("msolve")
    if arguments.backend == "msolve":
        assert msolve is not None

    expressions = [
        chart_expression(moment_terms(order, arguments.prime), 0, arguments.prime)
        for order in orders
    ]
    variables, polynomials = prepare_s0_branch_for_msolve(
        singular,
        expressions,
        arguments.prime,
        "s0-boundary",
        arguments.timeout,
    )
    assert variables == ("s1", "s2", "s3", "s5", "t0", "t1", "t2", "t4")
    reduced_orders = orders[1:]
    if arguments.adapted:
        base_polynomial, coefficients = adapted_coefficients(
            singular,
            arguments.prime,
            reduced_orders,
            polynomials,
            arguments.timeout,
        )
        if arguments.backend == "msolve":
            adapted_base_msolve(
                msolve,
                arguments.prime,
                reduced_orders,
                base_polynomial,
                coefficients,
                arguments.timeout,
                arguments.threads,
                arguments.linear_algebra,
                coefficient_counts,
                arguments.normalize_l,
            )
        else:
            adapted_base_ideal(
                singular,
                arguments.prime,
                reduced_orders,
                base_polynomial,
                coefficients,
                arguments.timeout,
            )
    else:
        coefficients = quotient_coefficients(
            singular,
            arguments.prime,
            reduced_orders,
            polynomials,
            arguments.timeout,
        )
        base_ideal(
            singular,
            arguments.prime,
            reduced_orders,
            polynomials[0],
            coefficients,
            arguments.timeout,
        )
    print(
        "EVIDENCE ONLY: modular L*Q*J-open coefficient ideal; "
        "no characteristic-zero theorem is promoted"
    )


if __name__ == "__main__":
    main()
