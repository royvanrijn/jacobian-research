#!/usr/bin/env python3
"""Modular full-anchor test for two-pair SIC in bidegree (3,3).

Normalize the non-null Sym^2 component to 2*X*T.  Its stabilizer is a
one-dimensional torus.  The locus with only torus-weight-zero higher
components is already closed by the diagonal-slice theorem.  The remaining
locus is covered, up to Weyl reflection, by five charts obtained by setting
one of s0,s1,s2,t0,t1 to one.

This script generates the restricted moments by their torus weights and
tests the five full Sym^6+Sym^4+Sym^2 charts over a finite field.  Its output
is evidence only unless a separate characteristic-zero certificate is made.
On the s0 chart it can also perform the exact mu_2 substitution and select
one of the two principal opens or their common mu_3-pivot boundary before
calling Singular or exporting the reduced variables to msolve.  On A=0 it
uses the constant t3 pivot in A; on A=B=0 it additionally uses the constant
s4 pivot in B.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from functools import lru_cache
from math import factorial
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time


ROOT = Path(__file__).resolve().parents[1]
PARAMETERS = (
    "s0", "s1", "s2", "s3", "s4", "s5", "s6",
    "t0", "t1", "t2", "t3", "t4",
)
WEIGHTS = (3, 2, 1, 0, -1, -2, -3, 2, 1, 0, -1, -2)

# After removing the unmatched X or Y power prescribed by WEIGHTS, each
# irreducible-basis coefficient is a polynomial in q=X*Y.
Q_POLYNOMIALS = (
    (1,),
    (-3, 3),
    (3, -9, 3),
    (-1, 9, -9, 1),
    (-3, 9, -3),
    (-3, 3),
    (-1,),
    (1, 1),
    (-2, 0, 2),
    (1, -3, -3, 1),
    (2, 0, -2),
    (1, 1),
)

# The normalized quadratic 2*X*T has biform coefficient polynomial
# -1-q+q^2+q^3 in the chosen divided-power basis.
QUADRATIC_Q = (-1, -1, 1, 1)
REPRESENTATIVE_CHARTS = (0, 1, 2, 7, 8)
OPPOSITE_PIVOTS = {
    0: 6,
    1: 5,
    2: 4,
    7: 11,
    8: 10,
}
S0_MU3_LINEAR_COEFFICIENTS = (
    (
        "6*s1^2*t1-3*s1*s2*t0-3*s1*t2-3*s2*t1"
        "+2*s3*t0-3*t0+t3"
    ),
    (
        "12*s1*s3+28*s1*t0*t1-18*s1-9*s2^2-14*s2*t0^2"
        "-3*s4-2*t0*t2-12*t1^2"
    ),
)


def convolve(left: tuple[int, ...], right: tuple[int, ...], prime: int) -> tuple[int, ...]:
    answer = [0] * (len(left) + len(right) - 1)
    for left_index, left_coefficient in enumerate(left):
        if left_coefficient % prime == 0:
            continue
        for right_index, right_coefficient in enumerate(right):
            answer[left_index + right_index] = (
                answer[left_index + right_index]
                + left_coefficient * right_coefficient
            ) % prime
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def polynomial_powers(
    polynomial: tuple[int, ...], maximum: int, prime: int
) -> tuple[tuple[int, ...], ...]:
    powers = [(1,)]
    for _ in range(maximum):
        powers.append(convolve(powers[-1], polynomial, prime))
    return tuple(powers)


def moment_terms(order: int, prime: int) -> dict[tuple[int, ...], int]:
    """Return mu_order after setting the quadratic scale c to one."""

    factorials = [factorial(index) % prime for index in range(3 * order + 1)]
    inverse_factorials = [
        pow(factorial(index) % prime, -1, prime)
        for index in range(order + 1)
    ]
    basis_powers = tuple(
        polynomial_powers(polynomial, order, prime)
        for polynomial in Q_POLYNOMIALS
    )
    quadratic_powers = polynomial_powers(QUADRATIC_Q, order, prime)
    parameter_order = (0, 6, 1, 5, 7, 11, 2, 4, 8, 10, 3, 9)
    exponents = [0] * len(PARAMETERS)
    answer: dict[tuple[int, ...], int] = defaultdict(int)
    order_factorial = factorials[order]

    @lru_cache(maxsize=None)
    def remaining_weight_bounds(position: int, degree_left: int) -> tuple[int, int]:
        remaining_weights = [WEIGHTS[index] for index in parameter_order[position:]]
        if not remaining_weights or degree_left == 0:
            return 0, 0
        return (
            min(0, degree_left * min(remaining_weights)),
            max(0, degree_left * max(remaining_weights)),
        )

    def visit(
        position: int,
        used_degree: int,
        weight: int,
        shift: int,
        inverse_denominator: int,
        q_polynomial: tuple[int, ...],
    ) -> None:
        if position == len(parameter_order):
            if weight != 0:
                return
            quadratic_exponent = order - used_degree
            product = convolve(
                q_polynomial,
                quadratic_powers[quadratic_exponent],
                prime,
            )
            scalar = (
                order_factorial
                * inverse_denominator
                * inverse_factorials[quadratic_exponent]
            ) % prime
            contraction = 0
            for q_degree, coefficient in enumerate(product):
                diagonal = shift + q_degree
                if 0 <= diagonal <= 3 * order:
                    contraction += (
                        coefficient
                        * factorials[3 * order - diagonal]
                        * factorials[diagonal]
                    )
            coefficient = scalar * contraction % prime
            if coefficient:
                exponent_tuple = tuple(exponents)
                answer[exponent_tuple] = (
                    answer[exponent_tuple] + coefficient
                ) % prime
            return

        parameter_index = parameter_order[position]
        parameter_weight = WEIGHTS[parameter_index]
        available = order - used_degree
        for exponent in range(available + 1):
            new_weight = weight + exponent * parameter_weight
            degree_left = available - exponent
            minimum, maximum = remaining_weight_bounds(position + 1, degree_left)
            if not minimum <= -new_weight <= maximum:
                continue
            exponents[parameter_index] = exponent
            visit(
                position + 1,
                used_degree + exponent,
                new_weight,
                shift + max(parameter_weight, 0) * exponent,
                inverse_denominator * inverse_factorials[exponent] % prime,
                convolve(
                    q_polynomial,
                    basis_powers[parameter_index][exponent],
                    prime,
                ),
            )
        exponents[parameter_index] = 0

    visit(0, 0, 0, 0, 1, (1,))
    return {
        exponents: coefficient
        for exponents, coefficient in answer.items()
        if coefficient % prime
    }


def chart_expression(
    terms: dict[tuple[int, ...], int],
    fixed_index: int,
    prime: int,
) -> str:
    combined: dict[tuple[int, ...], int] = defaultdict(int)
    for exponents, coefficient in terms.items():
        reduced = exponents[:fixed_index] + exponents[fixed_index + 1 :]
        combined[reduced] = (combined[reduced] + coefficient) % prime

    serialized: list[str] = []
    variable_names = PARAMETERS[:fixed_index] + PARAMETERS[fixed_index + 1 :]
    for exponents, coefficient in sorted(combined.items()):
        coefficient %= prime
        if not coefficient:
            continue
        factors: list[str] = []
        for variable, exponent in zip(variable_names, exponents):
            if exponent == 1:
                factors.append(variable)
            elif exponent > 1:
                factors.append(f"{variable}^{exponent}")
        monomial = "*".join(factors)
        if not monomial:
            serialized.append(str(coefficient))
        elif coefficient == 1:
            serialized.append(monomial)
        else:
            serialized.append(f"{coefficient}*{monomial}")
    return "+".join(serialized) or "0"


def s0_branch_setup(expressions: list[str], branch: str) -> str:
    """Return Singular code for one exact reduced s0 branch."""

    assert branch != "none"
    coefficient_a, coefficient_b = S0_MU3_LINEAR_COEFFICIENTS
    common_setup = f"""
ideal sourceI={",".join(expressions)};
poly pivotCoefficient=diff(sourceI[1],s6);
poly pivotValue=-(sourceI[1]-pivotCoefficient*s6)/pivotCoefficient;
ideal reducedI;
int equation;
for (equation=2;equation<=size(sourceI);equation++)
{{
  reducedI[size(reducedI)+1]=subst(sourceI[equation],s6,pivotValue);
}}
poly A={coefficient_a};
poly B={coefficient_b};
"""
    if branch == "s0-A-open":
        branch_setup = """
poly thirdRest=reducedI[1]+103680*A*s5;
poly s5Value=thirdRest*ainv/103680;
ideal I=ainv*A-1;
for (equation=2;equation<=size(reducedI);equation++)
{
  I[size(I)+1]=subst(reducedI[equation],s5,s5Value);
}
"""
    elif branch == "s0-A-open-sparse":
        branch_setup = """
ideal I=ainv*A-1;
for (equation=1;equation<=size(reducedI);equation++)
{
  I[size(I)+1]=reducedI[equation];
}
"""
    elif branch == "s0-B-open":
        branch_setup = """
poly t3Value=-(A-t3);
poly thirdOnBoundary=subst(reducedI[1],t3,t3Value);
poly thirdRest=thirdOnBoundary+17280*B*t4;
poly t4Value=thirdRest*binv/17280;
ideal I=binv*B-1;
for (equation=2;equation<=size(reducedI);equation++)
{
  I[size(I)+1]=subst(
    subst(reducedI[equation],t3,t3Value),t4,t4Value
  );
}
"""
    elif branch == "s0-B-open-sparse":
        branch_setup = """
poly t3Value=-(A-t3);
ideal I=binv*B-1;
for (equation=1;equation<=size(reducedI);equation++)
{
  I[size(I)+1]=subst(reducedI[equation],t3,t3Value);
}
"""
    else:
        assert branch == "s0-boundary"
        branch_setup = """
poly t3Value=-(A-t3);
poly s4Value=(B+3*s4)/3;
ideal I;
for (equation=1;equation<=size(reducedI);equation++)
{
  I[size(I)+1]=subst(
    subst(reducedI[equation],t3,t3Value),s4,s4Value
  );
}
"""
    return common_setup + branch_setup


def run_chart(
    singular: str,
    fixed_index: int,
    expressions: list[str],
    prime: int,
    timeout: int,
    algorithm: str,
    ordering: str,
    branch: str,
) -> tuple[int, int, bool, float]:
    variables = list(
        PARAMETERS[:fixed_index] + PARAMETERS[fixed_index + 1 :]
    )
    pivot = PARAMETERS[OPPOSITE_PIVOTS[fixed_index]]
    if branch in ("s0-A-open", "s0-A-open-sparse"):
        variables.append("ainv")
    elif branch in ("s0-B-open", "s0-B-open-sparse"):
        variables.append("binv")
    if ordering == "pivot-block":
        variables.remove(pivot)
        variables.insert(0, pivot)
        ring_ordering = f"(lp(1),dp({len(variables) - 1}))"
    else:
        ring_ordering = "dp"
    if branch != "none":
        assert fixed_index == 0
        ideal_setup = s0_branch_setup(expressions, branch)
    elif ordering == "pivot-substitute":
        ideal_setup = f"""
ideal sourceI={",".join(expressions)};
poly pivotCoefficient=diff(sourceI[1],{pivot});
if (deg(pivotCoefficient)>0)
{{
  print("PIVOT_ERROR");
  exit;
}}
poly pivotValue=-(sourceI[1]-pivotCoefficient*{pivot})/pivotCoefficient;
ideal I;
int equation;
for (equation=2;equation<=size(sourceI);equation++)
{{
  I[size(I)+1]=subst(sourceI[equation],{pivot},pivotValue);
}}
"""
    else:
        ideal_setup = f"ideal I={','.join(expressions)};"
    if algorithm == "incremental":
        basis_setup = """
ideal G=std(I[1]);
int basisEquation;
for (basisEquation=2;basisEquation<=size(I);basisEquation++)
{
  G=std(G+I[basisEquation]);
  if (G[1]==1) { break; }
}
"""
    else:
        basis_setup = f"ideal G={algorithm}(I);"
    started = time.monotonic()
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring anchor={prime},({",".join(variables)}),{ring_ordering};
option(redSB);
{ideal_setup}
{basis_setup}
print("ANCHOR "+string(dim(G))+" "+string(size(G))+" "+string(G[1]==1));
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    elapsed = time.monotonic() - started
    marker = re.search(r"(?m)^ANCHOR (-?\d+) (\d+) ([01])$", completed.stdout)
    if marker is None:
        raise AssertionError(completed.stdout[-2000:])
    dimension, basis_size, unit = marker.groups()
    return int(dimension), int(basis_size), unit == "1", elapsed


def prepare_s0_branch_for_msolve(
    singular: str,
    expressions: list[str],
    prime: int,
    branch: str,
    timeout: int,
) -> tuple[tuple[str, ...], list[str]]:
    """Use Singular only to perform the two exact branch substitutions."""

    ring_variables = list(PARAMETERS[1:])
    if branch == "s0-A-open":
        ring_variables.append("ainv")
        output_variables = tuple(
            variable
            for variable in ring_variables
            if variable not in ("s5", "s6")
        )
    elif branch == "s0-A-open-sparse":
        ring_variables.append("ainv")
        output_variables = tuple(
            variable for variable in ring_variables if variable != "s6"
        )
    elif branch == "s0-B-open":
        ring_variables.append("binv")
        output_variables = tuple(
            variable
            for variable in ring_variables
            if variable not in ("s6", "t3", "t4")
        )
    elif branch == "s0-B-open-sparse":
        ring_variables.append("binv")
        output_variables = tuple(
            variable
            for variable in ring_variables
            if variable not in ("s6", "t3")
        )
    else:
        assert branch == "s0-boundary"
        output_variables = tuple(
            variable
            for variable in ring_variables
            if variable not in ("s4", "s6", "t3")
        )
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring anchor={prime},({",".join(ring_variables)}),dp;
{s0_branch_setup(expressions, branch)}
print("BRANCH_BEGIN "+string(size(I)));
int outputEquation;
for (outputEquation=1;outputEquation<=size(I);outputEquation++)
{{
  print("BRANCH_POLY "+string(I[outputEquation]));
}}
print("BRANCH_END");
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    size_marker = re.search(
        r"(?m)^BRANCH_BEGIN (\d+)$",
        completed.stdout,
    )
    assert size_marker is not None, completed.stdout[-2000:]
    polynomials = re.findall(
        r"(?m)^BRANCH_POLY (.*)$",
        completed.stdout,
    )
    assert len(polynomials) == int(size_marker.group(1))
    assert "BRANCH_END" in completed.stdout
    return output_variables, polynomials


def run_chart_msolve(
    msolve: str,
    fixed_index: int,
    expressions: list[str],
    prime: int,
    timeout: int,
    threads: int,
    variables: tuple[str, ...] | None = None,
    linear_algebra: int = 2,
) -> tuple[str, float]:
    """Classify a modular chart with msolve."""

    if variables is None:
        variables = PARAMETERS[:fixed_index] + PARAMETERS[fixed_index + 1 :]
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="sic33-msolve-") as directory:
        input_path = Path(directory) / "system.ms"
        output_path = Path(directory) / "result.ms"
        polynomial_lines = [
            expression + ("," if index + 1 < len(expressions) else "")
            for index, expression in enumerate(expressions)
        ]
        input_path.write_text(
            ",".join(variables)
            + "\n"
            + str(prime)
            + "\n"
            + "\n".join(polynomial_lines)
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
            ],
            text=True,
            capture_output=True,
            check=False,
            timeout=timeout,
        )
        if completed.returncode != 0:
            return f"solver-error-{completed.returncode}", (
                time.monotonic() - started
            )
        result = output_path.read_text().strip()
    elapsed = time.monotonic() - started
    if result == "[-1]:" or result == "[-1]":
        return "unit", elapsed
    if result.startswith("[1,") and ",-1,[]" in result.replace(" ", ""):
        return "positive-dimensional", elapsed
    return "finite-nonempty", elapsed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=101)
    parser.add_argument("--max-order", type=int, default=10)
    parser.add_argument(
        "--orders",
        default="",
        help=(
            "comma-separated moment orders; when omitted, use every order "
            "from 2 through --max-order"
        ),
    )
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument(
        "--algorithm",
        choices=("std", "slimgb", "incremental"),
        default="slimgb",
    )
    parser.add_argument(
        "--ordering",
        choices=("dp", "pivot-block", "pivot-substitute"),
        default="pivot-substitute",
        help=(
            "put the coefficient opposite the normalized chart variable "
            "in a one-variable elimination block"
        ),
    )
    parser.add_argument(
        "--backend",
        choices=("singular", "msolve"),
        default="singular",
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--msolve-linear-algebra",
        choices=(1, 2, 42, 44),
        type=int,
        default=2,
    )
    parser.add_argument(
        "--extra-equations",
        choices=("none", "s0-mu3-boundary"),
        default="none",
        help=(
            "optionally impose the exact boundary where both coefficients "
            "used to solve the reduced s0-chart third moment vanish"
        ),
    )
    parser.add_argument(
        "--branch",
        choices=(
            "none",
            "s0-A-open",
            "s0-A-open-sparse",
            "s0-B-open",
            "s0-B-open-sparse",
            "s0-boundary",
        ),
        default="none",
        help=(
            "perform the exact mu_2 elimination and one branch of the "
            "reduced s0-chart mu_3 split before Gröbner computation"
        ),
    )
    parser.add_argument(
        "--charts",
        default="s0,s1,s2,t0,t1",
        help="comma-separated representative charts",
    )
    arguments = parser.parse_args()
    moment_orders = (
        tuple(int(order) for order in arguments.orders.split(",") if order)
        if arguments.orders
        else tuple(range(2, arguments.max_order + 1))
    )
    assert moment_orders
    assert len(set(moment_orders)) == len(moment_orders)
    assert min(moment_orders) >= 2
    assert 3 * max(moment_orders) < arguments.prime
    if arguments.branch != "none":
        assert moment_orders[:2] == (2, 3)
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    msolve = shutil.which("msolve")
    if arguments.backend == "msolve":
        assert msolve is not None, "msolve is required"

    requested = tuple(
        PARAMETERS.index(name)
        for name in arguments.charts.split(",")
        if name
    )
    assert set(requested) <= set(REPRESENTATIVE_CHARTS)

    all_terms: dict[int, dict[tuple[int, ...], int]] = {}
    for order in moment_orders:
        started = time.monotonic()
        terms = moment_terms(order, arguments.prime)
        all_terms[order] = terms
        print(
            f"MOMENT {order} terms={len(terms)} "
            f"seconds={time.monotonic()-started:.2f}",
            flush=True,
        )

    for fixed_index in requested:
        expressions = [
            chart_expression(all_terms[order], fixed_index, arguments.prime)
            for order in moment_orders
        ]
        if arguments.extra_equations == "s0-mu3-boundary":
            assert fixed_index == 0
            expressions.extend(S0_MU3_LINEAR_COEFFICIENTS)
        if arguments.backend == "msolve":
            try:
                msolve_variables = None
                if arguments.branch != "none":
                    msolve_variables, expressions = (
                        prepare_s0_branch_for_msolve(
                            singular,
                            expressions,
                            arguments.prime,
                            arguments.branch,
                            arguments.timeout,
                        )
                    )
                    print(
                        f"BRANCH {arguments.branch} "
                        f"variables={len(msolve_variables)} "
                        f"equations={len(expressions)}",
                        flush=True,
                    )
                status, elapsed = run_chart_msolve(
                    msolve,
                    fixed_index,
                    expressions,
                    arguments.prime,
                    arguments.timeout,
                    arguments.threads,
                    msolve_variables,
                    arguments.msolve_linear_algebra,
                )
            except subprocess.TimeoutExpired:
                print(
                    f"CHART {PARAMETERS[fixed_index]} TIMEOUT "
                    f"seconds={arguments.timeout}",
                    flush=True,
                )
                continue
            print(
                f"CHART {PARAMETERS[fixed_index]} status={status} "
                f"seconds={elapsed:.2f}",
                flush=True,
            )
            continue
        try:
            dimension, basis_size, unit, elapsed = run_chart(
                singular,
                fixed_index,
                expressions,
                arguments.prime,
                arguments.timeout,
                arguments.algorithm,
                arguments.ordering,
                arguments.branch,
            )
        except subprocess.TimeoutExpired:
            print(
                f"CHART {PARAMETERS[fixed_index]} TIMEOUT "
                f"seconds={arguments.timeout}",
                flush=True,
            )
            continue
        print(
            f"CHART {PARAMETERS[fixed_index]} dimension={dimension} "
            f"basis={basis_size} unit={int(unit)} seconds={elapsed:.2f}",
            flush=True,
        )

    print(
        "EVIDENCE ONLY: modular full-anchor charts; "
        "no characteristic-zero theorem is promoted"
    )


if __name__ == "__main__":
    main()
