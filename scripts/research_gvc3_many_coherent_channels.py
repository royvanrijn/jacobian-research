#!/usr/bin/env python3
"""Research four/five coherent harmonic channels below Delta^6.

For channel degrees ``ell_i`` and isotropic directions represented by points
of P1, this script compiles scalar Reynolds moments from the Gaussian/Wick
coefficient recurrence

    s_i W(s) = sum_{j != i} g_ij W(s-e_i-e_j),

where ``g_ij`` is the squared bracket of the two P1 points.  Three direction
points are normalized to infinity, zero, and one; later points are cross-ratio
parameters.  Modular Singular elimination is discovery.  An optional exact
msolve unit basis is accepted only after the same cutoff is found modularly.

The default case is the four-channel degree-eight tuple (2,4,6,8) on the
pairwise-distinct direction and nonzero-coefficient chart.  Direction
collision partitions can be supplied with ``--groups``.
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import itertools
import json
import math
import subprocess
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "gvc3_many_coherent_channels_research.json"
)
PRIMES = (101, 103, 107)


def compositions(total: int, length: int, prefix: tuple[int, ...] = ()):
    if length == 1:
        yield prefix + (total,)
        return
    for entry in range(total + 1):
        yield from compositions(total - entry, length - 1, prefix + (entry,))


def parse_degrees(text: str) -> tuple[int, ...]:
    degrees = tuple(int(entry) for entry in text.split(","))
    if len(degrees) not in (4, 5):
        raise argparse.ArgumentTypeError("degrees must contain four or five entries")
    if any(degree <= 0 or degree % 2 for degree in degrees):
        raise argparse.ArgumentTypeError("harmonic degrees must be positive and even")
    return degrees


def parse_groups(text: str, channels: int) -> tuple[tuple[int, ...], ...]:
    try:
        blocks = tuple(
            tuple(int(entry) for entry in block.split(","))
            for block in text.split("|")
        )
    except ValueError as error:
        raise argparse.ArgumentTypeError("groups must look like 0,1|2|3") from error
    flattened = tuple(itertools.chain.from_iterable(blocks))
    if sorted(flattened) != list(range(channels)) or any(not block for block in blocks):
        raise argparse.ArgumentTypeError("groups must partition every channel index once")
    return tuple(sorted((tuple(sorted(block)) for block in blocks), key=lambda block: block[0]))


def normalized_points(directions: int, parameters: tuple[sp.Symbol, ...]):
    if directions < 1 or directions > 5:
        raise ValueError("this compiler supports one through five directions")
    points = ((1, 0), (0, 1), (1, 1))[:directions]
    if directions >= 4:
        points += ((parameters[0], 1),)
    if directions >= 5:
        points += ((parameters[1], 1),)
    return points


def gram_matrix(directions: int, parameters: tuple[sp.Symbol, ...]):
    points = normalized_points(directions, parameters)
    matrix = [[sp.Integer(0) for _ in range(directions)] for _ in range(directions)]
    for left in range(directions):
        for right in range(left + 1, directions):
            u_left, v_left = points[left]
            u_right, v_right = points[right]
            bracket = u_left * v_right - v_left * u_right
            matrix[left][right] = matrix[right][left] = sp.expand(bracket**2)
    return tuple(tuple(row) for row in matrix)


def configuration_discriminant(directions: int, parameters: tuple[sp.Symbol, ...]):
    if directions <= 3:
        return sp.Integer(1)
    lam = parameters[0]
    answer = lam * (lam - 1)
    if directions == 5:
        mu = parameters[1]
        answer *= mu * (mu - 1) * (lam - mu)
    return sp.expand(answer)


def wick_compiler(grams: tuple[tuple[sp.Expr, ...], ...]):
    directions = len(grams)

    @functools.lru_cache(maxsize=None)
    def wick(stubs: tuple[int, ...]) -> sp.Expr:
        if not any(stubs):
            return sp.Integer(1)
        if sum(stubs) % 2 or 2 * max(stubs) > sum(stubs):
            return sp.Integer(0)
        pivot = next(index for index, degree in enumerate(stubs) if degree)
        reduced_pivot = list(stubs)
        reduced_pivot[pivot] -= 1
        answer = sp.Integer(0)
        for other in range(directions):
            if other == pivot or reduced_pivot[other] == 0:
                continue
            reduced = list(reduced_pivot)
            reduced[other] -= 1
            answer += grams[pivot][other] * wick(tuple(reduced))
        return sp.cancel(answer / stubs[pivot])

    return wick


def odd_double_factorial(value: int) -> int:
    answer = 1
    for factor in range(1, value + 1, 2):
        answer *= factor
    return answer


def moment(
    degrees: tuple[int, ...],
    groups: tuple[tuple[int, ...], ...],
    order: int,
    coefficients: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
) -> sp.Expr:
    channel_group = {}
    for group_index, block in enumerate(groups):
        for channel in block:
            channel_group[channel] = group_index
    grams = gram_matrix(len(groups), parameters)
    wick = wick_compiler(grams)
    answer = sp.Integer(0)
    for counts in compositions(order, len(degrees)):
        stubs = [0] * len(groups)
        for channel, count in enumerate(counts):
            stubs[channel_group[channel]] += degrees[channel] * count
        contraction = wick(tuple(stubs))
        if contraction == 0:
            continue
        multinomial = math.factorial(order)
        for count in counts:
            multinomial //= math.factorial(count)
        numerator = multinomial
        for stub_count in stubs:
            numerator *= math.factorial(stub_count)
        denominator = odd_double_factorial(sum(stubs) + 1)
        monomial = sp.prod(
            coefficient**count
            for coefficient, count in zip(coefficients, counts, strict=True)
        )
        answer += sp.Rational(numerator, denominator) * contraction * monomial
    return sp.expand(answer)


def primitive_polynomial(expression: sp.Expr, variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
    return polynomial.clear_denoms()[1].primitive()[1].as_expr()


def linear_pivot_localization(
    equations: dict[int, sp.Expr],
    pivot: sp.Symbol,
    variables: tuple[sp.Symbol, ...],
    saturation: sp.Expr,
) -> tuple[dict[int, sp.Expr], tuple[sp.Symbol, ...], sp.Expr, dict[str, object]]:
    """Solve the first equation when it is linear in ``pivot``.

    If ``A+B*pivot=0``, work on ``A*B != 0`` so that both the pivot and its
    coefficient remain nonzero.  Later equations are cleared by the smallest
    required power of B.
    """
    first_order = min(equations)
    first = sp.Poly(equations[first_order], pivot)
    if first.degree() != 1:
        raise ValueError(f"moment {first_order} is not linear in {pivot}")
    coefficient = sp.expand(first.coeff_monomial(pivot))
    remainder = sp.expand(first.coeff_monomial(1))
    base_variables = tuple(variable for variable in variables if variable != pivot)
    transformed: dict[int, sp.Expr] = {}
    transformed_counts: dict[str, int] = {}
    for order, equation in equations.items():
        if order == first_order:
            continue
        polynomial = sp.Poly(equation, pivot)
        degree = polynomial.degree()
        cleared = sp.Integer(0)
        for power in range(degree + 1):
            part = polynomial.coeff_monomial(pivot**power)
            if part:
                cleared += part * (-remainder) ** power * coefficient ** (degree - power)
        primitive = primitive_polynomial(sp.expand(cleared), base_variables)
        if primitive != 0:
            transformed[order] = primitive
            transformed_counts[str(order)] = len(
                sp.Poly(primitive, *base_variables).terms()
            )
    reduced_saturation = sp.expand(
        saturation.subs(pivot, 1) * remainder * coefficient
    )
    metadata = {
        "first_order": first_order,
        "pivot": str(pivot),
        "equation": "A+B*pivot=0",
        "A_sha256": hashlib.sha256(str(remainder).encode()).hexdigest(),
        "B_sha256": hashlib.sha256(str(coefficient).encode()).hexdigest(),
        "A_terms": len(sp.Poly(remainder, *base_variables).terms()),
        "B_terms": len(sp.Poly(coefficient, *base_variables).terms()),
        "transformed_term_counts": transformed_counts,
        "chart": "A*B!=0 (equivalently pivot!=0 and B!=0)",
    }
    return transformed, base_variables, reduced_saturation, metadata


def singular_expression(expression: sp.Expr, variables: tuple[sp.Symbol, ...], prime: int):
    polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
    terms: list[str] = []
    for exponent, coefficient in polynomial.terms():
        modular = (
            int(coefficient.p) * pow(int(coefficient.q), -1, prime)
        ) % prime
        if modular == 0:
            continue
        factors = []
        for variable, power in zip(variables, exponent, strict=True):
            if power == 1:
                factors.append(str(variable))
            elif power > 1:
                factors.append(f"{variable}^{power}")
        monomial = "*".join(factors)
        if not monomial:
            terms.append(str(modular))
        elif modular == 1:
            terms.append(monomial)
        else:
            terms.append(f"{modular}*{monomial}")
    return "+".join(terms) or "0"


def modular_cutoff(
    equations: dict[int, sp.Expr],
    saturation: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    prime: int,
    singular: str,
    timeout: int,
) -> dict[str, object]:
    z = sp.Symbol("z")
    ring_variables = (z,) + variables
    saturation_expression = singular_expression(z * saturation - 1, ring_variables, prime)
    expressions = [saturation_expression]
    expressions.extend(
        singular_expression(equation, variables, prime)
        for equation in equations.values()
    )
    cutoff = max(equations)
    blocks = [
        "option(redSB);",
        f"ring R={prime},({','.join(map(str, ring_variables))}),dp;",
        f"ideal I={','.join(expressions)};",
        "timer=1;",
        "ideal G=slimgb(I);",
        "int elapsed=timer;",
        "int unit=0;",
        "if ((size(G)==1)&&(G[1]==1)) { unit=1; }",
        'print("RESULT_BEGIN");',
        f'if (unit==1) {{ print("cutoff={cutoff}"); }} else {{ print("cutoff=0"); }}',
        'print("unit="+string(unit));',
        'print("basis_size="+string(size(G)));',
        'print("dimension="+string(dim(G)));',
        'print("elapsed_ticks="+string(elapsed));',
        'print("BASIS_BEGIN"); print(G); print("BASIS_END");',
        'print("RESULT_END"); quit;',
    ]
    source = "\n".join(blocks) + "\n"
    record: dict[str, object] = {
        "prime": prime,
        "input_sha256": hashlib.sha256(source.encode()).hexdigest(),
    }
    started = time.monotonic()
    try:
        completed = subprocess.run(
            [singular, "-q"],
            input=source,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        record.update(
            {
                "status": "timeout",
                "seconds": time.monotonic() - started,
                "stdout_tail": (error.stdout or "")[-2000:],
                "stderr": (error.stderr or "")[-2000:],
            }
        )
        return record
    record["seconds"] = time.monotonic() - started
    record["returncode"] = completed.returncode
    has_result_markers = (
        "RESULT_BEGIN" in completed.stdout and "RESULT_END" in completed.stdout
    )
    record["status"] = (
        "completed"
        if completed.returncode == 0 and has_result_markers
        else "failed"
    )
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        for key in ("cutoff", "unit", "basis_size", "dimension", "elapsed_ticks"):
            if stripped.startswith(key + "="):
                record[key] = int(stripped.split("=", 1)[1])
    if "BASIS_BEGIN" in completed.stdout and "BASIS_END" in completed.stdout:
        basis = completed.stdout.split("BASIS_BEGIN", 1)[1].split("BASIS_END", 1)[0].strip()
        record["basis_sha256"] = hashlib.sha256(basis.encode()).hexdigest()
        record["basis_preview"] = basis[:4000]
    record["stdout_tail"] = completed.stdout[-2000:]
    record["stderr"] = completed.stderr[-2000:]
    return record


def exact_unit(
    equations: dict[int, sp.Expr],
    saturation: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    cutoff: int,
    timeout: int,
    threads: int,
) -> dict[str, object]:
    from jcsearch import msolve

    z = sp.Symbol("z")
    selected = [equation for order, equation in equations.items() if order <= cutoff]
    started = time.monotonic()
    result = msolve.run(
        selected + [z * saturation - 1],
        (z,) + variables,
        prime=0,
        threads=threads,
        groebner=True,
        timeout=timeout,
    )
    exact = (
        result.returncode == 0
        and result.empty
        and "#field characteristic: 0" in result.output
        and "#length of basis:      1 element" in result.output
        and result.output.rstrip().endswith("[1]:")
    )
    return {
        "cutoff": cutoff,
        "unit": int(exact),
        "returncode": result.returncode,
        "seconds": time.monotonic() - started,
        "output_sha256": hashlib.sha256(result.output.encode()).hexdigest(),
        "output": result.output[-1000:],
        "stderr": result.stderr[-1000:],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--degrees", type=parse_degrees, default=parse_degrees("2,4,6,8"))
    parser.add_argument("--groups", help="direction partition, for example 0,1|2|3")
    parser.add_argument("--max-order", type=int, default=10)
    parser.add_argument("--primes", nargs="+", type=int, default=list(PRIMES))
    parser.add_argument("--singular", default="Singular")
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--exact", action="store_true")
    parser.add_argument(
        "--linear-pivot",
        help="solve the first nonzero moment for this chart variable",
    )
    parser.add_argument(
        "--specialize",
        nargs="*",
        default=[],
        metavar="NAME=VALUE",
        help="specialize chart parameters or coefficients before elimination",
    )
    parser.add_argument("--msolve-threads", type=int, default=4)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()

    degrees = arguments.degrees
    groups = (
        parse_groups(arguments.groups, len(degrees))
        if arguments.groups
        else tuple((index,) for index in range(len(degrees)))
    )
    directions = len(groups)
    coefficients = sp.symbols(f"a0:{len(degrees)}")
    parameters = sp.symbols("lam mu")
    active_parameters = parameters[: max(directions - 3, 0)]
    chart_variables = coefficients[1:] + active_parameters
    substitutions = {coefficients[0]: 1}

    equations: dict[int, sp.Expr] = {}
    term_counts: dict[str, int] = {}
    compile_seconds: dict[str, float] = {}
    for order in range(2, arguments.max_order + 1):
        print(f"COMPILE moment {order}", flush=True)
        started = time.monotonic()
        compiled = moment(degrees, groups, order, coefficients, active_parameters)
        specialized = primitive_polynomial(compiled.subs(substitutions), chart_variables)
        compile_seconds[str(order)] = time.monotonic() - started
        if specialized != 0:
            equations[order] = specialized
            term_counts[str(order)] = len(sp.Poly(specialized, *chart_variables).terms())
        else:
            term_counts[str(order)] = 0
        print(
            f"COMPILED moment {order} terms={term_counts[str(order)]} "
            f"seconds={compile_seconds[str(order)]:.3f}",
            flush=True,
        )

    saturation = sp.prod(coefficients[1:]) * configuration_discriminant(
        directions, active_parameters
    )
    specialization: dict[sp.Symbol, sp.Rational] = {}
    for assignment in arguments.specialize:
        if "=" not in assignment:
            raise SystemExit(f"invalid specialization {assignment!r}")
        name, value_text = assignment.split("=", 1)
        matching = [variable for variable in chart_variables if str(variable) == name]
        if len(matching) != 1:
            raise SystemExit(f"unknown or ambiguous specialization variable {name!r}")
        specialization[matching[0]] = sp.Rational(value_text)
    if specialization:
        equations = {
            order: primitive_polynomial(
                equation.subs(specialization),
                tuple(variable for variable in chart_variables if variable not in specialization),
            )
            for order, equation in equations.items()
            if equation.subs(specialization) != 0
        }
        saturation = sp.expand(saturation.subs(specialization))
        chart_variables = tuple(
            variable for variable in chart_variables if variable not in specialization
        )
        print(
            "SPECIALIZED",
            {str(variable): str(value) for variable, value in specialization.items()},
            flush=True,
        )
    pivot_metadata = None
    if arguments.linear_pivot:
        matching = [
            variable
            for variable in chart_variables
            if str(variable) == arguments.linear_pivot
        ]
        if len(matching) != 1:
            raise SystemExit(f"unknown or ambiguous pivot {arguments.linear_pivot!r}")
        equations, chart_variables, saturation, pivot_metadata = linear_pivot_localization(
            equations,
            matching[0],
            chart_variables,
            saturation,
        )
        print(
            "LOCALIZED",
            arguments.linear_pivot,
            pivot_metadata["transformed_term_counts"],
            flush=True,
        )
    modular_results = []
    for prime in arguments.primes:
        print(f"RUN modular p={prime}", flush=True)
        result = modular_cutoff(
            equations,
            saturation,
            chart_variables,
            prime,
            arguments.singular,
            arguments.timeout,
        )
        modular_results.append(result)
        print(
            f"DONE p={prime} status={result.get('status')} "
            f"cutoff={result.get('cutoff')} dim={result.get('dimension')} "
            f"basis={result.get('basis_size')} seconds={result.get('seconds', 0):.3f}",
            flush=True,
        )

    positive_cutoffs = [
        int(result["cutoff"])
        for result in modular_results
        if result.get("status") == "completed" and int(result.get("cutoff", 0)) > 0
    ]
    exact_result = None
    if arguments.exact:
        if len(positive_cutoffs) != len(arguments.primes) or len(set(positive_cutoffs)) != 1:
            raise RuntimeError("exact promotion requires one common modular unit cutoff")
        cutoff = positive_cutoffs[0]
        print(f"RUN exact QQ cutoff={cutoff}", flush=True)
        exact_result = exact_unit(
            equations,
            saturation,
            chart_variables,
            cutoff,
            arguments.timeout,
            arguments.msolve_threads,
        )
        print(
            f"DONE exact unit={exact_result['unit']} seconds={exact_result['seconds']:.3f}",
            flush=True,
        )
        if exact_result["unit"] != 1:
            raise RuntimeError("exact characteristic-zero unit certificate failed")

    artifact = {
        "format": "gvc3-many-coherent-channels-research-v1",
        "status": (
            "exact characteristic-zero chart exclusion"
            if exact_result and exact_result["unit"] == 1
            else "bounded modular discovery; not a characteristic-zero theorem"
        ),
        "degrees": list(degrees),
        "groups": [list(block) for block in groups],
        "directions": directions,
        "normalization": "a0=1; first three direction points are infinity, zero, one",
        "configuration_parameters": [str(parameter) for parameter in active_parameters],
        "configuration_discriminant": str(
            configuration_discriminant(directions, active_parameters)
        ),
        "coefficient_saturation": [str(coefficient) for coefficient in coefficients[1:]],
        "linear_pivot_localization": pivot_metadata,
        "specialization": {
            str(variable): str(value) for variable, value in specialization.items()
        },
        "max_order": arguments.max_order,
        "moment_term_counts": term_counts,
        "compile_seconds": compile_seconds,
        "modular_results": modular_results,
        "exact_result": exact_result,
        "scope": "one declared direction-partition and nonzero-coefficient chart",
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
