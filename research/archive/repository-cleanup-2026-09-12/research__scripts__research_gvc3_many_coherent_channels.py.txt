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


def rational_polynomial_mod(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    prime: int,
) -> sp.Poly:
    """Reduce a rational-coefficient polynomial without rescaling it."""
    polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
    coefficients: dict[tuple[int, ...], int] = {}
    for monomial, coefficient in polynomial.terms():
        rational = sp.Rational(coefficient)
        denominator = int(rational.q) % prime
        if denominator == 0:
            raise ZeroDivisionError(
                f"coefficient denominator vanishes modulo {prime}"
            )
        coefficients[monomial] = (
            int(rational.p) * pow(denominator, -1, prime)
        ) % prime
    return sp.Poly.from_dict(coefficients, *variables, modulus=prime)


def modular_linear_specialization(
    polynomial: sp.Poly,
    pivot: sp.Symbol,
    substitution: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    prime: int,
) -> sp.Poly:
    """Substitute a linear boundary coordinate directly in a finite field."""
    boundary_variables = tuple(variable for variable in variables if variable != pivot)
    pivot_index = variables.index(pivot)
    substitution_polynomial = rational_polynomial_mod(
        substitution, boundary_variables, prime
    )
    pivot_degree = polynomial.degree(pivot)
    powers = [sp.Poly(1, *boundary_variables, modulus=prime)]
    for _ in range(pivot_degree):
        powers.append(powers[-1] * substitution_polynomial)
    answer = sp.Poly(0, *boundary_variables, modulus=prime)
    for monomial, coefficient in polynomial.terms():
        pivot_power = monomial[pivot_index]
        boundary_monomial = monomial[:pivot_index] + monomial[pivot_index + 1 :]
        coefficient_monomial = sp.Poly.from_dict(
            {boundary_monomial: int(coefficient)},
            *boundary_variables,
            modulus=prime,
        )
        answer += coefficient_monomial * powers[pivot_power]
    return answer


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
    coefficient_polynomial = sp.Poly(coefficient, *base_variables, domain=sp.QQ)
    remainder_polynomial = sp.Poly(remainder, *base_variables, domain=sp.QQ)
    transformed: dict[int, sp.Expr] = {}
    transformed_counts: dict[str, int] = {}
    for order, equation in equations.items():
        if order == first_order:
            continue
        polynomial = sp.Poly(equation, pivot)
        degree = polynomial.degree()
        cleared = sp.Poly(0, *base_variables, domain=sp.QQ)
        for power in range(degree + 1):
            part = polynomial.coeff_monomial(pivot**power)
            if part:
                cleared += (
                    sp.Poly(part, *base_variables, domain=sp.QQ)
                    * (-remainder_polynomial) ** power
                    * coefficient_polynomial ** (degree - power)
                )
        primitive_polynomial_data = cleared.clear_denoms()[1].primitive()[1]
        primitive = primitive_polynomial_data.as_expr()
        if primitive != 0:
            transformed[order] = primitive
            transformed_counts[str(order)] = len(primitive_polynomial_data.terms())
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


def singular_expression(
    expression: sp.Expr | sp.Poly,
    variables: tuple[sp.Symbol, ...],
    prime: int,
):
    if isinstance(expression, sp.Poly) and expression.gens == variables:
        polynomial = expression
    else:
        polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
    terms: list[str] = []
    for exponent, coefficient in polynomial.terms():
        if polynomial.domain.is_FiniteField:
            modular = int(coefficient) % prime
        else:
            rational = sp.Rational(coefficient)
            modular = (
                int(rational.p) * pow(int(rational.q), -1, prime)
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
        for key in (
            "cutoff",
            "unit",
            "basis_size",
            "dimension",
            "elapsed_ticks",
        ):
            if stripped.startswith(key + "="):
                record[key] = int(stripped.split("=", 1)[1])
    if "BASIS_BEGIN" in completed.stdout and "BASIS_END" in completed.stdout:
        basis = completed.stdout.split("BASIS_BEGIN", 1)[1].split("BASIS_END", 1)[0].strip()
        record["basis_sha256"] = hashlib.sha256(basis.encode()).hexdigest()
        record["basis_preview"] = basis[:4000]
    record["stdout_tail"] = completed.stdout[-2000:]
    record["stderr"] = completed.stderr[-2000:]
    return record


def modular_saturation_cutoff(
    equations: dict[int, sp.Expr],
    saturation: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    prime: int,
    singular: str,
    timeout: int,
    report_quotient_dimension: bool = False,
) -> dict[str, object]:
    """Compute the same localized ideal using quotient saturation.

    This avoids adjoining a Rabinowitsch variable.  On the larger coherent
    channel systems it can be substantially faster than ``modular_cutoff``.
    The output is still discovery in characteristic ``prime``; it is not an
    exact characteristic-zero certificate.
    """
    names: list[str] = []
    declarations: list[str] = []
    for order, equation in equations.items():
        name = f"f{order}"
        names.append(name)
        declarations.append(
            f"poly {name}={singular_expression(equation, variables, prime)};"
        )
    cutoff = max(equations)
    quotient_dimension_setup = (
        [
            "int quotient_dimension=-1;",
            "if (dim(G)==0) { quotient_dimension=vdim(G); }",
        ]
        if report_quotient_dimension
        else []
    )
    quotient_dimension_output = (
        ['print("quotient_dimension="+string(quotient_dimension));']
        if report_quotient_dimension
        else []
    )
    source = "\n".join(
        [
            'LIB "elim.lib";',
            "option(redSB);",
            f"ring R={prime},({','.join(map(str, variables))}),dp;",
            *declarations,
            f"poly jsat={singular_expression(saturation, variables, prime)};",
            f"ideal I={','.join(names)};",
            "ideal J=jsat;",
            "timer=1;",
            "list L=sat_with_exp(I,J);",
            "ideal G=L[1];",
            "int elapsed=timer;",
            "int unit=0;",
            "if ((size(G)==1)&&(G[1]==1)) { unit=1; }",
            *quotient_dimension_setup,
            'print("RESULT_BEGIN");',
            f'if (unit==1) {{ print("cutoff={cutoff}"); }} else {{ print("cutoff=0"); }}',
            'print("unit="+string(unit));',
            'print("basis_size="+string(size(G)));',
            'print("dimension="+string(dim(G)));',
            *quotient_dimension_output,
            'print("saturation_exponent="+string(L[2]));',
            'print("elapsed_ticks="+string(elapsed));',
            'print("BASIS_BEGIN"); print(G); print("BASIS_END");',
            'print("RESULT_END"); quit;',
        ]
    ) + "\n"
    record: dict[str, object] = {
        "prime": prime,
        "method": "quotient_saturation",
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
        for key in (
            "cutoff",
            "unit",
            "basis_size",
            "dimension",
            "quotient_dimension",
            "saturation_exponent",
            "elapsed_ticks",
        ):
            if stripped.startswith(key + "="):
                record[key] = int(stripped.split("=", 1)[1])
    if "BASIS_BEGIN" in completed.stdout and "BASIS_END" in completed.stdout:
        basis = (
            completed.stdout.split("BASIS_BEGIN", 1)[1]
            .split("BASIS_END", 1)[0]
            .strip()
        )
        record["basis_sha256"] = hashlib.sha256(basis.encode()).hexdigest()
        record["basis_preview"] = basis[:4000]
    record["stdout_tail"] = completed.stdout[-2000:]
    record["stderr"] = completed.stderr[-2000:]
    return record


def saturation_factor_records(
    saturation: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> list[dict[str, object]]:
    """Factor a localization product over QQ in a stable, sparse-first order."""
    _, factor_data = sp.factor_list(
        sp.Poly(saturation, *variables, domain=sp.QQ)
    )
    variable_names = {str(variable) for variable in variables}
    records: list[dict[str, object]] = []
    for factor_polynomial, multiplicity in factor_data:
        expression = primitive_polynomial(factor_polynomial.as_expr(), variables)
        polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
        terms = len(polynomial.terms())
        degree = polynomial.total_degree()
        text = str(expression)
        if text in variable_names:
            category = 0
        elif degree == 1:
            category = 1
        else:
            category = 2
        records.append(
            {
                "expression": expression,
                "expression_text": text,
                "sha256": hashlib.sha256(text.encode()).hexdigest(),
                "multiplicity": multiplicity,
                "total_degree": degree,
                "terms": terms,
                "sort_key": (category, degree, terms, text),
            }
        )
    records.sort(key=lambda record: record["sort_key"])
    return records


def modular_principal_gcd_step(
    polynomial: sp.Poly,
    factor: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    prime: int,
    singular: str,
    timeout: int,
) -> tuple[dict[str, object], str]:
    """Remove all common factors from one principal generator in Singular."""
    source = "\n".join(
        [
            f"ring R={prime},({','.join(map(str, variables))}),dp;",
            (
                "poly f="
                + singular_expression(polynomial, variables, prime)
                + ";"
            ),
            f"poly g={singular_expression(factor, variables, prime)};",
            "timer=1;",
            "poly h=gcd(f,g);",
            "int saturation_exponent=0;",
            "while (deg(h)>0)",
            "{",
            "  f=f/h;",
            "  saturation_exponent=saturation_exponent+1;",
            "  h=gcd(f,g);",
            "}",
            "int elapsed=timer;",
            "int unit=0;",
            "if (deg(f)==0) { unit=1; }",
            'print("RESULT_BEGIN");',
            'print("unit="+string(unit));',
            'print("saturation_exponent="+string(saturation_exponent));',
            'print("elapsed_ticks="+string(elapsed));',
            'print("BASIS_BEGIN"); print(f); print("BASIS_END");',
            'print("RESULT_END"); quit;',
        ]
    ) + "\n"
    record: dict[str, object] = {
        "engine": "principal_singular_gcd",
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
        return record, ""
    record["seconds"] = time.monotonic() - started
    record["returncode"] = completed.returncode
    markers = "RESULT_BEGIN" in completed.stdout and "RESULT_END" in completed.stdout
    record["status"] = (
        "completed" if completed.returncode == 0 and markers else "failed"
    )
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        for key in ("unit", "saturation_exponent", "elapsed_ticks"):
            if stripped.startswith(key + "="):
                record[key] = int(stripped.split("=", 1)[1])
    basis = ""
    if "BASIS_BEGIN" in completed.stdout and "BASIS_END" in completed.stdout:
        basis = (
            completed.stdout.split("BASIS_BEGIN", 1)[1]
            .split("BASIS_END", 1)[0]
            .strip()
        )
        record["basis_sha256"] = hashlib.sha256(basis.encode()).hexdigest()
        record["basis_preview"] = basis[:4000]
    record["stdout_tail"] = completed.stdout[-2000:]
    record["stderr"] = completed.stderr[-2000:]
    return record, basis


def modular_pair_linear_saturation_test(
    polynomials: list[sp.Poly],
    factor: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    prime: int,
    singular: str,
    timeout: int,
) -> dict[str, object]:
    """Test saturation of a two-generator complete intersection by a line.

    If the two generators are coprime and remain coprime on ``factor=0``,
    their codimension-two complete intersection has no component supported on
    that line, hence is already saturated by the factor.
    """
    support = [variable for variable in variables if variable in factor.free_symbols]
    if not support:
        return {"status": "failed", "reason": "constant factor"}
    pivot = support[-1]
    factor_in_pivot = sp.Poly(factor, pivot)
    if factor_in_pivot.degree() != 1:
        return {"status": "failed", "reason": "factor is not linear in pivot"}
    coefficient = factor_in_pivot.coeff_monomial(pivot)
    remainder = factor_in_pivot.coeff_monomial(1)
    if coefficient.free_symbols:
        return {"status": "failed", "reason": "nonconstant pivot coefficient"}
    substitution = sp.cancel(-remainder / coefficient)
    boundary_variables = tuple(variable for variable in variables if variable != pivot)
    boundary_polynomials = [
        modular_linear_specialization(
            polynomial, pivot, substitution, variables, prime
        )
        for polynomial in polynomials
    ]
    source = "\n".join(
        [
            f"ring R={prime},({','.join(map(str, variables))}),dp;",
            "short=0;",
            f"poly f={singular_expression(polynomials[0], variables, prime)};",
            f"poly g={singular_expression(polynomials[1], variables, prime)};",
            "timer=1;",
            "poly h=gcd(f,g);",
            "int full_gcd_degree=deg(h);",
            f"ring S={prime},({','.join(map(str, boundary_variables))}),dp;",
            "short=0;",
            (
                "poly fb="
                + singular_expression(
                    boundary_polynomials[0], boundary_variables, prime
                )
                + ";"
            ),
            (
                "poly gb="
                + singular_expression(
                    boundary_polynomials[1], boundary_variables, prime
                )
                + ";"
            ),
            "poly hb=gcd(fb,gb);",
            "int boundary_gcd_degree=deg(hb);",
            "int elapsed=timer;",
            'print("RESULT_BEGIN");',
            'print("full_gcd_degree="+string(full_gcd_degree));',
            'print("boundary_gcd_degree="+string(boundary_gcd_degree));',
            'print("elapsed_ticks="+string(elapsed));',
            'print("BOUNDARY_GCD_BEGIN"); print(hb); print("BOUNDARY_GCD_END");',
            'print("RESULT_END"); quit;',
        ]
    ) + "\n"
    record: dict[str, object] = {
        "engine": "two_generator_linear_boundary_gcd",
        "pivot": str(pivot),
        "substitution": str(substitution),
        "boundary_term_counts": [
            len(polynomial.terms())
            for polynomial in boundary_polynomials
        ],
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
    markers = "RESULT_BEGIN" in completed.stdout and "RESULT_END" in completed.stdout
    record["status"] = (
        "completed" if completed.returncode == 0 and markers else "failed"
    )
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        for key in ("full_gcd_degree", "boundary_gcd_degree", "elapsed_ticks"):
            if stripped.startswith(key + "="):
                record[key] = int(stripped.split("=", 1)[1])
    if "BOUNDARY_GCD_BEGIN" in completed.stdout and "BOUNDARY_GCD_END" in completed.stdout:
        boundary_gcd = (
            completed.stdout.split("BOUNDARY_GCD_BEGIN", 1)[1]
            .split("BOUNDARY_GCD_END", 1)[0]
            .strip()
        )
        record["boundary_gcd_sha256"] = hashlib.sha256(
            boundary_gcd.encode()
        ).hexdigest()
        record["boundary_gcd_preview"] = boundary_gcd[:4000]
        record["_boundary_gcd_full"] = boundary_gcd
    record["saturated"] = int(
        record.get("status") == "completed"
        and record.get("full_gcd_degree") == 0
        and record.get("boundary_gcd_degree") == 0
    )
    record["stdout_tail"] = completed.stdout[-2000:]
    record["stderr"] = completed.stderr[-2000:]
    return record


def two_generator_linear_colon_generator(
    polynomials: list[sp.Poly],
    factor: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    prime: int,
    boundary_test: dict[str, object],
) -> sp.Poly:
    """Construct (g*F-f*G)/L from the common boundary factor H."""
    pivot = next(
        variable for variable in variables if str(variable) == boundary_test["pivot"]
    )
    substitution = sp.sympify(
        str(boundary_test["substitution"]),
        locals={str(variable): variable for variable in variables},
    )
    boundary_variables = tuple(variable for variable in variables if variable != pivot)
    boundary_polynomials = [
        modular_linear_specialization(
            polynomial, pivot, substitution, variables, prime
        )
        for polynomial in polynomials
    ]
    gcd_text = str(boundary_test["_boundary_gcd_full"])
    common = sp.Poly(
        sp.sympify(
            gcd_text.replace("^", "**"),
            locals={str(variable): variable for variable in boundary_variables},
        ),
        *boundary_variables,
        modulus=prime,
    ).monic()
    first_residual = boundary_polynomials[0].exquo(common)
    second_residual = boundary_polynomials[1].exquo(common)
    first_lift = sp.Poly(
        first_residual.as_expr(), *variables, modulus=prime
    )
    second_lift = sp.Poly(
        second_residual.as_expr(), *variables, modulus=prime
    )
    numerator = second_lift * polynomials[0] - first_lift * polynomials[1]
    factor_polynomial = sp.Poly(factor, *variables, modulus=prime)
    quotient, remainder = numerator.div(factor_polynomial)
    if not remainder.is_zero:
        raise RuntimeError("boundary syzygy numerator is not divisible by its line")
    return quotient.monic()


def modular_residual_boundary_gcd_test(
    colon_generator: sp.Poly,
    variables: tuple[sp.Symbol, ...],
    prime: int,
    singular: str,
    timeout: int,
    boundary_test: dict[str, object],
) -> dict[str, object]:
    """Test whether the first linear colon removed the whole boundary factor."""
    pivot = next(
        variable for variable in variables if str(variable) == boundary_test["pivot"]
    )
    substitution = sp.sympify(
        str(boundary_test["substitution"]),
        locals={str(variable): variable for variable in variables},
    )
    boundary_variables = tuple(variable for variable in variables if variable != pivot)
    colon_boundary = modular_linear_specialization(
        colon_generator, pivot, substitution, variables, prime
    )
    common_text = str(boundary_test["_boundary_gcd_full"])
    source = "\n".join(
        [
            f"ring R={prime},({','.join(map(str, boundary_variables))}),dp;",
            "short=0;",
            f"poly h={common_text};",
            (
                "poly k="
                + singular_expression(
                    colon_boundary, boundary_variables, prime
                )
                + ";"
            ),
            "timer=1;",
            "poly residual=gcd(h,k);",
            "int residual_gcd_degree=deg(residual);",
            "int elapsed=timer;",
            'print("RESULT_BEGIN");',
            'print("residual_gcd_degree="+string(residual_gcd_degree));',
            'print("elapsed_ticks="+string(elapsed));',
            'print("RESIDUAL_GCD_BEGIN"); print(residual); print("RESIDUAL_GCD_END");',
            'print("RESULT_END"); quit;',
        ]
    ) + "\n"
    record: dict[str, object] = {
        "engine": "first_colon_residual_boundary_gcd",
        "colon_boundary_terms": len(colon_boundary.terms()),
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
    markers = "RESULT_BEGIN" in completed.stdout and "RESULT_END" in completed.stdout
    record["status"] = (
        "completed" if completed.returncode == 0 and markers else "failed"
    )
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        for key in ("residual_gcd_degree", "elapsed_ticks"):
            if stripped.startswith(key + "="):
                record[key] = int(stripped.split("=", 1)[1])
    if "RESIDUAL_GCD_BEGIN" in completed.stdout and "RESIDUAL_GCD_END" in completed.stdout:
        residual = (
            completed.stdout.split("RESIDUAL_GCD_BEGIN", 1)[1]
            .split("RESIDUAL_GCD_END", 1)[0]
            .strip()
        )
        record["residual_gcd_sha256"] = hashlib.sha256(residual.encode()).hexdigest()
        record["residual_gcd_preview"] = residual[:4000]
    record["saturated_after_first_colon"] = int(
        record.get("status") == "completed"
        and record.get("residual_gcd_degree") == 0
    )
    record["stdout_tail"] = completed.stdout[-2000:]
    record["stderr"] = completed.stderr[-2000:]
    return record


def modular_linear_boundary_gcd_test(
    polynomials: list[sp.Poly],
    factor: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    prime: int,
    singular: str,
    timeout: int,
) -> dict[str, object]:
    """Compute the common boundary divisor for a linear factor.

    The successive saturation ideals start as a codimension-two complete
    intersection and every colon remains unmixed of the same codimension.
    Consequently a unit boundary gcd certifies that the linear factor avoids
    every associated component.
    """
    support = [variable for variable in variables if variable in factor.free_symbols]
    if not support:
        return {"status": "failed", "reason": "constant factor"}
    pivot = support[-1]
    factor_in_pivot = sp.Poly(factor, pivot)
    if factor_in_pivot.degree() != 1:
        return {"status": "failed", "reason": "factor is not linear in pivot"}
    coefficient = factor_in_pivot.coeff_monomial(pivot)
    remainder = factor_in_pivot.coeff_monomial(1)
    if coefficient.free_symbols:
        return {"status": "failed", "reason": "nonconstant pivot coefficient"}
    substitution = sp.cancel(-remainder / coefficient)
    boundary_variables = tuple(variable for variable in variables if variable != pivot)
    boundary_polynomials = [
        modular_linear_specialization(
            polynomial, pivot, substitution, variables, prime
        )
        for polynomial in polynomials
    ]
    source_lines = [
        f"ring R={prime},({','.join(map(str, boundary_variables))}),dp;",
        "short=0;",
    ]
    for index, polynomial in enumerate(boundary_polynomials, start=1):
        source_lines.append(
                f"poly f{index}="
            + singular_expression(polynomial, boundary_variables, prime)
            + ";"
        )
    source_lines.extend(
        [
            "timer=1;",
            "poly h=f1;",
            *[
                f"h=gcd(h,f{index});"
                for index in range(2, len(boundary_polynomials) + 1)
            ],
            "int boundary_gcd_degree=deg(h);",
            "int elapsed=timer;",
            'print("RESULT_BEGIN");',
            'print("boundary_gcd_degree="+string(boundary_gcd_degree));',
            'print("elapsed_ticks="+string(elapsed));',
            'print("BOUNDARY_GCD_BEGIN"); print(h); print("BOUNDARY_GCD_END");',
            'print("RESULT_END"); quit;',
        ]
    )
    source = "\n".join(source_lines) + "\n"
    record: dict[str, object] = {
        "engine": "unmixed_linear_boundary_gcd",
        "pivot": str(pivot),
        "substitution": str(substitution),
        "boundary_term_counts": [
            len(polynomial.terms())
            for polynomial in boundary_polynomials
        ],
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
    markers = "RESULT_BEGIN" in completed.stdout and "RESULT_END" in completed.stdout
    record["status"] = (
        "completed" if completed.returncode == 0 and markers else "failed"
    )
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        for key in ("boundary_gcd_degree", "elapsed_ticks"):
            if stripped.startswith(key + "="):
                record[key] = int(stripped.split("=", 1)[1])
    if "BOUNDARY_GCD_BEGIN" in completed.stdout and "BOUNDARY_GCD_END" in completed.stdout:
        boundary_gcd = (
            completed.stdout.split("BOUNDARY_GCD_BEGIN", 1)[1]
            .split("BOUNDARY_GCD_END", 1)[0]
            .strip()
        )
        record["boundary_gcd_sha256"] = hashlib.sha256(
            boundary_gcd.encode()
        ).hexdigest()
        record["boundary_gcd_preview"] = boundary_gcd[:4000]
        record["_boundary_gcd_full"] = boundary_gcd
    record["saturated"] = int(
        record.get("status") == "completed"
        and record.get("boundary_gcd_degree") == 0
    )
    record["stdout_tail"] = completed.stdout[-2000:]
    record["stderr"] = completed.stderr[-2000:]
    return record


def boundary_gcd_supported_on_factors(
    boundary_test: dict[str, object],
    supporting_factors: list[dict[str, object]],
    variables: tuple[sp.Symbol, ...],
    prime: int,
) -> tuple[bool, dict[str, int]]:
    """Test whether a boundary gcd uses only declared localization factors."""
    gcd_text = boundary_test.get("_boundary_gcd_full")
    if not gcd_text or not supporting_factors:
        return False, {}
    pivot = next(
        variable for variable in variables if str(variable) == boundary_test["pivot"]
    )
    substitution = sp.sympify(
        str(boundary_test["substitution"]),
        locals={str(variable): variable for variable in variables},
    )
    boundary_variables = tuple(variable for variable in variables if variable != pivot)
    remaining = sp.Poly(
        sp.sympify(
            str(gcd_text).replace("^", "**"),
            locals={str(variable): variable for variable in boundary_variables},
        ),
        *boundary_variables,
        modulus=prime,
    ).monic()
    exponents: dict[str, int] = {}
    restricted_factors: list[tuple[int, dict[str, object], sp.Poly]] = []
    for supporting_factor in supporting_factors:
        divisor = rational_polynomial_mod(
            supporting_factor["expression"].subs(pivot, substitution),
            boundary_variables,
            prime,
        )
        if divisor.is_zero or divisor.total_degree() == 0:
            continue
        divisor = divisor.monic()
        restricted_factors.append(
            (divisor.total_degree(), supporting_factor, divisor)
        )
    for _, supporting_factor, divisor in sorted(
        restricted_factors, key=lambda item: item[0], reverse=True
    ):
        exponent = 0
        quotient, remainder = remaining.div(divisor)
        while remainder.is_zero:
            remaining = quotient
            exponent += 1
            quotient, remainder = remaining.div(divisor)
        if exponent:
            exponents[str(supporting_factor["expression_text"])] = exponent
    return remaining.total_degree() == 0, exponents


def modular_unmixed_hypersurface_dimension_test(
    polynomials: list[sp.Poly],
    factor: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    prime: int,
    singular: str,
    timeout: int,
) -> dict[str, object]:
    """Test whether a factor contains a component of an unmixed ideal."""
    support = set(factor.free_symbols)
    ring_variables = tuple(
        variable for variable in variables if variable not in support
    ) + tuple(variable for variable in variables if variable in support)
    source = "\n".join(
        [
            "option(redSB);",
            f"ring R={prime},({','.join(map(str, ring_variables))}),dp;",
            "ideal I="
            + ",".join(
                singular_expression(polynomial, variables, prime)
                for polynomial in polynomials
            )
            + ";",
            "poly h=" + singular_expression(factor, variables, prime) + ";",
            "timer=1;",
            "ideal G=std(I+ideal(h));",
            "int boundary_dimension=dim(G);",
            "int elapsed=timer;",
            'print("RESULT_BEGIN");',
            'print("boundary_dimension="+string(boundary_dimension));',
            'print("elapsed_ticks="+string(elapsed));',
            'print("RESULT_END"); quit;',
        ]
    ) + "\n"
    record: dict[str, object] = {
        "engine": "unmixed_hypersurface_dimension",
        "input_generators": len(polynomials),
        "ring_variable_order": [str(variable) for variable in ring_variables],
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
    markers = "RESULT_BEGIN" in completed.stdout and "RESULT_END" in completed.stdout
    record["status"] = (
        "completed" if completed.returncode == 0 and markers else "failed"
    )
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        for key in ("boundary_dimension", "elapsed_ticks"):
            if stripped.startswith(key + "="):
                record[key] = int(stripped.split("=", 1)[1])
    record["avoids_components"] = int(
        record.get("status") == "completed"
        and int(record.get("boundary_dimension", len(variables)))
        < len(variables) - 2
    )
    record["stdout_tail"] = completed.stdout[-2000:]
    record["stderr"] = completed.stderr[-2000:]
    return record


def modular_linear_colon_syzygy_step(
    polynomials: list[sp.Poly],
    factor: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    prime: int,
    singular: str,
    timeout: int,
    boundary_test: dict[str, object],
) -> tuple[dict[str, object], list[sp.Poly]]:
    """Generate an exact linear colon from boundary syzygies.

    If ``P_i|_{L=0}=H r_i``, every syzygy of the residual tuple ``(r_i)``
    lifts to ``sum(s_i P_i)/L``.  Lifting a generating set of syzygies gives
    the full ideal quotient ``(P_i):L``.
    """
    pivot = next(
        variable for variable in variables if str(variable) == boundary_test["pivot"]
    )
    substitution = sp.sympify(
        str(boundary_test["substitution"]),
        locals={str(variable): variable for variable in variables},
    )
    boundary_variables = tuple(variable for variable in variables if variable != pivot)
    boundary_polynomials = [
        modular_linear_specialization(
            polynomial, pivot, substitution, variables, prime
        )
        for polynomial in polynomials
    ]
    common = sp.Poly(
        sp.sympify(
            str(boundary_test["_boundary_gcd_full"]).replace("^", "**"),
            locals={str(variable): variable for variable in boundary_variables},
        ),
        *boundary_variables,
        modulus=prime,
    ).monic()
    residuals = [polynomial.exquo(common) for polynomial in boundary_polynomials]
    ring_variables = tuple(
        variable for variable in variables if variable != pivot
    ) + (pivot,)
    source = "\n".join(
        [
            f"ring R={prime},({','.join(map(str, ring_variables))}),dp;",
            "short=0;",
            "ideal I="
            + ",".join(
                singular_expression(polynomial, variables, prime)
                for polynomial in polynomials
            )
            + ";",
            "ideal B="
            + ",".join(
                singular_expression(residual, boundary_variables, prime)
                for residual in residuals
            )
            + ";",
            "poly ell=" + singular_expression(factor, variables, prime) + ";",
            "timer=1;",
            "module Z=syz(B);",
            "int nonzero_generators=0;",
            "for (int i=1;i<=size(Z);i++) {",
            "  poly q=0;",
            "  for (int j=1;j<=size(I);j++) { q=q+Z[i][j]*I[j]; }",
            "  q=q/ell;",
            "  if (q!=0) {",
            "    nonzero_generators=nonzero_generators+1;",
            '    print("COLON_GENERATOR_BEGIN"); print(q); print("COLON_GENERATOR_END");',
            "  }",
            "}",
            "int elapsed=timer;",
            'print("RESULT_BEGIN");',
            'print("syzygies="+string(size(Z)));',
            'print("nonzero_generators="+string(nonzero_generators));',
            'print("elapsed_ticks="+string(elapsed));',
            'print("RESULT_END"); quit;',
        ]
    ) + "\n"
    record: dict[str, object] = {
        "engine": "linear_boundary_syzygy_colon",
        "input_generators": len(polynomials),
        "boundary_gcd_degree": boundary_test.get("boundary_gcd_degree"),
        "boundary_gcd_sha256": boundary_test.get("boundary_gcd_sha256"),
        "residual_term_counts": [len(residual.terms()) for residual in residuals],
        "ring_variable_order": [str(variable) for variable in ring_variables],
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
        return record, []
    record["seconds"] = time.monotonic() - started
    record["returncode"] = completed.returncode
    markers = "RESULT_BEGIN" in completed.stdout and "RESULT_END" in completed.stdout
    record["status"] = (
        "completed" if completed.returncode == 0 and markers else "failed"
    )
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        for key in ("syzygies", "nonzero_generators", "elapsed_ticks"):
            if stripped.startswith(key + "="):
                record[key] = int(stripped.split("=", 1)[1])
    generator_texts: list[str] = []
    remaining = completed.stdout
    while "COLON_GENERATOR_BEGIN" in remaining:
        remaining = remaining.split("COLON_GENERATOR_BEGIN", 1)[1]
        generator_text, remaining = remaining.split("COLON_GENERATOR_END", 1)
        generator_texts.append(generator_text.strip())
    generators = [
        sp.Poly(
            sp.sympify(
                generator.replace("^", "**"),
                locals={str(variable): variable for variable in variables},
            ),
            *variables,
            modulus=prime,
        ).monic()
        for generator in generator_texts
    ]
    record["generator_term_counts"] = [len(generator.terms()) for generator in generators]
    record["generator_total_degrees"] = [
        generator.total_degree() for generator in generators
    ]
    record["generator_sha256"] = [
        hashlib.sha256(
            singular_expression(generator, variables, prime).encode()
        ).hexdigest()
        for generator in generators
    ]
    record["stdout_tail"] = completed.stdout[-2000:]
    record["stderr"] = completed.stderr[-2000:]
    return record, generators


def modular_structured_triple_linear_colon_step(
    polynomials: list[sp.Poly],
    factor: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    prime: int,
    singular: str,
    timeout: int,
    boundary_test: dict[str, object],
) -> tuple[dict[str, object], list[sp.Poly]]:
    """Certify a three-generator linear colon by a regular residual triple.

    Write the boundary residuals after their common gcd as ``S*f,S*g,c``.
    When ``gcd(S,c)=1`` and ``f,g,c`` is a regular sequence, the three
    displayed pairwise syzygies generate the full residual syzygy module:
    ``(g,-f,0)``, ``(c,0,-S*f)``, and ``(0,c,-S*g)``.
    """
    if len(polynomials) != 3:
        return {"status": "failed", "reason": "requires three generators"}, []
    pivot = next(
        variable for variable in variables if str(variable) == boundary_test["pivot"]
    )
    substitution = sp.sympify(
        str(boundary_test["substitution"]),
        locals={str(variable): variable for variable in variables},
    )
    boundary_variables = tuple(variable for variable in variables if variable != pivot)
    boundary_polynomials = [
        modular_linear_specialization(
            polynomial, pivot, substitution, variables, prime
        )
        for polynomial in polynomials
    ]
    common = sp.Poly(
        sp.sympify(
            str(boundary_test["_boundary_gcd_full"]).replace("^", "**"),
            locals={str(variable): variable for variable in boundary_variables},
        ),
        *boundary_variables,
        modulus=prime,
    ).monic()
    residuals = [polynomial.exquo(common) for polynomial in boundary_polynomials]
    source = "\n".join(
        [
            f"ring R={prime},({','.join(map(str, boundary_variables))}),dp;",
            "short=0;",
            "poly b1="
            + singular_expression(residuals[0], boundary_variables, prime)
            + ";",
            "poly b2="
            + singular_expression(residuals[1], boundary_variables, prime)
            + ";",
            "poly c="
            + singular_expression(residuals[2], boundary_variables, prime)
            + ";",
            "timer=1;",
            "poly s=gcd(b1,b2);",
            "poly f=b1/s;",
            "poly g=b2/s;",
            "poly sc=gcd(s,c);",
            "ideal FG=std(ideal(f,g));",
            "poly membership_remainder=reduce(c,FG);",
            "int third_in_first_pair=0;",
            "if (membership_remainder==0) { third_in_first_pair=1; }",
            "ideal T=f,g,c;",
            "ideal G=std(T);",
            "int regular_dimension=dim(G);",
            "int elapsed=timer;",
            'print("RESULT_BEGIN");',
            'print("shared_degree="+string(deg(s)));',
            'print("shared_third_gcd_degree="+string(deg(sc)));',
            'print("regular_dimension="+string(regular_dimension));',
            'print("third_in_first_pair="+string(third_in_first_pair));',
            'print("elapsed_ticks="+string(elapsed));',
            'print("SHARED_BEGIN"); print(s); print("SHARED_END");',
            'print("FIRST_BEGIN"); print(f); print("FIRST_END");',
            'print("SECOND_BEGIN"); print(g); print("SECOND_END");',
            'print("REMAINDER_BEGIN"); print(membership_remainder); print("REMAINDER_END");',
            'print("RESULT_END"); quit;',
        ]
    ) + "\n"
    record: dict[str, object] = {
        "engine": "structured_triple_regular_sequence_colon",
        "boundary_gcd_degree": boundary_test.get("boundary_gcd_degree"),
        "boundary_gcd_sha256": boundary_test.get("boundary_gcd_sha256"),
        "residual_term_counts": [len(residual.terms()) for residual in residuals],
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
        return record, []
    record["seconds"] = time.monotonic() - started
    record["returncode"] = completed.returncode
    markers = "RESULT_BEGIN" in completed.stdout and "RESULT_END" in completed.stdout
    record["status"] = (
        "completed" if completed.returncode == 0 and markers else "failed"
    )
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        for key in (
            "shared_degree",
            "shared_third_gcd_degree",
            "regular_dimension",
            "third_in_first_pair",
            "elapsed_ticks",
        ):
            if stripped.startswith(key + "="):
                record[key] = int(stripped.split("=", 1)[1])
    marked_polynomials: dict[str, sp.Poly] = {}
    for label in ("SHARED", "FIRST", "SECOND"):
        begin = label + "_BEGIN"
        end = label + "_END"
        if begin in completed.stdout and end in completed.stdout:
            polynomial_text = (
                completed.stdout.split(begin, 1)[1].split(end, 1)[0].strip()
            )
            marked_polynomials[label] = sp.Poly(
                sp.sympify(
                    polynomial_text.replace("^", "**"),
                    locals={
                        str(variable): variable for variable in boundary_variables
                    },
                ),
                *boundary_variables,
                modulus=prime,
            ).monic()
    for label, polynomial in marked_polynomials.items():
        expression_text = singular_expression(
            polynomial, boundary_variables, prime
        )
        record[label.lower() + "_terms"] = len(polynomial.terms())
        record[label.lower() + "_total_degree"] = polynomial.total_degree()
        record[label.lower() + "_sha256"] = hashlib.sha256(
            expression_text.encode()
        ).hexdigest()
        record[label.lower() + "_preview"] = expression_text[:4000]
    if "REMAINDER_BEGIN" in completed.stdout and "REMAINDER_END" in completed.stdout:
        remainder_text = (
            completed.stdout.split("REMAINDER_BEGIN", 1)[1]
            .split("REMAINDER_END", 1)[0]
            .strip()
        )
        record["membership_remainder_sha256"] = hashlib.sha256(
            remainder_text.encode()
        ).hexdigest()
        record["membership_remainder_preview"] = remainder_text[:4000]
    expected_dimension = len(boundary_variables) - 3
    certified = int(
        record.get("status") == "completed"
        and int(record.get("shared_degree", 0)) > 0
        and record.get("shared_third_gcd_degree") == 0
        and record.get("regular_dimension") == expected_dimension
        and set(marked_polynomials) == {"SHARED", "FIRST", "SECOND"}
    )
    record["expected_regular_dimension"] = expected_dimension
    record["certified"] = certified
    generators: list[sp.Poly] = []
    if certified:
        shared = sp.Poly(
            marked_polynomials["SHARED"].as_expr(), *variables, modulus=prime
        )
        first = sp.Poly(
            marked_polynomials["FIRST"].as_expr(), *variables, modulus=prime
        )
        second = sp.Poly(
            marked_polynomials["SECOND"].as_expr(), *variables, modulus=prime
        )
        third = sp.Poly(residuals[2].as_expr(), *variables, modulus=prime)
        first_residual = shared * first
        second_residual = shared * second
        syzygies = [
            (second, -first, sp.Poly(0, *variables, modulus=prime)),
            (third, sp.Poly(0, *variables, modulus=prime), -first_residual),
            (sp.Poly(0, *variables, modulus=prime), third, -second_residual),
        ]
        factor_polynomial = sp.Poly(factor, *variables, modulus=prime)
        for syzygy in syzygies:
            numerator = sum(
                (coefficient * polynomial for coefficient, polynomial in zip(syzygy, polynomials)),
                sp.Poly(0, *variables, modulus=prime),
            )
            quotient, remainder = numerator.div(factor_polynomial)
            if not remainder.is_zero:
                raise RuntimeError("structured boundary syzygy is not divisible")
            if not quotient.is_zero:
                generators.append(quotient.monic())
    record["generator_term_counts"] = [len(generator.terms()) for generator in generators]
    record["generator_total_degrees"] = [
        generator.total_degree() for generator in generators
    ]
    record["generator_sha256"] = [
        hashlib.sha256(
            singular_expression(generator, variables, prime).encode()
        ).hexdigest()
        for generator in generators
    ]
    record["stdout_tail"] = completed.stdout[-2000:]
    record["stderr"] = completed.stderr[-2000:]
    return record, generators


def modular_component_quotient_step(
    polynomials: list[sp.Poly],
    component_generators: tuple[sp.Expr, sp.Expr],
    variables: tuple[sp.Symbol, ...],
    prime: int,
    singular: str,
    timeout: int,
) -> tuple[dict[str, object], list[sp.Poly], str]:
    """Take one exact quotient by a two-generator supported component."""
    support = set().union(
        *(generator.free_symbols for generator in component_generators)
    )
    ring_variables = tuple(
        variable for variable in variables if variable not in support
    ) + tuple(variable for variable in variables if variable in support)
    source = "\n".join(
        [
            "option(redSB);",
            f"ring R={prime},({','.join(map(str, ring_variables))}),dp;",
            "short=0;",
            "ideal I="
            + ",".join(
                singular_expression(polynomial, variables, prime)
                for polynomial in polynomials
            )
            + ";",
            "ideal C="
            + ",".join(
                singular_expression(generator, variables, prime)
                for generator in component_generators
            )
            + ";",
            "timer=1;",
            "ideal G=std(quotient(I,C));",
            "int elapsed=timer;",
            "int unit=0;",
            "if ((size(G)==1)&&(G[1]==1)) { unit=1; }",
            'print("RESULT_BEGIN");',
            'print("unit="+string(unit));',
            'print("basis_size="+string(size(G)));',
            'print("dimension="+string(dim(G)));',
            'print("elapsed_ticks="+string(elapsed));',
            'print("BASIS_BEGIN"); print(G); print("BASIS_END");',
            'print("RESULT_END"); quit;',
        ]
    ) + "\n"
    record: dict[str, object] = {
        "engine": "two_generator_component_quotient",
        "component_generators": [str(generator) for generator in component_generators],
        "input_generators": len(polynomials),
        "ring_variable_order": [str(variable) for variable in ring_variables],
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
        return record, [], ""
    record["seconds"] = time.monotonic() - started
    record["returncode"] = completed.returncode
    markers = "RESULT_BEGIN" in completed.stdout and "RESULT_END" in completed.stdout
    record["status"] = (
        "completed" if completed.returncode == 0 and markers else "failed"
    )
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        for key in ("unit", "basis_size", "dimension", "elapsed_ticks"):
            if stripped.startswith(key + "="):
                record[key] = int(stripped.split("=", 1)[1])
    basis = ""
    if "BASIS_BEGIN" in completed.stdout and "BASIS_END" in completed.stdout:
        basis = (
            completed.stdout.split("BASIS_BEGIN", 1)[1]
            .split("BASIS_END", 1)[0]
            .strip()
        )
        record["basis_sha256"] = hashlib.sha256(basis.encode()).hexdigest()
        record["basis_preview"] = basis[:4000]
    basis_entries = [
        entry.strip()
        for entry in basis.replace(",\r\n", ",\n").split(",\n")
        if entry.strip()
    ]
    basis_polynomials = [
        sp.Poly(
            sp.sympify(
                entry.replace("^", "**"),
                locals={str(variable): variable for variable in variables},
            ),
            *variables,
            modulus=prime,
        )
        for entry in basis_entries
    ]
    record["stdout_tail"] = completed.stdout[-2000:]
    record["stderr"] = completed.stderr[-2000:]
    return record, basis_polynomials, basis


def modular_component_minor_generators(
    polynomials: list[sp.Poly],
    linear_factor: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    prime: int,
    boundary_test: dict[str, object],
    selected_indices_override: list[int] | None = None,
    excluded_pairs: set[tuple[int, int]] | None = None,
    singular: str = "Singular",
    timeout: int = 60,
) -> tuple[dict[str, object], list[sp.Poly]]:
    """Construct certified residual-intersection minors for ``(L,H)``.

    For ``P_i=L*A_i+H*B_i``, each minor ``A_i*B_j-A_j*B_i`` is killed into
    the input ideal by both ``L`` and ``H``.  Thus adjoining the minors is a
    certified enlargement inside the desired localization saturation.
    """
    started = time.monotonic()
    pivot = next(
        variable for variable in variables if str(variable) == boundary_test["pivot"]
    )
    selected_indices = (
        list(selected_indices_override)
        if selected_indices_override is not None
        else sorted(
            range(len(polynomials)),
            key=lambda index: len(polynomials[index].terms()),
        )[:3]
    )
    selected_polynomials = [polynomials[index] for index in selected_indices]
    substitution = sp.sympify(
        str(boundary_test["substitution"]),
        locals={str(variable): variable for variable in variables},
    )
    boundary_variables = tuple(variable for variable in variables if variable != pivot)
    boundary_polynomials = [
        modular_linear_specialization(
            polynomial, pivot, substitution, variables, prime
        )
        for polynomial in selected_polynomials
    ]
    common_boundary = sp.Poly(
        sp.sympify(
            str(boundary_test["_boundary_gcd_full"]).replace("^", "**"),
            locals={str(variable): variable for variable in boundary_variables},
        ),
        *boundary_variables,
        modulus=prime,
    ).monic()
    boundary_coefficients = [
        boundary_polynomial.exquo(common_boundary)
        for boundary_polynomial in boundary_polynomials
    ]
    common_lift = sp.Poly(
        common_boundary.as_expr(), *variables, modulus=prime
    )
    lifted_boundary_coefficients = [
        sp.Poly(coefficient.as_expr(), *variables, modulus=prime)
        for coefficient in boundary_coefficients
    ]
    factor_polynomial = sp.Poly(linear_factor, *variables, modulus=prime)
    linear_coefficients: list[sp.Poly] = []
    for polynomial, boundary_coefficient in zip(
        selected_polynomials, lifted_boundary_coefficients, strict=True
    ):
        numerator = polynomial - common_lift * boundary_coefficient
        quotient, remainder = numerator.div(factor_polynomial)
        if not remainder.is_zero:
            raise RuntimeError("component decomposition is not exact")
        linear_coefficients.append(quotient)
    minors: list[sp.Poly] = []
    pairs = [
        (left, right)
        for left in range(len(selected_polynomials))
        for right in range(left + 1, len(selected_polynomials))
        if (left, right) not in (excluded_pairs or set())
    ]
    pairs.sort(
        key=lambda pair: (
            len(linear_coefficients[pair[0]].terms())
            * len(lifted_boundary_coefficients[pair[1]].terms())
            + len(linear_coefficients[pair[1]].terms())
            * len(lifted_boundary_coefficients[pair[0]].terms())
        )
    )
    selected_pair: tuple[int, int] | None = None
    for left, right in pairs:
        multiplication_started = time.monotonic()
        minor = (
            linear_coefficients[left] * lifted_boundary_coefficients[right]
            - linear_coefficients[right] * lifted_boundary_coefficients[left]
        )
        if minor.is_zero:
            continue
        minor = minor.monic()
        minors.append(minor)
        selected_pair = (left, right)
        multiplication_seconds = time.monotonic() - multiplication_started
        break
    unique_minors: dict[str, sp.Poly] = {}
    for minor in minors:
        expression_text = singular_expression(minor, variables, prime)
        unique_minors.setdefault(
            hashlib.sha256(expression_text.encode()).hexdigest(), minor
        )
    minors = list(unique_minors.values())
    record: dict[str, object] = {
        "engine": "two_generator_component_minors",
        "input_generators": len(polynomials),
        "selected_generator_indices": selected_indices,
        "selected_generator_term_counts": [
            len(polynomial.terms()) for polynomial in selected_polynomials
        ],
        "selected_generator_sha256": [
            hashlib.sha256(
                singular_expression(polynomial, variables, prime).encode()
            ).hexdigest()
            for polynomial in selected_polynomials
        ],
        "common_boundary_degree": common_boundary.total_degree(),
        "common_boundary_sha256": boundary_test.get("boundary_gcd_sha256"),
        "linear_coefficient_term_counts": [
            len(coefficient.terms()) for coefficient in linear_coefficients
        ],
        "boundary_coefficient_term_counts": [
            len(coefficient.terms()) for coefficient in boundary_coefficients
        ],
        "minor_count": len(minors),
        "selected_minor_pair": selected_pair,
        "multiplication_seconds": (
            multiplication_seconds if selected_pair is not None else 0.0
        ),
        "minor_term_counts": [len(minor.terms()) for minor in minors],
        "minor_total_degrees": [minor.total_degree() for minor in minors],
        "minor_sha256": list(unique_minors),
        "multiplication_identities_by_construction": 2 * len(minors),
        "status": "completed" if minors else "failed",
        "seconds": time.monotonic() - started,
    }
    return record, minors


def modular_component_minor_boundary_screen(
    polynomials: list[sp.Poly],
    linear_factor: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    prime: int,
    singular: str,
    timeout: int,
    boundary_test: dict[str, object],
) -> dict[str, object]:
    """Screen all residual-intersection minors on a linear boundary.

    The calculation is a cheap diagnostic for a possible transverse
    component, but it is promoted only after comparison with the reciprocal
    solved-hypersurface valuations.  No full minor is expanded in the ambient
    ring.
    """
    started = time.monotonic()
    pivot = next(
        variable for variable in variables if str(variable) == boundary_test["pivot"]
    )
    substitution = sp.sympify(
        str(boundary_test["substitution"]),
        locals={str(variable): variable for variable in variables},
    )
    boundary_variables = tuple(variable for variable in variables if variable != pivot)
    selected_indices = sorted(
        range(len(polynomials)), key=lambda index: len(polynomials[index].terms())
    )
    selected_polynomials = [polynomials[index] for index in selected_indices]
    boundary_polynomials = [
        modular_linear_specialization(
            polynomial, pivot, substitution, variables, prime
        )
        for polynomial in selected_polynomials
    ]
    common_boundary = sp.Poly(
        sp.sympify(
            str(boundary_test["_boundary_gcd_full"]).replace("^", "**"),
            locals={str(variable): variable for variable in boundary_variables},
        ),
        *boundary_variables,
        modulus=prime,
    ).monic()
    boundary_coefficients = [
        boundary_polynomial.exquo(common_boundary)
        for boundary_polynomial in boundary_polynomials
    ]
    factor_in_pivot = sp.Poly(linear_factor, pivot)
    pivot_coefficient = sp.Rational(factor_in_pivot.coeff_monomial(pivot))
    pivot_coefficient_mod = (
        int(pivot_coefficient.p) * pow(int(pivot_coefficient.q), -1, prime)
    ) % prime
    inverse_pivot_coefficient = pow(pivot_coefficient_mod, -1, prime)

    # If P=L*A+H*B and H,B are lifted independently of the pivot, then
    # (dP/dpivot)|_{L=0}=(dL/dpivot)*A|_{L=0}.  This obtains exactly the
    # boundary value needed by the minor test without expanding the often
    # enormous ambient quotient (P-H*B)/L.
    boundary_linear_coefficients = [
        modular_linear_specialization(
            polynomial.diff(pivot), pivot, substitution, variables, prime
        ).mul_ground(inverse_pivot_coefficient)
        for polynomial in selected_polynomials
    ]
    pairs = [
        (left, right)
        for left in range(len(selected_polynomials))
        for right in range(left + 1, len(selected_polynomials))
    ]
    pairs.sort(
        key=lambda pair: (
            len(boundary_linear_coefficients[pair[0]].terms())
            * len(boundary_coefficients[pair[1]].terms())
            + len(boundary_linear_coefficients[pair[1]].terms())
            * len(boundary_coefficients[pair[0]].terms())
        )
    )
    common_lift = sp.Poly(common_boundary.as_expr(), *variables, modulus=prime)
    source_lines = [
        f"ring R={prime},({','.join(map(str, variables))}),dp;",
        "short=0;",
        "poly h="
        + singular_expression(common_lift, variables, prime)
        + ";",
    ]
    substitution_text = singular_expression(substitution, variables, prime)
    for index, polynomial in enumerate(selected_polynomials):
        source_lines.extend(
            [
                f"poly p{index}="
                + singular_expression(polynomial, variables, prime)
                + ";",
                f"poly a{index}=subst(diff(p{index},{pivot}),{pivot},{substitution_text})*{inverse_pivot_coefficient};",
                f"poly b{index}=subst(p{index},{pivot},{substitution_text})/h;",
            ]
        )
    source_lines.extend(
        [
            "timer=1;",
            *[
                f"poly ha{index}=gcd(h,a{index});"
                for index in range(len(boundary_linear_coefficients))
            ],
            *[
                f"poly hb{index}=gcd(h,b{index});"
                for index in range(len(boundary_coefficients))
            ],
            "int all_a_divisible=1;",
            "int all_b_divisible=1;",
            *[
                f"if (deg(ha{index})<deg(h)) {{ all_a_divisible=0; }}"
                for index in range(len(boundary_coefficients))
            ],
            *[
                f"if (deg(hb{index})<deg(h)) {{ all_b_divisible=0; }}"
                for index in range(len(boundary_coefficients))
            ],
            "int forced_common_factor=0;",
            "if ((all_a_divisible==1)||(all_b_divisible==1)) { forced_common_factor=1; }",
            "int found=0;",
            "int pairs_examined=0;",
            "poly m=0;",
            "poly g=0;",
        ]
    )
    for left, right in pairs:
        source_lines.extend(
            [
                "if ((found==0)&&(forced_common_factor==0)) {",
                f"  m=a{left}*b{right}-a{right}*b{left};",
                "  g=gcd(h,m);",
                "  pairs_examined=pairs_examined+1;",
                f'  print("PAIR={left},{right},"+string(deg(g))+","+string(size(m)));',
                "  if (deg(g)==0) { found=1; }",
                "}",
            ]
        )
    source_lines.extend(
        [
            "int elapsed=timer;",
            'print("RESULT_BEGIN");',
            'print("elapsed_ticks="+string(elapsed));',
            'print("pairs_examined="+string(pairs_examined));',
            'print("forced_common_factor="+string(forced_common_factor));',
            *[
                f'print("A_GCD_{index}="+string(deg(ha{index})));'
                for index in range(len(boundary_linear_coefficients))
            ],
            *[
                f'print("B_GCD_{index}="+string(deg(hb{index})));'
                for index in range(len(boundary_coefficients))
            ],
            'print("RESULT_END"); quit;',
        ]
    )
    source = "\n".join(source_lines) + "\n"
    record: dict[str, object] = {
        "engine": "component_minor_boundary_screen",
        "selected_generator_indices": selected_indices,
        "selected_generator_term_counts": [
            len(polynomial.terms()) for polynomial in selected_polynomials
        ],
        "selected_generator_sha256": [
            hashlib.sha256(
                singular_expression(polynomial, variables, prime).encode()
            ).hexdigest()
            for polynomial in selected_polynomials
        ],
        "common_boundary_degree": common_boundary.total_degree(),
        "common_boundary_sha256": boundary_test.get("boundary_gcd_sha256"),
        "boundary_linear_coefficient_term_counts": [
            len(coefficient.terms()) for coefficient in boundary_linear_coefficients
        ],
        "boundary_coefficient_term_counts": [
            len(coefficient.terms()) for coefficient in boundary_coefficients
        ],
        "input_sha256": hashlib.sha256(source.encode()).hexdigest(),
    }
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
    markers = "RESULT_BEGIN" in completed.stdout and "RESULT_END" in completed.stdout
    record["status"] = (
        "completed" if completed.returncode == 0 and markers else "failed"
    )
    pair_records: list[dict[str, object]] = []
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        if stripped.startswith("PAIR="):
            left, right, gcd_degree, terms = map(
                int, stripped.split("=", 1)[1].split(",")
            )
            pair_records.append(
                {
                    "local_pair": [left, right],
                    "generator_pair": [
                        selected_indices[left],
                        selected_indices[right],
                    ],
                    "boundary_gcd_degree": gcd_degree,
                    "boundary_minor_terms": terms,
                }
            )
        else:
            for key in ("elapsed_ticks", "pairs_examined", "forced_common_factor"):
                if stripped.startswith(key + "="):
                    record[key] = int(stripped.split("=", 1)[1])
            for key, output_key in (
                ("A_GCD_", "boundary_linear_common_gcd_degrees"),
                ("B_GCD_", "boundary_residual_common_gcd_degrees"),
            ):
                if stripped.startswith(key):
                    index_text, value = stripped.removeprefix(key).split("=", 1)
                    values = record.setdefault(
                        output_key, [0] * len(selected_polynomials)
                    )
                    values[int(index_text)] = int(value)
    pair_records.sort(
        key=lambda item: (
            int(item["boundary_gcd_degree"]),
            int(item["boundary_minor_terms"]),
        )
    )
    record["pairs"] = pair_records
    if pair_records:
        record["best_pair"] = pair_records[0]
    record["coprime_boundary_minor_found"] = int(
        bool(pair_records) and pair_records[0]["boundary_gcd_degree"] == 0
    )
    record["stdout_tail"] = completed.stdout[-2000:]
    record["stderr"] = completed.stderr[-2000:]
    return record


def cleared_linear_hypersurface_expression(
    polynomial: sp.Poly,
    pivot: sp.Symbol,
    leading_coefficient: sp.Poly,
    remainder: sp.Poly,
    variables: tuple[sp.Symbol, ...],
    prime: int,
) -> str:
    """Serialize a polynomial restricted to C*pivot+R=0.

    If the pivot degree is n, the returned expression is
    C^n P(-R/C).  It is polynomial and differs from the restriction only by
    a power of C, which is harmless after C has been localized.
    """
    boundary_variables = tuple(variable for variable in variables if variable != pivot)
    pivot_index = variables.index(pivot)
    degree = polynomial.degree(pivot)
    coefficient_dicts: dict[int, dict[tuple[int, ...], int]] = {}
    for monomial, coefficient in polynomial.terms():
        pivot_power = monomial[pivot_index]
        boundary_monomial = monomial[:pivot_index] + monomial[pivot_index + 1 :]
        coefficient_dicts.setdefault(pivot_power, {})[boundary_monomial] = (
            int(coefficient) % prime
        )
    c_text = singular_expression(leading_coefficient, boundary_variables, prime)
    r_text = singular_expression(remainder, boundary_variables, prime)
    summands: list[str] = []
    for pivot_power in sorted(coefficient_dicts):
        coefficient_polynomial = sp.Poly.from_dict(
            coefficient_dicts[pivot_power],
            boundary_variables,
            modulus=prime,
        )
        coefficient_text = singular_expression(
            coefficient_polynomial, boundary_variables, prime
        )
        factors = [f"({coefficient_text})"]
        if pivot_power:
            factors.append(f"(-({r_text}))^{pivot_power}")
        c_power = degree - pivot_power
        if c_power:
            factors.append(f"({c_text})^{c_power}")
        summands.append("*".join(factors))
    return "+".join(summands) if summands else "0"


def full_linear_hypersurface_coordinate_expression(
    polynomial: sp.Poly,
    pivot: sp.Symbol,
    leading_coefficient: sp.Poly,
    remainder: sp.Poly,
    variables: tuple[sp.Symbol, ...],
    normal_coordinate: sp.Symbol,
    prime: int,
) -> str:
    """Serialize ``C^n P((d-R)/C)`` without expanding it in Python.

    Here the hypersurface equation is ``D=C*pivot+R`` and ``d`` is a new
    coordinate replacing ``D`` on the chart ``C != 0``.  Specializing
    ``d=0`` recovers :func:`cleared_linear_hypersurface_expression`.
    """
    boundary_variables = tuple(variable for variable in variables if variable != pivot)
    pivot_index = variables.index(pivot)
    degree = polynomial.degree(pivot)
    coefficient_dicts: dict[int, dict[tuple[int, ...], int]] = {}
    for monomial, coefficient in polynomial.terms():
        pivot_power = monomial[pivot_index]
        boundary_monomial = monomial[:pivot_index] + monomial[pivot_index + 1 :]
        coefficient_dicts.setdefault(pivot_power, {})[boundary_monomial] = (
            int(coefficient) % prime
        )
    c_text = singular_expression(leading_coefficient, boundary_variables, prime)
    r_text = singular_expression(remainder, boundary_variables, prime)
    summands: list[str] = []
    for pivot_power in sorted(coefficient_dicts):
        coefficient_polynomial = sp.Poly.from_dict(
            coefficient_dicts[pivot_power],
            boundary_variables,
            modulus=prime,
        )
        coefficient_text = singular_expression(
            coefficient_polynomial, boundary_variables, prime
        )
        factors = [f"({coefficient_text})"]
        if pivot_power:
            factors.append(
                f"({normal_coordinate}-({r_text}))^{pivot_power}"
            )
        c_power = degree - pivot_power
        if c_power:
            factors.append(f"({c_text})^{c_power}")
        summands.append("*".join(factors))
    return "+".join(summands) if summands else "0"


def modular_solved_hypersurface_coordinate_colon(
    polynomials: list[sp.Poly],
    hypersurface_factor: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    prime: int,
    singular: str,
    timeout: int,
    boundary_only: bool = False,
    component_quotient: bool = False,
    component_minor: bool = False,
    minor_boundary_only: bool = False,
) -> dict[str, object]:
    """Discover the colon by a solved hypersurface in its own coordinate.

    Direct saturation by the degree-17 factor is expensive.  If
    ``D=C*pivot+R``, the chart ``C != 0`` admits the polynomial coordinate
    ``d=D`` after clearing the predictable power of ``C`` in each input.
    Saturation is then by the single variable ``d``.  This remains a modular
    discovery calculation until the resulting basis identities are promoted
    over ``QQ``.
    """
    started = time.monotonic()
    pivot_candidates: list[
        tuple[tuple[int, int, str], sp.Symbol, sp.Expr, sp.Expr]
    ] = []
    for candidate in variables:
        factor_in_candidate = sp.Poly(hypersurface_factor, candidate)
        if factor_in_candidate.degree() != 1:
            continue
        coefficient = factor_in_candidate.coeff_monomial(candidate)
        remainder_expression = factor_in_candidate.coeff_monomial(1)
        boundary_variables = tuple(
            variable for variable in variables if variable != candidate
        )
        coefficient_polynomial = rational_polynomial_mod(
            coefficient, boundary_variables, prime
        )
        pivot_candidates.append(
            (
                (
                    len(coefficient_polynomial.terms()),
                    coefficient_polynomial.total_degree(),
                    str(candidate),
                ),
                candidate,
                coefficient,
                remainder_expression,
            )
        )
    if not pivot_candidates:
        return {
            "engine": "solved_hypersurface_coordinate_colon",
            "status": "failed",
            "reason": "hypersurface has no linear pivot",
        }
    _, pivot, coefficient_expression, remainder_expression = min(pivot_candidates)
    boundary_variables = tuple(variable for variable in variables if variable != pivot)
    leading_coefficient = rational_polynomial_mod(
        coefficient_expression, boundary_variables, prime
    )
    remainder = rational_polynomial_mod(
        remainder_expression, boundary_variables, prime
    )
    normal_coordinate = sp.Symbol("delta")
    if normal_coordinate in variables:
        normal_coordinate = sp.Symbol("delta_D")
    coordinate_variables = (normal_coordinate,) + boundary_variables
    transformed_expressions = [
        full_linear_hypersurface_coordinate_expression(
            polynomial,
            pivot,
            leading_coefficient,
            remainder,
            variables,
            normal_coordinate,
            prime,
        )
        for polynomial in polynomials
    ]
    declarations = [
        f"poly q{index}={expression};"
        for index, expression in enumerate(transformed_expressions)
    ]
    boundary_source = "\n".join(
        [
            f"ring R={prime},({','.join(map(str, coordinate_variables))}),dp;",
            "short=0;",
            *declarations,
            f"poly h=subst(q0,{normal_coordinate},0);",
            *[
                f"h=gcd(h,subst(q{index},{normal_coordinate},0));"
                for index in range(1, len(polynomials))
            ],
            "int boundary_gcd_degree=deg(h);",
            'print("RESULT_BEGIN");',
            'print("boundary_gcd_degree="+string(boundary_gcd_degree));',
            'print("BOUNDARY_GCD_BEGIN"); print(h); print("BOUNDARY_GCD_END");',
            'print("RESULT_END"); quit;',
        ]
    ) + "\n"
    source_lines = [
        'LIB "elim.lib";',
        "option(redSB);",
        f"ring R={prime},({','.join(map(str, coordinate_variables))}),dp;",
        "short=0;",
        *declarations,
        "timer=1;",
    ]
    if component_minor:
        if len(polynomials) < 2:
            return {
                "engine": "solved_hypersurface_coordinate_colon",
                "status": "failed",
                "reason": "component minor requires at least two generators",
            }
        source_lines.extend(
            [
                f"poly component_h=subst(q0,{normal_coordinate},0);",
                f"component_h=gcd(component_h,subst(q1,{normal_coordinate},0));",
                f"poly b0=subst(q0,{normal_coordinate},0)/component_h;",
                f"poly b1=subst(q1,{normal_coordinate},0)/component_h;",
                f"poly u0=(q0-component_h*b0)/{normal_coordinate};",
                f"poly u1=(q1-component_h*b1)/{normal_coordinate};",
                "poly component_minor=u0*b1-u1*b0;",
                "int component_minor_terms=size(component_minor);",
                "int component_minor_degree=deg(component_minor);",
                "ideal I=ideal(q0,q1,component_minor)"
                + (
                    "+ideal("
                    + ",".join(
                        f"q{index}" for index in range(2, len(polynomials))
                    )
                    + ")"
                    if len(polynomials) > 2
                    else ""
                )
                + ";",
            ]
        )
    elif component_quotient:
        if len(polynomials) < 2:
            return {
                "engine": "solved_hypersurface_coordinate_colon",
                "status": "failed",
                "reason": "component quotient requires at least two generators",
            }
        source_lines.extend(
            [
                f"poly component_h=subst(q0,{normal_coordinate},0);",
                f"component_h=gcd(component_h,subst(q1,{normal_coordinate},0));",
                "ideal K=std(quotient(ideal(q0,q1),"
                f"ideal({normal_coordinate},component_h)));",
                "int component_basis_size=size(K);",
                "int component_dimension=dim(K);",
                "ideal I=K"
                + (
                    "+ideal("
                    + ",".join(
                        f"q{index}" for index in range(2, len(polynomials))
                    )
                    + ")"
                    if len(polynomials) > 2
                    else ""
                )
                + ";",
            ]
        )
    else:
        source_lines.append(
            "ideal I="
            + ",".join(f"q{index}" for index in range(len(polynomials)))
            + ";"
        )
    if component_minor and minor_boundary_only:
        source_lines.extend(
            [
                f"poly post_h=subst(q0,{normal_coordinate},0);",
                f"post_h=gcd(post_h,subst(q1,{normal_coordinate},0));",
                *[
                    f"post_h=gcd(post_h,subst(q{index},{normal_coordinate},0));"
                    for index in range(2, len(polynomials))
                ],
                f"post_h=gcd(post_h,subst(component_minor,{normal_coordinate},0));",
                "int post_minor_boundary_gcd_degree=deg(post_h);",
                "int elapsed=timer;",
                'print("RESULT_BEGIN");',
                'print("component_gcd_degree="+string(deg(component_h)));',
                'print("component_minor_terms="+string(component_minor_terms));',
                'print("component_minor_degree="+string(component_minor_degree));',
                'print("post_minor_boundary_gcd_degree="+string(post_minor_boundary_gcd_degree));',
                'print("POST_MINOR_BOUNDARY_GCD_BEGIN"); print(post_h); print("POST_MINOR_BOUNDARY_GCD_END");',
                'print("elapsed_ticks="+string(elapsed));',
                'print("RESULT_END"); quit;',
            ]
        )
    else:
        source_lines.extend(
            [
                f"ideal J={normal_coordinate};",
                "list L=sat_with_exp(I,J);",
                "ideal G=std(L[1]);",
                "int elapsed=timer;",
                "int unit=0;",
                "if ((size(G)==1)&&(G[1]==1)) { unit=1; }",
                'print("RESULT_BEGIN");',
                'print("unit="+string(unit));',
                'print("basis_size="+string(size(G)));',
                'print("dimension="+string(dim(G)));',
                'print("saturation_exponent="+string(L[2]));',
            *(
                [
                    'print("component_gcd_degree="+string(deg(component_h)));',
                    'print("component_basis_size="+string(component_basis_size));',
                    'print("component_dimension="+string(component_dimension));',
                ]
                if component_quotient
                else []
            ),
            *(
                [
                    'print("component_gcd_degree="+string(deg(component_h)));',
                    'print("component_minor_terms="+string(component_minor_terms));',
                    'print("component_minor_degree="+string(component_minor_degree));',
                ]
                if component_minor
                else []
            ),
                'print("elapsed_ticks="+string(elapsed));',
                'print("BASIS_BEGIN"); print(G); print("BASIS_END");',
                'print("RESULT_END"); quit;',
            ]
        )
    source = "\n".join(source_lines) + "\n"
    record: dict[str, object] = {
        "engine": "solved_hypersurface_coordinate_colon",
        "component_quotient": int(component_quotient),
        "component_minor": int(component_minor),
        "minor_boundary_only": int(minor_boundary_only),
        "pivot": str(pivot),
        "normal_coordinate": str(normal_coordinate),
        "coordinate_variables": [str(variable) for variable in coordinate_variables],
        "leading_coefficient_degree": leading_coefficient.total_degree(),
        "leading_coefficient_terms": len(leading_coefficient.terms()),
        "remainder_degree": remainder.total_degree(),
        "remainder_terms": len(remainder.terms()),
        "input_generator_count": len(polynomials),
        "input_generator_term_counts": [
            len(polynomial.terms()) for polynomial in polynomials
        ],
        "input_generator_pivot_degrees": [
            polynomial.degree(pivot) for polynomial in polynomials
        ],
        "transformed_compact_expression_lengths": [
            len(expression) for expression in transformed_expressions
        ],
        "boundary_input_sha256": hashlib.sha256(
            boundary_source.encode()
        ).hexdigest(),
        "input_sha256": hashlib.sha256(source.encode()).hexdigest(),
    }
    try:
        boundary_completed = subprocess.run(
            [singular, "-q"],
            input=boundary_source,
            text=True,
            capture_output=True,
            timeout=min(timeout, 120),
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        record.update(
            {
                "status": "timeout",
                "stage": "boundary_gcd",
                "seconds": time.monotonic() - started,
                "stdout_tail": (error.stdout or "")[-2000:],
                "stderr": (error.stderr or "")[-2000:],
            }
        )
        return record
    boundary_markers = (
        "RESULT_BEGIN" in boundary_completed.stdout
        and "RESULT_END" in boundary_completed.stdout
    )
    if boundary_completed.returncode != 0 or not boundary_markers:
        record.update(
            {
                "status": "failed",
                "stage": "boundary_gcd",
                "seconds": time.monotonic() - started,
                "stdout_tail": boundary_completed.stdout[-2000:],
                "stderr": boundary_completed.stderr[-2000:],
            }
        )
        return record
    for line in boundary_completed.stdout.splitlines():
        stripped = line.strip()
        if stripped.startswith("boundary_gcd_degree="):
            record["boundary_gcd_degree"] = int(stripped.split("=", 1)[1])
    if (
        "BOUNDARY_GCD_BEGIN" in boundary_completed.stdout
        and "BOUNDARY_GCD_END" in boundary_completed.stdout
    ):
        boundary_gcd = (
            boundary_completed.stdout.split("BOUNDARY_GCD_BEGIN", 1)[1]
            .split("BOUNDARY_GCD_END", 1)[0]
            .strip()
        )
        record["boundary_gcd_sha256"] = hashlib.sha256(
            boundary_gcd.encode()
        ).hexdigest()
        record["boundary_gcd_preview"] = boundary_gcd[:4000]
    record["divisorial_boundary_component_absent"] = int(
        record.get("boundary_gcd_degree") == 0
    )
    if boundary_only:
        record.update(
            {
                "status": "completed",
                "colon_skipped": 1,
                "seconds": time.monotonic() - started,
                "stdout_tail": boundary_completed.stdout[-2000:],
                "stderr": boundary_completed.stderr[-2000:],
            }
        )
        return record
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
    markers = "RESULT_BEGIN" in completed.stdout and "RESULT_END" in completed.stdout
    record["status"] = (
        "completed" if completed.returncode == 0 and markers else "failed"
    )
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        for key in (
            "unit",
            "basis_size",
            "dimension",
            "saturation_exponent",
            "elapsed_ticks",
            "component_gcd_degree",
            "component_basis_size",
            "component_dimension",
            "component_minor_terms",
            "component_minor_degree",
            "post_minor_boundary_gcd_degree",
        ):
            if stripped.startswith(key + "="):
                record[key] = int(stripped.split("=", 1)[1])
    if "BASIS_BEGIN" in completed.stdout and "BASIS_END" in completed.stdout:
        basis = (
            completed.stdout.split("BASIS_BEGIN", 1)[1]
            .split("BASIS_END", 1)[0]
            .strip()
        )
        record["basis_sha256"] = hashlib.sha256(basis.encode()).hexdigest()
        record["basis_preview"] = basis[:4000]
    if (
        "POST_MINOR_BOUNDARY_GCD_BEGIN" in completed.stdout
        and "POST_MINOR_BOUNDARY_GCD_END" in completed.stdout
    ):
        post_minor_boundary_gcd = (
            completed.stdout.split("POST_MINOR_BOUNDARY_GCD_BEGIN", 1)[1]
            .split("POST_MINOR_BOUNDARY_GCD_END", 1)[0]
            .strip()
        )
        record["post_minor_boundary_gcd_sha256"] = hashlib.sha256(
            post_minor_boundary_gcd.encode()
        ).hexdigest()
        record["post_minor_boundary_gcd_preview"] = (
            post_minor_boundary_gcd[:4000]
        )
    record["stdout_tail"] = completed.stdout[-2000:]
    record["stderr"] = completed.stderr[-2000:]
    return record


def modular_solved_hypersurface_gcd_test(
    polynomials: list[sp.Poly],
    hypersurface_factor: sp.Expr,
    cycle_supporting_factors: list[dict[str, object]],
    already_localized_factors: list[dict[str, object]],
    variables: tuple[sp.Symbol, ...],
    prime: int,
    singular: str,
    timeout: int,
) -> dict[str, object]:
    """Classify height-one divisors after solving a hypersurface linearly.

    The input ideal is known to be unmixed of codimension two.  A gcd of its
    restrictions to H=0 therefore records exactly the codimension-two
    components contained in H, on the chart where the chosen leading
    coefficient is invertible.
    """
    started = time.monotonic()
    pivot_candidates: list[tuple[tuple[int, int, str], sp.Symbol, sp.Expr, sp.Expr]] = []
    for candidate in variables:
        factor_in_candidate = sp.Poly(hypersurface_factor, candidate)
        if factor_in_candidate.degree() != 1:
            continue
        coefficient = factor_in_candidate.coeff_monomial(candidate)
        remainder_expression = factor_in_candidate.coeff_monomial(1)
        boundary_variables = tuple(
            variable for variable in variables if variable != candidate
        )
        coefficient_polynomial = rational_polynomial_mod(
            coefficient, boundary_variables, prime
        )
        pivot_candidates.append(
            (
                (
                    len(coefficient_polynomial.terms()),
                    coefficient_polynomial.total_degree(),
                    str(candidate),
                ),
                candidate,
                coefficient,
                remainder_expression,
            )
        )
    if not pivot_candidates:
        return {
            "engine": "solved_hypersurface_boundary_gcd",
            "status": "failed",
            "reason": "hypersurface has no linear pivot",
        }
    _, pivot, coefficient_expression, remainder_expression = min(pivot_candidates)
    boundary_variables = tuple(variable for variable in variables if variable != pivot)
    leading_coefficient = rational_polynomial_mod(
        coefficient_expression, boundary_variables, prime
    )
    remainder = rational_polynomial_mod(
        remainder_expression, boundary_variables, prime
    )
    restricted_expressions = [
        cleared_linear_hypersurface_expression(
            polynomial,
            pivot,
            leading_coefficient,
            remainder,
            variables,
            prime,
        )
        for polynomial in polynomials
    ]

    def restricted_factor_record(
        factor: dict[str, object], prefix: str, index: int
    ) -> dict[str, object] | None:
        factor_polynomial = rational_polynomial_mod(
            factor["expression"], variables, prime
        )
        expression = cleared_linear_hypersurface_expression(
            factor_polynomial,
            pivot,
            leading_coefficient,
            remainder,
            variables,
            prime,
        )
        restricted = rational_polynomial_mod(
            sp.sympify(
                expression.replace("^", "**"),
                locals={str(variable): variable for variable in boundary_variables},
            ),
            boundary_variables,
            prime,
        )
        if restricted.is_zero or restricted.total_degree() == 0:
            return None
        restricted = restricted.monic()
        restricted_text = singular_expression(
            restricted, boundary_variables, prime
        )
        return {
            "name": f"{prefix}{index}",
            "factor": factor["expression_text"],
            "factor_sha256": factor["sha256"],
            "expression": restricted_text,
            "degree": restricted.total_degree(),
            "terms": len(restricted.terms()),
            "sha256": hashlib.sha256(restricted_text.encode()).hexdigest(),
        }

    localized_restrictions = [
        record
        for index, factor in enumerate(already_localized_factors)
        if (
            record := restricted_factor_record(factor, "u", index)
        )
        is not None
    ]
    support_restrictions = [
        record
        for index, factor in enumerate(cycle_supporting_factors)
        if (
            record := restricted_factor_record(factor, "s", index)
        )
        is not None
    ]
    source_lines = [
        f"ring R={prime},({','.join(map(str, boundary_variables))}),dp;",
        "short=0;",
        "timer=1;",
    ]
    for index, expression in enumerate(restricted_expressions):
        source_lines.append(f"poly q{index}={expression};")
    source_lines.extend(
        [
            "poly h=q0;",
            *[
                f"h=gcd(h,q{index});"
                for index in range(1, len(restricted_expressions))
            ],
            "poly raw_h=h;",
        ]
    )
    for record in localized_restrictions + support_restrictions:
        name = str(record["name"])
        source_lines.extend(
            [
                f"poly {name}={record['expression']};",
                f"ideal I{name}=std(ideal({name}));",
                f"int e{name}=0;",
                f"while ((h!=0)&&(reduce(h,I{name})==0)) {{",
                f"  h=h/{name}; e{name}=e{name}+1;",
                "}",
            ]
        )
    valuation_names: list[tuple[str, int, str]] = []
    for support_record in support_restrictions:
        support_name = str(support_record["name"])
        for polynomial_index in range(len(restricted_expressions)):
            valuation_name = f"v{support_name}_{polynomial_index}"
            valuation_names.append(
                (valuation_name, polynomial_index, str(support_record["factor"]))
            )
            source_lines.extend(
                [
                    f"poly {valuation_name}=q{polynomial_index};",
                    f"int e{valuation_name}=0;",
                    f"while (({valuation_name}!=0)&&(reduce({valuation_name},I{support_name})==0)) {{",
                    f"  {valuation_name}={valuation_name}/{support_name};",
                    f"  e{valuation_name}=e{valuation_name}+1;",
                    "}",
                ]
            )
    source_lines.extend(
        [
            "int elapsed=timer;",
            'print("RESULT_BEGIN");',
            'print("raw_gcd_degree="+string(deg(raw_h)));',
            'print("residual_gcd_degree="+string(deg(h)));',
            'print("elapsed_ticks="+string(elapsed));',
            *[
                f'print("EXPONENT_{record["name"]}="+string(e{record["name"]}));'
                for record in localized_restrictions + support_restrictions
            ],
            *[
                f'print("VALUATION_{name}="+string(e{name}));'
                for name, _, _ in valuation_names
            ],
            'print("RAW_GCD_BEGIN"); print(raw_h); print("RAW_GCD_END");',
            'print("RESIDUAL_GCD_BEGIN"); print(h); print("RESIDUAL_GCD_END");',
            'print("RESULT_END"); quit;',
        ]
    )
    source = "\n".join(source_lines) + "\n"
    record: dict[str, object] = {
        "engine": "solved_hypersurface_boundary_gcd",
        "pivot": str(pivot),
        "leading_coefficient_degree": leading_coefficient.total_degree(),
        "leading_coefficient_terms": len(leading_coefficient.terms()),
        "leading_coefficient_sha256": hashlib.sha256(
            singular_expression(
                leading_coefficient, boundary_variables, prime
            ).encode()
        ).hexdigest(),
        "remainder_degree": remainder.total_degree(),
        "remainder_terms": len(remainder.terms()),
        "restricted_input_count": len(restricted_expressions),
        "input_generator_term_counts": [
            len(polynomial.terms()) for polynomial in polynomials
        ],
        "input_generator_sha256": [
            hashlib.sha256(
                singular_expression(polynomial, variables, prime).encode()
            ).hexdigest()
            for polynomial in polynomials
        ],
        "localized_restrictions": localized_restrictions,
        "cycle_support_restrictions": support_restrictions,
        "input_sha256": hashlib.sha256(source.encode()).hexdigest(),
    }
    identity_source_lines = [
        f"ring V={prime},({','.join(map(str, variables))}),dp;",
        "short=0;",
        "poly d="
        + singular_expression(hypersurface_factor, variables, prime)
        + ";",
        "poly c="
        + singular_expression(leading_coefficient, boundary_variables, prime)
        + ";",
        "ideal GD=std(ideal(d));",
    ]
    for index, (polynomial, expression) in enumerate(
        zip(polynomials, restricted_expressions, strict=True)
    ):
        identity_source_lines.extend(
            [
                f"poly p{index}="
                + singular_expression(polynomial, variables, prime)
                + ";",
                f"poly qq{index}={expression};",
                f"poly rr{index}=reduce(c^{polynomial.degree(pivot)}*p{index}-qq{index},GD);",
            ]
        )
    conormal_names: list[tuple[str, int, str]] = []
    for support_index, support_factor in enumerate(cycle_supporting_factors):
        support_variables = [
            variable
            for variable in variables
            if variable in support_factor["expression"].free_symbols
        ]
        if not support_variables:
            continue
        support_pivot = support_variables[-1]
        support_in_pivot = sp.Poly(support_factor["expression"], support_pivot)
        if (
            support_in_pivot.degree() != 1
            or support_in_pivot.coeff_monomial(support_pivot).free_symbols
        ):
            continue
        support_name = f"l{support_index}"
        identity_source_lines.extend(
            [
                f"poly {support_name}="
                + singular_expression(
                    support_factor["expression"], variables, prime
                )
                + ";",
                f"ideal K{support_name}=std(ideal(d,{support_name}));",
            ]
        )
        for polynomial_index in range(len(polynomials)):
            conormal_name = f"n{support_index}_{polynomial_index}"
            conormal_names.append(
                (
                    conormal_name,
                    polynomial_index,
                    str(support_factor["expression_text"]),
                )
            )
            identity_source_lines.append(
                f"poly {conormal_name}=reduce(diff(p{polynomial_index},{support_pivot}),K{support_name});"
            )
    identity_source_lines.extend(
        [
            'print("RESULT_BEGIN");',
            *[
                f'print("REMAINDER_{index}="+string(size(rr{index})));'
                for index in range(len(polynomials))
            ],
            *[
                f'print("CONORMAL_{name}="+string(size({name})));'
                for name, _, _ in conormal_names
            ],
            'print("RESULT_END"); quit;',
        ]
    )
    identity_source = "\n".join(identity_source_lines) + "\n"
    record["restriction_identity_input_sha256"] = hashlib.sha256(
        identity_source.encode()
    ).hexdigest()
    try:
        identity_completed = subprocess.run(
            [singular, "-q"],
            input=identity_source,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        record.update(
            {
                "status": "timeout",
                "stage": "restriction_identity",
                "seconds": time.monotonic() - started,
                "stdout_tail": (error.stdout or "")[-2000:],
                "stderr": (error.stderr or "")[-2000:],
            }
        )
        return record
    identity_remainders: list[int] = [1] * len(polynomials)
    conormal_remainder_terms: dict[str, list[int]] = {
        str(factor["expression_text"]): [0] * len(polynomials)
        for factor in cycle_supporting_factors
    }
    for line in identity_completed.stdout.splitlines():
        stripped = line.strip()
        if stripped.startswith("REMAINDER_"):
            index_text, value = stripped.removeprefix("REMAINDER_").split("=", 1)
            identity_remainders[int(index_text)] = int(value)
        if stripped.startswith("CONORMAL_"):
            name, value = stripped.split("=", 1)
            conormal_name = name.removeprefix("CONORMAL_")
            _, polynomial_index, factor_name = next(
                item for item in conormal_names if item[0] == conormal_name
            )
            conormal_remainder_terms[factor_name][polynomial_index] = int(value)
    record["restriction_identity_remainder_terms"] = identity_remainders
    record["conormal_remainder_terms"] = conormal_remainder_terms
    record["restriction_identities_verified"] = int(
        identity_completed.returncode == 0
        and "RESULT_BEGIN" in identity_completed.stdout
        and "RESULT_END" in identity_completed.stdout
        and all(value == 0 for value in identity_remainders)
    )
    if record["restriction_identities_verified"] != 1:
        record.update(
            {
                "status": "failed",
                "stage": "restriction_identity",
                "seconds": time.monotonic() - started,
                "stdout_tail": identity_completed.stdout[-2000:],
                "stderr": identity_completed.stderr[-2000:],
            }
        )
        return record
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
    markers = "RESULT_BEGIN" in completed.stdout and "RESULT_END" in completed.stdout
    record["status"] = (
        "completed" if completed.returncode == 0 and markers else "failed"
    )
    exponents: dict[str, int] = {}
    individual_valuations: dict[str, list[int]] = {
        str(item["factor"]): [0] * len(restricted_expressions)
        for item in support_restrictions
    }
    restriction_lookup = {
        str(item["name"]): item
        for item in localized_restrictions + support_restrictions
    }
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        for key in ("raw_gcd_degree", "residual_gcd_degree", "elapsed_ticks"):
            if stripped.startswith(key + "="):
                record[key] = int(stripped.split("=", 1)[1])
        if stripped.startswith("EXPONENT_"):
            name, value = stripped.split("=", 1)
            restriction = restriction_lookup[name.removeprefix("EXPONENT_")]
            exponents[str(restriction["factor"])] = int(value)
        if stripped.startswith("VALUATION_"):
            name, value = stripped.split("=", 1)
            valuation_name = name.removeprefix("VALUATION_")
            _, polynomial_index, factor_name = next(
                item for item in valuation_names if item[0] == valuation_name
            )
            individual_valuations[factor_name][polynomial_index] = int(value)
    record["factor_exponents"] = exponents
    record["individual_support_valuations"] = individual_valuations
    for label in ("RAW_GCD", "RESIDUAL_GCD"):
        begin = label + "_BEGIN"
        end = label + "_END"
        if begin in completed.stdout and end in completed.stdout:
            polynomial_text = (
                completed.stdout.split(begin, 1)[1].split(end, 1)[0].strip()
            )
            key = label.lower()
            record[key + "_sha256"] = hashlib.sha256(
                polynomial_text.encode()
            ).hexdigest()
            record[key + "_preview"] = polynomial_text[:4000]
    support_names = [str(item["factor"]) for item in support_restrictions]
    record["support_complete"] = int(
        record.get("status") == "completed"
        and record.get("restriction_identities_verified") == 1
        and record.get("residual_gcd_degree") == 0
        and bool(support_names)
        and all(exponents.get(name, 0) > 0 for name in support_names)
    )
    record["stdout_tail"] = completed.stdout[-2000:]
    record["stderr"] = completed.stderr[-2000:]
    return record


def modular_two_normal_local_primary_test(
    polynomials: list[sp.Poly],
    hypersurface_factor: sp.Expr,
    hypersurface_pivot: sp.Symbol,
    collision_factor: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    prime: int,
    singular: str,
    timeout: int,
    parameter_values_override: dict[sp.Symbol, int] | None = None,
) -> dict[str, object]:
    """Compute a deterministic fiber in the two normal coordinates.

    By default the non-normal variables are specialized to the small integer
    tuple ``(2, 3, 5)``; a complete override may select another fiber.  A
    local ``ds`` standard basis then measures the primary length on that
    fiber of the intersection of the two declared divisors.
    """
    started = time.monotonic()
    collision_support = [
        variable for variable in variables if variable in collision_factor.free_symbols
    ]
    if not collision_support:
        return {"status": "failed", "reason": "constant collision factor"}
    collision_pivot = collision_support[-1]
    if collision_pivot == hypersurface_pivot:
        return {"status": "failed", "reason": "normal pivots coincide"}
    parameters = tuple(
        variable
        for variable in variables
        if variable not in {hypersurface_pivot, collision_pivot}
    )
    hypersurface_in_pivot = sp.Poly(hypersurface_factor, hypersurface_pivot)
    collision_in_pivot = sp.Poly(collision_factor, collision_pivot)
    if hypersurface_in_pivot.degree() != 1 or collision_in_pivot.degree() != 1:
        return {"status": "failed", "reason": "both normals must be linear"}
    hypersurface_coefficient = rational_polynomial_mod(
        hypersurface_in_pivot.coeff_monomial(hypersurface_pivot),
        parameters,
        prime,
    )
    hypersurface_remainder = rational_polynomial_mod(
        hypersurface_in_pivot.coeff_monomial(1), parameters, prime
    )
    collision_coefficient = sp.Rational(
        collision_in_pivot.coeff_monomial(collision_pivot)
    )
    if collision_coefficient.free_symbols:
        return {"status": "failed", "reason": "collision pivot coefficient varies"}
    collision_remainder = rational_polynomial_mod(
        collision_in_pivot.coeff_monomial(1), parameters, prime
    )
    collision_coefficient_mod = (
        int(collision_coefficient.p)
        * pow(int(collision_coefficient.q), -1, prime)
    ) % prime
    inverse_collision_coefficient = pow(collision_coefficient_mod, -1, prime)
    if parameter_values_override is None:
        parameter_values = {
            parameter: value % prime
            for parameter, value in zip(parameters, (2, 3, 5))
        }
    else:
        if set(parameter_values_override) != set(parameters):
            return {
                "status": "failed",
                "reason": "parameter override must specify every non-normal variable",
            }
        parameter_values = {
            parameter: int(parameter_values_override[parameter]) % prime
            for parameter in parameters
        }

    def evaluate_parameter_polynomial(polynomial: sp.Poly) -> int:
        answer = 0
        for monomial, coefficient in polynomial.terms():
            term = int(coefficient) % prime
            for parameter, power in zip(parameters, monomial, strict=True):
                term = term * pow(parameter_values[parameter], power, prime) % prime
            answer = (answer + term) % prime
        return answer

    c_value = evaluate_parameter_polynomial(hypersurface_coefficient)
    r_value = evaluate_parameter_polynomial(hypersurface_remainder)
    collision_remainder_value = evaluate_parameter_polynomial(
        collision_remainder
    )
    if c_value == 0:
        return {
            "status": "failed",
            "reason": "chosen parameter fiber kills the hypersurface pivot",
        }
    hypersurface_index = variables.index(hypersurface_pivot)
    collision_index = variables.index(collision_pivot)
    parameter_indices = [variables.index(parameter) for parameter in parameters]

    dd_symbol, ee_symbol = sp.symbols("dd ee")

    def local_polynomial(polynomial: sp.Poly) -> sp.Poly:
        hypersurface_degree = polynomial.degree(hypersurface_pivot)
        coefficients: dict[tuple[int, int], int] = {}
        for monomial, coefficient in polynomial.terms():
            modular_coefficient = int(coefficient) % prime
            parameter_coefficient = modular_coefficient
            for parameter, parameter_index in zip(
                parameters, parameter_indices, strict=True
            ):
                power = monomial[parameter_index]
                parameter_coefficient = (
                    parameter_coefficient
                    * pow(parameter_values[parameter], power, prime)
                    % prime
                )
            hypersurface_power = monomial[hypersurface_index]
            clearing_power = hypersurface_degree - hypersurface_power
            collision_power = monomial[collision_index]
            common_coefficient = (
                parameter_coefficient
                * pow(c_value, clearing_power, prime)
                * pow(inverse_collision_coefficient, collision_power, prime)
            ) % prime
            for dd_power in range(hypersurface_power + 1):
                dd_coefficient = (
                    math.comb(hypersurface_power, dd_power)
                    * pow(
                        -r_value % prime,
                        hypersurface_power - dd_power,
                        prime,
                    )
                ) % prime
                for ee_power in range(collision_power + 1):
                    ee_coefficient = (
                        math.comb(collision_power, ee_power)
                        * pow(
                            -collision_remainder_value % prime,
                            collision_power - ee_power,
                            prime,
                        )
                    ) % prime
                    local_monomial = (dd_power, ee_power)
                    coefficients[local_monomial] = (
                        coefficients.get(local_monomial, 0)
                        + common_coefficient * dd_coefficient * ee_coefficient
                    ) % prime
        return sp.Poly.from_dict(
            {
                monomial: coefficient
                for monomial, coefficient in coefficients.items()
                if coefficient
            },
            dd_symbol,
            ee_symbol,
            modulus=prime,
        )

    unique_polynomials: list[sp.Poly] = []
    input_hashes: list[str] = []
    seen_hashes: set[str] = set()
    for polynomial in polynomials:
        polynomial_hash = hashlib.sha256(
            singular_expression(polynomial, variables, prime).encode()
        ).hexdigest()
        if polynomial_hash in seen_hashes:
            continue
        seen_hashes.add(polynomial_hash)
        unique_polynomials.append(polynomial)
        input_hashes.append(polynomial_hash)
    selected = sorted(
        zip(unique_polynomials, input_hashes, strict=True),
        key=lambda item: len(item[0].terms()),
    )[:2]
    unique_polynomials = [item[0] for item in selected]
    input_hashes = [item[1] for item in selected]
    hypersurface_origin = (-r_value * pow(c_value, -1, prime)) % prime
    collision_origin = (
        -collision_remainder_value
        * inverse_collision_coefficient
        % prime
    )
    origin_values = {
        **parameter_values,
        hypersurface_pivot: hypersurface_origin,
        collision_pivot: collision_origin,
    }

    def evaluate_input_at_origin(polynomial: sp.Poly) -> int:
        answer = 0
        for monomial, coefficient in polynomial.terms():
            term = int(coefficient) % prime
            for variable, power in zip(variables, monomial, strict=True):
                term = term * pow(origin_values[variable], power, prime) % prime
            answer = (answer + term) % prime
        return answer

    direct_origin_values = [
        evaluate_input_at_origin(polynomial) for polynomial in unique_polynomials
    ]
    local_polynomials = [
        local_polynomial(polynomial) for polynomial in unique_polynomials
    ]
    expressions = [
        singular_expression(
            polynomial, (dd_symbol, ee_symbol), prime
        )
        for polynomial in local_polynomials
    ]
    summand_lists = [expression.split("+") for expression in expressions]
    polynomial_source_lines: list[str] = []
    for index, summands in enumerate(summand_lists):
        polynomial_source_lines.append(f"poly f{index}=0;")
        for chunk_start in range(0, len(summands), 100):
            chunk = "+".join(summands[chunk_start : chunk_start + 100])
            polynomial_source_lines.append(f"f{index}=f{index}+{chunk};")
    source = "\n".join(
        [
            f"ring L={prime},(dd,ee),ds;",
            "short=0;",
            *polynomial_source_lines,
            *[
                f"poly c{index}=subst(subst(f{index},dd,0),ee,0);"
                for index in range(len(expressions))
            ],
            "timer=1;",
            "ideal G=std(ideal("
            + ",".join(f"f{index}" for index in range(len(expressions)))
            + "));",
            "int elapsed=timer;",
            "int local_dimension=dim(G);",
            "int local_length=-1;",
            "if (local_dimension==0) { local_length=vdim(G); }",
            'print("RESULT_BEGIN");',
            'print("basis_size="+string(size(G)));',
            'print("local_dimension="+string(local_dimension));',
            'print("local_length="+string(local_length));',
            'print("elapsed_ticks="+string(elapsed));',
            *[
                f'print("CONSTANT_{index}="+string(size(c{index})));'
                for index in range(len(expressions))
            ],
            'print("BASIS_BEGIN"); print(G); print("BASIS_END");',
            'print("RESULT_END"); quit;',
        ]
    ) + "\n"
    record: dict[str, object] = {
        "engine": "two_normal_specialized_local_primary",
        "hypersurface_pivot": str(hypersurface_pivot),
        "collision_pivot": str(collision_pivot),
        "parameters": [str(parameter) for parameter in parameters],
        "parameter_values": {
            str(parameter): value for parameter, value in parameter_values.items()
        },
        "hypersurface_pivot_origin": hypersurface_origin,
        "collision_pivot_origin": collision_origin,
        "direct_origin_values": direct_origin_values,
        "input_generators": len(polynomials),
        "selected_input_generators": len(unique_polynomials),
        "selection_rule": "two sparsest distinct generators",
        "input_generator_sha256": input_hashes,
        "coordinate_expression_sha256": [
            hashlib.sha256(expression.encode()).hexdigest()
            for expression in expressions
        ],
        "input_sha256": hashlib.sha256(source.encode()).hexdigest(),
        "input_preview": source[:4000],
    }
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
    markers = "RESULT_BEGIN" in completed.stdout and "RESULT_END" in completed.stdout
    record["status"] = (
        "completed" if completed.returncode == 0 and markers else "failed"
    )
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        for key in ("basis_size", "local_dimension", "local_length", "elapsed_ticks"):
            if stripped.startswith(key + "="):
                record[key] = int(stripped.split("=", 1)[1])
        if stripped.startswith("CONSTANT_"):
            index_text, value = stripped.removeprefix("CONSTANT_").split("=", 1)
            constants = record.setdefault(
                "coordinate_origin_constant_terms", [1] * len(expressions)
            )
            constants[int(index_text)] = int(value)
    if "BASIS_BEGIN" in completed.stdout and "BASIS_END" in completed.stdout:
        basis = (
            completed.stdout.split("BASIS_BEGIN", 1)[1]
            .split("BASIS_END", 1)[0]
            .strip()
        )
        record["basis_sha256"] = hashlib.sha256(basis.encode()).hexdigest()
        record["basis_preview"] = basis[:4000]
    record["curvilinear_length_twelve"] = int(
        record.get("status") == "completed"
        and record.get("local_dimension") == 0
        and record.get("local_length") == 12
    )
    record["stdout_tail"] = completed.stdout[-2000:]
    record["stderr"] = completed.stderr[-2000:]
    return record


def modular_component_minor_sequence(
    polynomials: list[sp.Poly],
    linear_factor: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    prime: int,
    singular: str,
    timeout: int,
    initial_boundary_test: dict[str, object],
    max_steps: int = 8,
) -> tuple[dict[str, object], list[sp.Poly]]:
    """Adjoin sparse component minors until the linear boundary gcd is one."""
    current_polynomials = list(polynomials)
    current_boundary_test = initial_boundary_test
    layers: list[dict[str, object]] = []
    completed = 0
    candidate_order = sorted(
        range(len(polynomials)), key=lambda index: len(polynomials[index].terms())
    )[:3]
    used_pairs: set[tuple[int, int]] = set()
    for index in range(1, max_steps + 1):
        minor_step, minors = modular_component_minor_generators(
            current_polynomials,
            linear_factor,
            variables,
            prime,
            current_boundary_test,
            candidate_order,
            used_pairs,
            singular,
            max(timeout, 60),
        )
        layer: dict[str, object] = {
            "index": index,
            "boundary_gcd_degree_before": current_boundary_test.get(
                "boundary_gcd_degree"
            ),
            "minor_step": minor_step,
        }
        layers.append(layer)
        selected_pair = minor_step.get("selected_minor_pair")
        if isinstance(selected_pair, tuple):
            used_pairs.add(selected_pair)
        if minor_step.get("status") != "completed" or not minors:
            break
        current_polynomials.extend(minors)
        current_boundary_test = modular_linear_boundary_gcd_test(
            current_polynomials,
            linear_factor,
            variables,
            prime,
            singular,
            timeout,
        )
        layer["post_boundary_test"] = {
            key: value
            for key, value in current_boundary_test.items()
            if not key.startswith("_")
        }
        if current_boundary_test.get("saturated") == 1:
            completed = 1
            break
        if not (
            current_boundary_test.get("status") == "completed"
            and int(current_boundary_test.get("boundary_gcd_degree", 0)) > 0
        ):
            break
    return {
        "engine": "iterated_two_generator_component_minors",
        "status": "completed" if completed else "partial",
        "steps": len(layers),
        "generators_added": len(current_polynomials) - len(polynomials),
        "seconds": sum(
            float(layer["minor_step"].get("seconds", 0.0))
            + float(layer.get("post_boundary_test", {}).get("seconds", 0.0))
            for layer in layers
        ),
        "layers": layers,
    }, current_polynomials


def modular_successive_saturation_cutoff(
    equations: dict[int, sp.Expr],
    saturation: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    prime: int,
    singular: str,
    timeout: int,
    report_quotient_dimension: bool = False,
) -> dict[str, object]:
    """Saturate by irreducible factors one at a time with checkpoints.

    Each factor gets its own Singular process and timeout.  The complete
    reduced basis from one step is fed into the next step, while the artifact
    retains only hashes and previews.  Multiplicity in the localization
    product is recorded but one quotient saturation per irreducible factor is
    sufficient.
    """
    factors = saturation_factor_records(saturation, variables)
    factor_count = len(factors)
    processing_factors = list(factors)
    deferred_factor_hashes: set[str] = set()
    component_dependency_records: dict[str, list[dict[str, object]]] = {}
    component_cycles: list[dict[str, object]] = []
    cutoff = max(equations)
    ideal_text = ",".join(
        singular_expression(equation, variables, prime)
        for equation in equations.values()
    )
    factor_order = [
        {
            key: value
            for key, value in factor.items()
            if key not in {"expression", "sort_key"}
        }
        for factor in factors
    ]
    steps: list[dict[str, object]] = []
    final_basis = ideal_text
    final_dimension = None
    final_basis_size = None
    unit = 0
    total_seconds = 0.0
    status = "completed"
    principal_polynomial = (
        sp.Poly(next(iter(equations.values())), *variables, modulus=prime)
        if len(equations) == 1
        else None
    )
    ideal_polynomials = [
        sp.Poly(equation, *variables, modulus=prime)
        for equation in equations.values()
    ]
    for factor_index, factor in enumerate(processing_factors, start=1):
        if principal_polynomial is not None:
            step_started = time.monotonic()
            factor_polynomial = sp.Poly(
                factor["expression"], *variables, modulus=prime
            )
            if int(factor["total_degree"]) == 1:
                removed_exponent = 0
                quotient, remainder = principal_polynomial.div(factor_polynomial)
                while remainder.is_zero:
                    principal_polynomial = quotient
                    removed_exponent += 1
                    quotient, remainder = principal_polynomial.div(factor_polynomial)
                principal_polynomial = principal_polynomial.monic()
                basis = singular_expression(
                    principal_polynomial, variables, prime
                )
                step = {
                    "engine": "principal_linear_division",
                    "status": "completed",
                    "seconds": time.monotonic() - step_started,
                    "unit": int(principal_polynomial.total_degree() == 0),
                    "saturation_exponent": removed_exponent,
                    "basis_sha256": hashlib.sha256(basis.encode()).hexdigest(),
                    "basis_preview": basis[:4000],
                }
            else:
                step, basis = modular_principal_gcd_step(
                    principal_polynomial,
                    factor["expression"],
                    variables,
                    prime,
                    singular,
                    timeout,
                )
                if step.get("status") == "completed" and basis:
                    parsed_basis = sp.sympify(
                        basis.replace("^", "**"),
                        locals={str(variable): variable for variable in variables},
                    )
                    principal_polynomial = sp.Poly(
                        parsed_basis, *variables, modulus=prime
                    ).monic()
            if step.get("status") != "completed" or not basis:
                step.update(
                    {
                        "factor_index": factor_index,
                        "factor_sha256": factor["sha256"],
                        "factor": factor["expression_text"],
                    }
                )
                total_seconds += float(step.get("seconds", 0.0))
                steps.append(step)
                status = str(step.get("status", "failed"))
                print(
                    f"SAT_FACTOR_{status.upper()} p={prime} index={factor_index}",
                    flush=True,
                )
                break
            unit = int(principal_polynomial.total_degree() == 0)
            final_dimension = -1 if unit else len(variables) - 1
            final_basis_size = 1
            final_basis = basis
            ideal_text = basis
            step.update(
                {
                "factor_index": factor_index,
                "factor_sha256": factor["sha256"],
                "factor": factor["expression_text"],
                "unit": unit,
                "basis_size": 1,
                "dimension": final_dimension,
                "quotient_dimension": -1,
                }
            )
            total_seconds += float(step["seconds"])
            steps.append(step)
            print(
                f"SAT_FACTOR_PRINCIPAL p={prime} index={factor_index}/{len(factors)} "
                f"dim={final_dimension} exp={step.get('saturation_exponent')} "
                f"factor={factor['expression_text']}",
                flush=True,
            )
            if unit == 1:
                break
            continue
        predivision_exponents: list[int] = []
        if int(factor["total_degree"]) == 1:
            factor_polynomial = sp.Poly(
                factor["expression"], *variables, modulus=prime
            )
            reduced_polynomials: list[sp.Poly] = []
            for polynomial in ideal_polynomials:
                exponent = 0
                quotient, remainder = polynomial.div(factor_polynomial)
                while remainder.is_zero:
                    polynomial = quotient
                    exponent += 1
                    quotient, remainder = polynomial.div(factor_polynomial)
                reduced_polynomials.append(polynomial)
                predivision_exponents.append(exponent)
            ideal_polynomials = reduced_polynomials
            ideal_text = ",".join(
                singular_expression(polynomial, variables, prime)
                for polynomial in ideal_polynomials
            )
        hypersurface_test = None
        if int(factor["total_degree"]) > 1:
            cycle_dependencies = component_dependency_records.get(
                str(factor["sha256"]), []
            )
            if cycle_dependencies:
                cycle_factor_hashes = {
                    str(dependency["linear_factor_sha256"])
                    for dependency in cycle_dependencies
                }
                cycle_supporting_factors = [
                    supporting_factor
                    for supporting_factor in factors
                    if str(supporting_factor["sha256"]) in cycle_factor_hashes
                ]
                already_localized_factors = [
                    localized_factor
                    for localized_factor in factors
                    if str(localized_factor["sha256"])
                    not in cycle_factor_hashes | {str(factor["sha256"])}
                ]
                cycle_test = modular_solved_hypersurface_gcd_test(
                    ideal_polynomials,
                    factor["expression"],
                    cycle_supporting_factors,
                    already_localized_factors,
                    variables,
                    prime,
                    singular,
                    timeout,
                )
                support_closed = int(
                    cycle_test.get("support_complete") == 1
                    and len(cycle_supporting_factors)
                    == len(cycle_dependencies)
                )
                conormal_inconsistencies: list[dict[str, object]] = []
                scheme_multiplicity_obstructions: list[dict[str, object]] = []
                individual_valuations = cycle_test.get(
                    "individual_support_valuations", {}
                )
                conormal_remainders = cycle_test.get(
                    "conormal_remainder_terms", {}
                )
                for dependency in cycle_dependencies:
                    linear_name = str(dependency["linear_factor"])
                    valuations = individual_valuations.get(linear_name, [])
                    remainders = conormal_remainders.get(linear_name, [])
                    if valuations and min(valuations) > 1:
                        scheme_multiplicity_obstructions.append(
                            {
                                "linear_factor": linear_name,
                                "minimum_hypersurface_valuation": min(valuations),
                                "individual_hypersurface_valuations": valuations,
                            }
                        )
                    if (
                        dependency.get("coprime_boundary_minor_found") == 1
                        and valuations
                        and min(valuations) > 1
                        and remainders
                        and max(remainders) == 0
                    ):
                        conormal_inconsistencies.append(
                            {
                                "linear_factor": linear_name,
                                "minimum_hypersurface_valuation": min(valuations),
                                "direct_conormal_remainder_terms": remainders,
                                "reason": (
                                    "coprime boundary-minor diagnostic conflicts "
                                    "with the verified solved-hypersurface valuation"
                                ),
                            }
                        )
                local_primary_test: dict[str, object] | None = None
                thick_factor_names = {
                    str(obstruction["linear_factor"])
                    for obstruction in scheme_multiplicity_obstructions
                }
                if len(thick_factor_names) == 1:
                    thick_factor = next(
                        supporting_factor
                        for supporting_factor in cycle_supporting_factors
                        if str(supporting_factor["expression_text"])
                        in thick_factor_names
                    )
                    hypersurface_pivot = next(
                        variable
                        for variable in variables
                        if str(variable) == cycle_test.get("pivot")
                    )
                    local_primary_test = modular_two_normal_local_primary_test(
                        ideal_polynomials,
                        factor["expression"],
                        hypersurface_pivot,
                        thick_factor["expression"],
                        variables,
                        prime,
                        singular,
                        timeout,
                    )
                cycle_record = {
                    "hypersurface_factor": factor["expression_text"],
                    "hypersurface_factor_sha256": factor["sha256"],
                    "linear_dependencies": cycle_dependencies,
                    "hypersurface_test": cycle_test,
                    "support_closed": support_closed,
                    "scheme_closed": int(
                        support_closed == 1
                        and not conormal_inconsistencies
                        and not scheme_multiplicity_obstructions
                        and all(
                            dependency.get("coprime_boundary_minor_found") == 1
                            for dependency in cycle_dependencies
                        )
                    ),
                    "conormal_inconsistencies": conormal_inconsistencies,
                    "scheme_multiplicity_obstructions": (
                        scheme_multiplicity_obstructions
                    ),
                    "local_primary_test": local_primary_test,
                    "primary_profile_closed": int(
                        local_primary_test is not None
                        and local_primary_test.get(
                            "curvilinear_length_twelve"
                        )
                        == 1
                    ),
                }
                component_cycles.append(cycle_record)
                if cycle_record["support_closed"]:
                    basis = ideal_text
                    step = {
                        "factor_index": factor_index,
                        "factor_sha256": factor["sha256"],
                        "factor": factor["expression_text"],
                        "engine": "unmixed_localization_component_support_cycle",
                        "status": "classified",
                        "seconds": float(cycle_test.get("seconds", 0.0))
                        + float(
                            (local_primary_test or {}).get("seconds", 0.0)
                        ),
                        "unit": None,
                        "basis_size": len(ideal_polynomials),
                        "dimension": len(variables) - 2,
                        "quotient_dimension": -1,
                        "basis_kind": "pre_saturation_generators",
                        "basis_sha256": hashlib.sha256(basis.encode()).hexdigest(),
                        "basis_preview": basis[:4000],
                        "component_cycle": cycle_record,
                    }
                    steps.append(step)
                    total_seconds += float(step["seconds"])
                    final_basis = basis
                    final_dimension = len(variables) - 2
                    final_basis_size = len(ideal_polynomials)
                    status = "component_support_classified"
                    print(
                        f"SAT_FACTOR_COMPONENT_SUPPORT_CYCLE p={prime} "
                        f"index={factor_index}/{len(processing_factors)} "
                        f"support={len(cycle_supporting_factors)} "
                        f"factor={factor['expression_text']}",
                        flush=True,
                    )
                    break
            hypersurface_test = modular_unmixed_hypersurface_dimension_test(
                ideal_polynomials,
                factor["expression"],
                variables,
                prime,
                singular,
                timeout,
            )
            if hypersurface_test.get("avoids_components") == 1:
                basis = ideal_text
                step = {
                    "factor_index": factor_index,
                    "factor_sha256": factor["sha256"],
                    "factor": factor["expression_text"],
                    "engine": "unmixed_hypersurface_dimension",
                    "status": "completed",
                    "seconds": hypersurface_test.get("seconds", 0.0),
                    "unit": 0,
                    "basis_size": len(ideal_polynomials),
                    "dimension": len(variables) - 2,
                    "quotient_dimension": -1,
                    "saturation_exponent": 0,
                    "basis_sha256": hashlib.sha256(basis.encode()).hexdigest(),
                    "basis_preview": basis[:4000],
                    "hypersurface_test": hypersurface_test,
                }
                steps.append(step)
                total_seconds += float(step["seconds"])
                final_basis = basis
                final_dimension = len(variables) - 2
                final_basis_size = len(ideal_polynomials)
                print(
                    f"SAT_FACTOR_HYPERSURFACE_DIMENSION p={prime} "
                    f"index={factor_index}/{len(processing_factors)} "
                    f"dim={final_dimension} factor={factor['expression_text']}",
                    flush=True,
                )
                continue
        linear_boundary_test = None
        if int(factor["total_degree"]) == 1 and len(ideal_polynomials) > 2:
            linear_boundary_test = modular_linear_boundary_gcd_test(
                ideal_polynomials,
                factor["expression"],
                variables,
                prime,
                singular,
                timeout,
            )
            public_boundary_test = {
                key: value
                for key, value in linear_boundary_test.items()
                if not key.startswith("_")
            }
            if linear_boundary_test.get("saturated") == 1:
                basis = ideal_text
                step = {
                    "factor_index": factor_index,
                    "factor_sha256": factor["sha256"],
                    "factor": factor["expression_text"],
                    "engine": "unmixed_linear_boundary_gcd",
                    "status": "completed",
                    "seconds": linear_boundary_test.get("seconds", 0.0),
                    "unit": 0,
                    "basis_size": len(ideal_polynomials),
                    "dimension": len(variables) - 2,
                    "quotient_dimension": -1,
                    "saturation_exponent": 0,
                    "basis_sha256": hashlib.sha256(basis.encode()).hexdigest(),
                    "basis_preview": basis[:4000],
                    "boundary_test": public_boundary_test,
                }
                if any(predivision_exponents):
                    step["generator_predivision_exponents"] = predivision_exponents
                steps.append(step)
                total_seconds += float(step["seconds"])
                final_basis = basis
                final_dimension = len(variables) - 2
                final_basis_size = len(ideal_polynomials)
                print(
                    f"SAT_FACTOR_UNMIXED_BOUNDARY p={prime} "
                    f"index={factor_index}/{len(processing_factors)} "
                    f"dim={final_dimension} factor={factor['expression_text']}",
                    flush=True,
                )
                continue
            deferred_support, deferred_support_exponents = (
                boundary_gcd_supported_on_factors(
                    linear_boundary_test,
                    [
                        supporting_factor
                        for supporting_factor in factors
                        if supporting_factor["sha256"] != factor["sha256"]
                    ],
                    variables,
                    prime,
                )
            )
            if deferred_support and len(deferred_support_exponents) == 1:
                if int(linear_boundary_test.get("boundary_gcd_degree", 0)) <= 4:
                    component_minor_boundary_screen = (
                        modular_component_minor_boundary_screen(
                            ideal_polynomials,
                            factor["expression"],
                            variables,
                            prime,
                            singular,
                            timeout,
                            linear_boundary_test,
                        )
                    )
                else:
                    component_minor_boundary_screen = {
                        "engine": "component_minor_boundary_screen",
                        "status": "skipped",
                        "reason": (
                            "high-degree reciprocal solved-hypersurface "
                            "valuation is the authoritative certificate"
                        ),
                        "common_boundary_degree": linear_boundary_test.get(
                            "boundary_gcd_degree"
                        ),
                        "coprime_boundary_minor_found": 0,
                    }
                public_boundary_test["component_support_exponents"] = (
                    deferred_support_exponents
                )
                public_boundary_test["component_minor_boundary_screen"] = (
                    component_minor_boundary_screen
                )
                supporting_expression = next(
                    iter(deferred_support_exponents)
                )
                supporting_factor = next(
                    candidate
                    for candidate in factors
                    if candidate["expression_text"] == supporting_expression
                )
                component_dependency_records.setdefault(
                    str(supporting_factor["sha256"]), []
                ).append(
                    {
                        "linear_factor": factor["expression_text"],
                        "linear_factor_sha256": factor["sha256"],
                        "supporting_factor": supporting_expression,
                        "supporting_factor_sha256": supporting_factor["sha256"],
                        "support_exponent": deferred_support_exponents[
                            supporting_expression
                        ],
                        "coprime_boundary_minor_found": component_minor_boundary_screen.get(
                            "coprime_boundary_minor_found", 0
                        ),
                        "boundary_minor_pair": component_minor_boundary_screen.get(
                            "best_pair"
                        ),
                        "screen_input_sha256": component_minor_boundary_screen.get(
                            "input_sha256"
                        ),
                    }
                )
                print(
                    f"SAT_FACTOR_COMPONENT_SCREEN p={prime} "
                    f"index={factor_index}/{len(processing_factors)} "
                    f"status={component_minor_boundary_screen.get('status')} "
                    f"pairs={component_minor_boundary_screen.get('pairs_examined', 0)} "
                    f"coprime_minor={component_minor_boundary_screen.get('coprime_boundary_minor_found', 0)}",
                    flush=True,
                )
                if (
                    component_minor_boundary_screen.get("status") == "completed"
                    and int(linear_boundary_test.get("boundary_gcd_degree", 0)) <= 4
                    and not component_minor_boundary_screen.get(
                        "coprime_boundary_minor_found", 0
                    )
                ):
                    component_minor_sequence, component_polynomials = (
                        modular_component_minor_sequence(
                            ideal_polynomials,
                            factor["expression"],
                            variables,
                            prime,
                            singular,
                            timeout,
                            linear_boundary_test,
                        )
                    )
                    public_boundary_test["component_minor_sequence"] = (
                        component_minor_sequence
                    )
                    if component_minor_sequence.get("generators_added", 0):
                        ideal_polynomials = component_polynomials
                        ideal_text = ",".join(
                            singular_expression(polynomial, variables, prime)
                            for polynomial in ideal_polynomials
                        )
                        if component_minor_sequence.get("status") == "completed":
                            basis = ideal_text
                            step = {
                                "factor_index": factor_index,
                                "factor_sha256": factor["sha256"],
                                "factor": factor["expression_text"],
                                "engine": "linear_supported_component_minors",
                                "basis_kind": "colon_containment_generators",
                                "status": "completed",
                                "seconds": (
                                    float(
                                        linear_boundary_test.get("seconds", 0.0)
                                    )
                                    + float(
                                        component_minor_sequence.get("seconds", 0.0)
                                    )
                                ),
                                "unit": 0,
                                "basis_size": len(ideal_polynomials),
                                "dimension": len(variables) - 2,
                                "quotient_dimension": -1,
                                "saturation_exponent": 1,
                                "basis_sha256": hashlib.sha256(
                                    basis.encode()
                                ).hexdigest(),
                                "basis_preview": basis[:4000],
                                "boundary_test": public_boundary_test,
                            }
                            steps.append(step)
                            total_seconds += float(step["seconds"])
                            final_basis = basis
                            final_dimension = len(variables) - 2
                            final_basis_size = len(ideal_polynomials)
                            print(
                                f"SAT_FACTOR_COMPONENT_MINORS p={prime} "
                                f"index={factor_index}/{len(processing_factors)} "
                                f"dim={final_dimension} "
                                f"factor={factor['expression_text']}",
                                flush=True,
                            )
                            continue
            if (
                deferred_support
                and factor["sha256"] not in deferred_factor_hashes
            ):
                deferred_factor_hashes.add(str(factor["sha256"]))
                processing_factors.append(factor)
                basis = ideal_text
                step = {
                    "factor_index": factor_index,
                    "factor_sha256": factor["sha256"],
                    "factor": factor["expression_text"],
                    "engine": "boundary_gcd_supported_on_localization_factors",
                    "status": "deferred",
                    "seconds": linear_boundary_test.get("seconds", 0.0),
                    "unit": 0,
                    "basis_size": len(ideal_polynomials),
                    "dimension": len(variables) - 2,
                    "quotient_dimension": -1,
                    "saturation_exponent_lower_bound": 0,
                    "deferred_support_exponents": deferred_support_exponents,
                    "basis_sha256": hashlib.sha256(basis.encode()).hexdigest(),
                    "basis_preview": basis[:4000],
                    "boundary_test": public_boundary_test,
                }
                steps.append(step)
                total_seconds += float(step["seconds"])
                final_basis = basis
                final_dimension = len(variables) - 2
                final_basis_size = len(ideal_polynomials)
                print(
                    f"SAT_FACTOR_DEPENDENCY_DEFERRED p={prime} "
                    f"index={factor_index}/{len(processing_factors)} "
                    f"support={deferred_support_exponents} "
                    f"factor={factor['expression_text']}",
                    flush=True,
                )
                continue
            generic_colon_exponent = 0
            generic_colon_layers: list[dict[str, object]] = []
            current_boundary_test = linear_boundary_test
            generic_linear_complete = 0
            while (
                current_boundary_test.get("status") == "completed"
                and int(current_boundary_test.get("boundary_gcd_degree", 0)) > 0
                and generic_colon_exponent < 8
            ):
                if len(ideal_polynomials) == 3:
                    colon_step, new_generators = (
                        modular_structured_triple_linear_colon_step(
                            ideal_polynomials,
                            factor["expression"],
                            variables,
                            prime,
                            singular,
                            timeout,
                            current_boundary_test,
                        )
                    )
                else:
                    colon_step, new_generators = modular_linear_colon_syzygy_step(
                        ideal_polynomials,
                        factor["expression"],
                        variables,
                        prime,
                        singular,
                        timeout,
                        current_boundary_test,
                    )
                layer: dict[str, object] = {
                    "colon_exponent_before": generic_colon_exponent,
                    "boundary_test": {
                        key: value
                        for key, value in current_boundary_test.items()
                        if not key.startswith("_")
                    },
                    "colon_step": colon_step,
                }
                generic_colon_layers.append(layer)
                if colon_step.get("status") != "completed" or not new_generators:
                    break
                ideal_polynomials.extend(new_generators)
                ideal_text = ",".join(
                    singular_expression(polynomial, variables, prime)
                    for polynomial in ideal_polynomials
                )
                generic_colon_exponent += 1
                current_boundary_test = modular_linear_boundary_gcd_test(
                    ideal_polynomials,
                    factor["expression"],
                    variables,
                    prime,
                    singular,
                    timeout,
                )
                layer["post_boundary_test"] = {
                    key: value
                    for key, value in current_boundary_test.items()
                    if not key.startswith("_")
                }
                if current_boundary_test.get("saturated") == 1:
                    generic_linear_complete = 1
                    break
            if generic_colon_layers:
                public_boundary_test["colon_layers"] = generic_colon_layers
            generic_linear_seconds = float(
                linear_boundary_test.get("seconds", 0.0)
            ) + sum(
                float(layer["colon_step"].get("seconds", 0.0))
                + float(layer.get("post_boundary_test", {}).get("seconds", 0.0))
                for layer in generic_colon_layers
            )
            if generic_linear_complete:
                basis = ideal_text
                step = {
                    "factor_index": factor_index,
                    "factor_sha256": factor["sha256"],
                    "factor": factor["expression_text"],
                    "engine": "unmixed_linear_boundary_syzygy_colon",
                    "basis_kind": "colon_generators",
                    "status": "completed",
                    "seconds": generic_linear_seconds,
                    "unit": 0,
                    "basis_size": len(ideal_polynomials),
                    "dimension": len(variables) - 2,
                    "quotient_dimension": -1,
                    "saturation_exponent": generic_colon_exponent,
                    "basis_sha256": hashlib.sha256(basis.encode()).hexdigest(),
                    "basis_preview": basis[:4000],
                    "boundary_test": public_boundary_test,
                }
                steps.append(step)
                total_seconds += float(step["seconds"])
                final_basis = basis
                final_dimension = len(variables) - 2
                final_basis_size = len(ideal_polynomials)
                print(
                    f"SAT_FACTOR_UNMIXED_SYZYGY p={prime} "
                    f"index={factor_index}/{len(processing_factors)} "
                    f"dim={final_dimension} exp={generic_colon_exponent} "
                    f"factor={factor['expression_text']}",
                    flush=True,
                )
                continue
            if (
                generic_colon_exponent > 0
                and factor["sha256"] not in deferred_factor_hashes
            ):
                deferred_factor_hashes.add(str(factor["sha256"]))
                processing_factors.append(factor)
                basis = ideal_text
                step = {
                    "factor_index": factor_index,
                    "factor_sha256": factor["sha256"],
                    "factor": factor["expression_text"],
                    "engine": "unmixed_linear_partial_colon_deferred",
                    "basis_kind": "colon_generators",
                    "status": "deferred",
                    "seconds": generic_linear_seconds,
                    "unit": 0,
                    "basis_size": len(ideal_polynomials),
                    "dimension": len(variables) - 2,
                    "quotient_dimension": -1,
                    "saturation_exponent_lower_bound": generic_colon_exponent,
                    "basis_sha256": hashlib.sha256(basis.encode()).hexdigest(),
                    "basis_preview": basis[:4000],
                    "boundary_test": public_boundary_test,
                }
                steps.append(step)
                total_seconds += float(step["seconds"])
                final_basis = basis
                final_dimension = len(variables) - 2
                final_basis_size = len(ideal_polynomials)
                print(
                    f"SAT_FACTOR_UNMIXED_DEFERRED p={prime} "
                    f"index={factor_index}/{len(processing_factors)} "
                    f"exp>={generic_colon_exponent} "
                    f"factor={factor['expression_text']}",
                    flush=True,
                )
                continue
        if int(factor["total_degree"]) == 1 and len(ideal_polynomials) == 2:
            linear_boundary_test = modular_pair_linear_saturation_test(
                ideal_polynomials,
                factor["expression"],
                variables,
                prime,
                singular,
                timeout,
            )
            public_boundary_test = {
                key: value
                for key, value in linear_boundary_test.items()
                if not key.startswith("_")
            }
            if linear_boundary_test.get("saturated") == 1:
                basis = ideal_text
                step = {
                    "factor_index": factor_index,
                    "factor_sha256": factor["sha256"],
                    "factor": factor["expression_text"],
                    "engine": "two_generator_linear_boundary_gcd",
                    "status": "completed",
                    "seconds": linear_boundary_test.get("seconds", 0.0),
                    "unit": 0,
                    "basis_size": len(ideal_polynomials),
                    "dimension": len(variables) - 2,
                    "quotient_dimension": -1,
                    "saturation_exponent": 0,
                    "basis_sha256": hashlib.sha256(basis.encode()).hexdigest(),
                    "basis_preview": basis[:4000],
                    "boundary_test": public_boundary_test,
                }
                if any(predivision_exponents):
                    step["generator_predivision_exponents"] = predivision_exponents
                steps.append(step)
                total_seconds += float(step["seconds"])
                final_basis = basis
                final_dimension = len(variables) - 2
                final_basis_size = len(ideal_polynomials)
                print(
                    f"SAT_FACTOR_BOUNDARY_GCD p={prime} "
                    f"index={factor_index}/{len(factors)} "
                    f"dim={final_dimension} factor={factor['expression_text']}",
                    flush=True,
                )
                continue
            if (
                linear_boundary_test.get("status") == "completed"
                and linear_boundary_test.get("full_gcd_degree") == 0
                and int(linear_boundary_test.get("boundary_gcd_degree", 0)) > 0
            ):
                colon_generator = two_generator_linear_colon_generator(
                    ideal_polynomials,
                    factor["expression"],
                    variables,
                    prime,
                    linear_boundary_test,
                )
                ideal_polynomials.append(colon_generator)
                ideal_text = ",".join(
                    singular_expression(polynomial, variables, prime)
                    for polynomial in ideal_polynomials
                )
                public_boundary_test["colon_generator_terms"] = len(
                    colon_generator.terms()
                )
                public_boundary_test["colon_generator_total_degree"] = (
                    colon_generator.total_degree()
                )
                public_boundary_test["colon_generator_sha256"] = hashlib.sha256(
                    singular_expression(
                        colon_generator, variables, prime
                    ).encode()
                ).hexdigest()
                public_boundary_test["precomputed_colon_exponent"] = 1
                residual_boundary_test = modular_residual_boundary_gcd_test(
                    colon_generator,
                    variables,
                    prime,
                    singular,
                    timeout,
                    linear_boundary_test,
                )
                public_boundary_test["residual_boundary_test"] = (
                    residual_boundary_test
                )
                colon_exponent = 1
                layered_saturation_complete = int(
                    residual_boundary_test.get("saturated_after_first_colon") == 1
                )
                colon_layers: list[dict[str, object]] = []
                component_used_pairs: set[tuple[int, int]] = set()
                while not layered_saturation_complete and colon_exponent < 8:
                    next_boundary_test = modular_linear_boundary_gcd_test(
                        ideal_polynomials,
                        factor["expression"],
                        variables,
                        prime,
                        singular,
                        timeout,
                    )
                    public_next_boundary = {
                        key: value
                        for key, value in next_boundary_test.items()
                        if not key.startswith("_")
                    }
                    layer: dict[str, object] = {
                        "colon_exponent_before": colon_exponent,
                        "boundary_test": public_next_boundary,
                    }
                    colon_layers.append(layer)
                    if next_boundary_test.get("saturated") == 1:
                        layered_saturation_complete = 1
                        break
                    if not (
                        next_boundary_test.get("status") == "completed"
                        and int(next_boundary_test.get("boundary_gcd_degree", 0))
                        > 0
                    ):
                        break
                    component_support, component_support_exponents = (
                        boundary_gcd_supported_on_factors(
                            next_boundary_test,
                            [
                                supporting_factor
                                for supporting_factor in factors
                                if supporting_factor["sha256"] != factor["sha256"]
                            ],
                            variables,
                            prime,
                        )
                    )
                    if component_support and len(component_support_exponents) == 1:
                        support_text = next(iter(component_support_exponents))
                        supporting_factor = next(
                            supporting_factor
                            for supporting_factor in factors
                            if supporting_factor["expression_text"] == support_text
                        )
                        component_minor_step, component_minors = (
                            modular_component_minor_generators(
                                ideal_polynomials,
                                factor["expression"],
                                variables,
                                prime,
                                next_boundary_test,
                                excluded_pairs=component_used_pairs,
                                singular=singular,
                                timeout=max(timeout, 60),
                            )
                        )
                        layer["component_support_exponents"] = (
                            component_support_exponents
                        )
                        layer["component_minor_step"] = component_minor_step
                        selected_pair = component_minor_step.get(
                            "selected_minor_pair"
                        )
                        if isinstance(selected_pair, tuple):
                            component_used_pairs.add(selected_pair)
                        if (
                            component_minor_step.get("status") == "completed"
                            and component_minors
                        ):
                            ideal_polynomials.extend(component_minors)
                            ideal_text = ",".join(
                                singular_expression(
                                    polynomial, variables, prime
                                )
                                for polynomial in ideal_polynomials
                            )
                            colon_exponent = max(colon_exponent, 2)
                            layer["colon_step"] = component_minor_step
                            continue
                        component_step, component_basis, component_basis_text = (
                            modular_component_quotient_step(
                                ideal_polynomials,
                                (
                                    factor["expression"],
                                    supporting_factor["expression"],
                                ),
                                variables,
                                prime,
                                singular,
                                timeout,
                            )
                        )
                        layer["colon_step"] = component_step
                        if (
                            component_step.get("status") == "completed"
                            and component_basis
                            and component_basis_text
                        ):
                            ideal_polynomials = component_basis
                            ideal_text = component_basis_text
                            colon_exponent += 1
                            continue
                        break
                    colon_step: dict[str, object]
                    new_generators: list[sp.Poly]
                    if len(ideal_polynomials) == 3:
                        colon_step, new_generators = (
                            modular_structured_triple_linear_colon_step(
                                ideal_polynomials,
                                factor["expression"],
                                variables,
                                prime,
                                singular,
                                timeout,
                                next_boundary_test,
                            )
                        )
                        layer["structured_colon_step"] = colon_step
                    else:
                        colon_step, new_generators = (
                            modular_linear_colon_syzygy_step(
                                ideal_polynomials,
                                factor["expression"],
                                variables,
                                prime,
                                singular,
                                timeout,
                                next_boundary_test,
                            )
                        )
                    layer["colon_step"] = colon_step
                    if colon_step.get("status") != "completed" or not new_generators:
                        break
                    ideal_polynomials.extend(new_generators)
                    ideal_text = ",".join(
                        singular_expression(polynomial, variables, prime)
                        for polynomial in ideal_polynomials
                    )
                    colon_exponent += 1
                if colon_layers:
                    public_boundary_test["colon_layers"] = colon_layers
                if layered_saturation_complete:
                    basis = ideal_text
                    step = {
                        "factor_index": factor_index,
                        "factor_sha256": factor["sha256"],
                        "factor": factor["expression_text"],
                        "engine": "linear_boundary_syzygy_colon_certificate",
                        "basis_kind": "colon_generators",
                        "status": "completed",
                        "seconds": (
                            float(linear_boundary_test.get("seconds", 0.0))
                            + float(residual_boundary_test.get("seconds", 0.0))
                            + sum(
                                float(
                                    layer["boundary_test"].get("seconds", 0.0)
                                )
                                + float(
                                    layer.get("colon_step", {}).get("seconds", 0.0)
                                )
                                for layer in colon_layers
                            )
                        ),
                        "unit": 0,
                        "basis_size": len(ideal_polynomials),
                        "dimension": len(variables) - 2,
                        "quotient_dimension": -1,
                        "saturation_exponent": colon_exponent,
                        "basis_sha256": hashlib.sha256(basis.encode()).hexdigest(),
                        "basis_preview": basis[:4000],
                        "boundary_test": public_boundary_test,
                    }
                    if any(predivision_exponents):
                        step["generator_predivision_exponents"] = (
                            predivision_exponents
                        )
                    steps.append(step)
                    total_seconds += float(step["seconds"])
                    final_basis = basis
                    final_dimension = len(variables) - 2
                    final_basis_size = len(ideal_polynomials)
                    print(
                        f"SAT_FACTOR_SYZYGY_COLON p={prime} "
                        f"index={factor_index}/{len(factors)} "
                        f"dim={final_dimension} exp={colon_exponent} "
                        f"factor={factor['expression_text']}",
                        flush=True,
                    )
                    continue
                if factor["sha256"] not in deferred_factor_hashes:
                    deferred_factor_hashes.add(str(factor["sha256"]))
                    processing_factors.append(factor)
                    basis = ideal_text
                    step = {
                        "factor_index": factor_index,
                        "factor_sha256": factor["sha256"],
                        "factor": factor["expression_text"],
                        "engine": "linear_boundary_partial_colon_deferred",
                        "basis_kind": "colon_generators",
                        "status": "deferred",
                        "seconds": (
                            float(linear_boundary_test.get("seconds", 0.0))
                            + float(residual_boundary_test.get("seconds", 0.0))
                            + sum(
                                float(
                                    layer["boundary_test"].get("seconds", 0.0)
                                )
                                + float(
                                    layer.get("colon_step", {}).get("seconds", 0.0)
                                )
                                for layer in colon_layers
                            )
                        ),
                        "unit": 0,
                        "basis_size": len(ideal_polynomials),
                        "dimension": len(variables) - 2,
                        "quotient_dimension": -1,
                        "saturation_exponent_lower_bound": colon_exponent,
                        "basis_sha256": hashlib.sha256(basis.encode()).hexdigest(),
                        "basis_preview": basis[:4000],
                        "boundary_test": public_boundary_test,
                    }
                    if any(predivision_exponents):
                        step["generator_predivision_exponents"] = (
                            predivision_exponents
                        )
                    steps.append(step)
                    total_seconds += float(step["seconds"])
                    final_basis = basis
                    final_dimension = len(variables) - 2
                    final_basis_size = len(ideal_polynomials)
                    print(
                        f"SAT_FACTOR_DEFERRED p={prime} "
                        f"index={factor_index}/{len(processing_factors)} "
                        f"exp>={colon_exponent} factor={factor['expression_text']}",
                        flush=True,
                    )
                    continue
        quotient_dimension_setup = (
            [
                "int quotient_dimension=-1;",
                "if (dim(G)==0) { quotient_dimension=vdim(G); }",
            ]
            if report_quotient_dimension
            else []
        )
        quotient_dimension_output = (
            ['print("quotient_dimension="+string(quotient_dimension));']
            if report_quotient_dimension
            else []
        )
        factor_support = set(factor["expression"].free_symbols)
        step_ring_variables = tuple(
            variable for variable in variables if variable not in factor_support
        ) + tuple(variable for variable in variables if variable in factor_support)
        source = "\n".join(
            [
                'LIB "elim.lib";',
                "option(redSB);",
                f"ring R={prime},({','.join(map(str, step_ring_variables))}),dp;",
                f"ideal I={ideal_text};",
                (
                    "poly jsat="
                    + singular_expression(factor["expression"], variables, prime)
                    + ";"
                ),
                "ideal J=jsat;",
                "timer=1;",
                "list L=sat_with_exp(I,J);",
                "ideal G=L[1];",
                "int elapsed=timer;",
                "int unit=0;",
                "if ((size(G)==1)&&(G[1]==1)) { unit=1; }",
                *quotient_dimension_setup,
                'print("RESULT_BEGIN");',
                'print("unit="+string(unit));',
                'print("basis_size="+string(size(G)));',
                'print("dimension="+string(dim(G)));',
                *quotient_dimension_output,
                'print("saturation_exponent="+string(L[2]));',
                'print("elapsed_ticks="+string(elapsed));',
                'print("BASIS_BEGIN"); print(G); print("BASIS_END");',
                'print("RESULT_END"); quit;',
            ]
        ) + "\n"
        step: dict[str, object] = {
            "factor_index": factor_index,
            "factor_sha256": factor["sha256"],
            "factor": factor["expression_text"],
            "input_sha256": hashlib.sha256(source.encode()).hexdigest(),
            "ring_variable_order": [str(variable) for variable in step_ring_variables],
        }
        if any(predivision_exponents):
            step["generator_predivision_exponents"] = predivision_exponents
        if linear_boundary_test is not None:
            step["boundary_test"] = public_boundary_test
        if hypersurface_test is not None:
            step["hypersurface_test"] = hypersurface_test
        print(
            f"SAT_FACTOR p={prime} index={factor_index}/{len(factors)} "
            f"terms={factor['terms']} degree={factor['total_degree']} "
            f"factor={factor['expression_text']}",
            flush=True,
        )
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
            step.update(
                {
                    "status": "timeout",
                    "seconds": time.monotonic() - started,
                    "stdout_tail": (error.stdout or "")[-2000:],
                    "stderr": (error.stderr or "")[-2000:],
                }
            )
            total_seconds += float(step["seconds"])
            steps.append(step)
            status = "timeout"
            print(
                f"SAT_FACTOR_TIMEOUT p={prime} index={factor_index}",
                flush=True,
            )
            break
        step["seconds"] = time.monotonic() - started
        total_seconds += float(step["seconds"])
        step["returncode"] = completed.returncode
        has_result_markers = (
            "RESULT_BEGIN" in completed.stdout
            and "RESULT_END" in completed.stdout
        )
        step["status"] = (
            "completed"
            if completed.returncode == 0 and has_result_markers
            else "failed"
        )
        for line in completed.stdout.splitlines():
            stripped = line.strip()
            for key in (
                "unit",
                "basis_size",
                "dimension",
                "quotient_dimension",
                "saturation_exponent",
                "elapsed_ticks",
            ):
                if stripped.startswith(key + "="):
                    step[key] = int(stripped.split("=", 1)[1])
        basis = ""
        if "BASIS_BEGIN" in completed.stdout and "BASIS_END" in completed.stdout:
            basis = (
                completed.stdout.split("BASIS_BEGIN", 1)[1]
                .split("BASIS_END", 1)[0]
                .strip()
            )
            step["basis_sha256"] = hashlib.sha256(basis.encode()).hexdigest()
            step["basis_preview"] = basis[:4000]
        step["stdout_tail"] = completed.stdout[-2000:]
        step["stderr"] = completed.stderr[-2000:]
        steps.append(step)
        if step["status"] != "completed" or not basis:
            status = "failed"
            break
        ideal_text = basis
        basis_entries = [
            entry.strip()
            for entry in basis.replace(",\r\n", ",\n").split(",\n")
            if entry.strip()
        ]
        ideal_polynomials = [
            sp.Poly(
                sp.sympify(
                    entry.replace("^", "**"),
                    locals={str(variable): variable for variable in variables},
                ),
                *variables,
                modulus=prime,
            )
            for entry in basis_entries
        ]
        final_basis = basis
        final_dimension = step.get("dimension")
        final_basis_size = step.get("basis_size")
        unit = int(step.get("unit", 0))
        print(
            f"SAT_FACTOR_DONE p={prime} index={factor_index} "
            f"dim={final_dimension} basis={final_basis_size} "
            f"exp={step.get('saturation_exponent')}",
            flush=True,
        )
        if unit == 1:
            break

    if not factors:
        status = "failed"
    return {
        "prime": prime,
        "method": "successive_factor_quotients",
        "status": status,
        "cutoff": cutoff if status == "completed" and unit == 1 else 0,
        "attempted_cutoff": cutoff,
        "unit": unit,
        "dimension": final_dimension,
        "basis_size": final_basis_size,
        "basis_sha256": hashlib.sha256(final_basis.encode()).hexdigest(),
        "basis_preview": final_basis[:4000],
        "seconds": total_seconds,
        "factors_total": factor_count,
        "factors_completed": len(
            {
                str(step.get("factor_sha256"))
                for step in steps
                if step.get("status") == "completed"
            }
        ),
        "factors_component_classified": len(
            {
                str(step.get("factor_sha256"))
                for step in steps
                if step.get("status") == "completed"
            }
            | {
                str(dependency["linear_factor_sha256"])
                for cycle in component_cycles
                if cycle.get("support_closed") == 1
                for dependency in cycle["linear_dependencies"]
            }
            | {
                str(cycle["hypersurface_factor_sha256"])
                for cycle in component_cycles
                if cycle.get("support_closed") == 1
            }
        ),
        "factor_order": factor_order,
        "component_cycles": component_cycles,
        "steps": steps,
    }


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


def exact_two_normal_specialized_primary_test(
    equations: dict[int, sp.Expr],
    hypersurface_factor: sp.Expr,
    collision_factor: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    singular: str,
    timeout: int,
) -> dict[str, object]:
    """Promote the deterministic two-normal fiber to characteristic zero."""
    started = time.monotonic()
    hypersurface_pivot = next(
        variable for variable in variables if str(variable) == "a0"
    )
    collision_pivot = next(
        variable for variable in variables if str(variable) == "mu"
    )
    parameters = tuple(
        variable
        for variable in variables
        if variable not in {hypersurface_pivot, collision_pivot}
    )
    parameter_values = {
        parameter: sp.Rational(value)
        for parameter, value in zip(parameters, (2, 3, 5))
    }
    hypersurface_in_pivot = sp.Poly(hypersurface_factor, hypersurface_pivot)
    collision_in_pivot = sp.Poly(collision_factor, collision_pivot)
    c_value = sp.Rational(
        hypersurface_in_pivot.coeff_monomial(hypersurface_pivot).subs(
            parameter_values
        )
    )
    r_value = sp.Rational(
        hypersurface_in_pivot.coeff_monomial(1).subs(parameter_values)
    )
    collision_coefficient = sp.Rational(
        collision_in_pivot.coeff_monomial(collision_pivot)
    )
    collision_remainder_value = sp.Rational(
        collision_in_pivot.coeff_monomial(1).subs(parameter_values)
    )
    if c_value == 0 or collision_coefficient == 0:
        return {"status": "failed", "reason": "normal coordinate degenerates"}
    selected_polynomials = sorted(
        [sp.Poly(equation, *variables, domain=sp.QQ) for equation in equations.values()],
        key=lambda polynomial: len(polynomial.terms()),
    )[:2]
    hypersurface_index = variables.index(hypersurface_pivot)
    collision_index = variables.index(collision_pivot)
    parameter_indices = [variables.index(parameter) for parameter in parameters]
    dd_symbol, ee_symbol = sp.symbols("dd ee")

    def local_polynomial(polynomial: sp.Poly) -> sp.Poly:
        hypersurface_degree = polynomial.degree(hypersurface_pivot)
        coefficients: dict[tuple[int, int], sp.Rational] = {}
        for monomial, coefficient in polynomial.terms():
            parameter_coefficient = sp.Rational(coefficient)
            for parameter, parameter_index in zip(
                parameters, parameter_indices, strict=True
            ):
                parameter_coefficient *= parameter_values[parameter] ** monomial[
                    parameter_index
                ]
            hypersurface_power = monomial[hypersurface_index]
            collision_power = monomial[collision_index]
            common_coefficient = (
                parameter_coefficient
                * c_value ** (hypersurface_degree - hypersurface_power)
                * (1 / collision_coefficient) ** collision_power
            )
            for dd_power in range(hypersurface_power + 1):
                dd_coefficient = (
                    math.comb(hypersurface_power, dd_power)
                    * (-r_value) ** (hypersurface_power - dd_power)
                )
                for ee_power in range(collision_power + 1):
                    ee_coefficient = (
                        math.comb(collision_power, ee_power)
                        * (-collision_remainder_value)
                        ** (collision_power - ee_power)
                    )
                    local_monomial = (dd_power, ee_power)
                    coefficients[local_monomial] = sp.Rational(
                        coefficients.get(local_monomial, 0)
                        + common_coefficient * dd_coefficient * ee_coefficient
                    )
        return sp.Poly.from_dict(
            {
                monomial: coefficient
                for monomial, coefficient in coefficients.items()
                if coefficient
            },
            dd_symbol,
            ee_symbol,
            domain=sp.QQ,
        )

    local_polynomials = [local_polynomial(polynomial) for polynomial in selected_polynomials]

    def qq_expression(polynomial: sp.Poly) -> str:
        terms: list[str] = []
        for monomial, coefficient in polynomial.terms():
            rational = sp.Rational(coefficient)
            coefficient_text = (
                str(rational.p)
                if rational.q == 1
                else f"({rational.p}/{rational.q})"
            )
            factors = [coefficient_text]
            for variable, power in zip(
                (dd_symbol, ee_symbol), monomial, strict=True
            ):
                if power == 1:
                    factors.append(str(variable))
                elif power > 1:
                    factors.append(f"{variable}^{power}")
            terms.append("*".join(factors))
        return "+".join(f"({term})" for term in terms) if terms else "0"

    expressions = [qq_expression(polynomial) for polynomial in local_polynomials]
    source = "\n".join(
        [
            "ring L=0,(dd,ee),ds;",
            "short=0;",
            *[
                f"poly f{index}={expression};"
                for index, expression in enumerate(expressions)
            ],
            "timer=1;",
            "ideal G=std(ideal(f0,f1));",
            "int elapsed=timer;",
            "int local_dimension=dim(G);",
            "int local_length=-1;",
            "if (local_dimension==0) { local_length=vdim(G); }",
            'print("RESULT_BEGIN");',
            'print("basis_size="+string(size(G)));',
            'print("local_dimension="+string(local_dimension));',
            'print("local_length="+string(local_length));',
            'print("elapsed_ticks="+string(elapsed));',
            'print("BASIS_BEGIN"); print(G); print("BASIS_END");',
            'print("RESULT_END"); quit;',
        ]
    ) + "\n"
    record: dict[str, object] = {
        "engine": "two_normal_exact_specialized_local_primary",
        "field": "Q",
        "parameter_values": {
            str(parameter): str(value)
            for parameter, value in parameter_values.items()
        },
        "selected_input_term_counts": [
            len(polynomial.terms()) for polynomial in selected_polynomials
        ],
        "local_polynomial_term_counts": [
            len(polynomial.terms()) for polynomial in local_polynomials
        ],
        "input_sha256": hashlib.sha256(source.encode()).hexdigest(),
    }
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
    markers = "RESULT_BEGIN" in completed.stdout and "RESULT_END" in completed.stdout
    record["status"] = (
        "completed" if completed.returncode == 0 and markers else "failed"
    )
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        for key in ("basis_size", "local_dimension", "local_length", "elapsed_ticks"):
            if stripped.startswith(key + "="):
                record[key] = int(stripped.split("=", 1)[1])
    if "BASIS_BEGIN" in completed.stdout and "BASIS_END" in completed.stdout:
        basis = (
            completed.stdout.split("BASIS_BEGIN", 1)[1]
            .split("BASIS_END", 1)[0]
            .strip()
        )
        record["basis_sha256"] = hashlib.sha256(basis.encode()).hexdigest()
        record["basis_preview"] = basis[:4000]
    record["exact_curvilinear_length_twelve"] = int(
        record.get("status") == "completed"
        and record.get("local_dimension") == 0
        and record.get("local_length") == 12
    )
    record["stdout_tail"] = completed.stdout[-2000:]
    record["stderr"] = completed.stderr[-2000:]
    return record


def exact_cleared_linear_hypersurface_polynomial(
    polynomial: sp.Poly,
    pivot: sp.Symbol,
    leading_coefficient: sp.Poly,
    remainder: sp.Poly,
    variables: tuple[sp.Symbol, ...],
) -> sp.Poly:
    """Return ``C^n P(-R/C)`` over ``QQ`` for ``C*pivot+R=0``."""
    boundary_variables = tuple(variable for variable in variables if variable != pivot)
    univariate = sp.Poly(polynomial.as_expr(), pivot)
    degree = univariate.degree()
    coefficient_powers = [sp.Poly(1, *boundary_variables, domain=sp.QQ)]
    remainder_powers = [sp.Poly(1, *boundary_variables, domain=sp.QQ)]
    for _ in range(degree):
        coefficient_powers.append(coefficient_powers[-1] * leading_coefficient)
        remainder_powers.append(remainder_powers[-1] * (-remainder))
    restricted = sp.Poly(0, *boundary_variables, domain=sp.QQ)
    for power in range(degree + 1):
        coefficient = univariate.coeff_monomial(pivot**power)
        if coefficient == 0:
            continue
        restricted += (
            sp.Poly(coefficient, *boundary_variables, domain=sp.QQ)
            * remainder_powers[power]
            * coefficient_powers[degree - power]
        )
    return restricted


def exact_polynomial_valuation(
    polynomial: sp.Poly,
    factor: sp.Poly,
) -> tuple[int, sp.Poly]:
    """Return the exact factor valuation and the resulting quotient."""
    if factor.is_zero:
        raise ZeroDivisionError("cannot compute valuation at the zero polynomial")
    if factor.total_degree() == 0:
        return 0, polynomial
    valuation = 0
    quotient = polynomial
    while not quotient.is_zero:
        candidate, remainder = quotient.div(factor)
        if not remainder.is_zero:
            break
        quotient = candidate
        valuation += 1
    return valuation, quotient


def exact_solved_hypersurface_boundary_support(
    polynomials: list[sp.Poly],
    hypersurface_factor: sp.Expr,
    localization_factors: list[dict[str, object]],
    variables: tuple[sp.Symbol, ...],
) -> dict[str, object]:
    """Compute the exact common divisor on a solved hypersurface.

    This is the characteristic-zero promotion of the cheap coordinate-
    boundary screen.  It records which already inverted factors account for
    the complete common divisor, but it does not by itself compute the
    ambient colon or rule out lower-codimension torsion on the hypersurface.
    """
    started = time.monotonic()

    def digest(polynomial: sp.Poly) -> str:
        return hashlib.sha256(str(polynomial.monic().as_expr()).encode()).hexdigest()

    pivot_candidates: list[tuple[int, int, str, sp.Symbol]] = []
    for variable in variables:
        univariate = sp.Poly(hypersurface_factor, variable)
        if univariate.degree() != 1:
            continue
        boundary_variables = tuple(item for item in variables if item != variable)
        coefficient = sp.Poly(
            univariate.coeff_monomial(variable),
            *boundary_variables,
            domain=sp.QQ,
        )
        pivot_candidates.append(
            (
                len(coefficient.terms()),
                coefficient.total_degree(),
                str(variable),
                variable,
            )
        )
    if not pivot_candidates:
        return {
            "engine": "exact_solved_hypersurface_boundary_support",
            "status": "failed",
            "reason": "hypersurface has no linear pivot",
        }
    _, _, _, pivot = min(pivot_candidates)
    boundary_variables = tuple(variable for variable in variables if variable != pivot)
    hypersurface_in_pivot = sp.Poly(hypersurface_factor, pivot)
    leading_coefficient = sp.Poly(
        hypersurface_in_pivot.coeff_monomial(pivot),
        *boundary_variables,
        domain=sp.QQ,
    )
    remainder = sp.Poly(
        hypersurface_in_pivot.coeff_monomial(1),
        *boundary_variables,
        domain=sp.QQ,
    )
    restrictions = [
        exact_cleared_linear_hypersurface_polynomial(
            polynomial,
            pivot,
            leading_coefficient,
            remainder,
            variables,
        )
        for polynomial in polynomials
    ]
    exact_gcd = restrictions[0]
    for restriction in restrictions[1:]:
        exact_gcd = sp.gcd(exact_gcd, restriction)
    exact_gcd = exact_gcd.monic()

    restricted_factors: list[tuple[int, int, str, sp.Poly]] = []
    for factor in localization_factors:
        factor_polynomial = sp.Poly(
            factor["expression"], *variables, domain=sp.QQ
        )
        restricted = exact_cleared_linear_hypersurface_polynomial(
            factor_polynomial,
            pivot,
            leading_coefficient,
            remainder,
            variables,
        )
        if restricted.is_zero or restricted.total_degree() == 0:
            continue
        restricted = restricted.monic()
        restricted_factors.append(
            (
                restricted.total_degree(),
                len(restricted.terms()),
                str(factor["expression_text"]),
                restricted,
            )
        )
    support_residual = exact_gcd
    support_valuations: dict[str, int] = {}
    support_restrictions: dict[str, dict[str, object]] = {}
    for _, _, factor_text, restricted in sorted(
        restricted_factors,
        key=lambda item: (item[0], item[1], item[2]),
        reverse=True,
    ):
        valuation, support_residual = exact_polynomial_valuation(
            support_residual, restricted
        )
        if valuation:
            support_valuations[factor_text] = valuation
            support_restrictions[factor_text] = {
                "term_count": len(restricted.terms()),
                "total_degree": restricted.total_degree(),
                "sha256": digest(restricted),
            }
    support_residual = support_residual.monic()
    ordered_restricted_factors = sorted(
        restricted_factors,
        key=lambda item: (item[0], item[1], item[2]),
        reverse=True,
    )
    individual_chart_valuations: list[dict[str, int]] = []
    individual_chart_residuals: list[sp.Poly] = []
    for restriction in restrictions:
        residual = restriction.monic()
        valuations: dict[str, int] = {}
        for _, _, factor_text, restricted in ordered_restricted_factors:
            valuation, residual = exact_polynomial_valuation(residual, restricted)
            if valuation:
                valuations[factor_text] = valuation
        individual_chart_valuations.append(valuations)
        individual_chart_residuals.append(residual.monic())

    collision_specializations: list[dict[str, object]] = []
    collision_records = [
        factor
        for factor in localization_factors
        if str(factor["expression_text"]) == "lam - mu"
    ]
    if len(collision_records) == 1:
        collision_expression = collision_records[0]["expression"]
        collision_support = [
            variable
            for variable in variables
            if variable in collision_expression.free_symbols
        ]
        collision_pivot = collision_support[-1]
        collision_in_pivot = sp.Poly(collision_expression, collision_pivot)
        collision_coefficient = collision_in_pivot.coeff_monomial(collision_pivot)
        collision_remainder = collision_in_pivot.coeff_monomial(1)
        if not collision_coefficient.free_symbols:
            collision_substitution = sp.cancel(
                -collision_remainder / collision_coefficient
            )
            collision_boundary_variables = tuple(
                variable
                for variable in boundary_variables
                if variable != collision_pivot
            )
            for index, residual in enumerate(individual_chart_residuals):
                if individual_chart_valuations[index].get("lam - mu", 0):
                    continue
                specialized = sp.Poly(
                    residual.as_expr().subs(
                        collision_pivot, collision_substitution
                    ),
                    *collision_boundary_variables,
                    domain=sp.QQ,
                )
                if specialized.is_zero:
                    continue
                specialized = specialized.monic()
                factorization = []
                _, factor_data = sp.factor_list(specialized)
                for factor, exponent in factor_data:
                    factor = factor.monic()
                    factorization.append(
                        {
                            "exponent": exponent,
                            "term_count": len(factor.terms()),
                            "total_degree": factor.total_degree(),
                            "variable_degrees": {
                                str(variable): factor.degree(variable)
                                for variable in collision_boundary_variables
                            },
                            "sha256": digest(factor),
                            "expression": str(factor.as_expr()),
                        }
                    )
                collision_specializations.append(
                    {
                        "generator_index": index,
                        "pivot": str(collision_pivot),
                        "substitution": str(collision_substitution),
                        "term_count": len(specialized.terms()),
                        "total_degree": specialized.total_degree(),
                        "variable_degrees": {
                            str(variable): specialized.degree(variable)
                            for variable in collision_boundary_variables
                        },
                        "sha256": digest(specialized),
                        "factorization": factorization,
                    }
                )
    localized_factor_texts = {
        str(factor["expression_text"]) for factor in localization_factors
    }
    collision_nilpotent_generators = [
        {
            "generator_index": index,
            "collision_order": valuations["lam - mu"],
        }
        for index, valuations in enumerate(individual_chart_valuations)
        if valuations.get("lam - mu", 0) > 0
        and individual_chart_residuals[index].total_degree() == 0
    ]
    collision_transverse_units = [
        specialization
        for specialization in collision_specializations
        if all(
            str(factor["expression"]) in localized_factor_texts
            for factor in specialization["factorization"]
        )
    ]
    localized_boundary_empty = int(
        bool(collision_nilpotent_generators)
        and bool(collision_transverse_units)
    )
    return {
        "engine": "exact_solved_hypersurface_boundary_support",
        "status": "completed",
        "field": "Q",
        "pivot": str(pivot),
        "leading_coefficient_term_count": len(leading_coefficient.terms()),
        "leading_coefficient_total_degree": leading_coefficient.total_degree(),
        "remainder_term_count": len(remainder.terms()),
        "remainder_total_degree": remainder.total_degree(),
        "input_generator_count": len(polynomials),
        "input_generator_term_counts": [
            len(polynomial.terms()) for polynomial in polynomials
        ],
        "input_generator_pivot_degrees": [
            polynomial.degree(pivot) for polynomial in polynomials
        ],
        "restriction_term_counts": [
            len(restriction.terms()) for restriction in restrictions
        ],
        "restriction_total_degrees": [
            restriction.total_degree() for restriction in restrictions
        ],
        "restriction_sha256": [digest(restriction) for restriction in restrictions],
        "cleared_restriction_identities_by_construction": len(restrictions),
        "exact_gcd_term_count": len(exact_gcd.terms()),
        "exact_gcd_total_degree": exact_gcd.total_degree(),
        "exact_gcd_sha256": digest(exact_gcd),
        "support_valuations": support_valuations,
        "support_restrictions": support_restrictions,
        "support_residual_term_count": len(support_residual.terms()),
        "support_residual_total_degree": support_residual.total_degree(),
        "support_residual_sha256": digest(support_residual),
        "support_complete": int(support_residual.total_degree() == 0),
        "individual_chart_factor_valuations": individual_chart_valuations,
        "individual_chart_residual_term_counts": [
            len(residual.terms()) for residual in individual_chart_residuals
        ],
        "individual_chart_residual_total_degrees": [
            residual.total_degree() for residual in individual_chart_residuals
        ],
        "individual_chart_residual_sha256": [
            digest(residual) for residual in individual_chart_residuals
        ],
        "collision_specializations": collision_specializations,
        "collision_nilpotent_generators": collision_nilpotent_generators,
        "collision_transverse_units": [
            {
                "generator_index": specialization["generator_index"],
                "sha256": specialization["sha256"],
                "factorization": specialization["factorization"],
            }
            for specialization in collision_transverse_units
        ],
        "localized_boundary_empty": localized_boundary_empty,
        "localized_boundary_identity": (
            "if f=e^n*v and g=u+e*w with u,v chart units, the finite "
            "geometric-series identity puts u^n in (f,g)"
            if localized_boundary_empty
            else None
        ),
        "scope_warning": (
            "boundary common-divisor support only; ambient colon and "
            "lower-codimension hypersurface torsion are not classified"
        ),
        "seconds": time.monotonic() - started,
    }


def exact_two_normal_generic_primary_test(
    equations: dict[int, sp.Expr],
    hypersurface_factor: sp.Expr,
    collision_factor: sp.Expr,
    companion_linear_factor: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    singular: str = "Singular",
    timeout: int = 300,
) -> dict[str, object]:
    """Certify the generic normal profile over a rational-function field.

    The hypersurface is solved linearly.  Exact gcds of the two cleared
    restrictions classify its remaining divisorial support over ``QQ``.
    A two-by-two normal-jet determinant then computes the intersection
    length at the generic ``(hypersurface, collision)`` point without a
    Groebner basis over a rational-function coefficient field.
    """
    started = time.monotonic()
    pivot_candidates: list[tuple[int, int, str, sp.Symbol]] = []
    for variable in variables:
        univariate = sp.Poly(hypersurface_factor, variable)
        if univariate.degree() != 1:
            continue
        boundary_variables = tuple(item for item in variables if item != variable)
        coefficient = sp.Poly(
            univariate.coeff_monomial(variable),
            *boundary_variables,
            domain=sp.QQ,
        )
        pivot_candidates.append(
            (
                len(coefficient.terms()),
                coefficient.total_degree(),
                str(variable),
                variable,
            )
        )
    if not pivot_candidates:
        return {"status": "failed", "reason": "hypersurface has no linear pivot"}
    pivot = min(pivot_candidates)[-1]
    boundary_variables = tuple(variable for variable in variables if variable != pivot)
    hypersurface_univariate = sp.Poly(hypersurface_factor, pivot)
    leading_coefficient = sp.Poly(
        hypersurface_univariate.coeff_monomial(pivot),
        *boundary_variables,
        domain=sp.QQ,
    )
    hypersurface_remainder = sp.Poly(
        hypersurface_univariate.coeff_monomial(1),
        *boundary_variables,
        domain=sp.QQ,
    )

    collision_support = [
        variable
        for variable in boundary_variables
        if variable in collision_factor.free_symbols
    ]
    if not collision_support:
        return {"status": "failed", "reason": "constant collision factor"}
    collision_pivot = collision_support[-1]
    collision_univariate = sp.Poly(collision_factor, collision_pivot)
    if collision_univariate.degree() != 1:
        return {"status": "failed", "reason": "collision factor is not linear"}
    collision_coefficient = sp.Rational(
        collision_univariate.coeff_monomial(collision_pivot)
    )
    if collision_coefficient == 0 or collision_coefficient.free_symbols:
        return {
            "status": "failed",
            "reason": "collision pivot coefficient is not a nonzero constant",
        }
    collision_origin = sp.expand(
        -collision_univariate.coeff_monomial(1) / collision_coefficient
    )
    parameters = tuple(
        variable for variable in boundary_variables if variable != collision_pivot
    )

    selected = sorted(
        [
            (order, sp.Poly(equation, *variables, domain=sp.QQ))
            for order, equation in equations.items()
        ],
        key=lambda item: len(item[1].terms()),
    )[:2]
    if len(selected) != 2:
        return {"status": "failed", "reason": "two equations are required"}
    restrictions = [
        exact_cleared_linear_hypersurface_polynomial(
            polynomial,
            pivot,
            leading_coefficient,
            hypersurface_remainder,
            variables,
        )
        for _, polynomial in selected
    ]
    collision_polynomial = sp.Poly(
        collision_factor, *boundary_variables, domain=sp.QQ
    )
    collision_data = [
        exact_polynomial_valuation(restriction, collision_polynomial)
        for restriction in restrictions
    ]
    collision_valuations = [valuation for valuation, _ in collision_data]
    minimum_collision_valuation = min(collision_valuations)

    exact_gcd = sp.gcd(restrictions[0], restrictions[1])
    pivot_coordinate_restriction = exact_cleared_linear_hypersurface_polynomial(
        sp.Poly(pivot, *variables, domain=sp.QQ),
        pivot,
        leading_coefficient,
        hypersurface_remainder,
        variables,
    )
    companion_restriction = exact_cleared_linear_hypersurface_polynomial(
        sp.Poly(companion_linear_factor, *variables, domain=sp.QQ),
        pivot,
        leading_coefficient,
        hypersurface_remainder,
        variables,
    )
    support_factors = [
        (str(pivot), pivot_coordinate_restriction),
        (str(companion_linear_factor), companion_restriction),
        (str(collision_factor), collision_polynomial),
    ]
    support_valuations: dict[str, int] = {}
    support_residual = exact_gcd
    for name, factor in support_factors:
        valuation, support_residual = exact_polynomial_valuation(
            support_residual, factor
        )
        support_valuations[name] = valuation

    specialized_leading_coefficient = sp.Poly(
        leading_coefficient.as_expr().subs(collision_pivot, collision_origin),
        *parameters,
        domain=sp.QQ,
    )
    specialized_hypersurface_remainder = sp.Poly(
        hypersurface_remainder.as_expr().subs(
            collision_pivot, collision_origin
        ),
        *parameters,
        domain=sp.QQ,
    )
    normal_derivatives: list[sp.Poly] = []
    collision_initial_coefficients: list[sp.Poly] = []
    for (_, polynomial), (valuation, quotient) in zip(
        selected, collision_data, strict=True
    ):
        derivative_expression = sp.expand(
            polynomial.diff(pivot).as_expr().subs(
                collision_pivot, collision_origin
            )
        )
        derivative_polynomial = sp.Poly(
            derivative_expression,
            pivot,
            *parameters,
            domain=sp.QQ,
        )
        normal_derivatives.append(
            exact_cleared_linear_hypersurface_polynomial(
                derivative_polynomial,
                pivot,
                specialized_leading_coefficient,
                specialized_hypersurface_remainder,
                (pivot, *parameters),
            )
        )
        if valuation == minimum_collision_valuation:
            initial = sp.Poly(
                quotient.as_expr().subs(collision_pivot, collision_origin),
                *parameters,
                domain=sp.QQ,
            )
        else:
            initial = sp.Poly(0, *parameters, domain=sp.QQ)
        collision_initial_coefficients.append(initial)

    normal_jet_determinant = (
        normal_derivatives[0] * collision_initial_coefficients[1]
        - normal_derivatives[1] * collision_initial_coefficients[0]
    )
    determinant_chart_factors = [
        ("hypersurface leading coefficient", specialized_leading_coefficient),
        (
            str(pivot),
            sp.Poly(
                pivot_coordinate_restriction.as_expr().subs(
                    collision_pivot, collision_origin
                ),
                *parameters,
                domain=sp.QQ,
            ),
        ),
        (
            str(companion_linear_factor),
            sp.Poly(
                companion_restriction.as_expr().subs(
                    collision_pivot, collision_origin
                ),
                *parameters,
                domain=sp.QQ,
            ),
        ),
    ]
    determinant_chart_factor_valuations: dict[str, int] = {}
    determinant_chart_residual = normal_jet_determinant
    for name, factor in determinant_chart_factors:
        valuation, determinant_chart_residual = exact_polynomial_valuation(
            determinant_chart_residual, factor
        )
        determinant_chart_factor_valuations[name] = valuation

    def normalized(polynomial: sp.Poly) -> sp.Poly:
        return polynomial.clear_denoms()[1].primitive()[1]

    def polynomial_sha256(polynomial: sp.Poly) -> str:
        return hashlib.sha256(str(normalized(polynomial).as_expr()).encode()).hexdigest()

    exceptional_factor = normalized(determinant_chart_residual)
    exceptional_factorization = sp.factor_list(exceptional_factor.as_expr())
    exceptional_irreducible = int(
        len(exceptional_factorization[1]) == 1
        and exceptional_factorization[1][0][1] == 1
        and normalized(
            sp.Poly(
                exceptional_factorization[1][0][0],
                *parameters,
                domain=sp.QQ,
            )
        ).total_degree()
        == exceptional_factor.total_degree()
    )
    exceptional_partial_gcd_degrees = {
        str(parameter): sp.gcd(
            exceptional_factor, exceptional_factor.diff(parameter)
        ).total_degree()
        for parameter in parameters
        if exceptional_factor.degree(parameter) > 0
    }
    lower_index = collision_valuations.index(minimum_collision_valuation)
    higher_index = collision_valuations.index(max(collision_valuations))
    higher_polynomial = selected[higher_index][1]
    higher_normal_derivative_restriction = (
        exact_cleared_linear_hypersurface_polynomial(
            higher_polynomial.diff(pivot),
            pivot,
            leading_coefficient,
            hypersurface_remainder,
            variables,
        )
    )
    higher_normal_jet_coefficients: list[sp.Poly] = []
    higher_normal_jet_remainders: list[sp.Poly] = []
    derivative_expression = higher_normal_derivative_restriction.as_expr()
    for jet_order in range(7):
        coefficient = sp.Poly(
            sp.expand(
                derivative_expression.subs(collision_pivot, collision_origin)
                * (1 / collision_coefficient) ** jet_order
                / math.factorial(jet_order)
            ),
            *parameters,
            domain=sp.QQ,
        )
        _, remainder = coefficient.div(exceptional_factor)
        higher_normal_jet_coefficients.append(coefficient)
        higher_normal_jet_remainders.append(remainder)
        derivative_expression = sp.diff(
            derivative_expression, collision_pivot
        )
    first_exceptional_nonzero_jet = next(
        (
            index
            for index, remainder in enumerate(higher_normal_jet_remainders)
            if not remainder.is_zero
        ),
        None,
    )
    _, lower_initial_exceptional_remainder = (
        collision_initial_coefficients[lower_index].div(exceptional_factor)
    )
    _, lower_normal_exceptional_remainder = normal_derivatives[lower_index].div(
        exceptional_factor
    )
    exceptional_next_factor = normalized(
        higher_normal_jet_remainders[-1]
    )
    exceptional_next_chart_factors: list[tuple[str, sp.Poly]] = [
        (
            str(parameter),
            sp.Poly(parameter, *parameters, domain=sp.QQ),
        )
        for parameter in parameters
    ]
    unit_offset_parameters = [
        parameter for parameter in parameters if str(parameter) == "lam"
    ]
    exceptional_next_chart_factors.extend(
        (
            str(parameter) + " - 1",
            sp.Poly(parameter - 1, *parameters, domain=sp.QQ),
        )
        for parameter in unit_offset_parameters
    )
    exceptional_next_chart_factors.extend(
        [
            (str(pivot), determinant_chart_factors[1][1]),
            (str(companion_linear_factor), determinant_chart_factors[2][1]),
        ]
    )
    exceptional_next_chart_factor_valuations: dict[str, int] = {}
    exceptional_next_chart_residual = exceptional_next_factor
    for name, factor in exceptional_next_chart_factors:
        valuation, exceptional_next_chart_residual = exact_polynomial_valuation(
            exceptional_next_chart_residual, factor
        )
        exceptional_next_chart_factor_valuations[name] = valuation
    exceptional_next_chart_residual = normalized(
        exceptional_next_chart_residual
    )
    exceptional_next_factorization = sp.factor_list(
        exceptional_next_chart_residual.as_expr()
    )
    exceptional_next_irreducible = int(
        len(exceptional_next_factorization[1]) == 1
        and exceptional_next_factorization[1][0][1] == 1
        and normalized(
            sp.Poly(
                exceptional_next_factorization[1][0][0],
                *parameters,
                domain=sp.QQ,
            )
        ).total_degree()
        == exceptional_next_chart_residual.total_degree()
    )

    deep_exceptional_result: dict[str, object] = {
        "status": "not_applicable"
    }
    common_parameters = [
        parameter
        for parameter in parameters
        if exceptional_factor.degree(parameter) > 0
        and exceptional_next_chart_residual.degree(parameter) > 0
    ]
    if len(common_parameters) >= 2:
        elimination_parameter = min(
            common_parameters,
            key=lambda parameter: (
                exceptional_factor.degree(parameter)
                + exceptional_next_chart_residual.degree(parameter),
                str(parameter),
            ),
        )
        resultant_parameters = tuple(
            parameter
            for parameter in parameters
            if parameter != elimination_parameter
            and (
                exceptional_factor.degree(parameter) > 0
                or exceptional_next_chart_residual.degree(parameter) > 0
            )
        )
        if len(resultant_parameters) == 1:
            resultant_parameter = resultant_parameters[0]
            exceptional_resultant = sp.Poly(
                sp.resultant(
                    exceptional_factor.as_expr(),
                    exceptional_next_chart_residual.as_expr(),
                    elimination_parameter,
                ),
                resultant_parameter,
                domain=sp.QQ,
            ).primitive()[1]
            resultant_residual = exceptional_resultant
            resultant_boundary_valuations: dict[str, int] = {}
            for name, factor in (
                (
                    str(resultant_parameter),
                    sp.Poly(resultant_parameter, resultant_parameter, domain=sp.QQ),
                ),
                (
                    str(resultant_parameter) + " - 1",
                    sp.Poly(
                        resultant_parameter - 1,
                        resultant_parameter,
                        domain=sp.QQ,
                    ),
                ),
            ):
                valuation, resultant_residual = exact_polynomial_valuation(
                    resultant_residual, factor
                )
                resultant_boundary_valuations[name] = valuation
            resultant_residual = normalized(resultant_residual)
            resultant_residual_factorization = sp.factor_list(
                resultant_residual.as_expr()
            )
            resultant_residual_irreducible = int(
                len(resultant_residual_factorization[1]) == 1
                and resultant_residual_factorization[1][0][1] == 1
            )
            companion_in_resultant_variables = sp.Poly(
                determinant_chart_factors[2][1].as_expr(),
                elimination_parameter,
                resultant_parameter,
                domain=sp.QQ,
            )
            exceptional_in_resultant_variables = sp.Poly(
                exceptional_factor.as_expr(),
                elimination_parameter,
                resultant_parameter,
                domain=sp.QQ,
            )
            companion_resultant = sp.Poly(
                sp.resultant(
                    exceptional_in_resultant_variables.as_expr(),
                    companion_in_resultant_variables.as_expr(),
                    elimination_parameter,
                ),
                resultant_parameter,
                domain=sp.QQ,
            ).primitive()[1]
            boundary_projection_gcd = sp.gcd(
                resultant_residual, companion_resultant
            )

            def algebraic_extension_expression(
                polynomial: sp.Poly,
                polynomial_variables: tuple[sp.Symbol, ...],
            ) -> str:
                converted = sp.Poly(
                    polynomial.as_expr(),
                    *polynomial_variables,
                    domain=sp.QQ,
                )
                terms: list[str] = []
                for monomial, coefficient in converted.terms():
                    rational = sp.Rational(coefficient)
                    factors = [
                        str(rational.p)
                        if rational.q == 1
                        else f"({rational.p}/{rational.q})"
                    ]
                    for variable, power in zip(
                        polynomial_variables, monomial, strict=True
                    ):
                        variable_name = (
                            "@l"
                            if variable == resultant_parameter
                            else str(variable)
                        )
                        if power == 1:
                            factors.append(variable_name)
                        elif power > 1:
                            factors.append(f"{variable_name}^{power}")
                    terms.append("*".join(factors))
                return "+".join(f"({term})" for term in terms) if terms else "0"

            extension_source = "\n".join(
                [
                    f"ring E=(0,@l),({elimination_parameter}),dp;",
                    "short=0;",
                    "minpoly=number("
                    + algebraic_extension_expression(
                        resultant_residual, (resultant_parameter,)
                    )
                    + ");",
                    "poly h="
                    + algebraic_extension_expression(
                        exceptional_in_resultant_variables,
                        (elimination_parameter, resultant_parameter),
                    )
                    + ";",
                    "poly k="
                    + algebraic_extension_expression(
                        sp.Poly(
                            exceptional_next_chart_residual.as_expr(),
                            elimination_parameter,
                            resultant_parameter,
                            domain=sp.QQ,
                        ),
                        (elimination_parameter, resultant_parameter),
                    )
                    + ";",
                    "poly t="
                    + algebraic_extension_expression(
                        companion_in_resultant_variables,
                        (elimination_parameter, resultant_parameter),
                    )
                    + ";",
                    "timer=1;",
                    "poly g=gcd(h,k);",
                    "poly rem=reduce(t,std(ideal(g)));",
                    "int elapsed=timer;",
                    'print("RESULT_BEGIN");',
                    'print("gcd_degree="+string(deg(g)));',
                    'print("gcd_terms="+string(size(g)));',
                    'print("remainder_terms="+string(size(rem)));',
                    'print("elapsed_ticks="+string(elapsed));',
                    'print("GCD_BEGIN"); print(g); print("GCD_END");',
                    'print("RESULT_END"); quit;',
                    "",
                ]
            )
            extension_record: dict[str, object] = {
                "input_sha256": hashlib.sha256(
                    extension_source.encode()
                ).hexdigest(),
                "input_size": len(extension_source),
            }
            extension_started = time.monotonic()
            try:
                extension_completed = subprocess.run(
                    [singular, "-q"],
                    input=extension_source,
                    text=True,
                    capture_output=True,
                    timeout=timeout,
                    check=False,
                )
            except subprocess.TimeoutExpired as error:
                extension_record.update(
                    {
                        "status": "timeout",
                        "seconds": time.monotonic() - extension_started,
                        "stdout_tail": (error.stdout or "")[-2000:],
                        "stderr": (error.stderr or "")[-2000:],
                    }
                )
            else:
                extension_record.update(
                    {
                        "status": (
                            "completed"
                            if extension_completed.returncode == 0
                            and "RESULT_BEGIN" in extension_completed.stdout
                            and "RESULT_END" in extension_completed.stdout
                            else "failed"
                        ),
                        "returncode": extension_completed.returncode,
                        "seconds": time.monotonic() - extension_started,
                        "stdout_tail": extension_completed.stdout[-2000:],
                        "stderr": extension_completed.stderr[-2000:],
                    }
                )
                for line in extension_completed.stdout.splitlines():
                    stripped = line.strip()
                    for key in (
                        "gcd_degree",
                        "gcd_terms",
                        "remainder_terms",
                        "elapsed_ticks",
                    ):
                        if stripped.startswith(key + "="):
                            extension_record[key] = int(
                                stripped.split("=", 1)[1]
                            )
                if (
                    "GCD_BEGIN" in extension_completed.stdout
                    and "GCD_END" in extension_completed.stdout
                ):
                    extension_gcd = (
                        extension_completed.stdout.split("GCD_BEGIN", 1)[1]
                        .split("GCD_END", 1)[0]
                        .strip()
                    )
                    extension_record["gcd_sha256"] = hashlib.sha256(
                        extension_gcd.encode()
                    ).hexdigest()
                    extension_record["gcd_preview"] = extension_gcd[:4000]
            deep_exceptional_empty = int(
                resultant_residual_irreducible == 1
                and boundary_projection_gcd.degree()
                == resultant_residual.degree()
                and extension_record.get("status") == "completed"
                and extension_record.get("gcd_degree") == 1
                and extension_record.get("remainder_terms") == 0
            )
            deep_exceptional_result = {
                "status": "completed",
                "elimination_parameter": str(elimination_parameter),
                "resultant_parameter": str(resultant_parameter),
                "resultant_term_count": len(exceptional_resultant.terms()),
                "resultant_total_degree": exceptional_resultant.total_degree(),
                "resultant_sha256": polynomial_sha256(exceptional_resultant),
                "resultant_boundary_valuations": resultant_boundary_valuations,
                "resultant_residual_term_count": len(
                    resultant_residual.terms()
                ),
                "resultant_residual_total_degree": (
                    resultant_residual.total_degree()
                ),
                "resultant_residual_sha256": polynomial_sha256(
                    resultant_residual
                ),
                "resultant_residual_irreducible": (
                    resultant_residual_irreducible
                ),
                "companion_resultant_total_degree": (
                    companion_resultant.total_degree()
                ),
                "boundary_projection_gcd_degree": (
                    boundary_projection_gcd.total_degree()
                ),
                "boundary_projection_gcd_sha256": polynomial_sha256(
                    boundary_projection_gcd
                ),
                "algebraic_extension_test": extension_record,
                "deep_exceptional_empty_on_chart": deep_exceptional_empty,
            }

    support_complete = int(
        not support_residual.is_zero
        and support_residual.total_degree() == 0
        and all(support_valuations[name] > 0 for name, _ in support_factors)
    )
    generic_profile_closed = int(
        support_complete == 1
        and minimum_collision_valuation == 12
        and not normal_jet_determinant.is_zero
    )
    exceptional_profile_closed = int(
        generic_profile_closed == 1
        and exceptional_irreducible == 1
        and all(degree == 0 for degree in exceptional_partial_gcd_degrees.values())
        and first_exceptional_nonzero_jet == 6
        and not lower_initial_exceptional_remainder.is_zero
        and not lower_normal_exceptional_remainder.is_zero
        and collision_valuations[higher_index]
        > minimum_collision_valuation + 6
        and 2 * minimum_collision_valuation
        > minimum_collision_valuation + 6
    )
    exceptional_stratification_closed = int(
        exceptional_profile_closed == 1
        and deep_exceptional_result.get(
            "deep_exceptional_empty_on_chart"
        )
        == 1
    )
    return {
        "engine": "two_normal_exact_generic_jet_certificate",
        "status": "completed",
        "field": "Q(" + ",".join(map(str, parameters)) + ")",
        "hypersurface_pivot": str(pivot),
        "collision_pivot": str(collision_pivot),
        "selected_orders": [order for order, _ in selected],
        "selected_input_term_counts": [
            len(polynomial.terms()) for _, polynomial in selected
        ],
        "hypersurface_pivot_degrees": [
            polynomial.degree(pivot) for _, polynomial in selected
        ],
        "restriction_term_counts": [
            len(restriction.terms()) for restriction in restrictions
        ],
        "restriction_total_degrees": [
            restriction.total_degree() for restriction in restrictions
        ],
        "restriction_sha256": [
            polynomial_sha256(restriction) for restriction in restrictions
        ],
        "cleared_restriction_identities_by_construction": len(restrictions),
        "collision_valuations": collision_valuations,
        "minimum_collision_valuation": minimum_collision_valuation,
        "exact_gcd_term_count": len(exact_gcd.terms()),
        "exact_gcd_total_degree": exact_gcd.total_degree(),
        "exact_gcd_sha256": polynomial_sha256(exact_gcd),
        "support_valuations": support_valuations,
        "support_residual_degree": support_residual.total_degree(),
        "support_residual_sha256": polynomial_sha256(support_residual),
        "support_complete": support_complete,
        "normal_derivative_term_counts": [
            len(polynomial.terms()) for polynomial in normal_derivatives
        ],
        "normal_derivative_sha256": [
            polynomial_sha256(polynomial) for polynomial in normal_derivatives
        ],
        "collision_initial_term_counts": [
            len(polynomial.terms()) for polynomial in collision_initial_coefficients
        ],
        "collision_initial_sha256": [
            polynomial_sha256(polynomial)
            for polynomial in collision_initial_coefficients
        ],
        "normal_jet_determinant_term_count": len(normal_jet_determinant.terms()),
        "normal_jet_determinant_total_degree": (
            normal_jet_determinant.total_degree()
        ),
        "normal_jet_determinant_sha256": polynomial_sha256(
            normal_jet_determinant
        ),
        "normal_jet_determinant_nonzero": int(
            not normal_jet_determinant.is_zero
        ),
        "normal_jet_chart_factor_valuations": (
            determinant_chart_factor_valuations
        ),
        "normal_jet_chart_residual_term_count": len(
            determinant_chart_residual.terms()
        ),
        "normal_jet_chart_residual_total_degree": (
            determinant_chart_residual.total_degree()
        ),
        "normal_jet_chart_residual_sha256": polynomial_sha256(
            determinant_chart_residual
        ),
        "normal_jet_chart_residual_expression": str(
            normalized(determinant_chart_residual).as_expr()
        ),
        "normal_jet_chart_residual_nonzero": int(
            not determinant_chart_residual.is_zero
        ),
        "exceptional_factor_irreducible": exceptional_irreducible,
        "exceptional_factor_partial_gcd_degrees": (
            exceptional_partial_gcd_degrees
        ),
        "exceptional_field": (
            "Frac(Q["
            + ",".join(map(str, parameters))
            + "]/H)"
        ),
        "exceptional_higher_order": selected[higher_index][0],
        "exceptional_higher_normal_derivative_term_count": len(
            higher_normal_derivative_restriction.terms()
        ),
        "exceptional_higher_normal_derivative_sha256": polynomial_sha256(
            higher_normal_derivative_restriction
        ),
        "exceptional_normal_jet_term_counts": [
            len(polynomial.terms())
            for polynomial in higher_normal_jet_coefficients
        ],
        "exceptional_normal_jet_zero_flags": [
            int(polynomial.is_zero)
            for polynomial in higher_normal_jet_coefficients
        ],
        "exceptional_normal_jet_H_divisible": [
            int(polynomial.is_zero)
            for polynomial in higher_normal_jet_remainders
        ],
        "exceptional_normal_jet_remainder_term_counts": [
            len(polynomial.terms())
            for polynomial in higher_normal_jet_remainders
        ],
        "exceptional_normal_jet_remainder_sha256": [
            polynomial_sha256(polynomial)
            for polynomial in higher_normal_jet_remainders
        ],
        "exceptional_first_nonzero_normal_jet": (
            first_exceptional_nonzero_jet
        ),
        "exceptional_next_factor_term_count": len(
            exceptional_next_factor.terms()
        ),
        "exceptional_next_factor_total_degree": (
            exceptional_next_factor.total_degree()
        ),
        "exceptional_next_factor_variable_degrees": {
            str(parameter): exceptional_next_factor.degree(parameter)
            for parameter in parameters
        },
        "exceptional_next_factor_sha256": polynomial_sha256(
            exceptional_next_factor
        ),
        "exceptional_next_factor_expression": str(
            exceptional_next_factor.as_expr()
        ),
        "exceptional_next_chart_factor_valuations": (
            exceptional_next_chart_factor_valuations
        ),
        "exceptional_next_chart_residual_term_count": len(
            exceptional_next_chart_residual.terms()
        ),
        "exceptional_next_chart_residual_total_degree": (
            exceptional_next_chart_residual.total_degree()
        ),
        "exceptional_next_chart_residual_variable_degrees": {
            str(parameter): exceptional_next_chart_residual.degree(parameter)
            for parameter in parameters
        },
        "exceptional_next_chart_residual_sha256": polynomial_sha256(
            exceptional_next_chart_residual
        ),
        "exceptional_next_chart_residual_expression": str(
            exceptional_next_chart_residual.as_expr()
        ),
        "exceptional_next_chart_residual_irreducible": (
            exceptional_next_irreducible
        ),
        "deep_exceptional_result": deep_exceptional_result,
        "exceptional_lower_collision_initial_nonzero": int(
            not lower_initial_exceptional_remainder.is_zero
        ),
        "exceptional_lower_normal_derivative_nonzero": int(
            not lower_normal_exceptional_remainder.is_zero
        ),
        "exceptional_quadratic_normal_terms_start_at": (
            2 * minimum_collision_valuation
        ),
        "exceptional_generic_normal_length": (
            minimum_collision_valuation + first_exceptional_nonzero_jet
            if exceptional_profile_closed
            and first_exceptional_nonzero_jet is not None
            else None
        ),
        "exceptional_primary_profile_closed": exceptional_profile_closed,
        "exceptional_stratification_closed": (
            exceptional_stratification_closed
        ),
        "generic_normal_length": (
            minimum_collision_valuation if generic_profile_closed else None
        ),
        "generic_primary_profile_closed": generic_profile_closed,
        "seconds": time.monotonic() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--degrees", type=parse_degrees, default=parse_degrees("2,4,6,8"))
    parser.add_argument("--groups", help="direction partition, for example 0,1|2|3")
    parser.add_argument("--max-order", type=int, default=10)
    parser.add_argument("--primes", nargs="+", type=int, default=list(PRIMES))
    parser.add_argument("--singular", default="Singular")
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument(
        "--normalize-coefficient",
        default="a0",
        help="nonzero coefficient normalized to one (default: a0)",
    )
    parser.add_argument("--exact", action="store_true")
    parser.add_argument(
        "--compile-only",
        action="store_true",
        help="compile and hash moments without starting elimination",
    )
    parser.add_argument(
        "--saturation-method",
        choices=("rabinowitsch", "quotient", "successive"),
        default="rabinowitsch",
        help="finite-field localization algorithm",
    )
    parser.add_argument(
        "--report-quotient-dimension",
        action="store_true",
        help="compute vdim when a quotient-saturated basis is zero-dimensional",
    )
    parser.add_argument(
        "--solved-hypersurface-coordinate-colon",
        action="store_true",
        help=(
            "also replace the unique nonlinear localization factor by its "
            "own coordinate and saturate by that coordinate"
        ),
    )
    parser.add_argument(
        "--coordinate-colon-only",
        action="store_true",
        help=(
            "skip the ordinary modular saturation and run only the solved-"
            "hypersurface coordinate colon (implies the preceding option)"
        ),
    )
    parser.add_argument(
        "--coordinate-boundary-only",
        action="store_true",
        help=(
            "compute only the common restriction divisor in the solved "
            "hypersurface coordinate"
        ),
    )
    parser.add_argument(
        "--coordinate-component-quotient-only",
        action="store_true",
        help=(
            "in the solved coordinate, quotient the first two moments by "
            "their certified boundary component before saturating by the "
            "normal coordinate"
        ),
    )
    parser.add_argument(
        "--coordinate-component-minor-only",
        action="store_true",
        help=(
            "adjoin the explicit residual-intersection minor of the first "
            "two moments before saturating by the solved normal coordinate"
        ),
    )
    parser.add_argument(
        "--coordinate-component-minor-boundary-only",
        action="store_true",
        help=(
            "construct the residual-intersection minor and report only its "
            "normal-boundary gcd, without a subsequent Groebner saturation"
        ),
    )
    parser.add_argument(
        "--exact-coordinate-boundary-support",
        action="store_true",
        help=(
            "promote the solved-hypersurface boundary gcd and its support "
            "on already inverted factors over Q"
        ),
    )
    parser.add_argument(
        "--omit-nonlinear-localization-factor",
        action="store_true",
        help=(
            "enlarge the chart by omitting its unique nonlinear localization "
            "factor; a resulting unit certificate is therefore stronger"
        ),
    )
    parser.add_argument(
        "--linear-pivot",
        help="solve the first nonzero moment for this chart variable",
    )
    parser.add_argument(
        "--linear-pivot-boundary",
        help=(
            "replace the first equation A+B*pivot by the exceptional boundary "
            "equations A=B=0"
        ),
    )
    parser.add_argument(
        "--boundary-linear-pivot",
        help=(
            "after --linear-pivot-boundary creates A=B=0, solve A=0 "
            "linearly for this second coefficient"
        ),
    )
    parser.add_argument(
        "--scan-linear-pivots",
        action="store_true",
        help=(
            "record transformed term counts for every chart variable that is "
            "linear in the first nonzero moment"
        ),
    )
    parser.add_argument(
        "--scan-pivot-max-order",
        type=int,
        default=0,
        help="limit pivot-scan transformations to this moment order (zero means all)",
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
    normalized_matches = [
        coefficient
        for coefficient in coefficients
        if str(coefficient) == arguments.normalize_coefficient
    ]
    if len(normalized_matches) != 1:
        raise SystemExit(
            f"unknown normalization coefficient {arguments.normalize_coefficient!r}"
        )
    normalized_coefficient = normalized_matches[0]
    active_coefficients = tuple(
        coefficient
        for coefficient in coefficients
        if coefficient != normalized_coefficient
    )
    chart_variables = active_coefficients + active_parameters
    substitutions = {normalized_coefficient: 1}

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

    saturation = sp.prod(active_coefficients) * configuration_discriminant(
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
    pivot_scan = None
    if arguments.scan_linear_pivots:
        first_order = min(equations)
        scan_equations = {
            order: equation
            for order, equation in equations.items()
            if not arguments.scan_pivot_max_order
            or order <= arguments.scan_pivot_max_order
        }
        pivot_scan = {}
        for variable in chart_variables:
            degree = sp.Poly(equations[first_order], variable).degree()
            if degree != 1:
                pivot_scan[str(variable)] = {
                    "eligible": False,
                    "first_order": first_order,
                    "degree": degree,
                }
                continue
            _, _, _, metadata = linear_pivot_localization(
                scan_equations,
                variable,
                chart_variables,
                saturation,
            )
            pivot_scan[str(variable)] = {
                "eligible": True,
                "degree": degree,
                **metadata,
            }
            print(
                "PIVOT_SCAN",
                variable,
                metadata["transformed_term_counts"],
                flush=True,
            )
    if arguments.linear_pivot and arguments.linear_pivot_boundary:
        raise SystemExit(
            "--linear-pivot and --linear-pivot-boundary are mutually exclusive"
        )
    if arguments.boundary_linear_pivot and not arguments.linear_pivot_boundary:
        raise SystemExit(
            "--boundary-linear-pivot requires --linear-pivot-boundary"
        )
    pivot_boundary_metadata = None
    if arguments.linear_pivot_boundary:
        matching = [
            variable
            for variable in chart_variables
            if str(variable) == arguments.linear_pivot_boundary
        ]
        if len(matching) != 1:
            raise SystemExit(
                f"unknown or ambiguous boundary pivot "
                f"{arguments.linear_pivot_boundary!r}"
            )
        boundary_pivot = matching[0]
        first_order = min(equations)
        first = sp.Poly(equations[first_order], boundary_pivot)
        if first.degree() != 1:
            raise SystemExit(
                f"moment {first_order} is not linear in {boundary_pivot}"
            )
        coefficient = first.coeff_monomial(boundary_pivot)
        remainder = first.coeff_monomial(1)
        coefficient = primitive_polynomial(coefficient, chart_variables)
        remainder = primitive_polynomial(remainder, chart_variables)
        equations = {
            first_order - 1: remainder,
            first_order: coefficient,
            **{
                order: equation
                for order, equation in equations.items()
                if order > first_order
            },
        }
        pivot_boundary_metadata = {
            "pivot": str(boundary_pivot),
            "source_order": first_order,
            "equations": "A=B=0 for A+B*pivot=0",
            "A_sha256": hashlib.sha256(str(remainder).encode()).hexdigest(),
            "B_sha256": hashlib.sha256(str(coefficient).encode()).hexdigest(),
            "A_terms": len(sp.Poly(remainder, *chart_variables).terms()),
            "B_terms": len(sp.Poly(coefficient, *chart_variables).terms()),
            "chart": "B=0; the first moment and pivot!=0 force A=0",
        }
        print(
            "PIVOT_BOUNDARY",
            boundary_pivot,
            f"A_terms={pivot_boundary_metadata['A_terms']}",
            f"B_terms={pivot_boundary_metadata['B_terms']}",
            flush=True,
        )
    boundary_nested_pivot_metadata = None
    if arguments.boundary_linear_pivot:
        matching = [
            variable
            for variable in chart_variables
            if str(variable) == arguments.boundary_linear_pivot
        ]
        if len(matching) != 1:
            raise SystemExit(
                f"unknown or ambiguous nested boundary pivot "
                f"{arguments.boundary_linear_pivot!r}"
            )
        (
            equations,
            chart_variables,
            saturation,
            boundary_nested_pivot_metadata,
        ) = linear_pivot_localization(
            equations,
            matching[0],
            chart_variables,
            saturation,
        )
        print(
            "BOUNDARY_LOCALIZED",
            arguments.boundary_linear_pivot,
            boundary_nested_pivot_metadata["transformed_term_counts"],
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
    omitted_nonlinear_factor = None
    if arguments.omit_nonlinear_localization_factor:
        factor_records_before_omission = saturation_factor_records(
            saturation, chart_variables
        )
        nonlinear_factors = [
            factor
            for factor in factor_records_before_omission
            if int(factor["total_degree"]) > 1
        ]
        if len(nonlinear_factors) != 1:
            raise RuntimeError(
                "nonlinear-factor omission requires exactly one nonlinear "
                "localization factor"
            )
        omitted_nonlinear_factor = nonlinear_factors[0]
        divisor = sp.Poly(
            omitted_nonlinear_factor["expression"]
            ** int(omitted_nonlinear_factor["multiplicity"]),
            *chart_variables,
            domain=sp.QQ,
        )
        saturation_polynomial = sp.Poly(
            saturation, *chart_variables, domain=sp.QQ
        )
        saturation_quotient, saturation_remainder = saturation_polynomial.div(
            divisor
        )
        if not saturation_remainder.is_zero:
            raise RuntimeError("nonlinear localization factor division failed")
        saturation = saturation_quotient.as_expr()
        print(
            "OMIT_NONLINEAR_LOCALIZATION_FACTOR",
            omitted_nonlinear_factor["expression_text"],
            flush=True,
        )
    modular_results = []
    if not arguments.compile_only:
        coordinate_factor = None
        if (
            arguments.solved_hypersurface_coordinate_colon
            or arguments.coordinate_colon_only
            or arguments.coordinate_boundary_only
            or arguments.coordinate_component_quotient_only
            or arguments.coordinate_component_minor_only
            or arguments.coordinate_component_minor_boundary_only
        ):
            nonlinear_factors = [
                factor
                for factor in saturation_factor_records(saturation, chart_variables)
                if int(factor["total_degree"]) > 1
            ]
            if len(nonlinear_factors) != 1:
                raise RuntimeError(
                    "solved-hypersurface coordinate colon requires exactly "
                    "one nonlinear localization factor"
                )
            coordinate_factor = nonlinear_factors[0]
        for prime in arguments.primes:
            print(f"RUN modular p={prime}", flush=True)
            if (
                arguments.coordinate_colon_only
                or arguments.coordinate_boundary_only
                or arguments.coordinate_component_quotient_only
                or arguments.coordinate_component_minor_only
                or arguments.coordinate_component_minor_boundary_only
            ):
                result = {
                    "prime": prime,
                    "method": (
                        "solved_hypersurface_coordinate_boundary_only"
                        if arguments.coordinate_boundary_only
                        else "solved_hypersurface_coordinate_component_quotient_only"
                        if arguments.coordinate_component_quotient_only
                        else "solved_hypersurface_coordinate_component_minor_only"
                        if arguments.coordinate_component_minor_only
                        else "solved_hypersurface_coordinate_component_minor_boundary_only"
                        if arguments.coordinate_component_minor_boundary_only
                        else "solved_hypersurface_coordinate_colon_only"
                    ),
                    "status": "discovery",
                    "cutoff": 0,
                }
            elif arguments.saturation_method == "quotient":
                result = modular_saturation_cutoff(
                    equations,
                    saturation,
                    chart_variables,
                    prime,
                    arguments.singular,
                    arguments.timeout,
                    report_quotient_dimension=arguments.report_quotient_dimension,
                )
            elif arguments.saturation_method == "successive":
                result = modular_successive_saturation_cutoff(
                    equations,
                    saturation,
                    chart_variables,
                    prime,
                    arguments.singular,
                    arguments.timeout,
                    report_quotient_dimension=arguments.report_quotient_dimension,
                )
            else:
                result = modular_cutoff(
                    equations,
                    saturation,
                    chart_variables,
                    prime,
                    arguments.singular,
                    arguments.timeout,
                )
            if coordinate_factor is not None:
                print(f"RUN coordinate colon p={prime}", flush=True)
                coordinate_colon = modular_solved_hypersurface_coordinate_colon(
                    [
                        rational_polynomial_mod(
                            equation, chart_variables, prime
                        )
                        for equation in equations.values()
                    ],
                    coordinate_factor["expression"],
                    chart_variables,
                    prime,
                    arguments.singular,
                    arguments.timeout,
                    boundary_only=arguments.coordinate_boundary_only,
                    component_quotient=(
                        arguments.coordinate_component_quotient_only
                    ),
                    component_minor=(
                        arguments.coordinate_component_minor_only
                        or arguments.coordinate_component_minor_boundary_only
                    ),
                    minor_boundary_only=(
                        arguments.coordinate_component_minor_boundary_only
                    ),
                )
                result["solved_hypersurface_coordinate_colon"] = coordinate_colon
                if (
                    arguments.coordinate_colon_only
                    or arguments.coordinate_boundary_only
                    or arguments.coordinate_component_quotient_only
                    or arguments.coordinate_component_minor_only
                    or arguments.coordinate_component_minor_boundary_only
                ):
                    result["status"] = str(coordinate_colon.get("status", "failed"))
                    result["unit"] = coordinate_colon.get("unit")
                    result["dimension"] = coordinate_colon.get("dimension")
                    result["basis_size"] = coordinate_colon.get("basis_size")
                    result["seconds"] = coordinate_colon.get("seconds", 0.0)
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
    exact_generic_primary_result = None
    exact_local_primary_result = None
    exact_coordinate_boundary_support = None
    if arguments.exact_coordinate_boundary_support:
        if arguments.compile_only:
            raise RuntimeError(
                "--exact-coordinate-boundary-support cannot be combined "
                "with --compile-only"
            )
        factor_records = saturation_factor_records(saturation, chart_variables)
        nonlinear_factors = [
            factor for factor in factor_records if int(factor["total_degree"]) > 1
        ]
        if len(nonlinear_factors) != 1:
            raise RuntimeError(
                "exact coordinate boundary support requires exactly one "
                "nonlinear localization factor"
            )
        exact_coordinate_boundary_support = (
            exact_solved_hypersurface_boundary_support(
                [
                    sp.Poly(equation, *chart_variables, domain=sp.QQ)
                    for equation in equations.values()
                ],
                nonlinear_factors[0]["expression"],
                [
                    factor
                    for factor in factor_records
                    if factor["sha256"] != nonlinear_factors[0]["sha256"]
                ],
                chart_variables,
            )
        )
        print(
            "DONE exact coordinate boundary support "
            f"degree={exact_coordinate_boundary_support.get('exact_gcd_total_degree')} "
            f"complete={exact_coordinate_boundary_support.get('support_complete')}",
            flush=True,
        )
    if (
        modular_results
        and all(
            any(
                cycle.get("primary_profile_closed") == 1
                for cycle in result.get("component_cycles", [])
            )
            for result in modular_results
        )
    ):
        factor_records = saturation_factor_records(saturation, chart_variables)
        hypersurface_candidates = [
            factor
            for factor in factor_records
            if int(factor["total_degree"]) == 17
        ]
        collision_candidates = [
            factor
            for factor in factor_records
            if factor["expression_text"] == "lam - mu"
        ]
        companion_candidates = [
            factor
            for factor in factor_records
            if factor["expression_text"] == "143*a0 + 60*a2"
        ]
        if (
            len(hypersurface_candidates) == 1
            and len(collision_candidates) == 1
            and len(companion_candidates) == 1
        ):
            exact_generic_primary_result = exact_two_normal_generic_primary_test(
                equations,
                hypersurface_candidates[0]["expression"],
                collision_candidates[0]["expression"],
                companion_candidates[0]["expression"],
                chart_variables,
                arguments.singular,
                max(arguments.timeout, 120),
            )
            print(
                "DONE exact generic primary "
                f"length={exact_generic_primary_result.get('generic_normal_length')} "
                f"status={exact_generic_primary_result.get('status')}",
                flush=True,
            )
            exact_local_primary_result = exact_two_normal_specialized_primary_test(
                equations,
                hypersurface_candidates[0]["expression"],
                collision_candidates[0]["expression"],
                chart_variables,
                arguments.singular,
                arguments.timeout,
            )
            print(
                "DONE exact local primary "
                f"length={exact_local_primary_result.get('local_length')} "
                f"status={exact_local_primary_result.get('status')}",
                flush=True,
            )
    exact_result = None
    if arguments.exact:
        if arguments.compile_only:
            raise RuntimeError("--exact cannot be combined with --compile-only")
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
            "moment compilation only; no elimination result"
            if arguments.compile_only
            else "exact characteristic-zero chart exclusion"
            if exact_result and exact_result["unit"] == 1
            else (
                "exact characteristic-zero generic normal certificate with "
                "bounded modular saturation discovery"
            )
            if exact_generic_primary_result
            and exact_generic_primary_result.get(
                "generic_primary_profile_closed"
            )
            == 1
            else "bounded modular discovery; not a characteristic-zero theorem"
        ),
        "degrees": list(degrees),
        "groups": [list(block) for block in groups],
        "directions": directions,
        "normalization": (
            f"{normalized_coefficient}=1; first three direction points are "
            "infinity, zero, one"
        ),
        "configuration_parameters": [str(parameter) for parameter in active_parameters],
        "configuration_discriminant": str(
            configuration_discriminant(directions, active_parameters)
        ),
        "coefficient_saturation": [
            str(coefficient) for coefficient in active_coefficients
        ],
        "linear_pivot_scan": pivot_scan,
        "linear_pivot_scan_max_order": (
            arguments.scan_pivot_max_order if arguments.scan_linear_pivots else None
        ),
        "linear_pivot_localization": pivot_metadata,
        "linear_pivot_boundary": pivot_boundary_metadata,
        "boundary_linear_pivot_localization": boundary_nested_pivot_metadata,
        "specialization": {
            str(variable): str(value) for variable, value in specialization.items()
        },
        "max_order": arguments.max_order,
        "compile_only": arguments.compile_only,
        "saturation_method": arguments.saturation_method,
        "solved_hypersurface_coordinate_colon": bool(
            arguments.solved_hypersurface_coordinate_colon
            or arguments.coordinate_colon_only
            or arguments.coordinate_boundary_only
            or arguments.coordinate_component_quotient_only
            or arguments.coordinate_component_minor_only
            or arguments.coordinate_component_minor_boundary_only
        ),
        "coordinate_colon_only": arguments.coordinate_colon_only,
        "coordinate_boundary_only": arguments.coordinate_boundary_only,
        "coordinate_component_quotient_only": (
            arguments.coordinate_component_quotient_only
        ),
        "coordinate_component_minor_only": (
            arguments.coordinate_component_minor_only
        ),
        "coordinate_component_minor_boundary_only": (
            arguments.coordinate_component_minor_boundary_only
        ),
        "omitted_nonlinear_localization_factor": (
            {
                key: value
                for key, value in omitted_nonlinear_factor.items()
                if key not in {"expression", "sort_key"}
            }
            if omitted_nonlinear_factor is not None
            else None
        ),
        "moment_term_counts": term_counts,
        "compile_seconds": compile_seconds,
        "modular_results": modular_results,
        "exact_generic_primary_result": exact_generic_primary_result,
        "exact_local_primary_result": exact_local_primary_result,
        "exact_coordinate_boundary_support": (
            exact_coordinate_boundary_support
        ),
        "exact_result": exact_result,
        "scope": "one declared direction-partition and nonzero-coefficient chart",
    }
    if arguments.report_quotient_dimension:
        artifact["report_quotient_dimension"] = True
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
