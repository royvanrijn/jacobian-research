#!/usr/bin/env python3
"""Modular interpolation of the degree-seven order-five zero-scheme chart.

The exact function-field presentation has strong rank 68, but solving its
68-by-68 pivot system symbolically causes severe expression swell.  This
script evaluates the fixed pivot chart over one finite field and reconstructs
four residual rational functions by nested univariate interpolation.

The common denominator has total degree 32.  The selected numerator degree
bounds are 34, 36, 33, and 33.  Several extra sigma samples and tau lines are
reserved for exact hold-out validation.  Singular computes the zero-dimensional
ideal of the reconstructed numerators.  This remains modular discovery until
the resulting Groebner data are reconstructed and verified over Q.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import shutil
import subprocess

import sympy as sp
from sympy.polys.domains import GF, QQ
from sympy.polys.matrices.sdm import sdm_irref

from explore_degree_five_quantum_residue import add, scale
from reconstruct_degree_seven_order_five_zero_scheme import PIVOT_COLUMNS
from verify_degree_seven_relative_quantization_obstruction import (
    family_presentation,
)


DENOMINATOR_DEGREE = 32
NUMERATOR_DEGREES = (34, 36, 33, 33)
TARGET_INDICES = (68, 69, 72, 81)


def ambient_data():
    presentation = family_presentation(QQ, QQ.one, QQ.zero)
    monomials = sorted(
        set(presentation["constant"]).union(
            *(set(column) for column in presentation["strong_columns"])
        )
    )
    return tuple(monomials[:68]), tuple(monomials[index] for index in TARGET_INDICES)


PIVOT_MONOMIALS, TARGET_MONOMIALS = ambient_data()


def residual_task(task: tuple[int, int, int]):
    """Evaluate four fixed residual coordinates at one parameter point."""

    prime, sigma_value, tau_value = task
    field = GF(prime)
    presentation = family_presentation(
        field,
        field(sigma_value),
        field(tau_value),
    )
    columns = presentation["strong_columns"]
    basis = [columns[index] for index in PIVOT_COLUMNS]
    rows = {}
    for row_index, monomial in enumerate(PIVOT_MONOMIALS):
        row = {
            column_index: coefficient
            for column_index, column in enumerate(basis)
            if (coefficient := column.get(monomial, field.zero))
        }
        coefficient = presentation["constant"].get(monomial, field.zero)
        if coefficient:
            row[68] = coefficient
        rows[row_index] = row
    reduced, pivots, _ = sdm_irref(rows)
    if pivots != list(range(68)):
        return sigma_value, tau_value, None
    coordinates = [
        reduced[index].get(68, field.zero) for index in range(68)
    ]
    residual = dict(presentation["constant"])
    for coefficient, column in zip(coordinates, basis):
        residual = add(residual, scale(column, -coefficient))
    values = [
        int(residual.get(monomial, field.zero)) % prime
        for monomial in TARGET_MONOMIALS
    ]
    return sigma_value, tau_value, values


def polynomial_evaluate(coefficients: list[int], value: int, prime: int) -> int:
    result = 0
    for coefficient in reversed(coefficients):
        result = (result * value + coefficient) % prime
    return result


def polynomial_interpolate(
    points: list[tuple[int, int]], degree: int, prime: int
) -> list[int]:
    """Newton interpolation converted to the standard coefficient basis."""

    selected = points[: degree + 1]
    if len(selected) < degree + 1:
        raise ValueError("not enough interpolation points")
    polynomial = [0]
    basis = [1]
    for value, target in selected:
        current = polynomial_evaluate(polynomial, value, prime)
        divisor = polynomial_evaluate(basis, value, prime)
        if divisor == 0:
            raise ValueError("repeated interpolation abscissa")
        multiplier = (target - current) * pow(divisor, -1, prime) % prime
        if len(polynomial) < len(basis):
            polynomial.extend([0] * (len(basis) - len(polynomial)))
        for index, coefficient in enumerate(basis):
            polynomial[index] = (
                polynomial[index] + multiplier * coefficient
            ) % prime
        new_basis = [0] * (len(basis) + 1)
        for index, coefficient in enumerate(basis):
            new_basis[index] = (
                new_basis[index] - value * coefficient
            ) % prime
            new_basis[index + 1] = (
                new_basis[index + 1] + coefficient
            ) % prime
        basis = new_basis
    polynomial.extend([0] * (degree + 1 - len(polynomial)))
    polynomial = polynomial[: degree + 1]
    if any(
        polynomial_evaluate(polynomial, value, prime) != target % prime
        for value, target in points
    ):
        raise ValueError("polynomial hold-out validation failed")
    return polynomial


def rational_line_fit(
    samples: list[tuple[int, int]],
    numerator_degree: int,
    denominator_degree: int,
    prime: int,
) -> tuple[list[int], list[int]]:
    """Fit N/D with the leading denominator coefficient normalized to one."""

    field = GF(prime)
    unknowns = numerator_degree + denominator_degree + 1
    rows = {}
    for row_index, (value, target) in enumerate(samples):
        row = {}
        power = 1
        for index in range(numerator_degree + 1):
            coefficient = field(power)
            if coefficient:
                row[index] = coefficient
            power = power * value % prime
        power = 1
        for index in range(denominator_degree):
            coefficient = -field(target * power % prime)
            if coefficient:
                row[numerator_degree + 1 + index] = coefficient
            power = power * value % prime
        rhs = field(target * pow(value, denominator_degree, prime) % prime)
        if rhs:
            row[unknowns] = rhs
        rows[row_index] = row
    reduced, pivots, _ = sdm_irref(rows)
    if pivots != list(range(unknowns)):
        raise ValueError(f"unexpected rational interpolation pivots: {pivots}")
    solution = [
        int(reduced[index].get(unknowns, field.zero)) % prime
        for index in range(unknowns)
    ]
    numerator = solution[: numerator_degree + 1]
    denominator = solution[numerator_degree + 1 :] + [1]
    for value, target in samples:
        denominator_value = polynomial_evaluate(denominator, value, prime)
        if denominator_value == 0:
            continue
        if (
            polynomial_evaluate(numerator, value, prime)
            != target * denominator_value % prime
        ):
            raise ValueError("rational hold-out validation failed")
    return numerator, denominator


def reconstruct_line(
    tau_value: int,
    records: list[tuple[int, list[int]]],
    prime: int,
):
    """Reconstruct D(sigma,tau) and all selected N_j on one tau line."""

    first_samples = [(sigma, values[0]) for sigma, values in records]
    first_numerator, denominator = rational_line_fit(
        first_samples,
        NUMERATOR_DEGREES[0],
        DENOMINATOR_DEGREE,
        prime,
    )
    numerators = [first_numerator]
    for target_index, degree in enumerate(NUMERATOR_DEGREES[1:], start=1):
        polynomial_samples = [
            (
                sigma,
                values[target_index]
                * polynomial_evaluate(denominator, sigma, prime)
                % prime,
            )
            for sigma, values in records
        ]
        numerators.append(
            polynomial_interpolate(polynomial_samples, degree, prime)
        )
    return tau_value, denominator, numerators


def nested_reconstruction(lines, prime: int):
    usable = []
    for tau_value, records in sorted(lines.items()):
        try:
            usable.append(reconstruct_line(tau_value, records, prime))
        except ValueError:
            continue
    required_lines = max(NUMERATOR_DEGREES) + 1
    if len(usable) < required_lines + 3:
        raise AssertionError(
            f"only {len(usable)} usable lines; need {required_lines + 3}"
        )

    denominator_terms = {}
    for sigma_degree in range(DENOMINATOR_DEGREE + 1):
        tau_degree = DENOMINATOR_DEGREE - sigma_degree
        points = [
            (tau_value, denominator[sigma_degree])
            for tau_value, denominator, _ in usable
        ]
        coefficients = polynomial_interpolate(points, tau_degree, prime)
        for degree, coefficient in enumerate(coefficients):
            if coefficient:
                denominator_terms[(sigma_degree, degree)] = coefficient

    numerator_terms = []
    for target_index, total_degree in enumerate(NUMERATOR_DEGREES):
        terms = {}
        for sigma_degree in range(total_degree + 1):
            tau_degree = total_degree - sigma_degree
            points = [
                (tau_value, numerators[target_index][sigma_degree])
                for tau_value, _, numerators in usable
            ]
            coefficients = polynomial_interpolate(points, tau_degree, prime)
            for degree, coefficient in enumerate(coefficients):
                if coefficient:
                    terms[(sigma_degree, degree)] = coefficient
        numerator_terms.append(terms)
    return usable, denominator_terms, numerator_terms


def evaluate_bivariate(terms, sigma: int, tau: int, prime: int) -> int:
    return sum(
        coefficient
        * pow(sigma, sigma_degree, prime)
        * pow(tau, tau_degree, prime)
        for (sigma_degree, tau_degree), coefficient in terms.items()
    ) % prime


def render_polynomial(terms, prime: int) -> str:
    rendered = []
    for (sigma_degree, tau_degree), coefficient in sorted(terms.items()):
        signed = coefficient if coefficient <= prime // 2 else coefficient - prime
        factors = [str(signed)]
        if sigma_degree == 1:
            factors.append("sigma")
        elif sigma_degree > 1:
            factors.append(f"sigma^{sigma_degree}")
        if tau_degree == 1:
            factors.append("tau")
        elif tau_degree > 1:
            factors.append(f"tau^{tau_degree}")
        rendered.append("*".join(factors))
    return "+".join(rendered).replace("+-", "-") if rendered else "0"


def parse_singular_polynomial(text: str, prime: int):
    sigma, tau = sp.symbols("sigma tau")
    expression = sp.sympify(
        text.replace("^", "**"),
        locals={"sigma": sigma, "tau": tau},
    )
    polynomial = sp.Poly(expression, sigma, tau, modulus=prime)
    return {
        monomial: int(coefficient) % prime
        for monomial, coefficient in polynomial.terms()
    }


def singular_groebner(numerators, denominator, prime: int):
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required for the modular chart reconstruction")
    ideal = ",".join(render_polynomial(terms, prime) for terms in numerators)
    denominator_source = render_polynomial(denominator, prime)
    program = f"""LIB "elim.lib";
ring r={prime},(sigma,tau),dp;
option(redSB);
ideal I={ideal};
ideal RAW=std(I);
list L=sat(I,ideal({denominator_source}));
ideal G=std(L[1]);
print("RAW_DIMENSION="+string(dim(RAW)));
print("RAW_VDIM="+string(vdim(RAW)));
print("GB_SIZE="+string(size(G)));
print("DIMENSION="+string(dim(G)));
print("VDIM="+string(vdim(G)));
for (int i=1; i<=size(G); i++) {{ print("GB_"+string(i)+"="+string(G[i])); }}
"""
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    if result.stderr.strip() or "   ?" in result.stdout:
        raise AssertionError(result.stdout + result.stderr)
    values = {}
    basis = []
    for line in result.stdout.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key.startswith("GB_") and key != "GB_SIZE":
            basis.append(value)
        else:
            values[key] = value
    basis_terms = [
        parse_singular_polynomial(polynomial, prime) for polynomial in basis
    ]
    return values, basis, basis_terms, hashlib.sha256(program.encode()).hexdigest()


def serialize_terms(terms):
    return [
        {
            "sigma_degree": sigma_degree,
            "tau_degree": tau_degree,
            "coefficient": coefficient,
        }
        for (sigma_degree, tau_degree), coefficient in sorted(terms.items())
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--sigma-samples", type=int, default=76)
    parser.add_argument("--tau-lines", type=int, default=44)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.prime <= max(args.sigma_samples, args.tau_lines):
        raise SystemExit("prime must exceed both grid dimensions")
    if args.sigma_samples < 70 or args.tau_lines < 40:
        raise SystemExit("grid is too small for the declared degree bounds")

    tasks = [
        (args.prime, sigma, tau)
        for tau in range(args.tau_lines)
        for sigma in range(args.sigma_samples)
    ]
    if args.jobs == 1:
        records = [residual_task(task) for task in tasks]
    else:
        context = mp.get_context("spawn")
        with context.Pool(processes=args.jobs) as pool:
            records = pool.map(residual_task, tasks, chunksize=4)
    lines = {}
    chart_poles = []
    for sigma, tau, values in records:
        if values is None:
            chart_poles.append((sigma, tau))
            continue
        lines.setdefault(tau, []).append((sigma, values))

    usable, denominator, numerators = nested_reconstruction(lines, args.prime)
    # Validate every sampled chart point against the reconstructed functions.
    validated = 0
    for sigma, tau, values in records:
        if values is None:
            continue
        denominator_value = evaluate_bivariate(
            denominator, sigma, tau, args.prime
        )
        if denominator_value == 0:
            continue
        for target, terms in zip(values, numerators):
            assert evaluate_bivariate(terms, sigma, tau, args.prime) == (
                target * denominator_value % args.prime
            )
        validated += 1

    values, basis, basis_terms, digest = singular_groebner(
        numerators, denominator, args.prime
    )
    assert int(values["DIMENSION"]) == 0
    certificate = {
        "scope": "modular pivot-chart interpolation; no characteristic-zero claim",
        "prime": args.prime,
        "grid": {
            "sigma_samples": args.sigma_samples,
            "tau_lines": args.tau_lines,
            "points": len(tasks),
            "chart_poles": len(chart_poles),
            "usable_lines": len(usable),
            "validated_points": validated,
        },
        "degree_bounds": {
            "denominator": DENOMINATOR_DEGREE,
            "numerators": list(NUMERATOR_DEGREES),
        },
        "target_monomials": [list(monomial) for monomial in TARGET_MONOMIALS],
        "denominator_terms": serialize_terms(denominator),
        "numerator_terms": [serialize_terms(terms) for terms in numerators],
        "zero_scheme": {
            "raw_dimension": int(values["RAW_DIMENSION"]),
            "raw_vector_space_dimension": int(values["RAW_VDIM"]),
            "dimension": int(values["DIMENSION"]),
            "vector_space_dimension": int(values["VDIM"]),
            "groebner_basis": basis,
            "groebner_basis_terms": [
                serialize_terms(terms) for terms in basis_terms
            ],
        },
        "singular_program_sha256": digest,
    }
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.output:
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered)

    print(
        "PASS: reconstructed four residual numerators over GF({}) from {} points".format(
            args.prime, len(tasks)
        )
    )
    print(
        "PASS: validated {} chart points on {} usable tau lines".format(
            validated, len(usable)
        )
    )
    print(
        "ZERO SCHEME: dimension={} vdim={} gb_size={}".format(
            values["DIMENSION"], values["VDIM"], values["GB_SIZE"]
        )
    )
    for index, polynomial in enumerate(basis, start=1):
        print(f"GB[{index}] = {polynomial}")
    print("SCOPE: modular interpolation only; reconstruct and verify over Q")


if __name__ == "__main__":
    main()
