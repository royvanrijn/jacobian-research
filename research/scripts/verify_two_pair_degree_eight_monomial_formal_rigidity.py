#!/usr/bin/env python3
"""Good-prime formal-rigidity certificates in balanced degrees eight to eleven."""

from __future__ import annotations

import json
import subprocess
import sys
from math import comb
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from verify_two_pair_higher_degree_monomial_formal_rigidity import (  # noqa: E402
    correction,
    linear_term,
    quadratic_term,
)


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_degree_eight_monomial_formal_rigidity.json"
)
DEGREES = (8, 9, 10, 11)
PRIME = 1000003


def coefficient_mod_prime(coefficient: sp.Rational) -> int:
    numerator, denominator = coefficient.as_numer_denom()
    return (
        int(numerator)
        * pow(int(denominator), -1, PRIME)
    ) % PRIME


def polynomial_mod_prime(
    polynomial: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> str:
    terms = []
    for exponents, coefficient in sp.Poly(
        polynomial,
        *variables,
        domain=sp.QQ,
    ).terms():
        scalar = coefficient_mod_prime(coefficient)
        if scalar == 0:
            continue
        factors = [str(scalar)]
        for variable, exponent in zip(variables, exponents):
            if exponent == 1:
                factors.append(str(variable))
            elif exponent > 1:
                factors.append(f"{variable}^{exponent}")
        terms.append("*".join(factors))
    return "+".join(terms) if terms else "0"


def projected_system(degree: int) -> tuple[
    tuple[sp.Symbol, ...],
    sp.Expr,
    list[int],
    list[int],
    list[sp.Expr],
]:
    even_phases = [
        phase
        for phase in range(1, degree + 1)
        if (degree - phase) % 2 == 0
    ]
    odd_phases = [
        phase
        for phase in range(1, degree + 1)
        if (degree - phase) % 2 == 1
    ]
    even_terms = [
        correction(side, phase, degree)
        for side in ("positive", "negative")
        for phase in even_phases
    ]
    odd_terms = [
        correction(side, phase, degree)
        for side in ("positive", "negative")
        for phase in odd_phases
    ]
    variables = sp.symbols(f"x0:{len(odd_terms)}")
    orders = list(range(2, 2 + len(even_terms) + len(odd_terms)))
    linear_rows = []
    quadratic_coefficients = []
    for order in orders:
        linear_rows.append(
            [
                sp.factor(order * linear_term(order, term))
                for term in even_terms
            ]
        )
        quadratic = 0
        for left in range(len(odd_terms)):
            for right in range(left, len(odd_terms)):
                quadratic += (
                    sp.binomial(order, 2)
                    * (1 if left == right else 2)
                    * quadratic_term(
                        order,
                        odd_terms[left],
                        odd_terms[right],
                    )
                    * variables[left]
                    * variables[right]
                )
        quadratic_coefficients.append(sp.factor(quadratic))

    linear_count = len(even_terms)
    pivot_matrix = sp.Matrix(linear_rows[:linear_count])
    pivot_determinant = sp.factor(pivot_matrix.det())
    assert pivot_determinant
    eliminated_correction = (
        -pivot_matrix.inv()
        * sp.Matrix(quadratic_coefficients[:linear_count])
    )
    projected_obstructions = []
    for index in range(linear_count, len(orders)):
        projected = sp.factor(
            quadratic_coefficients[index]
            + (
                sp.Matrix([linear_rows[index]])
                * eliminated_correction
            )[0]
        )
        _, primitive = sp.Poly(projected, *variables).primitive()
        projected_obstructions.append(primitive.as_expr())
    return (
        variables,
        pivot_determinant,
        orders[:linear_count],
        orders[linear_count:],
        projected_obstructions,
    )


def singular_certificate(
    variables: tuple[sp.Symbol, ...],
    projected_obstructions: list[sp.Expr],
) -> tuple[int, list[str], list[str]]:
    modular_polynomials = [
        polynomial_mod_prime(polynomial, variables)
        for polynomial in projected_obstructions
    ]
    variable_names = ",".join(str(variable) for variable in variables)
    ideal_generators = ",\n".join(modular_polynomials)
    power_reductions = "\n".join(
        f"reduce({variable}^{len(variables) + 1},G);"
        for variable in variables
    )
    singular_source = f"""
ring r={PRIME},({variable_names}),dp;
ideal I={ideal_generators};
ideal G=std(I);
print("BEGIN_VDIM");
vdim(G);
print("END_VDIM");
print("BEGIN_POWERS");
{power_reductions}
print("END_POWERS");
print("BEGIN_LEADS");
lead(G);
print("END_LEADS");
quit;
"""
    result = subprocess.run(
        ["Singular", "-q"],
        input=singular_source,
        text=True,
        capture_output=True,
        check=True,
    )
    lines = [line.strip() for line in result.stdout.splitlines()]
    vdim_start = lines.index("BEGIN_VDIM") + 1
    vdim_end = lines.index("END_VDIM")
    power_start = lines.index("BEGIN_POWERS") + 1
    power_end = lines.index("END_POWERS")
    lead_start = lines.index("BEGIN_LEADS") + 1
    lead_end = lines.index("END_LEADS")
    assert vdim_end == vdim_start + 1
    dimension = int(lines[vdim_start])
    power_values = lines[power_start:power_end]
    leading_terms = lines[lead_start:lead_end]
    return dimension, power_values, leading_terms


def main() -> None:
    certificates = {}
    for degree in DEGREES:
        (
            variables,
            pivot_determinant,
            linear_orders,
            projected_orders,
            projected_obstructions,
        ) = projected_system(degree)
        kernel_count = 2 * (degree // 2)
        assert len(variables) == kernel_count
        linear_count = 2 * ((degree + 1) // 2)
        assert linear_orders == list(range(2, 2 + linear_count))
        assert projected_orders == list(
            range(
                2 + linear_count,
                2 + linear_count + kernel_count,
            )
        )
        pivot_determinant_mod_prime = coefficient_mod_prime(
            sp.Rational(pivot_determinant)
        )
        assert pivot_determinant_mod_prime

        (
            quotient_dimension,
            power_values,
            leading_terms,
        ) = singular_certificate(
            variables,
            projected_obstructions,
        )
        assert quotient_dimension == 2**kernel_count
        assert power_values == ["0"] * kernel_count
        certificates[str(degree)] = {
            "degree": degree,
            "linear_pivot_orders": linear_orders,
            "linear_pivot_determinant": str(pivot_determinant),
            "linear_pivot_determinant_mod_prime": (
                pivot_determinant_mod_prime
            ),
            "projected_quadratic_orders": projected_orders,
            "projected_quadratics_mod_prime": [
                polynomial_mod_prime(polynomial, variables)
                for polynomial in projected_obstructions
            ],
            "quotient_dimension_mod_prime": quotient_dimension,
            "variable_power_reductions": {
                "exponent": kernel_count + 1,
                "remainders": power_values,
            },
            "leading_terms_of_modular_standard_basis": leading_terms,
            "hilbert_vector": [
                comb(kernel_count, index)
                for index in range(kernel_count + 1)
            ],
            "length": 2**kernel_count,
        }

    artifact = {
        "format": "two-pair-degrees-eight-eleven-monomial-formal-rigidity-v1",
        "field": "characteristic zero via good-prime certificate",
        "prime": PRIME,
        "degrees": certificates,
        "written_source": (
            "extended-geometry/"
            "TWO_PAIR_OPPOSITE_MONOMIAL_OBSTRUCTION.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print(
        "PASS degrees 8 through 11: consecutive moment blocks give "
        "invertible linear pivots followed by square projected systems"
    )
    print(
        "PASS good prime 1000003: quotient dimensions are "
        "2^8, 2^8, 2^10, 2^10 and predicted variable powers reduce to zero"
    )
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
