#!/usr/bin/env python3
"""Build the sparse A=0, B!=0 incidence system on the cubic Hurwitz chart.

After the mu_2 pivot eliminates a3, write the primitive reduced mu_3 as

    A*a2^2 + B*a2 + C.

On the declared channel-minor open, A=0 forces b0!=0.  Hence the birational
coordinates q=b1/b0 and z=b0*a2 are valid, while A eliminates b2.  Retaining
z and the resulting linear incidence equation is much sparser than expanding
a2=-C/B.  The chart also exposes and removes the invertible b0 monomial
content of every later moment.  This script emits that generic component
system either over QQ or over one prime field.  It does not claim that the
system is empty or lift a finite-field calculation to characteristic zero.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

from research_two_pair_sic_bidegree33_rank_two_hurwitz import (
    EXACT_LIFT_PRIME,
    RationalParameterPolynomial,
    active_parameter_polynomial_string,
    active_rational_polynomial_string,
    add_rational_polynomial,
    base_polynomials,
    channel_minor_polynomial,
    coefficient_groups,
    discriminant_polynomial,
    exact_moment_polynomials,
    exact_parameter_polynomial,
    linear_substitution_numerator,
    moment,
    multiply_parameter_polynomials,
    multiply_rational_polynomials,
    powers,
    rational_linear_substitution_numerator,
    rational_polynomial,
    split_linear_parameter_polynomial,
    split_linear_rational_polynomial,
    substitute_parameter_variable,
    substitute_rational_variable,
)


VARIABLES = ("r", "a1", "z", "b0", "q", "lambda")
ACTIVE_PARAMETERS = (1, 2, 4, 5, 0)


def normalized_orders(raw: str) -> tuple[int, ...]:
    orders = tuple(sorted({int(value) for value in raw.split(",")}))
    if not orders or orders[0] != 2 or 3 not in orders:
        raise ValueError("orders must contain 2 and 3")
    return orders


def rational_coefficient_groups(
    polynomial: RationalParameterPolynomial, variable: int
) -> dict[int, RationalParameterPolynomial]:
    groups: dict[int, RationalParameterPolynomial] = {}
    for exponent, coefficient in polynomial.items():
        reduced_exponent = list(exponent)
        power = reduced_exponent[variable]
        reduced_exponent[variable] = 0
        add_rational_polynomial(
            groups.setdefault(power, {}),
            {tuple(reduced_exponent): coefficient},
        )
    return groups


def rational_ratio_replacements() -> tuple[
    RationalParameterPolynomial, RationalParameterPolynomial
]:
    b1_replacement = {
        (0, 0, 0, 0, 1, 1, 0): Fraction(1),
    }
    b2_replacement = {
        (0, 0, 0, 0, 0, 0, 0): Fraction(-1),
        (2, 0, 0, 0, 1, 0, 0): Fraction(3, 8),
        (1, 0, 0, 0, 1, 0, 0): Fraction(1, 8),
        (0, 0, 0, 0, 1, 0, 0): Fraction(3, 4),
        (1, 0, 0, 0, 1, 1, 0): Fraction(1, 8),
        (0, 0, 0, 0, 1, 1, 0): Fraction(5, 8),
        (0, 0, 0, 0, 1, 2, 0): Fraction(1, 2),
    }
    return b1_replacement, b2_replacement


def finite_ratio_replacements(
    prime: int,
) -> tuple[dict[tuple[int, ...], int], dict[tuple[int, ...], int]]:
    inverse_eight = pow(8, -1, prime)
    b1_replacement = {
        (0, 0, 0, 0, 1, 1, 0): 1,
    }
    b2_replacement = {
        (0, 0, 0, 0, 0, 0, 0): -1 % prime,
        (2, 0, 0, 0, 1, 0, 0): 3 * inverse_eight % prime,
        (1, 0, 0, 0, 1, 0, 0): inverse_eight,
        (0, 0, 0, 0, 1, 0, 0): 6 * inverse_eight % prime,
        (1, 0, 0, 0, 1, 1, 0): inverse_eight,
        (0, 0, 0, 0, 1, 1, 0): 5 * inverse_eight % prime,
        (0, 0, 0, 0, 1, 2, 0): 4 * inverse_eight % prime,
    }
    return b1_replacement, b2_replacement


def rescale_a2_by_b0(
    polynomial: dict[tuple[int, ...], int | Fraction],
    prime: int | None = None,
) -> tuple[dict[tuple[int, ...], int | Fraction], int]:
    """Substitute a2=z/b0, clear minimally, and remove b0 content."""

    maximum = max((exponent[2] for exponent in polynomial), default=0)
    transformed: dict[tuple[int, ...], int | Fraction] = {}
    for exponent, coefficient in polynomial.items():
        target = list(exponent)
        target[4] += maximum - target[2]
        key = tuple(target)
        value = transformed.get(key, 0) + coefficient
        if prime is not None:
            value %= prime
        if value:
            transformed[key] = value
        else:
            transformed.pop(key, None)
    b0_content = min(
        (exponent[4] for exponent in transformed), default=0
    )
    if b0_content:
        transformed = {
            tuple(
                power - (b0_content if index == 4 else 0)
                for index, power in enumerate(exponent)
            ): coefficient
            for exponent, coefficient in transformed.items()
        }
    return transformed, b0_content


def finite_field_source(
    orders: tuple[int, ...],
    prime: int,
    f4sat: bool,
    coefficient_branch: str,
) -> tuple[str, dict[str, object]]:
    maximum = max(orders)
    b_base, d_base = base_polynomials(prime)
    b_powers = powers(b_base, maximum, prime)
    d_powers = powers(d_base, maximum, prime)
    moments = {
        order: moment(order, b_powers, d_powers, prime)
        for order in orders
    }
    first_pivot, first_rest = split_linear_parameter_polynomial(
        moments[2], 3, prime
    )
    eliminated = {
        order: linear_substitution_numerator(
            polynomial, 3, first_pivot, first_rest, prime
        )
        for order, polynomial in moments.items()
        if order > 2
    }
    b1_replacement, b2_replacement = finite_ratio_replacements(prime)

    def transform(polynomial):
        answer = substitute_parameter_variable(
            polynomial, 5, b1_replacement, prime
        )
        answer = substitute_parameter_variable(
            answer, 6, b2_replacement, prime
        )
        return rescale_a2_by_b0(answer, prime)

    transformed: dict[int, dict[tuple[int, ...], int]] = {}
    b0_contents: dict[str, int] = {}
    for order, polynomial in eliminated.items():
        transformed_polynomial, b0_content = transform(polynomial)
        transformed[order] = transformed_polynomial
        b0_contents[str(order)] = b0_content
    incidence = transformed[3]
    groups = coefficient_groups(incidence, 2, prime)
    if set(groups) != {0, 1}:
        raise ValueError("A did not vanish or transformed mu_3 is constant")
    constant_coefficient = groups[0]
    linear_coefficient = groups[1]
    later = {
        order: polynomial
        for order, polynomial in transformed.items()
        if order > 3
    }
    reduced_first_pivot, _ = transform(first_pivot)
    reduced_discriminant, _ = transform(
        discriminant_polynomial(prime)
    )
    reduced_minor, _ = transform(
        channel_minor_polynomial("01", prime)
    )
    b0_polynomial = {
        (0, 0, 0, 0, 1, 0, 0): 1,
    }
    open_product = reduced_first_pivot
    open_factors = [
        b0_polynomial,
        reduced_discriminant,
        reduced_minor,
    ]
    if coefficient_branch == "open":
        open_factors.append(linear_coefficient)
    for factor in open_factors:
        open_product = multiply_parameter_polynomials(
            open_product, factor, prime
        )
    open_inverse = active_parameter_polynomial_string(
        open_product,
        prime,
        VARIABLES,
        ACTIVE_PARAMETERS,
        multiply_by_r=True,
        subtract_one=True,
    )
    open_saturation = active_parameter_polynomial_string(
        open_product, prime, VARIABLES, ACTIVE_PARAMETERS
    )
    if coefficient_branch == "open":
        branch_equations = [incidence]
    else:
        branch_equations = [linear_coefficient, constant_coefficient]
    generators = [
        *(
            active_parameter_polynomial_string(
                polynomial, prime, VARIABLES, ACTIVE_PARAMETERS
            )
            for polynomial in branch_equations
        ),
        *(
            active_parameter_polynomial_string(
                polynomial, prime, VARIABLES, ACTIVE_PARAMETERS
            )
            for polynomial in later.values()
        ),
    ]
    if f4sat:
        source_variables = VARIABLES[1:]
        generators.append(open_saturation)
    else:
        source_variables = VARIABLES
        generators.insert(0, open_inverse)
    source = (
        ",".join(source_variables)
        + f"\n{prime}\n"
        + ",\n".join(generators)
        + "\n"
    )
    return source, {
        "field": f"GF({prime})",
        "coefficient_branch": coefficient_branch,
        "linear_incidence_terms": len(incidence),
        "mu3_coefficient_terms": {
            "constant": len(constant_coefficient),
            "linear": len(linear_coefficient),
        },
        "removed_b0_powers": b0_contents,
        "later_moment_terms": {
            str(order): len(polynomial)
            for order, polynomial in later.items()
        },
    }


def characteristic_zero_source(
    orders: tuple[int, ...], f4sat: bool, coefficient_branch: str
) -> tuple[str, dict[str, object]]:
    moments = {
        order: rational_polynomial(polynomial)
        for order, polynomial in exact_moment_polynomials(orders).items()
    }
    first_pivot, first_rest = split_linear_rational_polynomial(
        moments[2], 3
    )
    eliminated = {
        order: rational_linear_substitution_numerator(
            polynomial, 3, first_pivot, first_rest
        )
        for order, polynomial in moments.items()
        if order > 2
    }
    b1_replacement, b2_replacement = rational_ratio_replacements()

    def transform(polynomial):
        answer = substitute_rational_variable(
            polynomial, 5, b1_replacement
        )
        answer = substitute_rational_variable(
            answer, 6, b2_replacement
        )
        return rescale_a2_by_b0(answer)

    transformed: dict[int, RationalParameterPolynomial] = {}
    b0_contents: dict[str, int] = {}
    for order, polynomial in eliminated.items():
        transformed_polynomial, b0_content = transform(polynomial)
        transformed[order] = transformed_polynomial
        b0_contents[str(order)] = b0_content
    incidence = transformed[3]
    groups = rational_coefficient_groups(incidence, 2)
    if set(groups) != {0, 1}:
        raise ValueError("A did not vanish or transformed mu_3 is constant")
    constant_coefficient = groups[0]
    linear_coefficient = groups[1]
    later = {
        order: polynomial
        for order, polynomial in transformed.items()
        if order > 3
    }
    discriminant = rational_polynomial(
        exact_parameter_polynomial(
            discriminant_polynomial(EXACT_LIFT_PRIME)
        )
    )
    minor = rational_polynomial(
        exact_parameter_polynomial(
            channel_minor_polynomial("01", EXACT_LIFT_PRIME)
        )
    )
    reduced_first_pivot, _ = transform(first_pivot)
    reduced_discriminant, _ = transform(discriminant)
    reduced_minor, _ = transform(minor)
    b0_polynomial: RationalParameterPolynomial = {
        (0, 0, 0, 0, 1, 0, 0): Fraction(1),
    }
    open_product = reduced_first_pivot
    open_factors = [
        b0_polynomial,
        reduced_discriminant,
        reduced_minor,
    ]
    if coefficient_branch == "open":
        open_factors.append(linear_coefficient)
    for factor in open_factors:
        open_product = multiply_rational_polynomials(
            open_product, factor
        )
    open_inverse = active_rational_polynomial_string(
        open_product,
        VARIABLES,
        ACTIVE_PARAMETERS,
        multiply_by_r=True,
        subtract_one=True,
    )
    open_saturation = active_rational_polynomial_string(
        open_product, VARIABLES, ACTIVE_PARAMETERS
    )
    if coefficient_branch == "open":
        branch_equations = [incidence]
    else:
        branch_equations = [linear_coefficient, constant_coefficient]
    generators = [
        *(
            active_rational_polynomial_string(
                polynomial, VARIABLES, ACTIVE_PARAMETERS
            )
            for polynomial in branch_equations
        ),
        *(
            active_rational_polynomial_string(
                polynomial, VARIABLES, ACTIVE_PARAMETERS
            )
            for polynomial in later.values()
        ),
    ]
    if f4sat:
        source_variables = VARIABLES[1:]
        generators.append(open_saturation)
    else:
        source_variables = VARIABLES
        generators.insert(0, open_inverse)
    source = (
        ",".join(source_variables)
        + "\n0\n"
        + ",\n".join(generators)
        + "\n"
    )
    return source, {
        "field": "QQ",
        "coefficient_branch": coefficient_branch,
        "linear_incidence_terms": len(incidence),
        "mu3_coefficient_terms": {
            "constant": len(constant_coefficient),
            "linear": len(linear_coefficient),
        },
        "removed_b0_powers": b0_contents,
        "later_moment_terms": {
            str(order): len(polynomial)
            for order, polynomial in later.items()
        },
    }


def as_singular_source(
    msolve_source: str,
    variable_order: tuple[str, ...],
    ordering: str,
) -> str:
    lines = msolve_source.splitlines()
    source_variables = tuple(lines[0].split(","))
    if set(variable_order) != set(source_variables):
        raise ValueError(
            "the Singular variable order must be a permutation of "
            + ",".join(source_variables)
        )
    characteristic = lines[1]
    generators = "\n".join(lines[2:])
    return "\n".join(
        [
            f"ring R={characteristic},({','.join(variable_order)}),{ordering};",
            f"ideal I={generators};",
            "option(redSB);",
            "int started=timer;",
            "ideal G=slimgb(I);",
            'print(\"SIC33_LINEAR_INCIDENCE_BEGIN\");',
            'print(\"elapsed_ticks=\"+string(timer-started));',
            'print(\"basis_size=\"+string(size(G)));',
            'print(\"dimension=\"+string(dim(G)));',
            'print(\"contains_one=\"+string(reduce(1,G)==0));',
            'print(\"first_basis_element=\"+string(G[1]));',
            'print(\"SIC33_LINEAR_INCIDENCE_END\");',
            "$",
        ]
    ) + "\n"


def omit_inverse_localizer(msolve_source: str) -> str:
    lines = msolve_source.splitlines()
    variables = lines[0].split(",")
    if not variables or variables[0] != "r":
        raise ValueError("the source does not begin with an inverse variable")
    generators = "\n".join(lines[2:]).split(",\n")
    return (
        ",".join(variables[1:])
        + f"\n{lines[1]}\n"
        + ",\n".join(generators[1:])
        + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orders", default="2,3,4,5,6,7")
    parser.add_argument("--prime", type=int, default=1_073_741_827)
    parser.add_argument("--characteristic-zero", action="store_true")
    parser.add_argument(
        "--coefficient-branch",
        choices=("open", "boundary"),
        default="open",
    )
    parser.add_argument("--f4sat", action="store_true")
    parser.add_argument("--no-localizer", action="store_true")
    parser.add_argument(
        "--format", choices=("msolve", "singular"), default="msolve"
    )
    parser.add_argument(
        "--singular-order", default="z,r,a1,b0,q,lambda"
    )
    parser.add_argument("--singular-ordering", default="dp")
    parser.add_argument("--emit", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    orders = normalized_orders(args.orders)
    if args.characteristic_zero and args.f4sat:
        parser.error("msolve F4 saturation is restricted to prime fields")
    if args.format == "singular" and args.f4sat:
        parser.error("the Singular source uses the inverse localizer")
    if args.no_localizer and args.f4sat:
        parser.error("--no-localizer and --f4sat are mutually exclusive")
    if args.characteristic_zero:
        source, profiles = characteristic_zero_source(
            orders, args.f4sat, args.coefficient_branch
        )
    else:
        source, profiles = finite_field_source(
            orders, args.prime, args.f4sat, args.coefficient_branch
        )
    if args.no_localizer:
        source = omit_inverse_localizer(source)
    if args.format == "singular":
        source = as_singular_source(
            source,
            tuple(args.singular_order.split(",")),
            args.singular_ordering,
        )
    if args.emit:
        args.emit.parent.mkdir(parents=True, exist_ok=True)
        args.emit.write_text(source)

    record: dict[str, object] = {
        "status": "exploratory",
        "chart": "rank-two cubic Hurwitz 01-minor open",
        "branch": (
            "P1!=0, A=0, B!=0 linear incidence"
            if args.coefficient_branch == "open"
            else "P1!=0, A=B=C=0 degree-drop boundary"
        ),
        "orders": list(orders),
        "localization": [
            "P1",
            "b0",
            *(["B"] if args.coefficient_branch == "open" else []),
            "Delta",
            "M01",
        ],
        "exact_observation": (
            "on A=0 and M01!=0, b0!=0, so A eliminates b2"
        ),
        "f4sat": args.f4sat,
        "localized": not args.no_localizer,
        "format": args.format,
        "source_bytes": len(source.encode()),
        "source_sha256": sha256(source.encode()).hexdigest(),
        "profiles": profiles,
        "interpretation": (
            "generic component source only; no emptiness or lifting claim"
        ),
        "reproduction_command": " ".join(sys.argv),
    }
    print(json.dumps(record, indent=2, sort_keys=True))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(record, indent=2, sort_keys=True) + "\n"
        )


if __name__ == "__main__":
    main()
