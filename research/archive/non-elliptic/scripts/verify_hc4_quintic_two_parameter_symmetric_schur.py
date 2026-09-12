#!/usr/bin/env python3
"""Verify generic Schur rigidity on a two-parameter sextic surface.

The sextic is

  (x^6+y^6+z^6)/30
  + mu*x^2*y^2*z^2
  + nu*sum_{i != j} x_i^4*x_j^2.

On nu != 0, six boundary coefficients solve the quadratic Schur quotient.
After clearing the resulting nu^2 denominator, 114 intrinsic equations
remain in the fifteen quartic coefficients.  Over Q(mu,nu), their exact
Groebner basis makes every quartic coefficient nilpotent of exponent three.

This proves rigidity at the generic point of the parameter surface.  It
does not claim that the specialization locus inside nu != 0 is empty.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess

import sympy as sp


x, y, z, mu, nu = sp.symbols("x y z mu nu")
quartic_coefficients = sp.symbols("s0:15")

quartic_monomials = [
    x**i * y**j * z ** (4 - i - j)
    for i in range(5)
    for j in range(5 - i)
]
quartic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(
        quartic_coefficients, quartic_monomials
    )
)

mixed_42 = sum(
    left**4 * right**2
    for left in (x, y, z)
    for right in (x, y, z)
    if left != right
)
h6 = (
    (x**6 + y**6 + z**6) / 30
    + mu * x**2 * y**2 * z**2
    + nu * mixed_42
)
hessian = sp.hessian(h6, (x, y, z))
hessian_determinant = sp.expand(hessian.det())
gradient = sp.Matrix(
    [sp.diff(quartic, variable) for variable in (x, y, z)]
)
schur_numerator = sp.expand(
    (gradient.T * hessian.adjugate() * gradient)[0]
)

quadratic_monomials = (
    x**2,
    y**2,
    z**2,
    x * y,
    x * z,
    y * z,
)
degree_14_monomials = [
    x**i * y**j * z ** (14 - i - j)
    for i in range(15)
    for j in range(15 - i)
]
quotient_columns = tuple(
    sp.Poly(
        hessian_determinant * quadratic_monomial,
        x,
        y,
        z,
    )
    for quadratic_monomial in quadratic_monomials
)
quotient_matrix = sp.Matrix(
    [
        [
            quotient_column.coeff_monomial(degree_14_monomial)
            for quotient_column in quotient_columns
        ]
        for degree_14_monomial in degree_14_monomials
    ]
)
schur_numerator_polynomial = sp.Poly(schur_numerator, x, y, z)
numerator_vector = sp.Matrix(
    [
        schur_numerator_polynomial.coeff_monomial(monomial)
        for monomial in degree_14_monomials
    ]
)

# These are z^14, y*z^13, y^2*z^12, x*z^13, x*y*z^12,
# and x^2*z^12.  Their quotient matrix has determinant 4096*nu^12.
pivot_rows = (0, 1, 2, 15, 16, 29)
pivot_matrix = quotient_matrix[list(pivot_rows), :]
assert sp.factor(pivot_matrix.det()) == 4096 * nu**12

s = quartic_coefficients

# To avoid repeated rational-matrix simplification, store 2*nu^2 times
# the six exact quotient solutions.
quotient_numerators = (
    -mu * s[1] ** 2
    + 640 * nu**3 * s[0] ** 2
    - 96 * nu**2 * s[0] * s[9]
    + 2 * nu**2 * s[5] ** 2
    + 2 * nu * s[1] * s[10]
    + 6 * nu * s[12] * s[5]
    - 6 * nu * s[5] ** 2
    + nu * s[6] ** 2
    + 4 * nu * s[9] ** 2,
    -mu * s[5] ** 2
    + 640 * nu**3 * s[0] ** 2
    - 96 * nu**2 * s[0] * s[2]
    + 2 * nu**2 * s[1] ** 2
    - 6 * nu * s[1] ** 2
    + 6 * nu * s[1] * s[3]
    + 4 * nu * s[2] ** 2
    + 2 * nu * s[5] * s[7]
    + nu * s[6] ** 2,
    nu * (32 * nu * s[0] ** 2 + s[1] ** 2 + s[5] ** 2),
    -4
    * (
        mu * s[1] * s[5]
        + 24 * nu**2 * s[0] * s[6]
        - nu**2 * s[1] * s[5]
        - nu * s[1] * s[7]
        - nu * s[10] * s[5]
        - nu * s[2] * s[6]
        - nu * s[6] * s[9]
    ),
    -2
    * nu
    * (8 * nu * s[0] * s[5] - s[1] * s[6] - 2 * s[5] * s[9]),
    -2
    * nu
    * (8 * nu * s[0] * s[1] - 2 * s[1] * s[2] - s[5] * s[6]),
)

for pivot_position, row in enumerate(pivot_rows):
    pivot_remainder = sp.expand(
        2 * nu**2 * numerator_vector[row]
        - sum(
            quotient_matrix[row, column]
            * quotient_numerators[column]
            for column in range(6)
        )
    )
    assert pivot_remainder == 0, pivot_position

intrinsic_equations = []
for row in range(len(degree_14_monomials)):
    equation = sp.expand(
        2 * nu**2 * numerator_vector[row]
        - sum(
            quotient_matrix[row, column]
            * quotient_numerators[column]
            for column in range(6)
        )
    )
    if equation != 0:
        intrinsic_equations.append(equation)
assert len(intrinsic_equations) == 114


def singular_expression(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


singular = shutil.which("Singular")
if singular is None:
    raise RuntimeError("Singular is required for the exact generic check")

parser = argparse.ArgumentParser()
parser.add_argument(
    "--extract-denominators",
    action="store_true",
    help=(
        "print the factored coefficient denominator of the generic basis "
        "and its lift transformation"
    ),
)
parser.add_argument(
    "--extract-basis-denominators",
    action="store_true",
    help=(
        "print the factored coefficient denominator of the generic basis "
        "without computing a lift transformation"
    ),
)
parser.add_argument(
    "--basis-profile",
    action="store_true",
    help="print the generic quotient dimension and degree profile",
)
parser.add_argument(
    "--fitting-profile",
    action="store_true",
    help="compute the degree-two fraction-free Fitting pivot",
)
parser.add_argument(
    "--cube-torsion-profile",
    action="store_true",
    help=(
        "compute the parameter torsion generated by the fifteen coefficient "
        "cubes in the degree-three Schur quotient"
    ),
)
parser.add_argument(
    "--cube-index",
    type=int,
    default=None,
    help=(
        "with --cube-torsion-profile, restrict to one coefficient cube "
        "(0 through 14)"
    ),
)
parser.add_argument(
    "--cube-torsion-stage",
    choices=(
        "finite-field",
        "fiber",
        "specialize",
        "generic",
        "annihilator",
        "relations",
        "fitting",
        "all",
    ),
    default="annihilator",
    help="largest exact cube-torsion calculation to run",
)
parser.add_argument(
    "--cube-prime",
    type=int,
    default=11,
    help="prime for the finite-field cube-support scan",
)
parser.add_argument(
    "--cube-mu-value",
    default="0",
    help="rational mu value for the specialized cube-membership check",
)
parser.add_argument(
    "--cube-nu-value",
    default="0",
    help="rational nu value for the specialized cube-membership check",
)
parser.add_argument(
    "--fourth-power-profile",
    action="store_true",
    help=(
        "scan the full fifteen-coefficient reduced Schur fiber by testing "
        "all coefficient fourth powers in degree four"
    ),
)
parser.add_argument(
    "--fourth-prime",
    type=int,
    default=11,
    help="prime for the finite-field fourth-power scan on nu != 0",
)
parser.add_argument(
    "--fourth-stage",
    choices=("scan", "annihilator"),
    default="scan",
    help="pointwise scan or symbolic parameter-annihilator calculation",
)
parser.add_argument(
    "--fourth-timeout",
    type=int,
    default=900,
    help="Singular timeout for --fourth-stage annihilator",
)
arguments = parser.parse_args()

if arguments.fourth_power_profile:
    from itertools import combinations_with_replacement

    prime = arguments.fourth_prime
    if prime < 2 or not sp.isprime(prime):
        parser.error("--fourth-prime must be prime")

    def exponent_vector(
        indices: tuple[int, ...],
    ) -> tuple[int, ...]:
        return tuple(
            indices.count(coefficient_index)
            for coefficient_index in range(len(s))
        )

    coefficient_characters = tuple(
        tuple(
            exponent % 2
            for exponent in sp.Poly(monomial, x, y, z).monoms()[0]
        )
        for monomial in quartic_monomials
    )

    def coefficient_monomial_character(
        exponent: tuple[int, ...],
    ) -> tuple[int, int, int]:
        return tuple(
            sum(
                exponent[index]
                * coefficient_characters[index][coordinate]
                for index in range(len(s))
            )
            % 2
            for coordinate in range(3)
        )

    zero_character = (0, 0, 0)
    quadratic_exponents = tuple(
        exponent_vector(indices)
        for indices in combinations_with_replacement(
            range(len(s)), 2
        )
    )
    fourth_exponents = tuple(
        exponent_vector(indices)
        for indices in combinations_with_replacement(
            range(len(s)), 4
        )
        if coefficient_monomial_character(
            exponent_vector(indices)
        )
        == zero_character
    )
    fourth_index = {
        exponent: index for index, exponent in enumerate(fourth_exponents)
    }

    primitive_equations = []
    row_contents = []
    for equation in intrinsic_equations:
        polynomial = sp.Poly(equation, *s)
        content = sp.Integer(0)
        for coefficient in polynomial.coeffs():
            content = sp.gcd(content, coefficient)
        content = sp.factor(content)
        row_contents.append(content)
        primitive_equations.append(sp.cancel(equation / content))

    content_counts = {
        content: row_contents.count(content)
        for content in set(row_contents)
    }
    print(
        "FOURTH_POWER_ROW_CONTENTS:",
        ", ".join(
            f"{sp.sstr(content)}:{multiplicity}"
            for content, multiplicity in sorted(
                content_counts.items(), key=lambda item: sp.sstr(item[0])
            )
        ),
        flush=True,
    )

    multiplication_sparse_columns = []
    for equation in primitive_equations:
        polynomial = sp.Poly(equation, *s)
        coefficient_by_exponent = dict(polynomial.terms())
        equation_characters = {
            coefficient_monomial_character(exponent)
            for exponent in coefficient_by_exponent
        }
        assert len(equation_characters) == 1
        equation_character = next(iter(equation_characters))
        for multiplier_exponent in quadratic_exponents:
            if (
                coefficient_monomial_character(multiplier_exponent)
                != equation_character
            ):
                continue
            entries = []
            for exponent, coefficient in coefficient_by_exponent.items():
                fourth_exponent = tuple(
                    exponent[index] + multiplier_exponent[index]
                    for index in range(len(s))
                )
                entries.append(
                    (fourth_index[fourth_exponent], coefficient)
                )
            multiplication_sparse_columns.append(tuple(entries))

    fourth_power_components = tuple(
        fourth_index[
            tuple(
                4 if index == coefficient_index else 0
                for index in range(len(s))
            )
        ]
        for coefficient_index in range(len(s))
    )
    fourth_power_component_set = set(fourth_power_components)
    multiplication_sparse_columns.sort(
        key=lambda entries: not any(
            row in fourth_power_component_set for row, _ in entries
        )
    )
    print(
        "FOURTH_POWER_CHARACTER_BLOCK:",
        zero_character,
        len(fourth_exponents),
        len(multiplication_sparse_columns),
        flush=True,
    )

    if arguments.fourth_stage == "annihilator":
        multiplication_columns = tuple(
            "+".join(
                f"({singular_expression(coefficient)})*gen({row + 1})"
                for row, coefficient in entries
            )
            for entries in multiplication_sparse_columns
        )
        fourth_power_columns = tuple(
            f"gen({component + 1})"
            for component in fourth_power_components
        )
        annihilator_program = f"""
LIB "primdec.lib";
ring fourth_parameter_ring={prime},(mu,nu),dp;
module multiplication_map={",".join(multiplication_columns)};
module fourth_power_targets={",".join(fourth_power_columns)};
multiplication_map=slimgb(multiplication_map);
print(
  "FOURTH_POWER_INTEGRAL_BASIS "
  +string(ncols(multiplication_map))
);
ideal fourth_annihilator=std(
  quotient(multiplication_map,fourth_power_targets)
);
print("FOURTH_POWER_ANNIHILATOR_BEGIN");
fourth_annihilator;
print("FOURTH_POWER_ANNIHILATOR_END");
ideal fourth_radical=std(radical(fourth_annihilator));
print("FOURTH_POWER_RADICAL_BEGIN");
fourth_radical;
print("FOURTH_POWER_RADICAL_END");
list fourth_components=minAssGTZ(fourth_radical);
print("FOURTH_POWER_COMPONENTS_BEGIN");
fourth_components;
print("FOURTH_POWER_COMPONENTS_END");
"""
        try:
            completed = subprocess.run(
                [singular, "-q"],
                input=annihilator_program,
                text=True,
                capture_output=True,
                check=True,
                timeout=arguments.fourth_timeout,
            )
        except subprocess.TimeoutExpired:
            print(
                "TIMEOUT: symbolic fourth-power annihilator over F_"
                f"{prime} produced no certificate in "
                f"{arguments.fourth_timeout} seconds"
            )
            raise SystemExit(2)
        if completed.stderr.strip() or "?" in completed.stdout:
            raise RuntimeError(
                "Singular fourth-power annihilator failed:\n"
                f"{completed.stdout}\n{completed.stderr}"
            )
        print(completed.stdout.strip())
        raise SystemExit

    coefficient_evaluators = {}
    for column in multiplication_sparse_columns:
        for _, coefficient in column:
            if coefficient in coefficient_evaluators:
                continue
            terms = []
            for exponent, scalar in sp.Poly(
                coefficient, mu, nu
            ).terms():
                numerator, denominator = map(
                    int, scalar.as_numer_denom()
                )
                terms.append(
                    (
                        exponent,
                        numerator
                        * pow(denominator, -1, prime)
                        % prime,
                    )
                )
            coefficient_evaluators[coefficient] = tuple(terms)

    def evaluate_coefficient(
        coefficient: sp.Expr, mu_value: int, nu_value: int
    ) -> int:
        return sum(
            scalar
            * pow(mu_value, exponent[0], prime)
            * pow(nu_value, exponent[1], prime)
            for exponent, scalar in coefficient_evaluators[coefficient]
        ) % prime

    def reduce_column(
        vector: dict[int, int],
        basis: dict[int, dict[int, int]],
    ) -> dict[int, int]:
        while vector:
            pivot = min(vector)
            value = vector[pivot]
            if pivot not in basis:
                inverse = pow(value, -1, prime)
                return {
                    row: coefficient * inverse % prime
                    for row, coefficient in vector.items()
                    if coefficient % prime
                }
            pivot_vector = basis[pivot]
            for row, coefficient in pivot_vector.items():
                updated = (
                    vector.get(row, 0) - value * coefficient
                ) % prime
                if updated:
                    vector[row] = updated
                else:
                    vector.pop(row, None)
        return {}

    exceptional_points = []
    processed_column_histogram = {}
    certified_empty_points = 0
    for mu_value in range(prime):
        for nu_value in range(1, prime):
            evaluated_coefficients = {
                coefficient: evaluate_coefficient(
                    coefficient, mu_value, nu_value
                )
                for coefficient in coefficient_evaluators
            }
            basis = {}
            certified_empty = False
            processed_columns = 0
            for processed_columns, column in enumerate(
                multiplication_sparse_columns, start=1
            ):
                vector = {
                    row: evaluated_coefficients[coefficient]
                    for row, coefficient in column
                    if evaluated_coefficients[coefficient]
                }
                reduced = reduce_column(vector, basis)
                if reduced:
                    basis[min(reduced)] = reduced
                if (
                    processed_columns % 64 == 0
                    and all(
                        not reduce_column({component: 1}, basis)
                        for component in fourth_power_components
                    )
                ):
                    certified_empty = True
                    break
            processed_column_histogram[processed_columns] = (
                processed_column_histogram.get(processed_columns, 0) + 1
            )
            if certified_empty:
                certified_empty_points += 1
                continue
            failing_targets = tuple(
                coefficient_index
                for coefficient_index, component in enumerate(
                    fourth_power_components
                )
                if reduce_column({component: 1}, basis)
            )
            if failing_targets:
                exceptional_points.append(
                    (mu_value, nu_value, failing_targets, len(basis))
                )
    print(
        "FOURTH_POWER_CERTIFIED_EMPTY_POINTS:",
        certified_empty_points,
    )
    print(
        "FOURTH_POWER_PROCESSED_COLUMN_HISTOGRAM:",
        processed_column_histogram,
    )
    print("FOURTH_POWER_EXCEPTIONAL_POINTS:", exceptional_points)
    print(
        "SCOPE: exhaustive F_"
        f"{prime} scan on nu!=0; all fourth powers in the ideal certify "
        "an empty reduced projective fiber"
    )
    raise SystemExit

if arguments.cube_torsion_profile:
    from itertools import combinations_with_replacement

    quadratic_exponents = tuple(
        tuple(
            indices.count(coefficient_index)
            for coefficient_index in range(len(s))
        )
        for indices in combinations_with_replacement(
            range(len(s)), 2
        )
    )
    all_cubic_exponents = tuple(
        tuple(
            indices.count(coefficient_index)
            for coefficient_index in range(len(s))
        )
        for indices in combinations_with_replacement(
            range(len(s)), 3
        )
    )
    assert len(quadratic_exponents) == 120
    assert len(all_cubic_exponents) == 680

    coefficient_characters = tuple(
        tuple(exponent % 2 for exponent in sp.Poly(
            monomial, x, y, z
        ).monoms()[0])
        for monomial in quartic_monomials
    )

    def coefficient_monomial_character(
        exponent: tuple[int, ...],
    ) -> tuple[int, int, int]:
        return tuple(
            sum(
                exponent[index] * coefficient_characters[index][coordinate]
                for index in range(len(s))
            )
            % 2
            for coordinate in range(3)
        )

    if arguments.cube_index is not None:
        if not 0 <= arguments.cube_index < len(s):
            parser.error("--cube-index must be between 0 and 14")
        target_character = coefficient_characters[arguments.cube_index]
        cubic_exponents = tuple(
            exponent
            for exponent in all_cubic_exponents
            if coefficient_monomial_character(exponent) == target_character
        )
    else:
        target_character = None
        cubic_exponents = all_cubic_exponents
    cubic_index = {
        exponent: index + 1
        for index, exponent in enumerate(cubic_exponents)
    }

    primitive_equations = []
    row_contents = []
    for equation in intrinsic_equations:
        polynomial = sp.Poly(equation, *s)
        content = sp.Integer(0)
        for coefficient in polynomial.coeffs():
            content = sp.gcd(content, coefficient)
        content = sp.factor(content)
        row_contents.append(content)
        primitive_equations.append(sp.cancel(equation / content))

    content_counts = {
        content: row_contents.count(content)
        for content in set(row_contents)
    }
    print(
        "CUBE_TORSION_ROW_CONTENTS:",
        ", ".join(
            f"{sp.sstr(content)}:{multiplicity}"
            for content, multiplicity in sorted(
                content_counts.items(), key=lambda item: sp.sstr(item[0])
            )
        ),
    )

    multiplication_sparse_columns = []
    for equation in primitive_equations:
        polynomial = sp.Poly(equation, *s)
        coefficient_by_exponent = dict(polynomial.terms())
        equation_characters = {
            coefficient_monomial_character(exponent)
            for exponent in coefficient_by_exponent
        }
        assert len(equation_characters) == 1
        equation_character = next(iter(equation_characters))
        for multiplier_index in range(len(s)):
            column_character = tuple(
                equation_character[coordinate]
                ^ coefficient_characters[multiplier_index][coordinate]
                for coordinate in range(3)
            )
            if (
                target_character is not None
                and column_character != target_character
            ):
                continue
            entries = []
            for exponent, coefficient in coefficient_by_exponent.items():
                cubic_exponent = list(exponent)
                cubic_exponent[multiplier_index] += 1
                component = cubic_index[tuple(cubic_exponent)]
                entries.append((component - 1, coefficient))
            multiplication_sparse_columns.append(tuple(entries))
    multiplication_columns = tuple(
        "+".join(
            f"({singular_expression(coefficient)})*gen({component + 1})"
            for component, coefficient in entries
        )
        for entries in multiplication_sparse_columns
    )
    if target_character is None:
        assert len(multiplication_columns) == 1710
    print(
        "CUBE_TORSION_CHARACTER_BLOCK:",
        target_character,
        len(cubic_exponents),
        len(multiplication_columns),
    )

    active_cube_indices = (
        (arguments.cube_index,)
        if arguments.cube_index is not None
        else tuple(range(len(s)))
    )
    cube_components = tuple(
        cubic_index[
            tuple(
                3 if index == coefficient_index else 0
                for index in range(len(s))
            )
        ]
        for coefficient_index in active_cube_indices
    )
    cube_columns = tuple(
        f"gen({component})" for component in cube_components
    )
    if arguments.cube_index is not None:
        selected_cube_columns = cube_columns
    else:
        selected_cube_columns = cube_columns

    if arguments.cube_torsion_stage == "finite-field":
        prime = arguments.cube_prime
        if prime < 2 or not sp.isprime(prime):
            parser.error("--cube-prime must be prime")

        coefficient_evaluators = {}
        for column in multiplication_sparse_columns:
            for _, coefficient in column:
                if coefficient not in coefficient_evaluators:
                    terms = []
                    for exponent, scalar in sp.Poly(
                        coefficient, mu, nu
                    ).terms():
                        numerator, denominator = map(
                            int, scalar.as_numer_denom()
                        )
                        terms.append(
                            (
                                exponent,
                                numerator
                                * pow(denominator, -1, prime)
                                % prime,
                            )
                        )
                    coefficient_evaluators[coefficient] = tuple(terms)

        def evaluate_coefficient(
            coefficient: sp.Expr, mu_value: int, nu_value: int
        ) -> int:
            return sum(
                scalar
                * pow(mu_value, exponent[0], prime)
                * pow(nu_value, exponent[1], prime)
                for exponent, scalar in coefficient_evaluators[coefficient]
            ) % prime

        def reduce_column(
            vector: dict[int, int],
            basis: dict[int, dict[int, int]],
        ) -> dict[int, int]:
            while vector:
                pivot = min(vector)
                value = vector[pivot]
                if pivot not in basis:
                    inverse = pow(value, -1, prime)
                    return {
                        row: coefficient * inverse % prime
                        for row, coefficient in vector.items()
                        if coefficient % prime
                    }
                pivot_vector = basis[pivot]
                for row, coefficient in pivot_vector.items():
                    updated = (
                        vector.get(row, 0) - value * coefficient
                    ) % prime
                    if updated:
                        vector[row] = updated
                    else:
                        vector.pop(row, None)
            return {}

        target_components = tuple(component - 1 for component in cube_components)
        exceptional_points = []
        rank_histogram = {}
        for mu_value in range(prime):
            for nu_value in range(prime):
                basis = {}
                for column in multiplication_sparse_columns:
                    vector = {
                        row: value
                        for row, coefficient in column
                        if (
                            value := evaluate_coefficient(
                                coefficient, mu_value, nu_value
                            )
                        )
                    }
                    reduced = reduce_column(vector, basis)
                    if reduced:
                        basis[min(reduced)] = reduced
                rank = len(basis)
                rank_histogram[rank] = rank_histogram.get(rank, 0) + 1
                failing_targets = []
                for target_index, component in zip(
                    active_cube_indices, target_components
                ):
                    remainder = reduce_column(
                        {component: 1}, basis
                    )
                    if remainder:
                        failing_targets.append(target_index)
                if failing_targets:
                    exceptional_points.append(
                        (mu_value, nu_value, tuple(failing_targets), rank)
                    )
        print("FINITE_FIELD_RANK_HISTOGRAM:", rank_histogram)
        print("FINITE_FIELD_EXCEPTIONAL_POINTS:", exceptional_points)
        raise SystemExit

    if arguments.cube_torsion_stage == "specialize":
        mu_value = sp.Rational(arguments.cube_mu_value)
        nu_value = sp.Rational(arguments.cube_nu_value)
        specialized_columns = tuple(
            "+".join(
                f"({singular_expression(coefficient.subs({
                    mu: mu_value, nu: nu_value
                }))})*gen({component + 1})"
                for component, coefficient in entries
                if coefficient.subs({mu: mu_value, nu: nu_value}) != 0
            )
            or "0"
            for entries in multiplication_sparse_columns
        )
        specialization_program = f"""
ring rational_specialization=0,t,dp;
module multiplication_map={",".join(specialized_columns)};
module cube_targets={",".join(selected_cube_columns)};
module specialized_basis=slimgb(multiplication_map);
module specialized_remainder=reduce(
  cube_targets,specialized_basis
);
print("SPECIALIZED_MULTIPLICATION_RANK "
  +string(ncols(specialized_basis)));
print("SPECIALIZED_CUBE_REMAINDER_BEGIN");
specialized_remainder;
print("SPECIALIZED_CUBE_REMAINDER_END");
"""
        completed = subprocess.run(
            [singular, "-q"],
            input=specialization_program,
            text=True,
            capture_output=True,
            check=True,
            timeout=300,
        )
        if completed.stderr.strip() or "?" in completed.stdout:
            raise RuntimeError(
                "Singular specialized cube calculation failed:\n"
                f"{completed.stdout}\n{completed.stderr}"
            )
        print(completed.stdout.strip())
        raise SystemExit

    if arguments.cube_torsion_stage == "fiber":
        mu_value = sp.Rational(arguments.cube_mu_value)
        nu_value = sp.Rational(arguments.cube_nu_value)
        specialized_equations = []
        for equation in primitive_equations:
            specialized_polynomial = sp.Poly(
                equation.subs({mu: mu_value, nu: nu_value}),
                *s,
                domain=sp.QQ,
            )
            _, cleared_polynomial = specialized_polynomial.clear_denoms()
            specialized_equations.append(
                singular_expression(cleared_polynomial.as_expr())
            )
        fiber_program = f"""
ring coefficient_fiber=0,({",".join(map(str, s))}),dp;
option(redSB);
ideal fiber_ideal={",".join(specialized_equations)};
ideal fiber_basis=slimgb(fiber_ideal);
print("FIBER_BASIS_SIZE "+string(size(fiber_basis)));
print("FIBER_VECTOR_SPACE_DIMENSION "+string(vdim(fiber_basis)));
ideal coefficient_maximal={",".join(map(str, s))};
int coefficient_index;
for (
  coefficient_index=1;
  coefficient_index<=size(coefficient_maximal);
  coefficient_index++
)
{{
  print(
    "FIBER_CUBE "
    +string(coefficient_index-1)+" "
    +string(
      reduce(
        coefficient_maximal[coefficient_index]^3,
        fiber_basis
      )==0
    )
  );
  print(
    "FIBER_FOURTH_POWER "
    +string(coefficient_index-1)+" "
    +string(
      reduce(
        coefficient_maximal[coefficient_index]^4,
        fiber_basis
      )==0
    )
  );
}}
"""
        completed = subprocess.run(
            [singular, "-q"],
            input=fiber_program,
            text=True,
            capture_output=True,
            check=True,
            timeout=300,
        )
        if completed.stderr.strip() or "?" in completed.stdout:
            raise RuntimeError(
                "Singular specialized fiber calculation failed:\n"
                f"{completed.stdout}\n{completed.stderr}"
            )
        print(completed.stdout.strip())
        raise SystemExit

    relation_program = ""
    if arguments.cube_torsion_stage in ("relations", "fitting", "all"):
        relation_program = """
module cube_relations=modulo(cube_targets,multiplication_map);
print("CUBE_RELATION_COLUMNS "+string(ncols(cube_relations)));
"""
    fitting_program = ""
    if arguments.cube_torsion_stage in ("fitting", "all"):
        fitting_program = """
ideal cube_fitting_zero=fitting(cube_relations,0);
print("CUBE_FITTING_ZERO_BEGIN");
cube_fitting_zero;
print("CUBE_FITTING_ZERO_END");
"""
    associated_prime_program = ""
    if arguments.cube_torsion_stage in ("annihilator", "all"):
        associated_prime_program = """
list associated_primes=minAssGTZ(cube_annihilator);
print("CUBE_ASSOCIATED_PRIMES_BEGIN");
associated_primes;
print("CUBE_ASSOCIATED_PRIMES_END");
"""
    if arguments.cube_torsion_stage == "generic":
        torsion_program = f"""
ring parameter_function_field=(0,mu,nu),t,dp;
module multiplication_map={",".join(multiplication_columns)};
module cube_targets={",".join(selected_cube_columns)};
module generic_basis=slimgb(multiplication_map);
module generic_remainder=reduce(cube_targets,generic_basis);
print("GENERIC_MULTIPLICATION_BASIS "
  +string(ncols(generic_basis)));
print("GENERIC_CUBE_REMAINDER "
  +string(ncols(generic_remainder)));
matrix cube_certificate=lift(multiplication_map,cube_targets);
print("GENERIC_CERTIFICATE_SHAPE "
  +string(nrows(cube_certificate))+" "
  +string(ncols(cube_certificate)));
int certificate_row;
int certificate_column;
number certificate_denominator;
for (
  certificate_row=1;
  certificate_row<=nrows(cube_certificate);
  certificate_row++
)
{{
  for (
    certificate_column=1;
    certificate_column<=ncols(cube_certificate);
    certificate_column++
  )
  {{
    if (cube_certificate[certificate_row,certificate_column]!=0)
    {{
      certificate_denominator=denominator(
        leadcoef(cube_certificate[
          certificate_row,certificate_column
        ])
      );
      print("CUBE_CERTIFICATE_DENOMINATOR "
        +string(certificate_denominator));
    }}
  }}
}}
"""
    else:
        torsion_program = f"""
LIB "homolog.lib";
LIB "primdec.lib";
ring parameter_ring=0,(mu,nu),dp;
module multiplication_map={",".join(multiplication_columns)};
module cube_targets={",".join(selected_cube_columns)};
multiplication_map=slimgb(multiplication_map);
print("INTEGRAL_MULTIPLICATION_BASIS "
  +string(ncols(multiplication_map)));
ideal cube_annihilator=std(
  quotient(multiplication_map,cube_targets)
);
print("CUBE_ANNIHILATOR_BEGIN");
cube_annihilator;
print("CUBE_ANNIHILATOR_END");
{relation_program}
{fitting_program}
{associated_prime_program}
"""
    completed = subprocess.run(
        [singular, "-q"],
        input=torsion_program,
        text=True,
        capture_output=True,
        check=True,
        timeout=900,
    )
    if completed.stderr.strip() or "?" in completed.stdout:
        raise RuntimeError(
            "Singular cube-torsion calculation failed:\n"
            f"{completed.stdout}\n{completed.stderr}"
        )
    print(completed.stdout.strip())
    raise SystemExit

if arguments.fitting_profile:
    from sympy.polys.matrices import DomainMatrix

    quadratic_coefficient_monomials = [
        s[left] * s[right]
        for left in range(len(s))
        for right in range(left, len(s))
    ]
    quadratic_coefficient_matrix = sp.Matrix(
        [
            [
                sp.Poly(
                    equation,
                    *s,
                ).coeff_monomial(monomial)
                for monomial in quadratic_coefficient_monomials
            ]
            for equation in intrinsic_equations
        ]
    )
    row_contents = []
    for row in range(quadratic_coefficient_matrix.rows):
        content = sp.Integer(0)
        for entry in quadratic_coefficient_matrix.row(row):
            content = sp.gcd(content, entry)
        content = sp.factor(content)
        row_contents.append(content)
        if content != 1:
            quadratic_coefficient_matrix.row_op(
                row, lambda entry, _: sp.cancel(entry / content)
            )
    content_counts = {
        content: row_contents.count(content)
        for content in set(row_contents)
    }
    print(
        "ROW_CONTENTS:",
        ", ".join(
            f"{sp.sstr(content)}:{multiplicity}"
            for content, multiplicity in sorted(
                content_counts.items(), key=lambda item: sp.sstr(item[0])
            )
        ),
    )
    polynomial_domain = sp.QQ.poly_ring(mu, nu)
    fitting_matrix = DomainMatrix.from_Matrix(
        quadratic_coefficient_matrix
    ).convert_to(polynomial_domain)
    fitting_columns = ()
    fitting_rows = ()
    fitting_point = None
    for candidate_mu, candidate_nu in (
        (1, 1),
        (1, 2),
        (2, 1),
        (2, 3),
        (0, 1),
    ):
        numerical_matrix = DomainMatrix.from_Matrix(
            quadratic_coefficient_matrix.subs(
                {mu: candidate_mu, nu: candidate_nu}
            )
        ).convert_to(sp.QQ)
        _, candidate_columns = numerical_matrix.rref()
        print(
            "SPECIALIZED_RANK:",
            candidate_mu,
            candidate_nu,
            len(candidate_columns),
        )
        if len(candidate_columns) > len(fitting_columns):
            fitting_columns = candidate_columns
            fitting_point = (candidate_mu, candidate_nu)
            selected_columns = numerical_matrix.extract(
                range(numerical_matrix.shape[0]),
                candidate_columns,
            )
            _, fitting_rows = selected_columns.transpose().rref()
    assert len(fitting_columns) == 99
    assert len(fitting_rows) == 99
    fitting_minor = fitting_matrix.extract(fitting_rows, fitting_columns)
    fitting_determinant = fitting_minor.det().as_expr()
    print(
        "DEGREE_TWO_FITTING:",
        len(fitting_columns),
        sp.factor(fitting_determinant),
    )
    raise SystemExit

denominator_program = ""
standard_basis_program = "ideal G=slimgb(I);"
if arguments.extract_denominators:
    standard_basis_program = """
matrix transformation;
ideal G=liftstd(I,transformation,"slimgb");
"""
    all_quartic_variables = "*".join(map(str, quartic_coefficients))
    denominator_program = f"""
matrix coefficient_matrix;
int coefficient_index;
number coefficient_denominator;
for (i=1;i<=size(G);i++)
{{
  coefficient_matrix=coef(G[i],{all_quartic_variables});
  for (
    coefficient_index=1;
    coefficient_index<=ncols(coefficient_matrix);
    coefficient_index++
  )
  {{
    coefficient_denominator=denominator(
      leadcoef(coefficient_matrix[2,coefficient_index])
    );
    print("DENOMINATOR "+string(coefficient_denominator));
  }}
}}
int transformation_row;
int transformation_column;
poly transformation_entry;
for (
  transformation_row=1;
  transformation_row<=nrows(transformation);
  transformation_row++
)
{{
  for (
    transformation_column=1;
    transformation_column<=ncols(transformation);
    transformation_column++
  )
  {{
    transformation_entry=transformation[
      transformation_row,transformation_column
    ];
    if (transformation_entry!=0)
    {{
      coefficient_matrix=coef(
        transformation_entry,{all_quartic_variables}
      );
      for (
        coefficient_index=1;
        coefficient_index<=ncols(coefficient_matrix);
        coefficient_index++
      )
      {{
        coefficient_denominator=denominator(
          leadcoef(coefficient_matrix[2,coefficient_index])
        );
        print("DENOMINATOR "+string(coefficient_denominator));
      }}
    }}
  }}
}}
"""
elif arguments.extract_basis_denominators:
    denominator_program = """
poly cleared_basis_element;
for (i=1;i<=size(G);i++)
{
  cleared_basis_element=cleardenom(G[i]);
  print("DENOMINATOR "+string(leadcoef(cleared_basis_element)));
}
"""

profile_program = ""
if arguments.basis_profile:
    profile_program = """
ideal standard_monomials=kbase(G);
int maximum_standard_degree=0;
int standard_index;
for (standard_index=1;standard_index<=size(standard_monomials);standard_index++)
{
  if (deg(standard_monomials[standard_index])>maximum_standard_degree)
  {
    maximum_standard_degree=deg(standard_monomials[standard_index]);
  }
}
print(
  "PROFILE "
  +string(vdim(G))+" "
  +string(size(standard_monomials))+" "
  +string(maximum_standard_degree)
);
intvec hilbert_function=hilb(G,1);
print("HILBERT "+string(hilbert_function));
"""

program = f"""
ring rr=(0,mu,nu),({",".join(map(str, quartic_coefficients))}),dp;
option(redSB);
ideal I={",".join(map(singular_expression, intrinsic_equations))};
{standard_basis_program}
ideal M={",".join(map(str, quartic_coefficients))};
ideal GM=std(M);
print(
  "GENERIC "
  +string(size(I))+" "
  +string(size(G))+" "
  +string(size(reduce(I,GM)))
);
int i;
poly cube;
for (i=1;i<=size(M);i++)
{{
  cube=M[i]^3;
  print("CUBE "+string(i)+" "+string(reduce(cube,G)==0));
}}
{denominator_program}
{profile_program}
"""
completed = subprocess.run(
    [singular, "-q"],
    input=program,
    text=True,
    capture_output=True,
    check=True,
    timeout=(
        900
        if (
            arguments.extract_denominators
            or arguments.extract_basis_denominators
        )
        else 180
    ),
)
if completed.stderr.strip():
    raise RuntimeError(completed.stderr)

generic_marker = re.search(
    r"(?m)^GENERIC (\d+) (\d+) (\d+)$", completed.stdout
)
assert generic_marker is not None
assert tuple(map(int, generic_marker.groups())) == (114, 117, 0)
cube_markers = re.findall(
    r"(?m)^CUBE (\d+) ([01])$", completed.stdout
)
assert len(cube_markers) == 15
assert all(success == "1" for _, success in cube_markers)

if arguments.extract_denominators or arguments.extract_basis_denominators:
    denominator_strings = set(
        re.findall(r"(?m)^DENOMINATOR (.+)$", completed.stdout)
    )
    coefficient_denominator = sp.Integer(1)
    for denominator_string in denominator_strings:
        denominator = sp.sympify(
            denominator_string.replace("^", "**"),
            locals={"mu": mu, "nu": nu},
        )
        coefficient_denominator = sp.lcm(
            coefficient_denominator,
            sp.Poly(denominator, mu, nu),
        ).as_expr()
    print(
        "GENERIC_BASIS_DENOMINATOR:",
        sp.factor(coefficient_denominator),
    )
if arguments.basis_profile:
    profile_marker = re.search(
        r"(?m)^PROFILE (\d+) (\d+) (\d+)$", completed.stdout
    )
    assert profile_marker is not None
    print("GENERIC_QUOTIENT_PROFILE:", " ".join(profile_marker.groups()))
    hilbert_marker = re.search(r"(?m)^HILBERT (.+)$", completed.stdout)
    assert hilbert_marker is not None
    print("GENERIC_HILBERT_DATA:", hilbert_marker.group(1))

print("PASS: six quotient pivots have determinant 4096*nu^12")
print("PASS: quotient elimination leaves 114 intrinsic equations")
print("PASS: the generic Groebner basis has 117 elements")
print("PASS: all fifteen quartic coefficient cubes reduce to zero")
print("SCOPE: generic rigidity only; exceptional nu!=0 fibers remain open")
