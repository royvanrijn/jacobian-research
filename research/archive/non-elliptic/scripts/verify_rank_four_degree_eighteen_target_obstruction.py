#!/usr/bin/env python3
"""Exclude target degree at most eighteen for the framed quartic endpoint.

This is the specialized Singular continuation of R4NP1.  After composing
off the intrinsic mu_5 action, a target self-equivalence has the form

    T = id + H V,

where H is the prime ramified-discriminant equation.  Preservation of H
modulo H^2 makes V a logarithmic derivation.  The exact degree-five
logarithmic space has dimension seven.  The endpoint condition leaves four
parameters, and ten exact Jacobian evaluations generate the unit ideal over
QQ.  Hence no polynomial target automorphism of degree at most eighteen can
carry the two marked endpoint fibers.
"""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp


S, P, B, C = sp.symbols("S P B C")
SEED = sp.Integer(-124416)
START = (
    sp.Rational(-1, 144),
    sp.Rational(1, 24),
    sp.Integer(0),
)
TARGET = (
    sp.Rational(-1, 144),
    sp.Rational(1, 96),
    sp.Integer(-9),
)
DELTA = sp.Matrix(
    [target - start for start, target in zip(START, TARGET)]
)


def monomial_exponents(bound: int) -> list[tuple[int, int, int]]:
    """Exponent triples of total degree at most ``bound``."""

    if bound < 0:
        return []
    return [
        (first, second, third)
        for first in range(bound + 1)
        for second in range(bound + 1 - first)
        for third in range(bound + 1 - first - second)
    ]


def logarithmic_system(
    bound: int,
    discriminant: sp.Expr,
) -> tuple[
    sp.Matrix,
    tuple[sp.Symbol, ...],
    sp.Matrix,
    sp.Expr,
    sp.Symbol,
]:
    """Linear system for V(H)-QH-constant=0 at degree ``bound``."""

    vector_monomials = monomial_exponents(bound)
    quotient_monomials = monomial_exponents(bound - 1)
    variable_count = (
        3 * len(vector_monomials)
        + len(quotient_monomials)
        + 1
    )
    unknowns = sp.symbols(f"logarithmic_unknown_0:{variable_count}")

    def polynomial(offset: int, monomials: list[tuple[int, int, int]]) -> sp.Expr:
        return sum(
            unknowns[offset + index] * P**first * B**second * C**third
            for index, (first, second, third) in enumerate(monomials)
        )

    block_size = len(vector_monomials)
    vector = sp.Matrix(
        [
            polynomial(0, vector_monomials),
            polynomial(block_size, vector_monomials),
            polynomial(2 * block_size, vector_monomials),
        ]
    )
    quotient = polynomial(3 * block_size, quotient_monomials)
    constant = unknowns[-1]
    logarithmic_identity = sp.Poly(
        sp.expand(
            vector.dot(
                sp.Matrix(
                    [
                        sp.diff(discriminant, P),
                        sp.diff(discriminant, B),
                        sp.diff(discriminant, C),
                    ]
                )
            )
            - quotient * discriminant
            - constant
        ),
        P,
        B,
        C,
    )
    matrix, _ = sp.linear_eq_to_matrix(
        logarithmic_identity.coeffs(),
        unknowns,
    )
    return matrix, unknowns, vector, quotient, constant


inverse_relation = (
    SEED * P**4 * S**4
    + P * S**3
    + B * S**2
    + S
    - C / 2
)
inverse_discriminant = sp.factor(sp.discriminant(inverse_relation, S))
H = sp.factor(-4 * inverse_discriminant / P**2)
H_polynomial = sp.Poly(H, P, B, C, domain=sp.QQ)
assert H_polynomial.total_degree() == 13
assert len(H_polynomial.terms()) == 16


# No logarithmic correction exists through vector degree three.  Degree
# four has one direction; degree five has the seven-dimensional frontier.
expected_nullities = (0, 0, 0, 0, 1, 7)
degree_five_data = None
for degree_bound, expected_nullity in enumerate(expected_nullities):
    system_data = logarithmic_system(degree_bound, H)
    system_matrix = system_data[0]
    nullity = system_matrix.cols - system_matrix.rank()
    assert nullity == expected_nullity
    if degree_bound == 5:
        degree_five_data = system_data

assert degree_five_data is not None
(
    degree_five_matrix,
    degree_five_unknowns,
    generic_vector,
    generic_quotient,
    generic_constant,
) = degree_five_data
nullspace = degree_five_matrix.nullspace()
assert len(nullspace) == 7

vector_basis: list[sp.Matrix] = []
for basis_column in nullspace:
    substitution = dict(zip(degree_five_unknowns, basis_column))
    vector = sp.Matrix(
        [
            sp.factor(component.subs(substitution))
            for component in generic_vector
        ]
    )
    quotient = sp.factor(generic_quotient.subs(substitution))
    constant = sp.factor(generic_constant.subs(substitution))
    assert constant == 0
    assert sp.factor(
        vector.dot(sp.Matrix([sp.diff(H, variable) for variable in (P, B, C)]))
        - quotient * H
    ) == 0
    vector_basis.append(vector)


# T(START)=TARGET means H(START)*V(START)=DELTA.
start_substitution = dict(zip((P, B, C), START))
start_discriminant = sp.factor(H.subs(start_substitution))
assert start_discriminant == sp.Rational(-1, 16)
evaluation_matrix = sp.Matrix.hstack(
    *[vector.subs(start_substitution) for vector in vector_basis]
)
assert evaluation_matrix.rank() == 3

basis_coefficients = sp.symbols("lambda_0:7")
solution_set = sp.linsolve(
    (evaluation_matrix, DELTA / start_discriminant),
    basis_coefficients,
)
assert solution_set != sp.EmptySet
solution = next(iter(solution_set))
free_parameters = tuple(
    sorted(
        set().union(*(coordinate.free_symbols for coordinate in solution))
        & set(basis_coefficients),
        key=lambda symbol: symbol.name,
    )
)
assert len(free_parameters) == 4

endpoint_vector = sum(
    (
        solution[index] * vector_basis[index]
        for index in range(len(vector_basis))
    ),
    sp.zeros(3, 1),
)
assert all(
    sp.factor(component) == 0
    for component in (
        start_discriminant
        * endpoint_vector.subs(start_substitution)
        - DELTA
    )
)


# If T=id+H*V were a polynomial automorphism preserving H, its determinant
# would be the constant one.  Ten exact point evaluations already have no
# common parameter solution.
gradient_H = sp.Matrix([sp.diff(H, variable) for variable in (P, B, C)])
jacobian_vector = endpoint_vector.jacobian((P, B, C))
evaluation_points = (
    (1, 1, 0),
    (1, 0, 1),
    (1, 1, 1),
    (-1, 1, 0),
    (2, 0, 0),
    (1, -1, 2),
    (2, 1, -1),
    (-1, 2, 1),
    (2, -1, 1),
    (-2, 1, 2),
)
jacobian_equations: list[sp.Expr] = []
for point in evaluation_points:
    substitution = dict(zip((P, B, C), point))
    derivative = (
        sp.eye(3)
        + endpoint_vector.subs(substitution)
        * gradient_H.subs(substitution).T
        + H.subs(substitution) * jacobian_vector.subs(substitution)
    )
    numerator = sp.together(sp.expand(derivative.det() - 1)).as_numer_denom()[0]
    primitive = sp.primitive(
        sp.Poly(
            numerator,
            *free_parameters,
            domain=sp.QQ,
        )
    )[1].as_expr()
    assert primitive != 0
    jacobian_equations.append(primitive)


singular = shutil.which("Singular")
if singular is None:
    raise RuntimeError(
        "Singular is required for the exact degree-eighteen unit-ideal check"
    )


def singular_expression(expression: sp.Expr) -> str:
    """Render an expanded QQ polynomial for Singular."""

    return sp.sstr(sp.expand(expression)).replace("**", "^")


singular_code = (
    f"ring r=0,({','.join(str(parameter) for parameter in free_parameters)}),dp;\n"
    "option(redSB);\n"
    "ideal I="
    + ",\n".join(singular_expression(equation) for equation in jacobian_equations)
    + ";\n"
    "ideal G=std(I);\n"
    "print(size(G));\n"
    "print(reduce(1,G));\n"
)
singular_result = subprocess.run(
    [singular, "-q"],
    input=singular_code,
    text=True,
    capture_output=True,
    check=False,
    timeout=120,
)
assert singular_result.returncode == 0, singular_result.stderr
singular_lines = [
    line.strip()
    for line in singular_result.stdout.splitlines()
    if line.strip()
]
assert singular_lines[-2:] == ["1", "0"], singular_result.stdout


print("PASS: the ramified discriminant has degree thirteen")
print("PASS: logarithmic nullities through degree five are (0,0,0,0,1,7)")
print("PASS: the endpoint condition leaves exactly four parameters")
print("PASS: ten exact Jacobian evaluations generate the unit ideal over QQ")
print("PASS: no endpoint target symmetry has degree at most eighteen")
print("SCOPE: target degree nineteen and above remains open")
