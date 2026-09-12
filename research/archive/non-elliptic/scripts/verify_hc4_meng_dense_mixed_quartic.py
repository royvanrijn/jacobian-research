#!/usr/bin/env python3
"""Exclude dense mixed homogeneous quartics in the Meng HC(4) chart.

After the polynomial unit-pivot descent and a base gauge, the normalized
potential and collision points are

    psi_0 = 2*y*r + 4*x*s,
    +/- p = +/-(1, -3/2, 6, 81/8).

Write a homogeneous quartic as a sum of base--dual bidegrees, where the
base variables are x,y and the dual variables are r,s.  This checker treats
three dense classes:

* every genuinely mixed term, of bidegree (1,3), (2,2), or (3,1);
* those terms together with an arbitrary pure-base quartic;
* those terms together with an arbitrary pure-dual quartic.

These spaces have respectively 25, 30, and 30 coefficients.  They are not
bounded-support searches.

For a homogeneous quartic h, collision of the gradients at +/-p is

    grad(h)(p) = -H_0*p,

where H_0=Hess(psi_0).  Constant Hessian determinant must equal 64 because
Hess(h) vanishes at the origin.  The spatial-degree-two determinant layer is
linear in h:

    trace(adj(H_0)*Hess(h)) = 0.

Collision and this linear layer have combined rank 14 in each class, leaving
11, 16, and 16 parameters.  A sparse polynomial-ring determinant expansion
then extracts the remaining spatial coefficients without asking SymPy to
expand a large nested expression.  Singular over QQ proves that each exact
coefficient ideal is the unit ideal.

The union of the pure-base and pure-dual sectors is deliberately not tested
by this coefficient calculation.  The complete 35-term homogeneous quartic
system is closed separately by the de Bondt--van den Essen reduction checked
in verify_hc4_meng_full_quartic_reduction.py.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations
import shutil
import subprocess

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.rings import PolyElement, ring


SINGULAR = shutil.which("Singular")
if SINGULAR is None:
    raise SystemExit("ERROR: Singular is required for this exact checker")


x, y, r, s = sp.symbols("x y r s")
spatial_variables = (x, y, r, s)
zero_exponent = (0, 0, 0, 0)

base_potential = 2 * y * r + 4 * x * s
base_hessian = sp.hessian(base_potential, spatial_variables)
base_adjugate = base_hessian.adjugate()
assert base_hessian.det() == 64

collision_point = {
    x: sp.Rational(1),
    y: -sp.Rational(3, 2),
    r: sp.Rational(6),
    s: sp.Rational(81, 8),
}
collision_target = -base_hessian * sp.Matrix(spatial_variables)


@dataclass(frozen=True)
class Chart:
    name: str
    base_degrees: tuple[int, ...]
    coefficient_count: int
    free_count: int
    determinant_coefficient_count: int


charts = (
    Chart("mixed", (1, 2, 3), 25, 11, 262),
    Chart("mixed_plus_pure_base", (1, 2, 3, 4), 30, 16, 273),
    Chart("mixed_plus_pure_dual", (0, 1, 2, 3), 30, 16, 273),
)


def quartic_exponents(base_degrees: tuple[int, ...]):
    exponents = []
    for base_degree in base_degrees:
        dual_degree = 4 - base_degree
        for x_degree in range(base_degree + 1):
            for r_degree in range(dual_degree + 1):
                exponents.append(
                    (
                        x_degree,
                        base_degree - x_degree,
                        r_degree,
                        dual_degree - r_degree,
                    )
                )
    return tuple(exponents)


def monomial(exponents):
    return sp.prod(
        variable**exponent
        for variable, exponent in zip(
            spatial_variables,
            exponents,
            strict=True,
        )
    )


def permutation_sign(permutation):
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(4)
        for right in range(left + 1, 4)
    )
    return -1 if inversions % 2 else 1


def add_term(polynomial, exponent, coefficient, coefficient_ring):
    if not coefficient:
        return
    polynomial[exponent] = (
        polynomial.get(exponent, coefficient_ring.zero) + coefficient
    )
    if not polynomial[exponent]:
        del polynomial[exponent]


def multiply_spatial_polynomials(left, right, coefficient_ring):
    product = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = tuple(
                left_value + right_value
                for left_value, right_value in zip(
                    left_exponent,
                    right_exponent,
                    strict=True,
                )
            )
            add_term(
                product,
                exponent,
                left_coefficient * right_coefficient,
                coefficient_ring,
            )
    return product


def sparse_hessian_determinant(
    exponents,
    coefficient_maps,
    coefficient_ring,
):
    matrix = [[{} for _ in range(4)] for _ in range(4)]

    for row in range(4):
        for column in range(4):
            if base_hessian[row, column]:
                add_term(
                    matrix[row][column],
                    zero_exponent,
                    coefficient_ring.domain.convert(
                        base_hessian[row, column]
                    ),
                    coefficient_ring,
                )

    for coefficient, exponent in zip(
        coefficient_maps,
        exponents,
        strict=True,
    ):
        for row in range(4):
            for column in range(4):
                factor = exponent[row] * (
                    exponent[column] - (1 if row == column else 0)
                )
                if not factor:
                    continue
                derivative_exponent = list(exponent)
                derivative_exponent[row] -= 1
                derivative_exponent[column] -= 1
                add_term(
                    matrix[row][column],
                    tuple(derivative_exponent),
                    factor * coefficient,
                    coefficient_ring,
                )

    determinant = {}
    for permutation in permutations(range(4)):
        product = {zero_exponent: coefficient_ring.one}
        for row in range(4):
            product = multiply_spatial_polynomials(
                product,
                matrix[row][permutation[row]],
                coefficient_ring,
            )
        sign = permutation_sign(permutation)
        for exponent, coefficient in product.items():
            add_term(
                determinant,
                exponent,
                sign * coefficient,
                coefficient_ring,
            )

    add_term(
        determinant,
        zero_exponent,
        coefficient_ring.domain.convert(-64),
        coefficient_ring,
    )
    return determinant


def singular_unit_ideal(parameter_names, equations):
    def singular_polynomial(polynomial: PolyElement):
        return str(polynomial).replace("**", "^")

    source = "\n".join(
        [
            f"ring q=0,({','.join(parameter_names)}),dp;",
            "option(redSB);",
            "ideal I="
            + ",".join(singular_polynomial(equation) for equation in equations)
            + ";",
            "ideal G=std(I);",
            'if (reduce(1,G)==0) { print("UNIT"); }',
            'else { print("NONUNIT"); print(size(G)); }',
        ]
    )
    result = subprocess.run(
        [SINGULAR, "-q"],
        input=source,
        text=True,
        capture_output=True,
        check=False,
        timeout=300,
    )
    if result.returncode != 0:
        raise AssertionError(
            "Singular failed:\n" + result.stdout + "\n" + result.stderr
        )
    assert result.stdout.strip() == "UNIT", result.stdout


def verify_chart(chart: Chart):
    exponents = quartic_exponents(chart.base_degrees)
    assert len(exponents) == chart.coefficient_count

    coefficients = sp.symbols(f"c0:{len(exponents)}")
    correction = sum(
        coefficient * monomial(exponent)
        for coefficient, exponent in zip(
            coefficients,
            exponents,
            strict=True,
        )
    )

    collision_equations = [
        sp.expand(
            sp.diff(correction, variable).subs(collision_point)
            - collision_target[index].subs(collision_point)
        )
        for index, variable in enumerate(spatial_variables)
    ]
    linear_determinant_layer = sp.expand(
        sum(
            base_adjugate[row, column]
            * sp.diff(
                correction,
                spatial_variables[row],
                spatial_variables[column],
            )
            for row in range(4)
            for column in range(4)
        )
    )
    linear_equations = collision_equations + sp.Poly(
        linear_determinant_layer,
        *spatial_variables,
    ).coeffs()
    linear_matrix, linear_rhs = sp.linear_eq_to_matrix(
        linear_equations,
        coefficients,
    )

    assert linear_matrix.rank() == 14
    assert linear_matrix.row_join(linear_rhs).rank() == 14
    solution = next(iter(sp.linsolve((linear_matrix, linear_rhs), coefficients)))
    free_parameters = tuple(
        sorted(
            set().union(*(entry.free_symbols for entry in solution))
            & set(coefficients),
            key=lambda symbol: int(str(symbol)[1:]),
        )
    )
    assert len(free_parameters) == chart.free_count
    solution_substitution = dict(zip(coefficients, solution, strict=True))
    assert all(
        sp.expand(equation.subs(solution_substitution)) == 0
        for equation in linear_equations
    )

    coefficient_ring, *_ = ring(
        ",".join(str(parameter) for parameter in free_parameters),
        QQ,
    )
    coefficient_maps = tuple(
        coefficient_ring.from_expr(entry) for entry in solution
    )
    determinant = sparse_hessian_determinant(
        exponents,
        coefficient_maps,
        coefficient_ring,
    )

    assert len(determinant) == chart.determinant_coefficient_count
    assert {sum(exponent) for exponent in determinant} == {4, 6, 8}

    # Independent exact evaluation of the custom sparse determinant builder.
    parameter_values = tuple(range(1, len(free_parameters) + 1))
    point_values = (2, -1, 3, 1)
    sparse_value = sum(
        coefficient(*parameter_values)
        * sp.prod(
            value**degree
            for value, degree in zip(
                point_values,
                exponent,
                strict=True,
            )
        )
        for exponent, coefficient in determinant.items()
    )
    direct_substitution = {
        parameter: value
        for parameter, value in zip(
            free_parameters,
            parameter_values,
            strict=True,
        )
    }
    numeric_correction = correction.subs(solution_substitution).subs(
        direct_substitution
    )
    direct_value = (
        sp.hessian(
            base_potential + numeric_correction,
            spatial_variables,
        )
        .subs(dict(zip(spatial_variables, point_values, strict=True)))
        .det()
        - 64
    )
    assert sp.Rational(sparse_value) == direct_value

    singular_unit_ideal(
        tuple(str(parameter) for parameter in free_parameters),
        tuple(determinant.values()),
    )
    print(
        f"PASS {chart.name}: {chart.coefficient_count} coefficients, "
        f"linear rank 14, {chart.free_count} nonlinear parameters, "
        f"{chart.determinant_coefficient_count} exact determinant equations"
    )


for chart in charts:
    verify_chart(chart)

print(
    "PASS: dense mixed quartics, even with either pure quartic sector, "
    "cannot retain the Meng collision and constant Hessian determinant"
)
print(
    "SCOPE: the simultaneous pure sectors are not tested here but are "
    "closed by the separate full-quartic reduction; mixed homogeneous "
    "degrees and non-coordinate coisotropic embeddings remain open"
)
