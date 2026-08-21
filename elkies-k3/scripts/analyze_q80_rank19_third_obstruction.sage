#!/usr/bin/env sage
"""Derive the third-order obstruction at the q=80 CM24 seed over GF(11).

The quadratic obstruction radical is an affine five-space in the nine kernel
coordinates.  Working symbolically on that five-space avoids an expensive
54-variable integer search and identifies exactly which first corrections
can lift through 11^4.
"""

from itertools import product

from sage.all import GF, Matrix, PolynomialRing, ZZ, vector


load("elkies-k3/scripts/verify_q80_cm24_mod7_seed.sage")

prime = ZZ(11)
field = GF(prime)
integer_parameters = PolynomialRing(ZZ, names=tuple(f"u{i}" for i in range(5)))
u_integer = integer_parameters.gens()
finite_parameters = PolynomialRing(
    field, names=tuple(f"u{i}" for i in range(5))
)
u = finite_parameters.gens()


def integer_lift(polynomial):
    polynomial = finite_parameters(polynomial)
    return integer_parameters(
        {
            monomial: ZZ(coefficient)
            for monomial, coefficient in polynomial.dict().items()
        }
    )


def divide_coefficients(polynomial, divisor):
    polynomial = integer_parameters(polynomial)
    assert all(coefficient % divisor == 0 for coefficient in polynomial.coefficients())
    return integer_parameters(
        {
            monomial: coefficient // divisor
            for monomial, coefficient in polynomial.dict().items()
        }
    )


# Radical parametrization, with (u0,...,u4)=(c5,...,c9).
c5, c6, c7, c8, c9 = u
surface_s = c5 + 5*c6 + 3*c7 - 3*c8 - c9
kernel_coordinates = vector(
    finite_parameters,
    (
        (surface_s-3)/5,
        -3*c5 - 4*c6 + 2*c7 - c8 + 5*c9,
        (-surface_s-1)/4,
        5,
        c5,
        c6,
        c7,
        c8,
        c9,
    ),
)
obstruction_substitution = dict(
    zip(tangent_field_variables, kernel_coordinates)
)
assert all(
    polynomial.subs(obstruction_substitution) == 0
    for polynomial in obstruction_polynomials
)

first_correction = []
for column in range(len(names)):
    value = finite_parameters(first_delta[column])
    value += sum(
        finite_parameters(kernel_coordinates[row]) * finite_parameters(kernel[row, column])
        for row in range(kernel.nrows())
    )
    first_correction.append(value)

seed_vector = vector(ZZ, seed_values)
point_mod_121 = vector(
    integer_parameters,
    [
        integer_parameters(seed_vector[column])
        + prime*integer_lift(first_correction[column])
        for column in range(len(names))
    ],
)
substitution_mod_121 = dict(zip(parameters.gens(), point_mod_121))
residual_mod_121 = [
    integer_parameters(equation.subs(substitution_mod_121))
    for equation in equations
]
rhs2 = vector(
    finite_parameters,
    [
        -finite_parameters(divide_coefficients(value, prime**2))
        for value in residual_mod_121
    ],
)
assert not any(left_kernel.change_ring(finite_parameters) * rhs2)

# A fixed polynomial right inverse of J on its image.
pivot_columns = tuple(jacobian.pivots())
column_basis = jacobian.matrix_from_columns(pivot_columns)
pivot_rows = tuple(column_basis.transpose().pivots())
square = jacobian.matrix_from_rows_and_columns(pivot_rows, pivot_columns)
assert square.is_invertible() and len(pivot_columns) == jacobian.rank()
selected_rhs = vector(finite_parameters, [rhs2[row] for row in pivot_rows])
selected_delta = square.change_ring(finite_parameters).solve_right(selected_rhs)
second_correction = vector(finite_parameters, [0] * len(names))
for column, value in zip(pivot_columns, selected_delta):
    second_correction[column] = value
assert jacobian.change_ring(finite_parameters) * second_correction == rhs2

point_mod_1331 = vector(
    integer_parameters,
    [
        point_mod_121[column]
        + prime**2*integer_lift(second_correction[column])
        for column in range(len(names))
    ],
)
substitution_mod_1331 = dict(zip(parameters.gens(), point_mod_1331))
residual_mod_1331 = [
    integer_parameters(equation.subs(substitution_mod_1331))
    for equation in equations
]
rhs3 = vector(
    finite_parameters,
    [
        -finite_parameters(divide_coefficients(value, prime**3))
        for value in residual_mod_1331
    ],
)
third_obstructions = tuple(
    polynomial
    for polynomial in left_kernel.change_ring(finite_parameters) * rhs3
    if polynomial
)
third_ideal = finite_parameters.ideal(third_obstructions)
print(
    f"Q80CM24THIRD|variables=5|equations={len(third_obstructions)}|"
    f"dimension={third_ideal.dimension()}",
    flush=True,
)
print(
    "Q80CM24THIRD|groebner="
    + ";".join(map(str, third_ideal.groebner_basis())),
    flush=True,
)

points = []
for values in product(field, repeat=5):
    if all(not polynomial(*values) for polynomial in third_obstructions):
        points.append(tuple(map(ZZ, values)))
print(
    f"Q80CM24THIRD|tested={prime**5}|points={len(points)}|first="
    + ";".join(",".join(map(str, point)) for point in points[:32]),
    flush=True,
)
assert points

# In this canonical right-inverse slice the third obstruction vanishes
# identically.  Continue once more; this is the first order that can separate
# the five quadratic-gate parameters.
selected_rhs3 = vector(finite_parameters, [rhs3[row] for row in pivot_rows])
selected_delta3 = square.change_ring(finite_parameters).solve_right(selected_rhs3)
third_correction = vector(finite_parameters, [0] * len(names))
for column, value in zip(pivot_columns, selected_delta3):
    third_correction[column] = value
assert jacobian.change_ring(finite_parameters) * third_correction == rhs3

point_mod_14641 = vector(
    integer_parameters,
    [
        point_mod_1331[column]
        + prime**3*integer_lift(third_correction[column])
        for column in range(len(names))
    ],
)
substitution_mod_14641 = dict(zip(parameters.gens(), point_mod_14641))
residual_mod_14641 = [
    integer_parameters(equation.subs(substitution_mod_14641))
    for equation in equations
]
rhs4 = vector(
    finite_parameters,
    [
        -finite_parameters(divide_coefficients(value, prime**4))
        for value in residual_mod_14641
    ],
)
fourth_obstructions = tuple(
    polynomial
    for polynomial in left_kernel.change_ring(finite_parameters) * rhs4
    if polynomial
)
fourth_ideal = finite_parameters.ideal(fourth_obstructions)
print(
    f"Q80CM24FOURTH|variables=5|equations={len(fourth_obstructions)}|"
    f"dimension={fourth_ideal.dimension()}|"
    f"degrees={','.join(map(str, sorted(set(map(lambda f: f.total_degree(), fourth_obstructions)))))}",
    flush=True,
)

fourth_points = []
for values in product(field, repeat=5):
    if all(not polynomial(*values) for polynomial in fourth_obstructions):
        fourth_points.append(tuple(map(ZZ, values)))
print(
    f"Q80CM24FOURTH|tested={prime**5}|points={len(fourth_points)}|first="
    + ";".join(",".join(map(str, point)) for point in fourth_points[:32]),
    flush=True,
)
