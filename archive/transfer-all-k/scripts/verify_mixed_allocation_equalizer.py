"""Quadratic-cubic equalizer audit for transfer vector (2,-1,-1)."""

import sympy as sp

W = sp.symbols("W")


def monic_deformation(prefix, degree, base):
    coefficients = sp.symbols(f"{prefix}0:{degree}")
    polynomial = sp.expand(
        base+sum(coefficients[index]*W**index for index in range(degree))
    )
    return polynomial, coefficients


paired_roots = (W-1)*(W-2)
Q_left, q_left = monic_deformation("ql", 6, W**6)
R_left, r_left = monic_deformation("rl", 4, paired_roots**2)
Q_right, q_right = monic_deformation("qr", 6, paired_roots**3)
R_right, r_right = monic_deformation("rr", 4, W**4)
variables = q_left+r_left+q_right+r_right
origin = {variable: 0 for variable in variables}

difference = sp.Poly(
    sp.expand(Q_left**2*R_left**3-Q_right**2*R_right**3), W
)
equations = sp.Matrix([
    difference.coeff_monomial(W**degree) for degree in range(2, 24)
])
jacobian = equations.jacobian(variables).subs(origin)
assert jacobian.rank() == 12

# Twelve normal variables/equations for formal implicit elimination.
pivot_columns = jacobian.rref()[1]
normal_matrix = sp.zeros(20, 12)
for column, variable_index in enumerate(pivot_columns):
    normal_matrix[variable_index, column] = 1
pivot_rows = (jacobian[:, list(pivot_columns)].T.rref()[1])
linear_block = jacobian[list(pivot_rows), list(pivot_columns)]
assert linear_block.det() != 0

# Four reduced directions: two coefficients of the k=2 factor at zero and
# one root-position parameter for each elementary compensating transfer.
p, q, r1, r2 = sp.symbols("p q r1 r2")
S0 = W**2+p*W+q
T12 = (W-1+r1)*(W-2+r2)
reduced_polynomials = (S0**3, T12**2, T12**3, S0**2)
reduced_parameters = (p, q, r1, r2)
reduced_origin = {p: 0, q: 0, r1: 0, r2: 0}
reduced_tangent = sp.zeros(20, 4)
for column, parameter in enumerate(reduced_parameters):
    entries = []
    for polynomial, degree in zip(reduced_polynomials, (6, 4, 6, 4)):
        derivative = sp.diff(polynomial, parameter).subs(reduced_origin)
        entries.extend(
            sp.Poly(derivative, W).coeff_monomial(W**power)
            for power in range(degree)
        )
    reduced_tangent[:, column] = sp.Matrix(entries)
assert reduced_tangent.rank() == 4
assert (jacobian*reduced_tangent).rank() == 0

tangent_basis = reduced_tangent
for vector in jacobian.nullspace():
    if tangent_basis.row_join(vector).rank() > tangent_basis.rank():
        tangent_basis = tangent_basis.row_join(vector)
assert tangent_basis.shape == (20, 8)
assert tangent_basis.row_join(normal_matrix).rank() == 20

e0, e1, e2, e3 = sp.symbols("e0 e1 e2 e3")
transverse_variables = (e0, e1, e2, e3)
linear_transverse = tangent_basis[:, 4:]*sp.Matrix(transverse_variables)
linear_substitution = {
    variables[index]: linear_transverse[index] for index in range(20)
}


def homogeneous_part(polynomial, degree):
    poly = sp.Poly(polynomial, *variables)
    return sp.Add(*[
        coefficient*sp.prod(
            variable**exponent
            for variable, exponent in zip(variables, monomial)
        )
        for monomial, coefficient in poly.terms()
        if sum(monomial) == degree
    ])


entries = list(equations)
quadratic_parts = [homogeneous_part(entry, 2) for entry in entries]
cubic_parts = [homogeneous_part(entry, 3) for entry in entries]
quadratic_selected = sp.Matrix([
    sp.expand(quadratic_parts[index].subs(linear_substitution))
    for index in pivot_rows
])
second_correction = -linear_block.inv()*quadratic_selected
normal_second = normal_matrix*second_correction

quadratic_cross = []
for quadratic in quadratic_parts:
    gradient = [
        sp.diff(quadratic, variable).subs(linear_substitution)
        for variable in variables
    ]
    quadratic_cross.append(sp.expand(sum(
        gradient[index]*normal_second[index] for index in range(20)
    )))
cubic_on_tangent = [
    sp.expand(cubic.subs(linear_substitution)) for cubic in cubic_parts
]
cubic_selected = sp.Matrix([
    quadratic_cross[index]+cubic_on_tangent[index] for index in pivot_rows
])
third_correction = -linear_block.inv()*cubic_selected

quadratic_obstruction = (
    sp.Matrix([
        sp.expand(quadratic.subs(linear_substitution))
        for quadratic in quadratic_parts
    ])
    + jacobian*normal_matrix*second_correction
)
cubic_obstruction = (
    jacobian*normal_matrix*third_correction
    + sp.Matrix(quadratic_cross)
    + sp.Matrix(cubic_on_tangent)
)
obstructions = [
    sp.factor(quadratic+cubic)
    for quadratic, cubic in zip(quadratic_obstruction, cubic_obstruction)
    if sp.factor(quadratic+cubic) != 0
]

expected_initial = [
    e0**2-20*e1**2,
    2*e0*e1+9*e1**2,
    e1**3,
    e2**2,
    e2*e3,
    e3**3,
]
computed = sp.groebner(
    obstructions, *transverse_variables, order="grevlex", domain=sp.QQ
)
expected = sp.groebner(
    expected_initial, *transverse_variables, order="grevlex", domain=sp.QQ
)
assert all(expected.reduce(obstruction)[1] == 0 for obstruction in obstructions)
assert all(computed.reduce(relation)[1] == 0 for relation in expected_initial)
assert computed.is_zero_dimensional

first_basis = (sp.Integer(1), e0, e1, e1**2)
second_basis = (sp.Integer(1), e2, e3, e3**2)
standard_basis = tuple(left*right for left in first_basis for right in second_basis)
assert len(standard_basis) == 16
assert all(expected.reduce(element)[1] == element for element in standard_basis)

print("PASS: the mixed (2,-1,-1) affine tangent space has dimension 8")
print("PASS: its quadratic-cubic initial ideal has colength 16")
print("PASS: the affine equalizer equals Z_2 tensor Z_1 tensor Z_1")
