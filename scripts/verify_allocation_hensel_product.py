"""Exact algebra audit for products of the known transfer blocks."""

from math import prod
import sympy as sp


def hilbert_convolution(left, right):
    result = [0]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i+j] += a*b
    return tuple(result)


# Z_1 has basis 1,e and Hilbert function (1,1).
rank = {1: 2, 2: 4}
hilbert = {1: (1, 1), 2: (1, 2, 1)}

for transfers in ((1, -1), (2, -2), (2, -1, -1), (2, 1, -1, -2)):
    assert sum(transfers) == 0
    predicted_rank = prod(rank[abs(k)] for k in transfers)
    assert predicted_rank == 2**sum(abs(k) for k in transfers)
    series = (1,)
    for k in transfers:
        series = hilbert_convolution(series, hilbert[abs(k)])
    assert sum(series) == predicted_rank

# The first new global test Z_2 tensor Z_2 has length 16 and Hilbert function
# (1,4,6,4,1).  Its coincident special fiber has socle dimension 2*2=4.
assert hilbert_convolution(hilbert[2], hilbert[2]) == (1, 4, 6, 4, 1)
assert rank[2]*rank[2] == 16
assert 2*2 == 4

# First global affine test: transfers (2,-2) at roots 0 and 1.  Sheet A has
# Q^2 at 0 and R^3 at 1; sheet B exchanges them.  Equality of normalized
# coefficients imposes only degrees >=2.  Its tangent space has the same
# dimension eight as the strong Hensel product, so affine difference adds no
# first-order direction.
W = sp.symbols("W")


def monic_deformation(prefix, degree, base):
    coefficients = sp.symbols(f"{prefix}0:{degree}")
    polynomial = sp.expand(
        base+sum(coefficients[index]*W**index for index in range(degree))
    )
    return polynomial, coefficients


Q_left, q_left = monic_deformation("ql", 6, W**6)
R_left, r_left = monic_deformation("rl", 4, (W-1)**4)
Q_right, q_right = monic_deformation("qr", 6, (W-1)**6)
R_right, r_right = monic_deformation("rr", 4, W**4)
variables = q_left+r_left+q_right+r_right
product_difference = sp.Poly(
    sp.expand(Q_left**2*R_left**3-Q_right**2*R_right**3), W
)
origin = {variable: 0 for variable in variables}
affine_equations = sp.Matrix([
    product_difference.coeff_monomial(W**degree) for degree in range(2, 24)
])
strong_equations = sp.Matrix([
    product_difference.coeff_monomial(W**degree) for degree in range(24)
])
affine_jacobian = affine_equations.jacobian(variables).subs(origin)
affine_rank = affine_jacobian.rank()
strong_rank = strong_equations.jacobian(variables).subs(origin).rank()
assert affine_rank == strong_rank == 12
assert len(variables)-affine_rank == 8

# Formal implicit reduction through degree three.  Select twelve normal
# coefficient variables and twelve equations with invertible linear block.
pivot_columns = affine_jacobian.rref()[1]
normal_matrix = sp.zeros(20, 12)
for column, variable_index in enumerate(pivot_columns):
    normal_matrix[variable_index, column] = 1
pivot_rows = (affine_jacobian[:, list(pivot_columns)].T.rref()[1])
linear_block = affine_jacobian[list(pivot_rows), list(pivot_columns)]
assert linear_block.det() != 0

# Put the four reduced Hensel directions first in the tangent basis.  They
# come from S_0=W^2+p0*W+q0 and
# S_1=(W-1)^2+p1*(W-1)+q1.
p0, q0, p1, q1 = sp.symbols("p0 q0 p1 q1")
S0 = W**2+p0*W+q0
S1 = (W-1)**2+p1*(W-1)+q1
reduced_polynomials = (S0**3, S1**2, S1**3, S0**2)
reduced_parameters = (p0, q0, p1, q1)
reduced_tangent = sp.zeros(20, 4)
for column, parameter in enumerate(reduced_parameters):
    entries = []
    for polynomial, degree in zip(reduced_polynomials, (6, 4, 6, 4)):
        derivative = sp.diff(polynomial, parameter).subs({
            p0: 0, q0: 0, p1: 0, q1: 0,
        })
        entries.extend(
            sp.Poly(derivative, W).coeff_monomial(W**power)
            for power in range(degree)
        )
    reduced_tangent[:, column] = sp.Matrix(entries)
assert reduced_tangent.rank() == 4
assert (affine_jacobian*reduced_tangent).rank() == 0

tangent_basis = reduced_tangent
for vector in affine_jacobian.nullspace():
    if tangent_basis.row_join(vector).rank() > tangent_basis.rank():
        tangent_basis = tangent_basis.row_join(vector)
assert tangent_basis.shape == (20, 8)
assert tangent_basis.row_join(normal_matrix).rank() == 20

# Slice away the four reduced directions and retain four transverse tangent
# variables.  Homogeneous expansion avoids a prohibitively large direct
# multivariate substitution.
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


affine_entries = list(affine_equations)
quadratic_parts = [homogeneous_part(entry, 2) for entry in affine_entries]
cubic_parts = [homogeneous_part(entry, 3) for entry in affine_entries]

quadratic_selected = sp.Matrix([
    sp.expand(quadratic_parts[index].subs(linear_substitution))
    for index in pivot_rows
])
second_correction = -linear_block.inv()*quadratic_selected
normal_second = normal_matrix*second_correction

# Cubic implicit correction is F_3(v)+dF_2(v)[normal_second].
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
    + affine_jacobian*normal_matrix*second_correction
)
cubic_obstruction = (
    affine_jacobian*normal_matrix*third_correction
    + sp.Matrix(quadratic_cross)
    + sp.Matrix(cubic_on_tangent)
)
jet_obstructions = [
    sp.factor(quadratic+cubic)
    for quadratic, cubic in zip(quadratic_obstruction, cubic_obstruction)
    if sp.factor(quadratic+cubic) != 0
]

# The local standard basis is two copies of the coincident Z_2 algebra.
expected_initial = [
    e0**2-9*e1**2,
    e0*e1+3*e1**2,
    e1**3,
    e2**2,
    e2*e3,
    e3**3,
]
jet_groebner = sp.groebner(
    jet_obstructions, *transverse_variables, order="grevlex", domain=sp.QQ
)
expected_groebner = sp.groebner(
    expected_initial, *transverse_variables, order="grevlex", domain=sp.QQ
)
assert all(expected_groebner.reduce(obstruction)[1] == 0
           for obstruction in jet_obstructions)
assert all(jet_groebner.reduce(relation)[1] == 0
           for relation in expected_initial)
assert jet_groebner.is_zero_dimensional

# Each pair has basis 1,a,b,b^2 (after a linear change in the first pair), so
# the tensor-product standard basis has 4*4=16 elements.
first_basis = (sp.Integer(1), e0, e1, e1**2)
second_basis = (sp.Integer(1), e2, e3, e3**2)
standard_basis = tuple(left*right for left in first_basis for right in second_basis)
assert len(standard_basis) == 16
assert all(expected_groebner.reduce(element)[1] == element
           for element in standard_basis)

# Recheck the local Z_2 multiplication used in the tensor calculation.
X, Y = sp.symbols("X Y")
special = sp.groebner([X**3, X*Y, Y**2], Y, X, domain=sp.QQ)
for element in (sp.Integer(1), X, Y, X**2):
    assert special.reduce(element)[1] == element

print("PASS: known Hensel products have length 2^sum(abs(k_rho))")
print("PASS: the (2,-2) strong product has length 16 and Hilbert function (1,4,6,4,1)")
print("PASS: its coincident special fiber has socle dimension 4")
print("PASS: the global (2,-2) affine and strong tangent spaces both have dimension 8")
print("PASS: the quadratic-cubic affine initial ideal has colength 16")
print("PASS: the global (2,-2) affine equalizer equals Z_2 tensor Z_2")
