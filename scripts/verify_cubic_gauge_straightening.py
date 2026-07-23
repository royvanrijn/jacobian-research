#!/usr/bin/env python3
"""Exact regressions for nonlinear cubic-slice gauge straightening."""

from __future__ import annotations

import itertools

import sympy as sp


# A polynomial upper-triangular GL_2 gauge has constant diagonal.  This
# degree-one saturated model checks the unit-determinant mechanism exactly;
# the written proof uses degree additivity in all degrees.
borel_t = sp.symbols("borel_t")
alpha_0, alpha_1, delta_0, delta_1, determinant_inverse = sp.symbols(
    "alpha_0 alpha_1 delta_0 delta_1 determinant_inverse"
)
borel_diagonal_product = sp.expand(
    (alpha_0 + alpha_1 * borel_t)
    * (delta_0 + delta_1 * borel_t)
)
assert sp.Poly(
    borel_diagonal_product,
    borel_t,
).all_coeffs() == [
    alpha_1 * delta_1,
    alpha_0 * delta_1 + alpha_1 * delta_0,
    alpha_0 * delta_0,
]
borel_unit_ideal = sp.groebner(
    (
        alpha_1 * delta_1,
        alpha_0 * delta_1 + alpha_1 * delta_0,
        alpha_0 * delta_0 * determinant_inverse - 1,
    ),
    determinant_inverse,
    delta_1,
    alpha_1,
    delta_0,
    alpha_0,
)
assert borel_unit_ideal.reduce(alpha_1**2)[1] == 0
assert borel_unit_ideal.reduce(delta_1**2)[1] == 0


# Universal normalized linear-times-quadratic factorization.
a, b, c, d, e, h, q = sp.symbols("a b c d e h q")
resultant = a**2 * e - a * b * d + b**2 * c
coefficients = (
    a * c,
    a * d + b * c,
    a * e + b * d,
    b * e,
)
A, C1, B, C = coefficients

# The target-dependent determinant-one variable shear T -> T+hS.
transformed_factors = (
    a,
    b + a * h,
    c,
    d + 2 * c * h,
    e + d * h + c * h**2,
)
aa, bb, cc, dd, ee = transformed_factors
transformed_resultant = aa**2 * ee - aa * bb * dd + bb**2 * cc
transformed_coefficients = (
    aa * cc,
    aa * dd + bb * cc,
    aa * ee + bb * dd,
    bb * ee,
)

assert sp.expand(transformed_resultant - resultant) == 0
assert tuple(map(sp.expand, transformed_coefficients)) == tuple(
    map(
        sp.expand,
        (
            A,
            C1 + 3 * A * h,
            B + 2 * C1 * h + 3 * A * h**2,
            C + B * h + C1 * h**2 + A * h**3,
        ),
    )
)

# The shear parameter may depend on any invariant of the upper-unipotent
# action, not just A.  Two standard translation invariants are the depressed
# cubic coefficients P and Q.
A_target, C1_target, B_target, C_target = sp.symbols(
    "A_target C1_target B_target C_target"
)
target_transformed = (
    A_target,
    C1_target + 3 * A_target * h,
    B_target + 2 * C1_target * h + 3 * A_target * h**2,
    C_target + B_target * h + C1_target * h**2 + A_target * h**3,
)
P_invariant = 3 * A_target * B_target - C1_target**2
Q_invariant = (
    27 * A_target**2 * C_target
    - 9 * A_target * C1_target * B_target
    + 2 * C1_target**3
)
discriminant_invariant = (
    C1_target**2 * B_target**2
    - 4 * A_target * B_target**3
    - 4 * C1_target**3 * C_target
    - 27 * A_target**2 * C_target**2
    + 18 * A_target * C1_target * B_target * C_target
)
assert sp.expand(
    P_invariant.subs(
        {
            C1_target: target_transformed[1],
            B_target: target_transformed[2],
        },
        simultaneous=True,
    )
    - P_invariant
) == 0
assert sp.expand(
    Q_invariant.subs(
        {
            C1_target: target_transformed[1],
            B_target: target_transformed[2],
            C_target: target_transformed[3],
        },
        simultaneous=True,
    )
    - Q_invariant
) == 0
assert sp.expand(
    discriminant_invariant.subs(
        {
            C1_target: target_transformed[1],
            B_target: target_transformed[2],
            C_target: target_transformed[3],
        },
        simultaneous=True,
    )
    - discriminant_invariant
) == 0
assert sp.expand(
    Q_invariant**2
    + 4 * P_invariant**3
    + 27 * A_target**2 * discriminant_invariant
) == 0

# If the shear time is an arbitrary polynomial h(A,C1,B,C), the ambient
# Jacobian is 1+D_+(h), where
# D_+=3*A*d/dC1+2*C1*d/dB+B*d/dC.
h_A, h_C1, h_B, h_C = sp.symbols("h_A h_C1 h_B h_C")
target_variables = (A_target, C1_target, B_target, C_target)
h_derivatives = (h_A, h_C1, h_B, h_C)
upper_variable_time_jacobian = sp.Matrix(
    4,
    4,
    lambda row, column: (
        sp.diff(target_transformed[row], target_variables[column])
        + sp.diff(target_transformed[row], h) * h_derivatives[column]
    ),
).det()
assert sp.factor(
    upper_variable_time_jacobian
    - (1 + 3 * A_target * h_C1 + 2 * C1_target * h_B + B_target * h_C)
) == 0

# On C1=q-3Ah the nonlinear slice becomes C1'=q.  The remaining
# coefficient change is a triangular polynomial automorphism of A^3_(A,B,C).
slice_relation = {C1_target: q - 3 * A_target * h}
B_shift = 2 * q * h - 3 * A_target * h**2
C_shift = h * B_target + q * h**2 - 2 * A_target * h**3
assert sp.expand(target_transformed[1].subs(slice_relation) - q) == 0
assert sp.expand(
    target_transformed[2].subs(slice_relation) - (B_target + B_shift)
) == 0
assert sp.expand(
    target_transformed[3].subs(slice_relation) - (C_target + C_shift)
) == 0

B_prime, C_prime = sp.symbols("B_prime C_prime")
B_inverse = B_prime - B_shift
C_inverse = sp.expand(
    C_prime
    - h * B_inverse
    - q * h**2
    + 2 * A_target * h**3
)
assert sp.expand(
    B_inverse.subs(B_prime, B_target + B_shift) - B_target
) == 0
assert sp.expand(
    C_inverse.subs(
        {
            B_prime: B_target + B_shift,
            C_prime: C_target + C_shift,
        },
        simultaneous=True,
    )
    - C_target
) == 0

# The opposite determinant-one shear S -> S+kT gives the symmetric gauge
# family based at C2=q.
k = sp.symbols("k")
lower_transformed_factors = (
    a + b * k,
    b,
    c + d * k + e * k**2,
    d + 2 * e * k,
    e,
)
la, lb, lc, ld, le = lower_transformed_factors
lower_resultant = la**2 * le - la * lb * ld + lb**2 * lc
lower_coefficients = (
    la * lc,
    la * ld + lb * lc,
    la * le + lb * ld,
    lb * le,
)
assert sp.expand(lower_resultant - resultant) == 0
assert tuple(map(sp.expand, lower_coefficients)) == tuple(
    map(
        sp.expand,
        (
            A + C1 * k + B * k**2 + C * k**3,
            C1 + 2 * B * k + 3 * C * k**2,
            B + 3 * C * k,
            C,
        ),
    )
)

lower_target_transformed = (
    A_target + C1_target * k + B_target * k**2 + C_target * k**3,
    C1_target + 2 * B_target * k + 3 * C_target * k**2,
    B_target + 3 * C_target * k,
    C_target,
)
lower_I = 3 * C_target * C1_target - B_target**2
lower_J = (
    27 * C_target**2 * A_target
    - 9 * C_target * B_target * C1_target
    + 2 * B_target**3
)
assert sp.expand(
    lower_I.subs(
        {
            C1_target: lower_target_transformed[1],
            B_target: lower_target_transformed[2],
        },
        simultaneous=True,
    )
    - lower_I
) == 0
assert sp.expand(
    lower_J.subs(
        {
            A_target: lower_target_transformed[0],
            C1_target: lower_target_transformed[1],
            B_target: lower_target_transformed[2],
        },
        simultaneous=True,
    )
    - lower_J
) == 0
assert sp.expand(
    lower_target_transformed[2].subs(
        B_target,
        q - 3 * C_target * k,
    )
    - q
) == 0

k_A, k_C1, k_B, k_C = sp.symbols("k_A k_C1 k_B k_C")
k_derivatives = (k_A, k_C1, k_B, k_C)
lower_variable_time_jacobian = sp.Matrix(
    4,
    4,
    lambda row, column: (
        sp.diff(lower_target_transformed[row], target_variables[column])
        + sp.diff(lower_target_transformed[row], k) * k_derivatives[column]
    ),
).det()
assert sp.factor(
    lower_variable_time_jacobian
    - (
        1
        + C1_target * k_A
        + 2 * B_target * k_C1
        + 3 * C_target * k_B
    )
) == 0

# For an alternating gauge U(p)L(k), with both parameters evaluated in the
# original coordinates, the Jacobian is the exact rank-two update from
# Proposition 2.4.
p = h
q_lower = k
lower_then_upper = (
    lower_target_transformed[0],
    lower_target_transformed[1]
    + 3 * lower_target_transformed[0] * p,
    lower_target_transformed[2]
    + 2 * lower_target_transformed[1] * p
    + 3 * lower_target_transformed[0] * p**2,
    lower_target_transformed[3]
    + lower_target_transformed[2] * p
    + lower_target_transformed[1] * p**2
    + lower_target_transformed[0] * p**3,
)
alternating_jacobian = sp.Matrix(
    4,
    4,
    lambda row, column: (
        sp.diff(lower_then_upper[row], target_variables[column])
        + sp.diff(lower_then_upper[row], p) * h_derivatives[column]
        + sp.diff(lower_then_upper[row], q_lower) * k_derivatives[column]
    ),
).det()


def apply_gradient(vector_field, gradient):
    return sp.expand(
        sum(
            coefficient * derivative
            for coefficient, derivative in zip(vector_field, gradient)
        )
    )


d_plus_vector = (0, 3 * A_target, 2 * C1_target, B_target)
d_minus_vector = (C1_target, 2 * B_target, 3 * C_target, 0)
h_weight_vector = (
    3 * A_target,
    C1_target,
    -B_target,
    -3 * C_target,
)
x_q_vector = tuple(
    sp.expand(
        d_plus_vector[index]
        - q_lower * h_weight_vector[index]
        - q_lower**2 * d_minus_vector[index]
    )
    for index in range(4)
)
x_q_p = apply_gradient(x_q_vector, h_derivatives)
x_q_q = apply_gradient(x_q_vector, k_derivatives)
d_minus_p = apply_gradient(d_minus_vector, h_derivatives)
d_minus_q = apply_gradient(d_minus_vector, k_derivatives)
alternating_jacobian_formula = sp.expand(
    (1 + x_q_p) * (1 + d_minus_q) - d_minus_p * x_q_q
)
alternating_jacobian_discrepancy = sp.factor(
    alternating_jacobian - alternating_jacobian_formula
)
assert alternating_jacobian_discrepancy == 0

# Expanding the rank-two formula separates the linear divergence, the
# bilinear sl2 divergence, and the genuinely new cubic interaction.
d_plus_p = apply_gradient(d_plus_vector, h_derivatives)
d_plus_q = apply_gradient(d_plus_vector, k_derivatives)
h_p = apply_gradient(h_weight_vector, h_derivatives)
h_q = apply_gradient(h_weight_vector, k_derivatives)
bilinear_interaction = sp.expand(
    -q_lower * h_p
    + d_plus_p * d_minus_q
    - d_minus_p * d_plus_q
)
cubic_interaction = sp.expand(
    q_lower
    * (
        d_minus_p * h_q
        - h_p * d_minus_q
        - q_lower * d_minus_p
    )
)
assert sp.expand(
    alternating_jacobian_formula
    - 1
    - d_plus_p
    - d_minus_q
    - bilinear_interaction
    - cubic_interaction
) == 0

# A nontrivial transported parameter checks the conjugated-kernel criterion:
# q=C3 is lower-invariant, r=C0 is upper-invariant in the intermediate
# coordinates, and p=r o psi_q is not presented in either original kernel.
transported_q = C_target
transported_p = sp.expand(
    lower_target_transformed[0].subs(k, transported_q)
)
transported_q_gradient = tuple(
    sp.diff(transported_q, variable) for variable in target_variables
)
transported_p_gradient = tuple(
    sp.diff(transported_p, variable) for variable in target_variables
)
assert apply_gradient(d_minus_vector, transported_q_gradient) == 0
transported_x_vector = tuple(
    sp.expand(component.subs(k, transported_q))
    for component in x_q_vector
)
transported_x_q = apply_gradient(
    transported_x_vector,
    transported_q_gradient,
)
transported_y_vector = tuple(
    sp.expand(
        transported_x_vector[index]
        - transported_x_q * d_minus_vector[index]
    )
    for index in range(4)
)
assert apply_gradient(
    transported_y_vector,
    transported_p_gradient,
) == 0

# The coefficient equations for arbitrary homogeneous-linear times have no
# genuinely alternating solution.  Their exact Groebner consequences leave
# p=a0*C0, q=b3*C3, with a0*b3=0.
linear_p_coefficients = sp.symbols("linear_p_0:4")
linear_q_coefficients = sp.symbols("linear_q_0:4")
linear_p = sum(
    linear_p_coefficients[index] * target_variables[index]
    for index in range(4)
)
linear_q = sum(
    linear_q_coefficients[index] * target_variables[index]
    for index in range(4)
)
linear_p_gradient = tuple(
    sp.diff(linear_p, variable) for variable in target_variables
)
linear_q_gradient = tuple(
    sp.diff(linear_q, variable) for variable in target_variables
)
linear_x_vector = tuple(
    sp.expand(
        d_plus_vector[index]
        - linear_q * h_weight_vector[index]
        - linear_q**2 * d_minus_vector[index]
    )
    for index in range(4)
)
linear_jacobian_minus_one = sp.expand(
    (
        1 + apply_gradient(linear_x_vector, linear_p_gradient)
    )
    * (
        1 + apply_gradient(d_minus_vector, linear_q_gradient)
    )
    - apply_gradient(d_minus_vector, linear_p_gradient)
    * apply_gradient(linear_x_vector, linear_q_gradient)
    - 1
)
linear_coefficient_equations = sp.Poly(
    linear_jacobian_minus_one,
    *target_variables,
).coeffs()
linear_time_groebner = sp.groebner(
    linear_coefficient_equations,
    *linear_p_coefficients,
    *linear_q_coefficients,
    order="grevlex",
)
lp0, lp1, lp2, lp3 = linear_p_coefficients
lq0, lq1, lq2, lq3 = linear_q_coefficients
for necessary_relation in (
    lq0**2,
    lq1**2,
    lp1,
    lq2,
    2 * lp2 + lq0,
    lp3 + 2 * lq1,
    lp0 * lq3,
):
    assert linear_time_groebner.reduce(necessary_relation)[1] == 0
linear_classification = {
    lp1: 0,
    lp2: 0,
    lp3: 0,
    lq0: 0,
    lq1: 0,
    lq2: 0,
}
classified_linear_jacobian = sp.factor(
    linear_jacobian_minus_one.subs(linear_classification)
)
assert sp.factor(classified_linear_jacobian.subs(lq3, 0)) == 0
assert sp.factor(classified_linear_jacobian.subs(lp0, 0)) == 0

# Exhaust the first 1,156 nonlinear monomial pairs.  Saturating the
# Jacobian coefficient ideal by the product of the two scalar times proves
# that no pair through degree three has both times nonzero.
monomial_times = []
for monomial_degree in range(1, 4):
    for exponents in itertools.product(
        range(monomial_degree + 1),
        repeat=4,
    ):
        if sum(exponents) == monomial_degree:
            monomial_times.append(
                sp.prod(
                    target_variables[index] ** exponents[index]
                    for index in range(4)
                )
            )

monomial_p_scalar, monomial_q_scalar, saturation_inverse = sp.symbols(
    "monomial_p_scalar monomial_q_scalar saturation_inverse"
)
for monomial_p, monomial_q in itertools.product(
    monomial_times,
    repeat=2,
):
    monomial_p_gradient = tuple(
        sp.diff(monomial_p_scalar * monomial_p, variable)
        for variable in target_variables
    )
    scaled_monomial_q = monomial_q_scalar * monomial_q
    monomial_q_gradient = tuple(
        sp.diff(scaled_monomial_q, variable)
        for variable in target_variables
    )
    monomial_x_vector = tuple(
        sp.expand(
            d_plus_vector[index]
            - scaled_monomial_q * h_weight_vector[index]
            - scaled_monomial_q**2 * d_minus_vector[index]
        )
        for index in range(4)
    )
    monomial_jacobian_minus_one = sp.expand(
        (
            1
            + apply_gradient(
                monomial_x_vector,
                monomial_p_gradient,
            )
        )
        * (
            1
            + apply_gradient(
                d_minus_vector,
                monomial_q_gradient,
            )
        )
        - apply_gradient(
            d_minus_vector,
            monomial_p_gradient,
        )
        * apply_gradient(
            monomial_x_vector,
            monomial_q_gradient,
        )
        - 1
    )
    monomial_coefficient_equations = sp.Poly(
        monomial_jacobian_minus_one,
        *target_variables,
    ).coeffs()
    nonzero_time_saturation = sp.groebner(
        (
            *monomial_coefficient_equations,
            saturation_inverse
            * monomial_p_scalar
            * monomial_q_scalar
            - 1,
        ),
        saturation_inverse,
        monomial_p_scalar,
        monomial_q_scalar,
        order="lex",
    )
    assert any(
        polynomial.as_expr() == 1
        for polynomial in nonzero_time_saturation.polys
    )

# The linearized homogeneous gauge equation has cokernel equal to the
# binary-cubic invariant line: zero off degrees divisible by four.
graded_cokernel_dimensions = []
for homogeneous_degree in range(1, 9):
    homogeneous_exponents = [
        exponents
        for exponents in itertools.product(
            range(homogeneous_degree + 1),
            repeat=4,
        )
        if sum(exponents) == homogeneous_degree
    ]
    homogeneous_basis = [
        sp.prod(
            target_variables[index] ** exponents[index]
            for index in range(4)
        )
        for exponents in homogeneous_exponents
    ]
    derivation_columns = []
    for basis_monomial in homogeneous_basis:
        for derivation_vector in (d_plus_vector, d_minus_vector):
            derivative = sum(
                derivation_vector[index]
                * sp.diff(
                    basis_monomial,
                    target_variables[index],
                )
                for index in range(4)
            )
            derivative_polynomial = sp.Poly(
                sp.expand(derivative),
                *target_variables,
            )
            derivation_columns.append(
                [
                    derivative_polynomial.coeff_monomial(
                        target_monomial
                    )
                    for target_monomial in homogeneous_basis
                ]
            )
    graded_operator = sp.Matrix(
        len(homogeneous_basis),
        2 * len(homogeneous_basis),
        lambda row, column: derivation_columns[column][row],
    )
    graded_cokernel_dimensions.append(
        len(homogeneous_basis) - graded_operator.rank()
    )
assert graded_cokernel_dimensions == [0, 0, 0, 1, 0, 0, 0, 1]


# The first genuinely nonlinear-looking test is
# C1+t*A^2=1, corresponding to h=t*A/3 and q=1.  Its source is A^3:
# start with the tangent slice, then add -t*c^2*S*L to Q.
x, y, z, t = sp.symbols("x y z t")
factor_b = 1 + x * y
factor_c = 1 - sp.Rational(3, 2) * x * y + x**2 * z
factor_d_0 = (
    sp.Rational(1, 2) * y
    - x * z
    + sp.Rational(3, 2) * x * y**2
    - x**2 * y * z
)
factor_e_0 = (
    -2 * z
    + 4 * y**2
    - 4 * x * y * z
    + 3 * x * y**3
    - 2 * x**2 * y**2 * z
)
factor_d_t = sp.expand(factor_d_0 - t * x * factor_c**2)
factor_e_t = sp.expand(factor_e_0 - t * factor_b * factor_c**2)

family_substitution = {
    a: x,
    b: factor_b,
    c: factor_c,
    d: factor_d_t,
    e: factor_e_t,
}
assert sp.expand((resultant - 1).subs(family_substitution)) == 0
assert sp.expand((C1 + t * A**2 - 1).subs(family_substitution)) == 0

family_map = tuple(
    sp.expand(expression.subs(family_substitution))
    for expression in (A, B, C)
)
assert sp.factor(
    sp.Matrix(family_map).jacobian((x, y, z)).det() + 1
) == 0

# Both quotient coordinates remain polynomial on the nonlinear source.
y_inverse = 2 * b * d - a * e + t * a * b * c**2
standard_z_inverse = (
    2 * d**2
    + c * e
    + 6 * b * d**2
    + 3 * b * c * e
    - sp.Rational(9, 2) * e
)
z_inverse = standard_z_inverse.subs(
    {
        d: d + t * a * c**2,
        e: e + t * b * c**2,
    },
    simultaneous=True,
)
assert sp.expand(y_inverse.subs(family_substitution) - y) == 0
assert sp.expand(z_inverse.subs(family_substitution) - z) == 0

print("PASS: determinant-one cubic-basis shears preserve the resultant")
print("PASS: I, J, and the discriminant are shear-invariant")
print("PASS: a variable-time shear has Jacobian 1+D(h)")
print("PASS: every polynomial Borel gauge reduces to a constant diagonal and one shear")
print("PASS: the alternating two-shear Jacobian is the exact rank-two formula")
print("PASS: the transported two-shear kernel criterion is exact")
print("PASS: normalized linear times admit no alternating cancellation")
print("PASS: 1,156 monomial pairs through degree three admit no cancellation")
print("PASS: graded gauge cokernels occur only in discriminant degrees 4 and 8")
print("PASS: C1=q-3*C0*h with h in ker(D_+) straightens to C1=q")
print("PASS: the induced target change is a polynomial automorphism")
print("PASS: the opposite invariant shear straightens the symmetric C2 slice")
print("PASS: the nonlinear C1+t*C0^2 family has source A^3 and Jacobian -1")
print("PASS: the nonlinear family is polynomially left-right gauge-equivalent")
