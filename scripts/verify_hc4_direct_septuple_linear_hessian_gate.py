#!/usr/bin/env python3
"""Exact checks for the septuple-linear Hessian pole--defect ladder.

The line-bundle, binary zero-Hessian, and coordinate-normalization arguments
in HC4_DIRECT_SEPTUPLE_LINEAR_HESSIAN_GATE.md are written proofs. This script
replays the exact degree ladder, tangent-pencil determinant identity, conic
boundary calibration, and first-normal obstruction.
"""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp


# K^2=O(B-6), together with a nonzero section of K(b-1), leaves six rows.
ladder: list[tuple[int, int, int]] = []
for pole_order in (1, 2, 3):
    for primitive_kernel_degree in range(pole_order):
        defect_length = 6 - 2 * primitive_kernel_degree
        assert defect_length >= 2 * (4 - pole_order)
        ladder.append(
            (pole_order, primitive_kernel_degree, defect_length)
        )

assert ladder == [
    (1, 0, 6),
    (2, 0, 6),
    (2, 1, 4),
    (3, 0, 6),
    (3, 1, 4),
    (3, 2, 2),
]


u, v = sp.symbols("u v")


# On the tangent-pencil arrow matrix, the first normal coefficient is
# -(4/3)*G*det(Hess(G)) for every binary quartic G.
quartic_coefficients = sp.symbols("g0:5")
G_general = sum(
    quartic_coefficients[index] * u ** (4 - index) * v**index
    for index in range(5)
)
G_u = sp.diff(G_general, u)
G_v = sp.diff(G_general, v)
G_hessian = sp.hessian(G_general, (u, v))
tangent_first_normal = sp.expand(
    -G_v**2 * G_hessian[0, 0]
    + 2 * G_u * G_v * G_hessian[0, 1]
    - G_u**2 * G_hessian[1, 1]
)
assert sp.expand(
    tangent_first_normal
    + sp.Rational(4, 3) * G_general * G_hessian.det()
) == 0


# A normalized 2+1+1 quartic has a nonzero first normal coefficient.
G_tangent = u**2 * v * (u - v)
assert sp.factor(sp.hessian(G_tangent, (u, v)).det()) != 0
tangent_substitution = {
    coefficient: sp.Poly(G_tangent, u, v).coeff_monomial(
        u ** (4 - index) * v**index
    )
    for index, coefficient in enumerate(quartic_coefficients)
}
assert sp.factor(tangent_first_normal.subs(tangent_substitution)) != 0

# Primitive degree-one kernels have projective-line image. In the tangent
# position their quartic jet has root type 3+1 or 2+2, and the same first
# normal identity makes the determinant line simple.
for linear_kernel_quartic in (u**3 * v, u**2 * v**2):
    gradient_gcd = sp.gcd(
        sp.diff(linear_kernel_quartic, u),
        sp.diff(linear_kernel_quartic, v),
    )
    assert sp.Poly(gradient_gcd, u, v).total_degree() == 2
    assert sp.factor(
        sp.hessian(linear_kernel_quartic, (u, v)).det()
    ) != 0


# A root of multiplicity m=1,...,4 in a binary quintic has exact Hessian
# multiplicity 2m-2. These are the four automatic square-divisor partitions.
for multiplicity in range(1, 5):
    monomial_quintic = u**multiplicity * v ** (5 - multiplicity)
    monomial_hessian = sp.factor(
        sp.hessian(monomial_quintic, (u, v)).det()
    )
    assert sp.Poly(monomial_hessian, u, v).degree(u) == 2 * multiplicity - 2

automatic_repeated_partitions = (
    (4, 1),
    (3, 2),
    (3, 1, 1),
    (2, 2, 1),
)
assert all(
    sum(2 * multiplicity - 2 for multiplicity in partition) >= 4
    for partition in automatic_repeated_partitions
)


# The one-double-root stratum has one exceptional anharmonic orbit.
t = sp.symbols("t")
F_double = sp.expand(u**2 * v * (u - v) * (u - t * v))
H_double = sp.factor(sp.hessian(F_double, (u, v)).det())
assert sp.rem(H_double, u**2, u) == 0
Q_double = sp.cancel(H_double / u**2)
Q_double_affine = sp.Poly(Q_double.subs(v, 1), u)
Q_discriminant = sp.factor(sp.discriminant(Q_double_affine.as_expr(), u))
exceptional_sextic = (
    4 * t**6
    - 12 * t**5
    - t**4
    + 22 * t**3
    - t**2
    - 12 * t
    + 4
)
assert sp.cancel(
    Q_discriminant
    / (t**2 * (t - 1) ** 2 * exceptional_sextic)
).is_number
assert sp.expand(
    exceptional_sextic
    - (
        4 * (t**2 - t + 1) ** 3
        - 25 * t**2 * (1 - t) ** 2
    )
) == 0
assert sp.expand(exceptional_sextic.subs(t, 1 - t) - exceptional_sextic) == 0
assert sp.expand(
    t**6 * exceptional_sextic.subs(t, 1 / t) - exceptional_sextic
) == 0


# The squarefree quintic stratum is a saturated degree-twelve scheme and one
# PGL2 orbit. Singular certifies the elimination from the factorization
# Hess(F)=q^2*lambda; SymPy checks the small lex basis and orbit arithmetic.
s = sp.symbols("s")
F_squarefree = sp.expand(
    u * v * (u - v) * (u - s * v) * (u - t * v)
)
H_squarefree_affine = sp.Poly(
    sp.expand(sp.hessian(F_squarefree, (u, v)).det().subs(v, 1)),
    u,
)
A, B, L1, L0 = sp.symbols("A B L1 L0")
q_square = u**2 + A * u + B
lambda_square = -16 * u**2 + L1 * u + L0
factor_difference = sp.Poly(
    sp.expand(
        H_squarefree_affine.as_expr() - q_square**2 * lambda_square
    ),
    u,
)
factor_equations = [
    factor_difference.coeff_monomial(u**index)
    for index in range(6)
]
L1_solution = sp.solve(factor_equations[5], L1)[0]
L0_solution = sp.solve(
    factor_equations[4].subs(L1, L1_solution),
    L0,
)[0]
reduced_factor_equations = []
for equation in factor_equations[:4]:
    specialized = sp.expand(
        equation.subs({L1: L1_solution, L0: L0_solution})
    )
    reduced_factor_equations.append(
        sp.Poly(specialized, A, B, s, t).primitive()[1].as_expr()
    )

squarefree_elimination = [
    s**3
    - s**2 * t
    - s * t**2
    + t**3
    - s**2
    + 2 * s * t
    - t**2
    - s
    - t
    + 1,
    s**2 * t**2
    + 2 * s * t**3
    - t**4
    - s**2 * t
    - 4 * s * t**2
    + t**3
    + s**2
    - 4 * s * t
    + 2 * t**2
    + 2 * s
    + t
    - 1,
    3 * s * t**4
    - t**5
    - 6 * s * t**3
    + t**4
    + 4 * s**2 * t
    - 7 * s * t**2
    + 3 * t**3
    - 2 * s**2
    + 6 * s * t
    + t**2
    - s
    - 3 * t
    + 1,
    t**6 - 3 * t**5 - 2 * t**4 + 9 * t**3 - 2 * t**2 - 3 * t + 1,
]


def singular_expression(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


singular = shutil.which("Singular")
if singular is None:
    raise RuntimeError("Singular is required for the squarefree elimination")

singular_program = f"""
ring r=0,(A,B,w,s,t),(dp(3),dp(2));
poly f0={singular_expression(reduced_factor_equations[0])};
poly f1={singular_expression(reduced_factor_equations[1])};
poly f2={singular_expression(reduced_factor_equations[2])};
poly f3={singular_expression(reduced_factor_equations[3])};
poly D=s*t*(s-1)*(t-1)*(s-t);
ideal I=f0,f1,f2,f3,w*D-1;
option(redSB);
ideal G=std(I);
ideal E=eliminate(G,A*B*w);
ideal X={",".join(
    singular_expression(polynomial)
    for polynomial in squarefree_elimination
)};
ideal GE=std(E);
ideal GX=std(X);
if (size(reduce(E,GX))!=0 || size(reduce(X,GE))!=0)
{{
  print("SQUAREFREE_ELIMINATION_FAILURE");
  exit(1);
}}
print("SQUAREFREE_ELIMINATION_CERTIFICATE");
"""
singular_result = subprocess.run(
    [singular, "-q"],
    input=singular_program,
    text=True,
    capture_output=True,
    check=True,
    timeout=120,
)
if singular_result.stderr.strip():
    raise RuntimeError(singular_result.stderr)
assert "SQUAREFREE_ELIMINATION_CERTIFICATE" in singular_result.stdout

squarefree_lex = sp.groebner(
    squarefree_elimination,
    s,
    t,
    order="lex",
)
assert len(squarefree_lex.polys) == 2
terminal_t = sp.factor(squarefree_lex.polys[1].as_expr())
t_factors = (
    t**2 - 3 * t + 1,
    t**2 - t - 1,
    t**2 + t - 1,
)
assert sp.expand(terminal_t - sp.prod(t_factors)) == 0

lex_s = squarefree_lex.polys[0].as_expr()
expected_s_remainders = (
    2 * (s**2 - t),
    2 * (s**2 - 2 * s * t + t),
    2 * (s**2 - 2 * s + t),
)
for t_factor, expected_remainder in zip(
    t_factors,
    expected_s_remainders,
):
    remainder = sp.rem(lex_s, t_factor, t)
    assert sp.expand(remainder - expected_remainder) == 0
    assert sp.resultant(t_factor, sp.discriminant(remainder, s), t) != 0

assert sp.gcd(terminal_t, sp.diff(terminal_t, t)) == 1
anharmonic_numerator = (t**2 - t + 1) ** 3
anharmonic_denominator = t**2 * (1 - t) ** 2
for t_factor in t_factors:
    assert sp.rem(
        anharmonic_numerator - 8 * anharmonic_denominator,
        t_factor,
        t,
    ) == 0

# The two s-roots for t=phi^2 are projectively equivalent, and an order-five
# Möbius map cycles the representative five-point set.
phi = sp.symbols("phi")
phi_relation = phi**2 - phi - 1


def reduce_phi(expression: sp.Expr) -> sp.Expr:
    return sp.rem(sp.expand(expression), phi_relation, phi)


assert reduce_phi(phi * (phi - 1) - 1) == 0
assert reduce_phi((2 - phi) * phi**2 - 1) == 0
assert reduce_phi(1 + (phi - 2) * phi**2) == 0
assert reduce_phi(1 + phi * (1 - phi)) == 0
assert reduce_phi(phi**2 - phi - 1) == 0

F_fermat_binary = u**5 + v**5
assert sp.factor(sp.hessian(F_fermat_binary, (u, v)).det()) == (
    400 * u**3 * v**3
)

# The unique squarefree orbit cannot carry the primitive conic kernel.
# The square-Hessian equation forces a=u*v. Polynomiality of the two tangent
# kernel coordinates forces G=alpha*u^4+beta*v^4, after which every kernel
# entry has the common factor u*v.
fermat_quartic_coefficients = sp.symbols("k0:5")
G_fermat = sum(
    fermat_quartic_coefficients[index]
    * u ** (4 - index)
    * v**index
    for index in range(5)
)
G_fermat_u = sp.diff(G_fermat, u)
G_fermat_v = sp.diff(G_fermat, v)
u_divisibility_remainder = sp.rem(G_fermat_u, u**2, u)
v_divisibility_remainder = sp.rem(G_fermat_v, v**2, v)
fermat_divisibility_equations = (
    sp.Poly(u_divisibility_remainder, u, v).coeffs()
    + sp.Poly(v_divisibility_remainder, u, v).coeffs()
)
fermat_divisibility_solution = sp.solve(
    fermat_divisibility_equations,
    fermat_quartic_coefficients[1:4],
    dict=True,
)
assert fermat_divisibility_solution == [
    {
        fermat_quartic_coefficients[1]: 0,
        fermat_quartic_coefficients[2]: 0,
        fermat_quartic_coefficients[3]: 0,
    }
]
G_fermat_reduced = sp.expand(
    G_fermat.subs(fermat_divisibility_solution[0])
)
assert G_fermat_reduced == (
    fermat_quartic_coefficients[0] * u**4
    + fermat_quartic_coefficients[4] * v**4
)
fermat_kernel = sp.Matrix(
    [
        u * v,
        -fermat_quartic_coefficients[0] * u * v / 5,
        -fermat_quartic_coefficients[4] * u * v / 5,
    ]
)
assert all(sp.rem(entry, u * v, u, v) == 0 for entry in fermat_kernel)


# Complete the repeated-root conic-kernel rows. For each displayed family,
# verify the tangent equations, its dimension, and the claimed common-factor
# or first-normal obstruction.
tangent_g_coefficients = sp.symbols("tg0:5")
tangent_b_coefficients = sp.symbols("tb0:3")
tangent_c_coefficients = sp.symbols("tc0:3")
tangent_generic_G = sum(
    tangent_g_coefficients[index] * u ** (4 - index) * v**index
    for index in range(5)
)
tangent_generic_b = sum(
    tangent_b_coefficients[index] * u ** (2 - index) * v**index
    for index in range(3)
)
tangent_generic_c = sum(
    tangent_c_coefficients[index] * u ** (2 - index) * v**index
    for index in range(3)
)
tangent_unknowns = (
    *tangent_g_coefficients,
    *tangent_b_coefficients,
    *tangent_c_coefficients,
)


def homogeneous_coefficients(
    polynomial: sp.Expr,
    degree: int,
) -> list[sp.Expr]:
    expanded = sp.Poly(sp.expand(polynomial), u, v)
    return [
        expanded.coeff_monomial(u ** (degree - index) * v**index)
        for index in range(degree + 1)
    ]


def verify_complete_tangent_family(
    binary_quintic: sp.Expr,
    normal_coordinate: sp.Expr,
    family_G: sp.Expr,
    family_b: sp.Expr,
    family_c: sp.Expr,
    parameters: tuple[sp.Symbol, ...],
) -> None:
    generic_equations = sp.hessian(binary_quintic, (u, v)) * sp.Matrix(
        [tangent_generic_b, tangent_generic_c]
    ) + normal_coordinate * sp.Matrix(
        [
            sp.diff(tangent_generic_G, u),
            sp.diff(tangent_generic_G, v),
        ]
    )
    coefficient_equations = []
    for equation in generic_equations:
        coefficient_equations.extend(homogeneous_coefficients(equation, 5))
    coefficient_matrix, _ = sp.linear_eq_to_matrix(
        coefficient_equations,
        tangent_unknowns,
    )
    assert len(coefficient_matrix.nullspace()) == len(parameters)

    family_equations = sp.hessian(binary_quintic, (u, v)) * sp.Matrix(
        [family_b, family_c]
    ) + normal_coordinate * sp.Matrix(
        [sp.diff(family_G, u), sp.diff(family_G, v)]
    )
    assert all(sp.expand(equation) == 0 for equation in family_equations)

    family_vector = sp.Matrix(
        homogeneous_coefficients(family_G, 4)
        + homogeneous_coefficients(family_b, 2)
        + homogeneous_coefficients(family_c, 2)
    )
    parameter_matrix = sp.Matrix.hstack(
        *(family_vector.diff(parameter) for parameter in parameters)
    )
    assert parameter_matrix.rank() == len(parameters)


p, q_parameter, r = sp.symbols("p q_parameter r")

# Root type 3+1+1.
F_311 = u**3 * v * (u - v)
a_311 = u**2
G_311 = -u**2 * (
    4 * p * u * v
    - 3 * p * v**2
    + q_parameter * u**2
    - 2 * q_parameter * u * v
)
b_311 = p * u**2
c_311 = q_parameter * u**2
verify_complete_tangent_family(
    F_311,
    a_311,
    G_311,
    b_311,
    c_311,
    (p, q_parameter),
)
assert all(
    sp.cancel(entry / a_311).as_numer_denom()[1] == 1
    for entry in (a_311, b_311, c_311)
)

# Root type 2+2+1.
F_221 = u**2 * v**2 * (u - v)
a_221 = u * v
G_221 = -u * v * (
    3 * p * u * v
    - 2 * p * v**2
    + 2 * q_parameter * u**2
    - 3 * q_parameter * u * v
)
b_221 = p * u * v
c_221 = q_parameter * u * v
verify_complete_tangent_family(
    F_221,
    a_221,
    G_221,
    b_221,
    c_221,
    (p, q_parameter),
)
assert all(
    sp.cancel(entry / a_221).as_numer_denom()[1] == 1
    for entry in (a_221, b_221, c_221)
)

# Root type 3+2, first square divisor a=u*v.
F_32 = u**2 * v**3
a_32_first = u * v
G_32_first = -u * v**2 * (
    2 * p * v + 3 * q_parameter * u
)
b_32_first = p * u * v
c_32_first = q_parameter * u * v
verify_complete_tangent_family(
    F_32,
    a_32_first,
    G_32_first,
    b_32_first,
    c_32_first,
    (p, q_parameter),
)
assert all(
    sp.cancel(entry / a_32_first).as_numer_denom()[1] == 1
    for entry in (a_32_first, b_32_first, c_32_first)
)

# Root type 3+2, second square divisor a=v^2. Polynomiality of H forces
# q_parameter=0, after which the kernel has common factor v^2.
a_32_second = v**2
G_32_second = -u * v * (
    4 * p * v**2
    + 3 * q_parameter * u**2
    + 6 * r * u * v
) / 2
b_32_second = (4 * p * v**2 - 3 * q_parameter * u**2) / 4
c_32_second = v * (q_parameter * u + r * v)
verify_complete_tangent_family(
    F_32,
    a_32_second,
    G_32_second,
    b_32_second,
    c_32_second,
    (p, q_parameter, r),
)
normal_numerator_32 = sp.expand(
    sp.diff(G_32_second, u) * b_32_second
    + sp.diff(G_32_second, v) * c_32_second
)
assert sp.factor(sp.rem(normal_numerator_32, v**2, v)) == (
    sp.Rational(15, 8) * q_parameter**2 * u**4 * v
)
for entry in (
    a_32_second,
    b_32_second.subs(q_parameter, 0),
    c_32_second.subs(q_parameter, 0),
):
    assert sp.cancel(entry / v**2).as_numer_denom()[1] == 1

# Root type 4+1. The primitive locus is p1!=0, but the first-normal
# remainder is 60*p1^3*u^4*v^2 modulo v^4.
p0, p1, p2 = sp.symbols("p0 p1 p2")
F_41 = u * v**4
a_41 = v**2
G_41 = -v**2 * (p0 * v**2 + 2 * p1 * u**2 + 4 * p2 * u * v)
b_41 = p0 * v**2 - 2 * p1 * u**2
c_41 = v * (p1 * u + p2 * v)
verify_complete_tangent_family(
    F_41,
    a_41,
    G_41,
    b_41,
    c_41,
    (p0, p1, p2),
)
H_41 = -sp.cancel(
    (
        sp.diff(G_41, u) * b_41
        + sp.diff(G_41, v) * c_41
    )
    / a_41
)
assert sp.cancel(H_41).as_numer_denom()[1] == 1
first_normal_41 = sp.expand(
    2
    * a_41
    * (
        b_41 * sp.diff(H_41, u)
        + c_41 * sp.diff(H_41, v)
    )
    + (
        sp.Matrix([b_41, c_41]).T
        * sp.hessian(G_41, (u, v))
        * sp.Matrix([b_41, c_41])
    )[0]
)
assert sp.factor(sp.rem(first_normal_41, v**4, v)) == (
    60 * p1**3 * u**4 * v**2
)

# Exceptional 2+1+1+1 orbit.
F_exceptional = v**2 * (5 * u**3 + 30 * u * v**2 + 8 * v**3)
assert sp.factor(sp.hessian(F_exceptional, (u, v)).det()) == (
    -600 * v**2 * (u - 2 * v) ** 2 * (u**2 + 4 * u * v + 6 * v**2)
)
a_exceptional = v * (u - 2 * v)
G_exceptional = (
    5
    * v
    * (
        3 * p * u**2 * v
        + 6 * p * v**3
        + 2 * q_parameter * u**3
        + 24 * q_parameter * u * v**2
        + 8 * q_parameter * v**3
    )
    / 2
)
b_exceptional = -p * a_exceptional / 2
c_exceptional = -q_parameter * a_exceptional / 2
verify_complete_tangent_family(
    F_exceptional,
    a_exceptional,
    G_exceptional,
    b_exceptional,
    c_exceptional,
    (p, q_parameter),
)
assert all(
    not sp.cancel(entry / a_exceptional).as_numer_denom()[1].has(u, v)
    for entry in (a_exceptional, b_exceptional, c_exceptional)
)

# Constant-kernel completion.
normal_x = sp.symbols("normal_x")
transverse_f_coefficients = sp.symbols("f0:6")
F_transverse = sum(
    transverse_f_coefficients[index]
    * u ** (5 - index)
    * v**index
    for index in range(6)
)
tj0, tj1, tj2, tk0, tk1, tc = sp.symbols(
    "tj0 tj1 tj2 tk0 tk1 tc"
)
J_transverse = tj0 * u**2 + tj1 * u * v + tj2 * v**2
K_transverse = tk0 * u + tk1 * v
h_transverse = (
    F_transverse
    + normal_x**3 * J_transverse / 6
    + normal_x**4 * K_transverse / 24
    + normal_x**5 * tc / 120
)
det_transverse = sp.expand(
    sp.hessian(h_transverse, (normal_x, u, v)).det()
)
binary_hessian_transverse = sp.hessian(F_transverse, (u, v)).det()
assert sp.expand(
    det_transverse.coeff(normal_x, 1)
    - J_transverse * binary_hessian_transverse
) == 0
transverse_without_J = det_transverse.subs({tj0: 0, tj1: 0, tj2: 0})
assert sp.expand(
    transverse_without_J.coeff(normal_x, 2)
    - K_transverse * binary_hessian_transverse / 2
) == 0
transverse_without_JK = transverse_without_J.subs({tk0: 0, tk1: 0})
assert sp.expand(
    transverse_without_JK.coeff(normal_x, 3)
    - tc * binary_hessian_transverse / 6
) == 0

# Tangent constant kernel, normalized to partial_v along normal_x=0.
alpha, beta = sp.symbols("alpha beta")
hh0, hh1, hh2, hh3 = sp.symbols("hh0:4")
jj0, jj1, jj2 = sp.symbols("jj0:3")
kk0, kk1, scalar_t = sp.symbols("kk0 kk1 scalar_t")
H_tangent = (
    hh0 * u**3 + hh1 * u**2 * v + hh2 * u * v**2 + hh3 * v**3
)
J_tangent = jj0 * u**2 + jj1 * u * v + jj2 * v**2
K_tangent = kk0 * u + kk1 * v
h_tangent = (
    alpha * u**5
    + normal_x * beta * u**4
    + normal_x**2 * H_tangent / 2
    + normal_x**3 * J_tangent / 6
    + normal_x**4 * K_tangent / 24
    + normal_x**5 * scalar_t / 120
)
det_tangent = sp.expand(
    sp.hessian(h_tangent, (normal_x, u, v)).det()
)

# alpha!=0: the x^2 face kills hh3, hh2, hh1 successively.
tangent_x2 = det_tangent.coeff(normal_x, 2)
assert sp.Poly(tangent_x2, u, v).coeff_monomial(u**3 * v**4) == (
    -120 * alpha * hh3**2
)
assert sp.Poly(
    tangent_x2.subs(hh3, 0),
    u,
    v,
).coeff_monomial(u**5 * v**2) == -60 * alpha * hh2**2
assert sp.Poly(
    tangent_x2.subs({hh3: 0, hh2: 0}),
    u,
    v,
).coeff_monomial(u**7) == -20 * alpha * hh1**2

case_alpha = {hh1: 0, hh2: 0, hh3: 0}
assert sp.factor(
    det_tangent.coeff(normal_x, 3).subs(case_alpha)
) == (
    sp.Rational(4, 3)
    * jj2
    * u**6
    * (5 * alpha * hh0 - 4 * beta**2)
)
assert sp.factor(
    det_tangent.coeff(normal_x, 4).subs({**case_alpha, jj2: 0})
) == -5 * alpha * jj1**2 * u**5
assert sp.factor(
    det_tangent.coeff(normal_x, 6).subs(
        {**case_alpha, jj2: 0, jj1: 0}
    )
) == -sp.Rational(5, 9) * alpha * kk1**2 * u**3
assert det_tangent.coeff(normal_x, 7).subs(
    {**case_alpha, jj2: 0, jj1: 0, kk1: 0}
) == 0

# alpha=0: the exact-seven coefficient is either zero or a square residual.
case_beta_initial = {alpha: 0, hh2: 0, hh3: 0}
assert sp.factor(
    det_tangent.coeff(normal_x, 3).subs(case_beta_initial)
) == (
    -sp.Rational(4, 3)
    * beta
    * u**6
    * (4 * beta * jj2 + 3 * hh1**2)
)
jj2_solution = -3 * hh1**2 / (4 * beta)
case_beta_jj2 = {**case_beta_initial, jj2: jj2_solution}
assert sp.Poly(
    sp.expand(det_tangent.coeff(normal_x, 4).subs(case_beta_jj2)),
    u,
    v,
).coeff_monomial(u**4 * v) == 15 * hh1**3
case_beta = {
    alpha: 0,
    hh1: 0,
    hh2: 0,
    hh3: 0,
    jj2: 0,
}
assert sp.factor(
    det_tangent.coeff(normal_x, 5).subs(case_beta)
) == -sp.Rational(7, 3) * beta * jj1**2 * u**4
assert sp.factor(
    det_tangent.coeff(normal_x, 7).subs({**case_beta, jj1: 0})
) == -sp.Rational(1, 3) * beta * kk1**2 * u**2


# Conic calibration F=u*v^4, G=u^2*v^2, H=-u^3.
F = u * v**4
G = u**2 * v**2
H = -u**3
C0 = sp.Matrix(
    [
        [H, sp.diff(G, u), sp.diff(G, v)],
        [sp.diff(G, u), sp.diff(F, u, 2), sp.diff(F, u, v)],
        [sp.diff(G, v), sp.diff(F, u, v), sp.diff(F, v, 2)],
    ]
)
kernel = sp.Matrix([-2 * v**2, -2 * u**2, u * v])
pairing_factor = -4 * v**2

assert C0 * kernel == sp.zeros(3, 1)
assert C0.adjugate().applyfunc(sp.expand) == (
    pairing_factor * kernel * kernel.T
).applyfunc(sp.expand)
assert sp.factor(sp.hessian(F, (u, v)).det()) == sp.factor(
    pairing_factor * kernel[0] ** 2
)


# The kernel entries span all of H^0(P^1,O(2)).
kernel_coefficient_matrix = sp.Matrix(
    [
        [
            sp.Poly(entry, u, v).coeff_monomial(monomial)
            for monomial in (u**2, u * v, v**2)
        ]
        for entry in kernel
    ]
)
assert kernel_coefficient_matrix.rank() == 3


# This calibration fails the first-normal divisibility by a^2.
a = kernel[0]
tangent_kernel = sp.Matrix(kernel[1:3, 0])
first_normal_numerator = sp.factor(
    2
    * a
    * tangent_kernel.dot(sp.Matrix([sp.diff(H, u), sp.diff(H, v)]))
    + (tangent_kernel.T * sp.hessian(G, (u, v)) * tangent_kernel)[0]
)
assert first_normal_numerator == -30 * u**4 * v**2
assert sp.cancel(first_normal_numerator / a**2).as_numer_denom()[1] != 1


# Verify the Jacobi coefficient formula for an arbitrary quadratic J.
x = sp.symbols("x")
j0, j1, j2 = sp.symbols("j0 j1 j2")
J = j0 * u**2 + j1 * u * v + j2 * v**2
h_first_normal = F + x * G + x**2 * H / 2 + x**3 * J / 6
determinant = sp.expand(sp.hessian(h_first_normal, (x, u, v)).det())
predicted_coefficient = sp.expand(
    pairing_factor * (a**2 * J + first_normal_numerator)
)
assert sp.expand(determinant.coeff(x, 1) - predicted_coefficient) == 0


print("PASS: septuple-line pole/defect ladder has exactly six rows")
print("PASS: Hessian integrability excludes the extremal pencil kernel")
print("PASS: both defect-four moving-linear-kernel rows are empty")
print("PASS: repeated binary-quintic square-Hessian strata classified")
print("PASS: squarefree binary square-Hessian locus is the Fermat orbit")
print("PASS: primitive conic kernel excludes the squarefree Fermat orbit")
print("PASS: extremal pole-three/defect-two septuple-line row is empty")
print("PASS: all three constant-kernel defect-six rows are empty")
print("PASS: conic square-Hessian and first-normal gates verified exactly")
print(
    "SCOPE: exact septuple line with squarefree quadratic cofactor is empty; "
    "multiplicity eight/nine are continued in HC4NHM3; lower-Smith rows "
    "remain open"
)
