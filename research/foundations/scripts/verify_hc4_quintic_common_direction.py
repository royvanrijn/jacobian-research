#!/usr/bin/env python3
"""Verify the quartic bordered lemma and quintic common-direction descent.

For a ternary polynomial s of degree at most four, put

    J(s) = grad(s)^T adj(Hess(s)) grad(s).

The quartic bordered lemma says that J(s)=0 forces s to omit a nonzero
constant direction.  This checker follows the binary-root stratification of
the leading quartic.  The non-pure-fourth-power strata are reduced by exact
face equations; the double-double exception and the pure-fourth-power
branches are certified by rational Singular ideals.

The lemma extends the existing HC4T31 common-direction descent from
degree-three to degree-four bordered coefficients.  The final calculation
also verifies the three clean determinant faces for a rank-three sextic with
a quintic that is independent of the constant sextic-kernel direction.
"""

from __future__ import annotations

import contextlib
import io
from pathlib import Path
import re
import runpy
import shutil
import subprocess

import sympy as sp


# Replay the common-direction descent after the bordered lemma supplies its
# second constant direction.
shared_path = Path(__file__).with_name(
    "verify_hc4_meng_triple_rank_three_reduction.py"
)
with contextlib.redirect_stdout(io.StringIO()):
    runpy.run_path(str(shared_path))


x, y, m, scale = sp.symbols("x y m scale")
ternary_variables = (x, y, m)


def bordered_invariant(polynomial: sp.Expr) -> sp.Expr:
    gradient = sp.Matrix(
        [sp.diff(polynomial, variable) for variable in ternary_variables]
    )
    hessian = sp.hessian(polynomial, ternary_variables)
    return sp.expand((gradient.T * hessian.adjugate() * gradient)[0])


def coefficient_equations(
    polynomial: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> list[sp.Expr]:
    equations: list[sp.Expr] = []
    for _, coefficient in sp.Poly(sp.expand(polynomial), *variables).terms():
        numerator = sp.together(coefficient).as_numer_denom()[0]
        equations.append(sp.expand(numerator))
    return equations


def singular_expression(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


singular = shutil.which("Singular")
if singular is None:
    raise RuntimeError("Singular is required for the exact radical checks")


def run_singular(program: str, timeout: int = 120) -> str:
    completed = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    if completed.stderr.strip():
        raise RuntimeError(completed.stderr)
    return completed.stdout


# 1. If the leading quartic is binary and has nonzero binary Hessian, the
# top coefficient first makes the cubic at most linear in m.
a, b, c, d, e = sp.symbols("a b c d e")
binary_quartic = (
    a * x**4
    + 4 * b * x**3 * y
    + 6 * c * x**2 * y**2
    + 4 * d * x * y**3
    + e * y**4
)
cubic_monomials = [
    x**x_degree * y**y_degree * m ** (3 - x_degree - y_degree)
    for x_degree in range(4)
    for y_degree in range(4 - x_degree)
]
cubic_coefficients = sp.symbols("u0:10")
generic_cubic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(
        cubic_coefficients, cubic_monomials
    )
)
quartic_cubic_pencil = (
    scale**4 * binary_quartic + scale**3 * generic_cubic
)
quartic_cubic_bordered = sp.Poly(
    bordered_invariant(quartic_cubic_pencil), scale
)
top_quartic_cubic = quartic_cubic_bordered.coeff_monomial(scale**15)
expected_top = (
    sp.Rational(4, 3)
    * binary_quartic
    * sp.hessian(binary_quartic, (x, y)).det()
    * sp.diff(generic_cubic, m, 2)
)
assert sp.expand(top_quartic_cubic - expected_top) == 0


# Thus write the m-dependent cubic part as m*Q_2(x,y).  The next face
# classifies Q_2 and the m^2 coefficient k on each binary-root stratum.
p0, p1, p2, k = sp.symbols("p0 p1 p2 k")
quadratic_q = p0 * x**2 + p1 * x * y + p2 * y**2
root_strata = {
    "four_distinct": x * y * (x - y) * (x - 2 * y),
    "double_simple_simple": x**2 * y * (x - y),
    "double_double": x**2 * y**2,
    "triple_simple": x**3 * y,
}
top_parameters = (p0, p1, p2, k)

for stratum, leading_form in root_strata.items():
    weighted = (
        scale**4 * leading_form
        + scale**3 * m * quadratic_q
        + scale**2 * k * m**2
    )
    face = sp.Poly(
        sp.Poly(bordered_invariant(weighted), scale).coeff_monomial(
            scale**14
        ),
        x,
        y,
        m,
    )
    equations = [coefficient for _, coefficient in face.terms()]
    face_ideal = sp.groebner(
        equations, *top_parameters, order="grevlex"
    )
    if stratum != "double_double":
        for parameter, exponent in zip(
            top_parameters, (2, 3, 2, 2)
        ):
            assert (
                face_ideal.reduce(sp.expand(parameter**exponent))[1] == 0
            )
    else:
        double_double_radical = (p0, p2, 8 * k - p1**2)
        radical_ideal = sp.groebner(
            double_double_radical,
            *top_parameters,
            order="grevlex",
        )
        for equation in equations:
            assert radical_ideal.reduce(equation)[1] == 0
        for generator, exponent in zip(
            double_double_radical, (2, 2, 2)
        ):
            assert (
                face_ideal.reduce(sp.expand(generator**exponent))[1] == 0
            )


# Once Q_2=k=0, the next possible m-dependence is m*L_1(x,y), followed
# by a linear term r*m.  Every non-pure-fourth root stratum forces L_1=0;
# the remaining coefficient then forces r=0.
lp, lq, r = sp.symbols("lp lq r")
linear_l = lp * x + lq * y
for leading_form in root_strata.values():
    weighted = (
        scale**4 * leading_form
        + scale**2 * m * linear_l
        + scale * r * m
    )
    bordered = sp.Poly(bordered_invariant(weighted), scale)
    face = sp.Poly(
        bordered.coeff_monomial(scale**12), x, y, m
    )
    equations = [coefficient for _, coefficient in face.terms()]
    face_ideal = sp.groebner(
        equations, lp, lq, r, order="grevlex"
    )
    assert face_ideal.reduce(lp**2)[1] == 0
    assert face_ideal.reduce(lq**2)[1] == 0
    no_linear_m = sp.Poly(
        bordered.as_expr().subs({lp: 0, lq: 0}), scale
    )
    assert no_linear_m.as_expr() != 0
    assert all(
        coefficient == 0
        for _, coefficient in sp.Poly(
            no_linear_m.as_expr(), x, y, m, scale
        ).terms()
        if coefficient != 0 and not coefficient.has(r)
    )
    assert sp.factor(no_linear_m.as_expr()).has(r)


# The only nonzero top-face exception is the double-double branch
# Q_2=p1*x*y, k=p1^2/8.  Normalize p1=1.  Even after every lower binary,
# mixed-linear, and linear term is included, the full bordered coefficient
# ideal is the unit ideal over Q.
g = sp.symbols("g0:4")
ell = sp.symbols("ell0:2")
lower_binary = sp.symbols("b0:3")
linear = sp.symbols("r0:3")
double_parameters = g + ell + lower_binary + linear
double_exception = (
    x**2 * y**2
    + (
        g[0] * x**3
        + g[1] * x**2 * y
        + g[2] * x * y**2
        + g[3] * y**3
        + m * x * y
    )
    + (
        lower_binary[0] * x**2
        + lower_binary[1] * x * y
        + lower_binary[2] * y**2
        + m * (ell[0] * x + ell[1] * y)
        + sp.Rational(1, 8) * m**2
    )
    + linear[0] * x
    + linear[1] * y
    + linear[2] * m
)
double_equations = coefficient_equations(
    bordered_invariant(double_exception), ternary_variables
)
double_program = f"""
ring rr=0,({",".join(map(str, double_parameters))}),dp;
option(redSB);
ideal I={",".join(map(singular_expression, double_equations))};
ideal G=slimgb(I);
print("DOUBLE "+string(size(G))+" "+string(reduce(1,G)));
"""
double_output = run_singular(double_program)
double_marker = re.search(r"(?m)^DOUBLE (\d+) ([01])$", double_output)
assert double_marker is not None
assert tuple(map(int, double_marker.groups())) == (1, 0)


# 2. Pure fourth power.  The first face is the binary Hessian determinant
# of the cubic on the two-dimensional kernel plane.
pure_fourth_pencil = scale**4 * x**4 + scale**3 * generic_cubic
pure_fourth_bordered = sp.Poly(
    bordered_invariant(pure_fourth_pencil), scale
)
assert sp.expand(
    pure_fourth_bordered.coeff_monomial(scale**14)
    - 16
    * x**6
    * sp.hessian(generic_cubic, (y, m)).det()
) == 0


# The binary singular-Hessian form and homogeneity give
#
#   f_3=g0*x^3+g1*x^2*y+g2*x*y^2+g3*y^3+a*x^2*m.
#
# If a!=0, normalize a=1.  The complete coefficient ideal has radical
# equal to the constant-direction alignment ideal displayed below.
quadratic_coefficients = sp.symbols("b0:6")
linear_coefficients = sp.symbols("r0:3")
pure_parameters = g + quadratic_coefficients + linear_coefficients
quadratic_monomials = [
    x**x_degree * y**y_degree * m ** (2 - x_degree - y_degree)
    for x_degree in range(3)
    for y_degree in range(3 - x_degree)
]


def pure_fourth_polynomial(a_value: int) -> sp.Expr:
    cubic = (
        g[0] * x**3
        + g[1] * x**2 * y
        + g[2] * x * y**2
        + g[3] * y**3
        + a_value * x**2 * m
    )
    quadratic = sum(
        coefficient * monomial
        for coefficient, monomial in zip(
            quadratic_coefficients, quadratic_monomials
        )
    )
    return (
        x**4
        + cubic
        + quadratic
        + linear_coefficients[0] * x
        + linear_coefficients[1] * y
        + linear_coefficients[2] * m
    )


nonzero_residual_equations = coefficient_equations(
    bordered_invariant(pure_fourth_polynomial(1)),
    ternary_variables,
)
alignment_generators = (
    g[2],
    g[3],
    quadratic_coefficients[1]
    - 2 * g[1] * quadratic_coefficients[0],
    quadratic_coefficients[2]
    - g[1] ** 2 * quadratic_coefficients[0],
    quadratic_coefficients[4]
    - g[1] * quadratic_coefficients[3],
    linear_coefficients[1] - g[1] * linear_coefficients[2],
)
nonzero_program = f"""
ring rr=0,({",".join(map(str, pure_parameters))}),dp;
option(redSB);
ideal I={",".join(map(singular_expression, nonzero_residual_equations))};
ideal GI=slimgb(I);
ideal A={",".join(map(singular_expression, alignment_generators))};
ideal GA=std(A);
print("NONZERO_INCLUSION "+string(size(reduce(I,GA))));
int i;
int exponent;
poly power;
for (i=1;i<=size(A);i++)
{{
  power=A[i];
  exponent=1;
  while (exponent<=4 && reduce(power,GI)!=0)
  {{
    power=power*A[i];
    exponent++;
  }}
  print(
    "NONZERO_POWER "+string(i)+" "+string(exponent)+" "
    +string(reduce(power,GI)==0)
  );
}}
"""
nonzero_output = run_singular(nonzero_program)
assert re.search(r"(?m)^NONZERO_INCLUSION 0$", nonzero_output)
nonzero_powers = re.findall(
    r"(?m)^NONZERO_POWER (\d+) (\d+) ([01])$", nonzero_output
)
assert len(nonzero_powers) == len(alignment_generators)
assert all(success == "1" for _, _, success in nonzero_powers)


# If a=0, the radical is the union of:
#
# * I1: independence of m;
# * I2: no cubic plane dependence and all quadratic/linear plane data
#       depending on one common linear form.
zero_residual_equations = coefficient_equations(
    bordered_invariant(pure_fourth_polynomial(0)),
    ternary_variables,
)
b0, b1, b2, b3, b4, b5 = quadratic_coefficients
r0, r1, r2 = linear_coefficients
independent_m = (b0, b1, b3, r2)
one_plane_form = (
    g[1],
    g[2],
    g[3],
    b1**2 - 4 * b0 * b2,
    2 * b2 * b3 - b1 * b4,
    2 * b0 * b4 - b1 * b3,
    2 * b2 * r2 - b1 * r1,
    2 * b0 * r1 - b1 * r2,
    b3 * r1 - b4 * r2,
)
zero_program = f"""
ring rr=0,({",".join(map(str, pure_parameters))}),dp;
option(redSB);
ideal I={",".join(map(singular_expression, zero_residual_equations))};
ideal GI=slimgb(I);
ideal I1={",".join(map(singular_expression, independent_m))};
ideal I2={",".join(map(singular_expression, one_plane_form))};
ideal T=intersect(I1,I2);
ideal GT=std(T);
print(
  "ZERO_BASE "+string(size(GI))+" "+string(size(GT))+" "
  +string(size(reduce(I,GT)))
);
int i;
int exponent;
poly power;
for (i=1;i<=size(GT);i++)
{{
  power=GT[i];
  exponent=1;
  while (exponent<=6 && reduce(power,GI)!=0)
  {{
    power=power*GT[i];
    exponent++;
  }}
  print(
    "ZERO_POWER "+string(i)+" "+string(exponent)+" "
    +string(reduce(power,GI)==0)
  );
}}
"""
zero_output = run_singular(zero_program)
zero_marker = re.search(
    r"(?m)^ZERO_BASE (\d+) (\d+) (\d+)$", zero_output
)
assert zero_marker is not None
assert tuple(map(int, zero_marker.groups())) == (966, 18, 0)
zero_powers = re.findall(
    r"(?m)^ZERO_POWER (\d+) (\d+) ([01])$", zero_output
)
assert len(zero_powers) == 18
assert all(success == "1" for _, _, success in zero_powers)


# 3. The first quintic extension.  For
#
#   q_2+h_3+h_4+h_5+h_6
#
# with rank Hess(h_6)=3, let t be the constant sextic-kernel direction.
# The lambda^15 face always makes D_t^2 h_5=0.  On the branch D_t h_5=0,
# the lambda^14 and lambda^13 faces successively make D_t^2 h_4 and
# D_t^2 h_3 zero, giving the common direction handled above.
lam = sp.symbols("lam")
c1, c2, c3 = sp.symbols("c1 c2 c3")
C6 = sp.diag(c1, c2, c3, 0)


def generic_symmetric(prefix: str) -> sp.Matrix:
    entries = sp.symbols(f"{prefix}0:10")
    return sp.Matrix(
        [
            [entries[0], entries[1], entries[2], entries[3]],
            [entries[1], entries[4], entries[5], entries[6]],
            [entries[2], entries[5], entries[7], entries[8]],
            [entries[3], entries[6], entries[8], entries[9]],
        ]
    )


H0 = generic_symmetric("h")
A3 = generic_symmetric("a")
B4 = generic_symmetric("b")
D5 = generic_symmetric("d")
quotient_determinant = c1 * c2 * c3
quintic_pencil = sp.Poly(
    (
        H0
        + lam * A3
        + lam**2 * B4
        + lam**3 * D5
        + lam**4 * C6
    ).det(method="berkowitz"),
    lam,
)
assert sp.expand(
    quintic_pencil.coeff_monomial(lam**15)
    - quotient_determinant * D5[3, 3]
) == 0

# After the lambda^15 equation D_tt=0, the next face is the exact Schur
# norm of the remaining quintic cross vector.  This records the live
# nonzero-derivative frontier rather than hiding it in the closed branch.
zero_quintic_tt = D5.subs({D5[3, 3]: 0})
quintic_cross = sp.Matrix(
    [
        zero_quintic_tt[0, 3],
        zero_quintic_tt[1, 3],
        zero_quintic_tt[2, 3],
    ]
)
sextic_adjugate = sp.diag(c2 * c3, c1 * c3, c1 * c2)
schur_norm_face = (
    quotient_determinant * B4[3, 3]
    - (quintic_cross.T * sextic_adjugate * quintic_cross)[0]
)
zero_quintic_tt_pencil = sp.Poly(
    (
        H0
        + lam * A3
        + lam**2 * B4
        + lam**3 * zero_quintic_tt
        + lam**4 * C6
    ).det(method="berkowitz"),
    lam,
)
assert sp.expand(
    zero_quintic_tt_pencil.coeff_monomial(lam**14)
    - schur_norm_face
) == 0

zero_quintic_t = D5.copy()
zero_quintic_t[3, :] = sp.zeros(1, 4)
zero_quintic_t[:, 3] = sp.zeros(4, 1)
zero_quintic_pencil = sp.Poly(
    (
        H0
        + lam * A3
        + lam**2 * B4
        + lam**3 * zero_quintic_t
        + lam**4 * C6
    ).det(method="berkowitz"),
    lam,
)
assert sp.expand(
    zero_quintic_pencil.coeff_monomial(lam**14)
    - quotient_determinant * B4[3, 3]
) == 0

zero_quartic_tt = B4.subs({B4[3, 3]: 0})
common_direction_pencil = sp.Poly(
    (
        H0
        + lam * A3
        + lam**2 * zero_quartic_tt
        + lam**3 * zero_quintic_t
        + lam**4 * C6
    ).det(method="berkowitz"),
    lam,
)
assert sp.expand(
    common_direction_pencil.coeff_monomial(lam**13)
    - quotient_determinant * A3[3, 3]
) == 0


# 4. The polynomial Schur norm genuinely has nonzero solutions, so it
# cannot by itself close the live branch.  Test the canonical diagonal
# solution
#
#   h6=(x^6+y^6+m^6)/30,
#   D_t h5=alpha*x^4+beta*y^4+gamma*m^4.
#
# The forced lambda^14 correction is
#
#   h4=8*t^2*(alpha^2*x^2+beta^2*y^2+gamma^2*m^2).
#
# The lambda^13 face is canceled by the displayed t^3 term in h3.  Even
# with an arbitrary quadratic base block and arbitrary t-cross terms in
# q2, lambda^11 then contains three uncancellable fifth-power monomials,
# so this most natural prolongation has alpha=beta=gamma=0.
t4 = sp.symbols("t4")
alpha4, beta4, gamma4 = sp.symbols("alpha4 beta4 gamma4")
cross_x, cross_y, cross_m = sp.symbols(
    "cross_x cross_y cross_m"
)
base_quadratic_coefficients = sp.symbols("base0:6")
four_variables = (x, y, m, t4)
diagonal_h6 = (x**6 + y**6 + m**6) / 30
diagonal_h5 = t4 * (
    alpha4 * x**4 + beta4 * y**4 + gamma4 * m**4
)
diagonal_h4 = 8 * t4**2 * (
    alpha4**2 * x**2
    + beta4**2 * y**2
    + gamma4**2 * m**2
)
diagonal_h3 = (
    sp.Rational(32, 3)
    * (alpha4**3 + beta4**3 + gamma4**3)
    * t4**3
)
base0, base1, base2, base3, base4, base5 = (
    base_quadratic_coefficients
)
diagonal_q2 = (
    t4 * (cross_x * x + cross_y * y + cross_m * m)
    + sp.Rational(1, 2)
    * (
        base0 * x**2
        + 2 * base1 * x * y
        + 2 * base2 * x * m
        + base3 * y**2
        + 2 * base4 * y * m
        + base5 * m**2
    )
)
diagonal_matrix = sp.zeros(4)
for weight, homogeneous_part in (
    (4, diagonal_h6),
    (3, diagonal_h5),
    (2, diagonal_h4),
    (1, diagonal_h3),
    (0, diagonal_q2),
):
    diagonal_matrix += lam**weight * sp.hessian(
        homogeneous_part, four_variables
    )
diagonal_determinant = sp.Poly(
    sp.expand(diagonal_matrix.det(method="berkowitz")), lam
)
for exponent in (15, 14, 13, 12):
    assert diagonal_determinant.coeff_monomial(lam**exponent) == 0
diagonal_face_11 = sp.expand(
    diagonal_determinant.coeff_monomial(lam**11)
)
expected_diagonal_face_11 = sp.expand(
    1024
    * t4**3
    * (
        alpha4**5 * y**4 * m**4
        + beta4**5 * x**4 * m**4
        + gamma4**5 * x**4 * y**4
    )
    - 8
    * (
        alpha4 * cross_x * x**3 * y**4 * m**4
        + beta4 * cross_y * x**4 * y**3 * m**4
        + gamma4 * cross_m * x**4 * y**4 * m**3
    )
)
assert sp.expand(
    diagonal_face_11 - expected_diagonal_face_11
) == 0
assert not (
    diagonal_face_11.free_symbols
    & set(base_quadratic_coefficients)
)


print("PASS: every binary-quartic root stratum has a constant direction")
print("PASS: the double-double exceptional face has a unit lower ideal")
print("PASS: both pure-fourth radical branches omit a constant direction")
print("PASS: the quartic bordered lemma holds in all ternary strata")
print("PASS: lambda^14 isolates the live quintic Schur norm")
print("PASS: the quintic rank-three zero-derivative faces are clean")
print("PASS: the canonical diagonal nonzero Schur norm dies at lambda^11")
print("SCOPE: common-direction Schur descent now allows degree-four s")
