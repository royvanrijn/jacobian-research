#!/usr/bin/env python3
"""Exact and bounded audit of the sharp degree-sixteen A4 selector plane.

The exact part constructs the birational degree-ten strict model of

    P<q0, q1, q3>,  q0=A^3,

on the rational (U,V) root chart and computes the complete rational
fixed-infinity tangent-cone hierarchy.  Singular verifies generic absolute
irreducibility, the genera and local delta invariants of the terminal fixed
degenerations, the first two exact rational points on the moving
discriminant, and two bounded rational-parameter censuses.

The height censuses are experiments, not proofs about every rational point
of the parameter plane.  The degree-twenty-three moving critical curve is
not fully implicitized in parameter space, and this script does not
construct a Keller map.
"""

from __future__ import annotations

from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from math import gcd
import shutil
import subprocess

import sympy as sp


# ---------------------------------------------------------------------------
# 1. The A4 root algebra and the sharp selector plane
# ---------------------------------------------------------------------------

a, b, T = sp.symbols("a b T")

jly_B = (
    a**3
    - 3 * a * b**2
    + 2 * b**3
    - 9 * a * b
    + 9 * b**2
    - 27 * a
    + 27 * b
    + 27
)
jly_A = a**3 - b**3 - 9 * b**2 - 27 * b - 54
jly_C = a**3 - b**3 + 27
quartic = (
    T**4
    - 6 * jly_A * jly_B * T**2
    - 8 * jly_B**3 * T
    + jly_B**2 * (9 * jly_A**2 - 12 * jly_C * jly_B)
)

q0 = a**3
q1_integral = (
    4 * b * T
    + 81 * a * b**2
    + 243 * a * b
    + 729 * a
    - 72 * b**3
    - 324 * b**2
    - 972 * b
    - 972
)
q3_integral = (
    3 * a * T
    + 4 * T
    - 54 * a * b**2
    - 162 * a * b
    - 486 * a
    + 12 * b**3
    - 324
)

# The basis used below is (q0, q1, q3), where the integral polynomials are
# q1_integral=4*q1 and q3_integral=3*q3.
assert sp.Poly(q1_integral, T).degree() == 1
assert sp.Poly(q3_integral, T).degree() == 1


# ---------------------------------------------------------------------------
# 2. Pullback to the rational root chart
# ---------------------------------------------------------------------------

U, V, Z = sp.symbols("U V Z")

H = (
    8 * U**3
    - 6 * U * V**2
    - 18 * U * V
    - 54 * U
    - 2 * V**3
    - 9 * V**2
    - 27 * V
    - 27
)
K = 4 * U**2 + 4 * U * V + 6 * U + V**2 + 3 * V + 9
M = U**2 + 2 * V**2 + 6 * V + 18
L = (
    U**3
    - 3 * U * V**2
    - 9 * U * V
    - 27 * U
    + 2 * V**3
    + 9 * V**2
    + 27 * V
    + 27
)
N1 = sp.expand(M * K)
N2 = (
    8 * U**3 * V
    + 12 * U**2 * V**2
    + 36 * U**2 * V
    + 108 * U**2
    + 6 * U * V**3
    + 36 * U * V**2
    + 108 * U * V
    + 162 * U
    + V**4
    + 9 * V**3
    + 27 * V**2
    + 54 * V
)
source_A = U**3 - V**3 - 9 * V**2 - 27 * V - 54
root_numerator = 3 * source_A * K**3 * L

# On this chart a=N1/H, b=N2/H, and T=root_numerator/H^3.  Clearing the
# common H^4 denominator gives the three numerators below.
pullback_q0 = sp.expand(N1**3 * H)
pullback_q1 = sp.expand(
    4 * N2 * root_numerator
    + (
        81 * N1 * N2**2
        + 243 * N1 * N2 * H
        + 729 * N1 * H**2
        - 72 * N2**3
        - 324 * N2**2 * H
        - 972 * N2 * H**2
        - 972 * H**3
    )
    * H
)
pullback_q3 = sp.expand(
    (3 * N1 + 4 * H) * root_numerator
    + (
        -54 * N1 * N2**2
        - 162 * N1 * N2 * H
        - 486 * N1 * H**2
        + 12 * N2**3
        - 324 * H**3
    )
    * H
)


def exact_quotient(numerator: sp.Expr, denominator: sp.Expr) -> sp.Expr:
    quotient = sp.cancel(numerator / denominator)
    assert sp.denom(quotient) == 1
    assert sp.expand(numerator - denominator * quotient) == 0
    return sp.expand(quotient)


strict_q0 = exact_quotient(pullback_q0, K**3)
strict_q1 = exact_quotient(pullback_q1, K**3)
strict_q3 = exact_quotient(pullback_q3, K**3)

assert sp.expand(strict_q0 - M**3 * H) == 0
assert tuple(
    sp.Poly(polynomial, U, V).total_degree()
    for polynomial in (strict_q0, strict_q1, strict_q3)
) == (9, 10, 10)
assert sp.gcd(sp.gcd(strict_q0, strict_q1), strict_q3) == 1

# Thus a selector A*q0+B*q1+C*q3 has strict degree-ten equation
#
#     12*A*strict_q0 + 3*B*strict_q1 + 4*C*strict_q3 = 0.
#
# The discarded K^3 is fixed and independent of [A:B:C].  The established
# root-chart field recovery makes this strict curve birational to the
# corresponding horizontal norm component whenever it is integral.


# ---------------------------------------------------------------------------
# 3. Complete fixed-infinity tangent-cone hierarchy
# ---------------------------------------------------------------------------


def homogenize_to_degree(expression: sp.Expr, degree: int) -> sp.Expr:
    polynomial = sp.Poly(expression, U, V)
    return sp.expand(
        sum(
            coefficient
            * U**exponents[0]
            * V**exponents[1]
            * Z ** (degree - sum(exponents))
            for exponents, coefficient in polynomial.terms()
        )
    )


def homogeneous_piece(
    expression: sp.Expr,
    first: sp.Symbol,
    second: sp.Symbol,
    degree: int,
) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(expression), first, second)
    return sp.expand(
        sum(
            coefficient * first**exponents[0] * second**exponents[1]
            for exponents, coefficient in polynomial.terms()
            if sum(exponents) == degree
        )
    )


projective_q0, projective_q1, projective_q3 = (
    homogenize_to_degree(polynomial, 10)
    for polynomial in (strict_q0, strict_q1, strict_q3)
)

x, y = sp.symbols("x y")
at_P = tuple(
    sp.expand(polynomial.subs({U: 1 + x, V: 1, Z: y}))
    for polynomial in (projective_q0, projective_q1, projective_q3)
)
at_Q = tuple(
    sp.expand(polynomial.subs({U: -sp.Rational(1, 2) + x, V: 1, Z: y}))
    for polynomial in (projective_q0, projective_q1, projective_q3)
)

A_parameter, B_parameter, C_parameter = sp.symbols("A_parameter B_parameter C_parameter")
selector_at_P = sp.expand(
    12 * A_parameter * at_P[0]
    + 3 * B_parameter * at_P[1]
    + 4 * C_parameter * at_P[2]
)
selector_at_Q = sp.expand(
    12 * A_parameter * at_Q[0]
    + 3 * B_parameter * at_Q[1]
    + 4 * C_parameter * at_Q[2]
)

p_quadratic = homogeneous_piece(selector_at_P, x, y, 2)
assert sp.expand(
    p_quadratic
    - 729
    * (4 * A_parameter + 9 * B_parameter - 56 * C_parameter)
    * y
    * (2 * x - 3 * y)
) == 0

# On the exceptional parameter line, the cubic is always a double-tangent
# cubic.  It vanishes only at the single rational point [65:-4:4].
p_line_substitution = {
    A_parameter: (-9 * B_parameter + 56 * C_parameter) / 4,
}
p_cubic_on_line = sp.factor(
    homogeneous_piece(
        selector_at_P.subs(p_line_substitution),
        x,
        y,
        3,
    )
)
assert sp.expand(
    p_cubic_on_line
    - 8748 * (B_parameter + C_parameter) * x * (x - 3 * y) ** 2
) == 0

deep_selector_at_P = selector_at_P.subs(
    {
        A_parameter: 65,
        B_parameter: -4,
        C_parameter: 4,
    }
)
assert homogeneous_piece(deep_selector_at_P, x, y, 2) == 0
assert homogeneous_piece(deep_selector_at_P, x, y, 3) == 0
deep_quartic = sp.factor(
    homogeneous_piece(deep_selector_at_P, x, y, 4)
)
assert deep_quartic == (
    2916
    * y
    * (2 * x - 3 * y)
    * (14 * x**2 - 90 * x * y + 135 * y**2)
)
assert sp.discriminant(14 * x**2 - 90 * x + 135, x) != 0

# At Q=[-1/2:1:0], C!=0 gives a nondegenerate quadratic tangent cone.
q_quadratic = sp.factor(homogeneous_piece(selector_at_Q, x, y, 2))
assert sp.expand(
    q_quadratic
    + sp.Rational(19683, 64)
    * C_parameter
    * (4 * x**2 + 6 * x * y + 9 * y**2)
) == 0
assert sp.discriminant(4 * x**2 + 6 * x + 9, x) != 0

# On C=0 the first tangent cone is cubic.  Its discriminant never vanishes
# at a rational exact selector because A^2+3AB+9B^2 is anisotropic over Q.
q_cubic_on_line = homogeneous_piece(
    selector_at_Q.subs(C_parameter, 0),
    x,
    y,
    3,
)
q_cubic_discriminant = sp.factor(
    sp.discriminant(q_cubic_on_line.subs(y, 1), x)
)
anisotropic_form = A_parameter**2 + 3 * A_parameter * B_parameter + 9 * B_parameter**2
discriminant_quotient = sp.factor(
    q_cubic_discriminant / anisotropic_form**2
)
assert discriminant_quotient != 0
assert not discriminant_quotient.free_symbols
assert sp.discriminant(anisotropic_form.subs(B_parameter, 1), A_parameter) == -27


# ---------------------------------------------------------------------------
# 4. Moving critical curve and its first rational selector components
# ---------------------------------------------------------------------------

# A member of the strict net is singular at (U,V) precisely when its
# coefficient vector lies in the kernel of the value-and-gradient matrix.
# The determinant has just two factors.  The conic M=0 is caused by
# R0=M^3*H; the genuine residual critical curve has degree twenty three.
selector_matrix = sp.Matrix(
    [
        [strict_q0, strict_q1, strict_q3],
        [
            sp.diff(strict_q0, U),
            sp.diff(strict_q1, U),
            sp.diff(strict_q3, U),
        ],
        [
            sp.diff(strict_q0, V),
            sp.diff(strict_q1, V),
            sp.diff(strict_q3, V),
        ],
    ]
)
critical_determinant = sp.expand(selector_matrix.det())
moving_critical = exact_quotient(critical_determinant, -54 * M**2)
assert sp.Poly(moving_critical, U, V).total_degree() == 23
critical_coefficient, critical_factors = sp.factor_list(moving_critical)
assert critical_coefficient == 1
assert len(critical_factors) == 1
assert critical_factors[0][1] == 1
assert sp.expand(critical_factors[0][0] - moving_critical) == 0

critical_top = homogeneous_piece(moving_critical, U, V, 23)
assert sp.expand(
    critical_top
    - 8
    * (U - V) ** 8
    * (U + 2 * V) ** 3
    * (2 * U + V) ** 6
    * (U**2 + 2 * V**2)
    * (U**2 + U * V + V**2) ** 2
) == 0

rho_source = V**2 + 3 * V + 9
assert sp.expand(
    moving_critical.subs(U, 0)
    - 8 * (2 * V + 3) ** 2 * (4 * V + 87) * rho_source**10
) == 0

kernel_coordinates = (
    exact_quotient(
        strict_q1 * sp.diff(strict_q3, U)
        - strict_q3 * sp.diff(strict_q1, U),
        9,
    ),
    exact_quotient(
        strict_q3 * sp.diff(strict_q0, U)
        - strict_q0 * sp.diff(strict_q3, U),
        9,
    ),
    exact_quotient(
        strict_q0 * sp.diff(strict_q1, U)
        - strict_q1 * sp.diff(strict_q0, U),
        9,
    ),
)
assert tuple(
    sp.Poly(coordinate, U, V).total_degree()
    for coordinate in kernel_coordinates
) == (19, 18, 18)
assert sp.gcd(
    sp.gcd(kernel_coordinates[0], kernel_coordinates[1]),
    kernel_coordinates[2],
) == 1

# The two rational critical points visible on U=0 give the first exact
# rational selectors on the moving discriminant.  Kernel coordinates below
# multiply (R0,R1,R3); converting to the original selector coordinates uses
# [lambda:mu:nu]=[s0:4*s1:3*s3].
point_77 = {U: 0, V: -sp.Rational(3, 2)}
point_103 = {U: 0, V: -sp.Rational(87, 4)}
kernel_77 = sp.Matrix([-sp.Rational(231, 8), sp.Rational(3, 2), 1])
kernel_103 = sp.Matrix([sp.Rational(309, 8), -sp.Rational(3, 2), 1])
matrix_77 = selector_matrix.subs(point_77)
matrix_103 = selector_matrix.subs(point_103)
assert matrix_77.rank() == 2
assert matrix_103.rank() == 2
assert matrix_77 * kernel_77 == sp.zeros(3, 1)
assert matrix_103 * kernel_103 == sp.zeros(3, 1)

curve_77 = sp.expand(924 * strict_q0 - 48 * strict_q1 - 32 * strict_q3)
curve_103 = sp.expand(1236 * strict_q0 - 48 * strict_q1 + 32 * strict_q3)
residual_9 = exact_quotient(curve_103, 36 * U)
assert sp.Poly(residual_9, U, V).total_degree() == 9
assert sp.gcd(U, residual_9) == 1

# The coefficient-plane norm of the [103:-16:8] selector has a rational
# conic factor.  The residual degree-fourteen factor corresponds
# birationally to residual_9 because the selector is linear in T.
selector_103 = sp.expand(
    1236 * q0 - 48 * q1_integral + 32 * q3_integral
)
norm_103 = sp.expand(sp.resultant(quartic, selector_103, T))
norm_coefficient, norm_factors = sp.factor_list(norm_103)
assert norm_coefficient == -2304
assert sorted(
    (sp.Poly(factor, a, b).total_degree(), exponent)
    for factor, exponent in norm_factors
) == [(2, 1), (14, 1)]
norm_conic = next(
    factor
    for factor, _ in norm_factors
    if sp.Poly(factor, a, b).total_degree() == 2
)
norm_residual_14 = next(
    factor
    for factor, _ in norm_factors
    if sp.Poly(factor, a, b).total_degree() == 14
)
assert sp.expand(norm_conic - (a**2 - 4 * b**2 - 12 * b - 36)) == 0
assert norm_conic.subs({a: 6, b: 0}) == 0
assert sp.gcd(norm_conic, norm_residual_14) == 1

# On M=0 the first column of selector_matrix vanishes.  The remaining two
# columns have generic rank two.  Their only rank-drop loci are the fixed
# quadratic cluster and one conjugate quadratic pair, verified below in
# Singular.  At the latter pair the kernel ratio s1/s3 depends on V, so the
# two conjugates have no common rational exact selector: their rational
# kernel intersection is again q0.
exceptional_v = sp.symbols("exceptional_v")
exceptional_u = (-8 * exceptional_v - 12) / 7
exceptional_minpoly = sp.Poly(
    3 * exceptional_v**2 + 9 * exceptional_v + 19,
    exceptional_v,
    domain=sp.QQ,
)


def reduce_at_exceptional_pair(expression: sp.Expr) -> sp.Expr:
    substituted = sp.together(
        expression.subs({U: exceptional_u, V: exceptional_v})
    )
    numerator, denominator = sp.fraction(substituted)
    numerator_polynomial = sp.Poly(numerator, exceptional_v, domain=sp.QQ)
    denominator_polynomial = sp.Poly(
        denominator,
        exceptional_v,
        domain=sp.QQ,
    )
    inverse = sp.invert(denominator_polynomial, exceptional_minpoly)
    return sp.rem(
        numerator_polynomial * inverse,
        exceptional_minpoly,
    ).as_expr()


assert reduce_at_exceptional_pair(strict_q1) == (
    sp.Rational(17915904, 7) * (2 * exceptional_v + 3)
)
assert reduce_at_exceptional_pair(strict_q3) == (
    -sp.Rational(5971968, 7) * (exceptional_v - 2)
)
exceptional_ratio = exceptional_v / 7 + sp.Rational(8, 21)
assert reduce_at_exceptional_pair(
    strict_q1 * exceptional_ratio + strict_q3
) == 0
assert sp.Poly(exceptional_ratio, exceptional_v).degree() == 1


# ---------------------------------------------------------------------------
# 5. Singular verification helpers
# ---------------------------------------------------------------------------


def singular_expression(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


R0_SINGULAR = singular_expression(strict_q0)
R1_SINGULAR = singular_expression(strict_q1)
R3_SINGULAR = singular_expression(strict_q3)
M_SINGULAR = singular_expression(M)
RESIDUAL9_SINGULAR = singular_expression(residual_9)


def singular_base(*, libraries: tuple[str, ...] = ()) -> str:
    headers = "\n".join(f'LIB "{library}";' for library in libraries)
    return f"""{headers}
ring r=0,(U,V),dp;
poly R0={R0_SINGULAR};
poly R1={R1_SINGULAR};
poly R3={R3_SINGULAR};
"""


SINGULAR = shutil.which("Singular")
if SINGULAR is None:
    raise SystemExit("Singular is required for the sharp-selector audit")


def run_singular(source: str, *, timeout: int = 120) -> str:
    completed = subprocess.run(
        [SINGULAR, "-q"],
        input=source,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )
    if completed.returncode != 0:
        raise AssertionError(
            "Singular failed:\n"
            + completed.stdout[-4000:]
            + completed.stderr[-4000:]
        )
    return completed.stdout


# Absolute factorization over the two rational function fields proves that
# the generic plane member (C!=0) and the generic C=0 pencil member are
# geometrically integral.
generic_absolute_source = f"""LIB "absfact.lib";
ring r=(0,l,m),(U,V),dp;
poly R0={R0_SINGULAR};
poly R1={R1_SINGULAR};
poly R3={R3_SINGULAR};
poly plane=l*R0+m*R1+R3;
def plane_ring=absFactorize(plane);
setring plane_ring;
if (absolute_factors[4]!=1) {{ ERROR("generic plane is not absolutely irreducible"); }}
setring r;
poly line=l*R0+R1;
def line_ring=absFactorize(line);
setring line_ring;
if (absolute_factors[4]!=1) {{ ERROR("generic line is not absolutely irreducible"); }}
print("GENERIC_ABSOLUTE_PASS");
"""
assert "GENERIC_ABSOLUTE_PASS" in run_singular(generic_absolute_source)


# Exact primary decomposition of the rank-drop locus inside M=0.  The
# first component is the universal fixed cluster.  The second is the
# conjugate quadratic pair analyzed symbolically above.
rank_drop_source = singular_base(libraries=("primdec.lib",)) + f"""
poly M={M_SINGULAR};
ideal rankdrop=std(ideal(
    M,
    R1*diff(R3,U)-R3*diff(R1,U),
    R1*diff(R3,V)-R3*diff(R1,V),
    diff(R1,U)*diff(R3,V)-diff(R3,U)*diff(R1,V)
));
list components=minAssGTZ(rankdrop);
if (size(components)!=2) {{ ERROR("unexpected M rank-drop decomposition"); }}
string encoded=string(components);
if (find(encoded,"V2+3V+9")==0 || find(encoded,"3V2+9V+19")==0)
{{
    ERROR("missing expected M rank-drop component");
}}
print("M_RANK_DROP_PASS");
"""
assert "M_RANK_DROP_PASS" in run_singular(rank_drop_source)


# The full projective genera and the local delta packages of the two
# terminal fixed-infinity degenerations are exact, not bounded evidence.
special_source = singular_base(libraries=("normal.lib",)) + r"""
proc checkDelta(poly f, ideal point, int delta, int tau, int branches)
{
    list data=deltaLoc(f,point);
    if (data[1]!=delta || data[2]!=tau || data[3]!=branches)
    {
        ERROR("unexpected local delta package");
    }
}

poly curve14=12*R0+3*R1+4*R3;       // [1:1:1]
poly curve12c=3*R1;                  // [0:1:0]
poly curve12p=396*R0+12*R1+12*R3;   // [33:4:3]
poly curve10=108*R0-12*R1;           // [9:-4:0]
poly curve9=780*R0-12*R1+16*R3;     // [65:-4:4]

if (genus(ideal(curve14))!=14) { ERROR("wrong genus at [1:1:1]"); }
if (genus(ideal(curve12c))!=12) { ERROR("wrong genus at [0:1:0]"); }
if (genus(ideal(curve12p))!=12) { ERROR("wrong genus at [33:4:3]"); }
if (genus(ideal(curve10))!=10) { ERROR("wrong genus at [9:-4:0]"); }
if (genus(ideal(curve9))!=9) { ERROR("wrong genus at [65:-4:4]"); }

ideal cluster=U,V^2+3*V+9;
checkDelta(curve10,cluster,20,32,6);
checkDelta(curve9,cluster,20,32,6);

ring projective_ring=0,(U,V,Z),dp;
poly R0=imap(r,R0);
poly R1=imap(r,R1);
poly R3=imap(r,R3);
poly curve14=homog(12*R0+3*R1+4*R3,Z);
poly curve12c=homog(3*R1,Z);
poly curve12p=homog(396*R0+12*R1+12*R3,Z);
poly curve10=homog(108*R0-12*R1,Z);
poly curve9=homog(780*R0-12*R1+16*R3,Z);

def old=projective_ring;
ring infinity_ring=0,(X,Y),dp;
map infinity_chart=old,X,1,Y;
poly curve14=infinity_chart(curve14);
poly curve12c=infinity_chart(curve12c);
poly curve12p=infinity_chart(curve12p);
poly curve10=infinity_chart(curve10);
poly curve9=infinity_chart(curve9);

ideal pointP=X-1,Y;
ideal pointQ=2*X+1,Y;
checkDelta(curve14,pointP,1,1,2);
checkDelta(curve14,pointQ,1,1,2);
checkDelta(curve12c,pointP,1,1,2);
checkDelta(curve12c,pointQ,3,4,3);
checkDelta(curve12p,pointP,3,5,2);
checkDelta(curve12p,pointQ,1,1,2);
checkDelta(curve10,pointP,3,5,2);
checkDelta(curve10,pointQ,3,4,3);
checkDelta(curve9,pointP,6,9,4);
checkDelta(curve9,pointQ,1,1,2);
print("SPECIAL_GENUS_DELTA_PASS");
"""
assert "SPECIAL_GENUS_DELTA_PASS" in run_singular(special_source)


# The first two rational points of the moving critical curve give one
# irreducible nodal member and one reducible member with a rational line.
moving_special_source = singular_base(libraries=("normal.lib",)) + f"""
proc checkDelta(poly f, ideal point, int delta, int tau, int branches)
{{
    list data=deltaLoc(f,point);
    if (data[1]!=delta || data[2]!=tau || data[3]!=branches)
    {{
        ERROR("unexpected moving local delta package");
    }}
}}

poly curve77=924*R0-48*R1-32*R3;
poly curve103=1236*R0-48*R1+32*R3;
poly residual9={RESIDUAL9_SINGULAR};
if (genus(ideal(curve77))!=12) {{ ERROR("wrong genus at [77:-16:-8]"); }}
if (genus(ideal(residual9))!=10) {{ ERROR("wrong residual genus at [103:-16:8]"); }}
checkDelta(curve77,ideal(2*U,2*V+3),1,1,2);
checkDelta(curve103,ideal(4*U,4*V+87),1,1,2);
checkDelta(curve77,ideal(U,V^2+3*V+9),20,32,6);

ring projective_ring=0,(U,V,Z),dp;
poly R0=imap(r,R0);
poly R1=imap(r,R1);
poly R3=imap(r,R3);
poly curve77=homog(924*R0-48*R1-32*R3,Z);
def old=projective_ring;
ring infinity_ring=0,(X,Y),dp;
map infinity_chart=old,X,1,Y;
poly curve77=infinity_chart(curve77);
checkDelta(curve77,ideal(X-1,Y),1,1,2);
checkDelta(curve77,ideal(2*X+1,Y),1,1,2);
checkDelta(curve77,ideal(X+2,Y),1,1,2);
print("MOVING_GENUS_DELTA_PASS");
"""
assert "MOVING_GENUS_DELTA_PASS" in run_singular(moving_special_source)

moving_absolute_source = singular_base(libraries=("absfact.lib",)) + f"""
poly curve77=924*R0-48*R1-32*R3;
poly residual9={RESIDUAL9_SINGULAR};
def curve77_ring=absFactorize(curve77);
setring curve77_ring;
if (absolute_factors[4]!=1) {{ ERROR("[77:-16:-8] is not absolutely irreducible"); }}
setring r;
def residual9_ring=absFactorize(residual9);
setring residual9_ring;
if (absolute_factors[4]!=1) {{ ERROR("degree-nine residual is not absolutely irreducible"); }}
print("MOVING_ABSOLUTE_PASS");
"""
assert "MOVING_ABSOLUTE_PASS" in run_singular(moving_absolute_source)


# ---------------------------------------------------------------------------
# 6. Bounded rational-parameter censuses
# ---------------------------------------------------------------------------


def projective_parameters(bound: int) -> tuple[tuple[int, int, int], ...]:
    parameters: list[tuple[int, int, int]] = []
    for first in range(-bound, bound + 1):
        for second in range(-bound, bound + 1):
            for third in range(-bound, bound + 1):
                # B=C=0 is the non-exact selector q0 and is excluded.
                if second == 0 and third == 0:
                    continue
                if gcd(gcd(abs(first), abs(second)), abs(third)) != 1:
                    continue
                leading = next(
                    value for value in (first, second, third) if value
                )
                if leading < 0:
                    continue
                parameters.append((first, second, third))
    return tuple(parameters)


height_six = projective_parameters(6)
height_three = projective_parameters(3)
assert len(height_six) == 864
assert len(height_three) == 144

# Absolute factorization is much faster than normalization, so the larger
# height-six box is used for the factorization census.  Add the two exact
# terminal parameters, which lie outside that box.
absolute_parameters = height_six + ((9, -4, 0), (65, -4, 4))
absolute_lines = [singular_base(libraries=("absfact.lib",))]
for index, (first, second, third) in enumerate(absolute_parameters):
    absolute_lines.append(
        f"poly curve{index}=({12 * first})*R0"
        f"+({3 * second})*R1+({4 * third})*R3;"
    )
    absolute_lines.append(
        f"def absolute_ring{index}=absFactorize(curve{index});"
        f" setring absolute_ring{index};"
        f" if (absolute_factors[4]!=1)"
        f" {{ ERROR(\"absolute factorization failed at index {index}\"); }}"
        f" setring r;"
    )
absolute_lines.append(f'print("ABSOLUTE_CENSUS={len(height_six)}");')
absolute_output = run_singular("\n".join(absolute_lines))
assert f"ABSOLUTE_CENSUS={len(height_six)}" in absolute_output


def genus_at(parameter: tuple[int, int, int]) -> tuple[tuple[int, int, int], int]:
    first, second, third = parameter
    source = singular_base(libraries=("normal.lib",)) + f"""
poly curve=({12 * first})*R0+({3 * second})*R1+({4 * third})*R3;
print("GENUS="+string(genus(ideal(curve))));
"""
    output = run_singular(source, timeout=20)
    markers = [
        line for line in output.splitlines() if line.startswith("GENUS=")
    ]
    assert len(markers) == 1
    return parameter, int(markers[0].split("=", 1)[1])


genus_results: dict[tuple[int, int, int], int] = {}
with ThreadPoolExecutor(max_workers=3) as executor:
    futures = [executor.submit(genus_at, parameter) for parameter in height_three]
    for future in as_completed(futures):
        parameter, genus_value = future.result()
        genus_results[parameter] = genus_value

assert len(genus_results) == len(height_three)
assert Counter(genus_results.values()) == Counter({14: 129, 12: 15})
for (_, _, third), genus_value in genus_results.items():
    assert genus_value == (12 if third == 0 else 14)


print("PASS: the sharp selector plane has an exact degree-ten strict model")
print("PASS: its two fixed infinity points have the complete rational tangent hierarchy")
print("PASS: the generic plane and C=0 pencil members are absolutely irreducible")
print("PASS: [9:-4:0] and [65:-4:4] are absolutely irreducible of genera 10 and 9")
print("PASS: the moving determinant is M^2 times one irreducible degree-23 curve")
print("PASS: [77:-16:-8] is absolutely irreducible of genus 12 with two new nodes")
print("PASS: [103:-16:8] splits into a rational line and an absolute genus-10 residual")
print("PASS: its norm splits into a rational conic and a degree-14 residual")
print("PASS: all 864 primitive parameters of height at most six are absolutely irreducible")
print("PASS: the height-three genus census is 15 of genus 12 and 129 of genus 14")
print("SCOPE: the height censuses are bounded and the full parameter discriminant remains open")
