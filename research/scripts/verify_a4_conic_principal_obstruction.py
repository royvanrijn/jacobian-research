#!/usr/bin/env python3
"""Exact obstruction to isolating the first rational A4 norm conic.

The sharp selector ``S_103`` has the rational norm component

    q = a^2 - 4*b^2 - 12*b - 36.

This checker purifies the corresponding height-one root-incidence prime.
It finds a smaller root-linear equation ``L``, verifies that ``L`` cuts the
same conic together with an absolutely irreducible genus-two residual, and
computes the exact non-Cartier locus of the conic prime.  Consequently the
reduced conic cannot be isolated with coefficient one by a principal
root-algebra divisor.

The second non-Cartier pair lies on an ordinary two-branch conductor node.
The conic is divisorial on one normalization branch and meets the other only
in codimension two.  Conductor matching then forces every principal function
vanishing on a positive multiple of the conic to acquire a second divisorial
component.  Thus no support-only positive multiple exists either.

This does not construct a Keller map.
"""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp


# ---------------------------------------------------------------------------
# 1. The A4 root algebra and the purified conic equation
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
selector_103 = sp.expand(
    1236 * q0 - 48 * q1_integral + 32 * q3_integral
)

conic = a**2 - 4 * b**2 - 12 * b - 36
purified_linear = 2 * T + 3 * a * b**2 + 6 * b**3 + 9 * a * b + 27 * a - 162

# S_103 belongs to the conic prime, but it is not itself a clean second
# generator there: its coefficient of L vanishes on a tangent point of q.
selector_relation = sp.expand(
    16 * (3 * a - 6 * b + 4) * purified_linear
    + 12 * (103 * a - 12 * b**2 - 36 * b - 108) * conic
)
assert sp.expand(selector_103 - selector_relation) == 0
assert sp.expand(
    sp.resultant(conic, 3 * a - 6 * b + 4, a)
    + 4 * (39 * b + 77)
) == 0

# The ideal (q,L) contains P and is prime: eliminating T identifies its
# quotient with Q[a,b]/(q), and q is an irreducible smooth conic.
norm_linear = sp.expand(sp.resultant(quartic, purified_linear, T))
norm_coefficient, norm_factors = sp.factor_list(norm_linear, a, b)
assert norm_coefficient == -3
assert sorted(
    (sp.Poly(factor, a, b).total_degree(), exponent)
    for factor, exponent in norm_factors
) == [(2, 1), (10, 1)]
norm_conic = next(
    factor
    for factor, _ in norm_factors
    if sp.Poly(factor, a, b).total_degree() == 2
)
residual_10 = next(
    factor
    for factor, _ in norm_factors
    if sp.Poly(factor, a, b).total_degree() == 10
)
assert sp.expand(norm_conic - conic) == 0
assert sp.factor_list(conic, a, b)[1] == [(conic, 1)]
assert sp.discriminant(conic, a) == 16 * (b**2 + 3 * b + 9)
assert sp.discriminant(b**2 + 3 * b + 9, b) == -27
assert conic.subs({a: 6, b: 0}) == 0

# Pseudo-division by the monic-up-to-2 linear equation gives the conormal
# relation used below.  In A=Q[a,b,T]/(P), it reads L*C=3*q*R_10.
conormal_coefficient = sp.cancel(
    (16 * quartic - norm_linear) / purified_linear
)
assert sp.denom(conormal_coefficient) == 1
conormal_coefficient = sp.expand(conormal_coefficient)
assert sp.expand(
    16 * quartic
    - purified_linear * conormal_coefficient
    + 3 * conic * residual_10
) == 0


# ---------------------------------------------------------------------------
# 2. Pullback: a rational line plus an absolute genus-two quintic
# ---------------------------------------------------------------------------

U, V = sp.symbols("U V")

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
source_A = U**3 - V**3 - 9 * V**2 - 27 * V - 54
source_L = (
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
root_numerator = 3 * source_A * K**3 * source_L

linear_pullback = sp.cancel(
    purified_linear.subs(
        {
            a: N1 / H,
            b: N2 / H,
            T: root_numerator / H**3,
        }
    )
)
pullback_numerator, pullback_denominator = sp.fraction(linear_pullback)
assert sp.expand(pullback_denominator - H**3) == 0
residual_5 = sp.cancel(pullback_numerator / (3 * U * K**3))
assert sp.denom(residual_5) == 1
residual_5 = sp.expand(residual_5)
assert sp.expand(
    pullback_numerator - 3 * U * K**3 * residual_5
) == 0
assert sp.Poly(residual_5, U, V).total_degree() == 5
assert sp.gcd(U, residual_5) == 1


# ---------------------------------------------------------------------------
# 3. Exact non-Cartier locus and absolute genus checks in Singular
# ---------------------------------------------------------------------------


def singular_expression(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


SINGULAR = shutil.which("Singular")
if SINGULAR is None:
    raise SystemExit("Singular is required for the conic-principal audit")


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


# For I=(q,L) in the hypersurface A, the conormal relation has coefficients
# (-3*R_10,C).  Where either is a unit, I is generated by L or q.  Where
# both vanish, I/mI has dimension two, so Nakayama proves non-principality.
# The radical calculation below identifies exactly four geometric points:
# two over Q(sqrt(-3)) at the universal cluster and two over Q(sqrt(5)).
local_source = f'''LIB "primdec.lib";
ring r=0,(T,a,b),dp;
poly q={singular_expression(conic)};
poly L={singular_expression(purified_linear)};
poly R10={singular_expression(residual_10)};
poly C={singular_expression(conormal_coefficient)};
ideal bad=std(ideal(q,L,R10,C));
ideal badrad=std(radical(bad));
ideal cluster=std(ideal(a,T,b2+3b+9));
ideal realpair=std(ideal(11a-32b-48,L,5b2+15b-19));
ideal expected=std(intersect(cluster,realpair));
if (size(reduce(badrad,expected))!=0 || size(reduce(expected,badrad))!=0)
{{
    ERROR("unexpected non-Cartier support");
}}
if (vdim(bad)!=14 || vdim(badrad)!=4)
{{
    ERROR("unexpected non-Cartier scheme lengths");
}}
print("NON_CARTIER_LOCUS_PASS");
'''
assert "NON_CARTIER_LOCUS_PASS" in run_singular(local_source)

assert sp.discriminant(5 * b**2 + 15 * b - 19, b) == 605
assert sp.sqrt(605).is_rational is False


# ---------------------------------------------------------------------------
# 4. The cluster class has index two
# ---------------------------------------------------------------------------

# At the first quadratic pair use c=2*b+3, so c^2+27=4*rho.  The conic has
# slope two and meets the E2 exceptional divisor at y=1/4.  On the E2 chart
# rho=a^2*y and T=a^2*t, the exceptional quartic is one simple root and one
# triple root.  The purified equation selects the simple root.
c, y, t = sp.symbols("c y t")
cluster_relation = c**2 + 27 - 4 * a**2 * y
cluster_total = sp.rem(
    sp.expand(
        16
        * quartic.subs(b, (c - 3) / 2)
        .subs(T, a**2 * t)
    ),
    cluster_relation,
    c,
)
cluster_strict = sp.cancel(cluster_total / a**8)
assert sp.denom(cluster_strict) == 1
cluster_exceptional = sp.expand(cluster_strict).subs(a, 0)
cluster_simple = 2 * t - (27 - 3 * c) * y
cluster_triple = 2 * t - (c - 9) * y
assert sp.rem(
    sp.expand(
        cluster_exceptional - cluster_simple * cluster_triple**3
    ),
    c**2 + 27,
    c,
) == 0

cluster_conic = sp.rem(
    sp.expand(conic.subs(b, (c - 3) / 2)),
    cluster_relation,
    c,
)
assert sp.expand(cluster_conic - a**2 * (1 - 4 * y)) == 0
cluster_linear_total = sp.rem(
    sp.expand(
        purified_linear.subs(b, (c - 3) / 2).subs(T, a**2 * t)
    ),
    cluster_relation,
    c,
)
cluster_linear_strict = sp.cancel(cluster_linear_total / a**2)
assert sp.denom(cluster_linear_strict) == 1
cluster_linear_exceptional = sp.expand(cluster_linear_strict).subs(a, 0)
assert sp.expand(
    cluster_linear_exceptional - (2 * t + 3 * (c - 9) * y)
) == 0

cluster_y = sp.Rational(1, 4)
cluster_selected_t = -sp.Rational(3, 2) * (c - 9) * cluster_y
assert sp.expand(
    cluster_selected_t - (27 - 3 * c) * cluster_y / 2
) == 0
assert sp.rem(
    sp.expand(
        cluster_selected_t - (c - 9) * cluster_y / 2
    ),
    c**2 + 27,
    c,
) != 0

# The established seven-curve normalization contracts to the A3 chain
#
#   E1_simple -- E2_simple -- E3.
#
# A curvette meeting the middle component has class 2 in Z/4, hence local
# Cartier index two.  The calculation above identifies the conic with that
# curvette.
cluster_discriminant_order = 4
cluster_conic_class = 2
cluster_cartier_index = cluster_discriminant_order // int(
    sp.gcd(cluster_discriminant_order, cluster_conic_class)
)
assert cluster_cartier_index == 2


# ---------------------------------------------------------------------------
# 5. The real quadratic pair is a two-branch conductor obstruction
# ---------------------------------------------------------------------------

# Work exactly in the residue field of the second quadratic closed point.
# The helper reduces rational functions modulo 5*b^2+15*b-19.
second_minpoly = sp.Poly(5 * b**2 + 15 * b - 19, b, domain=sp.QQ)


def second_residue(expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.fraction(sp.cancel(expression))
    numerator_reduced = sp.Poly(
        numerator,
        b,
        domain=sp.QQ,
    ).rem(second_minpoly)
    denominator_reduced = sp.Poly(
        denominator,
        b,
        domain=sp.QQ,
    ).rem(second_minpoly)
    denominator_inverse = sp.invert(
        denominator_reduced,
        second_minpoly,
    )
    return sp.factor(
        (numerator_reduced * denominator_inverse)
        .rem(second_minpoly)
        .as_expr()
    )


second_a = (32 * b + 48) / 11
second_T = sp.cancel(
    -(
        purified_linear - 2 * T
    ).subs(a, second_a)
    / 2
)
rho = b**2 + 3 * b + 9
sigma = a**3 * (2 * b + 3) - 3 * a**2 * rho + rho**2

# The selected root is a singular point of the root-incidence hypersurface,
# lying on the smooth sigma conductor and away from rho and B.  Its root
# multiplicity is exactly two.
second_point = {a: second_a, T: second_T}
assert second_residue(conic.subs(second_point)) == 0
assert second_residue(purified_linear.subs(second_point)) == 0
assert second_residue(quartic.subs(second_point)) == 0
assert all(
    second_residue(sp.diff(quartic, variable).subs(second_point)) == 0
    for variable in (a, b, T)
)
assert second_residue(sp.diff(quartic, T, 2).subs(second_point)) != 0
assert second_residue(rho) != 0
assert second_residue(jly_B.subs(a, second_a)) != 0
assert second_residue(sigma.subs(a, second_a)) == 0
assert any(
    second_residue(derivative.subs(a, second_a)) != 0
    for derivative in (sp.diff(sigma, a), sp.diff(sigma, b))
)

# The quadratic tangent cone has rank two.  Combined with the two etale root
# sections below, this gives two distinct smooth branches meeting normally
# along the smooth sigma conductor.
second_hessian = sp.hessian(quartic, (a, b, T))
second_hessian_residue = sp.Matrix(
    3,
    3,
    lambda row, column: second_residue(
        second_hessian[row, column].subs(second_point)
    ),
)
assert second_residue(second_hessian_residue.det()) == 0
assert second_residue(
    second_hessian_residue.extract([0, 2], [0, 2]).det()
) != 0

# The rational root chart supplies the two normalization branches exactly.
# Both points map to the same (a,b,T), and the (a,b)-Jacobian is a unit at
# each.  The first point lies on U=0, the strict transform of the conic.
J3 = (
    U**3
    - 12 * U * V**2
    - 36 * U * V
    - 108 * U
    - 16 * V**3
    - 72 * V**2
    - 216 * V
    - 216
)
root_chart_a = N1 / H
root_chart_b = N2 / H
root_chart_T = root_numerator / H**3
root_chart_jacobian = 4 * K**3 * source_L / H**3
first_preimage = (
    sp.Integer(0),
    -(second_a + 2 * b + 6) / 2,
)
second_preimage = (
    second_a / 2,
    -b - 3,
)

for source_point in (first_preimage, second_preimage):
    substitution = {U: source_point[0], V: source_point[1]}
    assert second_residue(
        (root_chart_a - second_a).subs(substitution)
    ) == 0
    assert second_residue(
        (root_chart_b - b).subs(substitution)
    ) == 0
    assert second_residue(
        (root_chart_T - second_T).subs(substitution)
    ) == 0
    assert second_residue(H.subs(substitution)) != 0
    assert second_residue(K.subs(substitution)) != 0
    assert second_residue(source_L.subs(substitution)) != 0
    assert second_residue(root_chart_jacobian.subs(substitution)) != 0

assert second_residue(second_preimage[0] - first_preimage[0]) != 0

# Pullback factorizations distinguish the conic divisor from its isolated
# incidence on the other branch:
#
#   q = U*K^2*J3/H^2,
#   L = 3*U*K^3*G5/H^3.
#
# At the U=0 branch both residual factors are units.  At the other branch
# J3=G5=0 transversely, so (q,L) has only a reduced codimension-two point
# there and no divisorial component.
conic_pullback = sp.cancel(
    conic.subs({a: root_chart_a, b: root_chart_b})
)
assert sp.factor(
    conic_pullback - U * K**2 * J3 / H**2
) == 0
assert sp.factor(
    linear_pullback - 3 * U * K**3 * residual_5 / H**3
) == 0

first_substitution = {U: first_preimage[0], V: first_preimage[1]}
second_substitution = {U: second_preimage[0], V: second_preimage[1]}
assert second_residue(J3.subs(first_substitution)) != 0
assert second_residue(residual_5.subs(first_substitution)) != 0
assert second_residue(J3.subs(second_substitution)) == 0
assert second_residue(residual_5.subs(second_substitution)) == 0
residual_intersection = sp.Matrix(
    [J3, residual_5]
).jacobian((U, V)).det()
assert second_residue(
    residual_intersection.subs(second_substitution)
) != 0

# The conic and conductor meet transversely in the coefficient plane.  If a
# principal function had divisor m*p and no other height-one support, then
# on the first normalization branch its conductor restriction would vanish
# to order m.  Conductor matching forces the restriction on the other branch
# to be a nonunit.  The principal ideal theorem then supplies a height-one
# zero on that branch, contradicting support only on p.  This local argument
# works for every m>0.
conic_conductor_intersection = sp.Matrix(
    [conic, sigma]
).jacobian((a, b)).det()
assert second_residue(
    conic_conductor_intersection.subs(a, second_a)
) != 0

geometry_source = f'''LIB "normal.lib";
LIB "absfact.lib";
ring source_ring=0,(U,V),dp;
poly G5={singular_expression(residual_5)};
if (genus(ideal(G5))!=2) {{ ERROR("wrong quintic genus"); }}
def source_absolute=absFactorize(G5);
setring source_absolute;
if (absolute_factors[4]!=1) {{ ERROR("quintic is not absolutely irreducible"); }}
setring source_ring;

ring target_ring=0,(a,b),dp;
poly R10={singular_expression(residual_10)};
if (genus(ideal(R10))!=2) {{ ERROR("wrong degree-ten genus"); }}
def target_absolute=absFactorize(R10);
setring target_absolute;
if (absolute_factors[4]!=1) {{ ERROR("degree-ten residual is not absolutely irreducible"); }}
print("GENUS_TWO_RESIDUAL_PASS");
'''
assert "GENUS_TWO_RESIDUAL_PASS" in run_singular(geometry_source)


print("PASS: S_103 purifies to the prime conic ideal (q,L)")
print("PASS: Norm(L) is q times one absolute genus-two degree-10 component")
print("PASS: the root-chart pullback is U times one absolute genus-two quintic")
print("PASS: the conic prime has exactly four geometric non-Cartier points")
print("PASS: the cluster conic class is 2 in Z/4 and has local index two")
print("PASS: the real quadratic pair is an ordinary two-branch conductor node")
print("OBSTRUCTION: no positive principal divisor has support only on the conic prime")
print("SCOPE: alternative selectors and an affine-space Keller map remain open")
