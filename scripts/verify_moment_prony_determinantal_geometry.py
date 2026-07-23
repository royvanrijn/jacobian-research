#!/usr/bin/env python3
"""Exact tests for intrinsic Moment--Prony determinantal geometry.

The equations under test are generated from optimal Gaussian moments, not by
substituting those moments into previously known seed equations.  Old
degree-six formulas and the universal coefficient atlases are used only as
independent scheme-theoretic comparisons after the intrinsic ideals exist.
"""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.moment_prony import (
    _universal_power_root,
    fixed_weight_hankel_prony,
    fixed_weight_332_prony,
    hessenberg_determinant_identity,
    moment_ritt_krylov_data,
    moment_fixed_weight_hankel,
    optimal_moment_seed,
    power_prony_data,
    primitive_numerator,
    reverse_power_sums,
    ritt_krylov_data,
)


def normalized_moment(H: sp.Expr, variable: sp.Symbol, order: int) -> sp.Expr:
    """Return ``mu_order=[W^(order-1)](1+H)^order/order``."""

    return sp.expand((1 + H) ** order).coeff(variable, order - 1) / order


# ---------------------------------------------------------------------------
# The determinant identity and the failure of ordinary unweighted Hankel rank.

for size in range(1, 5):
    assert hessenberg_determinant_identity(size) == 0

w = sp.Symbol("W")
nodes = (sp.Integer(1), sp.Integer(2), sp.Integer(3))


def weighted_power_sum(weights: tuple[int, ...], order: int) -> sp.Expr:
    return sum(weight * node**order for weight, node in zip(weights, nodes))


square_hankel = sp.Matrix(
    4, 4, lambda i, j: weighted_power_sum((2, 2, 2), i + j)
)
nonsquare_hankel = sp.Matrix(
    4, 4, lambda i, j: weighted_power_sum((4, 1, 1), i + j)
)
assert square_hankel.det() == 0
assert nonsquare_hankel.det() == 0
assert square_hankel[:3, :3].det() != 0
assert nonsquare_hankel[:3, :3].det() != 0
nonsquare_polynomial = (w - 1) ** 4 * (w - 2) * (w - 3)
assert sp.Poly(nonsquare_polynomial, w).sqf_part().degree() == 3
assert any(
    exponent % 2
    for _, exponent in sp.factor_list(nonsquare_polynomial, w)[1]
)


# ---------------------------------------------------------------------------
# Degree six: all-double and 3+3 from consecutive leading minors alone.

seed6 = optimal_moment_seed(6, variable=w)
u, v, omega = seed6.moments
c2, c3, c4, c5, c6 = (seed6.coefficient(j) for j in range(2, 7))
assert c2 == u
assert c3 == v
assert sp.expand(c4 - (omega - 2 * u**2)) == 0

double6 = power_prony_data(seed6, 2)
triple6 = power_prony_data(seed6, 3)
assert double6.equation_indices == (4,)
assert triple6.equation_indices == (3, 4)

F6 = primitive_numerator(double6.root_coefficients[3], seed6.moments)
G6 = primitive_numerator(triple6.root_coefficients[2], seed6.moments)
T6 = primitive_numerator(triple6.root_coefficients[3], seed6.moments)
assert sp.Poly(F6, *seed6.moments).total_degree() == 6
assert sp.Poly(G6, *seed6.moments).total_degree() == 6

# Independent comparison with the old coefficient equations.  These formulas
# do not participate in construction of F6,G6,T6.
old_all_double = (
    32 * c3 * c5 * c6**2
    + 64 * c3 * c6**3
    + 16 * c4**2 * c6**2
    - 24 * c4 * c5**2 * c6
    + 64 * c4 * c6**3
    + 5 * c5**4
    + 64 * c5 * c6**3
    + 64 * c6**4
)
old_ritt_32 = 27 * c3 * c6**2 - 18 * c4 * c5 * c6 + 5 * c5**3
assert sp.expand(F6 + old_all_double) == 0
assert sp.expand(G6 + old_ritt_32) == 0

ritt23_6 = moment_ritt_krylov_data(seed6, 2, 3)
ritt32_6 = moment_ritt_krylov_data(seed6, 3, 2)
K23 = primitive_numerator(ritt23_6.selected_minors[0], seed6.moments)
K32 = primitive_numerator(ritt32_6.selected_minors[0], seed6.moments)
assert sp.expand(K23 - F6) == 0
assert sp.expand(K32 + G6) == 0

# The all-triple normalization lies on both consecutive cubic-root minors.
q = sp.Symbol("q")
A = 15 * q**3 - 15 * q**2 + 6 * q - 1
B = 5 * q**2 - 5 * q + 1
triple_moments = {
    u: A / (3 * q * B),
    v: -(
        (3 * q**2 - 3 * q + 1)
        * (45 * q**4 - 15 * q**2 + 6 * q - 1)
    )
    / (27 * q**4 * B),
    omega: A
    * (30 * q**5 - 30 * q**4 - 3 * q**3 + 18 * q**2 - 8 * q + 1)
    / (9 * q**4 * B**2),
}
assert sp.factor(G6.subs(triple_moments)) == 0
assert sp.factor(T6.subs(triple_moments)) == 0

# Scheme structure is retained: on the reduced all-triple normalization the
# all-double determinant has a square, giving a dual-number transverse algebra
# at each of the four squarefree collision points.
collision6 = 45 * q**4 - 30 * q**3 + 15 * q**2 - 6 * q + 1
F6_pullback = sp.factor(F6.subs(triple_moments))
expected_F6_pullback = sp.factor(
    (3 * q - 1) ** 8 * collision6**2 / (3**11 * q**16 * B**4)
)
assert sp.factor(F6_pullback - expected_F6_pullback) == 0
assert sp.gcd(collision6, sp.diff(collision6, q)) == 1


# ---------------------------------------------------------------------------
# Degree eight in moments, with no elimination back to seed coefficients.

seed8 = optimal_moment_seed(8, variable=w)
double8 = power_prony_data(seed8, 2)
assert double8.equation_indices == (5, 6)
double8_equations = tuple(
    primitive_numerator(double8.root_coefficients[index - 1], seed8.moments)
    for index in double8.equation_indices
)
assert tuple(
    len(sp.Poly(equation, *seed8.moments).terms())
    for equation in double8_equations
) == (686, 1317)

ritt24_8 = moment_ritt_krylov_data(seed8, 2, 4)
ritt42_8 = moment_ritt_krylov_data(seed8, 4, 2)
assert ritt24_8.selected_rows == (2, 3)
assert ritt42_8.selected_rows == (3, 5)
moment_ritt24_equations = tuple(
    primitive_numerator(equation, seed8.moments)
    for equation in ritt24_8.selected_minors
)
moment_ritt42_equations = tuple(
    primitive_numerator(equation, seed8.moments)
    for equation in ritt42_8.selected_minors
)
assert all(moment_ritt24_equations)
assert all(moment_ritt42_equations)

# Compare the all-double and 2-o-4 ideals universally before applying any
# coordinate system.  Equality of generators proves equality of ideals, not
# merely equality of zero sets, and therefore retains nilpotents after every
# base change to moment coordinates.
c = sp.symbols("c2:8")
generic8 = w**8 + sum(c[j - 2] * w**j for j in range(2, 8))
generic_24 = ritt_krylov_data(generic8, w, 2, 4)
root_a, root_b = _universal_power_root(6, 2)
reverse_substitution = {root_a[j - 1]: c[6 - j] for j in range(1, 7)}
root5 = sp.together(root_b[4].subs(reverse_substitution)).as_numer_denom()[0]
root6 = sp.together(root_b[5].subs(reverse_substitution)).as_numer_denom()[0]
minor2 = sp.together(generic_24.selected_minors[0]).as_numer_denom()[0]
minor3 = sp.together(generic_24.selected_minors[1]).as_numer_denom()[0]
assert sp.expand(minor3 - root5) == 0
assert sp.expand(root6 - minor2 + 2 * c[-1] * root5) == 0

# The opposite factor order agrees generator-by-generator with the independent
# canonical degree-eight atlas, again as an ideal rather than only on points.
generic_42 = ritt_krylov_data(generic8, w, 4, 2)
c2g, c3g, c4g, c5g, c6g, c7g = c
atlas42 = (
    -2 * (256 * c3g - 128 * c4g * c7g + 20 * c6g * c7g**3 - 7 * c7g**5),
    -32 * c5g + 24 * c6g * c7g - 7 * c7g**3,
)
assert sp.expand(generic_42.selected_minors[0] + atlas42[0] / 512) == 0
assert sp.expand(generic_42.selected_minors[1] + atlas42[1] / 32) == 0

# A clean 2-o-2-o-2 witness satisfies all four degree-eight moment-minor
# systems by direct evaluation.
H8 = -w**2 * (w - 1) * (
    w**5 - 11 * w**4 + 37 * w**3 - 17 * w**2 - 189 * w + 519
) / 340
witness8 = {
    seed8.moments[order - 3]: normalized_moment(H8, w, order)
    for order in range(3, 8)
}
for equation in (
    *double8_equations,
    *moment_ritt24_equations,
    *moment_ritt42_equations,
):
    assert sp.factor(equation.subs(witness8)) == 0

# The degree-eight Ritt intersection is scheme-theoretically reduced.  In the
# 2-o-4 chart the two 4-o-2 fixed-pivot minors generate exactly (E), because
# one of them is a unit multiple of E; the other is (u2-u3^2)E.
u1, u2, u3, a1 = sp.symbols("u1 u2 u3 a1")
inner4 = w**4 + u3 * w**3 + u2 * w**2 + u1 * w
chart24 = sp.expand(inner4**2 + a1 * inner4)
pulled42 = ritt_krylov_data(chart24, w, 4, 2).selected_minors
E222 = 8 * u1 - 4 * u2 * u3 + u3**3
assert sp.factor(pulled42[0] - (u2 - u3**2) * E222 / 4) == 0
assert sp.factor(pulled42[1] - E222 / 4) == 0
assert sp.diff(E222, u1) == 8


# ---------------------------------------------------------------------------
# Mixed degree-eight weights: exact marked subresultants and the obstruction
# to claiming that their naive Fitting minors give the reduced global scheme.

mixed332 = fixed_weight_332_prony()
x = mixed332.distinguished_node
p = mixed332.power_sums
y, znode = sp.symbols("y znode")
weighted_parameterization = {
    p[k - 1]: 2 * x**k + 3 * y**k + 3 * znode**k for k in range(1, 7)
}
for relation in mixed332.relations:
    assert sp.expand(relation.subs(weighted_parameterization)) == 0
assert mixed332.fitting_matrix.shape == (3, 9)

# Substitution of Newton power sums places the marked 3+3+2 incidence directly
# over the five optimal degree-eight moment coordinates.
moment_power_sums = reverse_power_sums(seed8)
marked_moment_relations = tuple(
    relation.subs(dict(zip(p, moment_power_sums)))
    for relation in mixed332.relations
)
assert all(
    relation.free_symbols <= set(seed8.moments) | {x}
    for relation in marked_moment_relations
)
assert set(seed8.moments) <= set().union(
    *(relation.free_symbols for relation in marked_moment_relations)
)

# At the total collision p_k=8*t^k, every marked relation contains the cubic
# fiber (x-t)^3.  This is the precise nilpotent obstruction: maximal Fitting
# minors have the right support, but treating them as the reduced component
# ideal would silently add collision thickness.
t = sp.Symbol("t")
total_collision = {p[k - 1]: 8 * t**k for k in range(1, 7)}
collision_relations = tuple(
    sp.factor(relation.subs(total_collision)) for relation in mixed332.relations
)
assert sp.expand(collision_relations[0] + 80 * (x - t) ** 3) == 0
assert all(sp.rem(relation, (x - t) ** 3, x) == 0 for relation in collision_relations)

# The unmarked replacement uses the Hankel pencil, Christoffel numerator, and
# Sylvester rank drops.  It has one weight-two support root, two weight-three
# roots, and one flat-extension equation for p_6.
unmarked332 = fixed_weight_hankel_prony((2, 3, 3), p)
assert unmarked332.hankel_matrix.shape == (3, 3)
assert len(unmarked332.extension_equations) == 1
assert tuple(
    (
        condition.weight,
        condition.required_gcd_degree,
        condition.sylvester_matrix.shape,
        condition.rank_bound,
    )
    for condition in unmarked332.weight_conditions
) == (
    (2, 1, (7, 7), 6),
    (3, 2, (7, 7), 5),
)

vandermonde_squared = (x - y) ** 2 * (x - znode) ** 2 * (y - znode) ** 2
assert sp.expand(
    unmarked332.hankel_determinant.subs(weighted_parameterization)
    - 18 * vandermonde_squared
) == 0
assert sp.expand(
    unmarked332.support_polynomial.subs(weighted_parameterization)
    + 18
    * (x - unmarked332.variable)
    * (y - unmarked332.variable)
    * (znode - unmarked332.variable)
    * vandermonde_squared
) == 0
assert sp.factor(
    unmarked332.extension_equations[0].subs(weighted_parameterization)
) == 0

numeric_nodes = {x: 1, y: 2, znode: 4}
for condition in unmarked332.weight_conditions:
    specialized = condition.sylvester_matrix.subs(
        {**weighted_parameterization, **numeric_nodes}
    )
    assert specialized.rank() == condition.rank_bound

# The same matrices are evaluated directly over the five optimal moments.
unmarked_moment332 = moment_fixed_weight_hankel(seed8, (2, 3, 3))
allowed_symbols = set(seed8.moments) | {unmarked_moment332.variable}
assert unmarked_moment332.support_polynomial.free_symbols <= allowed_symbols
assert unmarked_moment332.hankel_determinant.free_symbols <= set(seed8.moments)
assert all(
    condition.sylvester_matrix.free_symbols <= allowed_symbols
    for condition in unmarked_moment332.weight_conditions
)

# Degree five is the shortest mixed case.  It has only p_1,p_2,p_3, so no
# flat-extension equation is available or needed: the pencil and two
# fixed-weight resultants already give the whole 3+2 model.
p5 = sp.symbols("q1:4")
mixed32 = fixed_weight_hankel_prony((2, 3), p5)
assert not mixed32.extension_equations
q_parameterization = {
    p5[k - 1]: 2 * x**k + 3 * y**k for k in range(1, 4)
}
for condition in mixed32.weight_conditions:
    specialized = condition.sylvester_matrix.subs(
        {**q_parameterization, **numeric_nodes}
    )
    assert specialized.rank() == condition.rank_bound

# A one-node collision stratum needs no weight resultant: p_0 already fixes
# its sole weight.  The remaining equations are exactly the support
# recurrence.
all_six_power_sums = sp.symbols("s1:5")
all_six_prony = fixed_weight_hankel_prony((6,), all_six_power_sums)
assert not all_six_prony.weight_conditions
assert len(all_six_prony.extension_equations) == 3
all_six_parameterization = {
    all_six_power_sums[k - 1]: 6 * x**k for k in range(1, 5)
}
assert all(
    sp.expand(equation.subs(all_six_parameterization)) == 0
    for equation in all_six_prony.extension_equations
)

# Finally compare the degree-eight mixed component with the all-double
# component on its normalization.  The two all-double equations acquire the
# common factor (y-z)^4.  Since the normalization coordinate transverse to
# the unordered triple-root collision is delta=(y-z)^2, this is delta^2 and
# gives the required length-two intersection, not a reduced divisor.
mixed_polynomial = sp.Poly(
    (w - x) ** 2 * (w - y) ** 3 * (w - znode) ** 3, w
)
mixed_coefficients = {
    j: mixed_polynomial.nth(j) for j in range(2, 9)
}
root_substitution = {
    root_a[j - 1]: mixed_coefficients[8 - j] for j in range(1, 7)
}
mixed_double_pullbacks = tuple(
    sp.factor(
        sp.together(root_b[index - 1].subs(root_substitution))
        .as_numer_denom()[0]
    )
    for index in (5, 6)
)
assert sp.expand(
    mixed_double_pullbacks[0]
    + 3 * (y - znode) ** 4 * (2 * x - y - znode)
) == 0
assert sp.expand(
    mixed_double_pullbacks[1]
    + (y - znode) ** 4
    * (
        12 * x * y
        + 12 * x * znode
        - 7 * y**2
        - 10 * y * znode
        - 7 * znode**2
    )
) == 0


print("PASS Moment-Prony: ordinary Hankel rank forgets fixed multiplicity weights")
print("PASS Moment-Prony: degree-six all-double is one log-Hessenberg leading minor")
print("PASS Moment-Prony: degree-six 3+3 is two consecutive leading minors")
print("PASS Moment-Prony: the degree-six collision pullback retains length two")
print("PASS Moment-Prony: degree-eight all-double has two intrinsic moment minors")
print("PASS Moment-Prony: both degree-eight Ritt ideals come from Krylov minors")
print("PASS Moment-Prony: the degree-eight Ritt intersection ideal is reduced")
print("PASS Moment-Prony: mixed 3+3+2 has a marked subresultant presentation")
print("PASS Moment-Prony: naive mixed-weight Fitting minors thicken total collision")
print("PASS Moment-Prony: saturated Christoffel subresultants give the unmarked 3+3+2 closure")
print("PASS Moment-Prony: degree five 3+2 is covered at the minimal Hankel length")
print("PASS Moment-Prony: one-node collision strata reduce to the support recurrence")
print("PASS Moment-Prony: the degree-eight mixed/all-double intersection retains length two")
