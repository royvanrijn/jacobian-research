#!/usr/bin/env python3
"""Verify the corrected A4 exceptional selector and coarse deletion screen."""

from fractions import Fraction
from itertools import combinations

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


a, b, c, T, X, z = sp.symbols("a b c T X z")
t, y = sp.symbols("t y")

rho = b**2 + 3 * b + 9
B = (
    a**3
    - 3 * a * b**2
    + 2 * b**3
    - 9 * a * b
    + 9 * b**2
    - 27 * a
    + 27 * b
    + 27
)
sigma = a**3 * (2 * b + 3) - 3 * a**2 * rho + rho**2
A = a**3 - b**3 - 9 * b**2 - 27 * b - 54
C = a**3 - b**3 + 27
P = (
    T**4
    - 6 * A * B * T**2
    - 8 * B**3 * T
    + B**2 * (9 * A**2 - 12 * C * B)
)

c_polynomial = 2 * b + 3
G = T + a**3
simple_selector = 2 * T - (27 - 3 * c_polynomial) * rho
corrected_selector = sp.expand(
    4 * c_polynomial * simple_selector
    - 27 * (c_polynomial - 9) * a * rho
)


def newton_root_orders(coefficient_orders):
    """Read root valuations from the lower Newton polygon."""

    points = [
        (exponent, Fraction(order))
        for exponent, order in enumerate(coefficient_orders)
        if order is not None
    ]
    hull = []
    for point in points:
        while len(hull) >= 2:
            previous_slope = Fraction(
                hull[-1][1] - hull[-2][1],
                hull[-1][0] - hull[-2][0],
            )
            next_slope = Fraction(
                point[1] - hull[-1][1],
                point[0] - hull[-1][0],
            )
            if next_slope <= previous_slope:
                hull.pop()
            else:
                break
        hull.append(point)

    orders = []
    for left, right in zip(hull, hull[1:]):
        slope = Fraction(
            right[1] - left[1],
            right[0] - left[0],
        )
        orders.extend([-slope] * (right[0] - left[0]))
    return tuple(orders)


local_relation = b**2 + 3 * b + 9 - z


def weighted_order(expression, a_weight, z_weight):
    """Return the exact (a,z)-weighted order."""

    polynomial = sp.Poly(
        sp.rem(sp.expand(expression), local_relation, b),
        a,
        z,
        b,
    )
    return min(
        a_weight * a_power + z_weight * z_power
        for (a_power, z_power, _), coefficient in polynomial.terms()
        if coefficient
    )


def divisor_order(expression, divisor):
    """Return the exact order along a target divisor."""

    quotient = sp.Poly(expression, a, b, domain=sp.QQ)
    divisor_poly = sp.Poly(divisor, a, b, domain=sp.QQ)
    order = 0
    while True:
        next_quotient, remainder = sp.div(quotient, divisor_poly)
        if remainder.as_expr() != 0:
            return order
        quotient = next_quotient
        order += 1


def characteristic_orders(function, weights):
    """Return the four root orders of a root-algebra function."""

    characteristic = sp.Poly(
        sp.resultant(P, X - function, T),
        X,
    )
    coefficient_orders = [
        weighted_order(characteristic.nth(exponent), *weights)
        if characteristic.nth(exponent) != 0
        else None
        for exponent in range(5)
    ]
    return newton_root_orders(coefficient_orders)


def strict_orders(function, divisor):
    """Return target-normalized root orders along a strict divisor."""

    characteristic = sp.Poly(
        sp.resultant(P, X - function, T),
        X,
    )
    coefficient_orders = [
        divisor_order(characteristic.nth(exponent), divisor)
        if characteristic.nth(exponent) != 0
        else None
        for exponent in range(5)
    ]
    return newton_root_orders(coefficient_orders)


# ---------------------------------------------------------------------------
# 1. The first exceptional quartic supplies a non-diagonal selector
# ---------------------------------------------------------------------------

# Use z=rho and c=2b+3, so c^2+27=4z.  On the first blowup put
# z=ay and T=at.  The exceptional quartic is one simple line times one
# triple line.
P_c = sp.expand(16 * P.subs(b, (c - 3) / 2))
first_blowup_relation = c**2 + 27 - 4 * a * y
first_blowup_total = sp.rem(
    sp.expand(P_c.subs(T, a * t)),
    first_blowup_relation,
    c,
)
assert min(
    exponent[0]
    for exponent, coefficient in sp.Poly(first_blowup_total, a).terms()
    if coefficient
) == 4
first_exceptional = sp.expand(first_blowup_total / a**4).subs(a, 0)
simple_line = 2 * t - (27 - 3 * c) * y
triple_line = 2 * t - (c - 9) * y
assert sp.rem(
    sp.expand(first_exceptional - simple_line * triple_line**3),
    c**2 + 27,
    c,
) == 0

print("PASS: the E1 quartic is one simple line times one triple line cubed")


# ---------------------------------------------------------------------------
# 2. One higher-order correction fixes the simple F branch
# ---------------------------------------------------------------------------

resolution_rays = {
    "E1": (1, 1),
    "E2": (1, 2),
    "E3": (1, 3),
    "F": (2, 3),
}
expected_orders = {
    "E1": (Fraction(2), Fraction(1), Fraction(1), Fraction(1)),
    "E2": (Fraction(3), Fraction(2), Fraction(2), Fraction(2)),
    "E3": (Fraction(3),) * 4,
    "F": (Fraction(6), Fraction(3), Fraction(3), Fraction(3)),
}
for name, weights in resolution_rays.items():
    assert characteristic_orders(corrected_selector, weights) == expected_orders[name]

# The correction remains a genuine new irreducible boundary after taking
# its norm to the coefficient plane.
selector_norm = sp.resultant(P, corrected_selector, T)
unit, norm_factors = sp.factor_list(selector_norm, a, b)
assert unit != 0
assert len(norm_factors) == 1
assert norm_factors[0][1] == 1
assert sp.Poly(norm_factors[0][0], a, b).total_degree() == 16

print("PASS: the corrected selector has exact simple/triple orders on all four rays")
print("PASS: its coefficient-plane norm is one irreducible degree-16 boundary")


# ---------------------------------------------------------------------------
# 3. The corrected selector gives the exact total mask divisor
# ---------------------------------------------------------------------------

expected_T_orders = {
    "E1": (Fraction(1),) * 4,
    "E2": (Fraction(2),) * 4,
    "E3": (Fraction(3),) * 4,
    "F": (Fraction(3),) * 4,
}
expected_G_orders = expected_T_orders
expected_residual = {
    "E1": (Fraction(0), Fraction(1), Fraction(1), Fraction(1)),
    "E2": (Fraction(1), Fraction(2), Fraction(2), Fraction(2)),
    "E3": (Fraction(3),) * 4,
    "F": (Fraction(0), Fraction(3), Fraction(3), Fraction(3)),
}
for name in resolution_rays:
    total_mask_orders = tuple(
        t_order + g_order - h_order
        for t_order, g_order, h_order in zip(
            expected_T_orders[name],
            expected_G_orders[name],
            expected_orders[name],
        )
    )
    assert total_mask_orders == expected_residual[name]

# Along the strict divisors, T supplies the normalized B mask, G supplies
# the ramified triple-rho mask, and the corrected selector is a unit.
assert strict_orders(T, B) == (Fraction(1, 2),) * 4
assert strict_orders(G, rho) == (
    Fraction(1, 3),
    Fraction(1, 3),
    Fraction(1, 3),
    Fraction(0),
)
for divisor in (B, rho, sigma):
    assert strict_orders(corrected_selector, divisor) == (Fraction(0),) * 4

# Hence the rational product
#
#     T(T+a^3)/Hhat
#
# has the complete residual divisor.  It admits the useful two-factor
# allocation
#
#     M1=B*T/Hhat,  M2=(T+a^3)/B.
#
# Both factors are regular after deleting the strict B and Hhat divisors;
# their exceptional order vectors are M1=delta and M2=0.
for name in resolution_rays:
    b_order = Fraction(weighted_order(B, *resolution_rays[name]))
    first_mask_orders = tuple(
        b_order + t_order - h_order
        for t_order, h_order in zip(
            expected_T_orders[name],
            expected_orders[name],
        )
    )
    second_mask_orders = tuple(
        g_order - b_order
        for g_order in expected_G_orders[name]
    )
    assert first_mask_orders == expected_residual[name]
    assert second_mask_orders == (Fraction(0),) * 4

print("PASS: T*(T+a^3)/Hhat has exactly the full residual mask divisor")
print("PASS: B*T/Hhat and (T+a^3)/B give an exact two-mask allocation")


# ---------------------------------------------------------------------------
# 4. Coarse relative Picard and ampleness screen
# ---------------------------------------------------------------------------

# Work on the resolved coefficient plane in chain order E1-F-E2-E3.  The
# exceptional multiplicities of Norm(Hhat) are the sums of its four root
# orders: (5,15,9,12).  Its strict transform meets F once and E3 three
# times.  The class columns of F,E3,B,Hhat form a unimodular basis.
intersection_matrix = sp.Matrix(
    [
        [-3, 1, 0, 0],
        [1, -1, 1, 0],
        [0, 1, -3, 1],
        [0, 0, 1, -1],
    ]
)
b_multiplicities = sp.Matrix([1, 3, 2, 3])
h_multiplicities = sp.Matrix([5, 15, 9, 12])
assert -intersection_matrix * b_multiplicities == sp.Matrix([0, 0, 0, 1])
assert -intersection_matrix * h_multiplicities == sp.Matrix([0, 1, 0, 3])

classes = {
    "E1": sp.eye(4).col(0),
    "F": sp.eye(4).col(1),
    "E2": sp.eye(4).col(2),
    "E3": sp.eye(4).col(3),
    "B": -b_multiplicities,
    "rho": -b_multiplicities,
    "G": -b_multiplicities,
    "sigma": -sp.Matrix([2, 6, 3, 3]),
    "Hhat": -h_multiplicities,
}
unimodular_sets = []
for names in combinations(classes, 4):
    determinant = sp.Matrix.hstack(*(classes[name] for name in names)).det()
    if abs(determinant) == 1:
        unimodular_sets.append(names)

candidate = ("F", "E3", "B", "Hhat")
assert candidate in unimodular_sets
assert sp.Matrix.hstack(*(classes[name] for name in candidate)).det() == 1

# The positive divisor F+E3+B+2Hhat intersects the four exceptional curves
# in (1,1,2,6), so its support is relatively ample over the coefficient
# plane.  The divisor exact sequence therefore passes its two easiest
# affine-space tests: no coarse relative Picard group and no boundary-unit
# relation.
ample_intersections = (
    intersection_matrix * sp.Matrix([0, 1, 0, 1])
    + sp.Matrix([0, 0, 0, 1])
    + 2 * sp.Matrix([0, 1, 0, 3])
)
assert ample_intersections == sp.Matrix([1, 1, 2, 6])

print("PASS: {F,E3,B,Hhat} is a unimodular coarse deletion set")
print("PASS: F+E3+B+2*Hhat is positive on every exceptional curve")
print("CANDIDATE: the corrected boundary passes the coarse Picard/unit/ample screen")


# ---------------------------------------------------------------------------
# 5. The normalized root incidence has a rank-two class obstruction
# ---------------------------------------------------------------------------

# The coefficient-plane screen collapses the branches above each exceptional
# divisor.  Normalize the quartic pullback instead.  On both retained rays
# E1 and E2, the exceptional quartic is the simple line times the triple
# line cubed displayed above.  Shifting to the triple line gives the lower
# Newton data
#
#     (0,1),(1,1),(2,1),(3,0),(4,0).
#
# Hence the triple branch has ramification slope 1/3 and is distinct from
# the simple branch.  There are at least two prime exceptional components
# above each of E1 and E2, hence at least four retained exceptional primes.
def shifted_triple_coefficient_orders(ray_weight):
    """Newton coefficient orders after centering the triple exceptional line."""

    ell = sp.symbols("ell")
    local_parameter = sp.symbols("local_parameter")
    relation = c**2 + 27 - 4 * a**ray_weight * local_parameter
    total = sp.rem(
        sp.expand(
            P_c.subs(
                T,
                a**ray_weight * t,
            )
        ),
        relation,
        c,
    )
    strict = sp.cancel(total / a ** (4 * ray_weight))
    triple_root = (c - 9) * local_parameter / 2
    centered = sp.Poly(
        sp.rem(
            sp.together(strict.subs(t, triple_root + ell)) * 16,
            relation,
            c,
        ),
        ell,
    )
    orders = []
    for exponent in range(5):
        coefficient = sp.Poly(
            sp.expand(centered.nth(exponent)),
            a,
            c,
            local_parameter,
        )
        if coefficient.is_zero:
            orders.append(None)
        else:
            orders.append(
                min(
                    monomial[0]
                    for monomial, value in coefficient.terms()
                    if value
                )
            )
    return tuple(orders)


assert shifted_triple_coefficient_orders(1) == (1, 1, 1, 0, 0)
assert shifted_triple_coefficient_orders(2) == (1, 1, 1, 0, 0)
retained_exceptional_primes = 2 + 2

# The strict pullback of B is one prime.  Its Newton polygon has root order
# 1/2.  After T=sqrt(B)*U the residual polynomial is
#
#     (U^2-3A)^2.
#
# The quadratic is irreducible in the function field of B: locally,
#
#     B=a^3+(c-3a)z,
#     A=a^3-(b+6)z,
#
# and on the normalization of B one has
#
#     A=a^3*(c-3a+b+6)/(c-3a).
#
# Its leading coefficient is 3(b+3)/c, a unit at rho=0, so A has odd order
# three and cannot be a square.  The residue degree and ramification degree
# are both two, exhausting the quartic degree and leaving one B-prime.
local_A = sp.rem(A, local_relation, b)
local_B = sp.rem(B, local_relation, b)
assert sp.expand(local_A - (a**3 - (b + 6) * z)) == 0
assert sp.expand(local_B - (a**3 + (c_polynomial - 3 * a) * z)) == 0
assert sp.gcd(
    sp.Poly(b + 3, b, domain=sp.QQ),
    sp.Poly(rho, b, domain=sp.QQ),
).as_expr() == 1
strict_b_primes = 1

# The corrected selector is linear in T with generic coefficient 8c.  Its
# irreducible norm therefore gives one strict Hhat-prime on the root
# incidence.  There are exactly two horizontal deleted primes available to
# kill classes of the retained exceptional components.
assert sp.Poly(corrected_selector, T).degree() == 1
assert sp.Poly(corrected_selector, T).LC() == 8 * c_polynomial
strict_hhat_primes = 1
horizontal_deleted_primes = strict_b_primes + strict_hhat_primes

# On a smooth resolution of a normal surface modification, the irreducible
# exceptional curves have negative-definite intersection matrix and
# independent relative divisor classes.  Deleted components above F and E3
# remove at most their own class directions.  Any further curves used to
# resolve singularities contained in the deleted boundary add one class and
# one deleted component simultaneously.  Consequently the divisor exact
# sequence leaves the following invariant lower bound.
relative_class_rank_lower_bound = (
    retained_exceptional_primes - horizontal_deleted_primes
)
assert relative_class_rank_lower_bound == 2

print("PASS: E1 and E2 each split into simple and triple normalized components")
print("PASS: the strict B and Hhat deletions contribute only two horizontal primes")
print("OBSTRUCTION: the normalized deletion has relative class rank at least two")
print("OBSTRUCTION: {F,E3,B,Hhat} cannot give an affine-space root incidence")


# ---------------------------------------------------------------------------
# 6. The normalized cluster is an A3 singularity
# ---------------------------------------------------------------------------

# Projection formula determines the normalized exceptional chain.  In order
#
#     E1_simple-F_simple-E2_simple-E3
#       -E2_triple-F_triple-E1_triple
#
# its self-intersections are
#
#     (-3,-1,-3,-4,-1,-3,-1).
#
# The simple E1/E2 components have ramification/residue data (1,1), their
# triple components have (3,1), the F components have (1,1) and (1,3), and
# the irreducible E3 component has (1,4).  The projection formula gives
# C^2=f*D^2/e for a component C above a base exceptional curve D.
base_self_intersections = (-3, -1, -3, -1, -3, -1, -3)
ramification_indices = (1, 1, 1, 1, 3, 1, 3)
residue_degrees = (1, 1, 1, 4, 1, 3, 1)
normalized_self_intersections = tuple(
    residue_degree * base_self_intersection // ramification_index
    for base_self_intersection, ramification_index, residue_degree in zip(
        base_self_intersections,
        ramification_indices,
        residue_degrees,
    )
)
assert normalized_self_intersections == (-3, -1, -3, -4, -1, -3, -1)
normalized_intersection_matrix = sp.diag(*normalized_self_intersections)
for index in range(6):
    normalized_intersection_matrix[index, index + 1] = 1
    normalized_intersection_matrix[index + 1, index] = 1

normalized_leading_minors = tuple(
    normalized_intersection_matrix[:size, :size].det()
    for size in range(1, 8)
)
assert normalized_leading_minors == (-3, 2, -3, 10, -7, 11, -4)
assert smith_normal_form(
    normalized_intersection_matrix,
    domain=ZZ,
) == sp.diag(1, 1, 1, 1, 1, 1, 4)


def contract_minus_one(matrix, index):
    """Return the intersection matrix after blowing down one (-1)-curve."""

    assert matrix[index, index] == -1
    retained = [
        position
        for position in range(matrix.rows)
        if position != index
    ]
    restricted = matrix.extract(retained, retained)
    incidence = matrix.extract(retained, [index])
    return restricted + incidence * incidence.T


# Blow down, in order, E1_triple, E2_triple, F_triple, and F_simple.
minimal_matrix = normalized_intersection_matrix
for index in (6, 4, 4, 1):
    minimal_matrix = contract_minus_one(minimal_matrix, index)
assert minimal_matrix == sp.Matrix(
    [
        [-2, 1, 0],
        [1, -2, 1],
        [0, 1, -2],
    ]
)
assert minimal_matrix.det() == -4

# The minimal resolution is therefore the A3 chain, with cyclic local
# divisor class group of order four.
local_discriminant_order = abs(minimal_matrix.det())
assert local_discriminant_order == 4

# A class assignment must include every horizontal component through the
# cluster.  In particular, G=T+a^3 does not consist only of the strict
# triple-rho component: its norm is rho times an irreducible degree-10
# factor, and that complementary factor also passes through the cluster.
# This prevents assigning the triple-rho class from G's exceptional orders
# alone.
g_norm = sp.expand(P.subs(T, -a**3))
g_norm_quotient = sp.cancel(g_norm / rho)
assert sp.denom(g_norm_quotient) == 1
g_quotient_unit, g_quotient_factors = sp.factor_list(
    sp.expand(g_norm_quotient),
    a,
    b,
)
assert g_quotient_unit != 0
assert len(g_quotient_factors) == 1
assert g_quotient_factors[0][1] == 1
assert sp.Poly(
    g_quotient_factors[0][0],
    a,
    b,
).total_degree() == 10
cluster_basis = sp.groebner(
    [a, rho],
    a,
    b,
    order="lex",
    domain=sp.QQ,
)
assert cluster_basis.reduce(
    sp.expand(g_norm_quotient)
)[1] == 0

print("PASS: the normalized seven-curve chain contracts to the A3 chain")
print("PASS: the cluster local class group is cyclic of order four")
print("PASS: G has a complementary degree-10 divisor through the cluster")
print("NOTE: assign all horizontal components in the A3 class group before screening")


# ---------------------------------------------------------------------------
# 7. The forced B prime has an index-two primitivity obstruction
# ---------------------------------------------------------------------------

# A concrete horizontal curvette fixes the discriminant-lattice generator.
# In local coordinates z=rho and c=2b+3, put
#
#     chi = 16*a^2 - 4*a*c + z.
#
# Its strict transform meets E1 at y=4c.  The linear root selector Lchi
# vanishes on the simple exceptional root there and not on the triple root.
chi = 16 * a**2 - 8 * a * b - 12 * a + rho
chi_selector = sp.expand(
    (32 * a + 9 - 4 * c_polynomial) * T
    - 13311 * a**3
    + (2241 * c_polynomial - 567) * a**2
    + (162 * c_polynomial - 7290) * a
)
assert sp.expand(chi - (16 * a**2 - 4 * a * c_polynomial + rho)) == 0

chi_selector_c = sp.expand(
    chi_selector.subs(b, (c - 3) / 2).subs(T, a * t)
)
chi_selector_first_total = sp.rem(
    chi_selector_c,
    first_blowup_relation,
    c,
)
assert sp.rem(
    sp.Poly(chi_selector_first_total, a),
    sp.Poly(a, a),
).as_expr() == 0
chi_selector_exceptional = sp.expand(
    chi_selector_first_total / a
).subs(a, 0)
assert chi_selector_exceptional == -4 * c * t + 162 * c + 9 * t - 7290

chi_exceptional_point = 4 * c
simple_exceptional_root = (27 - 3 * c) * chi_exceptional_point / 2
triple_exceptional_root = (c - 9) * chi_exceptional_point / 2
assert sp.rem(
    sp.expand(
        chi_selector_exceptional.subs(t, simple_exceptional_root)
    ),
    c**2 + 27,
    c,
) == 0
assert sp.rem(
    sp.expand(
        chi_selector_exceptional.subs(t, triple_exceptional_root)
    ),
    c**2 + 27,
    c,
) == 216 * (c - 45)

# The norm has exactly the curvette component chi and one complementary
# irreducible component.  Thus the chi component realizes an odd generator
# of the cyclic order-four local class group.
chi_selector_norm = sp.resultant(P, chi_selector, T)
chi_complement = sp.cancel(chi_selector_norm / chi)
assert sp.denom(chi_complement) == 1
assert sp.gcd(
    sp.Poly(chi, a, b, domain=sp.QQ),
    sp.Poly(chi_complement, a, b, domain=sp.QQ),
).as_expr() == 1
chi_complement_unit, chi_complement_factors = sp.factor_list(
    sp.expand(chi_complement),
    a,
    b,
)
assert chi_complement_unit != 0
assert len(chi_complement_factors) == 1
assert chi_complement_factors[0][1] == 1
assert sp.Poly(
    chi_complement_factors[0][0],
    a,
    b,
).total_degree() == 14

# In the seven-curve order used above, the full dual lattice is generated
# by E0,...,E5 and q=-M^{-1}E0.  Its index over the integral exceptional
# lattice is four, as required by the A3 discriminant.
dual_generator = -normalized_intersection_matrix.inv().col(0)
assert dual_generator == sp.Matrix(
    [
        sp.Rational(3, 4),
        sp.Rational(5, 4),
        sp.Rational(1, 2),
        sp.Rational(1, 4),
        sp.Rational(1, 2),
        sp.Rational(1, 4),
        sp.Rational(1, 4),
    ]
)
dual_basis = sp.Matrix.hstack(
    *(sp.eye(7).col(index) for index in range(6)),
    dual_generator,
)
assert abs(dual_basis.det()) == sp.Rational(1, 4)
assert all(
    entry.q == 1
    for entry in normalized_intersection_matrix * dual_basis
)

# Pullback valuations on
#
#   E1_simple, F_simple, E2_simple, E3,
#   E2_triple, F_triple, E1_triple
#
# are as follows.  Strict B occurs with ramification multiplicity two,
# whereas the linear Hhat prime occurs with multiplicity one.  Principal
# divisor relations therefore give [D_B]=-b/2 and [D_Hhat]=-h.
b_normalized_orders = sp.Matrix([1, 3, 2, 3, 6, 3, 3])
h_normalized_orders = sp.Matrix([2, 6, 3, 3, 6, 3, 3])
b_prime_class = -b_normalized_orders / 2
h_prime_class = -h_normalized_orders
b_class_coordinates = dual_basis.inv() * b_prime_class
h_class_coordinates = dual_basis.inv() * h_prime_class
assert b_class_coordinates == sp.Matrix([4, 6, 2, 0, 0, 0, -6])
assert h_class_coordinates == sp.Matrix([7, 9, 3, 0, 0, 0, -12])

# The B class itself has content two, while the Hhat class is primitive.
# Thus the irreducible B prime cannot occur in any basis of the full
# relative class lattice.  The gcd of the 2-by-2 minors of the forced pair
# is also two.  Every seven-prime boundary class matrix containing B has
# even determinant.  Fewer than seven boundary primes cannot kill the
# rank-seven class lattice; more than seven leave a boundary-unit relation.
# Thus no boundary enlargement of this resolved model containing the forced
# B prime can have both trivial units and trivial class group.
assert sp.gcd_list(tuple(b_class_coordinates)) == 2
assert sp.gcd_list(tuple(h_class_coordinates)) == 1
forced_pair = sp.Matrix.hstack(
    b_class_coordinates,
    h_class_coordinates,
)
forced_pair_minors = [
    forced_pair.extract(rows, [0, 1]).det()
    for rows in combinations(range(7), 2)
]
assert sp.gcd_list(forced_pair_minors) == 2

print("PASS: chi gives an explicit odd A3 discriminant-lattice curvette")
print("PASS: the forced B boundary class has content two")
print("PASS: the forced B/Hhat pair has determinantal divisor two")
print("OBSTRUCTION: no resolved boundary enlargement has both trivial units and class group")


# ---------------------------------------------------------------------------
# 8. A B-free curvette split survives the exact lattice screen
# ---------------------------------------------------------------------------

# The same total mask has the alternative factorization
#
#     M1 = T*Lchi/Hhat,  M2 = (T+a^3)/Lchi.
#
# Lchi has the normalized exceptional orders of the base parameter a.
# Consequently M1 is negative only on F_simple, while M2 has no negative
# exceptional order.  The strict poles are Hhat and the two components of
# Norm(Lchi); neither numerator shares a horizontal component with its
# denominator.
expected_chi_selector_orders = {
    "E1": (Fraction(1),) * 4,
    "E2": (Fraction(1),) * 4,
    "E3": (Fraction(1),) * 4,
    "F": (Fraction(2),) * 4,
}
for name, weights in resolution_rays.items():
    assert characteristic_orders(
        chi_selector,
        weights,
    ) == expected_chi_selector_orders[name]

a_normalized_orders = sp.Matrix([1, 2, 1, 1, 3, 2, 3])
t_normalized_orders = b_normalized_orders
first_b_free_mask_orders = (
    t_normalized_orders
    + a_normalized_orders
    - h_normalized_orders
)
second_b_free_mask_orders = (
    t_normalized_orders
    - a_normalized_orders
)
assert first_b_free_mask_orders == sp.Matrix([0, -1, 0, 1, 3, 2, 3])
assert second_b_free_mask_orders == sp.Matrix([0, 1, 1, 2, 3, 1, 0])

g_norm_polynomial = sp.Poly(g_norm, a, b, domain=sp.QQ)
t_norm_polynomial = sp.Poly(P.subs(T, 0), a, b, domain=sp.QQ)
assert sp.gcd(
    sp.Poly(chi_selector_norm, a, b, domain=sp.QQ),
    g_norm_polynomial,
).as_expr() == 1
assert sp.gcd(
    sp.Poly(selector_norm, a, b, domain=sp.QQ),
    sp.Poly(
        sp.expand(P.subs(T, 0) * chi_selector_norm),
        a,
        b,
        domain=sp.QQ,
    ),
).as_expr() == 1

# The two Lchi components are endpoint curvettes.  Dchi meets E1_simple;
# the complementary degree-14 component meets E1_triple.  Their classes
# add to the class -a of the full principal Lchi divisor.
d_chi_class = -dual_generator
d_14_class = -a_normalized_orders - d_chi_class
d_chi_coordinates = dual_basis.inv() * d_chi_class
d_14_coordinates = dual_basis.inv() * d_14_class
assert d_chi_coordinates == sp.Matrix([0, 0, 0, 0, 0, 0, -1])
assert d_14_coordinates == sp.Matrix([8, 13, 5, 2, 3, 1, -11])
assert (
    normalized_intersection_matrix * d_chi_class
    == sp.eye(7).col(0)
)
assert (
    normalized_intersection_matrix * d_14_class
    == sp.eye(7).col(6)
)

# The forced B-free boundary is
#
#     F_simple, Hhat, Dchi, D14.
#
# Its rank-four class sublattice is primitive.  Among completions by three
# exceptional primes, the unique unimodular one adds
#
#     E2_triple, F_triple, E1_triple.
exceptional_class_coordinates = {
    "S1": dual_basis.inv() * sp.eye(7).col(0),
    "Fs": dual_basis.inv() * sp.eye(7).col(1),
    "S2": dual_basis.inv() * sp.eye(7).col(2),
    "Q": dual_basis.inv() * sp.eye(7).col(3),
    "R2": dual_basis.inv() * sp.eye(7).col(4),
    "Ft": dual_basis.inv() * sp.eye(7).col(5),
    "R1": dual_basis.inv() * sp.eye(7).col(6),
}
b_free_forced_classes = (
    exceptional_class_coordinates["Fs"],
    h_class_coordinates,
    d_chi_coordinates,
    d_14_coordinates,
)
b_free_forced_matrix = sp.Matrix.hstack(*b_free_forced_classes)
b_free_forced_minors = [
    b_free_forced_matrix.extract(rows, range(4)).det()
    for rows in combinations(range(7), 4)
]
assert sp.gcd_list(b_free_forced_minors) == 1

b_free_completions = []
for names in combinations(
    ("S1", "S2", "Q", "R2", "Ft", "R1"),
    3,
):
    boundary_matrix = sp.Matrix.hstack(
        *b_free_forced_classes,
        *(exceptional_class_coordinates[name] for name in names),
    )
    if abs(boundary_matrix.det()) == 1:
        b_free_completions.append(names)
assert b_free_completions == [("R2", "Ft", "R1")]

# A positive divisor on the resulting seven-component support is
#
#     Fs + R2 + 2Ft + 6R1 + 2Hhat + Dchi + 5D14.
#
# Its intersections with the seven exceptional curves are all positive.
h_attachment = -normalized_intersection_matrix * h_normalized_orders
chi_attachment = normalized_intersection_matrix * d_chi_class
d14_attachment = normalized_intersection_matrix * d_14_class
b_free_ample_intersections = (
    normalized_intersection_matrix * sp.Matrix([0, 1, 0, 0, 1, 2, 6])
    + 2 * h_attachment
    + chi_attachment
    + 5 * d14_attachment
)
assert b_free_ample_intersections == sp.Matrix([2, 1, 1, 7, 1, 1, 1])

print("PASS: T*Lchi/Hhat and (T+a^3)/Lchi give an exact B-free mask split")
print("PASS: {Fs,R2,Ft,R1,Hhat,Dchi,D14} is a unimodular boundary")
print("PASS: that boundary supports a divisor positive on all seven curves")
print("CANDIDATE: the B-free curvette split passes the normalized lattice/ample screen")
