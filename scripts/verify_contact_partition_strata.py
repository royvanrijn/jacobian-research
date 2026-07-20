#!/usr/bin/env python3
"""Exact milestones for arbitrary full-contact partition strata."""

import itertools
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.discriminant_geometry import contact_partition_incidence  # noqa: E402


def primitive_coefficients(data):
    polynomial = sp.Poly(data.primitive, data.variable)
    return tuple(
        sp.cancel(polynomial.coeff_monomial(data.variable**power))
        for power in range(3, data.degree)
    )


def irreducible(polynomial):
    factors = sp.factor_list(polynomial)[1]
    return len(factors) == 1 and factors[0][1] == 1


# The arbitrary API also covers partial contact, with a monic residual factor
# and resultants preventing it from meeting any marked contact point.
partial = contact_partition_incidence(7, (3, 2), full_contact=False)
assert sp.degree(partial.residual_polynomial, partial.variable) == 2
assert len(partial.residual_coefficients) == 2
assert partial.expected_seed_dimension == 3
assert partial.residual_nonvanishing_factor != 1

# Equal-part permutations are quotiented by elementary symmetric coordinates.
triple_double = contact_partition_incidence(6, (2, 2, 2))
assert triple_double.symmetry_order == 6
assert len(triple_double.quotient_coordinates) == 3
raw_substitution = dict(triple_double.root_to_quotient)
raw_product = sp.prod(
    (triple_double.variable - root) ** 2
    for root in triple_double.marked_roots
)
assert sp.expand(triple_double.M.subs(raw_substitution) - raw_product) == 0

# Every coefficient-space ideal is returned as an exact saturated elimination
# presentation, and its matching equations vanish under the rational map.
for degree, partitions in {
    5: ((3, 2),),
    6: ((2, 2, 2), (3, 3), (4, 2), (6,)),
    7: ((3, 2, 2), (4, 3), (5, 2), (7,)),
}.items():
    for partition in partitions:
        data = contact_partition_incidence(degree, partition)
        ideal = data.coefficient_space_elimination_ideal
        assert ideal.contact_partition == tuple(sorted(partition, reverse=True))
        assert ideal.parameters == data.coefficient_parameters
        assert ideal.marked_points == data.quotient_coordinates
        coefficient_map = dict(zip(data.coefficient_parameters, primitive_coefficients(data)))
        W = data.variable
        assert sp.cancel(data.H.subs(W, 0)) == 0
        assert sp.cancel(sp.diff(data.H, W).subs(W, 0)) == 0
        assert sp.cancel(sp.diff(data.H, W).subs(W, 1) + 1) == 0
        assert sp.cancel(data.scale_denominator * data.H.subs(W, 1) + data.Phi) == 0
        assert sp.cancel(data.H - data.s * W + data.t - data.a * data.M) == 0
        assert all(
            sp.factor(equation.subs(coefficient_map)) == 0
            for equation in ideal.equations[1:]
        )

# Degree five: derive Phi_(3,2), map to the old (R,P) seed coordinates, and
# recover the established exceptional polynomial exactly up to 5/64.
degree_five = contact_partition_incidence(5, (3, 2))
W = degree_five.variable
a, b = degree_five.quotient_coordinates
R, P = sp.symbols("R P")
monic = sp.Poly(degree_five.M, W)
root_to_R = sp.expand(-monic.coeff_monomial(W**4) - 1)
root_to_P = sp.expand(monic.coeff_monomial(W**3) - root_to_R)
assert root_to_R == 3 * a + 2 * b - 1
b_from_R = sp.solve(sp.Eq(R, root_to_R), b)[0]
phi_in_a = sp.factor(degree_five.phi.subs(b, b_from_R))
P_equation = sp.factor((P - root_to_P).subs(b, b_from_R))
degree_five_elimination = sp.factor(sp.resultant(phi_in_a, P_equation, a))
F = (
    20 * P**3
    + 219 * P**2 * R**2
    - 312 * P**2 * R
    + 84 * P**2
    - 120 * P * R**4
    + 174 * P * R**3
    + 54 * P * R**2
    - 204 * P * R
    + 96 * P
    + 16 * R**6
    - 24 * R**5
    - 21 * R**4
    + 58 * R**3
    - 21 * R**2
    - 24 * R
    + 16
)
assert sp.factor(degree_five_elimination / F - sp.Rational(5, 64)) == 0

# Degree six: derive the main (2,2,2) equation by rationally recovering the
# cubic Q from the top coefficients of M=Q^2.
main = triple_double
W = main.variable
x, y, z = main.coefficient_parameters
h6 = sp.Poly(main.normalized_universal_primitive, W).coeff_monomial(W**6)
e1, e2, e3 = main.quotient_coordinates
m3, m4, m5 = x / h6, y / h6, z / h6
inverse_e1 = -m5 / 2
inverse_e2 = (m4 - inverse_e1**2) / 2
inverse_e3 = -(m3 + 2 * inverse_e1 * inverse_e2) / 2
main_pullback = sp.cancel(
    main.phi.subs({e1: inverse_e1, e2: inverse_e2, e3: inverse_e3})
)
F6 = sp.factor(main_pullback.as_numer_denom()[0])
assert sp.Poly(F6, x, y, z).total_degree() == 4
assert irreducible(F6)
assert sp.degree(main.phi, e3) == 1
assert sp.factor(sp.diff(main.phi, e3) - 2 * (e1 - 1)) == 0

degree_six = {
    partition: contact_partition_incidence(6, partition)
    for partition in ((2, 2, 2), (3, 3), (4, 2), (6,))
}
assert {
    partition: data.expected_seed_dimension for partition, data in degree_six.items()
} == {(2, 2, 2): 2, (3, 3): 1, (4, 2): 1, (6,): 0}
for data in degree_six.values():
    assert irreducible(data.phi)
    assert sp.gcd(
        sp.Poly(data.phi, *data.quotient_coordinates),
        sp.Poly(data.weighted_admissibility_factor, *data.quotient_coordinates),
    ).total_degree() == 0

# The (4,2) curve is the double-root boundary of Q in M=Q^2.
four_two = degree_six[(4, 2)]
u, v = four_two.quotient_coordinates
four_two_to_main = {e1: 2 * u + v, e2: u**2 + 2 * u * v, e3: u**2 * v}
assert sp.expand(main.phi.subs(four_two_to_main) - four_two.phi) == 0

# It is not contained in the singular locus: one exact admissible point is
# smooth on the quartic hypersurface.
four_two_point = {u: -4, v: 1}
four_two_coefficients = [value.subs(four_two_point) for value in primitive_coefficients(four_two)]
four_two_seed = dict(zip((x, y, z), four_two_coefficients))
assert F6.subs(four_two_seed) == 0
assert any(sp.diff(F6, variable).subs(four_two_seed) != 0 for variable in (x, y, z))

# The generic (3,3) curve is not in the main closure.  This rational point is
# admissible, satisfies Phi_(3,3), and does not satisfy F6.
three_three = degree_six[(3, 3)]
q1, q2 = three_three.quotient_coordinates
three_three_point = {q1: sp.Rational(-1, 2), q2: sp.Rational(-3, 2)}
assert three_three.phi.subs(three_three_point) == 0
assert three_three.weighted_admissibility_factor.subs(three_three_point) != 0
three_three_seed = dict(
    zip((x, y, z), (value.subs(three_three_point) for value in primitive_coefficients(three_three)))
)
assert F6.subs(three_three_seed) == -sp.Rational(640000, 177147)

# Every (6) point is the common collision boundary of (4,2), (3,3), and the
# main square stratum.  All are admissible and smooth on the main quartic.
six = degree_six[(6,)]
r = six.quotient_coordinates[0]
assert sp.expand(four_two.phi.subs({u: r, v: r}) - six.phi) == 0
assert sp.expand(three_three.phi.subs({q1: 2 * r, q2: r**2}) - six.phi) == 0
assert sp.expand(main.phi.subs({e1: 3 * r, e2: 3 * r**2, e3: r**3}) - six.phi) == 0
six_seed = dict(zip((x, y, z), primitive_coefficients(six)))
first_gradient_on_six = sp.cancel(sp.diff(F6, x).subs(six_seed)).as_numer_denom()[0]
assert sp.gcd(sp.Poly(six.phi, r), sp.Poly(first_gradient_on_six, r)).degree() == 0

# Degree seven: the leading stratum is an irreducible twofold in the
# four-dimensional normalized seed space.  An exact rank-two point proves the
# image really has dimension two, hence codimension two rather than one.
degree_seven = {
    partition: contact_partition_incidence(7, partition)
    for partition in ((3, 2, 2), (4, 3), (5, 2), (7,))
}
assert {
    partition: data.expected_seed_dimension for partition, data in degree_seven.items()
} == {(3, 2, 2): 2, (4, 3): 1, (5, 2): 1, (7,): 0}
for data in degree_seven.values():
    assert irreducible(data.phi)
    assert sp.gcd(
        sp.Poly(data.phi, *data.quotient_coordinates),
        sp.Poly(data.weighted_admissibility_factor, *data.quotient_coordinates),
    ).total_degree() == 0

leading_seven = degree_seven[(3, 2, 2)]
coordinates = leading_seven.quotient_coordinates
rank_point = dict(zip(coordinates, (-6, -3, -4)))
assert leading_seven.phi.subs(rank_point) == 0
assert leading_seven.weighted_admissibility_factor.subs(rank_point) != 0
seed_functions = primitive_coefficients(leading_seven)
jacobian_rows = [
    [sp.cancel(sp.diff(function, coordinate)).subs(rank_point) for coordinate in coordinates]
    for function in (leading_seven.phi, *seed_functions[:2])
]
assert sp.Matrix(jacobian_rows).det() != 0
normalized_seed_dimension = 7 - 3
assert normalized_seed_dimension - leading_seven.expected_seed_dimension == 2

print("PASS: arbitrary partitions, residual factors, and equal-part quotients")
print("PASS: degree-five (3,2) elimination recovers 5*F(R,P)/64")
print("PASS: degree-six main locus is an irreducible rational quartic surface")
print("PASS: degree-six closures and (6) collision intersections are exact")
print("PASS: degree-seven nonsurjectivity has genuine codimension two")
