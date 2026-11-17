#!/usr/bin/env python3
"""Exact rank theorem and the degree-eight component-poset test."""

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.discriminant_geometry import (  # noqa: E402
    affine_difference_mason_defect,
    collision_partition_blocks,
    collision_precedes,
    contact_partition_incidence,
    exceptional_stratum_codimension,
    exceptional_stratum_dimension,
    maximal_two_three_partitions,
    multiple_omission_incidence,
    normalized_seed_space_dimension,
    two_omission_incidence,
    weighted_newton_top_coefficients,
    weighted_vandermonde_determinant,
)


def restricted_partitions(total, maximum=None):
    """Partitions of total with every part at least two."""
    if total == 0:
        yield ()
        return
    maximum = min(total, maximum or total)
    for first in range(maximum, 1, -1):
        for tail in restricted_partitions(total - first, first):
            yield (first,) + tail


# The all-degree coefficient-rank certificate.  Newton identities make the
# passage from weighted power sums to top monic coefficients triangular, and
# the remaining determinant is a weighted Vandermonde.
for length in range(1, 6):
    roots = sp.symbols(f"r0:{length}")
    multiplicities = tuple(range(2, length + 2))
    coefficients = weighted_newton_top_coefficients(multiplicities, roots)
    jacobian = sp.det(sp.Matrix(coefficients).jacobian(roots))
    assert sp.factor(
        jacobian - weighted_vandermonde_determinant(multiplicities, roots)
    ) == 0

# Nonemptiness starts on the maximally collided locus M=(W-r)^n.  Its Phi
# polynomial has n-2 distinct roots; D is nonzero at all of them, and the
# weighted forbidden factor A does not vanish at all of them.  The written
# proof derives this from (z^n-nz+n-1)/(z-1)^2 and then splits the root by the
# formal implicit-function theorem.  This loop is a bounded exact regression
# of those three polynomial assertions, not the all-degree proof.
W_nonempty, r_nonempty = sp.symbols("W_nonempty r_nonempty")
for degree in range(3, 41):
    collided = (W_nonempty - r_nonempty) ** degree
    phi = sp.Poly(
        collided.subs(W_nonempty, 1)
        - collided.subs(W_nonempty, 0)
        - sp.diff(collided, W_nonempty).subs(W_nonempty, 0),
        r_nonempty,
    )
    D_collided = sp.Poly(
        sp.diff(collided, W_nonempty).subs(W_nonempty, 1)
        - sp.diff(collided, W_nonempty).subs(W_nonempty, 0),
        r_nonempty,
    )
    A_collided = sp.Poly(
        sp.diff(collided, W_nonempty, 2).subs(W_nonempty, 1)
        - 2 * D_collided.as_expr(),
        r_nonempty,
    )
    assert phi.degree() == degree - 2
    assert sp.gcd(phi, phi.diff()).degree() == 0
    assert sp.gcd(phi, D_collided).degree() == 0
    assert sp.rem(A_collided, phi).as_expr() != 0

# The theorem's dimension and codimension formulas are uniform over all
# permitted partitions, rather than a table of low-degree rank witnesses.
for degree in range(3, 21):
    for partition in restricted_partitions(degree):
        assert exceptional_stratum_dimension(partition) == len(partition) - 1
        assert exceptional_stratum_codimension(degree, partition) == (
            degree - len(partition) - 2
        )
        assert exceptional_stratum_codimension(degree, partition) == (
            normalized_seed_space_dimension(degree)
            - exceptional_stratum_dimension(partition)
        )

# Collision is a partial order: lambda <= mu means that mu is obtained by
# merging blocks of parts of lambda.  Substituting one root for every block
# sends M_lambda and Phi_lambda exactly to M_mu and Phi_mu.  Tangent-chord
# normalization supplies the nearby off-diagonal deformation used in the
# closure proof.
for degree in range(3, 11):
    partitions = tuple(restricted_partitions(degree))
    assert all(collision_precedes(partition, partition) for partition in partitions)
    assert all(
        any(
            collision_precedes(maximal, partition)
            for maximal in maximal_two_three_partitions(degree)
        )
        for partition in partitions
    )
    for finer in partitions:
        for coarser in partitions:
            forward = collision_precedes(finer, coarser)
            backward = collision_precedes(coarser, finer)
            if forward and backward:
                assert finer == coarser
            if not forward:
                continue
            blocks = collision_partition_blocks(finer, coarser)
            W_collision = sp.symbols(f"W_collision_{degree}")
            fine_roots = sp.symbols(f"fine_{degree}_0:{len(finer)}")
            coarse_roots = sp.symbols(f"coarse_{degree}_0:{len(coarser)}")
            substitution = {
                fine_roots[index]: coarse_roots[block_index]
                for block_index, block in enumerate(blocks)
                for index in block
            }
            fine_polynomial = sp.prod(
                (W_collision - root) ** multiplicity
                for root, multiplicity in zip(fine_roots, finer)
            )
            coarse_polynomial = sp.prod(
                (W_collision - root) ** multiplicity
                for root, multiplicity in zip(coarse_roots, coarser)
            )
            assert sp.expand(fine_polynomial.subs(substitution) - coarse_polynomial) == 0
            fine_phi = (
                fine_polynomial.subs(W_collision, 1)
                - fine_polynomial.subs(W_collision, 0)
                - sp.diff(fine_polynomial, W_collision).subs(W_collision, 0)
            )
            coarse_phi = (
                coarse_polynomial.subs(W_collision, 1)
                - coarse_polynomial.subs(W_collision, 0)
                - sp.diff(coarse_polynomial, W_collision).subs(W_collision, 0)
            )
            assert sp.expand(fine_phi.subs(substitution) - coarse_phi) == 0
    for first in partitions:
        for second in partitions:
            for third in partitions:
                if collision_precedes(first, second) and collision_precedes(
                    second, third
                ):
                    assert collision_precedes(first, third)

# Parts two and three give the minimal elements in the collision order and
# therefore the maximal stratum closures: neither can be split into permitted
# parts, whereas every m>=4 can.  For two distinct maximal stratum types,
# Mason--Stothers gives
# n <= ell(lambda)+ell(mu), but the combinatorics gives a strict deficit.
component_counts = []
for degree in range(3, 51):
    maximal = maximal_two_three_partitions(degree)
    permitted = {
        partition
        for partition in restricted_partitions(degree)
        if set(partition) <= {2, 3}
    }
    assert set(maximal) == permitted
    component_counts.append(len(maximal))
    cutoff = degree // 3
    explicit_count = (
        cutoff // 2 + 1 if degree % 2 == 0 else (cutoff + 1) // 2
    )
    assert len(maximal) == explicit_count
    assert min(
        exceptional_stratum_codimension(degree, partition)
        for partition in maximal
    ) == (degree + 1) // 2 - 2
    if degree <= 14:
        all_partitions = tuple(restricted_partitions(degree))
        for maximal_partition in maximal:
            coarsenings = tuple(
                partition
                for partition in all_partitions
                if collision_precedes(maximal_partition, partition)
            )
            for outside_partition in all_partitions:
                if collision_precedes(maximal_partition, outside_partition):
                    continue
                assert all(
                    affine_difference_mason_defect(
                        outside_partition, coarsening
                    )
                    >= 1
                    for coarsening in coarsenings
                )
    for index, left_partition in enumerate(maximal):
        for right_partition in maximal[index + 1 :]:
            assert affine_difference_mason_defect(
                left_partition, right_partition
            ) >= 1
            if degree <= 14:
                for coarsening in restricted_partitions(degree):
                    if collision_precedes(right_partition, coarsening):
                        assert affine_difference_mason_defect(
                            left_partition, coarsening
                        ) >= 1

z = sp.symbols("z")
component_series = sp.series(1 / ((1 - z**2) * (1 - z**3)), z, 0, 51).removeO()
assert component_counts == [
    sp.expand(component_series).coeff(z, degree) for degree in range(3, 51)
]

# The uniform two-omission API exposes the equation requested by the theorem.
two_omission_six = two_omission_incidence(6, (2, 2, 2), (3, 3))
assert len(two_omission_six.high_coefficient_equations) == 5
assert sp.factor(two_omission_six.high_coefficient_equations[-1]) == (
    two_omission_six.left_scale - two_omission_six.right_scale
)
assert two_omission_six.collision_discriminant != 0

# Degree eight is a theorem test.  Its only partitions using parts 2 and 3 are
# the predicted maximal types, with the dimensions forced by the theorem.
degree_eight_partitions = tuple(restricted_partitions(8))
assert set(degree_eight_partitions) == {
    (8,),
    (6, 2),
    (5, 3),
    (4, 4),
    (4, 2, 2),
    (3, 3, 2),
    (2, 2, 2, 2),
}
assert {
    partition: (
        exceptional_stratum_dimension(partition),
        exceptional_stratum_codimension(8, partition),
    )
    for partition in ((2, 2, 2, 2), (3, 3, 2))
} == {(2, 2, 2, 2): (3, 2), (3, 3, 2): (2, 3)}

square = contact_partition_incidence(8, (2, 2, 2, 2))
triple_triple_double = contact_partition_incidence(8, (3, 3, 2))
six_two = contact_partition_incidence(8, (6, 2))
eight = contact_partition_incidence(8, (8,))

# Both predicted maximal closures contain E_(6,2).  For Q^2 this is
# Q=(W-u)^3(W-v); for A^3 B^2 it is the collision of the two roots of A.
u, v = six_two.quotient_coordinates
e1, e2, e3, e4 = square.quotient_coordinates
p1, p2, b = triple_triple_double.quotient_coordinates
six_two_to_square = {
    e1: 3 * u + v,
    e2: 3 * u**2 + 3 * u * v,
    e3: u**3 + 3 * u**2 * v,
    e4: u**3 * v,
}
six_two_to_triple = {p1: 2 * u, p2: u**2, b: v}
assert sp.expand(square.Phi.subs(six_two_to_square) - six_two.Phi) == 0
assert sp.expand(
    triple_triple_double.Phi.subs(six_two_to_triple) - six_two.Phi
) == 0

# Their common collision boundary itself specializes to E_(8).
r = eight.quotient_coordinates[0]
assert sp.expand(six_two.Phi.subs({u: r, v: r}) - eight.Phi) == 0

# Test for genuinely different omitted values.  If
# Q(W)^2-A(W)^3(W-b)^2 is affine-linear, comparison from W^7 down to W^4
# determines Q.  The W^3 and W^2 coefficients then have the factorizations
# below.  Away from disc(A)=0 they force a1=2b and a2=b^2, a contradiction.
W = sp.symbols("W")
a1, a2 = sp.symbols("a1 a2")
Q = W**4 - e1 * W**3 + e2 * W**2 - e3 * W + e4
A = W**2 - a1 * W + a2
B = W - b
left = sp.Poly(Q**2, W)
right = sp.Poly(A**3 * B**2, W)
coefficient_differences = {
    power: sp.expand(
        left.coeff_monomial(W**power) - right.coeff_monomial(W**power)
    )
    for power in range(2, 8)
}
solved = {}
for power, coefficient in ((7, e1), (6, e2), (5, e3), (4, e4)):
    solved[coefficient] = sp.factor(
        sp.solve(coefficient_differences[power].subs(solved), coefficient)[0]
    )
disc_A = a1**2 - 4 * a2
degree_three_obstruction = sp.factor(coefficient_differences[3].subs(solved))
degree_two_obstruction = sp.factor(coefficient_differences[2].subs(solved))
assert sp.factor(
    degree_three_obstruction
    + sp.Rational(3, 128) * (a1 - 2 * b) * disc_A**2
) == 0
assert sp.factor(
    degree_two_obstruction
    - sp.Rational(1, 512)
    * disc_A**2
    * (11 * a1**2 - 12 * a1 * b + 4 * a2 - 24 * b**2)
) == 0
assert sp.factor(
    (degree_two_obstruction / disc_A**2).subs(a1, 2 * b)
    - (a2 - b**2) / 128
) == 0
assert sp.factor(disc_A.subs({a1: 2 * b, a2: b**2})) == 0

# Constructing the general shared incidence also audits that this compact
# calculation is the degree-eight instance of the public API.
degree_eight_multiple = multiple_omission_incidence(
    degree=8, partitions=((2, 2, 2, 2), (3, 3, 2))
)
assert degree_eight_multiple.partitions == ((2, 2, 2, 2), (3, 3, 2))
assert len(degree_eight_multiple.common_value_ideal) == 1

print("PASS: weighted Newton coefficients give the all-degree rank theorem")
print("PASS: maximally collided Phi has admissible roots in the nonemptiness regression")
print("PASS: every full-contact stratum has dimension ell(lambda)-1")
print("PASS: collision merging is a partial order and gives every closure inclusion")
print("PASS: Mason excludes all off-collision pairs of maximal 2/3 partitions")
print("PASS: component counts and total codimension have the uniform closed forms")
print("PASS: degree eight has maximal candidates (2,2,2,2) and (3,3,2)")
print("PASS: their closures share (6,2), with deeper boundary (8)")
print("PASS: their exact degree-eight intersection has no off-diagonal points")
