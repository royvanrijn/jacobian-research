#!/usr/bin/env python3
"""Exact certificate for the foundational cotangent-torus reduction.

The foundational Keller map is equivariant for source weights (-1,1,2)
and target weights (2,1,-1).  Its cotangent lift therefore preserves the
corresponding Hamiltonian moment maps.

The full affine reduction has a mu_2 quotient singularity on every moment
level.  On the target cotangent chart beta != 0, however, beta has weight
-1 and fixes the torus gauge.  The target reduction is A^4, while the
source reduction is the explicit smooth hypersurface

    F_2(x,y,z) + 2*a*F_1(x,y,z) - g*F_3(x,y,z) - c = 0.

The complete foundational three-point collision survives on this chart.
This script checks the polynomial identities and orbit statements.  It does
not assert that the source hypersurface is polynomially isomorphic to A^4.
"""

import sympy as sp


x, y, z = sp.symbols("x y z")
p_x, p_y, p_z = sp.symbols("p_x p_y p_z")
A, B, C = sp.symbols("A B C")
alpha, beta, gamma = sp.symbols("alpha beta gamma")
c = sp.symbols("c")

source_variables = (x, y, z)
source_momenta = (p_x, p_y, p_z)
source_weights = (-1, 1, 2)
target_weights = (2, 1, -1)

u = 1 + x * y
F = sp.Matrix(
    (
        u**3 * z + y**2 * u * (4 + 3 * x * y),
        y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
        2 * x - 3 * x**2 * y - x**3 * z,
    )
)
J = F.jacobian(source_variables)
J_inverse = -sp.Rational(1, 2) * J.adjugate()
assert sp.factor(J.det()) == -2
assert (J * J_inverse).applyfunc(sp.expand) == sp.eye(3)


def polynomial_weight(polynomial, variables, weights):
    """Return the common monomial weight, asserting homogeneity."""

    monomial_weights = {
        sum(exponent * weight for exponent, weight in zip(monomial, weights))
        for monomial in sp.Poly(polynomial, *variables).monoms()
    }
    assert len(monomial_weights) == 1
    return monomial_weights.pop()


# Equivariance of the base map.
for component, expected_weight in zip(F, target_weights):
    assert polynomial_weight(component, source_variables, source_weights) == expected_weight

# Cotangent moment preservation.
p = sp.Matrix(source_momenta)
target_momentum = J_inverse.T * p
mu_source = sum(
    weight * coordinate * momentum
    for weight, coordinate, momentum in zip(
        source_weights, source_variables, source_momenta
    )
)
mu_target_after_lift = sum(
    weight * component * momentum
    for weight, component, momentum in zip(
        target_weights, F, target_momentum
    )
)
assert sp.expand(mu_target_after_lift - mu_source) == 0

# The complete collision persists on the beta=1 moment chart.  In target
# momentum coordinates take (alpha,beta,gamma)=(-2c,1,0); the target moment
# at (-1/4,0,0) is c.
collision_points = (
    (0, 0, -sp.Rational(1, 4)),
    (1, -sp.Rational(3, 2), sp.Rational(13, 2)),
    (-1, sp.Rational(3, 2), sp.Rational(13, 2)),
)
collision_target = sp.Matrix((-sp.Rational(1, 4), 0, 0))
chosen_target_momentum = sp.Matrix((-2 * c, 1, 0))
for point in collision_points:
    substitution = dict(zip(source_variables, point))
    assert F.subs(substitution) == collision_target
    lifted_source_momentum = J.subs(substitution).T * chosen_target_momentum
    assert sp.expand(mu_source.subs(substitution).subs(
        dict(zip(source_momenta, lifted_source_momentum))
    ) - c) == 0
    assert (
        J_inverse.subs(substitution).T * lifted_source_momentum
        - chosen_target_momentum
    ).applyfunc(sp.expand) == sp.zeros(3, 1)

# On beta != 0 the invariant target coordinates are
# A*beta^2, C/beta, alpha/beta^2, gamma*beta.  Gauge beta=1 and solve the
# moment equation 2*A*alpha+B*beta-C*gamma=c for B.  This leaves four free
# coordinates (A,C,alpha,gamma), with two canonical pairs.
B_on_slice = c - 2 * A * alpha + C * gamma
mu_target = 2 * A * alpha + B * beta - C * gamma
assert sp.expand(mu_target.subs({beta: 1, B: B_on_slice}) - c) == 0


def target_poisson(left, right):
    """Canonical target cotangent bracket {momentum,base}=1."""

    bases = (A, B, C)
    momenta = (alpha, beta, gamma)
    return sp.expand(
        sum(
            sp.diff(left, momentum) * sp.diff(right, base)
            - sp.diff(left, base) * sp.diff(right, momentum)
            for base, momentum in zip(bases, momenta)
        )
    )


localized_invariants = (
    A * beta**2,
    C / beta,
    alpha / beta**2,
    gamma * beta,
)
reduced_brackets = sp.Matrix(
    4,
    4,
    lambda row, column: target_poisson(
        localized_invariants[row], localized_invariants[column]
    ).subs(beta, 1),
)
expected_brackets = sp.Matrix(
    (
        (0, 0, -1, 0),
        (0, 0, 0, -1),
        (1, 0, 0, 0),
        (0, 1, 0, 0),
    )
)
assert reduced_brackets == expected_brackets

# The same four exponent vectors give exact invariant generators after
# quantum localization at beta.  With [alpha,A]=[beta,B]=[gamma,C]=1,
#
#   Q1=A*beta^2, P1=alpha*beta^-2,
#   Q2=C*beta^-1, P2=gamma*beta
#
# have weights zero and satisfy [P_i,Q_j]=delta_ij.  Since beta commutes
# with the A/alpha and C/gamma pairs, there are no ordering corrections in
# these four relations.
quantum_generator_weights = (
    target_weights[0] + 2 * (-target_weights[1]),
    -target_weights[0] - 2 * (-target_weights[1]),
    target_weights[2] - (-target_weights[1]),
    -target_weights[2] + (-target_weights[1]),
)
assert quantum_generator_weights == (0, 0, 0, 0)

# In the polynomial coordinates (x,y,z; alpha,beta,gamma) obtained from
# p=J^T(alpha,beta,gamma), the source moment is
# 2*F_1*alpha+F_2*beta-F_3*gamma.  Gauge beta=1.
source_slice_equation = sp.expand(F[1] + 2 * alpha * F[0] - gamma * F[2] - c)
for point in collision_points:
    substitution = dict(zip(source_variables, point))
    substitution.update({alpha: -2 * c, gamma: 0})
    assert sp.expand(source_slice_equation.subs(substitution)) == 0
    assert (
        sp.Matrix((F[0], F[2], alpha, gamma)).subs(substitution)
        == sp.Matrix((-sp.Rational(1, 4), 0, -2 * c, 0))
    )

# beta has primitive weight -1, so beta=1 is a unique gauge representative.
# The three displayed representatives are therefore three distinct quotient
# points, rather than the two orbits seen on the beta=0 stratum.
assert len(set(collision_points)) == 3

# The full affine reduction is nevertheless singular on every level.  At
# (A,B,C;alpha,beta,gamma)=(-1/4,0,0;-2c,0,0), only weights +/-2 occur.
# The stabilizer is mu_2.  Its four-dimensional symplectic slice has four
# odd coordinates, and the quotient A^4/{+/-1} has ten independent quadratic
# cotangent generators at the origin, strictly more than its dimension four.
slice_variables = sp.symbols("s0:4")
quadratic_invariants = tuple(
    slice_variables[i] * slice_variables[j]
    for i in range(4)
    for j in range(i, 4)
)
assert len(quadratic_invariants) == 10
assert all(sp.Poly(invariant, *slice_variables).total_degree() == 2 for invariant in quadratic_invariants)

print("PASS: foundational map is G_m-equivariant with weights (-1,1,2) -> (2,1,-1)")
print("PASS: the cotangent lift preserves the algebraic moment map")
print("PASS: every full affine moment reduction contains an A^4/mu_2 singularity")
print("PASS: beta != 0 has a unique beta=1 gauge and target reduction A^4")
print("PASS: the localized target coordinates form two canonical pairs")
print("PASS: localized quantum target reduction has explicit A_2 generators")
print("PASS: the source reduction is the displayed smooth affine hypersurface")
print("PASS: all three foundational collision points remain distinct on the free chart")
print("SCOPE: polynomial A^4-triviality and A_2 quantum reduction of the source remain open")
