#!/usr/bin/env python3
"""Exact Fitting, self-intersection, conductor, and quartic-island audit."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.discriminant_geometry import (  # noqa: E402
    bitangent_equations,
    decorated_normalization,
    deterministic_generic_primitive,
    paired_quadratic_divisor_invariants,
    plane_parametric_conductor_polynomial,
    pointed_quadratic_divisor_invariant,
    quadratic_zero_cluster_boundary_chart,
)
from jcsearch.boundary import (  # noqa: E402
    cancellation_boundary_intersection_profile,
    cancellation_boundary_cover_profile,
    weighted_boundary_cover_profile,
)


W = sp.symbols("W")

SEEDS = {
    "F4a": W**3 * (W - 1),
    "F4b": W**2 * (W - 1) * (W - 3),
    "F4c": -W**2 * (W - 1) * (2 * W + 1),
}

EXPECTED = {
    "F4a": {
        "fitting": W * (W - sp.Rational(1, 2)),
        "pair_sum": sp.Rational(1, 2),
        "pair_product": -sp.Rational(1, 8),
        "pointed": 1,
    },
    "F4b": {
        "fitting": W**2 - 2 * W + sp.Rational(1, 2),
        "pair_sum": 2,
        "pair_product": -sp.Rational(1, 2),
        "pointed": 2,
    },
    "F4c": {
        "fitting": W**2 - sp.Rational(1, 4) * W - sp.Rational(1, 12),
        "pair_sum": sp.Rational(1, 4),
        "pair_product": -sp.Rational(9, 32),
        "pointed": sp.Rational(3, 19),
    },
}

decorations = {}
for name, H in SEEDS.items():
    decoration = decorated_normalization(
        H,
        W,
        prefix=name.lower(),
    )
    decorations[name] = decoration
    expected = EXPECTED[name]

    assert decoration.ordinary
    assert decoration.node_scheme_transverse
    assert decoration.infinity_mark == "infinity"
    assert decoration.boundary_marks == (("zero_cluster", 0),)
    assert decoration.boundary_divisor_polynomial in (W, W**2)
    assert sp.factor(
        decoration.fitting_polynomial - expected["fitting"]
    ) == 0
    conductor_map = decoration.conductor_map
    assert conductor_map is not None
    assert sp.degree(conductor_map.upstairs_polynomial, W) == 6
    downstairs_s, downstairs_t = conductor_map.downstairs_variables
    assert len(conductor_map.downstairs_ideal) == 2
    assert any(sp.degree(equation, downstairs_s) == 1 for equation in conductor_map.downstairs_ideal)
    assert any(sp.degree(equation, downstairs_t) == 3 for equation in conductor_map.downstairs_ideal)
    for equation in conductor_map.downstairs_ideal:
        pulled = sp.expand(
            equation.subs(
                dict(zip(conductor_map.downstairs_variables, conductor_map.map_images))
            )
        )
        assert sp.rem(
            pulled, conductor_map.upstairs_polynomial, W
        ) == 0
    assert sp.degree(decoration.node_branch_polynomial, W) == 2
    assert sp.gcd(
        decoration.fitting_polynomial,
        decoration.node_branch_polynomial,
    ) == 1
    assert sp.factor(
        decoration.conductor_polynomial
        - decoration.fitting_polynomial**2
        * decoration.node_branch_polynomial
    ) == 0

    pair_sum, pair_product = decoration.pair_coordinates
    pair_solution = sp.solve(
        decoration.node_pair_ideal,
        (pair_sum, pair_product),
        dict=True,
    )
    assert pair_solution == [
        {
            pair_sum: expected["pair_sum"],
            pair_product: expected["pair_product"],
        }
    ]
    # The unordered equations really are the S_2 quotient of the ordered
    # saturated scheme, rather than a separately matched support calculation.
    ordered_r, ordered_u = decoration.ordered_variables
    ordered_basis = sp.groebner(
        decoration.off_diagonal_ideal,
        ordered_r,
        ordered_u,
        order="lex",
    )
    for equation in decoration.node_pair_ideal:
        pullback = sp.expand(
            equation.subs(
                {
                    pair_sum: ordered_r + ordered_u,
                    pair_product: ordered_r * ordered_u,
                }
            )
        )
        assert ordered_basis.reduce(pullback)[1] == 0
    ordered_solutions = sp.solve(
        decoration.off_diagonal_ideal,
        (ordered_r, ordered_u),
        dict=True,
    )
    assert len(ordered_solutions) == 2
    pointed = pointed_quadratic_divisor_invariant(
        decoration.fitting_polynomial, W, 0
    )
    assert pointed == expected["pointed"]

# The cusp divisor plus node pairing alone does not separate the two split
# examples: the squared node radius is three times the squared cusp radius in
# both cases.  The canonical zero-cluster mark does separate them.
def quadratic_radius_squared(polynomial):
    poly = sp.Poly(polynomial, W)
    a, _, _ = poly.all_coeffs()
    return sp.factor(sp.discriminant(poly.as_expr(), W) / (4 * a**2))


for name in ("F4b", "F4c"):
    decoration = decorations[name]
    pair_sum, pair_product = decoration.pair_coordinates
    solution = sp.solve(
        decoration.node_pair_ideal,
        (pair_sum, pair_product),
        dict=True,
    )[0]
    node_polynomial = sp.expand(
        W**2 - solution[pair_sum] * W + solution[pair_product]
    )
    ratio = sp.factor(
        quadratic_radius_squared(node_polynomial)
        / quadratic_radius_squared(decoration.fitting_polynomial)
    )
    assert ratio == 3
    assert paired_quadratic_divisor_invariants(
        decoration.fitting_polynomial, node_polynomial, W
    ) == (0, 3)

assert EXPECTED["F4b"]["pointed"] != EXPECTED["F4c"]["pointed"]

# The mark r=0 is not a coordinate convention: in the normalized
# discriminant boundary chart, W=C*k is forced and specializes to the full
# boundary conic.  Thus its center on the normalization line is W=0.
k = sp.symbols("k")
for name in ("F4b", "F4c"):
    H = SEEDS[name]
    h2 = sp.Poly(H, W).coeff_monomial(W**2)
    boundary_A, boundary_B = quadratic_zero_cluster_boundary_chart(
        H, W, parameter=k
    )
    assert boundary_A == h2 * k**2
    assert boundary_B == 2 * h2 * k
    assert sp.factor(boundary_B**2 - 4 * h2 * boundary_A) == 0

# Multiplicity regression: Fitt_0 retains the doubled Hessian point instead
# of silently passing to reduced support.
H_double = W**2 * (W - 1) * (2 * W**2 - 2 * W + 1)
double_decoration = decorated_normalization(H_double, W, prefix="double_hessian")
assert dict(double_decoration.fitting_factorization) == {
    2 * W - 1: 2,
    5 * W - 1: 1,
}
assert not double_decoration.ordinary
assert double_decoration.conductor_map is not None
assert sp.factor(double_decoration.conductor_polynomial) == sp.factor(
    (2 * W - 1) ** 6
    * (5 * W - 1) ** 2
    * (200 * W**4 - 280 * W**3 + 128 * W**2 - 20 * W + 1)
    / (2**6 * 5**2 * 200)
)

# No arbitrary polynomial receives a fictitious mark at zero.  Its entire
# intrinsic second-boundary center scheme is gcd(H,H'), and a point label is
# inferred only when that scheme has one certified support point.
unmarked = decorated_normalization(W**4 + W + 1, W, prefix="unmarked")
assert unmarked.boundary_divisor_polynomial == 1
assert unmarked.boundary_marks == ()

# A fully ordinary quintic exercises the complete higher-degree package:
# three cusp points, three unordered nodes, six paired node branches and the
# general conductor.  Its irreducible cubic Hessian also tests that residue
# fields are retained over Q rather than split numerically.
H_quintic = deterministic_generic_primitive(5, W)
quintic = decorated_normalization(H_quintic, W, prefix="generic_quintic")
assert quintic.ordinary and quintic.node_scheme_transverse
assert sp.degree(quintic.fitting_polynomial, W) == 3
assert len(quintic.fitting_factorization) == 1
assert sp.degree(quintic.fitting_factorization[0][0], W) == 3
assert sp.degree(quintic.node_branch_polynomial, W) == 6
assert sp.degree(quintic.conductor_polynomial, W) == 12
pair_sum, pair_product = quintic.pair_coordinates
pair_basis = sp.groebner(
    quintic.node_pair_ideal, pair_sum, pair_product, order="lex"
)
product_eliminants = [
    polynomial.as_expr()
    for polynomial in pair_basis.polys
    if polynomial.as_expr().free_symbols <= {pair_product}
]
assert len(product_eliminants) == 1
assert sp.degree(product_eliminants[0], pair_product) == 3
implicit_s, implicit_t = quintic.implicit_variables
quintic_slope = sp.diff(H_quintic, W)
quintic_intercept = sp.expand(W * quintic_slope - H_quintic)
assert sp.expand(
    quintic.implicit_equation.subs(
        {implicit_s: quintic_slope, implicit_t: quintic_intercept}
    )
) == 0

# Affine reparameterization transports every divisor ideal.  This is the
# normalization-line shadow of left-right equivalence and also checks that
# multiplicities, the conductor and the boundary center move together.
V = sp.symbols("V")
a, b = sp.Rational(3, 2), sp.Rational(-2, 3)
H_reparameterized = sp.expand(H_quintic.subs(W, a * V + b))
transported = decorated_normalization(
    H_reparameterized,
    V,
    boundary_marks={"transported_center": -b / a},
    prefix="transported_quintic",
)
expected_fitting = sp.Poly(
    quintic.fitting_polynomial.subs(W, a * V + b), V
).monic().as_expr()
expected_conductor = sp.Poly(
    quintic.conductor_polynomial.subs(W, a * V + b), V
).monic().as_expr()
assert sp.factor(transported.fitting_polynomial - expected_fitting) == 0
assert sp.factor(transported.conductor_polynomial - expected_conductor) == 0

# A nonlinear polynomial target automorphism (s,t)->(s,t+s^2) leaves the
# normalization morphism, its equal-image relation, Fitting divisor and
# conductor unchanged.  This is a direct positive left-right-invariance
# regression rather than only a comparison of closed formulas.
x = quintic_slope
y = quintic_intercept
_, _, nonlinear_conductor = plane_parametric_conductor_polynomial(
    x, y + x**2, W, prefix="nonlinear_target"
)
assert nonlinear_conductor == quintic.conductor_polynomial
r, u = sp.symbols("nonlinear_r nonlinear_u")
Ds, Dt = bitangent_equations(H_quintic, W, r, u)
x_r, x_u = x.subs(W, r), x.subs(W, u)
nonlinear_Dt = sp.expand(Dt + Ds * (x_r + x_u))
original_basis = sp.groebner((Ds, Dt), r, u, order="lex")
nonlinear_basis = sp.groebner((Ds, nonlinear_Dt), r, u, order="lex")
assert original_basis.reduce(nonlinear_Dt)[1] == 0
assert nonlinear_basis.reduce(Dt)[1] == 0

# Generic divisorial layer of the upstairs finite normalization cover.  The
# completed tame-DVR model t=unit*q^e has different exponent e-1 and inertia
# an e-cycle.  The degree sums include affine and boundary primes.
for name, H in SEEDS.items():
    cover = weighted_boundary_cover_profile(H, W)
    assert cover.degree == 4
    assert all(boundary + affine == total == 4 for _, boundary, affine, total in cover.degree_sums)
    delta = [prime for prime in cover.primes if prime.target_divisor == "Z_Delta"]
    assert len(delta) == 1
    assert (delta[0].ramification_index, delta[0].residue_degree) == (2, 1)
    assert delta[0].different_exponent == 1
    assert delta[0].inertia_cycle == (2, 1, 1)

    second = [prime for prime in cover.primes if prime.target_divisor == "Z_0"]
    if name == "F4a":
        assert [(prime.label, prime.ramification_index) for prime in second] == [
            ("E_zero_cluster", 2)
        ]
        assert second[0].different_exponent == 1
    else:
        assert len(second) == 1
        assert second[0].ramification_index == 1
        assert second[0].different_exponent == 0

repeated_cover = weighted_boundary_cover_profile(
    sp.Rational(1, 4) * W**2 * (1 - W) * (W + 1) ** 2,
    W,
)
repeated_extra = [
    prime for prime in repeated_cover.primes if prime.label == "E_extra_0"
][0]
assert repeated_extra.ramification_index == 2
assert repeated_extra.different_exponent == 1
assert repeated_extra.local_parameter_equation == "C = unit*q^2"

for m, r in ((1, 1), (1, 2), (2, 1), (2, 2)):
    cover = cancellation_boundary_cover_profile(m, r)
    degree = r * (m + 1) + 1
    assert cover.degree == degree
    assert all(boundary + affine == total == degree for _, boundary, affine, total in cover.degree_sums)
    delta = cover.primes[0]
    assert delta.ramification_index == r + 1
    assert delta.different_exponent == r
    assert delta.inertia_cycle == (r + 1,) + (1,) * (degree - r - 1)
    second = [prime for prime in cover.primes if prime.target_divisor == "Z_0"]
    assert len(second) == (0 if (m, r) == (1, 1) else m * r - 1)
    assert all(
        prime.ramification_index == 1 and prime.different_exponent == 0
        for prime in second
    )
    intersection = cancellation_boundary_intersection_profile(m, r)
    assert intersection.nilpotency_index == m * r * (m + 1)
    assert sum(
        branch.length_contribution for branch in intersection.branches
    ) == intersection.nilpotency_index
    assert [
        (branch.divisor_multiplicity, branch.contact_order)
        for branch in intersection.branches
    ] == [(m, m * r), (1, m * r)]

print("PASS: Fitt_0 Omega is extracted with its full Hessian-root multiplicities")
print("PASS: diagonal saturation recovers the exact ordered and unordered node schemes")
print("PASS: the ordinary conductor is cusp^2 times the node-branch divisor")
print("PASS: the implicit-equation conductor handles nonordinary singularities")
print("PASS: the finite conductor map retains the downstairs singular scheme and upstairs branches")
print("PASS: cusp/node data alone match F4b and F4c, but the boundary mark separates them")
print("PASS: intrinsic marks, a generic quintic, residue fields, and target conjugation are exact")
print("PASS: weighted and cancellation boundary profiles retain (e,f), differents, DVR models, and inertia")
print("PASS: cancellation higher-stratum contacts retain both completed normalization branches")
