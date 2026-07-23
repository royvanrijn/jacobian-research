#!/usr/bin/env python3
"""Exact counterexample: conductor data does not remove rerooting.

The generic quartic and quintic checks compare the complete decorations
computed by ``decorated_normalization``.  The first quartic Hessian collision
is computed directly over Q(sqrt(-2)) because conductor formation there must
not be inferred by specializing the ordinary cusp-node factorization.
"""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.discriminant_geometry import (  # noqa: E402
    bitangent_equations,
    decorated_normalization,
)


w = sp.symbols("w")


def monic_substitution(expression, scale):
    """Pull a divisor back by w -> scale*w and choose its monic generator."""
    return sp.Poly(sp.expand(expression.subs(w, scale * w)), w).monic().as_expr()


def polynomial_numerator(expression):
    """Clear coefficient-field denominators without changing an ideal."""
    return sp.cancel(expression).as_numer_denom()[0]


def ideal_maps_into(source, source_variables, target, target_variables, images):
    """Check that a polynomial substitution sends one ideal into another."""
    basis = sp.groebner(
        tuple(sp.expand(equation) for equation in target),
        *target_variables,
        order="lex",
    )
    substitution = dict(zip(source_variables, images))
    for equation in source:
        pulled = polynomial_numerator(sp.expand(equation.subs(substitution)))
        assert basis.reduce(pulled)[1] == 0


def assert_ideal_transport(
    left,
    left_variables,
    right,
    right_variables,
    left_in_right,
    right_in_left,
):
    """Check equality of two ideals under mutually inverse substitutions."""
    ideal_maps_into(
        left, left_variables, right, right_variables, left_in_right
    )
    ideal_maps_into(
        right, right_variables, left, left_variables, right_in_left
    )


def verify_decoration_transport(H, root, label, base_decoration=None):
    """Verify the full rerooted decorated-normalization square."""
    root = sp.sympify(root)
    kappa = sp.factor(-1 / (root * sp.diff(H, w).subs(w, root)))
    G = sp.factor(kappa * H.subs(w, root * w))

    assert sp.expand(G.subs(w, 0)) == 0
    assert sp.expand(sp.diff(G, w).subs(w, 0)) == 0
    assert sp.expand(G.subs(w, 1)) == 0
    assert sp.expand(sp.diff(G, w).subs(w, 1)) == -1

    left = base_decoration or decorated_normalization(
        H, w, prefix=f"{label}_left"
    )
    right = decorated_normalization(G, w, prefix=f"{label}_right")

    # Full Fitting multiplicities, node-branch projection, conductor, and the
    # intrinsic second-boundary center all pull back along w_H=root*w_G.
    assert sp.factor(
        right.fitting_polynomial
        - monic_substitution(left.fitting_polynomial, root)
    ) == 0
    assert sp.factor(
        right.node_branch_polynomial
        - monic_substitution(left.node_branch_polynomial, root)
    ) == 0
    assert sp.factor(
        right.conductor_polynomial
        - monic_substitution(left.conductor_polynomial, root)
    ) == 0
    assert sp.factor(
        right.boundary_divisor_polynomial
        - monic_substitution(left.boundary_divisor_polynomial, root)
    ) == 0

    # Ordered equal-image schemes, including diagonal saturation.
    lr, lu = left.ordered_variables
    rr, ru = right.ordered_variables
    assert_ideal_transport(
        left.off_diagonal_ideal,
        (lr, lu),
        right.off_diagonal_ideal,
        (rr, ru),
        (root * rr, root * ru),
        (lr / root, lu / root),
    )

    # Their unordered S_2 quotients.  Sum has weight one and product weight
    # two under normalization-line scaling.
    lsum, lproduct = left.pair_coordinates
    rsum, rproduct = right.pair_coordinates
    assert_ideal_transport(
        left.node_pair_ideal,
        (lsum, lproduct),
        right.node_pair_ideal,
        (rsum, rproduct),
        (root * rsum, root**2 * rproduct),
        (lsum / root, lproduct / root**2),
    )

    # Compare the actual finite conductor maps, not only their upstairs
    # generators.
    left_conductor = left.conductor_map
    right_conductor = right.conductor_map
    assert left_conductor is not None and right_conductor is not None
    ls, lt = left_conductor.downstairs_variables
    rs, rt = right_conductor.downstairs_variables
    assert_ideal_transport(
        left_conductor.downstairs_ideal,
        (ls, lt),
        right_conductor.downstairs_ideal,
        (rs, rt),
        (rs / (kappa * root), rt / kappa),
        (kappa * root * ls, kappa * lt),
    )

    left_slope, left_intercept = left_conductor.map_images
    right_slope, right_intercept = right_conductor.map_images
    assert sp.factor(
        left_slope.subs(w, root * w)
        - right_slope / (kappa * root)
    ) == 0
    assert sp.factor(
        left_intercept.subs(w, root * w)
        - right_intercept / kappa
    ) == 0

    return G, right


# Universal identities behind the isomorphism of decorated covers.
r, u, a, kappa = sp.symbols("r u a kappa", nonzero=True)
coefficients = sp.symbols("h0:7")
universal_H = sum(coefficients[index] * w**index for index in range(7))
universal_G = sp.expand(kappa * universal_H.subs(w, a * w))

assert sp.expand(
    sp.diff(universal_G, w, 2)
    - kappa * a**2 * sp.diff(universal_H, w, 2).subs(w, a * w)
) == 0

left_Ds, left_Dt = bitangent_equations(universal_H, w, a * r, a * u)
right_Ds, right_Dt = bitangent_equations(universal_G, w, r, u)
assert sp.expand(right_Ds - kappa * a**2 * left_Ds) == 0
assert sp.expand(right_Dt - kappa * a * left_Dt) == 0

s, t = sp.symbols("s t")
left_incidence = universal_H.subs(w, a * w) - (
    s / (kappa * a)
) * a * w + t / kappa
right_incidence = universal_G - s * w + t
assert sp.expand(right_incidence - kappa * left_incidence) == 0


# Generic quartic: q=2 reroots at a=1/3 to q=-2.
H4 = -sp.Rational(1, 2) * w**2 * (w - 1) * (3 * w - 1)
quartic = decorated_normalization(H4, w, prefix="quartic_base")
G4, quartic_rerooted = verify_decoration_transport(
    H4, sp.Rational(1, 3), "quartic", quartic
)
assert sp.factor(
    G4 - sp.Rational(1, 2) * w**2 * (w - 1) * (w - 3)
) == 0
assert quartic.ordinary and quartic_rerooted.ordinary
assert sp.degree(quartic.fitting_polynomial, w) == 2
assert sp.degree(quartic.node_branch_polynomial, w) == 2
assert sp.degree(quartic.conductor_polynomial, w) == 6


# Generic quintic: three distinct rerootings, with the complete ordinary
# three-cusp/three-node/six-branch conductor package.
H5 = -sp.Rational(1, 2) * w**2 * (w - 1) * (w - 2) * (w - 3)
quintic = decorated_normalization(H5, w, prefix="quintic_base")
assert quintic.ordinary
assert sp.discriminant(sp.diff(H5, w, 2), w) == 39960

quintic_rerootings = []
for root in (1, 2, 3):
    rerooted, decoration = verify_decoration_transport(
        H5, root, f"quintic_{root}", quintic
    )
    quintic_rerootings.append(sp.expand(rerooted))
    assert decoration.ordinary
    assert sp.degree(decoration.fitting_polynomial, w) == 3
    assert sp.degree(decoration.node_branch_polynomial, w) == 6
    assert sp.degree(decoration.conductor_polynomial, w) == 12

assert len(set(quintic_rerootings)) == 3
assert sp.factor(
    quintic_rerootings[1]
    - w**2 * (w - 1) * (2 * w - 1) * (2 * w - 3)
) == 0
assert sp.factor(
    quintic_rerootings[2]
    + sp.Rational(1, 2)
    * w**2
    * (w - 1)
    * (3 * w - 1)
    * (3 * w - 2)
) == 0


def saturated_off_diagonal(H, extension, prefix):
    """Compute (D_s,D_t):(r-u)^infinity over an algebraic field."""
    rr, uu, gate = sp.symbols(f"{prefix}_r {prefix}_u {prefix}_gate")
    Ds, Dt = bitangent_equations(H, w, rr, uu)
    basis = sp.groebner(
        (Ds, Dt, 1 - gate * (rr - uu)),
        gate,
        rr,
        uu,
        order="lex",
        extension=extension,
    )
    return tuple(
        polynomial.as_expr()
        for polynomial in basis.polys
        if gate not in polynomial.as_expr().free_symbols
    )


def exact_plane_conductor(H, extension, prefix):
    """Use the primitive implicit equation and adjunction over Q(extension)."""
    z, slope, intercept = sp.symbols(
        f"{prefix}_z {prefix}_s {prefix}_t"
    )
    Hz = H.subs(w, z)
    derivative = sp.diff(Hz, z)
    tangent_intercept = sp.expand(z * derivative - Hz)
    resultant = sp.resultant(
        slope - derivative, intercept - tangent_intercept, z
    )
    implicit = sp.Poly(
        resultant, slope, intercept, extension=extension
    ).primitive()[1].as_expr()

    H_derivative = sp.diff(H, w)
    parametrized_intercept = sp.expand(w * H_derivative - H)
    numerator = sp.Poly(
        sp.diff(implicit, intercept).subs(
            {
                slope: H_derivative,
                intercept: parametrized_intercept,
            }
        ),
        w,
        extension=extension,
    )
    hessian = sp.Poly(sp.diff(H, w, 2), w, extension=extension)
    conductor, remainder = numerator.div(hessian)
    assert remainder.is_zero
    return conductor.monic()


def exact_conductor_map(H, conductor, extension, prefix):
    """Present the finite conductor quotient over an algebraic field."""
    slope, intercept = sp.symbols(f"{prefix}_s {prefix}_t")
    derivative = sp.diff(H, w)
    tangent_intercept = sp.expand(w * derivative - H)
    basis = sp.groebner(
        (
            conductor.as_expr(),
            slope - derivative,
            intercept - tangent_intercept,
        ),
        w,
        slope,
        intercept,
        order="lex",
        extension=extension,
    )
    downstairs = tuple(
        polynomial.as_expr()
        for polynomial in basis.polys
        if w not in polynomial.as_expr().free_symbols
    )
    return (slope, intercept), downstairs, (derivative, tangent_intercept)


# First Hessian collision for normalized quartics.  Work over
# K=Q(q)/(q^2+2), represented by q=sqrt(2)*I.
q = sp.sqrt(2) * sp.I
H_collision = (
    -sp.Rational(1, 2)
    * w**2
    * (w - 1)
    * (q * w - q + w + 1)
)
collision_root = sp.cancel((q - 1) / (q + 1), extension=q)
collision_kappa = sp.cancel(
    -1 / (collision_root * sp.diff(H_collision, w).subs(w, collision_root)),
    extension=q,
)
G_collision = sp.expand(
    collision_kappa * H_collision.subs(w, collision_root * w)
)
expected_G_collision = (
    -sp.Rational(1, 2)
    * w**2
    * (w - 1)
    * ((-q) * w + q + w + 1)
)
assert sp.Poly(
    G_collision - expected_G_collision, w, extension=q
).is_zero

midpoint = sp.cancel(q / (2 * (q + 1)), extension=q)
collision_fitting = sp.Poly(
    sp.diff(H_collision, w, 2), w, extension=q
).monic()
assert collision_fitting == sp.Poly(
    (w - midpoint) ** 2, w, extension=q
)

# At the collision, saturation removes the diagonal higher-cusp point: the
# ordinary node pairing has disappeared rather than surviving spuriously.
assert saturated_off_diagonal(
    H_collision, q, "collision_left"
) == (1,)
assert saturated_off_diagonal(
    G_collision, q, "collision_right"
) == (1,)

collision_conductor = exact_plane_conductor(
    H_collision, q, "collision_left"
)
rerooted_collision_conductor = exact_plane_conductor(
    G_collision, q, "collision_right"
)
assert collision_conductor == sp.Poly(
    (w - midpoint) ** 6, w, extension=q
)
assert sp.Poly(
    collision_conductor.as_expr().subs(w, collision_root * w),
    w,
    extension=q,
).monic() == rerooted_collision_conductor
assert sp.Poly(
    collision_fitting.as_expr().subs(w, collision_root * w),
    w,
    extension=q,
).monic() == sp.Poly(
    sp.diff(G_collision, w, 2), w, extension=q
).monic()

# The nonreduced downstairs conductor schemes and their maps commute with the
# same target scaling.  Coefficient cancellation is explicitly performed in
# Q(sqrt(-2)) so no numerical splitting or support reduction is involved.
left_variables, left_ideal, left_images = exact_conductor_map(
    H_collision, collision_conductor, q, "collision_map_left"
)
right_variables, right_ideal, right_images = exact_conductor_map(
    G_collision,
    rerooted_collision_conductor,
    q,
    "collision_map_right",
)
left_basis = sp.groebner(
    left_ideal, *left_variables, order="lex", extension=q
)
right_basis = sp.groebner(
    right_ideal, *right_variables, order="lex", extension=q
)
ls, lt = left_variables
rs, rt = right_variables
for equation in left_ideal:
    pulled = equation.subs(
        {
            ls: rs / (collision_kappa * collision_root),
            lt: rt / collision_kappa,
        }
    )
    numerator = sp.cancel(pulled, extension=q).as_numer_denom()[0]
    assert right_basis.reduce(numerator)[1] == 0
for equation in right_ideal:
    pulled = equation.subs(
        {
            rs: collision_kappa * collision_root * ls,
            rt: collision_kappa * lt,
        }
    )
    numerator = sp.cancel(pulled, extension=q).as_numer_denom()[0]
    assert left_basis.reduce(numerator)[1] == 0

assert sp.Poly(
    left_images[0].subs(w, collision_root * w)
    - right_images[0] / (collision_kappa * collision_root),
    w,
    extension=q,
).is_zero
assert sp.Poly(
    left_images[1].subs(w, collision_root * w)
    - right_images[1] / collision_kappa,
    w,
    extension=q,
).is_zero


print("PASS rerooting: the universal incidence and equal-image ideals transport")
print("PASS quartic: two distinct ordinary rerootings have isomorphic full decorations")
print("PASS quintic: three distinct ordinary rerootings have isomorphic full decorations")
print("PASS conductor: upstairs ideals, downstairs quotients, and finite maps transport")
print("PASS collision: Fitt=(w-m)^2, node pairing is empty, conductor=(w-m)^6")
print("COUNTEREXAMPLE: conductor and full Fitting data do not eliminate rerooting")
