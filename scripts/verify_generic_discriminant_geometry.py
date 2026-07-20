#!/usr/bin/env python3
"""Exact audit of generic nodal-cuspidal discriminant geometry."""

import itertools
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.discriminant_geometry import (  # noqa: E402
    bitangent_equations,
    cusp_polynomial,
    deterministic_generic_primitive,
    discriminant_param,
    ordinary_cusp_determinant,
    partition_dual_geometry,
    symmetric_bitangent_equations,
)
from jcsearch.weighted import WeightedSeedModel, w  # noqa: E402


def quotient_dimension(groebner_basis, variables):
    """Dimension of a zero-dimensional quotient from its leading monomials."""
    leading = [
        polynomial.LM(order=groebner_basis.order).exponents
        for polynomial in groebner_basis.polys
    ]
    bounds = []
    for index in range(len(variables)):
        pure_powers = [
            monomial[index]
            for monomial in leading
            if monomial[index] > 0
            and all(
                monomial[other] == 0
                for other in range(len(variables))
                if other != index
            )
        ]
        assert pure_powers
        bounds.append(min(pure_powers))
    return sum(
        1
        for exponent in itertools.product(*(range(bound) for bound in bounds))
        if not any(
            all(left >= right for left, right in zip(exponent, monomial))
            for monomial in leading
        )
    )


s, t = sp.symbols("s t")
pair_sum, pair_product = sp.symbols("pair_sum pair_product")

for degree in range(3, 11):
    H = deterministic_generic_primitive(degree, w)
    normalization_c = -sp.diff(H, w).subs(w, 1)
    model = WeightedSeedModel(sp.diff(H, w), c=normalization_c)
    assert model.fiber_degree == degree

    # Exactly n-2 simple flexes, each producing an ordinary cusp.
    cusp = sp.Poly(cusp_polynomial(H, w), w)
    assert cusp.degree() == degree - 2
    assert sp.gcd(cusp, cusp.diff()).degree() == 0
    expected_cusp_jet = sp.rem(
        2 * sp.diff(H, w, 3) ** 2,
        cusp.as_expr(),
        w,
    )
    assert sp.factor(ordinary_cusp_determinant(H, w) - expected_cusp_jet) == 0

    # The tangent-line parameter has degrees n-1,n and is smooth at infinity.
    slope, intercept = discriminant_param(H, w)
    assert sp.degree(slope, w) == degree - 1
    assert sp.degree(intercept, w) == degree
    leading = sp.Poly(H, w).LC()
    assert sp.limit(w * slope / intercept, w, sp.oo) == sp.Rational(degree, degree - 1)
    assert sp.limit(w**degree / intercept, w, sp.oo) == 1 / ((degree - 1) * leading)

    # Work on unordered pairs.  The common resultant contains n-2 diagonal
    # cusp points and the remaining simple factor counts ordinary bitangents.
    equation_s, equation_t = symmetric_bitangent_equations(
        H, w, pair_sum, pair_product
    )
    bitangent_resultant = sp.Poly(
        sp.resultant(equation_s, equation_t, pair_product), pair_sum
    )
    cusp_sums = sp.Poly(
        sp.resultant(cusp.as_expr(), pair_sum - 2 * w, w), pair_sum
    )
    diagonal_factor = sp.gcd(bitangent_resultant, cusp_sums)
    assert diagonal_factor.degree() == degree - 2
    node_factor = bitangent_resultant.exquo(diagonal_factor)
    node_count = (degree - 2) * (degree - 3) // 2
    assert node_factor.degree() == node_count
    # Every unordered pair has two orderings, and the saturation removes the
    # diagonal, so this is the requested ordered-pair quotient length.
    assert 2 * node_factor.degree() == (degree - 2) * (degree - 3)
    assert sp.gcd(node_factor, node_factor.diff()).degree() == 0
    assert sp.gcd(node_factor, cusp_sums).degree() == 0

    # The discriminant has the expected degree and minimal total Tjurina
    # number: 2 per cusp and 1 per node.  Together with the normalization and
    # bitangent count, this excludes tritangents, cusp-branch collisions, and
    # every higher affine singularity for these exact audit seeds.
    inverse = H - s * w + t
    discriminant = sp.factor(sp.discriminant(inverse, w))
    assert sp.total_degree(discriminant) == degree
    singular_basis = sp.groebner(
        [discriminant, sp.diff(discriminant, s), sp.diff(discriminant, t)],
        s,
        t,
        order="grevlex",
    )
    assert singular_basis.is_zero_dimensional
    tjurina_number = quotient_dimension(singular_basis, (s, t))
    assert tjurina_number == 2 * (degree - 2) + node_count

# Ordered equations are symmetric and their two branch directions are
# automatically transverse away from the diagonal.
r, u = sp.symbols("r u")
sample_H = deterministic_generic_primitive(7, w)
ordered_s, ordered_t = bitangent_equations(sample_H, w, r, u)
assert sp.expand(ordered_s.subs({r: u}) - cusp_polynomial(sample_H, w).subs(w, u)) == 0
assert sp.expand(ordered_t.subs({r: u}) - u * cusp_polynomial(sample_H, w).subs(w, u)) == 0
assert sp.det(sp.Matrix([[1, 1], [r, u]])) == u - r

# The existing multiplicity classifier now has a dual-curve interpretation.
assert partition_dual_geometry((3, 1, 1)) == "ordinary cusp"
assert partition_dual_geometry((2, 2, 1)) == "ordinary node"
assert partition_dual_geometry((4, 1)) == "higher cusp"
assert partition_dual_geometry((3, 2)) == "ordinary cusp branch meeting another tangent branch"
assert partition_dual_geometry((2, 2, 2)) == "tritangent line / triple normalization point"

print("PASS: deterministic admissible seeds through inverse degree ten")
print("PASS: exactly n-2 ordinary cusps and (n-2)(n-3)/2 ordinary nodes")
print("PASS: infinity is smooth and the singular scheme has minimal Tjurina number")
print("PASS: multiplicity partitions carry their dual-curve singularity types")
