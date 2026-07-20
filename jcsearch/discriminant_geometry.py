"""Dual-curve geometry of the discriminant normalization."""

from __future__ import annotations

import sympy as sp


def discriminant_param(H, W):
    """Return the tangent-line coordinates (s,t) of the graph of H."""
    derivative = sp.diff(H, W)
    return sp.expand(derivative), sp.expand(W * derivative - H)


def cusp_polynomial(H, W):
    """Return the ramification polynomial of the discriminant normalization."""
    return sp.expand(sp.diff(H, W, 2))


def bitangent_equations(H, W, r, u):
    """Return divided equations for two graph points sharing a tangent line."""
    derivative = sp.diff(H, W)
    slope_equation = sp.cancel(
        (derivative.subs(W, r) - derivative.subs(W, u)) / (r - u)
    )
    intercept_equation = sp.cancel(
        (
            r * derivative.subs(W, r)
            - H.subs(W, r)
            - u * derivative.subs(W, u)
            + H.subs(W, u)
        )
        / (r - u)
    )
    return sp.expand(slope_equation), sp.expand(intercept_equation)


def symmetric_bitangent_equations(H, W, pair_sum, pair_product):
    """Return bitangent equations on unordered pairs of normalization points."""
    r, u = sp.symbols("_bitangent_r _bitangent_u")
    equations = bitangent_equations(H, W, r, u)
    result = []
    for equation in equations:
        symmetric, remainder, mapping = sp.symmetrize(
            equation, [r, u], formal=True
        )
        assert remainder == 0
        result.append(
            sp.expand(
                symmetric.subs(
                    {
                        mapping[0][0]: pair_sum,
                        mapping[1][0]: pair_product,
                    }
                )
            )
        )
    return tuple(result)


def ordinary_cusp_determinant(H, W):
    """Return det(nu''(W),nu'''(W)) modulo the cusp equation H''=0."""
    slope, intercept = discriminant_param(H, W)
    second = sp.Matrix([sp.diff(slope, W, 2), sp.diff(intercept, W, 2)])
    third = sp.Matrix([sp.diff(slope, W, 3), sp.diff(intercept, W, 3)])
    determinant = sp.expand(sp.det(sp.Matrix.hstack(second, third)))
    remainder = sp.rem(determinant, cusp_polynomial(H, W), W)
    return sp.factor(remainder)


def contact_incidence_dimension(degree: int, multiplicities):
    """Dimension bound for a common-tangent contact incidence modulo lines.

    A tangent line with distinct contact points of multiplicities ``m_i`` has
    ``H-line = product((W-r_i)**m_i) * Q``.  The returned dimension counts the
    contact points and the coefficients of Q.  ``None`` means that the total
    contact exceeds the polynomial degree, so the incidence is empty.
    """
    multiplicities = tuple(int(value) for value in multiplicities)
    if degree < 2 or not multiplicities or any(value < 2 for value in multiplicities):
        raise ValueError("contact multiplicities must all be at least two")
    total_contact = sum(multiplicities)
    if total_contact > degree:
        return None
    quotient_coefficients = degree - total_contact + 1
    return len(multiplicities) + quotient_coefficients


def tangent_chord_normalization(G, W, alpha, beta):
    """Normalize a tangent chord of G to the weighted endpoints zero and one.

    If the tangent at ``alpha`` also meets the graph at ``beta``, the result H
    satisfies H(0)=H'(0)=H(1)=0.  Its dual parameterization differs from that
    of G only by an affine source reparameterization and an invertible affine
    target transformation.
    """
    difference = sp.sympify(beta) - sp.sympify(alpha)
    if difference == 0:
        raise ValueError("the tangent-chord endpoints must be distinct")
    shifted = sp.sympify(alpha) + difference * W
    tangent = G.subs(W, alpha) + sp.diff(G, W).subs(W, alpha) * difference * W
    return sp.expand(G.subs(W, shifted) - tangent)


def deterministic_generic_primitive(degree: int, W):
    """Return the rational admissible audit seed used in degrees 3 and above."""
    if degree < 3:
        raise ValueError("inverse degree must be at least three")
    extra = sum((index + 1) * W**index for index in range(degree - 2))
    return sp.expand(W**2 * (1 - W) * extra)


def partition_dual_geometry(partition):
    """Annotate a full root-multiplicity partition by its dual-curve geometry."""
    multiplicities = tuple(sorted((int(value) for value in partition), reverse=True))
    nontrivial = tuple(value for value in multiplicities if value > 1)
    if len(nontrivial) == 1 and nontrivial[0] == sum(multiplicities):
        if nontrivial[0] == 3:
            return "ordinary cusp; maximally ramified"
        return "maximally ramified higher cusp"
    if nontrivial == (3, 2):
        return "ordinary cusp branch meeting another tangent branch"
    if nontrivial == (2, 2, 2):
        return "tritangent line / triple normalization point"
    if nontrivial == (3,):
        return "ordinary cusp"
    if nontrivial == (2, 2):
        return "ordinary node"
    if nontrivial == (4,):
        return "higher cusp"
    return "higher or multiple dual-curve singularity"
