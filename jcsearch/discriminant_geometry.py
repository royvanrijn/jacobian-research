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
