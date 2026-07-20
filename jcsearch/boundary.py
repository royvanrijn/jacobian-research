"""Boundary saturation data for weighted inverse primitives."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class BoundarySaturationProfile:
    """Multiplicity data controlling the pulled-back discriminant at C=0."""

    degree: int
    zero_multiplicity: int
    one_multiplicity: int
    extra_degree: int
    extra_repetition_excess: int
    saturation_exponent: int
    zero_leading_coefficient: sp.Expr

    def boundary_trace(self, A, B, c=sp.Integer(1)):
        """Return the boundary trace of Q up to a nonzero scalar."""
        if self.zero_multiplicity == 2:
            conic = B**2 - 4 * sp.sympify(c) * self.zero_leading_coefficient * A
            return sp.expand(B**self.extra_repetition_excess * conic)
        return B**self.saturation_exponent


def boundary_saturation_profile(primitive, variable) -> BoundarySaturationProfile:
    """Compute the exact C-saturation exponent from the roots of ``H``.

    This applies to admissible primitives: the zero at 0 has multiplicity at
    least two and the distinguished zero at 1 is simple.  Extra factors may be
    simple or repeated and need not split over the coefficient field.
    """
    polynomial = sp.Poly(sp.sympify(primitive), variable)
    if polynomial.domain.characteristic() != 0:
        raise ValueError("boundary saturation requires characteristic zero")

    reduced = polynomial
    zero_multiplicity = 0
    zero_factor = sp.Poly(variable, variable)
    while reduced.eval(0) == 0:
        reduced = reduced.exquo(zero_factor)
        zero_multiplicity += 1

    one_multiplicity = 0
    one_factor = sp.Poly(variable - 1, variable)
    while reduced.eval(1) == 0:
        reduced = reduced.exquo(one_factor)
        one_multiplicity += 1

    if zero_multiplicity < 2:
        raise ValueError("an admissible primitive must have a double zero at 0")
    if one_multiplicity != 1:
        raise ValueError("the distinguished zero at 1 must be simple")

    extra_repetition_excess = sum(
        factor.degree() * (multiplicity - 1)
        for factor, multiplicity in reduced.sqf_list()[1]
    )
    saturation_exponent = zero_multiplicity + extra_repetition_excess
    zero_leading_coefficient = polynomial.coeff_monomial(variable**zero_multiplicity)

    return BoundarySaturationProfile(
        degree=polynomial.degree(),
        zero_multiplicity=zero_multiplicity,
        one_multiplicity=one_multiplicity,
        extra_degree=reduced.degree(),
        extra_repetition_excess=extra_repetition_excess,
        saturation_exponent=saturation_exponent,
        zero_leading_coefficient=zero_leading_coefficient,
    )
