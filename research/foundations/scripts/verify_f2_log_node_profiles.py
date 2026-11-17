#!/usr/bin/env python3
"""Verify the local logarithmic node profiles of the F2 ``(75,125)`` row."""

from __future__ import annotations

from pathlib import Path
import sys

from sympy import Matrix, Poly, Rational, cancel, fraction, symbols

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.log_node_profiles import (
    close_smooth_profile_by_boundary_support,
    smooth_target_boundary_profile,
    target_node_profile,
)


def factor_order(expression, factor, variable) -> int:
    """Return the order of a rational expression along an irreducible factor."""

    numerator, denominator = fraction(cancel(expression))

    def multiplicity(polynomial) -> int:
        value = Poly(polynomial, variable, domain="QQ")
        divisor = Poly(factor, variable, domain="QQ")
        answer = 0
        while value.rem(divisor).is_zero:
            value = value.exquo(divisor)
            answer += 1
        return answer

    return multiplicity(numerator) - multiplicity(denominator)


def infinity_order(expression, variable) -> int:
    """Return the order in ``w=1/variable`` at infinity."""

    numerator, denominator = fraction(cancel(expression))
    return Poly(denominator, variable).degree() - Poly(numerator, variable).degree()


def main() -> None:
    s = symbols("s")
    A = 1 + 3 * s + Rational(9, 5) * s**2
    h = 125 * s * (1 + s) ** 5 / (5 + 15 * s + 9 * s**2) ** 3
    eta = cancel(1 / h)

    # At h=0 the target-node chart is
    #   pi_0=a/b^2=(-Q)/P^2=tau_0*A/(1+s)^2,  xi_0=h.
    pi_zero_coefficient = cancel(A / (1 + s) ** 2)
    assert factor_order(pi_zero_coefficient, s, s) == 0
    assert factor_order(h, s, s) == 1
    assert factor_order(pi_zero_coefficient, s + 1, s) == -2
    assert factor_order(h, s + 1, s) == 5

    # At h=infinity the opposite node chart is
    #   pi_inf=b^3/a=P^3/(-Q)^2=tau_inf*(1+s)^3/A^2,
    #   xi_inf=eta=1/h.
    pi_infinity_coefficient = cancel((1 + s) ** 3 / A**2)
    assert factor_order(pi_infinity_coefficient, A, s) == -2
    assert factor_order(eta, A, s) == 3

    # The other endpoint maps to a smooth point of the extracted target
    # component.  Its normal coefficient vanishes once and its residue map
    # has contact order three.
    eta_infinity = Rational(729, 125)
    assert infinity_order(pi_infinity_coefficient, s) == 1
    assert infinity_order(eta - eta_infinity, s) == 3

    profiles = (
        target_node_profile(
            "s=0", transverse_order=0, residue_index=1
        ),
        target_node_profile(
            "s=-1", transverse_order=-2, residue_index=5
        ),
        target_node_profile(
            "first root of A", transverse_order=-2, residue_index=3
        ),
        target_node_profile(
            "second root of A", transverse_order=-2, residue_index=3
        ),
    )

    assert profiles[0].source_blowups == 0
    assert all(profile.cokernel_model == "0" for profile in profiles)
    assert all(profile.fitting_one_is_unit for profile in profiles)
    assert all(profile.normalization_defect == "zero" for profile in profiles)
    assert [profile.residual_determinant for profile in profiles] == [1, 5, 3, 3]
    assert [profile.source_blowups for profile in profiles] == [0, 2, 2, 2]
    assert [profile.resolved_exponent_matrix for profile in profiles] == [
        ((1, 0), (0, 1)),
        ((1, 0), (0, 5)),
        ((1, 0), (0, 3)),
        ((1, 0), (0, 3)),
    ]

    smooth_endpoint = smooth_target_boundary_profile(
        "s=infinity", transverse_order=1, contact_index=3
    )
    assert smooth_endpoint.source_blowups == 0
    assert smooth_endpoint.resolved_exponent_matrix == ((1, 1), (0, 3))
    assert smooth_endpoint.fitting_one_is_unit
    assert smooth_endpoint.residual_determinant is None
    assert "not determined" in smooth_endpoint.normalization_defect

    # Keller etaleness on the affine open forces the logarithmic determinant
    # divisor to be supported on the source boundary.  Its generic order on
    # the terminal component is zero because the residue map is generically
    # etale.  Hence the order-three restriction at this endpoint is the whole
    # local determinant: no transverse deformation of its support is possible.
    closed_smooth_endpoint = close_smooth_profile_by_boundary_support(
        smooth_endpoint, terminal_generic_order=0
    )
    assert closed_smooth_endpoint.cokernel_model == "R/(w^3)"
    assert closed_smooth_endpoint.normalization_defect == (
        "zero on the smooth reduced support w=0"
    )

    # The three interior points each need the valuation ray (2,1), hence two
    # source blowups in the fixed target-node chart.  A one-blowup attachment
    # only reaches (1,1) and leaves transverse target order -1.
    assert sum(profile.source_blowups for profile in profiles[1:]) == 6
    assert all(profile.transverse_order + 1 == -1 for profile in profiles[1:])

    # The second node in each two-exceptional chain maps into the next cone
    # of the already extracted target fan.  In local target cocharacter
    # coordinates the source rays are terminal=(1,0), inner=(2,1), and
    # outer=(1,1).  Expressing the last two images in the relevant target
    # cone gives tame triangular matrices once again.
    def second_node_matrix(index: int, next_target_ray: tuple[int, int]) -> Matrix:
        valuation_map = Matrix([[1, -2], [0, index]])
        inner_image = valuation_map * Matrix([2, 1])
        outer_image = valuation_map * Matrix([1, 1])
        target_basis = Matrix.hstack(Matrix([0, 1]), Matrix(next_target_ray))
        return Matrix.hstack(
            target_basis.inv() * inner_image,
            target_basis.inv() * outer_image,
        )

    h_zero_second = second_node_matrix(5, (-1, 3))
    h_infinity_second = second_node_matrix(3, (-1, 2))
    assert h_zero_second == Matrix([[5, 2], [0, 1]])
    assert h_infinity_second == Matrix([[3, 1], [0, 1]])
    assert h_zero_second.det() == 5
    assert h_infinity_second.det() == 3

    print(
        "PASS: F2 terminal log profiles; target-node cokernels vanish and "
        "the smooth endpoint is the cyclic boundary thickening R/(w^3)"
    )


if __name__ == "__main__":
    main()
