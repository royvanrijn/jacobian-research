#!/usr/bin/env python3
"""Verify the upstream carrier-extraction log profile for F2 ``(75,125)``."""

from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "plane-jc" / "cas"))

from jcsearch.log_node_profiles import (  # noqa: E402
    cyclic_snc_matching_profile,
    determinant,
)
from verify_f2_75_125_carrier_wronskian import (  # noqa: E402
    carrier_target_fan_audit,
)


def allowed_support(degree: int, edge_order: int) -> tuple[tuple[int, int], ...]:
    """Polynomial exponents above one certified carrier edge."""

    return tuple(
        (x_degree, y_degree)
        for x_degree in range(degree + 1)
        for y_degree in range(degree - x_degree + 1)
        if -5 * x_degree + y_degree >= -edge_order
    )


def support_audit() -> None:
    p_support = allowed_support(75, 15)
    q_support = allowed_support(125, 25)

    for index in range(1, 5):
        p_minimum = min(-index * x_degree + y_degree for x_degree, y_degree in p_support)
        q_minimum = min(-index * x_degree + y_degree for x_degree, y_degree in q_support)
        p_minimizers = tuple(
            exponent
            for exponent in p_support
            if -index * exponent[0] + exponent[1] == p_minimum
        )
        q_minimizers = tuple(
            exponent
            for exponent in q_support
            if -index * exponent[0] + exponent[1] == q_minimum
        )
        assert (p_minimum, p_minimizers) == (-3 * index, ((3, 0),))
        assert (q_minimum, q_minimizers) == (-5 * index, ((5, 0),))

    assert max(x_degree for x_degree, _ in p_support) == 15
    assert max(x_degree for x_degree, _ in q_support) == 25
    assert tuple(exponent for exponent in p_support if exponent[0] == 15) == (
        (15, 60),
    )
    assert tuple(exponent for exponent in q_support if exponent[0] == 25) == (
        (25, 100),
    )
    # Other top-total-degree monomials may lie strictly above the carrier
    # edge.  In the (W,U) node chart they acquire positive W order; the two
    # displayed endpoints are the unique terms of W order zero and therefore
    # make the normalized local coefficients units.
    assert (15, 60) in p_support and sum((15, 60)) == 75
    assert (25, 100) in q_support and sum((25, 100)) == 125


def upstream_ladder_audit() -> None:
    # On the carrier-zero side use a=1/x and y.  The first target functions
    # have leading monomials pi=a and w=a^7*y.  Thus the exponent map is
    # unimodular and sends the complete source ladder to the target ladder.
    exponent_map = ((1, 0), (7, 1))
    assert determinant(exponent_map[0], exponent_map[1]) == 1
    source_rays = tuple((index, 1) for index in range(0, 6))
    target_rays = tuple(
        (
            exponent_map[0][0] * ray[0] + exponent_map[0][1] * ray[1],
            exponent_map[1][0] * ray[0] + exponent_map[1][1] * ray[1],
        )
        for ray in source_rays
    )
    assert target_rays == (
        (0, 1),
        (1, 8),
        (2, 15),
        (3, 22),
        (4, 29),
        (5, 36),
    )
    assert all(
        abs(determinant(left, right)) == 1
        for left, right in zip(source_rays, source_rays[1:])
    )
    assert all(
        abs(determinant(left, right)) == 1
        for left, right in zip(target_rays, target_rays[1:])
    )
    target_fan = {
        tuple(ray) for ray in carrier_target_fan_audit()["boundary_rays"]
    }
    assert set(target_rays[1:]).issubset(target_fan)


def infinity_node_audit() -> None:
    # At the node of the first exceptional W=0 and the strict line at
    # infinity U=0, x=(W*U)^-1 and y=U^-1.  The certified edge and degree
    # bounds make the normalized P and -Q coefficients units:
    #   P=W^-15*U^-75*A,  -Q=W^-25*U^-125*B.
    p_orders = (-15, -75)
    q_orders = (-25, -125)
    pi_orders = (
        3 * p_orders[0] - 2 * q_orders[0],
        3 * p_orders[1] - 2 * q_orders[1],
    )
    assert pi_orders == (5, 25)

    # d(pi) wedge d(w)=(unit)*P^7/(-Q)^6*dx wedge dy, while
    # dx wedge dy=(unit)*W^-2*U^-3*dW wedge dU.
    wedge_orders = (
        7 * p_orders[0] - 6 * q_orders[0] - 2,
        7 * p_orders[1] - 6 * q_orders[1] - 3,
    )
    assert wedge_orders == (43, 222)

    # PF2CLP1 fixes ord_W(w)=36.  After absorbing the unit in pi into W, the
    # monomial differential operator sends W^r*U^s to a nonzero multiple of
    # W^(r+4)*U^(s+24), except on its kernel s=5*r.  Since the exact wedge is
    # a unit times W^43*U^222, every r=36,37,38 term must be a kernel term.
    # The nonzero r=36 coefficient is therefore forced to have s=180.  The
    # first nonshear correction is forced by
    #   (r+4,s+24)=(43,222).
    w_leading = (36, 180)
    assert determinant(pi_orders, w_leading) == 0
    pre_correction_kernel = tuple((ray, 5 * ray) for ray in range(36, 40))
    assert pre_correction_kernel == (
        (36, 180),
        (37, 185),
        (38, 190),
        (39, 195),
    )
    assert all(
        determinant(pi_orders, exponent) == 0
        for exponent in pre_correction_kernel
    )
    correction = (wedge_orders[0] - 4, wedge_orders[1] - 24)
    assert correction == (39, 198)
    assert determinant(pi_orders, correction) == 15
    residue_correction = (
        correction[0] - w_leading[0],
        correction[1] - w_leading[1],
    )
    assert residue_correction == (3, 18)

    # Primitive target normal and tangential residue after a harmless target
    # shear.  The local logarithmic determinant is 3*W^3*U^18.
    target_normal = (
        29 * pi_orders[0] - 4 * w_leading[0],
        29 * pi_orders[1] - 4 * w_leading[1],
    )
    assert target_normal == (1, 5)
    local_exponent_determinant = determinant(target_normal, residue_correction)
    assert local_exponent_determinant == 3

    matching = cyclic_snc_matching_profile(
        3,
        18,
        left_parameter="W",
        right_parameter="U",
    )
    assert matching.cokernel_model == "R/(W^3*U^18)"
    assert matching.fitting_zero == "(W^3*U^18)"
    assert matching.fitting_one == "R"
    assert matching.matching_quotient == "R/(W^3,U^18)"
    assert matching.matching_length == 54
    assert matching.finite_support_torsion == "zero"

    # Independent monomial-basis count for the finite matching quotient.
    quotient_basis = tuple(
        (w_degree, u_degree)
        for w_degree in range(3)
        for u_degree in range(18)
    )
    assert len(quotient_basis) == 54


def main() -> None:
    support_audit()
    upstream_ladder_audit()
    infinity_node_audit()
    print(
        "PASS: F2 upstream carrier ladder is unimodular; the extraction-root "
        "node has coker R/(W^3*U^18) and branch-matching length 54"
    )


if __name__ == "__main__":
    main()
