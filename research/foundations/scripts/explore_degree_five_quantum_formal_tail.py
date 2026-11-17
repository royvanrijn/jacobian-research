#!/usr/bin/env python3
"""Canonical X-adic hbar^5 primitive at the known degree-five seed.

This complements the Laurent no-go certificate.  For the free-zero
polynomial hbar^3 lift at ``(kappa,tau)=(0,1)``, set ``T4=0`` and solve

    {S4,T} = -O5

recursively in ``Q[Q,Z][[X]]``.  Since

    T|_(X=0) = -3Q/2,

the leading operator is ``{-,T_0}=-3*d_Z``.  At each X-order the zero
Z-integration-constant primitive is therefore unique.

The resulting formal correction starts regularly at X^12; it has no
negative-X principal part.  Its first term outside the filtered S4 space is
the nonzero X^14 Z^2 coefficient displayed below.  Thus the formal-local
primitive fails by differential order (and later by an infinite positive-X
tail), not by a pole at X=0.
"""

from sympy.polys.domains import QQ

from explore_degree_five_quantum_residue import (
    GENERIC_PROFILE,
    add,
    degree_five_family,
    third_order_family,
)


def main() -> None:
    field = QQ
    profile = GENERIC_PROFILE
    S, T = degree_five_family(
        field,
        -field.one / field(2),
        field.one,
    )
    third = third_order_family(S, T, field)
    S2, T2 = third.base

    defect = profile.poisson(S2, T2)
    defect = add(
        defect,
        profile.pi_power(S2, T, 3),
        field.one / field(24),
    )
    defect = add(
        defect,
        profile.pi_power(S, T2, 3),
        field.one / field(24),
    )
    defect = add(
        defect,
        profile.pi_power(S, T, 5),
        field.one / field(1920),
    )
    assert min(x_degree for x_degree, _, _ in defect) == 12

    S4 = {}
    max_x_order = 30
    for x_order in range(max_x_order + 1):
        residual = add(profile.poisson(S4, T), defect)
        coefficient = {
            (0, q_degree, z_degree): value
            for (x_degree, q_degree, z_degree), value in residual.items()
            if x_degree == x_order
        }
        # D_0=-3*d_Z.  Adding u with u_Z=residual/3 kills the coefficient.
        for (_, q_degree, z_degree), value in coefficient.items():
            monomial = (x_order, q_degree, z_degree + 1)
            S4[monomial] = S4.get(monomial, field.zero) + (
                value / (3 * (z_degree + 1))
            )
        residual = add(profile.poisson(S4, T), defect)
        assert not any(
            x_degree <= x_order
            for x_degree, _, _ in residual
        )

    expected_first_terms = {
        (12, 0, 1): field(-17302696248868, 7203),
        (13, 1, 1): field(-137721444722906, 2401),
        (14, 0, 2): field(69014023267862, 2401),
        (14, 2, 1): field(33989945008104, 16807),
    }
    for monomial, expected in expected_first_terms.items():
        assert S4[monomial] == expected

    filtered_violations = [
        (monomial, value)
        for monomial, value in sorted(S4.items())
        if monomial[2] > 1
        or monomial[0] + monomial[1] + 3 * monomial[2] > 21
    ]
    assert filtered_violations[0] == (
        (14, 0, 2),
        field(69014023267862, 2401),
    )

    print("PASS: canonical X-adic hbar^5 primitive constructed through X^30")
    print("PASS: the primitive is regular and starts at X^12*Z")
    print(
        "FIRST FILTERED FAILURE: "
        "coeff(X^14*Z^2)=69014023267862/2401"
    )
    print("SCOPE: this primitive fixes T4=0 and zero Z-integration constants")


if __name__ == "__main__":
    main()
