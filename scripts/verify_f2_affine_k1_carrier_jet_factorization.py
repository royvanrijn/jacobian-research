#!/usr/bin/env python3
"""Verify the exact F2 k=1 carrier-jet and terminal factor patterns."""

from __future__ import annotations

from itertools import product

import sympy as sp

from verify_f2_affine_target_k1_implicit_conductor import (
    expected_implicit_quintic,
)


u, z, C, s, q, t = sp.symbols("u z C s q t")
P, Q = sp.symbols("P Q")
a, b, c, d = sp.symbols("a b c d")
P0, Q0, Gamma = sp.symbols("P0 Q0 Gamma")
h1, h2, h3, h4, h5, h6, h7 = sp.symbols("h1:8")
j1, j2, j3, j4, j5, j6, j7 = sp.symbols("j1:8")


def normalized_forced_relations() -> list[sp.Expr]:
    """Return H5,H6,H7 as polynomials in the first four normalized jets."""

    return [
        (
            10 * h1**3 * h2
            - 2025 * h1**2 * h3
            - 2223 * h1 * h2**2
            - 4050 * h1 * h4
            - 4320 * h2 * h3
        )
        / 2025,
        (
            -100 * h1**6
            + 1575 * h1**4 * h2
            + 118800 * h1**3 * h3
            + 118746 * h1**2 * h2**2
            + 170100 * h1**2 * h4
            + 63180 * h1 * h2 * h3
            - 43011 * h2**3
            - 149445 * h2 * h4
            - 72900 * h3**2
        )
        / 54675,
        (
            4375 * h1**7
            - 78750 * h1**5 * h2
            - 2868750 * h1**4 * h3
            - 2685960 * h1**3 * h2**2
            - 3543750 * h1**3 * h4
            + 2126250 * h1**2 * h2 * h3
            + 3736854 * h1 * h2**3
            + 6561000 * h1 * h2 * h4
            + 455625 * h1 * h3**2
            + 2646270 * h2**2 * h3
            - 2733750 * h3 * h4
        )
        / 820125,
    ]


def coprime_parametrization_degree_audit() -> None:
    """Replay the numerical input to the birationality argument."""

    p = t**3 + a * t
    target_q = t**5 + b * t**4 + c * t**2 + d * t
    assert sp.Poly(p, t).degree() == 3
    assert sp.Poly(target_q, t).degree() == 5
    assert sp.gcd(3, 5) == 1


def normalized_edge_factor_audit() -> None:
    """Record the complete weighted restriction to (P,Q)=(C^3,C^5)."""

    implicit = expected_implicit_quintic()
    edge = sp.expand(implicit.subs({P: C**3, Q: C**5}, simultaneous=True))
    edge_jet = sp.expand(z**15 * edge.subs(C, 1 / z))
    S = a**4 + a**3 * b**2 - 2 * a**2 * b * c + 2 * a**2 * d + a * c**2 + d**2
    expected = {
        1: 3 * b,
        2: -5 * a,
        3: a * b + b**3 + 3 * c,
        4: -5 * a**2 - 4 * a * b**2 + 3 * b * c + 3 * d,
        5: 2 * a * (a * b - c),
        6: -a * b * c + 4 * a * d + 3 * b**2 * d + 3 * c**2,
        7: a**3 * b - 3 * a**2 * c - 5 * a * b * d + 3 * c * d,
        8: 0,
        9: (
            a**3 * c
            + a**2 * b**2 * c
            - a**2 * b * d
            - 2 * a * b * c**2
            + 5 * a * c * d
            + 3 * b * d**2
            + c**3
        ),
        10: -a * S,
        11: 0,
        12: d * S,
    }
    assert sp.expand(edge_jet - sum(expected[j] * z**j for j in expected)) == 0
    assert sp.rem(sp.Poly(edge, C), sp.Poly(C**3, C)) == 0
    assert sp.expand(sp.Poly(edge / C**3, C).coeff_monomial(1) - d * S) == 0

    # The factor C^3 only says that the chosen normal form and the pure
    # common-power cusp both pass through the translated point (0,0).
    assert sp.factor(edge.subs({a: 1, b: 0, c: 0, d: 0})) == (
        -C**5 * (5 * C**8 + 5 * C**6 + 1)
    )


def carrier_graph_jet_audit() -> tuple[list[sp.Expr], list[sp.Expr]]:
    """Verify the normalized seven-jet and its triangular inverse."""

    pi = u * (1 + a * u**2) ** 3 / (1 + b * u + c * u**3 + d * u**4) ** 2
    h = (1 + a * u**2) ** 5 / (1 + b * u + c * u**3 + d * u**4) ** 3
    H = [
        -3 * b,
        5 * a,
        14 * a * b - b**3 - 3 * c,
        -20 * a**2 + 46 * a * b**2 - 3 * b**4 - 12 * b * c - 3 * d,
        -187 * a**2 * b
        + 160 * a * b**3
        + 32 * a * c
        - 9 * b**5
        - 45 * b**2 * c
        - 18 * b * d,
        175 * a**3
        - 1204 * a**2 * b**2
        + 574 * a * b**4
        + 328 * a * b * c
        + 41 * a * d
        - 28 * b**6
        - 168 * b**3 * c
        - 84 * b**2 * d
        - 12 * c**2,
        2754 * a**3 * b
        - 6630 * a**2 * b**3
        - 442 * a**2 * c
        + 2100 * a * b**5
        + 2250 * a * b**2 * c
        + 500 * a * b * d
        - 90 * b**7
        - 630 * b**4 * c
        - 360 * b**3 * d
        - 135 * b * c**2
        - 30 * c * d,
    ]
    graph = 1 + sum(H[j - 1] * pi**j for j in range(1, 8))
    assert sp.series(h - graph, u, 0, 8).removeO().expand() == 0

    inverse = [
        h2 / 5,
        -h1 / 3,
        (5 * h1**3 - 126 * h1 * h2 - 135 * h3) / 405,
        (
            5 * h1**4
            - 90 * h1**2 * h2
            - 540 * h1 * h3
            - 324 * h2**2
            - 405 * h4
        )
        / 1215,
    ]
    inverse_substitution = dict(zip((a, b, c, d), inverse))
    for expected_coordinate, actual_coordinate in zip(
        (h1, h2, h3, h4), H[:4]
    ):
        assert sp.cancel(actual_coordinate.subs(inverse_substitution) - expected_coordinate) == 0
    return H, inverse


def jet_complete_intersection_audit(H: list[sp.Expr], inverse: list[sp.Expr]) -> None:
    """Verify the three forced jets and their weighted homogeneity."""

    forced = normalized_forced_relations()
    inverse_substitution = dict(zip((a, b, c, d), inverse))
    for actual, expected in zip(H[4:], forced):
        assert sp.cancel(actual.subs(inverse_substitution) - expected) == 0

    variables = (h1, h2, h3, h4)
    weights = (1, 2, 3, 4)
    for degree, polynomial in zip((5, 6, 7), forced):
        numerator = sp.cancel(polynomial).as_numer_denom()[0]
        for monomial, coefficient in sp.Poly(numerator, *variables).terms():
            if coefficient:
                assert sum(power * weight for power, weight in zip(monomial, weights)) == degree

    # Each relation is linear in a new variable with a nonzero constant
    # coefficient.  Consequently the quotient is Q[h1,h2,h3,h4], proving
    # that the seven-jet locus is a prime codimension-three complete
    # intersection, not merely a set-theoretic parametrized image.
    for variable, polynomial in zip((h5, h6, h7), forced):
        relation = sp.together(variable - polynomial)
        relation_polynomial = sp.Poly(relation.as_numer_denom()[0], variable)
        assert relation_polynomial.degree() == 1
        assert relation_polynomial.LC() != 0


def fixed_coordinate_transport_audit() -> tuple[
    list[sp.Expr], list[sp.Expr], list[sp.Expr]
]:
    """Compile the fixed-target centers into normalized carrier jets."""

    x = sp.Symbol("x")
    mu, eta, nu = sp.symbols("mu eta nu")
    normalized_jets = (h1, h2, h3, h4, h5, h6, h7)
    fixed_jets = (j1, j2, j3, j4, j5, j6, j7)
    normalized_graph = 1 + sum(
        coefficient * x**index
        for index, coefficient in enumerate(normalized_jets, start=1)
    )

    numerator_unit = 1 + mu * x**3 / normalized_graph**2
    denominator_unit = (
        1
        + eta * x**2 / normalized_graph
        + nu * x**5 / normalized_graph**3
    )
    fixed_x = x * numerator_unit**3 / denominator_unit**2
    fixed_y = normalized_graph * numerator_unit**5 / denominator_unit**3

    forward = [
        h1,
        h2 - 3 * eta,
        2 * h1 * eta + h3 + 5 * mu,
        -2 * h1**2 * eta - 8 * h1 * mu + 4 * h2 * eta + h4 - 6 * eta**2,
        2 * h1**3 * eta
        + 11 * h1**2 * mu
        - 6 * h1 * h2 * eta
        + 15 * h1 * eta**2
        - 11 * h2 * mu
        + 6 * h3 * eta
        + h5
        + 33 * eta * mu
        - 3 * nu,
        -2 * h1**4 * eta
        - 14 * h1**3 * mu
        + 8 * h1**2 * h2 * eta
        - 28 * h1**2 * eta**2
        + 28 * h1 * h2 * mu
        - 8 * h1 * h3 * eta
        - 112 * h1 * eta * mu
        + 8 * h1 * nu
        - 4 * h2**2 * eta
        + 28 * h2 * eta**2
        - 14 * h3 * mu
        + 8 * h4 * eta
        + h6
        - 28 * eta**3
        - 35 * mu**2,
        2 * h1**5 * eta
        + 17 * h1**4 * mu
        - 10 * h1**3 * h2 * eta
        + 45 * h1**3 * eta**2
        - 51 * h1**2 * h2 * mu
        + 10 * h1**2 * h3 * eta
        + 255 * h1**2 * eta * mu
        - 15 * h1**2 * nu
        + 10 * h1 * h2**2 * eta
        - 90 * h1 * h2 * eta**2
        + 34 * h1 * h3 * mu
        - 10 * h1 * h4 * eta
        + 120 * h1 * eta**3
        + 153 * h1 * mu**2
        + 17 * h2**2 * mu
        - 10 * h2 * h3 * eta
        - 170 * h2 * eta * mu
        + 10 * h2 * nu
        + 45 * h3 * eta**2
        - 17 * h4 * mu
        + 10 * h5 * eta
        + h7
        + 255 * eta**2 * mu
        - 30 * eta * nu,
    ]
    fixed_candidate = 1 + sum(
        coefficient * fixed_x**index
        for index, coefficient in enumerate(forward, start=1)
    )
    assert sp.series(fixed_y - fixed_candidate, x, 0, 8).removeO().expand() == 0

    inverse = [
        j1,
        j2 + 3 * eta,
        -2 * j1 * eta + j3 - 5 * mu,
        2 * j1**2 * eta + 8 * j1 * mu - 4 * j2 * eta + j4 - 6 * eta**2,
        -2 * j1**3 * eta
        - 11 * j1**2 * mu
        + 6 * j1 * j2 * eta
        + 15 * j1 * eta**2
        + 11 * j2 * mu
        - 6 * j3 * eta
        + j5
        + 30 * eta * mu
        + 3 * nu,
        2 * j1**4 * eta
        + 14 * j1**3 * mu
        - 8 * j1**2 * j2 * eta
        - 28 * j1**2 * eta**2
        - 28 * j1 * j2 * mu
        + 8 * j1 * j3 * eta
        - 104 * j1 * eta * mu
        - 8 * j1 * nu
        + 4 * j2**2 * eta
        + 28 * j2 * eta**2
        + 14 * j3 * mu
        - 8 * j4 * eta
        + j6
        + 28 * eta**3
        - 35 * mu**2,
        -2 * j1**5 * eta
        - 17 * j1**4 * mu
        + 10 * j1**3 * j2 * eta
        + 45 * j1**3 * eta**2
        + 51 * j1**2 * j2 * mu
        - 10 * j1**2 * j3 * eta
        + 240 * j1**2 * eta * mu
        + 15 * j1**2 * nu
        - 10 * j1 * j2**2 * eta
        - 90 * j1 * j2 * eta**2
        - 34 * j1 * j3 * mu
        + 10 * j1 * j4 * eta
        - 120 * j1 * eta**3
        + 153 * j1 * mu**2
        - 17 * j2**2 * mu
        + 10 * j2 * j3 * eta
        - 160 * j2 * eta * mu
        - 10 * j2 * nu
        + 45 * j3 * eta**2
        + 17 * j4 * mu
        - 10 * j5 * eta
        + j7
        - 225 * eta**2 * mu
        - 30 * eta * nu,
    ]
    forward_substitution = dict(zip(fixed_jets, forward))
    for expected, actual in zip(normalized_jets, inverse):
        assert sp.expand(actual.subs(forward_substitution) - expected) == 0

    # The fixed compatibility equations are obtained without expansion:
    # untransport J to H and apply the three normalized relations.  This is
    # triangular, linear in J5,J6,J7, and preserves the extended weights
    # wt(J_i)=i, wt(eta,mu,nu)=(2,3,5).
    first_four_substitution = dict(zip((h1, h2, h3, h4), inverse[:4]))
    fixed_relations = [
        sp.together(inverse[index] - relation.subs(first_four_substitution))
        for index, relation in enumerate(normalized_forced_relations(), start=4)
    ]
    weighted_variables = fixed_jets + (eta, mu, nu)
    weights = (1, 2, 3, 4, 5, 6, 7, 2, 3, 5)
    for degree, variable, relation in zip((5, 6, 7), fixed_jets[4:], fixed_relations):
        numerator = relation.as_numer_denom()[0]
        polynomial = sp.Poly(numerator, *weighted_variables)
        assert polynomial.degree(variable) == 1
        assert polynomial.LC() != 0
        for monomial, coefficient in polynomial.terms():
            if coefficient:
                assert sum(power * weight for power, weight in zip(monomial, weights)) == degree
    return forward, inverse, fixed_relations


def fixed_jet_dominance_audit(
    normalized_curve_jets: list[sp.Expr], fixed_transport: list[sp.Expr]
) -> None:
    """Prove that free fixed-target translations absorb the three residuals."""

    mu, eta, nu = sp.symbols("mu eta nu")
    normalized_substitution = dict(
        zip((h1, h2, h3, h4, h5, h6, h7), normalized_curve_jets)
    )
    fixed_curve_jets = [
        sp.expand(coefficient.subs(normalized_substitution))
        for coefficient in fixed_transport
    ]
    parameters = (a, b, c, d, mu, eta, nu)
    jet_jacobian = sp.Matrix(fixed_curve_jets).jacobian(parameters)

    immersion_resultant = sp.factor(
        sp.resultant(
            sp.diff(t**3 + a * t, t),
            sp.diff(t**5 + b * t**4 + c * t**2 + d * t, t),
            t,
        )
    )
    expected_resultant = (
        25 * a**4
        + 48 * a**3 * b**2
        - 144 * a**2 * b * c
        + 90 * a**2 * d
        + 108 * a * c**2
        + 81 * d**2
    )
    assert sp.expand(immersion_resultant - expected_resultant) == 0
    assert sp.expand(jet_jacobian.det() - 3 * immersion_resultant) == 0

    witness = {a: 1, b: 0, c: 0, d: 0, mu: 0, eta: 0, nu: 0}
    assert jet_jacobian.subs(witness).rank() == 7


def e6_carrier_substratum_audit(
    normalized_curve_jets: list[sp.Expr], fixed_transport: list[sp.Expr]
) -> None:
    """Identify the exact carrier-jet curve of the E6+A1 escape stratum."""

    beta = sp.Symbol("beta")
    e6_substitution = {a: 0, b: beta, c: 0, d: 0}
    e6_jets = [
        sp.expand(coefficient.subs(e6_substitution))
        for coefficient in normalized_curve_jets
    ]
    closed_form = [
        sp.simplify(
            -sp.Rational(3, index)
            * sp.binomial(2 * index - 4, index - 1)
            * beta**index
        )
        for index in range(1, 8)
    ]
    assert e6_jets == closed_form
    assert e6_jets == [
        -3 * beta,
        0,
        -beta**3,
        -3 * beta**4,
        -9 * beta**5,
        -28 * beta**6,
        -90 * beta**7,
    ]

    # The closed form follows from z=beta*u and
    # beta*pi=z/(1+z)^2, h=(1+z)^(-3).  The first seven coefficients
    # therefore cut out this explicit normalized monomial jet curve.
    normalized_relations = (
        h2,
        27 * h3 - h1**3,
        27 * h4 + h1**4,
        27 * h5 - h1**5,
        729 * h6 + 28 * h1**6,
        243 * h7 - 10 * h1**7,
    )
    jet_substitution = dict(zip((h1, h2, h3, h4, h5, h6, h7), e6_jets))
    assert all(
        sp.expand(relation.subs(jet_substitution)) == 0
        for relation in normalized_relations
    )

    # After the three target transports, the raw fixed E6 jet locus has
    # four parameters.  Its Jacobian attains rank four, so it has local
    # codimension three in the seven-center carrier space.
    mu, eta, nu = sp.symbols("mu eta nu")
    normalized_substitution = dict(
        zip((h1, h2, h3, h4, h5, h6, h7), e6_jets)
    )
    fixed_e6_jets = [
        sp.expand(coefficient.subs(normalized_substitution))
        for coefficient in fixed_transport
    ]
    e6_jacobian = sp.Matrix(fixed_e6_jets).jacobian((beta, mu, eta, nu))
    assert e6_jacobian.subs({beta: 1, mu: 0, eta: 0, nu: 0}).rank() == 4

    # Eliminate all four E6 and target-transport parameters triangularly.
    # The remaining fixed-coordinate locus is a prime codimension-three
    # complete intersection in J1,...,J7.
    recovery = {
        beta: -j1 / 3,
        eta: -j2 / 3,
        mu: -(j1**3 - 18 * j1 * j2 - 27 * j3) / 135,
        nu: -(
            6 * j1**5
            - 109 * j1**3 * j2
            - 297 * j1**2 * j3
            - 27 * j1 * j2**2
            + 297 * j2 * j3
            + 135 * j5
        )
        / 405,
    }
    for index in (0, 1, 2, 4):
        assert sp.expand(
            fixed_e6_jets[index].subs(recovery)
            - (j1, j2, j3, j4, j5, j6, j7)[index]
        ) == 0

    fixed_relations = (
        j1**4 - 18 * j1**2 * j2 - 72 * j1 * j3 - 30 * j2**2 - 45 * j4,
        -187 * j1**6
        + 3186 * j1**4 * j2
        + 11178 * j1**3 * j3
        + 6480 * j1**2 * j2**2
        - 972 * j1 * j2 * j3
        - 9720 * j1 * j5
        + 3780 * j2**3
        - 5103 * j3**2
        - 3645 * j6,
        89 * j1**7
        - 1604 * j1**5 * j2
        - 5181 * j1**4 * j3
        - 1314 * j1**3 * j2**2
        + 6408 * j1**2 * j2 * j3
        + 3375 * j1**2 * j5
        + 4131 * j1 * j3**2
        - 1125 * j2**2 * j3
        - 2250 * j2 * j5
        - 675 * j7,
    )
    for index, denominator, relation in zip(
        (3, 5, 6), (45, 3645, 675), fixed_relations
    ):
        assert sp.expand(
            fixed_e6_jets[index].subs(recovery)
            - (j1, j2, j3, j4, j5, j6, j7)[index]
            - relation / denominator
        ) == 0
    for variable, relation in zip((j4, j6, j7), fixed_relations):
        polynomial = sp.Poly(relation, variable)
        assert polynomial.degree() == 1
        assert polynomial.LC() != 0

    # Undo J_j=kappa^j*zeta_j/lambda.  Weighted homogeneity cancels kappa,
    # leaving three scale-free equations in the actual target-shear centers
    # and the carrier residue lambda.
    kappa, carrier_residue = sp.symbols("kappa carrier_residue", nonzero=True)
    zeta = sp.symbols("zeta1:8")
    raw_substitution = {
        variable: kappa**index * zeta[index - 1] / carrier_residue
        for index, variable in enumerate((j1, j2, j3, j4, j5, j6, j7), start=1)
    }
    raw_relations = (
        zeta[0] ** 4
        - 18 * carrier_residue * zeta[0] ** 2 * zeta[1]
        - 72 * carrier_residue**2 * zeta[0] * zeta[2]
        - 30 * carrier_residue**2 * zeta[1] ** 2
        - 45 * carrier_residue**3 * zeta[3],
        -187 * zeta[0] ** 6
        + 3186 * carrier_residue * zeta[0] ** 4 * zeta[1]
        + 11178 * carrier_residue**2 * zeta[0] ** 3 * zeta[2]
        + 6480 * carrier_residue**2 * zeta[0] ** 2 * zeta[1] ** 2
        - 972 * carrier_residue**3 * zeta[0] * zeta[1] * zeta[2]
        - 9720 * carrier_residue**4 * zeta[0] * zeta[4]
        + 3780 * carrier_residue**3 * zeta[1] ** 3
        - 5103 * carrier_residue**4 * zeta[2] ** 2
        - 3645 * carrier_residue**5 * zeta[5],
        89 * zeta[0] ** 7
        - 1604 * carrier_residue * zeta[0] ** 5 * zeta[1]
        - 5181 * carrier_residue**2 * zeta[0] ** 4 * zeta[2]
        - 1314 * carrier_residue**2 * zeta[0] ** 3 * zeta[1] ** 2
        + 6408 * carrier_residue**3 * zeta[0] ** 2 * zeta[1] * zeta[2]
        + 3375 * carrier_residue**4 * zeta[0] ** 2 * zeta[4]
        + 4131 * carrier_residue**4 * zeta[0] * zeta[2] ** 2
        - 1125 * carrier_residue**4 * zeta[1] ** 2 * zeta[2]
        - 2250 * carrier_residue**5 * zeta[1] * zeta[4]
        - 675 * carrier_residue**6 * zeta[6],
    )
    for weight, fixed_relation, raw_relation in zip(
        (4, 6, 7), fixed_relations, raw_relations
    ):
        assert sp.expand(
            fixed_relation.subs(raw_substitution)
            * carrier_residue**weight
            / kappa**weight
            - raw_relation
        ) == 0


def e8_carrier_endpoint_audit(
    normalized_curve_jets: list[sp.Expr], fixed_transport: list[sp.Expr]
) -> None:
    """Specialize the E6 carrier gate to the (3,5) E8 endpoint."""

    # At beta=0 the E6+A1 parametrization coalesces to
    # p=t^3, q=t^5.  Its normalized carrier graph is h=1, so all seven
    # normalized jets vanish.  Only the three target transports remain.
    e8_substitution = {a: 0, b: 0, c: 0, d: 0}
    e8_jets = [
        sp.expand(coefficient.subs(e8_substitution))
        for coefficient in normalized_curve_jets
    ]
    assert e8_jets == [0] * 7

    mu, eta, nu = sp.symbols("mu eta nu")
    normalized_substitution = dict(
        zip((h1, h2, h3, h4, h5, h6, h7), e8_jets)
    )
    fixed_e8_jets = [
        sp.expand(coefficient.subs(normalized_substitution))
        for coefficient in fixed_transport
    ]
    e8_jacobian = sp.Matrix(fixed_e8_jets).jacobian((mu, eta, nu))
    assert e8_jacobian.subs({mu: 0, eta: 0, nu: 0}).rank() == 3

    # The transport parameters are recovered from J2,J3,J5.  The other
    # four coordinates give a prime codimension-four complete intersection
    # that is triangular in J1,J4,J6,J7.
    recovery = {
        eta: -j2 / 3,
        mu: j3 / 5,
        nu: -(11 * j2 * j3 + 5 * j5) / 15,
    }
    for index in (1, 2, 4):
        assert sp.expand(
            fixed_e8_jets[index].subs(recovery)
            - (j1, j2, j3, j4, j5, j6, j7)[index]
        ) == 0

    fixed_relations = (
        j1,
        2 * j2**2 + 3 * j4,
        140 * j2**3 - 189 * j3**2 - 135 * j6,
        5 * j2**2 * j3 + 10 * j2 * j5 + 3 * j7,
    )
    for index, denominator, relation in zip(
        (0, 3, 5, 6), (-1, -3, 135, -3), fixed_relations
    ):
        assert sp.expand(
            fixed_e8_jets[index].subs(recovery)
            - (j1, j2, j3, j4, j5, j6, j7)[index]
            - relation / denominator
        ) == 0
    for variable, relation in zip((j1, j4, j6, j7), fixed_relations):
        polynomial = sp.Poly(relation, variable)
        assert polynomial.degree() == 1
        assert polynomial.LC() != 0

    # Eliminate the leading carrier scale from J_j=kappa^j*zeta_j/lambda.
    # The four raw equations can therefore be inserted directly into the
    # fixed-coordinate Laurent compiler.
    kappa, carrier_residue = sp.symbols("kappa carrier_residue", nonzero=True)
    zeta = sp.symbols("zeta1:8")
    raw_substitution = {
        variable: kappa**index * zeta[index - 1] / carrier_residue
        for index, variable in enumerate((j1, j2, j3, j4, j5, j6, j7), start=1)
    }
    raw_relations = (
        zeta[0],
        2 * zeta[1] ** 2 + 3 * carrier_residue * zeta[3],
        140 * zeta[1] ** 3
        - 189 * carrier_residue * zeta[2] ** 2
        - 135 * carrier_residue**2 * zeta[5],
        5 * zeta[1] ** 2 * zeta[2]
        + 10 * carrier_residue * zeta[1] * zeta[4]
        + 3 * carrier_residue**2 * zeta[6],
    )
    for weight, residue_power, fixed_relation, raw_relation in zip(
        (1, 4, 6, 7), (1, 2, 3, 3), fixed_relations, raw_relations
    ):
        assert sp.expand(
            fixed_relation.subs(raw_substitution)
            * carrier_residue**residue_power
            / kappa**weight
            - raw_relation
        ) == 0


def terminal_and_translation_factor_audit() -> None:
    """Separate coordinate factors from the invariant terminal leading fiber."""

    implicit = expected_implicit_quintic()
    terminal_p = q**-3 * s**2 * (1 + s)
    terminal_normalized_q = (
        sp.Rational(5, 9)
        * q**-5
        * s**3
        * (1 + 3 * s + sp.Rational(9, 5) * s**2)
    )
    normalized_pullback = sp.expand(
        q**15
        * implicit.subs(
            {P: terminal_p, Q: terminal_normalized_q}, simultaneous=True
        )
    )
    terminal_cubic = 135 * s**3 + 405 * s**2 + 396 * s + 125
    assert sp.factor(terminal_cubic) == terminal_cubic
    assert sp.discriminant(terminal_cubic, s) == -98_415
    assert sp.factor(normalized_pullback.subs(q, 0)) == -s**9 * terminal_cubic / 729

    # In arbitrary fixed leading scales the q=0 term is
    # B^3*P_top^5-A^5*Q_top^3.  The special residue relation
    # A^5/(-B)^3=125/729 leaves the same cubic up to the nonzero scalar B^3.
    A5, B3 = sp.symbols("A5 B3")
    terminal_unit = 1 + 3 * s + sp.Rational(9, 5) * s**2
    general_leading = s**9 * (
        B3 * s * (1 + s) ** 5 + A5 * terminal_unit**3
    )
    special_leading = sp.expand(
        general_leading.subs(A5, -sp.Rational(125, 729) * B3)
    )
    assert sp.factor(special_leading) == -B3 * s**9 * terminal_cubic / 729
    assert min(
        monomial[0]
        for monomial, coefficient in sp.Poly(normalized_pullback, s, q).terms()
        if coefficient != 0
    ) == 2

    # Restore fixed F2 target coordinates with A=1 and B=-9/5.  Translation
    # terms give a nonzero constant on the common-power edge and a nonzero
    # s=0 term at the terminal chart, so neither C^3 nor s^2 is invariant.
    B = sp.Rational(-9, 5)
    normalized_P = P - P0
    normalized_Q = (Q - Q0 - Gamma * (P - P0)) / B
    fixed_implicit = sp.expand(
        B**3
        * implicit.subs(
            {P: normalized_P, Q: normalized_Q}, simultaneous=True
        )
    )
    fixed_at_origin = sp.expand(fixed_implicit.subs({P: 0, Q: 0}))
    witness = {P0: 1, Q0: 0, Gamma: 0, a: 0, b: 0, c: 0, d: 0}
    assert fixed_at_origin.subs(witness) == sp.Rational(729, 125)

    terminal_fixed_Q = (
        -q**-5 * s**3 * (1 + 3 * s + sp.Rational(9, 5) * s**2)
    )
    fixed_pullback = sp.expand(
        q**15
        * fixed_implicit.subs(
            {P: terminal_p, Q: terminal_fixed_Q}, simultaneous=True
        )
    )
    assert sp.expand(fixed_pullback.subs(s, 0) - q**15 * fixed_at_origin) == 0
    assert sp.factor(fixed_pullback.subs(q, 0)) == s**9 * terminal_cubic / 125


def bounded_irreducibility_regression() -> None:
    """Experimental scan only: this is not a universal irreducibility proof."""

    implicit = expected_implicit_quintic()
    for values in product(range(-2, 3), repeat=4):
        specialized = sp.Poly(
            implicit.subs(dict(zip((a, b, c, d), values))), P, Q, domain=sp.QQ
        )
        coefficient, factors = sp.factor_list(specialized)
        assert coefficient != 0
        assert len(factors) == 1
        factor, multiplicity = factors[0]
        assert multiplicity == 1
        assert factor.total_degree() == 5


def main() -> None:
    coprime_parametrization_degree_audit()
    normalized_edge_factor_audit()
    H, inverse = carrier_graph_jet_audit()
    jet_complete_intersection_audit(H, inverse)
    fixed_transport, _, _ = fixed_coordinate_transport_audit()
    fixed_jet_dominance_audit(H, fixed_transport)
    e6_carrier_substratum_audit(H, fixed_transport)
    e8_carrier_endpoint_audit(H, fixed_transport)
    terminal_and_translation_factor_audit()
    bounded_irreducibility_regression()
    print(
        "PASS: the normalized k=1 carrier seven-jet is a triangular "
        "codimension-three complete intersection with an exact weighted "
        "fixed-coordinate transport; its raw seven-jet Jacobian is three "
        "times the immersion resultant, so free target translations absorb "
        "the three normalized residuals; the E6+A1 escape is an explicit "
        "one-parameter normalized jet curve and a prime codimension-three "
        "raw fixed-jet complete intersection with three explicit equations; "
        "its E8 endpoint is a prime codimension-four fixed-jet complete "
        "intersection with four scale-free equations; "
        "target translations also remove the "
        "apparent C^3 and s^2 factors, while the terminal leading cubic is "
        "invariant (625 normal-form targets pass the bounded irreducibility "
        "regression)"
    )


if __name__ == "__main__":
    main()
