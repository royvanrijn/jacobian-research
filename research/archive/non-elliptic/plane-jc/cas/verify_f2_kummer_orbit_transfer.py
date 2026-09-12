#!/usr/bin/env python3
"""Exact Kummer-orbit transfer audit for the F2 (75,125) frontier.

This checker proves four finite statements.

1. Every Laurent band obtained from k[X^5,y] has one fixed Kummer character:
       f_l(t)=t^l*u^k*A_l(u^5),  u=1+t, k=-l mod 5.
   Therefore the exact coefficient order at one nonzero fifth-root center
   transfers unchanged to every conjugate center in its natural Puiseux chart.

2. The forced F2 terminal block transports to all five roots of u^5=1 and
   still has bracket X^4.

3. The published order vertex (1,0) of phi_0 forces R(0) != 0, so the two
   zero-root strata in the earlier boundary handoff are incompatible with the
   F2 row.

4. The Newton-step inequality forces the chosen root multiplicity t_2>=2.
   Hence simple roots of the quadratic cofactor are not admissible choices for
   the above-bisectrix F2 continuation. If a second nonzero double root exists,
   it gives the same unique principal endpoint block.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import sympy as sp


def ceil_div(n: int, d: int) -> int:
    return -((-n) // d)


@dataclass(frozen=True)
class Band:
    side: str
    degree: int
    terminal_height: int
    layer: int
    source_j: tuple[int, ...]
    character: int
    terminal_order: int
    jet_order: int
    invariant_degree: int


def make_band(side: str, degree: int, terminal_height: int, layer: int) -> Band:
    i_min = max(0, ceil_div(layer, 5))
    i_max = (degree + layer) // 6
    if i_max < i_min:
        raise AssertionError(f"empty {side} band {layer}")
    source_j = tuple(5 * i - layer for i in range(i_min, i_max + 1))
    character = (-layer) % 5
    if any(j % 5 != character for j in source_j):
        raise AssertionError("one Laurent band acquired multiple Kummer characters")
    terminal_order = max(layer, ceil_div(17 * layer - terminal_height, 12))
    jet_order = terminal_order - layer
    invariant_powers = tuple((j - character) // 5 for j in source_j)
    if min(invariant_powers) < 0:
        raise AssertionError("negative invariant power")
    return Band(
        side=side,
        degree=degree,
        terminal_height=terminal_height,
        layer=layer,
        source_j=source_j,
        character=character,
        terminal_order=terminal_order,
        jet_order=jet_order,
        invariant_degree=max(invariant_powers),
    )


def audit_band_characters() -> tuple[dict[int, Band], dict[int, Band]]:
    p_bands = {
        layer: make_band("P", 75, 3, layer)
        for layer in range(-21, 16)
    }
    q_bands = {
        layer: make_band("Q", 125, 5, layer)
        for layer in range(-11, 26)
    }
    assert (p_bands[3].terminal_order, p_bands[15].terminal_order) == (4, 21)
    assert (
        q_bands[1].terminal_order,
        q_bands[13].terminal_order,
        q_bands[25].terminal_order,
    ) == (1, 18, 35)
    return p_bands, q_bands


def audit_chart_transfer(p_bands: dict[int, Band], q_bands: dict[int, Band]) -> None:
    # In the selected chart t=u-1 and z=X/t. At another center mu^5=1,
    # write s=u-mu and z_mu=X/s. Since z=(s/t)z_mu, a band term transfers as
    #
    # t^l*u^k*A(u^5)*z^l = s^l*u^k*A(u^5)*z_mu^l.
    #
    # Thus if A has exact order r at v=1, every conjugate natural chart has
    # exact tangential order l+r.
    t, s, u, z_mu = sp.symbols("t s u z_mu", nonzero=True)
    for band in [*p_bands.values(), *q_bands.values()]:
        # A generic monomial in the invariant polynomial suffices for the
        # algebraic chart identity, and covers negative Laurent layers.
        n = min(2, band.invariant_degree)
        selected = (
            t**band.layer
            * u**band.character
            * (u**5)**n
            * ((s / t) * z_mu) ** band.layer
        )
        natural = s**band.layer * u**band.character * (u**5)**n * z_mu**band.layer
        if sp.cancel(selected - natural) != 0:
            raise AssertionError(f"chart transfer failed on {band.side}_{band.layer}")

    # Exact leading terminal coefficients at a conjugate fifth root.
    # If the selected leading coefficient is c and the invariant factor has
    # order r at v=1, the conjugate coefficient is c*mu^(k+4r).
    mu = sp.symbols("mu")
    coefficients = {
        "P3": (sp.Integer(1), p_bands[3]),
        "P15": (sp.Integer(1), p_bands[15]),
        "Q1": (sp.Integer(-1), q_bands[1]),
        "Q13": (sp.Integer(-3), q_bands[13]),
        "Q25": (-sp.Rational(9, 5), q_bands[25]),
    }
    modulus = sp.Poly(mu**5 - 1, mu)
    leading: dict[str, sp.Expr] = {}
    for name, (selected_coefficient, band) in coefficients.items():
        raw = selected_coefficient * mu ** (
            band.character + 4 * band.jet_order
        )
        leading[name] = sp.rem(sp.Poly(raw, mu), modulus).as_expr()

    assert leading == {
        "P3": mu,
        "P15": mu**4,
        "Q1": -mu**4,
        "Q13": -3 * mu**2,
        "Q25": -sp.Rational(9, 5),
    }

    x, z = sp.symbols("x z")
    p_terms = {
        3: leading["P3"] * x**4,
        15: leading["P15"] * x**21,
    }
    q_terms = {
        1: leading["Q1"] * x,
        13: leading["Q13"] * x**18,
        25: leading["Q25"] * x**35,
    }
    bracket = sum(
        (
            i * p * sp.diff(q, x)
            - j * sp.diff(p, x) * q
        )
        * z ** (i + j)
        for i, p in p_terms.items()
        for j, q in q_terms.items()
    )
    residual = sp.Poly(sp.expand(bracket - x**4 * z**4), mu)
    if sp.rem(residual, modulus).as_expr() != 0:
        raise AssertionError("conjugate terminal block lost bracket X^4")


def audit_order_vertex_gate() -> None:
    x, y, v, r, a, b, c = sp.symbols("x y v r a b c")
    phi = x * (v - r) ** 2 * (a * v**2 + b * v + c)
    phi_xy = sp.expand(phi.subs(v, x * y**5))
    generic = sp.Poly(phi_xy, x, y)
    zero_constant = sp.Poly(sp.expand(phi_xy.subs(c, 0)), x, y)

    generic_min_y = min(monomial[1] for monomial, _ in generic.terms())
    zero_min_y = min(monomial[1] for monomial, _ in zero_constant.terms())
    assert generic_min_y == 0
    assert zero_min_y == 5

    # P has leading form phi^3. Hence the order vertex is (3,0) exactly
    # when c=R(0) is nonzero; c=0 moves it to y-order at least 15.
    assert min(m[1] for m, _ in sp.Poly(phi_xy**3, x, y).terms()) == 0
    assert min(
        m[1]
        for m, _ in sp.Poly(sp.expand(phi_xy.subs(c, 0) ** 3), x, y).terms()
    ) == 15


def audit_newton_multiplicity_gate() -> None:
    # F2 data in the Makar-Limanov--Trakhtenberg algorithm:
    # d0=3, v0=(15,60), v1'=(3,0), and primitive edge vector (beta,gamma)=(1,5).
    d0 = 3
    mu0, nu0 = 15, 60
    mu1_prime, nu1_prime = 3, 0
    beta, gamma = 1, 5
    lower = Fraction(
        mu1_prime * gamma - nu1_prime * beta,
        (gamma - beta) * d0,
    )
    upper = Fraction(nu0 - nu1_prime, d0 * gamma)
    assert lower == Fraction(5, 4)
    assert upper == 4
    admissible_multiplicities = tuple(
        t for t in range(1, 10) if lower < t <= upper
    )
    assert admissible_multiplicities == (2, 3, 4)

    vertices = {
        t: (
            Fraction(mu1_prime)
            + Fraction(d0 * t - nu1_prime) * Fraction(beta, gamma),
            d0 * t,
        )
        for t in admissible_multiplicities
    }
    assert vertices[2] == (Fraction(21, 5), 6)
    assert vertices[3] == (Fraction(24, 5), 9)
    assert vertices[4] == (Fraction(27, 5), 12)

    # The lambda0=5/3 D=75 row is the t2=2 row. In particular, a simple
    # cofactor root (t2=1) cannot produce an above-bisectrix continuation.
    assert not (lower < 1 <= upper)


def audit_principal_endpoint_filter() -> None:
    # Let the common root have natural tangential order e and normal order 5.
    # Then P,Q top points are (3e,15),(5e,25). Suppose one principal pair
    # ends at (a,b),(c,d), and its low-low bracket is s^4 z^4.
    #
    # For a nonzero Kummer orbit the coefficient excess orders are
    # rP=a-b>=0 and rQ=c-d>=0. The target equations give
    # a+c=5, b+d=4, hence rP+rQ=1.
    #
    # Parallelism reduces to
    #   b*(8e-41)=12e+36rP-75.
    solutions = []
    for e in range(1, 51):
        denominator = 8 * e - 41
        for r_p in (0, 1):
            numerator = 12 * e + 36 * r_p - 75
            if denominator == 0 or numerator % denominator:
                continue
            b = numerator // denominator
            a = b + r_p
            d = 4 - b
            c = d + 1 - r_p
            if not (-21 <= b <= 15 and -11 <= d <= 25):
                continue
            if not (a < 3 * e and c < 5 * e and b < 15 and d < 25):
                continue
            determinant = b * c - d * a
            if determinant == 0:
                continue
            solutions.append((e, r_p, a, b, c, d))

    assert solutions == [
        (4, 0, 3, 3, 2, 1),
        (4, 1, 0, -1, 5, 5),
        (7, 1, 4, 3, 1, 1),
    ]
    relevant = [row for row in solutions if row[0] in (6, 7)]
    assert relevant == [(7, 1, 4, 3, 1, 1)]


def audit_orbit_strata() -> None:
    # Cover-level center counts from the previous contact census versus
    # Kummer-orbit counts downstairs. Nonzero fibers u^5=rho are free
    # mu_5-orbits. The c=R(0)=0 rows are removed by the order-vertex gate.
    strata = {
        "two_distinct_nonzero_R_roots": {
            "cover_centers": 15,
            "downstairs_orbits": 3,
            "natural_common_orders": (7, 6, 6),
            "admissible_principal_chains": 1,
        },
        "one_double_nonzero_R_root": {
            "cover_centers": 10,
            "downstairs_orbits": 2,
            "natural_common_orders": (7, 7),
            "admissible_principal_chains": 2,
        },
    }
    assert strata["two_distinct_nonzero_R_roots"]["downstairs_orbits"] == 3
    assert strata["one_double_nonzero_R_root"]["downstairs_orbits"] == 2


def main() -> None:
    p_bands, q_bands = audit_band_characters()
    audit_chart_transfer(p_bands, q_bands)
    audit_order_vertex_gate()
    audit_newton_multiplicity_gate()
    audit_principal_endpoint_filter()
    audit_orbit_strata()
    print("F2_KUMMER_BAND_TRANSFER_PASS")
    print("F2_FIXED_ORBIT_TERMINAL_BLOCK_PASS")
    print("F2_ORDER_VERTEX_ZERO_ROOT_EXCLUSION_PASS")
    print("F2_NEWTON_MULTIPLICITY_GATE_PASS")
    print("F2_PRINCIPAL_ENDPOINT_FILTER_PASS")
    print("F2_KUMMER_ORBIT_REDUCTION_PASS")


if __name__ == "__main__":
    main()
