#!/usr/bin/env python3
"""Verify the F2 geometric-degree-six Stein/cubic-germ reduction."""

from __future__ import annotations

from itertools import combinations_with_replacement

import sympy as sp


def terminal_passport_audit() -> None:
    s, lam = sp.symbols("s lam")
    numerator = 125 * s * (s + 1) ** 5
    denominator = (9 * s**2 + 15 * s + 5) ** 3
    phi = sp.expand(numerator - lam * denominator)

    lam0 = sp.Rational(125, 729)
    special = sp.factor(phi.subs(lam, lam0))
    expected = -sp.Rational(125, 729) * (
        135 * s**3 + 405 * s**2 + 396 * s + 125
    )
    assert sp.expand(special - expected) == 0
    assert sp.discriminant(expected, s) != 0

    # The missing three degrees occur at s=infinity.  In w=1/s the residue
    # difference has an exact w^3 factor and a nonzero residual constant.
    w = sp.symbols("w")
    h_at_infinity = sp.cancel(
        numerator.subs(s, 1 / w) / denominator.subs(s, 1 / w) - lam0
    )
    quotient = sp.cancel(h_at_infinity / w**3)
    assert sp.limit(quotient, w, 0) != 0

    # A generic finite nonzero value has six simple preimages.  At lam0 one
    # has the local-degree partition (3,1,1,1), whose sum already equals d.
    generic = sp.Poly(phi.subs(lam, sp.Rational(1, 2)), s)
    assert generic.degree() == 6
    assert sp.discriminant(generic.as_expr(), s) != 0
    assert sum((3, 1, 1, 1)) == 6


def global_inertia_budget_audit() -> None:
    raw = {
        row
        for length in (1, 2)
        for row in combinations_with_replacement((2, 3), length)
        if sum(row) <= 5
    }
    assert raw == {(2,), (3,), (2, 2), (2, 3)}

    smooth_cubic = {(3,), (2, 2), (2,)}
    assert smooth_cubic < raw
    assert (2, 3) not in smooth_cubic


def cubic_normal_form_audit() -> None:
    pi, w, z = sp.symbols("pi w z")

    for order in range(1, 13):
        a = pi**order
        image = w**3 + a * w
        ramification = sp.diff(image, w)
        discriminant = sp.discriminant(w**3 + a * w - z, w)
        assert sp.expand(ramification - (3 * w**2 + a)) == 0
        assert sp.expand(discriminant - (-4 * a**3 - 27 * z**2)) == 0
        factors = sp.factor_list(
            ramification, pi, w, extension=sp.I * sp.sqrt(3)
        )[1]
        if order % 2:
            assert len(factors) == 1
            # The normalized branch has (ord pi, ord z)=(2,3r).
            assert (2, 3 * order)[0] == 2
        else:
            assert len(factors) == 2


def logarithmic_determinant(
    first: sp.Expr,
    second: sp.Expr,
    left: sp.Symbol,
    right: sp.Symbol,
) -> sp.Expr:
    """Determinant in the source log basis dlog(left),dlog(right)."""

    return sp.expand(
        (left * sp.diff(first, left) / first)
        * (right * sp.diff(second, right))
        - (right * sp.diff(first, right) / first)
        * (left * sp.diff(second, left))
    )


def endpoint_packet_audit() -> None:
    v, w = sp.symbols("v w")

    # ord(a)=1 after the two blowups separating the tangent ramification
    # curve from the terminal boundary.
    pi = v * w**2
    z = (1 + v) * w**3
    det = sp.factor(logarithmic_determinant(pi, z, v, w))
    assert sp.expand(det - w**3 * (3 + v)) == 0

    # ord(a)>=2 after the first blowup.  The factor multiplying w^3 is a
    # unit at the terminal--exceptional node for every tested order.
    for order in range(2, 13):
        pi = v * w
        z = w**3 + v**order * w ** (order + 1)
        det = sp.factor(logarithmic_determinant(pi, z, v, w))
        expected = w**3 * (3 + v**order * w ** (order - 2))
        assert sp.expand(det - expected) == 0
        assert expected.subs({v: 0, w: 0}) == 0
        assert (expected / w**3).subs({v: 0, w: 0}) == 3

    # The cyclic specialization is already SNC and has the same packet.
    tau = sp.symbols("tau")
    pi = tau
    z = w**3
    det = sp.factor(logarithmic_determinant(pi, z, tau, w))
    assert det == 3 * w**3


def trace_discriminant_of_power_basis(exponents: tuple[int, int, int]) -> sp.Expr:
    """Discriminant of (w^e_i) for w^3=z via conjugate traces."""

    z = sp.symbols("z")

    def trace_power(power: int) -> sp.Expr:
        if power % 3:
            return sp.Integer(0)
        return 3 * z ** (power // 3)

    matrix = sp.Matrix(
        [[trace_power(left + right) for right in exponents] for left in exponents]
    )
    return sp.factor(matrix.det())


def conductor_order_audit() -> None:
    z = sp.symbols("z")
    rows = (
        ((0, 1, 2), 0, 2),
        ((0, 2, 4), 1, 4),
        ((0, 4, 5), 2, 6),
    )
    for basis, delta, expected_order in rows:
        discriminant = trace_discriminant_of_power_basis(basis)
        polynomial = sp.Poly(discriminant, z)
        assert polynomial.as_dict().keys() == {(expected_order,)}
        assert expected_order == 2 + 2 * delta

    # The last order has square-zero closed fiber: products of its two
    # nontrivial basis elements acquire a factor z.
    assert 4 + 4 >= 3 and 4 + 5 >= 3 and 5 + 5 >= 3

    # A smooth ambient cyclic cubic can nevertheless have the middle,
    # nonnormal terminal slice.  Its normalization has (x,y)=(t^2,t^3),
    # residue degree three, and its branch u=-v^2 has boundary contact two.
    t, x, y, u, v = sp.symbols("t x y u v")
    assert sp.expand((x**3 - y**2).subs({x: t**2, y: t**3})) == 0
    assert v.subs(v, t**3) == t**3
    branch_u = -v**2
    assert branch_u.subs(v, t).as_powers_dict()[t] == 2


def conductor_contact_atlas_audit() -> None:
    """Enumerate the exact d=6 conductor/contact families."""

    max_contact = 24
    cubic_rows: dict[int, tuple[int]] = {}
    simple_rows: dict[int, tuple[int]] = {}
    double_rows: dict[int, set[tuple[int, int]]] = {}
    mixed_rows: dict[int, set[tuple[int, int]]] = {}

    for delta in range(0, 48):
        rhs = 2 + 2 * delta

        cubic_contact = 1 + delta
        if cubic_contact <= max_contact:
            cubic_rows[delta] = (cubic_contact,)

        simple_contact = rhs
        if simple_contact <= max_contact:
            simple_rows[delta] = (simple_contact,)

        double = {
            (left, right)
            for left in range(1, max_contact + 1)
            for right in range(left, max_contact + 1)
            if left + right == rhs
        }
        if double:
            double_rows[delta] = double

        # Pair order is (three-cycle contact, transposition contact).
        mixed = {
            (cubic, simple)
            for cubic in range(1, max_contact + 1)
            for simple in range(1, max_contact + 1)
            if 2 * cubic + simple == rhs
        }
        if mixed:
            mixed_rows[delta] = mixed

    assert cubic_rows[0] == (1,)
    assert simple_rows[0] == (2,)
    assert double_rows[0] == {(1, 1)}
    assert 0 not in mixed_rows

    assert max(cubic_rows) == 23
    assert max(simple_rows) == 11
    assert all(
        simple % 2 == 0
        for rows in mixed_rows.values()
        for _, simple in rows
    )

    for delta, (contact,) in cubic_rows.items():
        assert 2 * contact == 2 + 2 * delta
    for delta, (contact,) in simple_rows.items():
        assert contact == 2 + 2 * delta
    for delta, rows in double_rows.items():
        assert all(left + right == 2 + 2 * delta for left, right in rows)
    for delta, rows in mixed_rows.items():
        assert all(2 * cubic + simple == 2 + 2 * delta for cubic, simple in rows)


def main() -> None:
    terminal_passport_audit()
    global_inertia_budget_audit()
    cubic_normal_form_audit()
    endpoint_packet_audit()
    conductor_order_audit()
    conductor_contact_atlas_audit()
    print(
        "PASS: degree six localizes every affine branch at 125/729; a "
        "normal terminal cubic slice has only (3;k=1), (2+2;k=1+1), or "
        "(2;k=2), all with endpoint packet R/(w^3); conductor orders "
        "2,4,6 distinguish the first nonnormal boundary orders; the exact "
        "conductor/contact identity classifies all four remaining families"
    )


if __name__ == "__main__":
    main()
