#!/usr/bin/env python3
"""Verify the side-packet compactifications and their vertex matching.

WQRS1 computes the complete generic base-point clusters on the three open
sides of the quartic base triangle, and WQCA1 computes the vertex-corner
atlas.  This checker supplies the first global gluing layer between them:

* it extracts the three nontrivial algebraic direction packets;
* it compactifies each packet over the side parameter P1;
* it computes the endpoint conductors and packet Hilbert corrections;
* it matches the L1 packets to the V12 P18 branches;
* it identifies the L1 and L3 quadratic residue fields at V13; and
* it matches the L3 first direction and cubic contact to q23 and h23 at V23.

The calculation concerns the normalized cluster-center curves.  It does not
yet assemble every rational center surface into the global normalized Rees
algebra or compute the latter's bigraded Hilbert polynomial.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from verify_quartic_biprojective_graph import (  # noqa: E402
    X0,
    X1,
    X2,
    X3,
    primitive_integral,
    quartic_projective_coordinates,
    run_singular,
    singular_polynomial,
)
from verify_quartic_rees_stratification import (  # noqa: E402
    compact_newton_facets,
    linear_quadratic_transform,
    surface_data,
    weight_initial,
    weight_order,
)


u, v, t, z = sp.symbols("u v t z")
a, b, c = sp.symbols("a b c")
d, e = sp.symbols("d e")
T0, T1, U, V = sp.symbols("T0 T1 U V")
p, q = sp.symbols("p q")


@dataclass(frozen=True)
class Packet:
    name: str
    polynomial: sp.Expr
    exponent: int
    endpoint_delta: int
    factor_degrees: tuple[int, ...]

    @property
    def degree(self) -> int:
        return sp.Poly(self.polynomial, z, domain=sp.QQ).degree()


def assert_associate(actual: sp.Expr, expected: sp.Expr) -> None:
    quotient = sp.cancel(actual / expected)
    assert quotient != 0 and not quotient.free_symbols, (actual, expected, quotient)


def extract_packets(
    coordinates: tuple[sp.Expr, ...],
) -> tuple[Packet, Packet, Packet]:
    line_ideals = {
        "L1": tuple(
            sp.expand(item.subs({X0: u, X1: v, X2: t, X3: 1}))
            for item in coordinates
        ),
        "L3": tuple(
            sp.expand(item.subs({X0: u, X1: t, X2: 1, X3: v}))
            for item in coordinates
        ),
    }

    # L1: the multiplicity-four child has one rational direction and one
    # irreducible cubic.  In z=t*v/u their product is constant in t.
    l1 = line_ideals["L1"]
    root_multiplicity, _, _, _ = surface_data(l1)
    assert root_multiplicity == 6
    first = linear_quadratic_transform(l1, 6, v)
    first_multiplicity, _, _, first_factors = surface_data(first)
    assert first_multiplicity == 4
    rational = next(
        factor.as_expr()
        for factor, _ in first_factors
        if factor.total_degree() == 1
    )
    cubic = next(
        factor.as_expr()
        for factor, _ in first_factors
        if factor.total_degree() == 3
    )
    rational_residual = sp.expand(rational.subs({u: 1, v: z / t}))
    cubic_residual = sp.expand(cubic.subs({u: 1, v: z / t}))
    assert rational_residual == z + 1
    assert cubic_residual == 25 * z**3 + 45 * z**2 + 54 * z + 46
    first_packet = Packet(
        "L1 first packet",
        sp.expand(rational_residual * cubic_residual),
        1,
        6,
        (1, 3),
    )

    # L1: the terminal quadratic is constant in z=t^3*v/u.
    second = linear_quadratic_transform(l1, 6, 5 * t * u - 3 * v)
    second_multiplicity, _, _, second_factors = surface_data(second)
    assert second_multiplicity == 2 and len(second_factors) == 1
    second_child = linear_quadratic_transform(
        second,
        2,
        second_factors[0][0].as_expr(),
    )
    _, _, _, terminal_factors = surface_data(second_child)
    assert len(terminal_factors) == 1
    terminal_l1 = terminal_factors[0][0].as_expr()
    l1_quadratic = sp.expand(terminal_l1.subs({u: 1, v: z / t**3}))
    assert l1_quadratic == -15625 * z**2 + 15750 * z + 756
    second_packet = Packet(
        "L1 second packet",
        l1_quadratic,
        3,
        3,
        (2,),
    )

    # L3: four rational transforms leave one quadratic, constant in the
    # same z=t^3*v/u coordinate.
    current = line_ideals["L3"]
    for direction in (-3 * t * v + 5 * u, v, t**2 * v + u, v):
        multiplicity, _, _, factors = surface_data(current)
        assert multiplicity == 2 and len(factors) == 1
        current = linear_quadratic_transform(current, multiplicity, direction)
    _, _, _, l3_factors = surface_data(current)
    assert len(l3_factors) == 1
    terminal_l3 = l3_factors[0][0].as_expr()
    l3_quadratic = sp.expand(terminal_l3.subs({u: 1, v: z / t**3}))
    assert l3_quadratic == 3 * z**2 - 12 * z + 5
    third_packet = Packet(
        "L3 terminal packet",
        l3_quadratic,
        3,
        3,
        (2,),
    )
    return first_packet, second_packet, third_packet


def packet_equation(packet: Packet) -> sp.Expr:
    polynomial = sp.Poly(packet.polynomial, z, domain=sp.QQ)
    degree = packet.degree
    exponent = packet.exponent
    return sp.expand(
        sum(
            coefficient
            * T0 ** (exponent * monomial[0])
            * T1 ** (exponent * (degree - monomial[0]))
            * V ** monomial[0]
            * U ** (degree - monomial[0])
            for monomial, coefficient in polynomial.terms()
        )
    )


def bihomogeneous_degree(expression: sp.Expr) -> tuple[int, int]:
    polynomial = sp.Poly(expression, T0, T1, U, V, domain=sp.QQ)
    degrees = {
        (monomial[0] + monomial[1], monomial[2] + monomial[3])
        for monomial, _ in polynomial.terms()
    }
    assert len(degrees) == 1
    return next(iter(degrees))


def endpoint_equations(packet: Packet) -> tuple[sp.Expr, sp.Expr]:
    polynomial = packet.polynomial
    degree = packet.degree
    exponent = packet.exponent
    zero = sp.expand(q**degree * polynomial.subs(z, p**exponent / q))
    infinity = sp.expand(
        p ** (exponent * degree) * polynomial.subs(z, q / p**exponent)
    )
    return (
        primitive_integral(zero, (p, q)),
        primitive_integral(infinity, (p, q)),
    )


def verify_packet_geometry(packets: tuple[Packet, ...]) -> None:
    for packet in packets:
        factor_degrees = tuple(
            factor.degree()
            for factor, power in sp.factor_list(
                sp.Poly(packet.polynomial, z, domain=sp.QQ)
            )[1]
            for _ in range(power)
        )
        assert factor_degrees == packet.factor_degrees
        assert sp.discriminant(packet.polynomial, z) != 0
        assert packet.endpoint_delta == math.comb(packet.degree, 2) * packet.exponent

        equation = packet_equation(packet)
        assert bihomogeneous_degree(equation) == (
            packet.degree * packet.exponent,
            packet.degree,
        )

        # All normalized branches meet [U:V]=[0:1] over T0=0 and
        # [U:V]=[1:0] over T1=0, with pairwise contact `exponent`.
        at_zero = sp.factor(equation.subs(T0, 0))
        at_infinity = sp.factor(equation.subs(T1, 0))
        packet_polynomial = sp.Poly(packet.polynomial, z, domain=sp.QQ)
        assert at_zero == (
            packet_polynomial.coeff_monomial(1)
            * T1 ** (packet.exponent * packet.degree)
            * U**packet.degree
        )
        assert at_infinity == (
            packet_polynomial.coeff_monomial(z**packet.degree)
            * T0 ** (packet.exponent * packet.degree)
            * V**packet.degree
        )

        base_degree = packet.degree * packet.exponent
        direction_degree = packet.degree
        image_hilbert = (
            direction_degree,
            base_degree,
            base_degree
            + direction_degree
            - base_degree * direction_degree,
        )
        normalized_hilbert = (
            packet.degree,
            packet.degree * packet.exponent,
            packet.degree,
        )
        assert image_hilbert[:2] == normalized_hilbert[:2]
        assert (
            normalized_hilbert[2] - image_hilbert[2]
            == 2 * packet.endpoint_delta
        )


def verify_packet_conductors(packets: tuple[Packet, ...]) -> None:
    program_parts = [
        'LIB "normal.lib";',
        "proc assertReductionZero(poly f,ideal G,string label)",
        "{ if(reduce(f,std(G))!=0){ ERROR(label); } }",
        "proc assertIdealEqual(ideal A,ideal B,string label)",
        "{ int i; for(i=1;i<=size(A);i++){assertReductionZero(A[i],B,label);}",
        "  for(i=1;i<=size(B);i++){assertReductionZero(B[i],A,label);} }",
    ]
    for index, packet in enumerate(packets):
        zero, infinity = endpoint_equations(packet)
        if packet.degree == 4:
            expected = "ideal(p3,p2q,pq2,q3)"
        else:
            expected = "ideal(p3,q)"
        for suffix, expression in (("zero", zero), ("infinity", infinity)):
            ring_name = f"packet{index}{suffix}Ring"
            polynomial_name = f"packet{index}{suffix}"
            conductor_name = f"packet{index}{suffix}Conductor"
            program_parts.extend(
                [
                    f"ring {ring_name}=0,(p,q),dp;",
                    f"poly {polynomial_name}={singular_polynomial(expression)};",
                    f"ideal {conductor_name}=std(normalConductor(ideal({polynomial_name})));",
                    f'assertIdealEqual({conductor_name},{expected},"{packet.name} {suffix} conductor");',
                    f"if(vdim(std(ideal({polynomial_name})+{conductor_name}))"
                    f"!={packet.endpoint_delta})"
                    f'{{ ERROR("{packet.name} {suffix} delta"); }}',
                ]
            )
    program_parts.append('"PASS quartic side-packet conductors";')
    output = run_singular("\n".join(program_parts), timeout=120)
    assert "PASS quartic side-packet conductors" in output


def projective_point_at(
    section: tuple[sp.Expr, sp.Expr],
    substitution: dict[sp.Symbol, int],
) -> tuple[int, int]:
    values = tuple(sp.expand(item.subs(substitution)) for item in section)
    assert values != (0, 0)
    integers = [int(value) for value in values]
    divisor = math.gcd(*(abs(value) for value in integers))
    integers = [value // divisor for value in integers]
    first_nonzero = next(value for value in integers if value)
    if first_nonzero < 0:
        integers = [-value for value in integers]
    return tuple(integers)  # type: ignore[return-value]


def verify_rational_center_sections() -> None:
    """Compactify all rational infinitely-near centers over the side P1s."""

    sections = {
        "L1 root v": (sp.Integer(1), sp.Integer(0)),
        "L1 root q": (3 * T1, 5 * T0),
        "L1 first rational": (T0, -T1),
        "L1 second rational": (5 * T0, -3 * T1),
        "L2 first": (sp.Integer(1), sp.Integer(0)),
        "L2 second": (T0, -T1),
        "L2 third": (sp.Integer(0), sp.Integer(1)),
        "L2 fourth": (sp.Integer(1), sp.Integer(0)),
        "L3 first": (3 * T0, 5 * T1),
        "L3 second": (sp.Integer(1), sp.Integer(0)),
        "L3 third": (T0**2, -T1**2),
        "L3 fourth": (sp.Integer(1), sp.Integer(0)),
    }
    expected_endpoints = {
        "L1 root v": ((1, 0), (1, 0)),
        "L1 root q": ((1, 0), (0, 1)),
        "L1 first rational": ((0, 1), (1, 0)),
        "L1 second rational": ((0, 1), (1, 0)),
        "L2 first": ((1, 0), (1, 0)),
        "L2 second": ((0, 1), (1, 0)),
        "L2 third": ((0, 1), (0, 1)),
        "L2 fourth": ((1, 0), (1, 0)),
        "L3 first": ((0, 1), (1, 0)),
        "L3 second": ((1, 0), (1, 0)),
        "L3 third": ((0, 1), (1, 0)),
        "L3 fourth": ((1, 0), (1, 0)),
    }
    for name, section in sections.items():
        first, second = section
        common = sp.gcd(
            sp.Poly(first, T0, T1, domain=sp.QQ),
            sp.Poly(second, T0, T1, domain=sp.QQ),
        )
        assert common.total_degree() == 0
        endpoints = (
            projective_point_at(section, {T0: 0, T1: 1}),
            projective_point_at(section, {T0: 1, T1: 0}),
        )
        assert endpoints == expected_endpoints[name]


def verify_kummer_endpoint_table(
    coordinates: tuple[sp.Expr, ...],
) -> None:
    vertices = {
        "V12": tuple(
            sp.expand(item.subs({X0: a, X1: b, X2: c, X3: 1}))
            for item in coordinates
        ),
        "V13": tuple(
            sp.expand(item.subs({X0: a, X1: b, X2: 1, X3: c}))
            for item in coordinates
        ),
        "V23": tuple(
            sp.expand(item.subs({X0: a, X1: 1, X2: b, X3: c}))
            for item in coordinates
        ),
    }
    facets = {
        name: tuple(sorted(compact_newton_facets(generators)))
        for name, generators in vertices.items()
    }
    assert facets["V13"] == ()

    # The side parameter is c on L1 at V12, b on L2 at V12, c on L2 at
    # V23, and b on L3 at V23.  A facet chart requires t=lambda^w_j.
    table = {
        ("L1", "V12"): ("c", tuple(weight[2] for weight in facets["V12"])),
        ("L2", "V12"): ("b", tuple(weight[1] for weight in facets["V12"])),
        ("L2", "V23"): ("c", tuple(weight[2] for weight in facets["V23"])),
        ("L3", "V23"): ("b", tuple(weight[1] for weight in facets["V23"])),
    }
    assert table == {
        ("L1", "V12"): ("c", (1, 1)),
        ("L2", "V12"): ("b", (1, 3)),
        ("L2", "V23"): ("c", (2, 1, 3, 5)),
        ("L3", "V23"): ("b", (1, 2, 2, 3)),
    }


def verify_vertex_packet_matching(
    coordinates: tuple[sp.Expr, ...],
    packets: tuple[Packet, Packet, Packet],
) -> None:
    first_packet, second_packet, third_packet = packets
    local_v12 = tuple(
        sp.expand(item.subs({X0: a, X1: b, X2: c, X3: 1}))
        for item in coordinates
    )
    order = weight_order(local_v12[1], (2, 3, 1))
    p18 = sp.factor(
        sp.cancel(
            weight_initial(local_v12[1], (2, 3, 1), order)
            / (a**2 + b * c)
        )
    )

    # V12, L1 endpoint c=1: the cubic and quadratic side packets are
    # exactly the two P18 tangent packets from the corner atlas.
    pc = sp.expand(p18.subs(c, 1))
    zero_tangent = sp.cancel(pc.subs(b, a * d) / a**5)
    zero_face = weight_initial_2(zero_tangent, (a, d), (1, 1))
    cubic_residual = sp.expand(zero_face.subs(d, z * a) / a**3)
    cubic_factor = sp.factor_list(first_packet.polynomial)[1][1][0].as_expr()
    assert_associate(cubic_residual, cubic_factor)

    line_tangent = sp.cancel(pc.subs(b, a * (sp.Rational(5, 3) + d)) / a**5)
    line_shift = sp.expand(line_tangent.subs(d, e - sp.Rational(3, 5) * a))
    line_face = weight_initial_2(line_shift, (a, e), (1, 2))
    quadratic_residual = sp.expand(line_face.subs(e, z * a**2) / a**4)
    assert_associate(quadratic_residual, second_packet.polynomial)

    # V12, L2 endpoint b=lambda^3: this is the only cubic Kummer chart.
    lam, A, C = sp.symbols("lam A C")
    q12 = a**2 + b * c
    assert sp.factor(
        q12.subs({a: lam**2 * A, b: lam**3, c: lam * C}) / lam**4
    ) == A**2 + C
    assert sp.factor(
        p18.subs({a: lam**2 * A, b: lam**3, c: lam * C}) / lam**18
        - p18.subs({a: A, b: 1, c: C})
    ) == 0

    # V13: the two degree-two packets have the same residue field and their
    # roots are identified by z_L1=(9/25)z_L3-27/125.
    transition = sp.Rational(9, 25) * z - sp.Rational(27, 125)
    assert sp.factor(
        second_packet.polynomial.subs(z, transition)
        + 675 * third_packet.polynomial
    ) == 0
    assert sp.factor(sp.discriminant(second_packet.polynomial, z)) > 0
    assert sp.factor(sp.discriminant(third_packet.polynomial, z)) == 84

    # V23, L3 endpoint: q23 is precisely the first generic L3 direction.
    # Along q23=0, h23 has cubic contact, explaining the next three centers.
    q23 = 5 * a * b - 3 * c
    h23 = 3 * a**3 - 5 * a * b + 3 * c
    side_substitution = {a: u / t, b: 1 / t, c: v / t}
    assert sp.factor(q23.subs(side_substitution) * t**2) == 5 * u - 3 * t * v
    pulled_h = sp.factor(h23.subs(side_substitution) * t**3)
    assert pulled_h == 3 * u**3 - 5 * t * u + 3 * t**2 * v
    assert sp.factor(pulled_h.subs(v, 5 * u / (3 * t))) == 3 * u**3

    # The L2 endpoint is the c=1 fixed direction and misses both V23 colors.
    assert q23.subs({a: 0, b: 0, c: 1}) == -3
    assert h23.subs({a: 0, b: 0, c: 1}) == 3


def weight_initial_2(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, sp.Symbol],
    weight: tuple[int, int],
) -> sp.Expr:
    polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
    order = min(
        weight[0] * monomial[0] + weight[1] * monomial[1]
        for monomial, _ in polynomial.terms()
    )
    return sp.expand(
        sum(
            coefficient * variables[0] ** monomial[0] * variables[1] ** monomial[1]
            for monomial, coefficient in polynomial.terms()
            if weight[0] * monomial[0] + weight[1] * monomial[1] == order
        )
    )


def main() -> None:
    coordinates = quartic_projective_coordinates()
    packets = extract_packets(coordinates)
    verify_packet_geometry(packets)
    verify_packet_conductors(packets)
    verify_rational_center_sections()
    verify_kummer_endpoint_table(coordinates)
    verify_vertex_packet_matching(coordinates, packets)
    print(
        "PASS quartic Rees side-corner matching: packet bidegrees "
        "(4,4),(6,2),(6,2), endpoint deltas 6,3,3, V12 P18 matches, "
        "V13 quadratic transition, and V23 q/h contact verified"
    )


if __name__ == "__main__":
    main()
