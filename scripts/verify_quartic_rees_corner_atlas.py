#!/usr/bin/env python3
"""Verify the finite corner atlas for the quartic normalized Rees problem.

The generic side clusters and generic facet colors are checked by
``verify_quartic_rees_stratification.py``.  This continuation treats the
places where those strata meet:

* it constructs the complete compact Newton-fan incidence at V12 and V23;
* it proves that all eight maximal toric orbits are base-point free;
* it resolves the two boundary singularities of the weighted-degree-18
  V12 color and computes their exact conductors;
* it checks the two V12 facet-overlap cusp charts and their explicit
  quadratic normalizations; and
* it verifies that the V23 q-colors, and then q and h, glue smoothly across
  the two colored facet overlaps.

This is a vertex-corner normalization result.  It does not yet glue these
charts to the closures of all three generic side clusters, compute the
global bigraded Hilbert polynomial, or prove simultaneous normalization in
the weighted parameter.
"""

from __future__ import annotations

import itertools
import math
import sys
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
    minimal_support_points,
    weight_initial,
    weight_order,
)


a, b, c = sp.symbols("a b c")
d, e, s, t, u, v, x, y, z = sp.symbols("d e s t u v x y z")


def affine_rank(points: tuple[tuple[int, int, int], ...]) -> int:
    if len(points) < 2:
        return 0
    return sp.Matrix(
        [
            [point[index] - points[0][index] for index in range(3)]
            for point in points[1:]
        ]
    ).rank()


def compact_facet_faces(
    generators: tuple[sp.Expr, ...],
) -> dict[tuple[int, int, int], tuple[tuple[int, int, int], ...]]:
    points = minimal_support_points(generators)
    return {
        weight: tuple(
            point
            for point in points
            if sum(weight[index] * point[index] for index in range(3))
            == order
        )
        for weight, order in compact_newton_facets(generators).items()
    }


def verify_facet_adjacency(
    local_ideals: dict[str, tuple[sp.Expr, ...]],
) -> None:
    expected = {
        "V12": {
            frozenset(((1, 1, 1), (2, 3, 1))),
        },
        "V23": {
            frozenset(((1, 1, 2), (2, 3, 5))),
            frozenset(((1, 2, 1), (2, 3, 5))),
            frozenset(((1, 2, 3), (2, 3, 5))),
        },
    }
    for name, generators in local_ideals.items():
        faces = compact_facet_faces(generators)
        adjacent = set()
        for first, second in itertools.combinations(faces, 2):
            intersection = tuple(sorted(set(faces[first]) & set(faces[second])))
            if affine_rank(intersection) == 1:
                adjacent.add(frozenset((first, second)))
        assert adjacent == expected[name], (name, adjacent)


def primitive_relation(
    rays: tuple[tuple[int, int, int], ...],
) -> tuple[int, ...]:
    kernel = sp.Matrix.hstack(*(sp.Matrix(ray) for ray in rays)).nullspace()
    assert len(kernel) == 1
    vector = kernel[0]
    denominator = math.lcm(*(int(entry.q) for entry in vector))
    entries = [int(entry * denominator) for entry in vector]
    divisor = math.gcd(*(abs(entry) for entry in entries))
    entries = [entry // divisor for entry in entries]
    first_nonzero = next(entry for entry in entries if entry)
    if first_nonzero < 0:
        entries = [-entry for entry in entries]
    return tuple(entries)


def maximal_newton_cones(
    generators: tuple[sp.Expr, ...],
) -> dict[
    tuple[int, int, int],
    tuple[tuple[tuple[int, int, int], ...], tuple[sp.Rational, ...]],
]:
    points = minimal_support_points(generators)
    compact_rays = tuple(sorted(compact_newton_facets(generators)))
    rays = ((1, 0, 0), (0, 1, 0), (0, 0, 1)) + compact_rays
    minima = {
        ray: min(
            sum(ray[index] * point[index] for index in range(3))
            for point in points
        )
        for ray in rays
    }
    polynomials = tuple(sp.Poly(item, a, b, c, domain=sp.QQ) for item in generators)
    result = {}
    for point in points:
        active = tuple(
            sorted(
                ray
                for ray in rays
                if sum(ray[index] * point[index] for index in range(3))
                == minima[ray]
            )
        )
        if sp.Matrix.hstack(*(sp.Matrix(ray) for ray in active)).rank() != 3:
            continue
        coefficient_vector = tuple(
            polynomial.coeff_monomial(point) for polynomial in polynomials
        )
        assert any(coefficient != 0 for coefficient in coefficient_vector)
        result[point] = (active, coefficient_vector)
    return result


def verify_maximal_cones(
    local_ideals: dict[str, tuple[sp.Expr, ...]],
) -> None:
    e1, e2, e3 = (1, 0, 0), (0, 1, 0), (0, 0, 1)
    w111, w231 = (1, 1, 1), (2, 3, 1)
    w112, w121, w123, w235 = (1, 1, 2), (1, 2, 1), (1, 2, 3), (2, 3, 5)
    expected = {
        "V12": {
            (0, 6, 4): (frozenset((e1, w111, w231)), (0, sp.Rational(3, 4), 0, 0), 2),
            (8, 2, 0): (frozenset((e3, w111, w231)), (0, sp.Rational(3, 4), 0, 0), 1),
            (11, 0, 0): (frozenset((e2, e3, w231)), (0, sp.Rational(-3, 2), 0, 0), 2),
        },
        "V23": {
            (1, 3, 2): (frozenset((w112, w121, w235)), (0, 0, 1, 0), 2),
            (3, 5, 0): (frozenset((e3, w112, w235)), (0, 0, sp.Rational(25, 9), 0), 1),
            (8, 0, 1): (frozenset((e2, w121, w123, w235)), (0, 0, 0, 1), None),
            (9, 1, 0): (frozenset((e3, w123, w235)), (0, 0, 0, sp.Rational(-5, 3)), 1),
            (11, 0, 0): (frozenset((e2, e3, w123)), (0, 0, 0, 1), 1),
        },
    }
    for name, generators in local_ideals.items():
        cones = maximal_newton_cones(generators)
        assert set(cones) == set(expected[name])
        for point, (expected_rays, expected_coefficients, expected_index) in expected[name].items():
            rays, coefficients = cones[point]
            assert frozenset(rays) == expected_rays
            assert coefficients == expected_coefficients
            if expected_index is not None:
                assert len(rays) == 3
                index = abs(int(sp.Matrix.hstack(*(sp.Matrix(ray) for ray in rays)).det()))
                assert index == expected_index

    nonsimplicial_rays = maximal_newton_cones(local_ideals["V23"])[(8, 0, 1)][0]
    assert len(nonsimplicial_rays) == 4
    relation = primitive_relation(nonsimplicial_rays)
    relation_by_ray = dict(zip(nonsimplicial_rays, relation, strict=True))
    assert relation_by_ray == {
        e2: 2,
        w121: -1,
        w123: -3,
        w235: 2,
    }
    assert all(
        sum(relation[index] * nonsimplicial_rays[index][coordinate] for index in range(4)) == 0
        for coordinate in range(3)
    )


def weighted_initial_2(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, sp.Symbol],
    weight: tuple[int, int],
) -> tuple[int, sp.Expr]:
    polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
    order = min(
        weight[0] * monomial[0] + weight[1] * monomial[1]
        for monomial, _ in polynomial.terms()
    )
    initial = sp.expand(
        sum(
            coefficient * variables[0] ** monomial[0] * variables[1] ** monomial[1]
            for monomial, coefficient in polynomial.terms()
            if weight[0] * monomial[0] + weight[1] * monomial[1] == order
        )
    )
    return order, initial


def assert_associate(actual: sp.Expr, expected: sp.Expr) -> None:
    quotient = sp.cancel(actual / expected)
    assert quotient != 0 and not quotient.free_symbols, (actual, expected, quotient)


def degree_signature(expression: sp.Expr, variable: sp.Symbol) -> tuple[tuple[int, int], ...]:
    return tuple(
        (factor.degree(), power)
        for factor, power in sp.factor_list(sp.Poly(expression, variable, domain=sp.QQ))[1]
    )


def verify_p18_branches_and_contacts(p18: sp.Expr) -> None:
    q12 = a**2 + b * c
    assert sp.factor(p18.subs({b: 1, c: -a**2})) == a**12
    assert sp.factor(p18.subs({c: 1, b: -a**2})) == a**8

    assert sp.factor(p18.subs(a, 0)) == sp.Rational(3, 4) * b**5 * c**3
    assert sp.factor(p18.subs(b, 0)) == -a**8 * (9 * a - 23 * c**2) / 6
    assert sp.factor(p18.subs(c, 0)) == -3 * a**6 * (2 * a**3 - b**2) / 4

    # At b=1 all three branches first follow c=-a^2, then split with
    # c+a^2=z*a^4 and residual polynomial 3*z^3+4.
    pb_shift = sp.expand(p18.subs({b: 1, c: d - a**2}))
    order, initial = weighted_initial_2(pb_shift, (a, d), (1, 4))
    assert order == 12
    assert_associate(initial, 4 * a**12 + 3 * d**3)
    pb_residual = sp.expand(pb_shift.subs(d, z * a**4) / a**12).subs(a, 0)
    assert_associate(pb_residual, 3 * z**3 + 4)
    assert degree_signature(pb_residual, z) == ((3, 1),)
    assert sp.discriminant(pb_residual, z) != 0

    # At c=1 there is a degree-three branch packet tangent to b=0.  After
    # two ordinary transforms its residual direction polynomial is cubic.
    pc = sp.expand(p18.subs(c, 1))
    pc_zero_tangent = sp.cancel(pc.subs(b, a * d) / a**5)
    order, initial = weighted_initial_2(pc_zero_tangent, (a, d), (1, 1))
    assert order == 3
    cubic_residual = sp.expand(initial.subs(d, z * a) / a**3)
    assert_associate(cubic_residual, 25 * z**3 + 45 * z**2 + 54 * z + 46)
    assert degree_signature(cubic_residual, z) == ((3, 1),)
    assert sp.discriminant(cubic_residual, z) != 0

    # The other two branches follow b/a=5/3 and then d=-3a/5; their final
    # quadratic directions differ first in b at order three.
    pc_line = sp.cancel(pc.subs(b, a * (sp.Rational(5, 3) + d)) / a**5)
    pc_line_shift = sp.expand(pc_line.subs(d, e - sp.Rational(3, 5) * a))
    order, initial = weighted_initial_2(pc_line_shift, (a, e), (1, 2))
    assert order == 4
    quadratic_residual = sp.expand(initial.subs(e, z * a**2) / a**4)
    assert_associate(
        quadratic_residual,
        -15625 * z**2 + 15750 * z + 756,
    )
    assert degree_signature(quadratic_residual, z) == ((2, 1),)
    primitive_quadratic = primitive_integral(quadratic_residual, (z,))
    assert sp.factor(sp.discriminant(primitive_quadratic, z)) == 295312500

    # Pairwise intersection ledgers agree with the conductor colengths:
    # 3 branches of mutual contact 4, and packets (3,2) with contacts
    # 2 internally, 3 internally, and 1 across the two tangent directions.
    assert math.comb(3, 2) * 4 == 12
    assert math.comb(3, 2) * 2 + 3 * 2 + 3 == 15
    assert sp.factor(q12.subs({b: 1, c: -a**2})) == 0


def quadratic_coefficients(
    expression: sp.Expr,
    variable: sp.Symbol,
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    polynomial = sp.Poly(expression, variable)
    assert polynomial.degree() == 2
    return (
        polynomial.coeff_monomial(variable**2),
        polynomial.coeff_monomial(variable),
        polynomial.coeff_monomial(1),
    )


def verify_colored_overlaps(p18: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    q12 = a**2 + b * c

    # The two affine covers of the V12 (111)-(231) toric overlap.  The
    # strict q-color is x+1 on the unimodular e3 chart and x^2+1 on the
    # index-two e1 chart.  The P-color has a cusp over each closed point.
    q_e3 = sp.factor(q12.subs({a: y, b: y, c: x * y}) / y**2)
    p_e3 = sp.expand(p18.subs({a: y, b: y, c: x * y}) / y**8)
    p_e3 = sp.expand(p_e3.subs(x, s - 1))
    assert q_e3 == x + 1
    assert sp.factor(p_e3.subs(y, 0)) == sp.Rational(3, 4) * s**3
    assert sp.rem(sp.Poly(p_e3, s), sp.Poly(s, s)).as_expr() == y**2
    order, initial = weighted_initial_2(p_e3, (y, s), (3, 2))
    assert order == 6
    assert_associate(initial, 3 * s**3 + 4 * y**2)

    a3, b3, c3 = quadratic_coefficients(p_e3, y)
    discriminant_e3 = sp.factor(b3**2 - 4 * a3 * c3)
    assert discriminant_e3 == sp.Rational(3, 4) * s**3 * (7 * s - 4)
    normal_rhs_e3 = sp.factor(discriminant_e3 / s**2)
    assert normal_rhs_e3 == sp.Rational(3, 4) * s * (7 * s - 4)
    assert sp.gcd(sp.Poly(normal_rhs_e3, s), sp.diff(sp.Poly(normal_rhs_e3, s), s)).degree() == 0

    q_e1 = sp.factor(q12.subs({a: x * y, b: y, c: y}) / y**2)
    p_e1 = sp.expand(p18.subs({a: x * y, b: y, c: y}) / y**8)
    assert q_e1 == x**2 + 1
    assert sp.factor(p_e1.subs(y, 0)) == sp.Rational(3, 4) * q_e1**3
    remainder_e1 = sp.rem(
        sp.Poly(p_e1, x, domain=sp.QQ.frac_field(y)),
        sp.Poly(q_e1, x, domain=sp.QQ.frac_field(y)),
    ).as_expr()
    assert remainder_e1 == y**2

    a1, b1, c1 = quadratic_coefficients(p_e1, y)
    discriminant_e1 = sp.factor(b1**2 - 4 * a1 * c1)
    expected_e1 = (
        sp.Rational(3, 4)
        * x**10
        * (x**2 + 1) ** 3
        * (3 * x**2 + 7)
    )
    assert discriminant_e1 == expected_e1
    normal_rhs_e1 = sp.factor(discriminant_e1 / (x**5 * (x**2 + 1)) ** 2)
    assert sp.factor(
        normal_rhs_e1
        - sp.Rational(3, 4) * (x**2 + 1) * (3 * x**2 + 7)
    ) == 0
    assert sp.gcd(sp.Poly(normal_rhs_e1, x), sp.diff(sp.Poly(normal_rhs_e1, x), x)).degree() == 0

    # V23: q glues across (112)-(235), q and h glue across (235)-(123),
    # and the q-color misses the (235)-(121) overlap.
    q23 = 5 * a * b - 3 * c
    h23 = 3 * a**3 - 5 * a * b + 3 * c
    q_112_235 = sp.factor(
        q23.subs({a: u * v**2, b: u * v**3, c: s * u**2 * v**5})
        / (u**2 * v**5)
    )
    assert q_112_235 == 5 - 3 * s

    q_235_123 = sp.factor(
        q23.subs({a: u**2 * v, b: u**3 * v**2, c: s * u**5 * v**3})
        / (u**5 * v**3)
    )
    h_235_123 = sp.factor(
        h23.subs({a: u**2 * v, b: u**3 * v**2, c: s * u**5 * v**3})
        / (u**5 * v**3)
    )
    assert q_235_123 == 5 - 3 * s
    assert h_235_123 == 3 * s + 3 * u - 5
    assert sp.factor(h_235_123.subs(u, 0) + q_235_123) == 0

    q_235_121 = sp.factor(
        q23.subs({a: s * u**2 * v, b: s * u**3 * v**2, c: s**2 * u**5 * v})
        / (s**2 * u**5 * v)
    )
    assert q_235_121 == 5 * v**2 - 3
    assert q_235_121.subs(v, 0) != 0
    return p_e3, p_e1


def verify_conductors(
    p18: sp.Expr,
    p_e3: sp.Expr,
    p_e1: sp.Expr,
) -> None:
    pb = primitive_integral(p18.subs(b, 1), (a, c))
    pc = primitive_integral(p18.subs(c, 1), (a, b))
    cusp_e3 = primitive_integral(p_e3, (y, s))
    cusp_e1 = primitive_integral(p_e1, (y, x))
    program = rf'''
LIB "normal.lib";

proc assertReductionZero(poly f, ideal G, string label)
{{
  if (reduce(f,std(G)) != 0)
  {{
    "FAIL: "+label;
    exit(1);
  }}
}}

proc assertIdealEqual(ideal A, ideal B, string label)
{{
  int i;
  for (i=1; i<=size(A); i++)
  {{
    assertReductionZero(A[i],B,label);
  }}
  for (i=1; i<=size(B); i++)
  {{
    assertReductionZero(B[i],A,label);
  }}
}}

ring pbRing=0,(a,c),dp;
poly pb={singular_polynomial(pb)};
ideal pbOrder=pb;
ideal pbConductor=std(normalConductor(pbOrder));
ideal pbExpected=c4,a2c2+c3,a4+2a2c+c2;
assertIdealEqual(pbConductor,pbExpected,"P18 b=1 conductor");
if (vdim(std(pbOrder+pbConductor)) != 12)
{{
  "FAIL: P18 b=1 conductor colength";
  exit(1);
}}
assertIdealEqual(radical(pbConductor),ideal(a,c),"P18 b=1 conductor support");

ring pcRing=0,(a,b),dp;
poly pc={singular_polynomial(pc)};
ideal pcOrder=pc;
ideal pcConductor=std(normalConductor(pcOrder));
ideal pcExpected=
  25a2b2-30ab3+9b4,
  81b5-625ab3+375b4,
  27ab4-125ab3+75b4,
  5a4b-5ab3+3b4,
  a6;
assertIdealEqual(pcConductor,pcExpected,"P18 c=1 conductor");
if (vdim(std(pcOrder+pcConductor)) != 15)
{{
  "FAIL: P18 c=1 conductor colength";
  exit(1);
}}
assertIdealEqual(radical(pcConductor),ideal(a,b),"P18 c=1 conductor support");

ring e3Ring=0,(y,s),dp;
poly cusp3={singular_polynomial(cusp_e3)};
ideal cusp3Order=cusp3;
ideal cusp3Conductor=std(normalConductor(cusp3Order));
assertIdealEqual(cusp3Conductor,ideal(y,s),"rational cusp conductor");
if (vdim(std(cusp3Order+cusp3Conductor)) != 1)
{{
  "FAIL: rational cusp delta";
  exit(1);
}}

ring e1Ring=0,(y,x),dp;
poly cusp1={singular_polynomial(cusp_e1)};
ideal cusp1Order=cusp1;
ideal cusp1Conductor=std(normalConductor(cusp1Order));
assertIdealEqual(cusp1Conductor,ideal(y,x2+1),"quadratic cusp conductor");
if (vdim(std(cusp1Order+cusp1Conductor)) != 2)
{{
  "FAIL: quadratic cusp delta";
  exit(1);
}}
assertIdealEqual(
  radical(cusp1Conductor),ideal(y,x2+1),"quadratic cusp support"
);

"PASS quartic corner conductors";
'''
    output = run_singular(program, timeout=120)
    assert "PASS quartic corner conductors" in output


def main() -> None:
    coordinates = quartic_projective_coordinates()
    local_ideals = {
        "V12": tuple(
            sp.expand(item.subs({X0: a, X1: b, X2: c, X3: 1}))
            for item in coordinates
        ),
        "V23": tuple(
            sp.expand(item.subs({X0: a, X1: 1, X2: b, X3: c}))
            for item in coordinates
        ),
    }
    verify_facet_adjacency(local_ideals)
    verify_maximal_cones(local_ideals)

    order = weight_order(local_ideals["V12"][1], (2, 3, 1))
    initial = weight_initial(local_ideals["V12"][1], (2, 3, 1), order)
    p18 = sp.factor(sp.cancel(initial / (a**2 + b * c)))
    assert weight_order(p18, (2, 3, 1)) == 18
    verify_p18_branches_and_contacts(p18)
    p_e3, p_e1 = verify_colored_overlaps(p18)
    verify_conductors(p18, p_e3, p_e1)

    print(
        "PASS quartic Rees corner atlas: 8 maximal Newton cones are "
        "base-point free; P18 boundary deltas (12,15), residual cusp "
        "deltas (1,2), and V23 colored overlaps verified"
    )


if __name__ == "__main__":
    main()
