#!/usr/bin/env python3
"""Exact regressions for the fixed-quintic Galois-incidence atlas.

This checker verifies the fresh simple-discriminant branch used to obtain
pairwise nonisomorphic S5 fields and independently enumerates the cycle
types controlling real signatures and unramified local factorizations.
The large pair-sum and Dummit identities remain in their canonical existing
checkers; this file deliberately does not duplicate them.
"""

from __future__ import annotations

import os
import warnings

os.environ.setdefault("SYMPY_GROUND_TYPES", "python")

import sympy as sp
from sympy.utilities.exceptions import SymPyDeprecationWarning

warnings.filterwarnings("ignore", category=SymPyDeprecationWarning)


T, Pi, B, C = sp.symbols("T Pi B C")

f = T**5 - 5 * T**3 - 2 * Pi * B * T**2 + 4 * Pi**3 * T - 2 * Pi**5 * C

Delta = (
    432 * B**5 * C * Pi**2
    - 432 * B**4 * Pi**2
    + 12600 * B**3 * C * Pi**3
    - 2000 * B**3 * C
    + 9000 * B**2 * C**2 * Pi**7
    + 20625 * B**2 * C**2 * Pi**4
    - 11520 * B**2 * Pi**3
    + 2000 * B**2
    + 18750 * B * C**3 * Pi**8
    - 25600 * B * C * Pi**7
    + 56000 * B * C * Pi**4
    - 45000 * B * C * Pi
    + 3125 * C**4 * Pi**12
    - 40000 * C**2 * Pi**8
    + 112500 * C**2 * Pi**5
    - 84375 * C**2 * Pi**2
    + 16384 * Pi**7
    - 51200 * Pi**4
    + 40000 * Pi
)

assert sp.expand(sp.discriminant(f, T) - 16 * Pi**8 * Delta) == 0

node = {Pi: 1, B: -sp.Rational(3, 2), C: sp.Rational(3, 2)}
assert sp.factor(f.subs(node)) == (T - 1) ** 2 * (T + 1) * (T**2 + T - 3)
assert [sp.diff(Delta, variable).subs(node) for variable in (Pi, B, C)] == [
    0,
    -468,
    -468,
]

p = sp.symbols("p")
perturbed_delta = sp.Poly(
    sp.expand(Delta.subs({**node, C: sp.Rational(3, 2) + p})),
    p,
    domain=sp.QQ,
)
assert perturbed_delta.eval(0) == 0
assert perturbed_delta.diff().eval(0) == -468


Permutation = tuple[int, ...]


def compose(left: Permutation, right: Permutation) -> Permutation:
    """Return left after right."""

    return tuple(left[right[index]] for index in range(5))


def generated_group(*generators: Permutation) -> set[Permutation]:
    identity = tuple(range(5))
    group = {identity}
    frontier = [identity]
    while frontier:
        element = frontier.pop()
        for generator in generators:
            product = compose(generator, element)
            if product not in group:
                group.add(product)
                frontier.append(product)
    return group


def cycle_type(permutation: Permutation) -> tuple[int, ...]:
    unseen = set(range(5))
    lengths = []
    while unseen:
        start = next(iter(unseen))
        current = start
        length = 0
        while current in unseen:
            unseen.remove(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


cycle_5 = (1, 2, 3, 4, 0)
reflection = (0, 4, 3, 2, 1)
multiplier_2 = (0, 2, 4, 1, 3)
cycle_3 = (1, 2, 0, 3, 4)
transposition = (1, 0, 2, 3, 4)

groups = {
    "C5": generated_group(cycle_5),
    "D5": generated_group(cycle_5, reflection),
    "F20": generated_group(cycle_5, multiplier_2),
    "A5": generated_group(cycle_5, cycle_3),
    "S5": generated_group(cycle_5, transposition),
}

assert {name: len(group) for name, group in groups.items()} == {
    "C5": 5,
    "D5": 10,
    "F20": 20,
    "A5": 60,
    "S5": 120,
}

expected_cycle_types = {
    "C5": {(1, 1, 1, 1, 1), (5,)},
    "D5": {(1, 1, 1, 1, 1), (5,), (2, 2, 1)},
    "F20": {(1, 1, 1, 1, 1), (5,), (2, 2, 1), (4, 1)},
    "A5": {(1, 1, 1, 1, 1), (5,), (3, 1, 1), (2, 2, 1)},
    "S5": {
        (1, 1, 1, 1, 1),
        (5,),
        (4, 1),
        (3, 2),
        (3, 1, 1),
        (2, 2, 1),
        (2, 1, 1, 1),
    },
}
assert {
    name: {cycle_type(element) for element in group}
    for name, group in groups.items()
} == expected_cycle_types

real_roots_for_involution_type = {
    (1, 1, 1, 1, 1): 5,
    (2, 1, 1, 1): 3,
    (2, 2, 1): 1,
}
assert {
    name: {
        real_roots_for_involution_type[kind]
        for kind in types
        if kind in real_roots_for_involution_type
    }
    for name, types in expected_cycle_types.items()
} == {
    "C5": {5},
    "D5": {1, 5},
    "F20": {1, 5},
    "A5": {1, 5},
    "S5": {1, 3, 5},
}

# Rational F20 surface from the De Moivre family.
parameter_t, parameter_pi, source_h, source_x = sp.symbols(
    "parameter_t parameter_pi source_h source_x"
)
surface_u = (1 - parameter_t**2) / (1 + parameter_t**2)
surface_v = 2 * parameter_t / (1 + parameter_t**2)
assert sp.factor(surface_u**2 + surface_v**2 - 1) == 0

surface_h = (
    sp.Rational(4, 5) * parameter_pi**3
    - 1
    + 3 * surface_u**2 * surface_v**2
) / (surface_u * surface_v)

de_moivre = source_x**5 - 5 * source_x**3 + 5 * source_x + source_h
assert sp.expand(
    sp.discriminant(de_moivre, source_x)
    - 3125 * (source_h**2 - 4) ** 2
) == 0

de_moivre_companion = sp.zeros(5)
for column in range(4):
    de_moivre_companion[column + 1, column] = 1
de_moivre_companion[:, 4] = sp.Matrix(
    (-source_h, -5, 0, 5, 0)
)

abstract_u, abstract_v = sp.symbols("abstract_u abstract_v")
abstract_eta = (
    abstract_u * de_moivre_companion
    + abstract_v * (de_moivre_companion**2 - 2 * sp.eye(5))
)
assert sp.trace(abstract_eta) == 0
assert sp.expand(
    sp.trace(abstract_eta**2)
    - 10 * (abstract_u**2 + abstract_v**2)
) == 0

eta_charpoly = sp.Poly(abstract_eta.charpoly(T).as_expr(), T)
eta_q = 5 * abstract_u * abstract_v * (
    abstract_v * source_h - 2 * abstract_u
)
eta_r = 5 * (
    (abstract_u**2 + abstract_v**2) ** 2
    - 3 * abstract_u**2 * abstract_v**2
    + abstract_u
    * abstract_v
    * (abstract_u**2 + abstract_v**2)
    * source_h
)
eta_s = (
    abstract_u**5 * source_h
    + 10 * abstract_u**4 * abstract_v
    + 5 * abstract_u**3 * abstract_v**2 * source_h
    - 10 * abstract_u**2 * abstract_v**3
    - 5 * abstract_u * abstract_v**4 * source_h
    + 2 * abstract_v**5
    - abstract_v**5 * source_h**2
)
expected_eta_charpoly = (
    T**5
    - 5 * (abstract_u**2 + abstract_v**2) * T**3
    + eta_q * T**2
    + eta_r * T
    + eta_s
)
assert sp.expand(eta_charpoly.as_expr() - expected_eta_charpoly) == 0
assert sp.factor(
    eta_r.subs(
        {
            abstract_u: surface_u,
            abstract_v: surface_v,
            source_h: surface_h,
        }
    )
    - 4 * parameter_pi**3
) == 0


def factor_degrees_mod_prime(poly: sp.Poly, prime: int) -> tuple[int, ...]:
    _, cleared = sp.Poly(
        poly.as_expr(), poly.gens[0], domain=sp.QQ
    ).clear_denoms(convert=True)
    _, poly = cleared.primitive()
    assert int(poly.LC()) % prime
    assert int(sp.discriminant(poly.as_expr(), poly.gens[0])) % prime
    factors = sp.factor_list(poly.as_expr(), modulus=prime)[1]
    return tuple(
        sorted(
            (
                int(sp.degree(factor, poly.gens[0]))
                for factor, exponent in factors
                for _ in range(exponent)
            ),
            reverse=True,
        )
    )


surface_point = {
    parameter_t: sp.Rational(1, 2),
    parameter_pi: 1,
}
surface_u_point = sp.factor(surface_u.subs(surface_point))
surface_v_point = sp.factor(surface_v.subs(surface_point))
surface_h_point = sp.factor(surface_h.subs(surface_point))
assert (surface_u_point, surface_v_point, surface_h_point) == (
    sp.Rational(3, 5),
    sp.Rational(4, 5),
    sp.Rational(307, 300),
)

source_witness = sp.Poly(
    300 * source_x**5 - 1500 * source_x**3 + 1500 * source_x + 307,
    source_x,
    domain=sp.ZZ,
)
assert factor_degrees_mod_prime(source_witness, 11) == (5,)

surface_q_point = sp.factor(
    eta_q.subs(
        {
            abstract_u: surface_u_point,
            abstract_v: surface_v_point,
            source_h: surface_h_point,
        }
    )
)
surface_s_point = sp.factor(
    eta_s.subs(
        {
            abstract_u: surface_u_point,
            abstract_v: surface_v_point,
            source_h: surface_h_point,
        }
    )
)
assert (surface_q_point, surface_s_point) == (
    -sp.Rational(572, 625),
    -sp.Rational(67834669, 70312500),
)
surface_target_b = sp.factor(
    -eta_q.subs(
        {
            abstract_u: surface_u,
            abstract_v: surface_v,
            source_h: surface_h,
        }
    )
    / (2 * parameter_pi)
)
assert sp.factor(
    sp.diff(surface_target_b, parameter_t).subs(surface_point)
) == -sp.Rational(1884, 625)

target_witness = sp.Poly(
    70312500 * T**5
    - 351562500 * T**3
    - 64350000 * T**2
    + 281250000 * T
    - 67834669,
    T,
    domain=sp.ZZ,
)
assert factor_degrees_mod_prime(target_witness, 11) == (5,)
assert {
    prime: factor_degrees_mod_prime(target_witness, prime)
    for prime in (7, 11, 29, 89)
} == {
    7: (4, 1),
    11: (5,),
    29: (2, 2, 1),
    89: (1, 1, 1, 1, 1),
}

assert sp.Poly(de_moivre.subs(source_h, sp.Rational(307, 300)), source_x).count_roots(
    -sp.oo, sp.oo
) == 5
assert sp.Poly(de_moivre.subs(source_h, sp.Rational(1269, 100)), source_x).count_roots(
    -sp.oo, sp.oo
) == 1

# The tame double-collision used for fresh ramification.
assert sp.factor(de_moivre.subs(source_h, 2)) == (
    (source_x + 2) * (source_x**2 - source_x - 1) ** 2
)
fixed_t_h = sp.factor(surface_h.subs(parameter_t, sp.Rational(1, 2)))
assert sp.expand(
    fixed_t_h
    - sp.Rational(5, 3) * parameter_pi**3
    + sp.Rational(193, 300)
) == 0

# Rational D5 curve from a Brumer specialization.
d5_r, pair_x = sp.symbols("d5_r pair_x", nonzero=True)
d5_pi = 1 / (40 * d5_r**4)
d5_target_b = (400 * d5_r**6 - 1) / (80 * d5_r**5)
d5_target_c = -64 * d5_r**5
d5_polynomial = sp.Poly(
    T**5
    - 5 * T**3
    - (400 * d5_r**6 - 1) / (1600 * d5_r**9) * T**2
    + 1 / (16000 * d5_r**12) * T
    + 1 / (800000 * d5_r**15),
    T,
)
assert sp.expand(
    d5_polynomial.as_expr()
    - (
        T**5
        - 5 * T**3
        - 2 * d5_pi * d5_target_b * T**2
        + 4 * d5_pi**3 * T
        - 2 * d5_pi**5 * d5_target_c
    )
) == 0

d5_discriminant_numerator = (
    256000000 * d5_r**18
    + 160000 * d5_r**12
    - 800 * d5_r**6
    - 7
)
assert sp.factor(
    sp.discriminant(d5_polynomial.as_expr(), T)
    - (
        d5_discriminant_numerator
        / (40960000000 * d5_r**30)
    )
    ** 2
) == 0

d5_a = d5_polynomial.coeff_monomial(T**2)
d5_b = d5_polynomial.coeff_monomial(T)
d5_c = d5_polynomial.coeff_monomial(1)
d5_pair_resolvent = sp.Poly(
    pair_x**10
    - 15 * pair_x**8
    + d5_a * pair_x**7
    + (75 - 3 * d5_b) * pair_x**6
    + (-10 * d5_a - 11 * d5_c) * pair_x**5
    + (-d5_a**2 + 10 * d5_b - 125) * pair_x**4
    + (-4 * d5_a * d5_b + 25 * d5_a + 20 * d5_c) * pair_x**3
    + (
        5 * d5_a**2
        + 7 * d5_a * d5_c
        - 4 * d5_b**2
        + 25 * d5_b
    )
    * pair_x**2
    + (-d5_a**3 + 4 * d5_b * d5_c - 25 * d5_c) * pair_x
    - d5_a**2 * d5_b
    - 5 * d5_a * d5_c
    - d5_c**2,
    pair_x,
)
d5_pair_u = (
    3200000 * pair_x**5 * d5_r**15
    - 16000000 * pair_x**3 * d5_r**15
    - 40000 * pair_x**3 * d5_r**9
    - 800000 * pair_x**2 * d5_r**12
    - 2000 * pair_x**2 * d5_r**6
    + 40000 * pair_x * d5_r**9
    + 2000 * d5_r**6
    - 7
)
d5_pair_v = (
    1600000 * pair_x**5 * d5_r**15
    - 16000000 * pair_x**3 * d5_r**15
    + 20000 * pair_x**3 * d5_r**9
    + 2000 * pair_x**2 * d5_r**6
    + 40000000 * pair_x * d5_r**15
    - 120000 * pair_x * d5_r**9
    - 50 * pair_x * d5_r**3
    - 6000 * d5_r**6
    + 19
)
assert sp.expand(
    d5_pair_resolvent.as_expr()
    - d5_pair_u
    * d5_pair_v
    / (5120000000000 * d5_r**30)
) == 0

d5_point = {d5_r: -sp.Rational(1, 2)}
assert (
    sp.factor(d5_pi.subs(d5_point)),
    sp.factor(d5_target_b.subs(d5_point)),
    sp.factor(d5_target_c.subs(d5_point)),
) == (
    sp.Rational(2, 5),
    -sp.Rational(21, 10),
    sp.Integer(2),
)
d5_witness = sp.Poly(
    3125 * T**5 - 15625 * T**3 + 5250 * T**2 + 800 * T - 128,
    T,
    domain=sp.ZZ,
)
assert {
    prime: factor_degrees_mod_prime(d5_witness, prime)
    for prime in (3, 11, 23)
} == {
    3: (5,),
    11: (2, 2, 1),
    23: (1, 1, 1, 1, 1),
}
d5_pair_u_point = sp.Poly(
    d5_pair_u.subs(d5_point), pair_x, domain=sp.QQ
)
d5_pair_v_point = sp.Poly(
    d5_pair_v.subs(d5_point), pair_x, domain=sp.QQ
)
assert factor_degrees_mod_prime(d5_pair_u_point, 3) == (5,)
assert factor_degrees_mod_prime(d5_pair_v_point, 3) == (5,)
assert sp.Poly(
    d5_polynomial.as_expr().subs(d5_r, sp.Rational(1, 2)),
    T,
).count_roots(-sp.oo, sp.oo) == 5
assert sp.Poly(
    d5_polynomial.as_expr().subs(d5_r, sp.Rational(1, 4)),
    T,
).count_roots(-sp.oo, sp.oo) == 1

d5_lambda = sp.symbols("d5_lambda")
d5_source = sp.Poly(
    source_x**5
    - 5 * source_x**4
    + 5 * (2 - d5_lambda**2) * source_x**3
    + 5 * (2 * d5_lambda**2 - 1) * source_x**2
    + 5 * (1 - d5_lambda**2) * source_x
    - 2,
    source_x,
)
d5_branch = 4 * d5_lambda**6 + d5_lambda**4 - 2 * d5_lambda**2 - 7
assert sp.expand(
    sp.discriminant(d5_source.as_expr(), source_x)
    - 62500 * d5_branch**2
) == 0
d5_fresh_polynomial = sp.Poly(
    d5_branch.subs(d5_lambda, 20 * d5_r**3),
    d5_r,
)
assert sp.gcd(d5_fresh_polynomial, d5_fresh_polynomial.diff()).degree() == 0
d5_branch_mod_13 = sp.Poly(
    d5_source.as_expr().subs(d5_lambda, 4),
    source_x,
    modulus=13,
)
assert sorted(
    (
        int(sp.degree(factor, source_x)),
        exponent,
    )
    for factor, exponent in sp.factor_list(
        d5_branch_mod_13.as_expr(), modulus=13
    )[1]
) == [(1, 1), (1, 2), (1, 2)]

# Smaller alternating descent surface.
a5_k, a5_u, a5_v, a5_pi = sp.symbols(
    "a5_k a5_u a5_v a5_pi", nonzero=True
)
a5_a = 5 * a5_k**2 - 1
a5_source = sp.Poly(
    source_x**5 + 5 * a5_a * source_x + 4 * a5_a,
    source_x,
)
assert sp.factor(
    sp.discriminant(a5_source.as_expr(), source_x)
    - (2000 * a5_k * a5_a**2) ** 2
) == 0
a5_specialization = sp.Poly(
    a5_source.as_expr().subs(a5_k, 1),
    source_x,
)
assert factor_degrees_mod_prime(a5_specialization, 3) == (5,)
assert factor_degrees_mod_prime(a5_specialization, 7) == (3, 1, 1)

a5_companion = sp.zeros(5)
for column in range(4):
    a5_companion[column + 1, column] = 1
a5_companion[:, 4] = sp.Matrix((-4 * a5_a, -5 * a5_a, 0, 0, 0))
a5_eta = a5_u * a5_companion + a5_v * a5_companion**3
a5_charpoly = sp.Poly(a5_eta.charpoly(T).as_expr(), T)
a5_expected = (
    T**5
    + 20 * a5_a * a5_u * a5_v * T**3
    + 20
    * a5_a
    * a5_v
    * (a5_u**2 - 3 * a5_a * a5_v**2)
    * T**2
    + 5
    * a5_a
    * (
        25 * a5_a**2 * a5_v**4
        + 10 * a5_a * a5_u**2 * a5_v**2
        - 16 * a5_a * a5_u * a5_v**3
        + a5_u**4
    )
    * T
    + 4
    * a5_a
    * (
        25 * a5_a**2 * a5_u * a5_v**4
        + 16 * a5_a**2 * a5_v**5
        + 10 * a5_a * a5_u**3 * a5_v**2
        + a5_u**5
    )
)
assert sp.expand(a5_charpoly.as_expr() - a5_expected) == 0
assert sp.trace(a5_eta) == 0
assert sp.factor(
    sp.trace(a5_eta**2) + 40 * a5_a * a5_u * a5_v
) == 0

# Eliminate the first A5 normalization equation.
a5_s = sp.symbols("a5_s", nonzero=True)
a5_u_eliminated = a5_s / (2 * a5_a)
a5_v_eliminated = -1 / (2 * a5_s)
assert sp.factor(
    a5_a * a5_u_eliminated * a5_v_eliminated
    + sp.Rational(1, 4)
) == 0
a5_pi_cube = sp.factor(
    sp.Rational(5, 4)
    * a5_a
    * (
        25 * a5_a**2 * a5_v**4
        + 10 * a5_a * a5_u**2 * a5_v**2
        - 16 * a5_a * a5_u * a5_v**3
        + a5_u**4
    )
)
a5_pi_cube_eliminated = (
    5
    * (
        25 * a5_a**6
        + 16 * a5_a**4 * a5_s**2
        + 10 * a5_a**3 * a5_s**4
        + a5_s**8
    )
    / (64 * a5_a**3 * a5_s**4)
)
assert sp.factor(
    a5_pi_cube.subs(
        {
            a5_u: a5_u_eliminated,
            a5_v: a5_v_eliminated,
        }
    )
    - a5_pi_cube_eliminated
) == 0

# Polynomial sections of the A5 cubic cover have an exact degree gate.
a5_curve_k = sp.symbols("a5_curve_k")
a5_curve_a = 5 * a5_curve_k**2 - 1
a5_curve_c0, a5_curve_c1, a5_curve_c2 = sp.symbols(
    "a5_curve_c0 a5_curve_c1 a5_curve_c2"
)


def a5_branch_for_s(a5_curve_s: sp.Expr) -> sp.Expr:
    return (
        25 * a5_curve_a**6
        + 16 * a5_curve_a**4 * a5_curve_s**2
        + 10 * a5_curve_a**3 * a5_curve_s**4
        + a5_curve_s**8
    )


a5_linear_s = a5_curve_c1 * a5_curve_k + a5_curve_c0
a5_quadratic_s = (
    a5_curve_c2 * a5_curve_k**2
    + a5_curve_c1 * a5_curve_k
    + a5_curve_c0
)
assert sp.degree(
    a5_branch_for_s(a5_linear_s),
    a5_curve_k,
) == 12
assert sp.degree(
    a5_curve_a**3 * a5_linear_s**4,
    a5_curve_k,
) == 10
assert sp.degree(
    a5_branch_for_s(a5_quadratic_s),
    a5_curve_k,
) == 16
assert sp.degree(
    a5_curve_a**3 * a5_quadratic_s**4,
    a5_curve_k,
) == 14

a5_constant_s = sp.symbols("a5_constant_s")
a5_cube_a, a5_cube_b, a5_cube_c = sp.symbols(
    "a5_cube_a a5_cube_b a5_cube_c"
)
a5_constant_branch = a5_branch_for_s(a5_constant_s)
a5_cube_polynomial = (
    a5_cube_a * a5_curve_k**4
    + a5_cube_b * a5_curve_k**2
    + a5_cube_c
)
a5_constant_identity = sp.Poly(
    sp.expand(
        a5_cube_polynomial**3
        - 5 * a5_constant_s**2 * a5_constant_branch
    ),
    a5_curve_k,
)
a5_constant_ideal = sp.groebner(
    [
        coefficient
        for coefficient in a5_constant_identity.all_coeffs()
        if coefficient != 0
    ],
    a5_cube_a,
    a5_cube_b,
    a5_cube_c,
    a5_constant_s,
    order="lex",
)
assert any(
    sp.factor(polynomial.as_expr()) == a5_constant_s**8
    for polynomial in a5_constant_ideal.polys
)

# Divisor multiplicities force s=c*h^3.  In degrees three and six,
# finite-place and infinity valuations of the cleared cube identity
# are incompatible.
a5_curve_r, a5_degree_three_c = sp.symbols(
    "a5_curve_r a5_degree_three_c",
    nonzero=True,
)
a5_degree_three_s = (
    a5_degree_three_c * (a5_curve_k - a5_curve_r) ** 3
)
a5_degree_three_branch = a5_branch_for_s(a5_degree_three_s)
assert sp.factor(
    a5_degree_three_branch.subs(a5_curve_k, a5_curve_r)
    - 25 * (5 * a5_curve_r**2 - 1) ** 6
) == 0
assert sp.Poly(
    a5_degree_three_branch,
    a5_curve_k,
).LC() == a5_degree_three_c**8
# At the factor a=5k^2-1, the branch valuation is six for every
# positive exponent e of a in s; the coefficient tie at e=1 is
# 25+16(s/a)^2 and cannot vanish in the real residue field Q(sqrt(5)).
assert [
    min(6, 4 + 2 * exponent, 3 + 4 * exponent, 8 * exponent)
    for exponent in range(1, 7)
] == [6, 6, 6, 6, 6, 6]
assert [
    exponent_class
    for exponent_class in range(3)
    if (2 * exponent_class) % 3 == 0
] == [0]
# In a field of degree at most two, ramification indices above 5 are
# one or two.  Either forces v_5(c)=0 mod 3 when c^2 is a cube.
for a5_ramification_index in (1, 2):
    assert [
        valuation_class
        for valuation_class in range(3)
        if (
            2 * a5_ramification_index * valuation_class
        )
        % 3
        == 0
    ] == [0]
assert [
    valuation_class
    for valuation_class in range(3)
    if valuation_class == 0
    and (1 + 10 * valuation_class) % 3 == 0
] == []

# Exact A5 anchors in both real chambers.
a5_total_integral = sp.Poly(
    5 * source_x**5 - 25 * source_x**3 + 20 * source_x + 4,
    source_x,
)
assert sp.discriminant(a5_total_integral.as_expr(), source_x) == 145000**2
assert a5_total_integral.count_roots(-sp.oo, sp.oo) == 5
assert {
    prime: factor_degrees_mod_prime(a5_total_integral, prime)
    for prime in (3, 17, 23, 211)
} == {
    3: (5,),
    17: (2, 2, 1),
    23: (3, 1, 1),
    211: (1, 1, 1, 1, 1),
}

a5_one_real = sp.Poly(
    source_x**5
    - 5 * source_x**3
    - 16 * source_x**2
    + 32 * source_x
    + 16,
    source_x,
)
assert sp.discriminant(a5_one_real.as_expr(), source_x) == 54176**2
assert a5_one_real.count_roots(-sp.oo, sp.oo) == 1
assert {
    prime: factor_degrees_mod_prime(a5_one_real, prime)
    for prime in (3, 7, 13, 389)
} == {
    3: (5,),
    7: (3, 1, 1),
    13: (2, 2, 1),
    389: (1, 1, 1, 1, 1),
}

# Mestre source pencils through the two A5 anchors.
a5_t = sp.symbols("a5_t")
a5_total_p = (
    source_x**5
    - 5 * source_x**3
    + 4 * source_x
    + sp.Rational(4, 5)
)
a5_total_q = (
    -source_x**4
    + source_x**3
    + 2 * source_x**2
    - source_x
    - sp.Rational(21, 25)
)
a5_total_r = (
    source_x**4
    - source_x**3
    - source_x**2
    + source_x
    + sp.Rational(8, 5)
)
assert sp.expand(
    a5_total_p * sp.diff(a5_total_q, source_x)
    - sp.diff(a5_total_p, source_x) * a5_total_q
    - a5_total_r**2
) == 0
a5_total_disc_factor = (
    427 * a5_t**4
    + 2335 * a5_t**3
    + 7925 * a5_t**2
    + 12125 * a5_t
    + 14500
)
assert sp.factor(
    sp.discriminant(a5_total_p - a5_t * a5_total_q, source_x)
    - sp.Rational(4, 15625) * a5_total_disc_factor**2
) == 0
assert sp.Poly(a5_total_disc_factor, a5_t).count_roots(
    -sp.oo, sp.oo
) == 0

a5_one_q = (
    -source_x**4
    + sp.Rational(9499, 1790) * source_x**3
    - sp.Rational(36290687, 3204100) * source_x**2
    + sp.Rational(23868087, 1602050) * source_x
    - sp.Rational(37781321, 3204100)
)
a5_one_r = (
    source_x**4
    - sp.Rational(9499, 1790) * source_x**3
    + sp.Rational(4841, 895) * source_x**2
    + sp.Rational(26683, 1790) * source_x
    - sp.Rational(22208, 895)
)
assert sp.expand(
    a5_one_real.as_expr() * sp.diff(a5_one_q, source_x)
    - sp.diff(a5_one_real.as_expr(), source_x) * a5_one_q
    - a5_one_r**2
) == 0
a5_one_disc_factor = (
    599061012975492710 * a5_t**4
    + 18908224953293263933 * a5_t**3
    + 193635834698041154800 * a5_t**2
    + 591172039750728710000 * a5_t
    - 526305815118736000000
)
a5_one_disc_root = (
    sp.Rational(1693, 16447056722460500000)
    * a5_one_disc_factor
)
assert sp.factor(
    sp.discriminant(
        a5_one_real.as_expr() - a5_t * a5_one_q,
        source_x,
    )
    - a5_one_disc_root**2
) == 0
assert sp.Poly(a5_one_disc_factor, a5_t).count_roots(
    -sp.Rational(1, 2), sp.Rational(1, 2)
) == 0

# The affine normalization cover of the total-real pencil has genus eight.
a5_z, a5_r_parameter = sp.symbols("a5_z a5_r_parameter")
a5_centered = sp.Poly(
    sp.expand(
        (a5_total_p - a5_t * a5_total_q).subs(
            source_x,
            a5_z - a5_t / 5,
        )
    ),
    a5_z,
)
a5_center_cubic = -(2 * a5_t**2 + 5 * a5_t + 25) / 5
a5_center_linear = -(
    3 * a5_t**4
    + 15 * a5_t**3
    - 25 * a5_t**2
    - 125 * a5_t
    - 500
) / 125
assert sp.factor(
    a5_centered.coeff_monomial(a5_z**3) - a5_center_cubic
) == 0
assert sp.factor(
    a5_centered.coeff_monomial(a5_z) - a5_center_linear
) == 0
a5_lift_t = 5 * (1 - 2 * a5_r_parameter) / (
    a5_r_parameter**2 - 2
)
a5_lift_lambda = -(
    a5_r_parameter**2 - a5_r_parameter + 2
) / (a5_r_parameter**2 - 2)
assert sp.factor(
    a5_lift_lambda**2
    - (
        2 * a5_lift_t**2 + 5 * a5_lift_t + 25
    )
    / 25
) == 0
a5_n8 = (
    4 * a5_r_parameter**8
    - 10 * a5_r_parameter**7
    - 7 * a5_r_parameter**6
    + 160 * a5_r_parameter**5
    - 429 * a5_r_parameter**4
    + 290 * a5_r_parameter**3
    - 23 * a5_r_parameter**2
    - 60 * a5_r_parameter
    + 59
)
a5_lift_q = a5_r_parameter**2 - a5_r_parameter + 2
a5_pi_cube_before_parametrization = -5 * (
    3 * a5_t**4
    + 15 * a5_t**3
    - 25 * a5_t**2
    - 125 * a5_t
    - 500
) / (4 * (2 * a5_t**2 + 5 * a5_t + 25) ** 2)
assert sp.factor(
    a5_pi_cube_before_parametrization.subs(a5_t, a5_lift_t)
    - a5_n8 / (4 * a5_lift_q**4)
) == 0
a5_n8_poly = sp.Poly(a5_n8, a5_r_parameter)
a5_lift_q_poly = sp.Poly(a5_lift_q, a5_r_parameter)
assert sp.gcd(a5_n8_poly, a5_n8_poly.diff()).degree() == 0
assert sp.gcd(a5_n8_poly, a5_lift_q_poly).degree() == 0
assert a5_n8_poly.degree() == 8
assert a5_lift_q_poly.degree() == 2
assert sp.factor(
    a5_n8.subs(a5_r_parameter, sp.Rational(1, 2))
    - 4
    * a5_lift_q.subs(a5_r_parameter, sp.Rational(1, 2)) ** 4
) == 0

# A quadratic non-affine generator has a trace-conic obstruction at t=infinity.
a5_centered_companion = sp.zeros(5)
for column in range(4):
    a5_centered_companion[column + 1, column] = 1
a5_centered_companion[:, 4] = sp.Matrix(
    tuple(
        -a5_centered.coeff_monomial(a5_z**exponent)
        for exponent in range(5)
    )
)
a5_centered_trace_2 = sp.factor(
    sp.trace(a5_centered_companion**2)
)
a5_h = sp.symbols("a5_h")
a5_quadratic_eta = (
    a5_centered_companion
    + a5_h
    * (
        a5_centered_companion**2
        - a5_centered_trace_2 * sp.eye(5) / 5
    )
)
a5_trace_conic = (
    18 * a5_h**2 * a5_t**4
    + 90 * a5_h**2 * a5_t**3
    + 325 * a5_h**2 * a5_t**2
    + 500 * a5_h**2 * a5_t
    + 875 * a5_h**2
    - 60 * a5_h * a5_t**3
    - 225 * a5_h * a5_t**2
    - 375 * a5_h * a5_t
    + 50 * a5_t**2
    + 125 * a5_t
    + 625
)
assert sp.factor(sp.trace(a5_quadratic_eta)) == 0
assert sp.factor(
    sp.trace(a5_quadratic_eta**2)
    - sp.Rational(2, 125) * a5_trace_conic
) == 0

a5_u, a5_h_infinity, a5_k_infinity = sp.symbols(
    "a5_u a5_h_infinity a5_k_infinity"
)
a5_trace_at_infinity = sp.factor(
    a5_u**2
    * a5_trace_conic.subs(
        {
            a5_t: 1 / a5_u,
            a5_h: a5_u * a5_h_infinity,
        }
    )
)
assert sp.factor(
    a5_trace_at_infinity.subs(a5_u, 0)
    - 2 * (3 * a5_h_infinity - 5) ** 2
) == 0
a5_trace_after_blowup = sp.Poly(
    sp.expand(
        a5_trace_at_infinity.subs(
            a5_h_infinity,
            sp.Rational(5, 3) + a5_k_infinity,
        )
    ),
    a5_u,
)
assert sp.factor(
    a5_trace_after_blowup.coeff_monomial(1)
    - 18 * a5_k_infinity**2
) == 0
assert sp.factor(
    a5_trace_after_blowup.coeff_monomial(a5_u)
    - 90 * a5_k_infinity**2
    - 75 * a5_k_infinity
) == 0
assert sp.factor(
    a5_trace_after_blowup.coeff_monomial(a5_u**2)
    - 325 * a5_k_infinity**2
    - sp.Rational(2125, 3) * a5_k_infinity
    - sp.Rational(8125, 9)
) == 0
a5_residue_c, a5_residue_d = sp.symbols(
    "a5_residue_c a5_residue_d"
)
a5_residue_form = (
    18 * a5_residue_c**2
    + 75 * a5_residue_c
    + sp.Rational(8125, 9)
)
assert sp.expand(
    72 * a5_residue_form
    - (36 * a5_residue_c + 75) ** 2
    - 59375
) == 0
assert 59375 == 5**5 * 19
assert pow(2, 2, 5) == 4

# Direct affine normalization of Lehmer's C5 family.
c5_n = sp.symbols("c5_n")
c5_coefficients = (
    c5_n**2,
    -(2 * c5_n**3 + 6 * c5_n**2 + 10 * c5_n + 10),
    c5_n**4 + 5 * c5_n**3 + 11 * c5_n**2 + 15 * c5_n + 5,
    c5_n**3 + 4 * c5_n**2 + 10 * c5_n + 10,
    sp.Integer(1),
)
c5_power_sums = {0: sp.Integer(5)}
for exponent in range(1, 5):
    c5_power_sums[exponent] = sp.factor(
        -sum(
            c5_coefficients[index - 1]
            * c5_power_sums[exponent - index]
            for index in range(1, exponent)
        )
        - exponent * c5_coefficients[exponent - 1]
    )
c5_center = -c5_power_sums[1] / 5
c5_trace_2 = sp.factor(
    c5_power_sums[2]
    + 2 * c5_center * c5_power_sums[1]
    + 5 * c5_center**2
)
c5_trace_4 = sp.factor(
    c5_power_sums[4]
    + 4 * c5_center * c5_power_sums[3]
    + 6 * c5_center**2 * c5_power_sums[2]
    + 4 * c5_center**3 * c5_power_sums[1]
    + 5 * c5_center**4
)
c5_q = c5_n**4 + 5 * c5_n**3 + 15 * c5_n**2 + 25 * c5_n + 25
assert sp.factor(c5_trace_2 - 4 * c5_q / 5) == 0
c5_e4 = sp.factor((c5_trace_2**2 - 2 * c5_trace_4) / 8)
assert sp.factor(
    c5_e4
    + (
        3 * c5_n**4 + 15 * c5_n**3 + 20 * c5_n**2 - 50
    )
    * c5_q
    / 125
) == 0
assert sp.factor(
    25 * c5_e4 / c5_trace_2**2
    + 5
    * (3 * c5_n**4 + 15 * c5_n**3 + 20 * c5_n**2 - 50)
    / (16 * c5_q)
) == 0

# A cyclic normal-basis section removes the Lehmer square obstruction.
c5_lehmer = sp.Poly(
    source_x**5
    + sum(
        coefficient * source_x ** (4 - index)
        for index, coefficient in enumerate(c5_coefficients)
    ),
    source_x,
)
c5_sigma_numerator = c5_n + 2 + c5_n * source_x - source_x**2
c5_sigma_denominator = 1 + (c5_n + 2) * source_x
c5_sigma = sp.factor(
    sp.rem(
        c5_sigma_numerator
        * sp.invert(
            c5_sigma_denominator,
            c5_lehmer.as_expr(),
            source_x,
        ),
        c5_lehmer.as_expr(),
        source_x,
    )
)
assert sp.factor(
    sp.rem(
        sp.together(
            c5_lehmer.as_expr().subs(source_x, c5_sigma)
        ).as_numer_denom()[0],
        c5_lehmer.as_expr(),
        source_x,
    )
) == 0
c5_orbit = [source_x]
for _ in range(1, 6):
    c5_next = sp.together(c5_orbit[-1].subs(source_x, c5_sigma))
    c5_next_numerator, c5_next_denominator = c5_next.as_numer_denom()
    c5_orbit.append(
        sp.factor(
            sp.rem(
                c5_next_numerator
                * sp.invert(
                    c5_next_denominator,
                    c5_lehmer.as_expr(),
                    source_x,
                ),
                c5_lehmer.as_expr(),
                source_x,
            )
        )
    )
assert sp.factor(c5_orbit[5] - source_x) == 0

c5_a_vector = (2, 1, 0, -1, -2)
c5_b_vector = (
    sp.Rational(5, 2),
    sp.Rational(5, 2),
    5,
    -sp.Rational(5, 2),
    -sp.Rational(15, 2),
)
c5_c_vector = (0, 5, 10, -10, -5)
c5_d_vector = tuple(
    c5_a_vector[index] * c5_n**2
    + c5_b_vector[index] * c5_n
    + c5_c_vector[index]
    for index in range(5)
)
assert sum(c5_d_vector) == 0
assert sp.factor(
    sum(entry**2 for entry in c5_d_vector) - 10 * c5_q
) == 0

# Explicit second-intersection parametrization of the full trace quadric.
c5_z_free = sp.symbols("c5_z0:4")
c5_z_vector = (
    *c5_z_free,
    -sum(c5_z_free),
)
c5_z_norm = sp.expand(sum(entry**2 for entry in c5_z_vector))
c5_dz = sp.expand(
    sum(
        c5_d_vector[index] * c5_z_vector[index]
        for index in range(5)
    )
)
c5_quadric_coefficients = tuple(
    c5_d_vector[index] / c5_q
    - 2 * c5_dz * c5_z_vector[index] / (c5_q * c5_z_norm)
    for index in range(5)
)
assert sp.factor(sum(c5_quadric_coefficients)) == 0
assert sp.factor(
    c5_q * sum(entry**2 for entry in c5_quadric_coefficients)
    - 10
) == 0

c5_companion = sp.zeros(5)
for column in range(4):
    c5_companion[column + 1, column] = 1
c5_companion[:, 4] = sp.Matrix(
    tuple(-coefficient for coefficient in reversed(c5_coefficients))
)


def c5_evaluate(element: sp.Expr) -> sp.Matrix:
    polynomial = sp.Poly(element, source_x)
    return sum(
        (
            polynomial.coeff_monomial(source_x**exponent)
            * c5_companion**exponent
            for exponent in range(5)
        ),
        sp.zeros(5),
    )


c5_eta_matrix = sum(
    (
        c5_d_vector[index] * c5_evaluate(c5_orbit[index])
        for index in range(5)
    ),
    sp.zeros(5),
) / c5_q
assert sp.factor(sp.trace(c5_eta_matrix)) == 0
assert sp.factor(sp.trace(c5_eta_matrix**2) - 10) == 0
c5_section_trace_4 = sp.factor(sp.trace(c5_eta_matrix**4))
c5_section_numerator = (
    64 * c5_n**12
    + 960 * c5_n**11
    + 6780 * c5_n**10
    + 31800 * c5_n**9
    + 119875 * c5_n**8
    + 385875 * c5_n**7
    + 1055000 * c5_n**6
    + 2370000 * c5_n**5
    + 4206250 * c5_n**4
    + 5600000 * c5_n**3
    + 5237500 * c5_n**2
    + 3125000 * c5_n
    + 1000000
)
c5_section_numerator_poly = sp.Poly(
    c5_section_numerator,
    c5_n,
)
assert c5_section_numerator_poly.degree() == 12
assert sp.gcd(
    c5_section_numerator_poly,
    c5_section_numerator_poly.diff(),
).degree() == 0
assert sp.gcd(
    c5_section_numerator_poly,
    sp.Poly(c5_q, c5_n),
).degree() == 0
assert sp.factor(
    (50 - c5_section_trace_4) / 16
    - c5_section_numerator / (64 * c5_q**3)
) == 0
c5_section_at_zero = sp.Poly(
    c5_eta_matrix.subs(c5_n, 0).charpoly(T).as_expr(),
    T,
)
assert c5_section_at_zero.as_expr() == T**5 - 5 * T**3 + 4 * T + sp.Rational(7, 5)
c5_compact_integral = sp.Poly(
    5 * T**5 - 25 * T**3 + 20 * T + 7,
    T,
)
assert sp.discriminant(c5_compact_integral.as_expr(), T) == 26875**2
assert c5_compact_integral.count_roots(-sp.oo, sp.oo) == 5
assert factor_degrees_mod_prime(c5_compact_integral, 2) == (5,)
assert factor_degrees_mod_prime(c5_compact_integral, 7) == (
    1,
    1,
    1,
    1,
    1,
)

print("PASS: discriminant pullback and simple branch are exact")
print("PASS: fresh primes outside {2,3,13} have discriminant valuation one")
print("PASS: quintic signature and unramified cycle-type tables are complete")
print("PASS: the rational De Moivre surface has generic F20 certificates")
print("PASS: both F20 signatures and all four unramified types occur")
print("PASS: the F20 surface has a tame fresh-ramification branch")
print("PASS: the rational Brumer curve has generic D5 certificates")
print("PASS: both D5 signatures and all three unramified types occur")
print("PASS: the D5 curve has a tame fresh-ramification branch")
print("PASS: fixed-map A5 anchors realize both signatures and all local types")
print("PASS: both Mestre A5 source-pencil identities are exact")
print("PASS: the total-real affine A5 lift is a genus-eight cubic cover")
print("PASS: the quadratic A5 trace conic has a 5-adic obstruction at infinity")
print("PASS: the smaller A5 descent surface is one explicit cubic cover")
print("PASS: A5 polynomial sections have no degrees zero through six")
print("PASS: the full C5 trace quadric has an explicit birational parametrization")
print("PASS: the displayed Lehmer C5 section lies on a genus-ten cubic cover")
print("PASS: the compact C5 anchor realizes both unramified local types")
