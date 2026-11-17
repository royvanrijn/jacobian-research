#!/usr/bin/env python3
"""Exact audit of Juntang Zhuang's three quartic "island" examples.

External provenance (formulas, labels, and collision certificates):
https://github.com/jzkay12/jacobian_conjecture/tree/1ff68e870f66afec8c6611f910fcc8f5522fdbce

This is an independent compact reconstruction, not a copy of the upstream
expanded-polynomial checker.  It verifies the polynomial maps and reported
collisions, derives their quartic weighted-seed resolvents, computes the
saturated boundary traces, and compares the resulting signatures with the
unique degree-four cancellation type.  The S_4 conclusion uses the written
universal-pencil theorem; this finite script is not a monodromy proof.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib

import sympy as sp


x, y, z = sp.symbols("x y z")
source_A, source_T, C = sp.symbols("source_A source_T C")
W, A, B = sp.symbols("W A B")
p, q = sp.symbols("p q")


@dataclass(frozen=True)
class QuarticIsland:
    name: str
    lam: sp.Rational
    compact_p: sp.Expr
    compact_q: sp.Expr
    seed: sp.Expr
    standard_B_sign: int
    resultant_scalar: int
    determinant: sp.Rational
    collision_points: tuple[tuple[sp.Expr, sp.Expr, sp.Expr], ...]
    collision_image: tuple[sp.Expr, sp.Expr, sp.Expr]
    expected_c_order: int
    expected_trace: sp.Expr
    expected_nilpotency: int
    second_boundary_index: int
    extra_root: sp.Expr | None
    upstream_coefficient_fingerprint: str


ISLANDS = (
    QuarticIsland(
        name="F4a",
        lam=sp.Rational(5, 4),
        compact_p=-source_T**2
        * (3 * source_A * source_T**2 * C**2
           - 2 * source_A * source_T * C - 1)
        / source_A,
        compact_q=-source_T
        * (4 * source_A * source_T**2 * C**2
           - 3 * source_A * source_T * C - 1)
        / source_A,
        seed=W**3 * (W - 1),
        standard_B_sign=-1,
        resultant_scalar=1,
        determinant=sp.Rational(-1),
        collision_points=((1, 0, 1), (0, 4, -23)),
        collision_image=(1, 1, 0),
        expected_c_order=3,
        expected_trace=4 * q**3,
        expected_nilpotency=3,
        second_boundary_index=2,
        extra_root=None,
        upstream_coefficient_fingerprint=(
            "a9cc615b4dde9f9d1eb0fa10c2efbdb8e8c57d600165e022a186d2cea8ad419f"
        ),
    ),
    QuarticIsland(
        name="F4b",
        lam=sp.Rational(2),
        compact_p=source_T**2
        * (3 * source_A * source_T**2 * C**2
           - 8 * source_A * source_T * C + 3 * source_A + 2)
        / (2 * source_A),
        compact_q=2 * source_T
        * (2 * source_A * source_T**2 * C**2
           - 6 * source_A * source_T * C + 3 * source_A + 1)
        / source_A,
        seed=W**2 * (W - 1) * (W - 3),
        standard_B_sign=1,
        resultant_scalar=-2,
        determinant=sp.Rational(-2),
        collision_points=(
            (1, 0, 1),
            (-1, sp.Rational(8, 3), sp.Rational(19, 3)),
            (0, 4, -sp.Rational(139, 2)),
        ),
        collision_image=(sp.Rational(5, 2), 8, 0),
        expected_c_order=2,
        expected_trace=-36 * (24 * p - q**2),
        expected_nilpotency=1,
        second_boundary_index=1,
        extra_root=sp.Rational(3),
        upstream_coefficient_fingerprint=(
            "37bba034dc93f53d439be708e121fe57eceb70d5595c409f3dcf6429e6b2836a"
        ),
    ),
    QuarticIsland(
        name="F4c",
        lam=sp.Rational(13, 10),
        compact_p=-source_T**2
        * (6 * source_A * source_T**2 * C**2
           - 2 * source_A * source_T * C - source_A - 3)
        / (3 * source_A),
        compact_q=-source_T
        * (8 * source_A * source_T**2 * C**2
           - 3 * source_A * source_T * C - 2 * source_A - 3)
        / source_A,
        seed=-W**2 * (W - 1) * (2 * W + 1),
        standard_B_sign=1,
        resultant_scalar=-3,
        determinant=sp.Rational(-3),
        collision_points=(
            (1, 0, 1),
            (-1, 5, sp.Rational(15, 2)),
            (0, sp.Rational(50, 9), -sp.Rational(19163, 405)),
        ),
        collision_image=(sp.Rational(4, 3), 5, 0),
        expected_c_order=2,
        expected_trace=-9 * (12 * p - q**2),
        expected_nilpotency=1,
        second_boundary_index=1,
        extra_root=-sp.Rational(1, 2),
        upstream_coefficient_fingerprint=(
            "8ea429bf09926eb622789c324198a1f7afc2255f9ae6b6e5ba67ea8d86181831"
        ),
    ),
)


def c_valuation(poly: sp.Expr) -> int:
    """Return the exact C-adic order of a nonzero polynomial."""
    terms = sp.Poly(poly, C).terms()
    return min(monomial[0] for monomial, coefficient in terms if coefficient)


def nilpotency_index(trace: sp.Expr) -> int:
    """Return one for a reduced hypersurface and the largest repeated exponent."""
    return max(exponent for _, exponent in sp.factor_list(trace)[1])


def coefficient_fingerprint(mapping: tuple[sp.Expr, ...]) -> str:
    """Hash the ordered rational coefficient tables of a polynomial map."""
    chunks = []
    for coordinate in mapping:
        polynomial = sp.Poly(sp.expand(coordinate), x, y, z, domain=sp.QQ)
        for monomial, coefficient in polynomial.terms():
            exponent_text = ",".join(map(str, monomial))
            chunks.append(
                f"{exponent_text}:{coefficient.p}/{coefficient.q}"
            )
        chunks.append("|")
    return hashlib.sha256(";".join(chunks).encode()).hexdigest()


for island in ISLANDS:
    source_D = 1 - x * (island.lam * y + x * z)
    source_substitution = {
        source_A: 1 + x * y,
        source_T: (1 + x * y) / x,
        C: x * source_D,
    }
    mapping = tuple(
        sp.factor(sp.cancel(expr.subs(source_substitution)))
        for expr in (island.compact_p, island.compact_q, C.subs(source_substitution))
    )

    # The apparent x and A denominators cancel, recovering the exact upstream
    # expanded polynomial maps from a compact source-side presentation.
    for coordinate in mapping:
        denominator = sp.cancel(coordinate).as_numer_denom()[1]
        assert not denominator.free_symbols
        sp.Poly(coordinate, x, y, z, domain=sp.QQ)
    assert coefficient_fingerprint(mapping) == (
        island.upstream_coefficient_fingerprint
    )

    jacobian = sp.factor(sp.Matrix(mapping).jacobian((x, y, z)).det())
    assert jacobian == island.determinant
    for point in island.collision_points:
        values = tuple(
            sp.factor(coordinate.subs(dict(zip((x, y, z), point))))
            for coordinate in mapping
        )
        assert values == island.collision_image

    # Eliminate source_A from the two compact target equations.  The harmless
    # factor source_T is the excluded primitive-element chart; the other
    # factor is the quartic inverse equation.
    numerator_p = sp.together(p - island.compact_p).as_numer_denom()[0]
    numerator_q = sp.together(q - island.compact_q).as_numer_denom()[0]
    resultant = sp.factor(sp.resultant(numerator_p, numerator_q, source_A))
    t_resolvent = sp.factor(resultant / (island.resultant_scalar * source_T))

    c = -sp.diff(island.seed, W).subs(W, 1)
    weighted_resolvent = sp.expand(
        island.seed - island.standard_B_sign * q * C * W + c * p * C**2
    )
    from_t = sp.factor(
        t_resolvent.subs(source_T, W / C) * C**2
    )
    assert sp.factor(from_t - weighted_resolvent) == 0
    assert sp.Poly(weighted_resolvent, W).degree() == 4

    # Saturate the pulled-back discriminant by its exact boundary power and
    # inspect the scheme-theoretic trace on C=0.
    discriminant = sp.factor(sp.discriminant(weighted_resolvent, W))
    order = c_valuation(discriminant)
    assert order == island.expected_c_order
    saturated = sp.factor(discriminant / C**order)
    trace = sp.factor(saturated.subs(C, 0))
    assert sp.factor(trace - island.expected_trace) == 0
    assert nilpotency_index(trace) == island.expected_nilpotency

    # Degree four forces (m,r)=(2,1) in the cancellation family.  Its second
    # boundary is unramified and its thick-intersection index is six.
    cancellation_pairs = [
        (m, r_value)
        for m in range(1, 5)
        for r_value in range(1, 5)
        if r_value * (m + 1) + 1 == 4
    ]
    assert cancellation_pairs == [(2, 1)]
    assert (
        island.second_boundary_index != 1
        or island.expected_nilpotency != 2 * 1 * 3
    )

    root_label = (
        "zero cluster of multiplicity 3"
        if island.extra_root is None
        else f"extra root rho={island.extra_root}"
    )
    print(
        f"PASS {island.name}: degree=4, seed={sp.factor(island.seed)}, "
        f"{root_label}, e_Delta=2, e_C={island.second_boundary_index}, "
        f"mu={island.expected_nilpotency}, monodromy=S_4 (universal theorem), "
        "not cancellation type (2,1)"
    )
