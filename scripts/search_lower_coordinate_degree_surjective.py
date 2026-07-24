#!/usr/bin/env python3
"""Exact first search for a surjective Keller counterexample of degree < 17.

The search compares the three established suspension families in geometric
degree five, applies the exact full-contact test before accepting a
candidate, and checks all single elementary source shears of degree at most
two on the two competitive surjective normal forms.
"""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from jcsearch.weighted import WeightedSeedModel, w, x, y, z  # noqa: E402
from master_cancellation import hensel_jet, parameter_polynomial  # noqa: E402


def profile(mapping: tuple[sp.Expr, ...]) -> tuple[int, ...]:
    return tuple(
        sp.Poly(component, x, y, z).total_degree() for component in mapping
    )


def cancellation_map(m: int, r: int) -> tuple[tuple[sp.Expr, ...], sp.Expr]:
    """Return the normalized cancellation map over QQ[q]/(M_(m,r))."""

    A, q, t = sp.symbols("A q t")
    modulus = parameter_polynomial(m, r, q)
    h = hensel_jet(m, r, A, q)
    actual_A = 1 + x * y**m
    B = actual_A ** (r + 1) * z + y ** (m + 1) * h.subs(A, actual_A)
    P = sp.expand(actual_A * B)
    Q = sp.expand(y + x * B)
    antiderivative = sp.integrate(
        sp.expand((1 - t * (Q - P * t) ** m) ** r), t
    )
    rational_R = sp.cancel(antiderivative.subs(t, x / actual_A))
    numerator, denominator = sp.fraction(rational_R)
    numerator = sp.rem(
        sp.Poly(numerator, q), sp.Poly(modulus, q)
    ).as_expr()
    R = sp.cancel(numerator / denominator)
    assert sp.denom(R) == 1
    return (P, Q, R), modulus


def no_single_quadratic_source_shear_improves(
    mapping: tuple[sp.Expr, ...],
    coefficient_modulus: sp.Expr | None = None,
) -> None:
    """Check elementary v -> v-lambda*m shears, deg(m)<=2.

    When ``coefficient_modulus`` is supplied, coefficients are reduced in
    its number-field quotient before the common exceptional-lambda test.
    """

    lam = sp.symbols("lambda")
    coefficient_variable = (
        next(iter(coefficient_modulus.free_symbols))
        if coefficient_modulus is not None
        else None
    )
    variables = (x, y, z)
    baseline = max(profile(mapping))

    for variable in variables:
        others = tuple(item for item in variables if item != variable)
        monomials = []
        for degree in (1, 2):
            monomials.extend(
                others[0] ** exponent * others[1] ** (degree - exponent)
                for exponent in range(degree + 1)
            )

        for monomial in monomials:
            transformed = tuple(
                sp.Poly(
                    sp.expand(
                        component.subs(
                            {variable: variable - lam * monomial},
                            simultaneous=True,
                        )
                    ),
                    x,
                    y,
                    z,
                )
                for component in mapping
            )
            generic_degree = max(item.total_degree() for item in transformed)
            top_coefficients = []
            for component in transformed:
                for exponents, coefficient in component.terms():
                    if sum(exponents) != generic_degree:
                        continue
                    if coefficient_modulus is not None:
                        coefficient = sp.rem(
                            sp.Poly(coefficient, coefficient_variable),
                            sp.Poly(
                                coefficient_modulus, coefficient_variable
                            ),
                        ).as_expr()
                        top_coefficients.extend(
                            sp.Poly(
                                coefficient, coefficient_variable
                            ).all_coeffs()
                        )
                    else:
                        top_coefficients.append(coefficient)

            common = None
            for coefficient in top_coefficients:
                polynomial = sp.Poly(coefficient, lam)
                common = (
                    polynomial
                    if common is None
                    else sp.gcd(common, polynomial)
                )

            exceptional = (
                sp.solve(common.as_expr(), lam)
                if common is not None and common.degree() > 0
                else ()
            )
            for value in exceptional:
                if value == 0:
                    continue
                specialized = tuple(
                    sp.Poly(
                        sp.expand(component.as_expr().subs(lam, value)),
                        x,
                        y,
                        z,
                    ).total_degree()
                    for component in transformed
                )
                assert max(specialized) >= baseline


# Weighted family: the canonical quintic is surjective by the exact
# full-contact criterion and has a rational two-point collision.
H = w**4 * (1 - w)
weighted = WeightedSeedModel(sp.diff(H, w))
weighted_map = weighted.mapping()
assert profile(weighted_map) == (17, 16, 4)
assert sp.factor(
    sp.Matrix(weighted_map).jacobian((x, y, z)).det()
) == 1

# A quintic without a simple root has type (3,2), with the quintuple-root
# case included by a=b.  The unaffected W^4,W^3,W^2 coefficients of
# W^5-W^4+sW-t make the full-contact ideal the unit ideal.
a, b = sp.symbols("a b")
full_contact = sp.Poly((w - a) ** 3 * (w - b) ** 2, w)
weighted_contact_ideal = (
    full_contact.nth(4) + 1,
    full_contact.nth(3),
    full_contact.nth(2),
)
assert sp.groebner(weighted_contact_ideal, a, b).polys == [
    sp.Poly(1, a, b)
]

weighted_target = (-16, -16, 1)
weighted_points = (
    (-sp.Rational(1, 15), 16, -3880),
    (sp.Rational(1, 32), -30, 30624),
)
for point in weighted_points:
    assert tuple(
        sp.factor(component.subs(dict(zip((x, y, z), point))))
        for component in weighted_map
    ) == weighted_target

no_single_quadratic_source_shear_improves(weighted_map)


# Cancellation family: geometric degree five has exactly two types.
cancellation_31, modulus_31 = cancellation_map(3, 1)
cancellation_12, modulus_12 = cancellation_map(1, 2)
assert profile(cancellation_31) == (13, 10, 24)
assert profile(cancellation_12) == (9, 8, 13)
for mapping, modulus in (
    (cancellation_31, modulus_31),
    (cancellation_12, modulus_12),
):
    q_parameter = next(iter(modulus.free_symbols))
    determinant_error = sp.fraction(
        sp.cancel(
            sp.Matrix(mapping).jacobian((x, y, z)).det() + 1
        )
    )[0]
    assert sp.rem(
        sp.Poly(determinant_error, q_parameter),
        sp.Poly(modulus, q_parameter),
    ).is_zero

# Type (1,2) is the attractive degree-13 false lead.  At this target its
# inverse polynomial is a fifth power, so all five roots lie on D=0 and no
# affine source point reconstructs.
T = sp.symbols("T")
Psi_12 = sp.integrate((1 - T * (2 - T)) ** 2, T) - sp.Rational(1, 5)
assert sp.factor(Psi_12 - (T - 1) ** 5 / 5) == 0
assert sp.factor(1 - T * (2 - T) - (T - 1) ** 2) == 0

# Type (3,1) is surjective: coefficient comparison with every possible
# (3,2) full-contact polynomial gives the unit ideal.  It is nevertheless
# too large in the normalized chart.
P, Q, scale = sp.symbols("P Q scale")
Psi_31 = sp.integrate(1 - T * (Q - P * T) ** 3, T)
contact_31 = sp.expand(scale * (T - a) ** 3 * (T - b) ** 2)
contact_equations_31 = tuple(
    sp.expand(Psi_31 - contact_31).coeff(T, degree)
    for degree in range(1, 6)
)
contact_basis_31 = sp.groebner(
    contact_equations_31, scale, P, Q, b, a, order="lex"
)
assert len(contact_basis_31.polys) == 1
assert contact_basis_31.polys[0].as_expr() == 1
collision_31 = sp.Poly(T + T**5 / 5, T)
assert collision_31.gcd(collision_31.diff()).degree() == 0

no_single_quadratic_source_shear_improves(
    cancellation_31, modulus_31
)


# Quadratic gauge: the sparsest genuine quintic G=S+mu*S^3+lambda*S^5
# still has degree profile (7,32,30), and in fact always has a full-contact
# omitted target.  The relations below give that target over the algebraic
# closure for every nonzero mu and lambda.
S, mu, nu = sp.symbols("S mu nu", nonzero=True)
t_source = 1 + x * y
q_source = t_source**2 * z + y**2 * (1 + 3 * t_source) / mu
quadratic_map = (
    t_source * q_source,
    y
    + 3 * mu * x * q_source
    + 5 * nu * t_source**2 * x**3 * q_source**5,
    x * (5 - 3 * t_source)
    - mu * x**3 * z
    - 3 * nu * (x * q_source) ** 5,
)
assert profile(quadratic_map) == (7, 32, 30)

p = sp.symbols("p")
quadratic_inverse = (
    S + mu * p * S**3 + nu * p**5 * S**5
)
quadratic_full_contact = sp.expand(
    nu * p**5 * (S - a) ** 3 * (S - b) ** 2
)
B_contact = -2 * sp.Poly(quadratic_full_contact, S).nth(2)
C_contact = -2 * sp.Poly(quadratic_full_contact, S).nth(0)
quadratic_inverse -= (B_contact * S**2 + C_contact) / 2
coefficient_field = sp.QQ.frac_field(mu, nu)
quadratic_contact_relations = sp.groebner(
    (
        15 * nu * p**3 - 4 * mu**2,
        mu * p * a**2 + 1,
        2 * b + 3 * a,
    ),
    b,
    a,
    p,
    domain=coefficient_field,
)
contact_difference = sp.Poly(
    sp.expand(quadratic_inverse - quadratic_full_contact), S
)
assert all(
    quadratic_contact_relations.reduce(coefficient)[1] == 0
    for coefficient in contact_difference.all_coeffs()
)


print("PASS weighted quintic: surjective, collision, degree profile (17,16,4)")
print("PASS cancellation (1,2): degree 13 but excluded by quintuple contact")
print("PASS cancellation (3,1): surjective but degree profile (13,10,24)")
print("PASS sparse quadratic quintics: degree 32 and explicit full contact")
print("PASS bounded source-shear search: no degree below 17")
