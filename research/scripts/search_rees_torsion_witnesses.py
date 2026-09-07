#!/usr/bin/env python3
"""Search low target degrees for opposite-weight quadratic Rees witnesses.

For the degree-five base map F_2, enumerate monomial target vector fields
Y_p,Y_-p.  Compute the second fundamental form

    II_F(Y_p,Y_-p) = -(DF)^(-1) D^2F[ell_F(Y_p),ell_F(Y_-p)]

and test its weight-zero class in the invariant-ring-saturated normal quotient
from TORUS_FILTERED_LR_MODULE.md.  The quotient is deliberately enlarged, so
a nonzero residue is a certificate that survives the actual target image.

This is a witness search, not a proof of strictness: zero residues in the
enlarged quotient do not prove that the corresponding normal symbol vanishes.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import product

import sympy as sp


x, y, z = sp.symbols("x y z")
A, B, C = sp.symbols("A B C")
u, gamma = sp.symbols("u gamma")
v = u - 1
S = gamma - 1 + sp.Rational(8, 7) * v
W = u * gamma

SOURCE_VARIABLES = (x, y, z)
TARGET_VARIABLES = (A, B, C)
SOURCE_WEIGHTS = (1, -1, -2)
TARGET_WEIGHTS = (-2, -1, 1)


def H(root: sp.Expr) -> sp.Expr:
    return sp.expand(root**2 * (root - 1) * (3 * root**2 - 11 * root + 6) / 60)


H_2 = H(W)
p_2 = sp.diff(H(sp.Symbol("root")), sp.Symbol("root")).subs(
    sp.Symbol("root"), W
)
c_2 = sp.Rational(1, 30)
q_2 = sp.cancel((W * p_2 - H_2) / c_2)
a = sp.expand(sp.cancel(u + q_2 / gamma**2))
b = sp.expand(sp.cancel(c_2 + p_2 / gamma))

# Polynomial form of F_2.  The Laurent-looking invariant presentation cancels.
substitute_invariants = {u: 1 + x * y, gamma: 1 - sp.Rational(8, 7) * x * y + x**2 * z}
F = sp.Matrix(
    [
        sp.cancel(x**-2 * a.subs(substitute_invariants)),
        sp.cancel(x**-1 * b.subs(substitute_invariants)),
        sp.cancel(x * gamma.subs(substitute_invariants)),
    ]
)
assert all(sp.denom(entry) == 1 for entry in F)
F = sp.Matrix([sp.expand(entry) for entry in F])
DF = F.jacobian(SOURCE_VARIABLES)
assert sp.factor(DF.det()) == c_2
DF_inverse = DF.adjugate().applyfunc(sp.expand) / c_2

# Logarithmic differential and the invariant-ring-saturated target quotient.
J = sp.Matrix(
    [
        [-2 * a, sp.diff(a, u), sp.diff(a, gamma)],
        [-b, sp.diff(b, u), sp.diff(b, gamma)],
        [gamma, 0, 1],
    ]
)
assert sp.factor(J.det()) == c_2
G_1 = sp.groebner([a, b**2], gamma, u, order="grevlex")
G_2 = sp.groebner([b, a * gamma], gamma, u, order="grevlex")


@dataclass(frozen=True)
class TargetMonomialField:
    component: int
    exponents: tuple[int, int, int]

    @property
    def degree(self) -> int:
        return sum(self.exponents)

    @property
    def monomial_weight(self) -> int:
        return sum(e * w for e, w in zip(self.exponents, TARGET_WEIGHTS))

    @property
    def field_weight(self) -> int:
        return self.monomial_weight - TARGET_WEIGHTS[self.component]

    def vector(self) -> sp.Matrix:
        monomial = sp.prod(variable**exponent for variable, exponent in zip(TARGET_VARIABLES, self.exponents))
        entries = [sp.Integer(0)] * 3
        entries[self.component] = monomial
        return sp.Matrix(entries)

    def label(self) -> str:
        monomial = sp.prod(variable**exponent for variable, exponent in zip(TARGET_VARIABLES, self.exponents))
        return f"({sp.sstr(monomial)}) e_{TARGET_VARIABLES[self.component]}"


def fields_through_degree(maximum_degree: int) -> list[TargetMonomialField]:
    fields = []
    for total_degree in range(maximum_degree + 1):
        for exponents in product(range(total_degree + 1), repeat=3):
            if sum(exponents) != total_degree:
                continue
            for component in range(3):
                fields.append(TargetMonomialField(component, exponents))
    return fields


LIFT_CACHE: dict[TargetMonomialField, sp.Matrix] = {}


def target_lift(field: TargetMonomialField) -> sp.Matrix:
    if field in LIFT_CACHE:
        return LIFT_CACHE[field]
    pulled_back = field.vector().subs(dict(zip(TARGET_VARIABLES, F)))
    result = DF_inverse * pulled_back
    lift = sp.Matrix([sp.expand(entry) for entry in result])
    LIFT_CACHE[field] = lift
    return lift


def second_fundamental_form(
    first: TargetMonomialField, second: TargetMonomialField
) -> sp.Matrix:
    first_lift = target_lift(first)
    second_lift = target_lift(second)
    entries = []
    for coordinate in F:
        hessian = sp.hessian(coordinate, SOURCE_VARIABLES)
        entries.append(sp.expand((first_lift.T * hessian * second_lift)[0]))
    result = -DF_inverse * sp.Matrix(entries)
    return sp.Matrix([sp.expand(entry) for entry in result])


def prelie_defect(
    first: TargetMonomialField, second: TargetMonomialField
) -> sp.Matrix:
    """Independent definition of II_F used to check the Hessian formula."""
    first_vector = first.vector()
    second_vector = second.vector()
    first_lift = target_lift(first)
    second_lift = target_lift(second)
    target_product = (
        first_vector.jacobian(TARGET_VARIABLES) * second_vector
        + second_vector.jacobian(TARGET_VARIABLES) * first_vector
    )
    pulled_back_product = target_product.subs(dict(zip(TARGET_VARIABLES, F)))
    lifted_product = DF_inverse * pulled_back_product
    result = (
        first_lift.jacobian(SOURCE_VARIABLES) * second_lift
        + second_lift.jacobian(SOURCE_VARIABLES) * first_lift
        - lifted_product
    ) / 2
    return sp.Matrix([sp.expand(entry) for entry in result])


def invariant_polynomial(expression: sp.Expr) -> sp.Expr | None:
    """Rewrite a source-torus invariant polynomial in R=Q[v,S]."""
    polynomial = sp.Poly(sp.expand(expression), x, y, z)
    invariant = sp.Integer(0)
    for (x_power, y_power, z_power), coefficient in polynomial.terms():
        # A weight-zero monomial equals v^y_power S^z_power.
        if x_power != y_power + 2 * z_power:
            return None
        invariant += coefficient * v**y_power * S**z_power
    return sp.expand(invariant)


def logarithmic_source_coordinates(field: sp.Matrix) -> sp.Matrix | None:
    # Use the polynomial logarithmic coordinates (delta x/x,delta u,delta gamma).
    # This includes the non-diagonal semi-invariant generators x*z in the y
    # component and y^2 in the z component.
    logarithmic_source = [
        sp.cancel(field[0] / x),
        sp.expand(y * field[0] + x * field[1]),
        sp.expand(
            (-sp.Rational(8, 7) * y + 2 * x * z) * field[0]
            - sp.Rational(8, 7) * x * field[1]
            + x**2 * field[2]
        ),
    ]
    if sp.denom(logarithmic_source[0]) != 1:
        return None
    source = sp.Matrix([invariant_polynomial(entry) for entry in logarithmic_source])
    if any(entry is None for entry in source):
        return None
    return source


def third_normal_residue(field: sp.Matrix) -> sp.Expr | None:
    """Certificate in the third saturated normal summand R/(gamma)."""
    source = logarithmic_source_coordinates(field)
    if source is None:
        return None
    return sp.expand((gamma * source[0] + source[2]).subs(gamma, 0))


def saturated_normal_residue(field: sp.Matrix) -> tuple[sp.Expr, sp.Expr, sp.Expr] | None:
    source = logarithmic_source_coordinates(field)
    if source is None:
        return None
    target = J * source
    return (
        sp.expand(G_1.reduce(sp.expand(target[0]))[1]),
        sp.expand(G_2.reduce(sp.expand(target[1]))[1]),
        sp.expand(target[2]).subs(gamma, 0),
    )


def source_degree(field: sp.Matrix) -> int:
    return max(sp.Poly(entry, x, y, z).total_degree() for entry in field if entry != 0)


def normal_symbol_mod_gamma(residue: sp.Expr) -> sp.Expr:
    """Weighted-leading symbol in gr(R/(gamma)); deg(u)=2."""
    polynomial = sp.Poly(sp.expand(residue), u)
    return polynomial.LC() * u ** polynomial.degree()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-target-degree", type=int, default=2)
    parser.add_argument("--all", action="store_true", help="print all nonzero witnesses")
    parser.add_argument("--verbose", action="store_true", help="print the full II_F vector")
    arguments = parser.parse_args()

    fields = fields_through_degree(arguments.max_target_degree)
    pairs = []
    for first in fields:
        if first.field_weight == 0:
            continue
        for second in fields:
            if first.field_weight != -second.field_weight:
                continue
            # Remove the ordering symmetry.
            if (first.degree, first.component, first.exponents) > (
                second.degree,
                second.component,
                second.exponents,
            ):
                continue
            pairs.append((first, second))
    pairs.sort(
        key=lambda pair: (
            max(pair[0].degree, pair[1].degree),
            pair[0].degree + pair[1].degree,
            abs(pair[0].field_weight),
            pair[0].degree,
            pair[0].component,
            pair[0].exponents,
            pair[1].component,
            pair[1].exponents,
        )
    )

    witnesses = []
    unsupported = 0
    for index, (first, second) in enumerate(pairs, start=1):
        form = second_fundamental_form(first, second)
        third_residue = third_normal_residue(form)
        if third_residue is None:
            unsupported += 1
            continue
        residue = None
        if third_residue == 0:
            residue = saturated_normal_residue(form)
        if third_residue != 0 or (residue is not None and any(entry != 0 for entry in residue)):
            witnesses.append((first, second, form, third_residue))
            first_lift_degree = source_degree(target_lift(first))
            second_lift_degree = source_degree(target_lift(second))
            print(
                f"WITNESS {len(witnesses)}: target degrees ({first.degree},{second.degree}), "
                f"weights ({first.field_weight},{second.field_weight}), "
                f"lift degrees ({first_lift_degree},{second_lift_degree}), "
                f"II source degree {source_degree(form)}"
            )
            print(f"  Y_p    = {first.label()}")
            print(f"  Y_-p   = {second.label()}")
            if third_residue != 0:
                print(f"  R/(gamma) residue = {sp.factor(third_residue)}")
                print(
                    "  normal symbol       = "
                    f"{normal_symbol_mod_gamma(third_residue)}"
                )
            else:
                print(f"  saturated residue  = {tuple(map(sp.sstr, residue))}")
            if arguments.verbose:
                print(f"  II_F = {tuple(map(sp.sstr, form))}")
            independent = prelie_defect(first, second)
            assert all(sp.expand(independent[i] - form[i]) == 0 for i in range(3))
            print("  cross-check         = Hessian and symmetric pre-Lie formulas agree")
            if not arguments.all:
                break

    print(f"searched_pairs={len(pairs)} supported_pairs={len(pairs)-unsupported}")
    if not witnesses:
        print("NO WITNESS in the supported diagonal weight-zero sector")


if __name__ == "__main__":
    main()
