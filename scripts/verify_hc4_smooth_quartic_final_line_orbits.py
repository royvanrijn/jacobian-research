#!/usr/bin/env python3
"""Verify the two remaining Fermat-symmetry line-fiber normal forms.

The four roots of

    tau^4 - 4*tau^3 + 10*tau^2 - 4*tau + 1

meet both remaining six-point Fermat-symmetry orbits.  On this quartic the
visible polar conic is the union of a resultant line and a residual-polar
line.  For each component this checker constructs the complete 81-equation
reciprocal-Hessian packet, solves its ten linear equations over the exact
function field, and checks that the reduced equations span the eight
remaining deformation variables.  The last check is coefficient-field
linear algebra, not a sampled or finite-field calculation.
"""

from __future__ import annotations

import argparse
import hashlib
from dataclasses import dataclass
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import PolyElement, ring

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "research_hc4_smooth_quartic_simple_line.py"
DRIVER_SHA256 = "99f3f94c9ee4bac0a489f25916ff290b076d33e7165e88b0a952754548c419ec"


@dataclass(frozen=True)
class ExactField:
    tau: sp.Symbol
    parameters: tuple[sp.Symbol, ...]
    field: object
    tau_element: object
    parameter_elements: dict[sp.Symbol, object]


def make_exact_field() -> ExactField:
    tau = sp.symbols("tau")
    quartic = tau**4 - 4 * tau**3 + 10 * tau**2 - 4 * tau + 1
    assert sp.factor(quartic) == quartic
    alpha = sp.CRootOf(quartic, 0)
    number_field = QQ.algebraic_field(alpha, alias="tau")
    parameters = sp.symbols("m c sigma b15 b16 b17")
    function_field = number_field.frac_field(*parameters)
    tau_element = function_field.new(number_field.from_sympy(alpha))
    parameter_elements = dict(zip(parameters, function_field.gens, strict=True))
    return ExactField(
        tau=tau,
        parameters=parameters,
        field=function_field,
        tau_element=tau_element,
        parameter_elements=parameter_elements,
    )


def to_field(expression: sp.Expr, exact: ExactField):
    """Map a polynomial in tau and the six parameters to the exact field."""

    polynomial = sp.Poly(
        sp.expand(expression), exact.tau, *exact.parameters, domain=QQ
    )
    value = exact.field.zero
    for exponents, coefficient in polynomial.terms():
        term = exact.field.convert(coefficient)
        term *= exact.tau_element ** exponents[0]
        for parameter, exponent in zip(
            exact.parameters, exponents[1:], strict=True
        ):
            term *= exact.parameter_elements[parameter] ** exponent
        value += term
    return value


def line_substitution(component: str, exact: ExactField) -> dict[sp.Symbol, sp.Expr]:
    tau = exact.tau
    m, c, *_ = exact.parameters
    p, q, r = sp.symbols("p q r")
    if component == "resultant":
        return {p: 3 * c, q: 3 * c * m, r: -c * (3 + m)}
    if component == "residual-polar":
        coefficient_p = 6 * tau**3 - 30 * tau**2 + 15 * tau - 4
        coefficient_q = -4 * tau**3 + 10 * tau**2 - 5 * tau + 1
        coefficient_r = 6 * tau**3 + 1
        return {
            p: coefficient_r * c,
            q: coefficient_r * c * m,
            r: -c * (coefficient_p + coefficient_q * m),
        }
    raise ValueError(component)


def solve_linear_layer(
    equations: list[sp.Expr],
    unknowns: tuple[sp.Symbol, ...],
    exact: ExactField,
) -> tuple[tuple[int, ...], tuple[int, ...], list[list[object]]]:
    from research_hc4_smooth_quartic_simple_line import unknown_degree

    linear_equations = [
        equation for equation in equations if unknown_degree(equation, unknowns) <= 1
    ]
    assert len(linear_equations) == 10
    matrix, right_hand_side = sp.linear_eq_to_matrix(linear_equations, unknowns)
    rows = [
        [to_field(matrix[row, column], exact) for column in range(len(unknowns))]
        + [to_field(right_hand_side[row, 0], exact)]
        for row in range(matrix.rows)
    ]
    reduced, pivots = DomainMatrix(
        rows, (len(rows), len(unknowns) + 1), exact.field
    ).rref(method="CD")
    assert len(pivots) == 10
    assert all(pivot < len(unknowns) for pivot in pivots)
    free = tuple(index for index in range(len(unknowns)) if index not in pivots)
    assert len(free) == 8
    return tuple(pivots), free, reduced.to_list()


def quotient_representatives(
    unknowns: tuple[sp.Symbol, ...],
    pivots: tuple[int, ...],
    free: tuple[int, ...],
    reduced_linear: list[list[object]],
    exact: ExactField,
):
    names = ",".join(str(unknowns[index]) for index in free)
    ring_data = ring(names, exact.field, order="grevlex")
    polynomial_ring, free_generators = ring_data[0], ring_data[1:]
    representatives: list[PolyElement | None] = [None] * len(unknowns)
    for index, generator in zip(free, free_generators, strict=True):
        representatives[index] = generator
    for row, pivot in enumerate(pivots):
        representative = polynomial_ring.ground_new(reduced_linear[row][-1])
        for index, generator in zip(free, free_generators, strict=True):
            representative -= reduced_linear[row][index] * generator
        representatives[pivot] = representative
    assert all(representative is not None for representative in representatives)
    return polynomial_ring, free_generators, tuple(representatives)


def reduce_equation(
    equation: sp.Expr,
    unknowns: tuple[sp.Symbol, ...],
    representatives: tuple[PolyElement, ...],
    polynomial_ring,
    exact: ExactField,
) -> PolyElement:
    polynomial = sp.Poly(sp.expand(equation), *unknowns, domain="EX")
    result = polynomial_ring.zero
    for exponents, coefficient in polynomial.terms():
        term = polynomial_ring.ground_new(to_field(coefficient, exact))
        for representative, exponent in zip(
            representatives, exponents, strict=True
        ):
            if exponent:
                term *= representative**exponent
        result += term
    return result


def row_space_contains_variables(
    polynomials: list[PolyElement], free_generators, exact: ExactField
) -> tuple[int, int, tuple[bool, ...], tuple[tuple[int, int, int], ...]]:
    monomials = {
        monomial
        for polynomial in polynomials
        for monomial, _ in polynomial.terms()
    }
    target_monomials = [generator.terms()[0][0] for generator in free_generators]
    monomials.update(target_monomials)
    ordered_monomials = sorted(monomials, reverse=True)
    column = {monomial: index for index, monomial in enumerate(ordered_monomials)}
    rows: list[list[object]] = []
    for polynomial in polynomials:
        row = [exact.field.zero] * len(ordered_monomials)
        for monomial, coefficient in polynomial.terms():
            row[column[monomial]] = coefficient
        rows.append(row)
    matrix = DomainMatrix(rows, (len(rows), len(ordered_monomials)), exact.field)
    _, denominator, pivots = matrix.rref_den(method="CD", keep_domain=False)
    factor_profiles: list[tuple[int, int, int]] = []
    for factor, exponent in denominator.factor_list()[1]:
        monomials_in_factor = [monomial for monomial, _ in factor.terms()]
        m_degree = max(monomial[0] for monomial in monomials_in_factor)
        c_degree = max(monomial[1] for monomial in monomials_in_factor)
        assert all(
            all(value == 0 for value in monomial[2:])
            for monomial in monomials_in_factor
        )
        factor_profiles.append((m_degree, c_degree, exponent))
    contained = (len(pivots) == len(ordered_monomials),) * len(target_monomials)
    return (
        len(ordered_monomials),
        len(pivots),
        tuple(contained),
        tuple(sorted(factor_profiles)),
    )


def verify_polar_factorization(exact: ExactField) -> None:
    tau = exact.tau
    p, q, r = sp.symbols("p q r")
    quartic = tau**4 - 4 * tau**3 + 10 * tau**2 - 4 * tau + 1
    pivot = (
        (3 * p**2 - q * r) * tau**5
        + (9 * p * r - q**2) * tau**4
        + (18 * r**2 - 6 * p * q) * tau**3
        + (18 * p**2 - 6 * q * r) * tau**2
        + (9 * p * r - q**2) * tau
        + (3 * r**2 - p * q)
    )
    resultant = 3 * p + q + 3 * r
    residual_polar = (
        (6 * tau**3 - 30 * tau**2 + 15 * tau - 4) * p
        + (-4 * tau**3 + 10 * tau**2 - 5 * tau + 1) * q
        + (6 * tau**3 + 1) * r
    )
    remainder = sp.rem(
        sp.Poly(pivot - resultant * residual_polar, tau), sp.Poly(quartic, tau)
    )
    assert remainder.as_expr() == 0

    fixed_root_evaluation = (
        tau**5 - 3 * tau**4 + 6 * tau**3 + 6 * tau**2 - 3 * tau + 1
    )
    assert sp.expand(fixed_root_evaluation - (tau + 1) * quartic) == 0

    # The two quadratic factors over Q(i) are the two six-point quotient
    # values.  Their product is the irreducible rational quartic used below.
    assert sp.expand(
        (tau**2 - (2 + 2 * sp.I) * tau + 1)
        * (tau**2 - (2 - 2 * sp.I) * tau + 1)
        - quartic
    ) == 0
    assert sp.expand((2 + 2 * sp.I) ** 3 - 3 * (2 + 2 * sp.I)) == (
        -22 + 10 * sp.I
    )
    assert sp.expand((2 - 2 * sp.I) ** 3 - 3 * (2 - 2 * sp.I)) == (
        -22 - 10 * sp.I
    )
    print("PASS: the quartic polar conic factors into the two declared lines")


def verify_zero_boundary_support() -> None:
    x, y, z = sp.symbols("x y z")
    p, q, r = sp.symbols("p q r")
    b15, b16, b17 = sp.symbols("b15 b16 b17")
    matrix = sp.Matrix(
        [
            [0, 0, -y**2],
            [0, 0, x**2],
            [
                -y**2,
                x**2,
                p * x**2
                + q * x * y
                + r * y**2
                + z * (b15 * x + b16 * y + b17 * z),
            ],
        ]
    )
    assert matrix.det() == 0


def verify_component(component: str, exact: ExactField) -> None:
    from research_hc4_smooth_quartic_simple_line import build_equations

    equations, unknowns, _ = build_equations("squarefree-line", False)
    assert len(equations) == 81
    substitution = line_substitution(component, exact)
    specialized = [
        sp.expand(equation.subs(substitution, simultaneous=True))
        for equation in equations
    ]
    pivots, free, reduced_linear = solve_linear_layer(specialized, unknowns, exact)
    polynomial_ring, free_generators, representatives = quotient_representatives(
        unknowns, pivots, free, reduced_linear, exact
    )
    # These eleven rows are a compact exact witness extracted from the full
    # packet (indices are one-based in the mathematical note).  After the
    # linear quotient they involve exactly eleven monomials.  Nonsingularity
    # of that coefficient matrix is enough; reducing the other sixty rows is
    # mathematically unnecessary and much slower.
    witness_indices = (28, 34, 46, 49, 52, 56, 62, 64, 67, 71, 76)
    reduced_equations = [
        reduce_equation(
            specialized[index - 1],
            unknowns,
            representatives,
            polynomial_ring,
            exact,
        )
        for index in witness_indices
    ]
    reduced_equations = [equation for equation in reduced_equations if equation]
    assert len(reduced_equations) == 11
    monomial_count, rank, contained, factor_profiles = row_space_contains_variables(
        reduced_equations, free_generators, exact
    )
    assert monomial_count == 11
    assert rank == 11
    assert all(contained)
    expected_profiles = {
        "resultant": ((0, 1, 10), (1, 0, 17), (3, 0, 1)),
        "residual-polar": ((0, 1, 10), (2, 0, 17), (7, 0, 1)),
    }
    assert factor_profiles == expected_profiles[component]
    pivot_names = tuple(str(unknowns[index]) for index in pivots)
    free_names = tuple(str(unknowns[index]) for index in free)
    print(
        f"PASS: {component}: linear pivots={pivot_names}, free={free_names}, "
        f"witness_rows={witness_indices}, monomials={monomial_count}, rank={rank}, "
        f"determinant_profiles={factor_profiles}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--audit-existing-only",
        action="store_true",
        help="verify committed equation-builder provenance without exact-field replay",
    )
    arguments = parser.parse_args()

    digest = hashlib.sha256(DRIVER.read_bytes()).hexdigest()
    assert digest == DRIVER_SHA256, (digest, DRIVER_SHA256)
    if arguments.audit_existing_only:
        print(
            "PASS committed HC4 smooth-quartic final-line provenance is intact; "
            "no symbolic or exact-field replay"
        )
        return
    exact = make_exact_field()
    verify_polar_factorization(exact)
    verify_zero_boundary_support()
    for component in ("resultant", "residual-polar"):
        verify_component(component, exact)
    print(
        "THEOREM: both remaining six-point line-fiber normal forms have "
        "zero generic deformation on both polar-line components"
    )


if __name__ == "__main__":
    main()
