#!/usr/bin/env python3
r"""Exact full-box tangent audits for the explicit maps F_4, F_5, and F_6.

For d=max_i deg(F_i), the determinant-one coefficient scheme is

    X(3,d) = {Phi in (QQ[x,y,z]_{\le d})^3 : det(D Phi)=1}.

This checker constructs the sparse rational matrix of its linearized
Jacobian operator at the all-degree rational-fiber map F_N,

    L_F(G) = trace(adj(D F) D G),

and computes its rank by exact sparse Gaussian elimination over QQ.  It also
differentiates the N-3-dimensional normalized weighted-seed family and checks
the canonical source trivializers adj(D F) G.

The calculation concerns the full bounded-degree coefficient box, not the
smaller weighted-support coefficient scheme.
"""

from __future__ import annotations

import argparse
import math
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import WeightedSeedModel, w, x, y, z  # noqa: E402


VARIABLES = (x, y, z)
EXPECTED = {
    4: {
        "degree": 12,
        "rank": 1307,
        "nullity": 58,
        "obstruction_exponent": 3,
    },
    5: {
        "degree": 17,
        "rank": 3332,
        "nullity": 88,
        "obstruction_exponent": 3,
    },
    6: {
        "degree": 22,
        "rank": 6777,
        "nullity": 123,
        "obstruction_exponent": 4,
    },
}


def extra_roots(degree: int) -> tuple[int, ...]:
    """The N-2 extra integer roots in the uniform construction."""
    if degree < 3:
        raise ValueError("inverse degree must be at least three")
    k = degree // 2
    if degree % 2:
        return (2,) + tuple(
            root for j in range(3, k + 2) for root in (j, 1 - j)
        )
    return (3, 4) + tuple(
        root for j in range(5, k + 3) for root in (j, 1 - j)
    )


def explicit_seed(degree: int) -> sp.Expr:
    """Return the normalized integer-root seed H_N."""
    roots = (0, -1) + extra_roots(degree)
    polynomial = sp.prod(w - root for root in roots)
    derivative = sp.diff(polynomial, w)
    scale = derivative.subs(w, 0) - derivative.subs(w, 1)
    return sp.cancel(
        (polynomial - derivative.subs(w, 0) * w) / scale
    )


def monomials_through(degree: int) -> list[tuple[int, int, int]]:
    """Three-variable exponent triples of total degree at most degree."""
    return [
        (first, second, total - first - second)
        for total in range(degree + 1)
        for first in range(total + 1)
        for second in range(total - first + 1)
    ]


def rational(value: sp.Expr) -> Fraction:
    """Convert a SymPy rational to the standard-library exact type."""
    value = sp.Rational(value)
    return Fraction(int(value.p), int(value.q))


def adjugate_terms(
    adjugate: sp.Matrix,
) -> list[list[list[tuple[tuple[int, int, int], Fraction]]]]:
    """Sparse monomial terms of every adjugate entry."""
    return [
        [
            [
                (exponent, rational(coefficient))
                for exponent, coefficient in sp.Poly(
                    adjugate[row, column], *VARIABLES
                ).terms()
            ]
            for column in range(3)
        ]
        for row in range(3)
    ]


def linearized_column(
    component: int,
    exponent: tuple[int, int, int],
    terms: list[
        list[list[tuple[tuple[int, int, int], Fraction]]]
    ],
) -> dict[tuple[int, int, int], Fraction]:
    """Sparse coefficient dictionary for L_F(x^exponent e_component)."""
    column: dict[tuple[int, int, int], Fraction] = {}
    for variable in range(3):
        if exponent[variable] == 0:
            continue
        base = list(exponent)
        base[variable] -= 1
        derivative_coefficient = exponent[variable]
        for adjugate_exponent, coefficient in terms[variable][component]:
            output_exponent = tuple(
                base[index] + adjugate_exponent[index]
                for index in range(3)
            )
            new_coefficient = (
                column.get(output_exponent, Fraction())
                + derivative_coefficient * coefficient
            )
            if new_coefficient:
                column[output_exponent] = new_coefficient
            else:
                column.pop(output_exponent, None)
    return column


def leading_exponent(
    column: dict[tuple[int, int, int], Fraction],
) -> tuple[int, int, int]:
    """Graded-lexicographic leading exponent."""
    return max(column, key=lambda exponent: (sum(exponent), exponent))


def exact_sparse_rank(
    adjugate: sp.Matrix, degree: int
) -> tuple[
    int,
    int,
    int,
    dict[tuple[int, int, int], dict[tuple[int, int, int], Fraction]],
]:
    """Rank L_F over QQ, number of nonzero rows, and maximum pivot bit size."""
    terms = adjugate_terms(adjugate)
    pivots: dict[
        tuple[int, int, int], dict[tuple[int, int, int], Fraction]
    ] = {}
    output_exponents: set[tuple[int, int, int]] = set()
    maximum_bit_size = 0

    for component in range(3):
        for exponent in monomials_through(degree):
            column = linearized_column(component, exponent, terms)
            output_exponents.update(column)
            while column:
                lead = leading_exponent(column)
                coefficient = column[lead]
                pivot = pivots.get(lead)
                if pivot is None:
                    normalized = {
                        monomial: value / coefficient
                        for monomial, value in column.items()
                    }
                    pivots[lead] = normalized
                    maximum_bit_size = max(
                        maximum_bit_size,
                        max(
                            max(
                                abs(value.numerator).bit_length(),
                                value.denominator.bit_length(),
                            )
                            for value in normalized.values()
                        ),
                    )
                    break
                for monomial, value in pivot.items():
                    new_value = (
                        column.get(monomial, Fraction())
                        - coefficient * value
                    )
                    if new_value:
                        column[monomial] = new_value
                    else:
                        column.pop(monomial, None)

    return (
        len(pivots),
        len(output_exponents),
        maximum_bit_size,
        pivots,
    )


def reduce_modulo_linearized_image(
    polynomial: sp.Expr,
    pivots: dict[
        tuple[int, int, int], dict[tuple[int, int, int], Fraction]
    ],
) -> dict[tuple[int, int, int], Fraction]:
    """Return the exact echelon remainder modulo im(L_F)."""
    working = {
        exponent: rational(coefficient)
        for exponent, coefficient in sp.Poly(
            polynomial, *VARIABLES
        ).terms()
    }
    remainder: dict[tuple[int, int, int], Fraction] = {}
    while working:
        lead = leading_exponent(working)
        coefficient = working.pop(lead)
        pivot = pivots.get(lead)
        if pivot is None:
            remainder[lead] = coefficient
            continue
        for monomial, value in pivot.items():
            if monomial == lead:
                continue
            new_value = (
                working.get(monomial, Fraction()) - coefficient * value
            )
            if new_value:
                working[monomial] = new_value
            else:
                working.pop(monomial, None)
    return remainder


def mapping_from_primitive(primitive: sp.Expr) -> sp.Matrix:
    """Weighted map formula, allowing rational dependence on a parameter."""
    seed = sp.diff(primitive, w)
    kappa = sp.diff(seed, w).subs(w, 1)
    parameter_a = sp.cancel(-(1 + kappa) / (2 + kappa))
    invariant_v = x * y
    invariant_s = x**2 * z
    invariant_u = 1 + invariant_v
    gamma = 1 + parameter_a * invariant_v + invariant_s
    capital_w = invariant_u * gamma
    lagrangian = sp.expand(w * seed - primitive)
    return sp.Matrix(
        [
            sp.cancel(
                (
                    invariant_u
                    + lagrangian.subs(w, capital_w) / gamma**2
                )
                / x**2
            ),
            sp.cancel((1 + seed.subs(w, capital_w) / gamma) / x),
            sp.expand(x * gamma),
        ]
    )


def polynomial_degree(polynomial: sp.Expr) -> int:
    """Total degree, using -1 for zero."""
    if polynomial == 0:
        return -1
    return sp.Poly(polynomial, *VARIABLES).total_degree()


def term_count(polynomial: sp.Expr) -> int:
    """Expanded term count, using zero for zero."""
    if polynomial == 0:
        return 0
    return len(sp.Poly(polynomial, *VARIABLES).terms())


def visible_seed_directions(
    primitive: sp.Expr, mapping: sp.Matrix, adjugate: sp.Matrix, degree: int
) -> list[dict[str, object]]:
    """Audit the N-3 normalized-seed tangents and source trivializers."""
    parameter = sp.symbols("seed_parameter")
    directions: list[sp.Matrix] = []
    records: list[dict[str, object]] = []

    for index in range(degree - 3):
        seed_tangent = w ** (index + 2) * (w - 1) ** 2
        family = mapping_from_primitive(
            primitive + parameter * seed_tangent
        )
        direction = sp.Matrix(
            [
                sp.cancel(sp.diff(component, parameter).subs(parameter, 0))
                for component in family
            ]
        )
        assert all(
            polynomial_degree(component) <= max(
                polynomial_degree(entry) for entry in mapping
            )
            for component in direction
        )
        trivializer = sp.Matrix(
            [sp.expand(component) for component in adjugate * direction]
        )
        divergence = sp.expand(
            sum(
                sp.diff(trivializer[index], VARIABLES[index])
                for index in range(3)
            )
        )
        assert divergence == 0
        reconstructed = (
            sp.Matrix(mapping).jacobian(VARIABLES) * trivializer - direction
        )
        assert all(sp.expand(component) == 0 for component in reconstructed)
        directions.append(direction)
        records.append(
            {
                "direction_degrees": tuple(
                    polynomial_degree(component) for component in direction
                ),
                "trivializer_degrees": tuple(
                    polynomial_degree(component) for component in trivializer
                ),
                "trivializer_terms": tuple(
                    term_count(component) for component in trivializer
                ),
            }
        )

    coefficient_rows = sorted(
        {
            (component_index, exponent)
            for direction in directions
            for component_index, component in enumerate(direction)
            for exponent, _coefficient in sp.Poly(
                component, *VARIABLES
            ).terms()
        }
    )
    coefficient_matrix = sp.Matrix(
        [
            [
                sp.Poly(
                    direction[component_index], *VARIABLES
                ).coeff_monomial(exponent)
                for direction in directions
            ]
            for component_index, exponent in coefficient_rows
        ]
    )
    assert coefficient_matrix.rank() == degree - 3
    return records


def quadratic_obstruction_witness(
    mapping: sp.Matrix,
    jacobian: sp.Matrix,
    coefficient_degree: int,
    exponent: int,
    pivots: dict[
        tuple[int, int, int], dict[tuple[int, int, int], Fraction]
    ],
) -> dict[str, object]:
    """Certify one nonzero raw bounded-box quadratic Kuranishi class."""
    gamma = sp.cancel(mapping[2] / x)
    scalar = sp.expand(x**exponent * gamma ** (exponent - 1))
    direction = sp.Matrix([sp.expand(y * scalar), scalar, 0])
    direction_degrees = tuple(
        polynomial_degree(component) for component in direction
    )
    assert max(direction_degrees) <= coefficient_degree

    linear_term = sp.expand(
        sp.trace(jacobian.adjugate() * direction.jacobian(VARIABLES))
    )
    assert linear_term == 0

    parameter = sp.symbols("obstruction_parameter")
    determinant = sp.Poly(
        sp.expand(
            (
                jacobian
                + parameter * direction.jacobian(VARIABLES)
            ).det()
        ),
        parameter,
    )
    quadratic_term = sp.factor(
        determinant.coeff_monomial(parameter**2)
    )
    expected_quadratic = sp.expand(
        -x ** (2 * exponent + 2) * gamma ** (2 * exponent - 2)
    )
    assert sp.expand(quadratic_term - expected_quadratic) == 0
    remainder = reduce_modulo_linearized_image(quadratic_term, pivots)
    assert remainder

    return {
        "exponent": exponent,
        "direction_degrees": direction_degrees,
        "quadratic_term": quadratic_term,
        "remainder_terms": len(remainder),
    }


def audit_degree(degree: int) -> None:
    """Run the complete exact audit for one inverse degree."""
    primitive = explicit_seed(degree)
    model = WeightedSeedModel(sp.diff(primitive, w))
    mapping = sp.Matrix(model.mapping())
    jacobian = mapping.jacobian(VARIABLES)
    assert sp.factor(jacobian.det()) == 1

    coordinate_degrees = tuple(
        polynomial_degree(component) for component in mapping
    )
    coefficient_degree = max(coordinate_degrees)
    coefficient_monomials = math.comb(coefficient_degree + 3, 3)
    coefficient_variables = 3 * coefficient_monomials
    equation_slots = math.comb(3 * coefficient_degree, 3)

    rank, nonzero_rows, maximum_bit_size, pivots = exact_sparse_rank(
        jacobian.adjugate(), coefficient_degree
    )
    tangent_dimension = coefficient_variables - rank
    expected = EXPECTED.get(degree)
    if expected is not None:
        assert coefficient_degree == expected["degree"]
        assert rank == expected["rank"]
        assert tangent_dimension == expected["nullity"]

    visible = visible_seed_directions(
        primitive, mapping, jacobian.adjugate(), degree
    )
    maximum_trivializer_degrees = tuple(
        max(record["trivializer_degrees"]) for record in visible
    )
    obstruction = None
    if expected is not None:
        obstruction = quadratic_obstruction_witness(
            mapping,
            jacobian,
            coefficient_degree,
            expected["obstruction_exponent"],
            pivots,
        )

    print(
        f"N={degree}: coordinate degrees={coordinate_degrees}, "
        f"d={coefficient_degree}"
    )
    print(
        f"  X(3,{coefficient_degree}): variables={coefficient_variables}, "
        f"determinant coefficient slots={equation_slots}"
    )
    print(
        f"  L_F: nonzero rows={nonzero_rows}, rank_QQ={rank}, "
        f"tangent dimension={tangent_dimension}, "
        f"maximum pivot bit size={maximum_bit_size}"
    )
    print(
        f"  visible seed rank={degree - 3}, "
        "source-trivializer maximum degrees="
        f"{maximum_trivializer_degrees}"
    )
    if obstruction is not None:
        print(
            "  raw quadratic obstruction: "
            f"m={obstruction['exponent']}, "
            f"direction degrees={obstruction['direction_degrees']}, "
            f"remainder terms={obstruction['remainder_terms']}"
        )
        print(f"    Q(G)={obstruction['quadratic_term']}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "degrees",
        nargs="*",
        type=int,
        default=(4, 5, 6),
        help="inverse degrees to audit (default: 4 5 6)",
    )
    args = parser.parse_args()
    if any(degree < 4 for degree in args.degrees):
        parser.error("this coefficient-tangent audit starts in degree four")

    for degree in args.degrees:
        audit_degree(degree)

    print("PASS: exact full-box tangent ranks certified over QQ")
    print("PASS: all visible normalized-seed tangents are independent")
    print(
        "PASS: every visible tangent has its exact divergence-free "
        "canonical source trivializer"
    )
    print(
        "PASS: explicit raw-box tangents have nonzero quadratic "
        "Kuranishi classes in degrees four through six"
    )


if __name__ == "__main__":
    main()
