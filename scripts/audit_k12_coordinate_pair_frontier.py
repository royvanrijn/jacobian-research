#!/usr/bin/env python3
"""Exact bounded audit of literal coordinate-pair reductions of K12.

For each component K_j that is triangular in z_j, restrict to the source
graph K_j=0 and ask whether every component of degree greater than three can
be repaired by one target shear

    y_i -> y_i - P_i(y_1,...,hat(y_i),...,hat(y_j),...,y_12).

The target polynomial P_i is bounded by degree three.  The two closest raw
restrictions, j=11 and j=12, are also checked through target degree four.
The script additionally classifies every linear target coordinate whose
pullback is a graph coordinate and proves by unit ideals that none gives a
raw degree-three restriction.  Inconsistency of the shear systems is
certified by a rank increase after reduction modulo a good prime.
"""

from __future__ import annotations

import hashlib
import json
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Iterator

import sympy as sp

from audit_macfarlane_g20_dimension_reduction import build_maps


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "k12_coordinate_pair_frontier.json"
)
PRIME = 1_000_003
DESIRED_DEGREE = 3

Exponent = tuple[int, ...]
SparsePolynomial = dict[Exponent, int]


def residue(value: sp.Expr) -> int:
    rational = sp.Rational(value)
    return (
        int(rational.p) % PRIME
    ) * pow(int(rational.q) % PRIME, PRIME - 2, PRIME) % PRIME


def as_sparse_mod(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> SparsePolynomial:
    return {
        exponents: coefficient
        for exponents, value in sp.Poly(
            expression, *variables, domain=sp.QQ
        ).terms()
        if (coefficient := residue(value))
    }


def multiply(
    left: SparsePolynomial,
    right: SparsePolynomial,
) -> SparsePolynomial:
    result: SparsePolynomial = {}
    for left_exponents, left_coefficient in left.items():
        for right_exponents, right_coefficient in right.items():
            exponents = tuple(
                left_value + right_value
                for left_value, right_value in zip(
                    left_exponents, right_exponents
                )
            )
            value = (
                result.get(exponents, 0)
                + left_coefficient * right_coefficient
            ) % PRIME
            if value:
                result[exponents] = value
            else:
                result.pop(exponents, None)
    return result


def high_degree_part(polynomial: SparsePolynomial) -> SparsePolynomial:
    return {
        exponents: coefficient
        for exponents, coefficient in polynomial.items()
        if sum(exponents) > DESIRED_DEGREE
    }


class SparseColumnSpace:
    """Incremental normalized column echelon form over F_PRIME."""

    def __init__(self) -> None:
        self.pivots: dict[Exponent, SparsePolynomial] = {}

    @property
    def rank(self) -> int:
        return len(self.pivots)

    def add(self, source: SparsePolynomial) -> bool:
        column = dict(source)
        while column:
            pivot = min(column)
            coefficient = column[pivot]
            if pivot not in self.pivots:
                inverse = pow(coefficient, PRIME - 2, PRIME)
                normalized = {
                    exponents: value * inverse % PRIME
                    for exponents, value in column.items()
                }
                self.pivots[pivot] = normalized
                return True
            existing = self.pivots[pivot]
            for exponents, value in existing.items():
                updated = (
                    column.get(exponents, 0) - coefficient * value
                ) % PRIME
                if updated:
                    column[exponents] = updated
                else:
                    column.pop(exponents, None)
        return False


def target_monomials(
    polynomials: list[SparsePolynomial],
    maximum_degree: int,
) -> Iterator[tuple[int, SparsePolynomial]]:
    """Yield every target monomial of degrees 1,...,maximum_degree."""

    one = {(0,) * len(next(iter(polynomials[0]))): 1}

    def extend(
        degree: int,
        start: int,
        current: SparsePolynomial,
    ) -> Iterator[SparsePolynomial]:
        if degree == 0:
            yield current
            return
        for index in range(start, len(polynomials)):
            yield from extend(
                degree - 1,
                index,
                multiply(current, polynomials[index]),
            )

    for degree in range(1, maximum_degree + 1):
        for polynomial in extend(degree, 0, one):
            yield degree, polynomial


def screen_component(
    restricted: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    component_index: int,
    maximum_target_degree: int,
) -> dict[str, object]:
    available = [
        as_sparse_mod(component, variables)
        for index, component in enumerate(restricted)
        if index != component_index
    ]
    space = SparseColumnSpace()
    basis_count = 0
    nonzero_high_columns = 0
    basis_by_degree: dict[int, int] = {}
    for degree, polynomial in target_monomials(
        available, maximum_target_degree
    ):
        basis_count += 1
        basis_by_degree[degree] = basis_by_degree.get(degree, 0) + 1
        column = high_degree_part(polynomial)
        if column:
            nonzero_high_columns += 1
        space.add(column)
    rank_before = space.rank
    target = high_degree_part(
        as_sparse_mod(restricted[component_index], variables)
    )
    target_adds_pivot = space.add(target)
    return {
        "component": component_index + 1,
        "maximum_target_degree": maximum_target_degree,
        "target_variable_count": len(available),
        "basis_count": basis_count,
        "basis_count_by_degree": {
            str(degree): basis_by_degree[degree]
            for degree in sorted(basis_by_degree)
        },
        "nonzero_high_degree_columns": nonzero_high_columns,
        "rank_mod_prime": rank_before,
        "augmented_rank_mod_prime": space.rank,
        "target_adds_pivot": target_adds_pivot,
    }


def build_k12() -> tuple[tuple[sp.Symbol, ...], list[sp.Expr]]:
    data = build_maps()
    x = data["x"]
    f13 = data["F13"]
    z = tuple(sp.symbols("z1:13"))
    s = sp.Symbol("s")
    inverse_source = {x[index]: z[index] for index in range(12)}
    inverse_source[x[12]] = s - z[1] ** 2
    transformed = [
        sp.expand(component.subs(inverse_source, simultaneous=True))
        for component in f13
    ]
    transformed[3] = sp.expand(transformed[3] - transformed[7] ** 2)
    assert transformed[12] == s
    return z, [
        sp.expand(component.subs(s, 0))
        for component in transformed[:12]
    ]


def graph_restriction(
    k12: list[sp.Expr],
    z: tuple[sp.Symbol, ...],
    omitted: int,
) -> tuple[tuple[sp.Symbol, ...], sp.Expr, list[sp.Expr]]:
    correction = sp.expand(k12[omitted] - z[omitted])
    assert z[omitted] not in correction.free_symbols
    graph = sp.expand(-correction)
    retained_variables = tuple(
        variable for index, variable in enumerate(z) if index != omitted
    )
    restricted = [
        sp.expand(component.subs(z[omitted], graph))
        for index, component in enumerate(k12)
        if index != omitted
    ]
    return retained_variables, graph, restricted


def audit_linear_graph_coordinates(
    k12: list[sp.Expr],
    z: tuple[sp.Symbol, ...],
) -> list[dict[str, object]]:
    """Classify linear g whose pullback g(K) is a graph coordinate."""

    nonlinear = [
        sp.expand(component - variable)
        for component, variable in zip(k12, z)
    ]
    nonlinear_polynomials = [
        sp.Poly(component, *z, domain=sp.QQ)
        for component in nonlinear
    ]
    monomials = sorted(
        set().union(
            *[
                {
                    exponents
                    for exponents, coefficient in polynomial.terms()
                    if coefficient
                }
                for polynomial in nonlinear_polynomials
            ]
        )
    )
    collision_image = (
        sp.Integer(0),
        sp.Integer(0),
        -sp.Rational(1, 4),
        *([sp.Integer(0)] * 9),
    )
    records: list[dict[str, object]] = []
    possible_pivots: list[int] = []
    for pivot_index in range(12):
        constraint_rows = [
            [
                polynomial.coeff_monomial(exponents)
                for polynomial in nonlinear_polynomials
            ]
            for exponents in monomials
            if exponents[pivot_index]
        ]
        constraint_matrix = sp.Matrix(constraint_rows)
        kernel = constraint_matrix.nullspace()
        pivot_possible = any(vector[pivot_index] for vector in kernel)
        record: dict[str, object] = {
            "source_pivot_variable": pivot_index + 1,
            "constraint_rank": constraint_matrix.rank(),
            "linear_coordinate_space_dimension": len(kernel),
            "nonzero_pivot_possible": pivot_possible,
        }
        if not pivot_possible:
            records.append(record)
            continue

        possible_pivots.append(pivot_index + 1)
        pivot_row = [sp.Integer(0)] * 12
        pivot_row[pivot_index] = sp.Integer(1)
        system = sp.Matrix.vstack(
            constraint_matrix, sp.Matrix([pivot_row])
        )
        right_hand_side = sp.Matrix(
            [sp.Integer(0)] * constraint_matrix.rows + [sp.Integer(1)]
        )
        normalized = next(
            iter(sp.linsolve((system, right_hand_side)))
        )
        parameters = sorted(
            set().union(
                *[coefficient.free_symbols for coefficient in normalized]
            ),
            key=str,
        )
        slice_value = sum(
            normalized[index] * collision_image[index]
            for index in range(12)
        )
        correction = sp.expand(
            sum(
                normalized[index] * nonlinear[index]
                for index in range(12)
            )
        )
        graph = sp.expand(
            slice_value
            - sum(
                normalized[index] * z[index]
                for index in range(12)
                if index != pivot_index
            )
            - correction
        )
        retained_variables = tuple(
            variable
            for index, variable in enumerate(z)
            if index != pivot_index
        )
        high_degree_equations: list[sp.Expr] = []
        for index, component in enumerate(k12):
            if index == pivot_index:
                continue
            restricted = sp.Poly(
                sp.expand(component.subs(z[pivot_index], graph)),
                *retained_variables,
                domain=sp.QQ.frac_field(*parameters),
            )
            high_degree_equations.extend(
                sp.factor(coefficient)
                for exponents, coefficient in restricted.terms()
                if sum(exponents) > DESIRED_DEGREE
            )
        high_degree_equations = list(dict.fromkeys(high_degree_equations))
        groebner = sp.groebner(
            high_degree_equations, *parameters, order="grevlex"
        )
        unit_ideal = groebner.contains(sp.Integer(1))
        assert unit_ideal
        record.update(
            {
                "normalized_parameter_count": len(parameters),
                "normalized_coefficients": [
                    str(coefficient) for coefficient in normalized
                ],
                "raw_degree_three_equation_count": len(
                    high_degree_equations
                ),
                "raw_degree_three_ideal_is_unit": unit_ideal,
                "reduced_groebner_basis": [
                    str(polynomial.as_expr())
                    for polynomial in groebner.polys
                ],
            }
        )
        records.append(record)
    assert possible_pivots == list(range(4, 13))
    return records


def main() -> None:
    z, k12 = build_k12()
    linear_graph_coordinates = audit_linear_graph_coordinates(k12, z)
    triangular = [
        index
        for index in range(12)
        if z[index] not in (k12[index] - z[index]).free_symbols
    ]
    assert triangular == list(range(3, 12))

    deletions: list[dict[str, object]] = []
    closest: dict[int, tuple[tuple[sp.Symbol, ...], list[sp.Expr], int]] = {}
    expected_bad = {
        4: [3, 10],
        5: [2, 3, 7, 9],
        6: [2, 3, 9],
        7: [3],
        8: [2, 3, 4],
        9: [2],
        10: [3],
        11: [1],
        12: [1],
    }

    for omitted in triangular:
        variables, graph, restricted = graph_restriction(k12, z, omitted)
        polynomials = [
            sp.Poly(component, *variables, domain=sp.QQ)
            for component in restricted
        ]
        component_degrees = [
            polynomial.total_degree() for polynomial in polynomials
        ]
        bad = [
            index
            for index, degree in enumerate(component_degrees)
            if degree > DESIRED_DEGREE
        ]
        bad_original_indices = [
            index + 1 if index < omitted else index + 2
            for index in bad
        ]
        assert bad_original_indices == expected_bad[omitted + 1]
        high_terms = [
            (
                index + 1 if index < omitted else index + 2,
                sum(
                    1
                    for exponents, coefficient in polynomial.terms()
                    if coefficient and sum(exponents) > DESIRED_DEGREE
                ),
                sp.expand(
                    sum(
                        coefficient
                        * sp.prod(
                            variable**exponent
                            for variable, exponent in zip(
                                variables, exponents
                            )
                        )
                        for exponents, coefficient in polynomial.terms()
                        if coefficient
                        and sum(exponents) > DESIRED_DEGREE
                    )
                ),
            )
            for index, polynomial in enumerate(polynomials)
            if component_degrees[index] > DESIRED_DEGREE
        ]
        screens = []
        for index in bad:
            screen = screen_component(restricted, variables, index, 3)
            screen["component"] = (
                index + 1 if index < omitted else index + 2
            )
            screens.append(screen)
        assert any(screen["target_adds_pivot"] for screen in screens)
        if omitted in (10, 11):
            assert len(bad) == 1
            closest[omitted] = (variables, restricted, bad[0])
        deletions.append(
            {
                "omitted_component": omitted + 1,
                "source_graph": f"z{omitted + 1}={graph}",
                "source_graph_degree": sp.Poly(
                    graph, *variables, domain=sp.QQ
                ).total_degree(),
                "restricted_component_degrees": component_degrees,
                "restricted_maximum_degree": max(component_degrees),
                "bad_original_components": bad_original_indices,
                "high_degree_term_counts": [
                    {
                        "component": index,
                        "count": count,
                        "defect": str(defect),
                    }
                    for index, count, defect in high_terms
                ],
                "degree_three_target_screens": screens,
            }
        )

    quartic_screens = []
    for omitted in (10, 11):
        variables, restricted, bad = closest[omitted]
        screen = screen_component(restricted, variables, bad, 4)
        assert screen["basis_count"] == 1_000
        # Ten linear columns have no high-degree terms in these one-bad-output
        # cases; the 990 nonlinear columns are independent modulo PRIME.
        assert screen["nonzero_high_degree_columns"] == 990
        assert screen["rank_mod_prime"] == 990
        assert screen["target_adds_pivot"]
        quartic_screens.append(
            {"omitted_component": omitted + 1, **screen}
        )

    assert all(
        any(
            screen["target_adds_pivot"]
            for screen in deletion["degree_three_target_screens"]
        )
        for deletion in deletions
    )
    artifact = {
        "format": "k12-coordinate-pair-frontier-v1",
        "status": "exact bounded obstruction over Q",
        "prime": PRIME,
        "desired_restricted_degree": DESIRED_DEGREE,
        "literal_triangular_components": [index + 1 for index in triangular],
        "nontriangular_components": [
            index + 1 for index in range(12) if index not in triangular
        ],
        "linear_target_coordinate_graph_audit": linear_graph_coordinates,
        "deletions": deletions,
        "quartic_target_screens_for_closest_cases": quartic_screens,
        "conclusion": (
            "No linear target coordinate g has a raw degree-at-most-three "
            "graph restriction g(K)=g(p). No literal triangular component "
            "K_j, j=4,...,12, gives such a restriction after one parallel "
            "target-shear stage of polynomial degree at most three. For the "
            "closest literal deletions j=11,12, the same statement holds "
            "through target degree four."
        ),
        "proof_logic": (
            "For at least one bad retained component in every deletion, the "
            "high-degree defect is outside the span of all high-degree parts "
            "of bounded target monomials in the other raw retained outputs. "
            "The modular augmented-rank increase exhibits a nonzero minor "
            "over F_1000003; the corresponding rational minor is nonzero, "
            "so the target-completion linear system is inconsistent over Q."
        ),
        "scope": (
            "This is not a lower bound for Keller counterexamples and does "
            "not exclude nonlinear source coordinates, nonliteral linear "
            "combinations of output coordinates, or multi-stage target "
            "automorphisms."
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    digest = hashlib.sha256(OUTPUT.read_bytes()).hexdigest()
    print("PASS K12: literal triangular components are exactly 4,...,12")
    print(
        "PASS K11: every linear target-coordinate graph family has a "
        "unit raw-degree-three obstruction ideal"
    )
    print(
        "PASS K11: all nine raw graph deletions exceed degree three"
    )
    print(
        "PASS K11: every degree-three one-stage completion has an exact "
        "modular-rank obstruction"
    )
    print(
        "PASS K11: closest deletions 11 and 12 remain obstructed through "
        "target degree four"
    )
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")
    print(f"SHA256 {digest}")


if __name__ == "__main__":
    main()
