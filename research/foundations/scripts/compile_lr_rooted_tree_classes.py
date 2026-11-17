#!/usr/bin/env python3
"""Compile constant-direction rooted-tree defects for the F_2 target lift.

For a rooted tree ``tau`` decorated by the constant target fields
``partial_A, partial_B, partial_C``, compare the target and source elementary
differentials

    Delta_F(tau) = E_source(tau; ell_F(Y_v))
                   - ell_F(E_target(tau; Y_v)).

The one-edge tree is the directed pre-Lie defect, hence the symmetric second
fundamental form from COMPLEXITY_FILTERED_CONTACT.md.  All larger target
elementary differentials vanish for constant decorations, but their lifted
source trees need not vanish in the filtered normal quotient.

The calculation is performed in torus semi-invariant coordinates

    x^p (q_r d_r + q_u d_u + q_gamma d_gamma),
    r = log(x),

over Q[u,gamma].  This avoids expanding the high-degree source polynomials in
x,y,z while retaining exact affine-coordinate elementary differentials.

The experiment certifies an all-order family of nonzero *individual tree
classes*.  It does not assert that these classes survive with nonzero total
coefficient in a BCH/LR forcing sum; that separate cancellation problem is
recorded in the accompanying note.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from pathlib import Path

import sympy as sp

from search_rees_torsion_witnesses import (
    G_1,
    G_2,
    J,
    TargetMonomialField,
    gamma,
    logarithmic_source_coordinates,
    second_fundamental_form,
    u,
)


LABELS = ("A", "B", "C")
LABEL_INDEX = {label: index for index, label in enumerate(LABELS)}
FIELD_WEIGHTS = (2, 1, -1)
COMPONENT_SHIFTS = (1, -1, -2)
PARTIAL_SHIFTS = (-1, 1, 2)
SEPARATOR = sp.Rational(1, 6)

v = u - 1
S = gamma - 1 + sp.Rational(8, 7) * v
J_INVERSE = J.inv().applyfunc(sp.expand)


@dataclass(frozen=True)
class SemiInvariantField:
    """A homogeneous source field in logarithmic invariant coordinates."""

    weight: int
    logarithmic: tuple[sp.Expr, sp.Expr, sp.Expr]

    def vector(self) -> sp.Matrix:
        return sp.Matrix(self.logarithmic)


@dataclass(frozen=True)
class RootedTree:
    """A non-planar rooted tree with a constant target direction at each vertex."""

    label: str
    children: tuple["RootedTree", ...] = ()

    def __post_init__(self) -> None:
        if self.label not in LABEL_INDEX:
            raise ValueError(f"unknown target direction {self.label!r}")
        canonical = tuple(sorted(self.children, key=lambda child: child.encoding()))
        object.__setattr__(self, "children", canonical)

    def order(self) -> int:
        return 1 + sum(child.order() for child in self.children)

    def weight(self) -> int:
        return FIELD_WEIGHTS[LABEL_INDEX[self.label]] + sum(
            child.weight() for child in self.children
        )

    def encoding(self) -> str:
        if not self.children:
            return self.label
        return f"{self.label}({','.join(child.encoding() for child in self.children)})"


def add_fields(
    first: SemiInvariantField, second: SemiInvariantField
) -> SemiInvariantField:
    """Add homogeneous fields of the same torus weight."""
    if first.weight != second.weight:
        raise ValueError("cannot add fields of different torus weights")
    return SemiInvariantField(
        first.weight,
        tuple(
            sp.expand(first.logarithmic[index] + second.logarithmic[index])
            for index in range(3)
        ),
    )


def scale_field(scalar: sp.Expr, field: SemiInvariantField) -> SemiInvariantField:
    """Scale a homogeneous field without changing its torus weight."""
    return SemiInvariantField(
        field.weight,
        tuple(sp.expand(scalar * entry) for entry in field.logarithmic),
    )


def ladder(*labels: str) -> RootedTree:
    """Return the rooted ladder whose labels are listed from root to leaf."""
    if not labels:
        raise ValueError("a ladder needs at least one vertex")
    tree = RootedTree(labels[-1])
    for label in reversed(labels[:-1]):
        tree = RootedTree(label, (tree,))
    return tree


def constant_source_lift(label: str) -> SemiInvariantField:
    index = LABEL_INDEX[label]
    return SemiInvariantField(
        FIELD_WEIGHTS[index],
        tuple(sp.expand(entry) for entry in J_INVERSE[:, index]),
    )


def cartesian_coefficients(field: SemiInvariantField) -> tuple[sp.Expr, ...]:
    """Strip the forced x powers from the Cartesian components of ``field``."""
    q_r, q_u, q_gamma = field.logarithmic
    return (
        sp.expand(q_r),
        sp.expand(q_u - v * q_r),
        sp.expand(q_gamma - 2 * S * q_r + sp.Rational(8, 7) * q_u),
    )


def logarithmic_from_cartesian(
    weight: int, coefficients: tuple[sp.Expr, sp.Expr, sp.Expr]
) -> SemiInvariantField:
    q_r = sp.expand(coefficients[0])
    q_u = sp.expand(coefficients[1] + v * q_r)
    q_gamma = sp.expand(
        coefficients[2] + 2 * S * q_r - sp.Rational(8, 7) * q_u
    )
    return SemiInvariantField(weight, (q_r, q_u, q_gamma))


def partial_coefficient(
    coefficient: sp.Expr, x_exponent: int, coordinate: int
) -> tuple[sp.Expr, int]:
    """Differentiate x^k*f(u,gamma) in x, y, or z, retaining its x power."""
    if coordinate == 0:
        differentiated = (
            x_exponent * coefficient
            + v * sp.diff(coefficient, u)
            + (2 * gamma - 2 + sp.Rational(8, 7) * v)
            * sp.diff(coefficient, gamma)
        )
    elif coordinate == 1:
        differentiated = sp.diff(coefficient, u) - sp.Rational(8, 7) * sp.diff(
            coefficient, gamma
        )
    elif coordinate == 2:
        differentiated = sp.diff(coefficient, gamma)
    else:
        raise ValueError(f"invalid Cartesian coordinate {coordinate}")
    return sp.expand(differentiated), x_exponent + PARTIAL_SHIFTS[coordinate]


def elementary_derivative(
    root: SemiInvariantField, children: tuple[SemiInvariantField, ...]
) -> SemiInvariantField:
    """Compute D^k(root)[children] in the original affine source connection."""
    if not children:
        return root

    child_cartesian = tuple(cartesian_coefficients(child) for child in children)
    total_weight = root.weight + sum(child.weight for child in children)
    result = []
    for component, root_coefficient in enumerate(cartesian_coefficients(root)):
        expected_exponent = total_weight + COMPONENT_SHIFTS[component]
        component_result = sp.Integer(0)
        for choices in product(range(3), repeat=len(children)):
            coefficient = root_coefficient
            exponent = root.weight + COMPONENT_SHIFTS[component]
            for coordinate in choices:
                coefficient, exponent = partial_coefficient(
                    coefficient, exponent, coordinate
                )
            for child, coordinate, cartesian in zip(
                children, choices, child_cartesian, strict=True
            ):
                coefficient *= cartesian[coordinate]
                exponent += child.weight + COMPONENT_SHIFTS[coordinate]
            assert exponent == expected_exponent
            component_result += coefficient
        result.append(sp.expand(component_result))
    return logarithmic_from_cartesian(total_weight, tuple(result))


def prelie(
    first: SemiInvariantField, second: SemiInvariantField
) -> SemiInvariantField:
    """The affine map-composition pre-Lie product ``(D first) second``."""
    return elementary_derivative(first, (second,))


def map_bracket(
    first: SemiInvariantField, second: SemiInvariantField
) -> SemiInvariantField:
    """The map-composition bracket, the antisymmetrization of ``prelie``."""
    return add_fields(prelie(first, second), scale_field(-1, prelie(second, first)))


@lru_cache(maxsize=None)
def source_elementary(tree: RootedTree) -> SemiInvariantField:
    root = constant_source_lift(tree.label)
    children = tuple(source_elementary(child) for child in tree.children)
    return elementary_derivative(root, children)


def target_elementary(tree: RootedTree) -> sp.Matrix:
    """Elementary differential on the target for constant decorations."""
    if tree.children:
        return sp.zeros(3, 1)
    result = sp.zeros(3, 1)
    result[LABEL_INDEX[tree.label]] = 1
    return result


def lift_constant_target_vector(vector: sp.Matrix) -> SemiInvariantField:
    logarithmic = J_INVERSE * vector
    nonzero = [index for index, entry in enumerate(vector) if entry != 0]
    if not nonzero:
        return SemiInvariantField(0, (sp.Integer(0),) * 3)
    weights = {FIELD_WEIGHTS[index] for index in nonzero}
    if len(weights) != 1:
        raise ValueError("target vector is not weight homogeneous")
    return SemiInvariantField(
        weights.pop(), tuple(sp.expand(entry) for entry in logarithmic)
    )


def tree_defect(tree: RootedTree) -> SemiInvariantField:
    source = source_elementary(tree)
    target_lift = lift_constant_target_vector(target_elementary(tree))
    if target_elementary(tree) == sp.zeros(3, 1):
        target_lift = SemiInvariantField(source.weight, (sp.Integer(0),) * 3)
    assert source.weight == target_lift.weight
    return SemiInvariantField(
        source.weight,
        tuple(
            sp.expand(source.logarithmic[index] - target_lift.logarithmic[index])
            for index in range(3)
        ),
    )


def saturated_normal_residue(
    field: SemiInvariantField,
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    if field.weight != 0:
        raise ValueError("the saturated normal residue expects total weight zero")
    target_coordinates = J * field.vector()
    return (
        sp.expand(G_1.reduce(sp.expand(target_coordinates[0]))[1]),
        sp.expand(G_2.reduce(sp.expand(target_coordinates[1]))[1]),
        sp.expand(target_coordinates[2].subs(gamma, 0)),
    )


def third_normal_residue(field: SemiInvariantField) -> sp.Expr:
    if field.weight != 0:
        raise ValueError("the third normal residue expects total weight zero")
    # The third row of J is (gamma,0,1), so modulo gamma this is q_gamma.
    return sp.expand(field.logarithmic[2].subs(gamma, 0))


def root_action_matrix(label: str, child_weight: int = 0) -> sp.Matrix:
    """Matrix for grafting ``label`` above a homogeneous child tree."""
    columns = []
    for index in range(3):
        basis = [sp.Integer(0)] * 3
        basis[index] = 1
        child = SemiInvariantField(child_weight, tuple(basis))
        columns.append(elementary_derivative(constant_source_lift(label), (child,)).vector())
    return sp.Matrix.hstack(*columns).applyfunc(sp.expand)


def weighted_symbol_mod_gamma(residue: sp.Expr) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(residue), u)
    if polynomial.is_zero:
        return sp.Integer(0)
    return polynomial.LC() * u ** polynomial.degree()


def matrix_to_strings(matrix: sp.Matrix) -> list[list[str]]:
    return [
        [sp.sstr(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def separator_sequence(matrix: sp.Matrix, seed: sp.Matrix, terms: int) -> list[sp.Expr]:
    result = []
    vector = seed
    for _ in range(terms):
        result.append(sp.factor(vector[2]))
        vector = matrix * vector
    return result


def family_tree(order: int) -> RootedTree:
    if order < 2:
        raise ValueError("the all-order family starts at order two")
    tree = ladder("B", "C") if order % 2 == 0 else ladder("A", "C", "C")
    for _ in range((order - tree.order()) // 2):
        tree = RootedTree("B", (RootedTree("C", (tree,)),))
    assert tree.order() == order
    assert tree.weight() == 0
    return tree


def compile_certificate(max_order: int) -> dict[str, object]:
    if max_order < 3:
        raise ValueError("max_order must be at least three")

    one_edge = ladder("B", "C")
    compiled_edge = tree_defect(one_edge)
    old_edge = second_fundamental_form(
        TargetMonomialField(LABEL_INDEX["B"], (0, 0, 0)),
        TargetMonomialField(LABEL_INDEX["C"], (0, 0, 0)),
    )
    old_logarithmic = logarithmic_source_coordinates(old_edge)
    assert all(
        sp.expand(compiled_edge.logarithmic[index] - old_logarithmic[index]) == 0
        for index in range(3)
    )

    matrices = {label: root_action_matrix(label) for label in LABELS}
    for label in LABELS:
        # The stripped coefficient matrix is independent of the child's torus
        # weight; only the factored x power changes.
        assert root_action_matrix(label, -2) == matrices[label]
        assert root_action_matrix(label, 3) == matrices[label]
    transfer = (matrices["B"] * matrices["C"]).applyfunc(sp.expand)
    transfer_mod_gamma = transfer.subs(gamma, 0).applyfunc(sp.expand)
    numeric_transfer = transfer_mod_gamma.subs(u, SEPARATOR)

    even_seed = tree_defect(ladder("B", "C")).vector().subs(gamma, 0)
    odd_seed = tree_defect(ladder("A", "C", "C")).vector().subs(gamma, 0)
    numeric_even_seed = even_seed.subs(u, SEPARATOR)
    numeric_odd_seed = odd_seed.subs(u, SEPARATOR)

    lambda_symbol = sp.Symbol("lambda")
    characteristic = sp.Poly(
        numeric_transfer.charpoly(lambda_symbol).as_expr(), lambda_symbol
    )
    coefficients = characteristic.all_coeffs()
    assert coefficients[0] == 1
    recurrence = tuple(sp.factor(-coefficient) for coefficient in coefficients[1:])
    assert all(coefficient > 0 for coefficient in recurrence)

    even_separator = separator_sequence(numeric_transfer, numeric_even_seed, 4)
    odd_separator = separator_sequence(numeric_transfer, numeric_odd_seed, 3)
    assert even_separator[0] > 0
    assert all(value < 0 for value in even_separator[1:])
    assert all(value < 0 for value in odd_separator)

    predicted = {
        2: tree_defect(ladder("B", "C")).vector(),
        3: tree_defect(ladder("A", "C", "C")).vector(),
    }
    rows = []
    for order in range(2, max_order + 1):
        tree = family_tree(order)
        defect = tree_defect(tree)
        if order >= 4:
            predicted[order] = (transfer * predicted[order - 2]).applyfunc(sp.expand)
        assert defect.vector() == predicted[order]
        residue = third_normal_residue(defect)
        polynomial = sp.Poly(residue, u)
        assert not polynomial.is_zero
        rows.append(
            {
                "order": order,
                "tree": tree.encoding(),
                "third_residue_degree_u": polynomial.degree(),
                "third_residue_terms": len(polynomial.terms()),
                "third_residue_leading_coefficient": sp.sstr(polynomial.LC()),
                "third_residue_at_u_1_over_6": sp.sstr(
                    sp.factor(residue.subs(u, SEPARATOR))
                ),
                "associated_graded_symbol": sp.sstr(
                    weighted_symbol_mod_gamma(residue)
                ),
            }
        )

    second_residue = saturated_normal_residue(compiled_edge)
    expected_third = -sp.Rational(30, 7) * (
        4896 * u**5
        - 25092 * u**4
        + 15232 * u**3
        - 1887 * u**2
        + 126 * u
        - 21
    )
    assert sp.expand(second_residue[2] - expected_third) == 0

    return {
        "description": "Constant-direction rooted-tree normal classes for the F_2 target lift",
        "scope": {
            "decorations": ["partial_A", "partial_B", "partial_C"],
            "target_weights": {"partial_A": 2, "partial_B": 1, "partial_C": -1},
            "normal_quotient_third_summand": "Q[u,gamma]/(gamma)",
            "filtration": "deg(u)=2, deg(gamma)=3",
        },
        "compiler_cross_checks": {
            "one_vertex_defect_zero": all(
                tree_defect(RootedTree(label)).vector() == sp.zeros(3, 1)
                for label in LABELS
            ),
            "one_edge_BC_equals_second_fundamental_form": True,
            "one_edge_third_residue": sp.sstr(second_residue[2]),
            "one_edge_associated_graded_symbol": sp.sstr(
                weighted_symbol_mod_gamma(second_residue[2])
            ),
        },
        "all_order_family": {
            "even_seed": "B(C)",
            "odd_seed": "A(C(C))",
            "recursion": "tau_(n+2)=B(C(tau_n))",
            "orders_covered": "every n>=2",
            "transfer_matrix_at_gamma_0_u_1_over_6": matrix_to_strings(
                numeric_transfer
            ),
            "even_seed_vector_at_gamma_0_u_1_over_6": [
                sp.sstr(entry) for entry in numeric_even_seed
            ],
            "odd_seed_vector_at_gamma_0_u_1_over_6": [
                sp.sstr(entry) for entry in numeric_odd_seed
            ],
            "characteristic_polynomial": sp.sstr(characteristic.as_expr()),
            "cayley_hamilton_recurrence_coefficients": [
                sp.sstr(coefficient) for coefficient in recurrence
            ],
            "even_initial_third_values_k_0_through_3": [
                sp.sstr(value) for value in even_separator
            ],
            "odd_initial_third_values_k_0_through_2": [
                sp.sstr(value) for value in odd_separator
            ],
            "proof": (
                "The recurrence coefficients are positive. The odd initial "
                "values are negative, and the even values from k=1 onward "
                "start with three negative values; Cayley-Hamilton therefore "
                "makes every later separator value negative. The exceptional "
                "even k=0 value is positive. Hence every third residue is a "
                "nonzero polynomial, so its weighted associated-graded symbol "
                "in gr(Q[u,gamma]/(gamma)) is nonzero."
            ),
        },
        "computed_orders": rows,
        "status_boundary": (
            "This proves all-order nonvanishing for the displayed individual "
            "rooted-tree normal classes. It does not prove that their total "
            "coefficient survives cancellations in the mixed BCH forcing "
            "series or uniformly over all lower LR jets."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-order", type=int, default=12)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "artifacts/generated-results/lr_rooted_tree_normal_classes.json"
        ),
    )
    parser.add_argument("--no-write", action="store_true")
    arguments = parser.parse_args()

    certificate = compile_certificate(arguments.max_order)
    payload = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    if not arguments.no_write:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(payload, encoding="utf-8")

    family = certificate["all_order_family"]
    print("PASS: rooted-tree compiler reproduces the quadratic pre-Lie defect")
    print(
        "PASS: tau_2=B(C), tau_3=A(C(C)), "
        "tau_(n+2)=B(C(tau_n)) have nonzero third normal residue for every n>=2"
    )
    print(
        "PASS: Cayley-Hamilton recurrence coefficients = "
        + ", ".join(family["cayley_hamilton_recurrence_coefficients"])
    )
    print(
        "LIMIT: individual tree nonvanishing does not exclude cancellation "
        "in a BCH/LR forcing sum"
    )
    print(f"certificate_sha256={digest}")


if __name__ == "__main__":
    main()
