#!/usr/bin/env python3
"""Exact fifth obstruction for the explicit nonradial SIC2C4 four-jet."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import sympy as sp
from flint import fmpq_mat


ROOT = Path(__file__).resolve().parents[1]
FOURTH_SCRIPT = ROOT / "scripts" / "verify_two_pair_counterexample_fourth_order.py"
LOCAL_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_counterexample_local_moduli.json"
)
FOURTH_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_counterexample_fourth_order.json"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_counterexample_fifth_order.json"
)
MAX_ORDER = 67
TAIL_COUNT = MAX_ORDER - 11

spec = importlib.util.spec_from_file_location("sic2_fourth", FOURTH_SCRIPT)
assert spec and spec.loader
fourth = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fourth)

local = fourth.local
seed = fourth.seed
BASIS = fourth.BASIS


def parse_vector(values: dict[str, str]) -> sp.Matrix:
    by_exponent = {
        tuple(map(int, label.split("^"))): sp.sympify(value)
        for label, value in values.items()
    }
    return sp.Matrix([by_exponent.get(exponent, 0) for exponent in BASIS])


def vector_polynomial(vector: sp.Matrix) -> seed.Polynomial:
    return {
        BASIS[index]: sp.factor(value)
        for index, value in enumerate(vector)
        if value
    }


def contraction(
    polynomial: seed.Polynomial,
    degree: int,
    f_order: int,
) -> sp.Expr:
    weights = fourth.generating_weights(degree, f_order)
    return sp.factor(
        sum(
            coefficient * weights.get(exponent, 0)
            for exponent, coefficient in polynomial.items()
        )
    )


def master_row(order: int) -> list[str]:
    blocks: list[sp.Expr] = []
    for degree, f_order, scale in (
        (4, order, sp.Rational(1)),
        (8, order - 1, sp.Rational(order)),
        (12, order - 2, sp.Rational(order * (order - 1), 2)),
        (
            16,
            order - 3,
            sp.Rational(order * (order - 1) * (order - 2), 6),
        ),
        (
            20,
            order - 4,
            sp.Rational(
                order * (order - 1) * (order - 2) * (order - 3),
                120,
            ),
        ),
    ):
        weights = fourth.generating_weights(degree, f_order)
        blocks.extend(
            scale
            * weights.get(
                (xi1, degree - xi1, z1, degree - z1),
                0,
            )
            for xi1 in range(degree + 1)
            for z1 in range(degree + 1)
        )
    return [str(value) for value in blocks]


def main() -> None:
    local_data = json.loads(LOCAL_ARTIFACT.read_text())
    fourth_data = json.loads(FOURTH_ARTIFACT.read_text())
    jet = fourth_data[
        "generic_direction_full_cubic_freedom"
    ]["explicit_component_point"]["jet"]
    h, k, ell, m = (parse_vector(jet[label]) for label in ("H", "K", "L", "M"))

    tangent_basis = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in local_data["tangent_basis_columns"]
        ]
    )
    rows = [
        [
            fourth.generating_weights(4, order).get(exponent, 0)
            for exponent in BASIS
        ]
        for order in range(MAX_ORDER + 1)
    ]
    first_rows = sp.Matrix(rows[:12])
    pivot_columns = list(first_rows.rref()[1])
    pivot_inverse = first_rows[:, pivot_columns].inv()
    row_coordinates = [
        sp.Matrix([[row[index] for index in pivot_columns]]) * pivot_inverse
        for row in rows
    ]

    h_poly, k_poly, ell_poly, m_poly = (
        vector_polynomial(vector) for vector in (h, k, ell, m)
    )
    hk2 = seed.multiply(h_poly, seed.multiply(k_poly, k_poly))
    h2ell = seed.multiply(seed.multiply(h_poly, h_poly), ell_poly)
    h3k = seed.multiply(seed.power(h_poly, 3), k_poly)
    h5 = seed.power(h_poly, 5)
    hm = seed.multiply(h_poly, m_poly)
    kell = seed.multiply(k_poly, ell_poly)

    known = []
    tangent_effect = []
    for order in range(MAX_ORDER + 1):
        value = sp.Rational(0)
        effect = sp.zeros(1, 13)
        if order >= 1:
            value += order * (
                contraction(hm, 8, order - 1)
                + contraction(kell, 8, order - 1)
            )
            weights8 = fourth.generating_weights(8, order - 1)
            bilinear = sp.zeros(25)
            for left, left_exponent in enumerate(BASIS):
                for right, right_exponent in enumerate(BASIS):
                    exponent = tuple(
                        left_exponent[index] + right_exponent[index]
                        for index in range(4)
                    )
                    bilinear[left, right] = weights8.get(exponent, 0)
            effect = order * h.T * bilinear * tangent_basis
        if order >= 2:
            value += sp.Rational(order * (order - 1), 2) * (
                contraction(h2ell, 12, order - 2)
                + contraction(hk2, 12, order - 2)
            )
        if order >= 3:
            value += (
                sp.Rational(order * (order - 1) * (order - 2), 6)
                * contraction(h3k, 16, order - 3)
            )
        if order >= 4:
            value += (
                sp.Rational(
                    order * (order - 1) * (order - 2) * (order - 3),
                    120,
                )
                * contraction(h5, 20, order - 4)
            )
        known.append(sp.factor(value))
        tangent_effect.append(effect)

    tail_constant = []
    tail_matrix_rows = []
    for order in range(12, MAX_ORDER + 1):
        constant = known[order]
        effect = tangent_effect[order].copy()
        for prefix in range(12):
            constant -= row_coordinates[order][prefix] * known[prefix]
            effect -= row_coordinates[order][prefix] * tangent_effect[prefix]
        tail_constant.append(sp.factor(constant))
        tail_matrix_rows.append([sp.factor(value) for value in effect])
    tail_matrix = sp.Matrix(tail_matrix_rows)
    tail_constant_vector = sp.Matrix(tail_constant)

    coefficient_rank = tail_matrix.rank()
    augmented_rank = tail_matrix.row_join(-tail_constant_vector).rank()

    master_rows = [master_row(order) for order in range(12, 81)]
    master_basis_rank = fmpq_mat(master_rows[:56]).rank()
    master_replay_rank = fmpq_mat(master_rows).rank()
    assert master_basis_rank == master_replay_rank == 56

    result = {
        "format": "two-pair-counterexample-fifth-order-v1",
        "field": "Q(sqrt(41))",
        "source_fourth_artifact": str(FOURTH_ARTIFACT.relative_to(ROOT)),
        "combined_all_order_tail": {
            "degrees": [4, 8, 12, 16, 20],
            "symbolic_beta_numerator_degree_bound": 55,
            "basis_orders": list(range(12, 68)),
            "basis_rank": master_basis_rank,
            "replay_orders": list(range(12, 81)),
            "replay_rank": master_replay_rank,
        },
        "free_fourth_tangent_parameters": 13,
        "fifth_compatibility_equations": TAIL_COUNT,
        "coefficient_rank": coefficient_rank,
        "augmented_rank": augmented_rank,
        "consistent": coefficient_rank == augmented_rank,
    }

    if coefficient_rank == augmented_rank:
        solution = sp.linsolve(
            (tail_matrix, -tail_constant_vector),
            *sp.symbols("v0:13"),
        )
        point = next(iter(solution))
        zero_point = [
            sp.factor(value.subs({symbol: 0 for symbol in set().union(
                *(entry.free_symbols for entry in point)
            ) if str(symbol).startswith("tau")}))
            for value in point
        ]
        result["one_fourth_tangent_solution"] = [
            str(value) for value in zero_point
        ]
    else:
        left_kernel = tail_matrix.T.nullspace()
        certificate = next(
            vector
            for vector in left_kernel
            if sp.factor((vector.T * tail_constant_vector)[0]) != 0
        )
        pairing = sp.factor((certificate.T * tail_constant_vector)[0])
        conjugate_pairing = sp.factor(
            pairing.xreplace({sp.sqrt(41): -sp.sqrt(41)})
        )
        pairing_norm = sp.factor(sp.expand(pairing * conjugate_pairing))
        assert pairing_norm != 0
        assert not pairing_norm.has(sp.sqrt(41))
        common_denominator = sp.ilcm(
            *[
                sp.denom(value)
                for value in certificate
                if value
            ]
        )
        result["obstruction_certificate"] = {
            "support": {
                str(12 + index): str(sp.factor(value))
                for index, value in enumerate(certificate)
                if value
            },
            "primitive_integer_support": {
                str(12 + index): str(sp.factor(common_denominator * value))
                for index, value in enumerate(certificate)
                if value
            },
            "pairing": str(pairing),
            "conjugate_pairing": str(conjugate_pairing),
            "pairing_norm": str(pairing_norm),
        }

    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(
        "PASS SIC2C4 fifth: combined all-order tail rank "
        f"{master_basis_rank}"
    )
    print(
        "PASS SIC2C4 fifth: tangent/augmented ranks "
        f"{coefficient_rank}/{augmented_rank}"
    )
    print(f"PASS SIC2C4 fifth: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
