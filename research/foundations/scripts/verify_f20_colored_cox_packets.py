#!/usr/bin/env python3
"""Audit the first genuinely colored Cox packet architecture for F20.

The calculation has four deliberately separate levels.

1. Six primitive packet divisors break the six proportional-row witnesses.
2. The complete positive derivative support has sixteen rational packets.
   Their indicator columns give the Hilbert basis of the packet-supported
   nonnegative divisor monoid and contain the derivative target uniquely.
3. The packet basis compresses to three different-factor Cartier candidates
   with ``div(P_X)=3*D_d+D_q+D_r`` and unique model ``(3,1,1)``.
4. Orders and the known q-conductor unit lattice pass, but no residue cocycle
   for the individual Cox sections is supplied.  Thus conductor descent is
   still uncertified.
5. Since there is no conductor-certified Cox algebra, entrywise adjugate and
   affine-space tests are not run.  An independent conditional class-lattice
   screen records why the compressed packet architecture could not itself be
   affine space under the usual factorial-core hypotheses.

This is a certificate for a scoped combinatorial frontier, not a construction
of a global Cox ring or a polynomial Keller map.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

from boundary_package_compiler import (
    colored_proportionality_witnesses,
    determinant,
    f20_toroidal_ledger_datum,
    matrix_rank,
    smith_invariant_factors,
)


PACKETS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("d_ramification", ("d_ramified_4",)),
    (
        "q_crossing",
        ("q_collision_plus", "q_collision_minus"),
    ),
    (
        "q_node_slopes",
        tuple(f"q_node_slope_{index}" for index in range(1, 5)),
    ),
    (
        "r_ramification",
        ("r_ramified_2_plus", "r_ramified_2_minus"),
    ),
    ("r_cusp_E1", ("r_cusp_E1_total_5",)),
    ("r_cusp_E2", ("r_cusp_E2_total_5",)),
    ("r_cusp_E3_unramified", ("r_cusp_E3_unramified",)),
    (
        "r_cusp_E3_ramification",
        ("r_cusp_E3_ramified_plus", "r_cusp_E3_ramified_minus"),
    ),
    (
        "r_cusp_E4",
        tuple(f"r_cusp_E4_sheet_{index}" for index in range(1, 6)),
    ),
    (
        "triple_E1_ramification",
        (
            "triple_plus_E1_ramified_4",
            "triple_minus_E1_ramified_4",
        ),
    ),
    (
        "triple_E2_cluster",
        tuple(
            f"triple_{sign}_E2_cluster_{index}"
            for sign in ("plus", "minus")
            for index in range(1, 5)
        ),
    ),
    (
        "qr_E1_A_ramification",
        tuple(
            f"qr_tangent_{center}_E1_A_ramified_2"
            for center in range(1, 4)
        ),
    ),
    (
        "qr_E1_A_unramified",
        tuple(
            f"qr_tangent_{center}_E1_A_unramified"
            for center in range(1, 4)
        ),
    ),
    (
        "qr_E1_B_ramification",
        tuple(
            f"qr_tangent_{center}_E1_B_ramified_2"
            for center in range(1, 4)
        ),
    ),
    (
        "qr_E2_A",
        tuple(
            f"qr_tangent_{center}_E2_A_sheet_{index}"
            for center in range(1, 4)
            for index in range(1, 4)
        ),
    ),
    (
        "qr_E2_B",
        tuple(
            f"qr_tangent_{center}_E2_B_sheet_{index}"
            for center in range(1, 4)
            for index in range(1, 3)
        ),
    ),
)


PRIMARY_BREAKERS = (
    "d_ramification",
    "q_crossing",
    "r_ramification",
    "triple_E1_ramification",
    "triple_E2_cluster",
    "qr_E1_A_ramification",
)


# Coefficients of three compact Galois-stable Cartier candidates in the
# primitive packet basis.  They allocate the local d^3*q^2*r^2 different:
# the q-node belongs to D_q, the r-cusp to D_r, the first triple exceptional
# has allocation 3+2+2, the second triple exceptional belongs to D_d, and
# the q-r residual discriminants split as q^2*r on A and r on B.
DIFFERENT_FACTOR_PACKETS: dict[str, dict[str, int]] = {
    "D_d": {
        "d_ramification": 1,
        "triple_E1_ramification": 1,
        "triple_E2_cluster": 1,
    },
    "D_q": {
        "q_crossing": 1,
        "q_node_slopes": 1,
        "triple_E1_ramification": 2,
        "qr_E1_A_ramification": 1,
        "qr_E1_A_unramified": 1,
        "qr_E2_A": 1,
    },
    "D_r": {
        "r_ramification": 1,
        "r_cusp_E1": 4,
        "r_cusp_E2": 8,
        "r_cusp_E3_unramified": 2,
        "r_cusp_E3_ramification": 4,
        "r_cusp_E4": 4,
        "triple_E1_ramification": 2,
        "qr_E1_A_ramification": 1,
        "qr_E1_B_ramification": 1,
        "qr_E2_A": 1,
        "qr_E2_B": 1,
    },
}


def transpose_columns(
    columns: tuple[tuple[int, ...], ...],
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(column[row] for column in columns)
        for row in range(len(columns[0]))
    )


def exact_f20_polynomial() -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    s, t, X = sp.symbols("s t X")
    d = s**2 + 4
    q = 4 * s**2 * t**2 + 4 * s**2 * t + 8 * s * t + 6 * s - 8 * t - 5
    P = sp.expand(
        X**5
        + (t**2 * d - 2 * s - sp.Rational(17, 4)) * X**4
        + (3 * t * d + d + sp.Rational(13, 2) * s + 1) * X**3
        - (t * d + sp.Rational(11, 2) * s - 8) * X**2
        + (s - 6) * X
        + 1
    )
    return P, d, q, X


def build_certificate() -> dict[str, object]:
    datum = f20_toroidal_ledger_datum()
    color_names = tuple(color.name for color in datum.boundary_colors)
    color_index = {name: index for index, name in enumerate(color_names)}
    target = tuple(row[3] for row in datum.valuation_matrix)

    packet_names = tuple(name for name, _colors in PACKETS)
    packet_colors = dict(PACKETS)
    packet_columns = tuple(
        tuple(int(color in colors) for color in color_names)
        for _name, colors in PACKETS
    )
    packet_matrix = transpose_columns(packet_columns)

    flattened = tuple(color for _name, colors in PACKETS for color in colors)
    positive_colors = tuple(
        color_names[index] for index, order in enumerate(target) if order > 0
    )
    assert len(color_names) == 63
    assert len(PACKETS) == 16
    assert len(flattened) == len(set(flattened))
    assert set(flattened) == set(positive_colors)
    assert matrix_rank(packet_matrix) == 16

    packet_orders: dict[str, int] = {}
    for name, colors in PACKETS:
        orders = {target[color_index[color]] for color in colors}
        assert len(orders) == 1
        order = orders.pop()
        assert order > 0
        packet_orders[name] = order

    reconstructed_target = tuple(
        sum(
            packet_orders[name] * packet_columns[column][row]
            for column, name in enumerate(packet_names)
        )
        for row in range(len(color_names))
    )
    assert reconstructed_target == target

    different_names = tuple(DIFFERENT_FACTOR_PACKETS)
    different_columns = tuple(
        tuple(
            sum(
                coefficient
                * packet_columns[packet_names.index(packet)][row]
                for packet, coefficient in packet_coefficients.items()
            )
            for row in range(len(color_names))
        )
        for packet_coefficients in DIFFERENT_FACTOR_PACKETS.values()
    )
    different_matrix = transpose_columns(different_columns)
    assert matrix_rank(different_matrix) == 3
    assert tuple(
        3 * different_columns[0][row]
        + different_columns[1][row]
        + different_columns[2][row]
        for row in range(len(color_names))
    ) == target
    assert matrix_rank(
        tuple(
            row + (target[index],)
            for index, row in enumerate(different_matrix)
        )
    ) == 3

    # The old base-factor architecture has columns d,q,r and a finite-order
    # q selector whose finite colored orders vanish.
    base_columns = tuple(
        tuple(row[column] for row in datum.valuation_matrix)
        for column in (4, 5, 6, 7)
    )
    base_matrix = transpose_columns(base_columns)
    old_witnesses = colored_proportionality_witnesses(
        base_matrix, target, color_names
    )
    assert len(old_witnesses) == 6

    breaker_columns = tuple(
        packet_columns[packet_names.index(name)] for name in PRIMARY_BREAKERS
    )
    breaker_matrix = transpose_columns(base_columns + breaker_columns)
    breaker_rank = matrix_rank(breaker_matrix)
    breaker_augmented_rank = matrix_rank(
        tuple(row + (target[index],) for index, row in enumerate(breaker_matrix))
    )
    assert breaker_augmented_rank > breaker_rank

    broken_relations: list[dict[str, object]] = []
    for witness in old_witnesses:
        left_name, right_name = witness["colors"]
        left = color_index[left_name]
        right = color_index[right_name]
        left_scale, right_scale = witness["row_scales"]
        breakers: list[dict[str, object]] = []
        for packet_name, column in zip(packet_names, packet_columns):
            mismatch = right_scale * column[left] - left_scale * column[right]
            if mismatch:
                breakers.append({"packet": packet_name, "mismatch": mismatch})
        assert breakers
        broken_relations.append(
            {
                "colors": witness["colors"],
                "old_mismatch": witness["mismatch"],
                "packet_breakers": tuple(breakers),
            }
        )

    different_broken_relations: list[dict[str, object]] = []
    for witness in old_witnesses:
        left_name, right_name = witness["colors"]
        left = color_index[left_name]
        right = color_index[right_name]
        left_scale, right_scale = witness["row_scales"]
        breakers = tuple(
            {
                "column": name,
                "mismatch": (
                    right_scale * column[left]
                    - left_scale * column[right]
                ),
            }
            for name, column in zip(different_names, different_columns)
            if (
                right_scale * column[left]
                - left_scale * column[right]
            )
        )
        assert breakers
        different_broken_relations.append(
            {"colors": witness["colors"], "breakers": breakers}
        )

    compact_columns = base_columns + different_columns
    compact_matrix = transpose_columns(compact_columns)
    compact_rank = matrix_rank(compact_matrix)
    compact_augmented_rank = matrix_rank(
        tuple(
            row + (target[index],)
            for index, row in enumerate(compact_matrix)
        )
    )
    assert (compact_rank, compact_augmented_rank) == (6, 6)

    compact_models = tuple(
        (exponent_d, exponent_q, exponent_r)
        for exponent_d in range(9)
        for exponent_q in range(9)
        for exponent_r in range(9)
        if all(
            exponent_d * different_columns[0][row]
            + exponent_q * different_columns[1][row]
            + exponent_r * different_columns[2][row]
            == target[row]
            for row in range(len(color_names))
        )
    )
    assert compact_models == ((3, 1, 1),)

    full_columns = base_columns + packet_columns
    full_matrix = transpose_columns(full_columns)
    full_rank = matrix_rank(full_matrix)
    full_augmented_rank = matrix_rank(
        tuple(row + (target[index],) for index, row in enumerate(full_matrix))
    )
    assert (full_rank, full_augmented_rank) == (19, 19)

    # A visibly unimodular minor proves saturation without enumerating the
    # enormous family of 19-minors of a 63-by-20 matrix.  The first three
    # rows kill the base d,q,r columns; one pivot from every disjoint packet
    # then gives an identity block after integral row operations.
    pivot_colors = (
        "d_unramified",
        "q_residual_1",
        "r_unramified",
    ) + tuple(colors[0] for _name, colors in PACKETS)
    nonzero_column_indices = (0, 1, 2) + tuple(range(4, 20))
    saturated_minor = tuple(
        tuple(
            full_matrix[color_index[color]][column]
            for column in nonzero_column_indices
        )
        for color in pivot_colors
    )
    saturated_determinant = determinant(saturated_minor)
    assert abs(saturated_determinant) == 1

    # The nonnegative solution is forced.  Target-zero rows force all three
    # base exponents to zero.  Disjoint packet supports then force one packet
    # exponent at a time.  The zero finite-order selector is normalized to
    # exponent zero modulo its descended norm.
    semigroup_model = {
        "mask_d": 0,
        "mask_q": 0,
        "mask_r": 0,
        "q_selector_w_minus_1": 0,
        **{f"cox_{name}": order for name, order in packet_orders.items()},
    }
    assert tuple(
        sum(
            semigroup_model[name] * full_matrix[row][column]
            for column, name in enumerate(semigroup_model)
        )
        for row in range(len(color_names))
    ) == target

    # The q-conductor lattice from the exact normalization calculation.
    q_pullback = (
        (-1, -1, -2),
        (-1, -1, -2),
        (1, 0, 0),
        (0, 0, 1),
    )
    q_completion = tuple(
        row + (int(index == 0),) for index, row in enumerate(q_pullback)
    )
    assert smith_invariant_factors(q_pullback) == (1, 1, 1)
    q_completion_determinant = determinant(q_completion)
    assert q_completion_determinant == -1
    assert target[color_index["q_collision_plus"]] == target[
        color_index["q_collision_minus"]
    ]

    # Natural principal sheet selectors are not boundary units: their norms
    # introduce interior divisors.  This proves that the packet columns are
    # genuinely Cox/line-bundle data rather than already-known functions.
    P, d, q, X = exact_f20_polynomial()
    s, t = sp.symbols("s t")
    h_d = 4 * X - 1
    h_q = 2 * (s - 1) * X - (2 * s**2 * t + 2 * s**2 + 3 * s - 4)
    norm_d = sp.factor(sp.resultant(P, h_d, X))
    norm_q = sp.factor(sp.resultant(P, h_q, X))
    assert sp.factor(norm_d + 4 * d * (t - 2) ** 2) == 0
    q_cofactor = sp.cancel(-2 * norm_q / q**2)
    assert sp.Poly(q_cofactor, s, t).total_degree() > 0
    assert sp.rem(sp.Poly(q_cofactor, t), sp.Poly(q, t)) != 0

    # If the packet variables were the complete unit basis of a factorial
    # core, the saturated rank-19 presentation would have one residual unit
    # (the zero-order q selector) and free class rank 63-19=44.  This is a
    # conditional affine-space obstruction, not a claim that the missing Cox
    # algebra satisfies the factorial-core hypotheses.
    packet_affine_screen = {
        "conditional_on_factorial_core_and_complete_boundary": True,
        "boundary_count": 63,
        "character_count": 20,
        "matrix_rank": 19,
        "unit_rank": 1,
        "class_group_free_rank": 44,
        "class_group_torsion": (),
        "saturated_minor_determinant": saturated_determinant,
        "passes_affine_space_necessary_conditions": False,
    }
    compact_affine_screen = {
        "conditional_on_factorial_core_and_complete_boundary": True,
        "boundary_count": 63,
        "character_count": 7,
        "matrix_rank": 6,
        "unit_rank": 1,
        "class_group_free_rank": 57,
        "passes_affine_space_necessary_conditions": False,
    }

    certificate: dict[str, object] = {
        "status": "conductor_residue_uncertified",
        "scope": (
            "the sixteen Galois-stable positive-derivative packets on the "
            "certified finite F20 color atlas; global regular SNC resolution "
            "and Cox-ring realization remain assumptions"
        ),
        "color_count": len(color_names),
        "packet_count": len(PACKETS),
        "packets": {
            name: {
                "colors": colors,
                "primitive_column": {
                    color: 1 for color in colors
                },
                "target_exponent": packet_orders[name],
            }
            for name, colors in PACKETS
        },
        "cartier_gate": {
            "status": "conditional",
            "local_statement": (
                "on a regular resolved colored surface, every packet orbit "
                "sum is an effective Cartier divisor"
            ),
            "global_regular_snc_resolution_certified": False,
            "global_cox_sections_constructed": False,
            "compact_different_factor_columns": {
                name: packet_coefficients
                for name, packet_coefficients in DIFFERENT_FACTOR_PACKETS.items()
            },
            "different_identity": "div(P_X)=3*D_d+D_q+D_r",
        },
        "compact_three_column_gate": {
            "status": "passes_combinatorially",
            "columns": different_names,
            "generator_rank": compact_rank,
            "augmented_rank": compact_augmented_rank,
            "unique_nonnegative_model": {
                "D_d": compact_models[0][0],
                "D_q": compact_models[0][1],
                "D_r": compact_models[0][2],
            },
            "broken_relations": tuple(different_broken_relations),
            "scope": (
                "Cartier on a global regular colored resolution; the "
                "individual Cox sections and their descent data are not "
                "constructed"
            ),
        },
        "six_witness_gate": {
            "status": "passes",
            "primary_breaker_packets": PRIMARY_BREAKERS,
            "old_witnesses": tuple(broken_relations),
            "six_breaker_generator_rank": breaker_rank,
            "six_breaker_augmented_rank": breaker_augmented_rank,
            "full_target_spanned_by_six_breakers": False,
        },
        "integral_span_gate": {
            "status": "passes",
            "generator_rank": full_rank,
            "augmented_rank": full_augmented_rank,
            "target_class_order": 1,
            "saturated_minor_determinant": saturated_determinant,
        },
        "nonnegative_semigroup_gate": {
            "status": "passes",
            "model": semigroup_model,
            "model_is_unique_after_selector_normalization": True,
            "proof": (
                "target-zero d/q/r rows force base exponents zero; the "
                "sixteen packet columns have pairwise disjoint support"
            ),
        },
        "conductor_gluing_gate": {
            "status": "uncertified",
            "packet_orders_are_galois_orbit_constant": True,
            "q_crossing_orders_match": True,
            "q_pullback_smith_diagonal": (1, 1, 1),
            "q_selector_completion_determinant": q_completion_determinant,
            "integral_unit_lattice_obstruction": False,
            "residue_cocycle_certificate": None,
            "reason": (
                "equal valuations and a unimodular unit lattice do not "
                "determine the conductor fiber identifications of the "
                "individual Cox sections"
            ),
        },
        "principal_selector_screen": {
            "status": "fails_for_tested_natural_selectors",
            "d_simple_selector_norm": str(norm_d),
            "d_interior_factor": "(t - 2)^2",
            "q_collision_selector_has_nonconstant_cofactor": True,
            "conclusion": (
                "the primitive packet columns are not supplied by these "
                "obvious principal boundary units"
            ),
        },
        "nonlinear_pipeline": {
            "tropical_packet_survivor": True,
            "conductor_certified_survivor": False,
            "divisorial_inverse_adjugate_cancellation": "passes_formally",
            "entrywise_inverse_adjugate_polynomiality": "not_reached",
            "affine_space_recognition": "not_reached",
            "reason": (
                "the pipeline stops at the missing conductor residue cocycle "
                "and global Cox algebra"
            ),
        },
        "independent_packet_affine_screen": packet_affine_screen,
        "independent_compact_affine_screen": compact_affine_screen,
        "full_geometric_color_basis": {
            "column_count": 63,
            "valuation_matrix": "I_63",
            "unit_rank": 0,
            "class_group_free_rank": 0,
            "class_group_torsion": (),
            "status": "lattice_neutral_only",
            "unresolved": (
                "Galois/conductor descent, Cox relations and irrelevant "
                "locus, dimension-preserving quotient, and affine-space "
                "recognition of the resulting source"
            ),
        },
    }
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n")
    print(rendered)
    print("PASS: six primitive F20 Cox packets break all six row witnesses")
    print("PASS: three different-factor Cartier columns give model (3,1,1)")
    print("PASS: sixteen packets give the unique nonnegative derivative model")
    print("PASS: q-conductor orders and unit lattice pass integrally")
    print("SCOPE: conductor residues and the global Cox algebra remain uncertified")


if __name__ == "__main__":
    main()
