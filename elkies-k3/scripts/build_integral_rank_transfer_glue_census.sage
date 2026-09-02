#!/usr/bin/env python3
"""Build the equation-free integral rank-transfer/glue census.

The output deliberately separates three operations which look similar in an
equation calculation but are different integrally:

* replacing the primitive hyperbolic plane ``U`` inside a fixed NS lattice;
* graph gluing an auxiliary and its orthogonal complement in a Niemeier
  lattice; and
* saturating invariant and anti-invariant Mordell--Weil character lattices.

Run with Sage's Python because the finite quadratic-form normalizations are
computed from the pinned integral Gram matrices, rather than copied as prose.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import Genus, QQ, ZZ, block_diagonal_matrix, matrix


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-integral-rank-transfer-glue-census-v1.json"
)


def load_json(relative_path):
    return json.loads((ROOT / relative_path).read_text())


def load_matrix(relative_path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in (ROOT / relative_path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def digest(relative_path):
    return hashlib.sha256((ROOT / relative_path).read_bytes()).hexdigest()


def rows_as_strings(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def discriminant_form_key(gram):
    normal = Genus(gram).discriminant_form().normal_form()
    return {
        "invariants": list(map(int, normal.invariants())),
        "quadratic_gram": rows_as_strings(normal.gram_matrix_quadratic()),
        "value_module": str(normal.value_module_qf()),
    }


def route_edges(source, moves):
    edges = []
    current_root, current_rank, current_mw = source
    assert current_rank + current_mw == 17
    for edge_index, (q, target_root, target_rank, target_mw) in enumerate(moves, 1):
        assert target_rank + target_mw == 17
        edges.append(
            {
                "edge_index": edge_index,
                "operation": "primitive_U_change",
                "old_fibre_degree": 2,
                "q": q,
                "source_root_type": current_root,
                "source_root_rank": current_rank,
                "source_mw_rank": current_mw,
                "target_root_type": target_root,
                "target_root_rank": target_rank,
                "target_mw_rank": target_mw,
                "root_rank_change": target_rank - current_rank,
                "mw_rank_change": target_mw - current_mw,
                "ambient_discriminant_form_change": "none",
                "integral_transport_index": 1,
            }
        )
        current_root, current_rank, current_mw = target_root, target_rank, target_mw
    return edges


def ambient_tuple(name, gram, group, eigensublattices, glue_subgroups):
    smith = gram.smith_form()[0].diagonal()
    return {
        "L": {
            "name": name,
            "rank": gram.nrows(),
            "determinant_absolute": abs(int(gram.det())),
        },
        "G": group,
        "L_chi": eigensublattices,
        "A_L": [abs(int(value)) for value in smith if abs(int(value)) > 1],
        "q_L": discriminant_form_key(gram),
        "glue_subgroups": glue_subgroups,
    }


def assert_character_partition(record):
    ranks = [piece["rank"] for piece in record["L_chi"]]
    assert all(isinstance(rank, int) for rank in ranks)
    assert sum(ranks) == record["L"]["rank"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    j2_path = "artifacts/generated-results/elkies-k3-rootless-j2-completeness-track.json"
    foundry_path = "artifacts/generated-results/elkies-k3-lattice-foundry-v1.json"
    ns0024_route_path = (
        "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-r13-nef-route.json"
    )
    golay_design_path = (
        "artifacts/generated-results/elkies-k3-golay-octad-rank17-det720.json"
    )
    golay_ns_path = (
        "artifacts/generated-results/elkies-k3-golay-det720-ns-saturation-v1.json"
    )
    golay_route_path = (
        "artifacts/generated-results/elkies-k3-golay-det720-3a5-to-mw17-physical-corridor-v1.json"
    )
    golay_source_path = (
        "artifacts/generated-results/elkies-k3-golay-octad-det720-source-niemeier.json"
    )
    e6_21_path = (
        "artifacts/generated-results/elkies-k3-e6-ii-rank3-quadratic-base-change-v1.json"
    )
    e6_22_path = (
        "artifacts/generated-results/elkies-k3-e6-rank4-linear-chord-incidence-v1.json"
    )
    r17_genus1_path = (
        "artifacts/generated-results/elkies-k3-r17-rank28-genus-one-bisection-pilot-v1.json"
    )
    r17_pair_path = (
        "artifacts/generated-results/elkies-2026-bisection-pair-cover-geometry-full.json"
    )
    r17_split_path = (
        "artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json"
    )
    golay_rejection_path = (
        "artifacts/generated-results/elkies-k3-golay-det720-3a5-saturation-rejection-v1.json"
    )

    j2 = load_json(j2_path)
    foundry = load_json(foundry_path)
    ns0024_route = load_json(ns0024_route_path)
    golay_design = load_json(golay_design_path)
    golay_ns = load_json(golay_ns_path)
    golay_route = load_json(golay_route_path)
    golay_source = load_json(golay_source_path)
    e6_21 = load_json(e6_21_path)
    e6_22 = load_json(e6_22_path)
    r17_genus1 = load_json(r17_genus1_path)
    r17_pair = load_json(r17_pair_path)
    r17_split = load_json(r17_split_path)
    golay_rejection = load_json(golay_rejection_path)

    assert j2["status"] == "PASS_EXACT_ROOTLESS_J2_CONTROLS_AND_OFF_GENUS_REJECTION_NOT_COMPLETE"
    assert ns0024_route["status"] == "PASS_EXACT_NEW_K3_ROOTFUL_TO_ROOTLESS_NEF_ROUTE"
    assert golay_design["status"] == "PASS_EXACT_GOLAY_OCTAD_RANK17_DET720_LATTICE_DESIGN"
    assert golay_ns["status"] == "PASS_EXACT_PRIMITIVE_K3_NS_DET720_AND_DISCRIMINANT_FORM"
    assert golay_source["status"] == "PASS_EXACT_DET720_NONCYCLIC_GLUE_AND_NIEMEIER_SOURCE"
    assert golay_source["discriminant_glue"]["invariants"] == [2, 6, 60]
    assert golay_source["discriminant_glue"]["graph_glue_order"] == 720
    assert (
        golay_source["discriminant_glue"][
            "anti_isometry_count_in_smith_generator_enumeration"
        ]
        == 96
    )
    assert e6_21["status"] == "PASS_EXACT_E6_II_RANK_SUM_3_RHO19_ROOTLESS_IMPOSSIBLE"
    assert e6_22["status"] == "PASS_EXACT_E6_RANK4_INCIDENCE_DESCENT"
    assert r17_genus1["status"] == "PASS_EXACT_R17_RANK28_GENUS_ONE_BISECTION_PILOT"
    assert r17_pair["status"] == "PASS_COMPLETE_CONIC_AND_GENUS_ONE_PAIR_CLASSIFICATION"
    assert r17_split["status"] == "PASS_EXACT_SIMULTANEOUS_SPLIT_HITS_QUOTIENTED"
    assert golay_rejection["status"] == (
        "PASS_EXACT_RATIONAL_POINT_REJECTED_NS_DET20_TORSION3_HALF_SECTION"
    )

    U = matrix(ZZ, [[0, 1], [1, 0]])
    r17 = load_matrix("elkies-k3/data/lattice/rank17_gram.txt")
    h3_ns = block_diagonal_matrix(U, -r17)
    assert abs(h3_ns.det()) == 948

    ns0024 = next(row for row in foundry["ns_classes"] if row["ns_id"] == "NS0024")
    ns0024_gram = matrix(ZZ, ns0024["ns_gram_representative"])
    assert abs(ns0024_gram.det()) == 950

    golay_frame = matrix(ZZ, golay_design["frame"]["gram"])
    golay_ns_gram = block_diagonal_matrix(U, -golay_frame)
    assert abs(golay_ns_gram.det()) == 720
    # The saturation artifact stores the literal sign-negation of the positive
    # frame key.  Re-normalizing that equivalent form need not print the same
    # rational representatives in Q/2Z, so compare group invariants here and
    # compute the canonical key directly from the NS Gram below.
    assert discriminant_form_key(golay_ns_gram)["invariants"] == golay_ns[
        "neron_severi"
    ]["discriminant_group_invariants"]

    e6_21_ns = matrix(ZZ, e6_21["k3"]["neron_severi_gram"])
    e6_22_ns = matrix(ZZ, e6_22["neron_severi"]["integral_gram"])

    # The height pairings are rational.  Multiplying by 12 gives even
    # integral lattices on which A_L and q_L are literal finite quadratic
    # forms.  This is only a bookkeeping model; the unscaled height
    # determinants are retained below.
    e6_21_character_gram = block_diagonal_matrix(
        matrix(ZZ, [[16, -8], [-8, 16]]), matrix(ZZ, [[8]])
    )
    e6_22_character_gram = block_diagonal_matrix(
        matrix(ZZ, [[16, -8], [-8, 16]]),
        matrix(ZZ, [[88, 16], [16, 88]]),
    )
    e6_22_saturated_gram = matrix(
        ZZ,
        [[16, -8, 8, -4], [-8, 16, -4, 8], [8, -4, 26, 2], [-4, 8, 2, 26]],
    )
    e6_22_half_sum_basis = matrix(
        QQ,
        [
            [1, 0, QQ(1) / 2, 0],
            [0, 1, 0, QQ(1) / 2],
            [0, 0, QQ(1) / 2, 0],
            [0, 0, 0, QQ(1) / 2],
        ],
    )
    assert e6_21_character_gram.det() == ZZ(12) ** 3 * ZZ(8) / 9
    assert e6_22_character_gram.det() == ZZ(12) ** 4 * ZZ(208) / 3
    assert (
        e6_22_half_sum_basis.transpose()
        * e6_22_character_gram
        * e6_22_half_sum_basis
        == e6_22_saturated_gram
    )
    assert e6_22_saturated_gram.det() == e6_22_character_gram.det() / 16

    h3_edges = route_edges(
        ("E8+E7", 15, 2),
        [
            (6, "E8+E6", 14, 3),
            (8, "D13", 13, 4),
            (24, "D12", 12, 5),
            (6, "A11", 11, 6),
            (8, "2A5", 10, 7),
            (4, "3A3", 9, 8),
            (4, "A3+2A2", 7, 10),
            (4, "5A1", 5, 12),
            (4, "4A1", 4, 13),
            (4, "3A1", 3, 14),
            (4, "2A1", 2, 15),
            (4, "A1", 1, 16),
            (6, "rootless", 0, 17),
        ],
    )
    q80_edges = route_edges(
        ("E6+D5+A3", 14, 3),
        [
            (4, "D9+A4", 13, 4),
            (4, "D7+D5", 12, 5),
            (6, "D7+D4", 11, 6),
            (4, "A6+A4", 10, 7),
            (4, "A6+A3", 9, 8),
            (6, "A4+A2+A1", 7, 10),
            (4, "A3+A2", 5, 12),
            (4, "4A1", 4, 13),
            (4, "A1", 1, 16),
            (6, "rootless", 0, 17),
        ],
    )
    ns0024_edges = route_edges(
        ("A3+A4+A6", 13, 4),
        [
            (4, "A1+A2+A4+D5", 12, 5),
            (4, "A1+A5+D5", 11, 6),
            (4, "2A1+A8", 10, 7),
            (4, "A1+A2+A6", 9, 8),
            (4, "A1+A7", 8, 9),
            (4, "A3+A4", 7, 10),
            (4, "A6", 6, 11),
            (6, "A5", 5, 12),
            (6, "A4", 4, 13),
            (6, "A3", 3, 14),
            (6, "A2", 2, 15),
            (6, "A1", 1, 16),
            (6, "rootless", 0, 17),
        ],
    )
    assert [edge["q"] for edge in ns0024_edges] == [
        edge["q"] for edge in ns0024_route["edges"]
    ]

    golay_edges = route_edges(
        ("3A5", 15, 2),
        [
            (4, "4A2+A5", 13, 4),
            (4, "3A1+2A2+A3", 10, 7),
            (4, "4A1+A2", 6, 11),
            (4, "3A1", 3, 14),
            (4, "2A1", 2, 15),
            (4, "rootless", 0, 17),
        ],
    )
    assert [edge["target_root_rank"] for edge in golay_edges] == [
        edge["target_root_rank"] for edge in golay_route["route"]["edges"]
    ]

    h3_tuple = ambient_tuple(
        "NS_H3",
        h3_ns,
        {"name": "trivial_on_same_NS_corridor", "order": 1},
        [
            {
                "character": "trivial",
                "description": "the whole NS_H3 lattice",
                "rank": 19,
                "discriminant_form": discriminant_form_key(h3_ns),
            }
        ],
        [
            {
                "context": "primitive_U_split",
                "subgroup": "0",
                "index": 1,
                "reason": "U is unimodular",
            },
            {
                "context": "J2_Niemeier_complement",
                "subgroup": "graph(phi) in A_K plus A_W",
                "structure": "Z/948",
                "order": 948,
                "anti_isometry_count_per_control": 8,
                "ambient": "N(2A7+2D5)",
                "auxiliary_q_generator": "1267/948 in Q/2Z",
            },
        ],
    )

    ns0024_tuple = ambient_tuple(
        "NS0024",
        ns0024_gram,
        {"name": "trivial_on_same_NS_corridor", "order": 1},
        [
            {
                "character": "trivial",
                "description": "the whole NS0024 lattice",
                "rank": 19,
                "discriminant_form": discriminant_form_key(ns0024_gram),
            }
        ],
        [
            {
                "context": "primitive_U_split",
                "subgroup": "0",
                "index": 1,
                "reason": "U is unimodular",
            },
            {
                "context": "source_Niemeier_complement",
                "subgroup": "graph of a cyclic discriminant anti-isometry",
                "order": 950,
                "ambient": "N(A15+D9)",
            },
        ],
    )

    golay_tuple = ambient_tuple(
        "NS_Golay720",
        golay_ns_gram,
        {"name": "trivial_on_same_NS_corridor", "order": 1},
        [
            {
                "character": "trivial",
                "description": "the whole NS_Golay720 lattice",
                "rank": 19,
                "discriminant_form": discriminant_form_key(golay_ns_gram),
            }
        ],
        [
            {
                "context": "Golay_octad_complement",
                "subgroup": "graph(phi) in A_K plus A_W",
                "structure": "Z/2 + Z/6 + Z/60",
                "order": 720,
                "ambient": "N(24A1)",
            },
            {
                "context": "rootful_source_complement",
                "subgroup": "selected three-generator graph glue",
                "structure": "Z/2 + Z/6 + Z/60",
                "order": 720,
                "anti_isometry_count": golay_source["discriminant_glue"][
                    "anti_isometry_count_in_smith_generator_enumeration"
                ],
                "source_smith_generator_images": golay_source[
                    "discriminant_glue"
                ]["selected_images_in_source_smith_coordinates"],
                "source_smith_generator_orders": golay_source[
                    "discriminant_glue"
                ]["selected_image_orders"],
                "ambient": "N(4A5+D4)",
            },
            {
                "context": "primitive_U_split",
                "subgroup": "0",
                "index": 1,
                "reason": "U is unimodular",
            },
        ],
    )

    e6_21_tuple = ambient_tuple(
        "12_times_E6_2_plus_1_MW_character_lattice",
        e6_21_character_gram,
        {"name": "quadratic_base_change_deck", "order": 2},
        [
            {
                "character": "+",
                "role": "12 times the Mordell-Weil height lattice",
                "rank": 2,
                "gram": [[16, -8], [-8, 16]],
                "unscaled_determinant": "4/3",
            },
            {
                "character": "-",
                "role": "12 times the Mordell-Weil height lattice",
                "rank": 1,
                "gram": [[8]],
                "unscaled_determinant": "2/3",
            },
        ],
        [
            {
                "context": "MW_character_saturation",
                "subgroup": "0",
                "index": 1,
                "determinant_check": "(4/3)*(2/3)=8/9",
                "kind": "integral_isotypic_saturation",
            }
        ],
    )

    e6_22_tuple = ambient_tuple(
        "12_times_E6_2_plus_2_pure_MW_character_lattice",
        e6_22_character_gram,
        {"name": "quadratic_base_change_deck", "order": 2},
        [
            {
                "character": "+",
                "role": "12 times the pure Mordell-Weil height lattice",
                "basis": ["P", "Q"],
                "rank": 2,
                "gram": [[16, -8], [-8, 16]],
                "unscaled_determinant": "4/3",
            },
            {
                "character": "-",
                "role": "12 times the pure Mordell-Weil height lattice",
                "basis": ["T1", "T2"],
                "rank": 2,
                "gram": [[88, 16], [16, 88]],
                "unscaled_determinant": "52",
            },
        ],
        [
            {
                "context": "MW_character_saturation",
                "subgroup": "<(P+T1)/2, (Q+T2)/2>",
                "structure": "(Z/2)^2",
                "index": 4,
                "determinant_before": "208/3",
                "determinant_after": "13/3",
                "scaled_determinant_before": int(e6_22_character_gram.det()),
                "scaled_determinant_after": int(e6_22_saturated_gram.det()),
                "kind": "even_discriminant_glue_after_scaling_height_by_12",
                "scaled_saturated_gram": rows_as_strings(e6_22_saturated_gram),
            }
        ],
    )

    genus1_visible = block_diagonal_matrix(2 * r17, matrix(ZZ, [[16]]))
    genus1_tuple = ambient_tuple(
        "orthogonal_character_sublattice_on_R17_genus_one_bisection_cover",
        genus1_visible,
        {"name": "quadratic_base_change_deck", "order": 2},
        [
            {
                "character": "+",
                "role": "pulled-back invariant Mordell-Weil lattice",
                "rank": 17,
                "gram": "R17(2)",
            },
            {
                "character": "-",
                "role": "certified anti-invariant line",
                "rank": 1,
                "gram": [[16]],
            },
        ],
        [
            {
                "context": "certified_half_sum_extension",
                "subgroup": "<R mod (R17(2)+<T>)>",
                "structure": "Z/2",
                "index": 2,
                "determinant_after_known_extension": int(genus1_visible.det() / 4),
                "known_relation": "for conjugate lifts R,Rprime, R+Rprime=tau",
                "full_saturation_index": "not determined",
                "warning": "This known index-2 extension need not be the full MW lattice or the full anti-invariant space.",
            }
        ],
    )

    paired_visible = block_diagonal_matrix(
        4 * r17, matrix(ZZ, [[24]]), matrix(ZZ, [[24]])
    )
    paired_tuple = ambient_tuple(
        "orthogonal_character_sublattice_on_R17_paired_rational_bisection_compositum",
        paired_visible,
        {"name": "biquadratic_deck", "structure": "V4", "order": 4},
        [
            {
                "character": "1",
                "role": "pulled-back invariant Mordell-Weil lattice",
                "rank": 17,
                "gram": "R17(4)",
            },
            {
                "character": "chi_1",
                "role": "first pulled anti-invariant line",
                "rank": 1,
                "gram": [[24]],
            },
            {
                "character": "chi_2",
                "role": "second pulled anti-invariant line",
                "rank": 1,
                "gram": [[24]],
            },
            {
                "character": "chi_1_chi_2",
                "role": "product-twist eigenspace inside this declared orthogonal character sublattice",
                "rank": 0,
                "full_MW_rank_warning": "The product-character rank in the full MW lattice is not determined.",
            },
        ],
        [
            {
                "context": "two_certified_half_sum_extensions",
                "subgroup": "<R1,R2> modulo (R17(4)+<T1>+<T2>)",
                "structure": "(Z/2)^2",
                "index": 4,
                "determinant_after_known_extension": int(paired_visible.det() / 16),
                "known_relations": ["2R1=tau1+T1", "2R2=tau2+T2"],
                "full_saturation_index": "not determined",
                "warning": "The separate anti-character coordinates make these classes independent, but further saturation and product-character rank are unknown.",
            }
        ],
    )

    e6_21_tuple["bookkeeping"] = {
        "minimum_in_scaled_character_lattice": 8,
        "norm_two_roots_in_scaled_character_lattice": 0,
        "saturation_index": 1,
        "ambient_NS": {
            "rank": 19,
            "determinant_absolute": abs(int(e6_21_ns.det())),
            "A_NS": [
                abs(int(value))
                for value in e6_21_ns.smith_form()[0].diagonal()
                if abs(int(value)) > 1
            ],
            "q_NS": discriminant_form_key(e6_21_ns),
            "root_system": "2E6+A2",
            "root_rank": 14,
        },
    }
    e6_22_tuple["bookkeeping"] = {
        "minimum_in_scaled_pure_character_lattice": 16,
        "norm_two_roots_in_scaled_pure_character_lattice": 0,
        "saturation_index": 4,
        "ambient_NS": {
            "rank": 19,
            "determinant_absolute": abs(int(e6_22_ns.det())),
            "A_NS": [
                abs(int(value))
                for value in e6_22_ns.smith_form()[0].diagonal()
                if abs(int(value)) > 1
            ],
            "q_NS": discriminant_form_key(e6_22_ns),
            "root_system": "2E6+A1",
            "root_rank": 13,
        },
    }
    genus1_tuple["bookkeeping"] = {
        "minimum_of_orthogonal_character_sublattice": 8,
        "norm_two_roots": 0,
        "certified_extension_index": 2,
        "full_saturation_index": "not determined",
    }
    paired_tuple["bookkeeping"] = {
        "minimum_of_orthogonal_character_sublattice": 16,
        "norm_two_roots": 0,
        "certified_extension_index": 4,
        "full_saturation_index": "not determined",
    }
    for record in [
        h3_tuple,
        ns0024_tuple,
        golay_tuple,
        e6_21_tuple,
        e6_22_tuple,
        genus1_tuple,
        paired_tuple,
    ]:
        assert_character_partition(record)

    all_edges = h3_edges + q80_edges + ns0024_edges + golay_edges
    q_histogram = {}
    for edge in all_edges:
        q_histogram[str(edge["q"])] = q_histogram.get(str(edge["q"]), 0) + 1

    payload = {
        "schema": "elkies-k3.integral-rank-transfer-glue-census.v1",
        "status": "PASS_EXACT_EXTRACTION_WITH_EXPLICIT_OPEN_GLUE_FIELDS",
        "reproduce": "sage -python elkies-k3/scripts/build_integral_rank_transfer_glue_census.sage --check",
        "inputs": {
            path: digest(path)
            for path in [
                j2_path,
                foundry_path,
                ns0024_route_path,
                golay_design_path,
                golay_ns_path,
                golay_route_path,
                golay_source_path,
                e6_21_path,
                e6_22_path,
                r17_genus1_path,
                r17_pair_path,
                r17_split_path,
                golay_rejection_path,
                "elkies-k3/data/lattice/rank17_gram.txt",
            ]
        },
        "tuple_order": ["L", "G", "L_chi", "A_L", "q_L", "glue_subgroups"],
        "operation_classes": {
            "primitive_U_change": {
                "definition": "Replace primitive U by primitive Uprime in fixed NS and take Wprime=Uprime^perp.",
                "finite_form_effect": "None: U is unimodular, so q_Wprime is fixed by q_NS.",
                "root_effect": "Roots may leave or enter the frame because the rational orthogonal complement changes.",
                "reversible": True,
            },
            "maximal_graph_glue": {
                "definition": "Glue primitive orthogonal lattices K,W along the graph of an anti-isometry q_K to -q_W.",
                "finite_form_effect": "The maximal isotropic graph kills the full discriminant and produces a unimodular ambient.",
                "root_effect": "The complement changes when the embedding or anti-isometry changes; roots must be tested in the actual complement.",
                "reversible": True,
            },
            "character_saturation": {
                "definition": "Adjoin |G|-primary finite-index classes to full integral rational-isotypic lattices, or record the quotient first when only rational height blocks are available.",
                "finite_form_effect": "For the full integral C2 isotypic lattices the saturation quotient is killed by two; a visible or rational-block quotient requires a separate proof. Any actual index changes determinant by its square.",
                "integrality_warning": "A finite quotient between rational Mordell-Weil height blocks is not by itself a Nikulin discriminant-form glue subgroup.",
                "root_effect": "An overlattice alone cannot remove roots already present.",
                "reversible": True,
            },
            "Kneser_reglue_candidate": {
                "definition": "Replace, rather than merely enlarge, isotropic glue across a common index-p sublattice.",
                "finite_form_effect": "Preserves the genus for ordinary p-neighbors while removing old cosets and adding new ones.",
                "root_effect": "Can annihilate roots and is the appropriate atomic candidate for an inverted calculus.",
                "observed_edge_identification": "not yet computed for the four U-change corridors",
                "reversible": True,
            },
        },
        "experiments": [
            {
                "id": "historical_H3_to_J2_to_R17",
                "tuple": h3_tuple,
                "transitions": h3_edges,
                "frame_family": {
                    "definition": "W_i=U_i^perp(-1)",
                    "rank": 17,
                    "determinant_absolute": 948,
                    "discriminant_form": discriminant_form_key(r17),
                    "terminal_minimum": 4,
                    "terminal_roots": 0,
                    "saturation_index_in_NS": 1,
                },
                "J2_result": {
                    "meaning": "frame integral-isometry class, not a geometric arrow",
                    "rootless_class_count": 2,
                    "classes": ["published R17", "alternate Q80-derived"],
                    "minimum": 4,
                    "norm_four_pairs": [1311, 1313],
                    "automorphism_orders": [2, 4],
                },
                "proof_boundary": "The H3 route and J2 classification are exact; J1 surface-automorphism classification is not claimed.",
            },
            {
                "id": "Q80_low_q_corridor",
                "tuple": h3_tuple,
                "transitions": q80_edges,
                "frame_family": {
                    "definition": "W_i=U_i^perp(-1)",
                    "rank": 17,
                    "determinant_absolute": 948,
                    "discriminant_form": discriminant_form_key(r17),
                    "terminal_minimum": 4,
                    "terminal_roots": 0,
                    "saturation_index_in_NS": 1,
                },
                "endpoint": "alternate Q80-derived rootless J2 class",
                "proof_boundary": "Exact generic lattice corridor; the alternate generic characteristic-zero equation is open.",
            },
            {
                "id": "NS0024_corridor",
                "tuple": ns0024_tuple,
                "transitions": ns0024_edges,
                "frame_family": {
                    "definition": "W_i=U_i^perp(-1)",
                    "rank": 17,
                    "determinant_absolute": 950,
                    "discriminant_form": ns0024[
                        "discriminant_form_normal_key"
                    ],
                    "terminal_minimum": 4,
                    "terminal_roots": 0,
                    "saturation_index_in_NS": 1,
                },
                "endpoint": "NS0024-F005",
                "proof_boundary": "Exact lattice/arithmetic class and marked route; equation and rational descent are open.",
            },
            {
                "id": "Golay_720_construction_and_corridor",
                "tuple": golay_tuple,
                "transitions": golay_edges,
                "frame_family": {
                    "definition": "W_i=U_i^perp(-1)",
                    "rank": 17,
                    "determinant_absolute": 720,
                    "discriminant_form": golay_design["frame"][
                        "discriminant_form_normal_key"
                    ],
                    "terminal_roots": 0,
                    "saturation_index_in_NS": 1,
                },
                "rootless_minimum": 4,
                "rootless_norm_four_vectors": 3064,
                "proof_boundary": "Exact abstract K3 NS, complement glue, and lattice corridor; a rational source equation is open.",
            },
            {
                "id": "E6_quadratic_base_change_2_plus_1",
                "tuple": e6_21_tuple,
                "transition": {
                    "operation": "quadratic_pullback_then_character_split",
                    "rank_split": "2+1",
                    "character_glue_index": 1,
                },
                "proof_boundary": "Exact generic ranks, saturation, torsion and determinant 24; same-NS rootless MW17 is impossible by the determinant bound.",
            },
            {
                "id": "E6_quadratic_base_change_2_plus_2",
                "tuple": e6_22_tuple,
                "transition": {
                    "operation": "quadratic_pullback_then_character_split_and_2_glue",
                    "geometric_rank_split": "2+2",
                    "descent_warning": "The ordered incidence base has genus one; the rational unordered quotient has arithmetic rank two, not four.",
                },
                "proof_boundary": "Exact over the ordered genus-one incidence field; no rational rank-four P1 family is claimed.",
            },
            {
                "id": "R17_genus_one_bisections",
                "tuple": genus1_tuple,
                "transition": {
                    "operation": "quadratic_pullback_then_visible_anti_invariant_line",
                    "visible_rank_lower_bound": 18,
                    "anti_invariant_height": 16,
                    "base_genus": 1,
                },
                "proof_boundary": "The orthogonal character sum and its index-two half-sum extension are exact; further saturation, total anti-rank, and the full MW lattice are not determined.",
            },
            {
                "id": "R17_paired_rational_bisection_V4_composita",
                "tuple": paired_tuple,
                "transition": {
                    "operation": "two_independent_quadratic_pullbacks_then_V4_character_split",
                    "visible_rank_lower_bound": 19,
                    "new_height_block": [[24, 0], [0, 24]],
                    "base_genus": 1,
                    "pair_count": 765167640,
                },
                "proof_boundary": "These are pairs from the 39120-class rational norm-ten bisection atlas. All are geometrically connected genus-one V4 covers; rational points, product-character rank, and full saturation are not uniform consequences.",
            },
            {
                "id": "R17_norm12_simultaneous_split_specialization",
                "tuple": {
                    "L": {
                        "name": "specialized_section_lattice_at_t=1/25",
                        "rank": ">=18",
                        "determinant_absolute": "not determined",
                    },
                    "G": {
                        "name": "not a computed deck-action lattice record",
                        "order": "not applicable to the specialization quotient certificate",
                    },
                    "L_chi": "not determined",
                    "A_L": "not determined",
                    "q_L": "not determined",
                    "glue_subgroups": "not determined",
                },
                "transition": {
                    "operation": "simultaneous_specialization_split",
                    "parameter": "1/25",
                    "split_classes": [
                        "norm8-orbit-0f6b1",
                        "norm12-orbit-103b2",
                    ],
                    "certified_rank_lower_bound": 18,
                    "independent_escape_basis_labels": ["norm12-orbit-103b2"],
                },
                "classification_status": "success recorded, but no integral Gram, finite form, or glue mutation has yet been extracted",
                "proof_boundary": "Finite-quotient escape proves one new specialization direction modulo the generic MW17 subgroup; it is not a generic rank transfer or a completed lattice/glue tuple.",
            },
        ],
        "scope_boundary": {
            "Golay": "Golay here means the binary-octad determinant-720 construction and its source corridor. The ternary-Golay N(12A2) foundry backend is a large frame-family census, not one of the named transition corridors, and is not silently summarized as a single tuple.",
        },
        "aggregate_same_NS_clicks": {
            "edge_count": len(all_edges),
            "old_fibre_degree_histogram": {"2": len(all_edges)},
            "q_histogram": q_histogram,
            "q_4_or_6_count": q_histogram.get("4", 0) + q_histogram.get("6", 0),
            "discriminant_form_changes": 0,
            "nonunimodular_transports": 0,
        },
        "main_conclusion": {
            "positive": "Three operations recur: primitive-U change, maximal graph glue/complement, and character-primary saturation.",
            "negative": "No same-NS corridor click is detected by a change of finite discriminant form; a pure overlattice extension cannot annihilate an existing root.",
            "inversion_target": "Compute common-substructure p-neighbor/reglue data for frame classes, then enumerate only reversible regluings whose actual new cosets pass minimum at least four.",
            "missing_exact_experiment": "For each consecutive frame pair compute a common embedded lattice K, both quotient glue subgroups H_old,H_new, prime-factor indices, and root survival by glue coset.",
        },
        "negative_controls": [
            {
                "id": "wrong_H3_local_genus",
                "observation": "65 even rootless determinant-948 seeds form 19 classes but have the wrong 2-adic and 79-adic symbols.",
                "lesson": "Rank, determinant, parity, minimum and roots do not determine q_L or the genus.",
            },
            {
                "id": "Golay_rational_3I6_overglue",
                "observation": "An index-2 half-section and order-3 torsion give maximal index 6, changing determinant 720 to 20.",
                "glue_structure": "2-primary half-section plus 3-primary torsion",
                "lesson": "A simple rational equation can lie on the wrong integral NS overlattice.",
            },
            {
                "id": "E6_rank4_descent",
                "observation": "The ordered genus-one field has geometric split 2+2, while the rational unordered quotient has rank split 1+1.",
                "lesson": "Character rank must be computed after descent; unordered incidence does not preserve individual eigensections.",
            },
            {
                "id": "R17_single_C2_collision",
                "observation": "All 39120 rational bisection squareclasses are distinct.",
                "lesson": "Two anti-invariant gains require a V4 compositum or a new twist-rank mechanism, not two bisections on one quadratic cover.",
            },
        ],
    }

    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = args.output if args.output.is_absolute() else ROOT / args.output
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise SystemExit(f"stale or missing artifact: {output}")
        print("PASS integral rank-transfer glue census")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(serialized)
    print(output.relative_to(ROOT))


if __name__ == "__main__":
    main()
