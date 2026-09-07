#!/usr/bin/env sage-python
"""Certify maximal-graph core-genus inversion and its first core gate.

status: ACTIVE_PROOF
claim: A full rank-two graph completion W of K+C forces and is forced at
  finite-form level by q_K = q_W orthogonal_sum (-q_C).  The exact recorded
  bridge presentations satisfy this identity, and the minimum-two gate rejects
  the complete held-out E6 core shell before any graph-glue enumeration.
inputs: artifacts/generated-results/
  elkies-k3-integral-rank-transfer-bridge-reglue-v1.json,
  artifacts/generated-results/
  elkies-k3-integral-rank-transfer-theta-convolution-v1.json,
  artifacts/generated-results/
  elkies-k3-e6-rank4-det78-prospective-bridge-predictor-v1.json
outputs: artifacts/generated-results/
  elkies-k3-integral-rank-transfer-core-generation-v1.json
supersedes/superseded-by: none
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from sage.all import (
    Genus,
    QQ,
    ZZ,
    block_diagonal_matrix,
    identity_matrix,
    lcm,
    matrix,
    pari,
    vector,
)
from sage.quadratic_forms.genera.genus import genera


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
BRIDGES = GENERATED / "elkies-k3-integral-rank-transfer-bridge-reglue-v1.json"
THETA = GENERATED / "elkies-k3-integral-rank-transfer-theta-convolution-v1.json"
E6 = GENERATED / "elkies-k3-e6-rank4-det78-prospective-bridge-predictor-v1.json"
OUTPUT = GENERATED / "elkies-k3-integral-rank-transfer-core-generation-v1.json"


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rational_rows(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def finite_form_key(discriminant_form):
    normal = discriminant_form.normal_form()
    return {
        "invariants": list(map(int, normal.invariants())),
        "quadratic_gram": rational_rows(normal.gram_matrix_quadratic()),
        "value_module": str(normal.value_module_qf()),
    }


def discriminant_form_key(gram):
    """Canonical key for the finite quadratic form of an even Gram matrix."""

    return finite_form_key(Genus(gram).discriminant_form())


def glued_frame(core, bridge, glue_vector):
    """Construct the integral overlattice selected by one glue generator."""

    split = block_diagonal_matrix(core, bridge)
    denominator = lcm(value.denominator() for value in glue_vector)
    generators = (
        denominator * identity_matrix(QQ, split.nrows())
    ).stack(matrix(QQ, [denominator * glue_vector])).change_ring(ZZ)
    basis = generators.row_module(ZZ).basis_matrix().change_ring(QQ) / denominator
    candidate = basis * split * basis.transpose()
    assert all(value in ZZ for value in candidate.list())
    candidate = candidate.change_ring(ZZ)
    assert not any(value % 2 for value in candidate.diagonal())
    return candidate


def minimum_norm(gram):
    """Exact minimum for the few rank-15 cores in the positive corpus."""

    bound = 2
    while True:
        result = pari(gram).qfminim(bound)
        if int(result[0]):
            vectors = matrix(ZZ, result[2].sage()).columns()
            return min(
                int(vector(ZZ, item) * gram * vector(ZZ, item))
                for item in vectors
            )
        bound += 2


def split_record(edge, side):
    core = matrix(ZZ, edge["core"]["gram"])
    presentation = edge[side]
    bridge = matrix(ZZ, presentation["bridge_gram"])
    generators = presentation["glue_generators"]
    assert len(generators) == 1
    glue = vector(
        QQ,
        [QQ(value) for value in generators[0]["K_plus_C_dual_coordinates"]],
    )
    frame = glued_frame(core, bridge, glue)

    bridge_determinant = abs(int(bridge.det()))
    core_determinant = abs(int(core.det()))
    frame_determinant = abs(int(frame.det()))
    glue_order = int(generators[0]["order"])
    assert glue_order == bridge_determinant
    assert int(presentation["K_plus_C_index_in_W"]) == bridge_determinant
    assert core_determinant == frame_determinant * bridge_determinant

    # The right side has discriminant form q_W orthogonal_sum (-q_C).
    generated_core_form = discriminant_form_key(
        block_diagonal_matrix(frame, -bridge)
    )
    observed_core_form = discriminant_form_key(core)
    assert generated_core_form == observed_core_form

    return {
        "side": side.removesuffix("_frame"),
        "core_determinant": core_determinant,
        "frame_determinant": frame_determinant,
        "bridge_determinant": bridge_determinant,
        "glue_order": glue_order,
        "determinant_identity": True,
        "finite_quadratic_form_identity": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    bridge_data = json.loads(BRIDGES.read_text())
    theta_data = json.loads(THETA.read_text())
    e6_data = json.loads(E6.read_text())
    assert bridge_data["status"] == "PASS_EXACT_BRIDGE_REGLUE_CERTIFICATES"
    assert theta_data["status"] == (
        "PASS_EXACT_THETA_CONVOLUTION_ZERO_SUPPORT_ENUMERATOR"
    )
    assert e6_data["status"] == (
        "PASS_BLIND_DET78_PROSPECTIVE_BRIDGE_PREDICTOR_NEGATIVE_CONTROL"
    )

    presentations = []
    for edge in bridge_data["edges"]:
        for side in ("old_frame", "new_frame"):
            presentations.append(
                {
                    "corridor": edge["corridor"],
                    "edge_index": int(edge["edge_index"]),
                    **split_record(edge, side),
                }
            )
    assert len(presentations) == 84

    theta_by_corridor = {
        row["corridor"]: row for row in theta_data["corridors"]
    }
    terminal_edges = {
        row["corridor"]: row
        for row in bridge_data["edges"]
        if int(row["target_root_rank"]) == 0
    }
    assert set(theta_by_corridor) == set(terminal_edges)
    positive_cores = []
    for corridor in sorted(terminal_edges):
        edge = terminal_edges[corridor]
        theta = theta_by_corridor[corridor]
        core = matrix(ZZ, edge["core"]["gram"])
        bridge = matrix(ZZ, edge["new_frame"]["bridge_gram"])
        glue = vector(
            QQ,
            [
                QQ(value)
                for value in edge["new_frame"]["glue_generators"][0][
                    "K_plus_C_dual_coordinates"
                ]
            ],
        )
        frame = glued_frame(core, bridge, glue)
        generated_form = discriminant_form_key(
            block_diagonal_matrix(frame, -bridge)
        )
        determinant_genera = genera(
            (15, 0), abs(int(core.det())), even=True
        )
        matching_genera = [
            genus
            for genus in determinant_genera
            if finite_form_key(genus.discriminant_form()) == generated_form
        ]
        assert len(matching_genera) == 1
        assert [str(value) for value in matching_genera[0].local_symbols()] == [
            str(value) for value in Genus(core).local_symbols()
        ]
        core_minimum = minimum_norm(core)
        assert core_minimum > 2
        class_count = len(theta["classes"])
        rootless_count = sum(row["predicted_rootless"] for row in theta["classes"])
        assert rootless_count > 0
        positive_cores.append(
            {
                "corridor": corridor,
                "rank": int(core.nrows()),
                "determinant": abs(int(core.det())),
                "minimum": core_minimum,
                "discriminant_form": discriminant_form_key(core),
                "even_rank15_genera_at_determinant": len(determinant_genera),
                "genera_matching_generated_discriminant_form": len(matching_genera),
                "generated_core_genus_local_symbols": [
                    str(value) for value in matching_genera[0].local_symbols()
                ],
                "theta_signature_through_norm_two": theta["core_theta_profile"],
                "bridge_determinant": int(theta["cyclic_glue_order"]),
                "enumerated_binary_bridge_classes": class_count,
                "zero_support_binary_bridge_classes": rootless_count,
                "admits_zero_support_rank_two_completion": True,
            }
        )
    assert len({row["determinant"] for row in positive_cores}) == 4

    ranked = e6_data["ranked_candidates"]
    assert len(ranked) == 277
    rejected = [row for row in ranked if int(row["core_minimum"]) == 2]
    assert len(rejected) == len(ranked)
    rejected_by_order = Counter(
        int(row["glue_group_invariants"][0]) for row in rejected
    )

    payload = {
        "schema": "elkies-k3.integral-rank-transfer-core-generation.v1",
        "status": "PASS_EXACT_MAXIMAL_GRAPH_CORE_GENUS_INVERSION",
        "inputs": {
            relative(BRIDGES): digest(BRIDGES),
            relative(THETA): digest(THETA),
            relative(E6): digest(E6),
        },
        "theorem": {
            "hypotheses": (
                "K, C and W are even nondegenerate lattices; W is the "
                "overlattice of K+C defined by an isotropic graph H whose "
                "projection H->A_C is an isomorphism."
            ),
            "forward": (
                "The image B in A_K is nondegenerate and anti-isometric to "
                "A_C. Hence q_K is q_W orthogonal_sum (-q_C), and "
                "det(K)=det(W)*det(C)."
            ),
            "converse": (
                "A splitting q_K = q_W orthogonal_sum (-q_C) gives the "
                "diagonal isotropic graph between the -q_C summand and A_C; "
                "its overlattice has discriminant form q_W."
            ),
            "rootless_criterion": (
                "For a fixed binary bridge C, K admits a rootless maximal-graph "
                "completion exactly when K and C are rootless and some "
                "anti-isometry A_C->B in A_K has zero theta convolution on "
                "its entire graph."
            ),
            "complete_core_descriptor": (
                "For any declared finite bridge universe, rank plus the finite "
                "quadratic form and all discriminant-coset theta coefficients "
                "through norm two determine the complete zero-support "
                "completion spectrum of K."
            ),
        },
        "core_first_algorithm": [
            "Set q_W=-q_NS for the desired Picard-rank-19 frame genus.",
            "Enumerate a declared finite universe of positive even binary bridges C.",
            "Generate the rank-15 core genus by q_K=q_W orthogonal_sum (-q_C) and det(K)=det(W)*det(C).",
            "Enumerate and mass-close lattice classes K in that genus; reject theta_K(0,2)>0 before any glue graph.",
            "Deduplicate survivors by isomorphism of their theta-decorated discriminant forms through norm two.",
            "Enumerate anti-isometric bridge summands and retain exactly the zero-support graphs; construct rank-17 frames only for survivors.",
        ],
        "recorded_click_replay": {
            "edges": int(bridge_data["aggregate"]["edge_count"]),
            "maximal_graph_presentations_checked": len(presentations),
            "all_determinant_identities": True,
            "all_finite_quadratic_form_identities": True,
            "corridor_presentation_counts": {
                key: value
                for key, value in sorted(
                    Counter(row["corridor"] for row in presentations).items()
                )
            },
            "presentations": presentations,
        },
        "observed_positive_core_classes": positive_cores,
        "held_out_e6_early_gate": {
            "source_shell_complete_by_declared_weyl_enumeration": (
                int(e6_data["shell"]["dominant_orbits"]) == 280
                and int(e6_data["shell"]["nonprimitive_orbits"]) == 3
            ),
            "target_j2_genus_mass_closed": bool(
                e6_data["truth_set"]["mass_closed"]
            ),
            "primitive_rank15_cores": len(ranked),
            "rejected_by_core_minimum_two": len(rejected),
            "survivors_before_graph_enumeration": len(ranked) - len(rejected),
            "rejected_by_cyclic_order": {
                str(key): value for key, value in sorted(rejected_by_order.items())
            },
            "interpretation": (
                "The necessary core-rootlessness gate rejects the full declared "
                "E6 source shell without computing a bridge theta table, an "
                "isotropic subgroup, or a rank-17 child root system."
            ),
        },
        "proof_boundary": {
            "proved": (
                "The maximal-graph finite-form equivalence, the exact decorated-"
                "theta completion criterion, all 84 recorded splitting identities, "
                "and the 277/277 E6 minimum-two rejection."
            ),
            "not_proved": (
                "No unbounded bridge-determinant cutoff, uniform algorithm for "
                "enumerating or mass-closing every rank-15 genus, complexity "
                "bound, speedup theorem, or equation-level fibration lift is claimed."
            ),
        },
        "reproduce": (
            "sage -python elkies-k3/scripts/"
            "certify_integral_rank_transfer_core_generation.sage --check"
        ),
    }

    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    if arguments.check:
        if not output.exists():
            raise SystemExit(f"missing artifact: {output}")
        if output.read_text() != encoded:
            raise SystemExit(f"stale artifact: {output}")
        print("PASS maximal-graph core-genus inversion")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
    print(relative(output))


if __name__ == "__main__":
    main()
