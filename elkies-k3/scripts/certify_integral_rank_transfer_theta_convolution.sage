#!/usr/bin/env sage-python
"""Certify the inverse theta-convolution root test on terminal bridge glues.

status: ACTIVE_PROOF
claim: On the complete fourteen-class terminal binary-bridge census for H3,
  Q80, NS0024, and Golay720, cached low-norm discriminant-coset theta tables
  for the rank-15 core and rank-2 bridge reproduce the exact root count of
  every admissible graph glue.  Their zero-support test selects exactly the
  five rootless bridge classes without constructing a rank-17 child lattice.
inputs: artifacts/generated-results/
  elkies-k3-integral-rank-transfer-bridge-reglue-v1.json,
  artifacts/generated-results/
  elkies-k3-integral-rank-transfer-bridge-predictor-benchmark-v1.json
outputs: artifacts/generated-results/
  elkies-k3-integral-rank-transfer-theta-convolution-v1.json
supersedes/superseded-by: none
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, lcm, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
BRIDGES = GENERATED / "elkies-k3-integral-rank-transfer-bridge-reglue-v1.json"
TERMINAL = (
    GENERATED
    / "elkies-k3-integral-rank-transfer-bridge-predictor-benchmark-v1.json"
)
OUTPUT = (
    GENERATED / "elkies-k3-integral-rank-transfer-theta-convolution-v1.json"
)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def fractional_part(value):
    value = QQ(value)
    return value - value.floor()


def discriminant_class(coordinates):
    return tuple(fractional_part(value) for value in coordinates)


def theta_profile(gram):
    """Return theta coefficients in all dual cosets through norm two."""

    inverse = gram.inverse()
    denominator = lcm(value.denominator() for value in inverse.list())
    scaled_inverse = (denominator * inverse).change_ring(ZZ)
    enumeration = pari(scaled_inverse).qfminim(2 * denominator)
    representatives = matrix(ZZ, enumeration[2].sage()).columns()

    profile = defaultdict(Counter)
    zero = vector(QQ, [0] * gram.nrows())
    profile[discriminant_class(zero)][QQ(0)] = 1
    for representative in representatives:
        for dual_pairing in (
            vector(ZZ, representative),
            -vector(ZZ, representative),
        ):
            dual_vector = dual_pairing * inverse
            norm = dual_pairing * inverse * dual_pairing
            assert 0 < norm <= 2
            profile[discriminant_class(dual_vector)][norm] += 1

    vector_count = sum(sum(row.values()) for row in profile.values())
    assert vector_count == int(enumeration[0]) + 1
    return profile


def canonical_profile(profile):
    return [
        {
            "discriminant_class": [str(value) for value in residue],
            "norm_counts": {
                str(norm): count for norm, count in sorted(counts.items())
            },
        }
        for residue, counts in sorted(profile.items())
    ]


def profile_summary(profile):
    encoded = json.dumps(
        canonical_profile(profile), sort_keys=True, separators=(",", ":")
    ).encode()
    return {
        "represented_discriminant_classes": len(profile),
        "dual_vectors_of_norm_at_most_two_including_zero": sum(
            sum(row.values()) for row in profile.values()
        ),
        "theta_table_sha256": hashlib.sha256(encoded).hexdigest(),
    }


def bridge_generator(gram):
    generator = vector(ZZ, [1, 0]) * gram.inverse()
    if generator in ZZ**2:
        generator = vector(ZZ, [0, 1]) * gram.inverse()
    return generator


def reduced_even_binary_forms(determinant):
    """Enumerate every Minkowski-reduced positive even binary form."""

    result = []
    for half_left in range(1, determinant + 1):
        for half_right in range(half_left, determinant + 1):
            for off_diagonal in range(half_left + 1):
                if 4 * half_left * half_right - off_diagonal**2 != determinant:
                    continue
                result.append(
                    matrix(
                        ZZ,
                        [
                            [2 * half_left, off_diagonal],
                            [off_diagonal, 2 * half_right],
                        ],
                    )
                )
    return result


def gram_key(gram):
    return tuple(tuple(int(value) for value in row) for row in gram.rows())


def convolved_root_profile(
    core_profile,
    bridge_profile,
    core_generator,
    bridge_generator_value,
    core_multiplier,
    order,
):
    """Count norm-two vectors in every class of one cyclic graph glue."""

    result = {}
    for label in range(order):
        core_class = discriminant_class(
            label * core_multiplier * core_generator
        )
        bridge_class = discriminant_class(label * bridge_generator_value)
        count = 0
        for norm, multiplicity in core_profile.get(core_class, {}).items():
            count += multiplicity * bridge_profile.get(bridge_class, {}).get(
                QQ(2) - norm, 0
            )
        if count:
            result[label] = count
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    bridge_data = json.loads(BRIDGES.read_text())
    terminal_data = json.loads(TERMINAL.read_text())
    assert bridge_data["status"] == "PASS_EXACT_BRIDGE_REGLUE_CERTIFICATES"
    assert terminal_data["status"].startswith("PASS_EXACT_RETROSPECTIVE")

    terminal_edges = {
        row["corridor"]: row
        for row in bridge_data["edges"]
        if int(row["target_root_rank"]) == 0
    }
    declared = terminal_data["terminal_binary_bridge_census"]["corridors"]
    assert set(terminal_edges) == {row["corridor"] for row in declared}

    corridors = []
    total_bridge_classes = 0
    total_graphs = 0
    total_rootless = 0
    for terminal in declared:
        corridor = terminal["corridor"]
        edge = terminal_edges[corridor]
        order = int(terminal["cyclic_glue_order"])
        assert order == int(edge["bridge_replacement"]["common_cyclic_glue_order"])
        assert ZZ(order).is_prime()

        core = matrix(ZZ, edge["core"]["gram"])
        assert int(pari(core).qfminim(2)[0]) == 0
        stored_generator = edge["new_frame"]["glue_generators"]
        assert len(stored_generator) == 1
        split_generator = vector(
            QQ,
            [
                QQ(value)
                for value in stored_generator[0]["K_plus_C_dual_coordinates"]
            ],
        )
        core_generator = split_generator[:-2]
        assert order * core_generator in ZZ ** core.nrows()
        assert core_generator not in ZZ ** core.nrows()

        core_profile = theta_profile(core)
        # Enumerate the complete bridge universe from the cyclic order.  The
        # earlier census is retained only as an independent outcome label.
        generated_bridges = reduced_even_binary_forms(order)
        declared_by_gram = {
            gram_key(matrix(ZZ, row["bridge_gram"])): row
            for row in terminal["classes"]
        }
        assert {gram_key(bridge) for bridge in generated_bridges} == set(
            declared_by_gram
        )

        classes = []
        for bridge_index, bridge in enumerate(generated_bridges, 1):
            generator = bridge_generator(bridge)
            assert order * generator in ZZ**2
            assert generator not in ZZ**2
            bridge_profile = theta_profile(bridge)

            # This is the constructive step: derive graph labels only from
            # the two finite quadratic forms, before reading child outcomes.
            admissible = []
            for multiplier in range(1, order):
                glue_norm = (
                    (multiplier * core_generator)
                    * core
                    * (multiplier * core_generator)
                    + generator * bridge * generator
                )
                if glue_norm in ZZ and ZZ(glue_norm) % 2 == 0:
                    admissible.append((multiplier, int(glue_norm)))
            assert len(admissible) == 2

            predictions = []
            for multiplier, glue_norm in admissible:
                root_profile = convolved_root_profile(
                    core_profile,
                    bridge_profile,
                    core_generator,
                    generator,
                    multiplier,
                    order,
                )
                predicted = sum(root_profile.values())
                predictions.append(
                    {
                        "core_glue_multiplier": multiplier,
                        "isotropic_generator_norm": glue_norm,
                        "nonzero_root_cosets": {
                            str(label): count
                            for label, count in sorted(root_profile.items())
                        },
                        "predicted_signed_root_count": predicted,
                        "zero_support_accepts": predicted == 0,
                    }
                )

            # Only now use the earlier child-construction census as a label.
            declared_class = declared_by_gram[gram_key(bridge)]
            stored = {
                int(row["core_glue_multiplier"]): int(row["signed_root_count"])
                for row in declared_class["admissible_oriented_graph_labels"]
            }
            assert set(stored) == {
                row["core_glue_multiplier"] for row in predictions
            }
            graphs = []
            for prediction in predictions:
                expected = stored[prediction["core_glue_multiplier"]]
                assert prediction["predicted_signed_root_count"] == expected
                graphs.append(
                    {**prediction, "stored_child_signed_root_count": expected}
                )

            assert len(graphs) == 2
            assert graphs[0]["predicted_signed_root_count"] == graphs[1][
                "predicted_signed_root_count"
            ]
            predicted_rootless = graphs[0]["zero_support_accepts"]
            assert predicted_rootless == bool(declared_class["rootless"])
            classes.append(
                {
                    "bridge_class_index": bridge_index,
                    "bridge_gram": [
                        [int(value) for value in row] for row in bridge.rows()
                    ],
                    "bridge_theta_profile": profile_summary(bridge_profile),
                    "admissible_oriented_graphs": graphs,
                    "predicted_rootless": predicted_rootless,
                }
            )

        corridor_rootless = sum(row["predicted_rootless"] for row in classes)
        assert corridor_rootless == int(terminal["rootless_binary_bridge_classes"])
        total_bridge_classes += len(classes)
        total_graphs += sum(len(row["admissible_oriented_graphs"]) for row in classes)
        total_rootless += corridor_rootless
        corridors.append(
            {
                "corridor": corridor,
                "terminal_edge_index": int(edge["edge_index"]),
                "cyclic_glue_order": order,
                "core_theta_profile": profile_summary(core_profile),
                "bridge_class_count": len(classes),
                "rootless_bridge_class_count": corridor_rootless,
                "classes": classes,
            }
        )

    assert total_bridge_classes == 14
    assert total_graphs == 28
    assert total_rootless == 5
    payload = {
        "schema": "elkies-k3.integral-rank-transfer-theta-convolution.v1",
        "status": "PASS_EXACT_THETA_CONVOLUTION_ZERO_SUPPORT_ENUMERATOR",
        "inputs": {
            relative(BRIDGES): digest(BRIDGES),
            relative(TERMINAL): digest(TERMINAL),
        },
        "theorem": {
            "theta_coefficient": (
                "theta_K(a,nu)=#{x in K^dual: x mod K=a and x^2=nu}"
            ),
            "convolution": (
                "rho_KC(a,b)=sum_nu theta_K(a,nu)*theta_C(b,2-nu)"
            ),
            "zero_support_gate": (
                "The overlattice selected by isotropic H is rootless iff "
                "rho_KC vanishes on every element of H."
            ),
        },
        "algorithm_boundary": {
            "does_before_child_construction": (
                "Enumerates every Minkowski-reduced positive even binary bridge "
                "of the required determinant, derives every oriented graph "
                "multiplier from finite-form isotropy, caches every dual-coset "
                "theta coefficient through norm two, convolves the tables, and "
                "retains precisely the zero-support graphs."
            ),
            "comparison_only": (
                "Stored child root counts are read only after each convolution "
                "and are used as exact regression labels."
            ),
            "not_proved": (
                "A speedup, an automatic core-generation rule, or completeness "
                "outside the fourteen declared fixed-core bridge classes."
            ),
        },
        "aggregate": {
            "corridor_count": len(corridors),
            "bridge_class_count": total_bridge_classes,
            "oriented_graph_count": total_graphs,
            "exact_root_count_matches": total_graphs,
            "rootless_bridge_classes_selected": total_rootless,
            "constructed_rank17_children_during_prediction": 0,
        },
        "corridors": corridors,
        "reproduce": (
            "sage -python elkies-k3/scripts/"
            "certify_integral_rank_transfer_theta_convolution.sage --check"
        ),
    }

    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit(f"stale or missing artifact: {output}")
        print("PASS exact theta-convolution zero-support enumerator")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
    print(relative(output))


if __name__ == "__main__":
    main()
