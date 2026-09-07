#!/usr/bin/env python3
"""Construct the exact point set for the bounded wgxli elementary mutation.

The declared search is complete for one common elementary shear from the
retained signed permutation.  It fixes the stable denominator anchors, uses
columns of l1 norm at most two, requires determinant one, bounds the changed
canonical height by 1.25 times its displayed value on every fibre, and retains
only a joint Gram-objective improvement of at least one percent.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import time


ROOT = Path(__file__).resolve().parents[2]
LINEAGE = (
    ROOT / "artifacts/generated-results/elliptic-curves"
    / "icarm_wgxli_rank17_lineage_v1.json"
)
REBASING = (
    ROOT / "artifacts/generated-results/elliptic-curves"
    / "icarm_wgxli_rank17_signed_permutation_rebasing_v1.json"
)
OUTPUT = (
    ROOT / "artifacts/generated-results/elliptic-curves"
    / "icarm_wgxli_rank17_mutation_p4_minus_p1_v1.json"
)
TARGETS = (351, 356, 376, 377, 385)
HEIGHT_INFLATION_BOUND = 1.25
MINIMUM_JOINT_IMPROVEMENT = 0.01


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=LINEAGE)
    parser.add_argument("--rebasing", type=Path, default=REBASING)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    return parser.parse_args()


def sha256_bytes(raw):
    return hashlib.sha256(raw).hexdigest()


def rational_text(value):
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def height_matrices(fibres):
    program = ["default(realprecision,80);"]
    for fibre in fibres:
        curve_id = int(fibre["curve_id"])
        A, B = fibre["short_model"]
        points = ",".join(
            f"[{x_value},{y_value}]"
            for x_value, y_value in fibre["short_points_first_17"]
        )
        program.extend((
            f"E=ellinit([0,0,0,{A},{B}]);P=[{points}];H=ellheightmatrix(E,P);",
            f'print("BEGIN|{curve_id}");',
            "for(i=1,17,for(j=1,17,if(j>1,print1(\"|\"));print1(H[i,j]));print());",
        ))
    program.append('print("PARI|",version());')
    completed = subprocess.run(
        ["gp", "-q"], input="\n".join(program) + "\n", text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, timeout=120,
    )
    if completed.stderr.strip():
        raise RuntimeError(completed.stderr.strip())
    matrices = {}
    current = None
    pari_version = None
    for line in completed.stdout.splitlines():
        line = line.strip()
        if line.startswith("BEGIN|"):
            current = int(line.split("|", 1)[1])
            matrices[current] = []
        elif line.startswith("PARI|"):
            pari_version = line.split("|", 1)[1]
        elif line:
            matrices[current].append([float(value) for value in line.split("|")])
    return matrices, pari_version


def matmul(left, right):
    return [
        [sum(left[i][k] * right[k][j] for k in range(len(right))) for j in range(len(right[0]))]
        for i in range(len(left))
    ]


def transpose(matrix):
    return [list(column) for column in zip(*matrix)]


def transform_gram(gram, transform):
    return matmul(transpose(transform), matmul(gram, transform))


def fit_residual_squared(left, right):
    dot = sum(x * y for left_row, right_row in zip(left, right) for x, y in zip(left_row, right_row))
    left_norm = sum(x * x for row in left for x in row)
    right_norm = sum(x * x for row in right for x in row)
    return max(0.0, 1.0 - dot * dot / (left_norm * right_norm))


def joint_objective(matrices):
    answer = 0.0
    for left_index, left in enumerate(TARGETS):
        for right in TARGETS[left_index + 1 :]:
            answer += fit_residual_squared(matrices[left], matrices[right])
    return answer


def identity_matrix(size):
    return [[int(row == column) for column in range(size)] for row in range(size)]


def add_short(left, right, A):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = map(Fraction, left)
    x2, y2 = map(Fraction, right)
    A = Fraction(A)
    if x1 == x2:
        if y2 == -y1:
            return None
        slope = (3 * x1 * x1 + A) / (2 * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope * slope - x1 - x2
    y3 = -y1 + slope * (x1 - x3)
    return x3, y3


def main():
    arguments = parse_args()
    started = time.monotonic()
    lineage_raw = arguments.input.read_bytes()
    rebasing_raw = arguments.rebasing.read_bytes()
    lineage = json.loads(lineage_raw)
    rebasing = json.loads(rebasing_raw)
    fibres = lineage["rootless_k3_interpolation_input"]["fibres"]
    if [int(fibre["curve_id"]) for fibre in fibres] != list(TARGETS):
        raise AssertionError("unexpected target order")
    retained = rebasing["retained_signed_permutation_candidates"]
    if len(retained) != 1 or any(
        row["signs"] != "+" * 17
        or row["permutation_new_label_to_old_point"] != list(range(1, 18))
        for row in retained[0].values()
    ):
        raise AssertionError("the bounded mutation must start from the literal retained basis")

    matrices, pari_version = height_matrices(fibres)
    base_objective = joint_objective(matrices)
    anchors = {
        int(record["label"]) - 1
        for record in rebasing["bounded_permutation_search"]["stable_denominator_anchors"]
    }
    search_records = []
    for changed in range(17):
        if changed in anchors:
            continue
        for added in range(17):
            if changed == added:
                continue
            for sign in (-1, 1):
                transform = identity_matrix(17)
                transform[added][changed] = sign
                transformed = {
                    curve_id: transform_gram(matrices[curve_id], transform)
                    for curve_id in TARGETS
                }
                ratios = [
                    transformed[curve_id][changed][changed]
                    / matrices[curve_id][changed][changed]
                    for curve_id in TARGETS
                ]
                objective = joint_objective(transformed)
                record = {
                    "changed_label": changed + 1,
                    "new_point_word": [
                        [changed + 1, 1],
                        [added + 1, sign],
                    ],
                    "determinant": 1,
                    "maximum_column_l1_norm": 2,
                    "height_ratios_by_curve": {
                        str(curve_id): ratio for curve_id, ratio in zip(TARGETS, ratios)
                    },
                    "maximum_height_ratio": max(ratios),
                    "joint_objective": objective,
                    "relative_improvement": (base_objective - objective) / base_objective,
                }
                if (
                    max(ratios) <= HEIGHT_INFLATION_BOUND
                    and record["relative_improvement"] >= MINIMUM_JOINT_IMPROVEMENT
                ):
                    search_records.append(record)
    search_records.sort(key=lambda row: (row["joint_objective"], row["changed_label"], row["new_point_word"]))
    if len(search_records) != 1 or search_records[0]["new_point_word"] != [[4, 1], [1, -1]]:
        raise AssertionError(f"bounded mutation proposal changed: {search_records}")

    proposal = search_records[0]
    transformed_fibres = []
    for fibre in fibres:
        A, B = map(Fraction, fibre["short_model"])
        points = [tuple(map(Fraction, point)) for point in fibre["short_points_first_17"]]
        changed = add_short(points[3], (points[0][0], -points[0][1]), A)
        if changed is None:
            raise AssertionError("the proposed exact mutation became the origin")
        if changed[1] ** 2 != changed[0] ** 3 + A * changed[0] + B:
            raise AssertionError("exact mutated point missed its short curve")
        points[3] = changed
        transformed = dict(fibre)
        transformed["short_points_first_17"] = [
            [rational_text(x_value), rational_text(y_value)]
            for x_value, y_value in points
        ]
        transformed_fibres.append(transformed)

    interpolation = dict(lineage["rootless_k3_interpolation_input"])
    interpolation["fibres"] = transformed_fibres
    interpolation["hypothesis_only"] = True
    clean_test_primes = {}
    for prime in (17, 53):
        for fibre in transformed_fibres:
            A, B = map(Fraction, fibre["short_model"])
            discriminant = -16 * (4 * A**3 + 27 * B**2)
            if (
                discriminant.numerator % prime == 0
                or discriminant.denominator % prime == 0
            ):
                raise AssertionError(
                    f"mutated fibre {fibre['curve_id']} is not clean at {prime}"
                )
            if any(
                Fraction(value).denominator % prime == 0
                for point in fibre["short_points_first_17"]
                for value in point
            ):
                raise AssertionError(
                    f"mutated point coordinate has a denominator at {prime}"
                )
        clean_test_primes[str(prime)] = {
            "all_five_fibres_nonsingular": True,
            "all_85_point_coordinates_defined": True,
        }
    payload = {
        "schema": "icarm.wgxli-rank17-bounded-elementary-mutation.v1",
        "status": "PASS_EXACT_POINT_CONSTRUCTION_FOR_ONE_BOUNDED_MUTATION",
        "inputs": {
            str(arguments.input.relative_to(ROOT)): sha256_bytes(lineage_raw),
            str(arguments.rebasing.relative_to(ROOT)): sha256_bytes(rebasing_raw),
        },
        "software": {"pari_gp": pari_version},
        "declared_search_bound": {
            "starting_bases": "the unique retained signed permutation",
            "common_transform_across_fibres": True,
            "elementary_mutation_count": 1,
            "allowed_column": "e_i plus or minus e_j, i != j",
            "maximum_column_l1_norm": 2,
            "unimodular_determinant": 1,
            "stable_denominator_anchor_labels_fixed": sorted(index + 1 for index in anchors),
            "maximum_changed_height_ratio_on_every_fibre": HEIGHT_INFLATION_BOUND,
            "minimum_joint_gram_objective_improvement": MINIMUM_JOINT_IMPROVEMENT,
            "enumerated_transform_count": (17 - len(anchors)) * 16 * 2,
            "retained_transform_count": len(search_records),
        },
        "displayed_basis_joint_objective": base_objective,
        "retained_proposal": proposal,
        "exact_group_law": {
            "word": "P4-P1",
            "all_five_mutated_points_verified_on_their_short_curves": True,
            "transformation_determinant": 1,
        },
        "clean_modular_test_primes": clean_test_primes,
        "rootless_k3_interpolation_input": interpolation,
        "proof_boundary": (
            "The finite shear enumeration and rational group-law construction are exact "
            "inside the declared one-mutation bound. The Gram objective is numerical and "
            "only proposes this candidate; modular first-jet elimination is a separate gate."
        ),
        "runtime_seconds": time.monotonic() - started,
    }
    if arguments.check:
        expected = json.loads(arguments.output.read_text())
        expected.pop("runtime_seconds", None)
        payload.pop("runtime_seconds", None)
        if expected != payload:
            raise SystemExit("stale bounded-mutation artifact")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
        arguments.output.write_text(rendered)
        print(f"WGXLIMUTATION|output={arguments.output}|sha256={sha256_bytes(rendered.encode())}")
    print(
        "WGXLIMUTATION|enumerated=352|retained=1|word=P4-P1|"
        "status=PASS_EXACT_POINT_CONSTRUCTION_FOR_ONE_BOUNDED_MUTATION"
    )


if __name__ == "__main__":
    main()
