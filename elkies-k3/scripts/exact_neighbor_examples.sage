#!/usr/bin/env sage
"""Pinned section/component inputs for exact-neighbor regression certificates."""

import csv
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, matrix, vector


HERE = Path(__file__).resolve().parent
BASE = HERE.parent
load(str(HERE / "exact_neighbor_engine.sage"))


def load_frame(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def q80_first_q4_example():
    """Return q80's first globally-nef q=4 datum and its supplied walls."""
    frame = load_frame(BASE / "data/fibrations/kumar_q80_e6_d5_a3_mw3_frame.txt")
    ns = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -frame)
    old_fiber = vector(ZZ, [1, 0] + [0] * 17)
    zero = vector(ZZ, [-1, 1] + [0] * 17)
    with (BASE / "data/fibrations/kumar_q80_to_rootless_path.tsv").open() as handle:
        first = next(csv.DictReader(handle, delimiter="\t"))
    if tuple(first[key] for key in ("step", "q", "a", "b")) != ("1", "4", "2", "2"):
        raise AssertionError("pinned q80 first neighbor changed")
    divisor = vector(ZZ, [2, 2] + [ZZ(value) for value in first["v"].split(",")])
    walls, _ = component_walls(frame, old_fiber, include_zero=zero)
    return {
        "name": "q80-first-q4",
        "ns": ns,
        "old_fiber": old_fiber,
        "divisor": divisor,
        "curves": walls,
        "expected_root_data": (13, 164, 20),
        "proof_metadata": {
            "global_nef_certificate": "elkies-k3/scripts/analyze_q80_rootless_first_neighbor.sage",
            "scope": "The companion performs the exact MW/CVP and degree-two bisection checks.",
        },
    }


def h3_d13_q24_example():
    """Return the first rank-growing H3 q=24 degree-two chamber datum."""
    frame = load_frame(BASE / "data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt")
    ns = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -frame)
    old_fiber = vector(ZZ, [1, 0] + [0] * 17)
    zero = vector(ZZ, [-1, 1] + [0] * 17)
    divisor = vector(ZZ, (
        12, 2,
        0, 5, 0, 1, 2, 1, 2, 2, 2, 2, 4, 8, 2, 0, -1, 1, 1,
    ))
    highest = (2, 2, 1, 1, 2, 1, 2, 2, 2, 2, 2, 2, 2)
    curves = [("O", zero)] + [
        (f"D13_{index + 1}", vector(ZZ, [0, 0] + [
            -ZZ(index == node) for node in range(17)
        ]))
        for index in range(13)
    ]
    curves.append((
        "D13_affine",
        old_fiber + vector(ZZ, [0, 0] + list(highest) + [0] * 4),
    ))
    return {
        "name": "h3-d13-q24",
        "ns": ns,
        "old_fiber": old_fiber,
        "divisor": divisor,
        "curves": tuple(curves),
        "expected_root_data": (12, 264, 4),
        "proof_metadata": {
            "global_nef_certificate": "elkies-k3/scripts/analyze_h3_d13_q4_chamber.sage",
            "scope": "The companion completes the section quotient and bisection parity checks.",
        },
    }


def run_example(example):
    """Run the engine and package the result with its exact proof metadata."""
    result = degree_two_neighbor(
        example["ns"], example["divisor"], example["old_fiber"], example["curves"]
    )
    if result["child_root_data"] != example["expected_root_data"]:
        raise AssertionError(f"{example['name']} child root data changed")
    certificate = neighbor_certificate(
        example["ns"], example["old_fiber"], example["curves"], result,
        proof_metadata=example["proof_metadata"],
    )
    return result, certificate
