#!/usr/bin/env sage-python
"""Catalogue all paired bisection bases with an evident rational point.

The complete equation batch contains 104 quadratic covers with square leading
coefficient and 21 with square constant coefficient, with one cover in both
sets.  Pairing within either set gives 5,566 distinct genus-one V4 bases with
a rational point over t=infinity or t=0.

For every such pair this script computes the third-quotient binary quartic,
its classical invariants, a global minimal Jacobian, conductor and global root
number.  The output is a complete exact arithmetic catalogue; Mordell--Weil
point searches and independence certificates are a separate resumable stage.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import reduce
from hashlib import sha256
from itertools import combinations
import json
from math import gcd
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, lcm


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-2026-immediate-point-pair-catalogue-full.json"
EXPECTED_INPUT_SCHEMA = "elkies-k3.bisection-extension-input.v1"
OUTPUT_SCHEMA = "elkies-k3.elkies-2026-immediate-point-pair-catalogue.v1"


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def rational_text(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def primitive_quadratic(coefficients):
    values = tuple(QQ(value) for value in coefficients)
    denominator = lcm([value.denominator() for value in values])
    integers = [ZZ(value * denominator) for value in values]
    content = reduce(gcd, (abs(int(value)) for value in integers))
    integers = tuple(value // content for value in integers)
    if integers[2] < 0:
        integers = tuple(-value for value in integers)
    return integers


def integer_bits(value) -> int:
    return max(1, abs(int(ZZ(value))).bit_length())


def pair_key(left, right):
    return tuple(sorted((int(left["record"]["lattice_orbit_mask"]), int(right["record"]["lattice_orbit_mask"]))))


def cover_entry(record):
    actual = tuple(QQ(value) for value in record["residual_chord"]["q_coefficients"])
    primitive = primitive_quadratic(actual)
    scale = actual[0] / primitive[0]
    if not scale.is_square():
        raise ArithmeticError("an immediate-point cover did not normalize by a rational square")
    square_rescaling = QQ(scale.sqrt())
    if actual != tuple(scale * value for value in primitive):
        raise ArithmeticError("quadratic normalization is not scalar")
    return {
        "record": record,
        "q": primitive,
        "square_rescaling": square_rescaling,
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

batch = json.loads(args.input.read_text())
if batch.get("schema") != EXPECTED_INPUT_SCHEMA:
    raise ValueError(f"unexpected input schema: {batch.get('schema')!r}")

at_zero = []
at_infinity = []
for record in batch["bisections"]:
    actual = tuple(QQ(value) for value in record["residual_chord"]["q_coefficients"])
    if actual[0].is_square() or actual[2].is_square():
        entry = cover_entry(record)
        if actual[0].is_square():
            at_zero.append(entry)
        if actual[2].is_square():
            at_infinity.append(entry)

if len(at_zero) != 21 or len(at_infinity) != 104:
    raise ArithmeticError(
        f"immediate-point cover counts changed: zero={len(at_zero)}, infinity={len(at_infinity)}"
    )

pairs = {}
for point_label, group in (("zero", at_zero), ("infinity", at_infinity)):
    for left, right in combinations(group, 2):
        key = pair_key(left, right)
        if key not in pairs:
            pairs[key] = {"left": left, "right": right, "points": []}
        pairs[key]["points"].append(point_label)
if len(pairs) != 5566:
    raise ArithmeticError(f"immediate-point pair count changed: {len(pairs)}")

R = PolynomialRing(QQ, "t")
rows = []
root_numbers = Counter()
point_sources = Counter()
for key, pair in pairs.items():
    left, right = sorted(
        (pair["left"], pair["right"]),
        key=lambda entry: int(entry["record"]["lattice_orbit_mask"]),
    )
    left_record = left["record"]
    right_record = right["record"]
    q_left = R(list(left["q"]))
    q_right = R(list(right["q"]))
    if not q_left.is_irreducible() or not q_right.is_irreducible() or q_left.gcd(q_right) != 1:
        raise ArithmeticError(f"pair {key} does not have four distinct branch points")

    quartic = q_left * q_right
    e, d, c, b, a = (quartic[index] for index in range(5))
    invariant_i = 12 * a * e - 3 * b * d + c**2
    invariant_j = (
        72 * a * c * e
        + 9 * b * c * d
        - 27 * a * d**2
        - 27 * b**2 * e
        - 2 * c**3
    )
    raw_jacobian = EllipticCurve(QQ, [0, 0, 0, -27 * invariant_i, -27 * invariant_j])
    minimal = raw_jacobian.global_minimal_model()
    minimal_model = tuple(ZZ(value) for value in minimal.a_invariants())
    conductor = ZZ(minimal.conductor())
    root_number = int(minimal.root_number())
    root_numbers[root_number] += 1
    point_label = "+".join(sorted(pair["points"]))
    point_sources[point_label] += 1

    common_points = {}
    if "zero" in pair["points"]:
        common_points["zero"] = {
            "t": "0",
            "u": rational_text(QQ(q_left[0]).sqrt()),
            "v": rational_text(QQ(q_right[0]).sqrt()),
        }
    if "infinity" in pair["points"]:
        common_points["infinity"] = {
            "u_over_t": rational_text(QQ(q_left[2]).sqrt()),
            "v_over_t": rational_text(QQ(q_right[2]).sqrt()),
        }

    quartic_coefficients = tuple(ZZ(quartic[index]) for index in range(5))
    minimal_coefficient_bits = max(integer_bits(value) for value in minimal_model)
    quartic_coefficient_bits = max(integer_bits(value) for value in quartic_coefficients)
    priority_ranks = tuple(
        sorted((int(left_record["priority_rank"]), int(right_record["priority_rank"])))
    )
    rows.append(
        {
            "pair_key": f"{key[0]}:{key[1]}",
            "orbit_masks": list(key),
            "orbit_hex": [f"0x{mask:05x}" for mask in key],
            "priority_ranks": list(priority_ranks),
            "equation_ranks": sorted(
                (int(left_record["equation_rank"]), int(right_record["equation_rank"]))
            ),
            "rational_point_sources": sorted(pair["points"]),
            "common_points": common_points,
            "covers": [
                {
                    "lattice_orbit_mask": int(entry["record"]["lattice_orbit_mask"]),
                    "orbit_hex": f"0x{int(entry['record']['lattice_orbit_mask']):05x}",
                    "priority_rank": int(entry["record"]["priority_rank"]),
                    "equation_rank": int(entry["record"]["equation_rank"]),
                    "q_coefficients_low_to_high": [int(value) for value in entry["q"]],
                    "batch_q_square_rescaling": rational_text(entry["square_rescaling"]),
                }
                for entry in (left, right)
            ],
            "q_coefficients_low_to_high": [
                [int(value) for value in left["q"]],
                [int(value) for value in right["q"]],
            ],
            "batch_q_square_rescalings": [
                rational_text(left["square_rescaling"]),
                rational_text(right["square_rescaling"]),
            ],
            "quartic_coefficients_low_to_high": [int(value) for value in quartic_coefficients],
            "binary_quartic_I": rational_text(invariant_i),
            "binary_quartic_J": rational_text(invariant_j),
            "minimal_jacobian_a1_a2_a3_a4_a6": [int(value) for value in minimal_model],
            "conductor": int(conductor),
            "global_root_number": root_number,
            "complexity": {
                "minimal_model_max_coefficient_bits": minimal_coefficient_bits,
                "quartic_max_coefficient_bits": quartic_coefficient_bits,
                "priority_rank_max": max(priority_ranks),
                "priority_rank_sum": sum(priority_ranks),
            },
        }
    )

rows.sort(
    key=lambda row: (
        row["complexity"]["minimal_model_max_coefficient_bits"],
        row["complexity"]["quartic_max_coefficient_bits"],
        row["complexity"]["priority_rank_max"],
        row["complexity"]["priority_rank_sum"],
        row["orbit_masks"],
    )
)
for index, row in enumerate(rows, start=1):
    row["arithmetic_complexity_rank"] = index

published_key = "7713:42110"
new_key = "28257:81769"
by_key = {row["pair_key"]: row for row in rows}
if published_key not in by_key or new_key not in by_key:
    raise ArithmeticError("published or new positive-control pair is absent")

result = {
    "schema": OUTPUT_SCHEMA,
    "status": "PASS_COMPLETE_5566_IMMEDIATE_POINT_PAIR_CATALOGUE",
    "input": {
        "path": display_path(args.input),
        "sha256": digest(args.input),
    },
    "counts": {
        "square_constant_covers": len(at_zero),
        "square_leading_covers": len(at_infinity),
        "immediate_point_pairs": len(rows),
        "rational_point_sources": dict(sorted(point_sources.items())),
        "global_root_numbers": {str(key): value for key, value in sorted(root_numbers.items())},
    },
    "ranking": {
        "key": [
            "minimal_model_max_coefficient_bits",
            "quartic_max_coefficient_bits",
            "priority_rank_max",
            "priority_rank_sum",
            "orbit_masks",
        ],
        "interpretation": (
            "This is an exact equation-size ordering, not a theorem that point-search time or "
            "Mordell--Weil rank is monotone in the score."
        ),
    },
    "positive_control_rows": {
        "published_rank19_pair": by_key[published_key],
        "new_rank_at_least_3_pair": by_key[new_key],
    },
    "pairs": rows,
    "proof_boundary": (
        "Every row is an exact minimal Jacobian/root-number computation for a genus-one paired "
        "base with a displayed rational point. No Mordell--Weil rank follows from the root number; "
        "independent rational points require a separate exact certificate."
    ),
}

args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(
    "ELKIES2026IMMEDIATEPAIRS|"
    f"pairs={len(rows)}|root_plus={root_numbers[1]}|root_minus={root_numbers[-1]}|"
    f"published_complexity_rank={by_key[published_key]['arithmetic_complexity_rank']}|"
    f"new_complexity_rank={by_key[new_key]['arithmetic_complexity_rank']}|"
    f"status={result['status']}|output={display_path(args.output)}"
)
