#!/usr/bin/env sage
"""Search the frozen rank-28 quartics for a second simultaneous split.

status: ACTIVE_SEARCH
claim: exhaustive compact-height and bounded canonical-subgroup negative search
inputs: frozen R17 rank-28 genus-one bisection pilot
outputs: generated simultaneous-splitting search artifact

The compact scan is exhaustive for primitive ``t=a/b`` with
``|a|, b <= H``.  A C++ hot loop intersects exact quadratic-residue masks at
good primes; every surviving parameter is then tested by an exact integer
square criterion.  A complementary search enumerates multiples of the
canonical opposite-ordinate point on each pointed quartic Jacobian.

This is a bounded experiment.  An empty result is not a global nonexistence
theorem and says nothing about rational points outside the two stated search
regions.
"""

from __future__ import annotations

import argparse
from array import array
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import shutil
import struct
import subprocess

from sage.all import EllipticCurve, PolynomialRing, QQ


SCRIPTS = Path(__file__).resolve().parent
ROOT = SCRIPTS.parents[1]
INPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-rank28-genus-one-bisection-pilot-v1.json"
)
CPP = SCRIPTS / "scan_elkies_2026_rank28_genus_one_splitting.cpp"
OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-rank28-simultaneous-splitting-h10000-v1.json"
)
LOCAL = ROOT / "artifacts/local/elkies-k3/r17-rank28-genus-one-splitting"
PRIMES = (19, 23, 31, 37, 43, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 107)
MAGIC = 0x47315331


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def rational_text(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else str(value)


def normalized_integer_quartic(record: dict) -> tuple[int, tuple[int, ...]]:
    coefficients = [
        Fraction(value)
        for value in record["branch_polynomial_q_coefficients_low_to_high"]
    ]
    denominator = 1
    for value in coefficients:
        denominator = math.lcm(denominator, value.denominator)
    integral = tuple(
        value.numerator * (denominator // value.denominator)
        for value in coefficients
    )
    if len(integral) != 5:
        raise ArithmeticError("expected a quartic with five coefficients")
    return denominator, integral


def homogeneous_value(curve: tuple[int, tuple[int, ...]], a: int, b: int) -> int:
    _denominator, coefficients = curve
    b2 = b * b
    b3 = b2 * b
    b4 = b2 * b2
    return (
        (((coefficients[4] * a + coefficients[3] * b) * a + coefficients[2] * b2) * a
         + coefficients[1] * b3)
        * a
        + coefficients[0] * b4
    )


def exact_square_root(
    curve: tuple[int, tuple[int, ...]], a: int, b: int
) -> int | None:
    denominator, _coefficients = curve
    radicand = denominator * homogeneous_value(curve, a, b)
    if radicand < 0:
        return None
    root = math.isqrt(radicand)
    return root if root * root == radicand else None


def export_residue_tables(
    path: Path, curves: list[tuple[int, tuple[int, ...]]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as stream:
        stream.write(struct.pack("<III", MAGIC, len(curves), len(PRIMES)))
        for prime in PRIMES:
            if any(denominator % prime == 0 for denominator, _ in curves):
                raise ArithmeticError(f"bad residue prime {prime}")
            squares = {value * value % prime for value in range(prime)}
            masks = array("H")
            for a in range(prime):
                for b in range(prime):
                    mask = 0
                    for index, (denominator, coefficients) in enumerate(curves):
                        value = denominator * homogeneous_value(
                            (1, tuple(coefficient % prime for coefficient in coefficients)),
                            a,
                            b,
                        )
                        if value % prime in squares:
                            mask |= 1 << index
                    masks.append(mask)
            stream.write(struct.pack("<I", prime))
            stream.write(masks.tobytes())


def compile_and_scan(
    table_path: Path, candidate_path: Path, binary_path: Path, height: int
) -> tuple[list[str], str]:
    compiler = shutil.which("g++")
    if compiler is None:
        raise SystemExit("g++ is required for the exhaustive compact scan")
    compile_command = [
        compiler,
        "-O3",
        "-std=c++17",
        "-Wall",
        "-Wextra",
        "-pedantic",
        str(CPP),
        "-o",
        str(binary_path),
    ]
    subprocess.run(compile_command, check=True)
    scan_command = [str(binary_path), str(table_path), str(height), str(candidate_path)]
    completed = subprocess.run(scan_command, check=True, text=True, capture_output=True)
    summary = completed.stdout.strip()
    if "|status=PASS" not in summary:
        raise ArithmeticError("compact scanner did not return PASS")
    return compile_command, summary


def exact_compact_results(
    candidate_path: Path,
    curves: list[tuple[int, tuple[int, ...]]],
    t0: Fraction,
) -> dict:
    candidates = []
    exact_tests = 0
    for line in candidate_path.read_text().splitlines():
        numerator, denominator, mask = map(int, line.split())
        split = []
        roots = {}
        for index, curve in enumerate(curves):
            if not ((mask >> index) & 1):
                continue
            exact_tests += 1
            root = exact_square_root(curve, numerator, denominator)
            if root is not None:
                split.append(index + 1)
                roots[str(index + 1)] = str(root)
        if len(split) >= 2:
            candidates.append(
                {
                    "t": rational_text(QQ(numerator) / denominator),
                    "projective_pair": [numerator, denominator],
                    "split_target_indices_one_based": split,
                    "integral_square_roots": roots,
                    "is_original_target_parameter": Fraction(numerator, denominator) == t0,
                }
            )
    original = [row for row in candidates if row["is_original_target_parameter"]]
    if len(original) != 1 or len(original[0]["split_target_indices_one_based"]) != len(curves):
        raise ArithmeticError("the eleven-fold t0 positive control was not recovered")
    return {
        "modular_candidate_count": sum(1 for _ in candidate_path.open()),
        "exact_square_tests": exact_tests,
        "simultaneous_splits": candidates,
        "new_simultaneous_splits_away_from_t0": [
            row for row in candidates if not row["is_original_target_parameter"]
        ],
    }


def residue_survives(
    curve: tuple[int, tuple[int, ...]], numerator: int, denominator: int
) -> bool:
    scalar, coefficients = curve
    for prime in PRIMES:
        value = scalar * homogeneous_value(
            (1, tuple(coefficient % prime for coefficient in coefficients)),
            numerator % prime,
            denominator % prime,
        )
        if value % prime != 0 and pow(value % prime, (prime - 1) // 2, prime) != 1:
            return False
    return True


def jacobian_multiple_search(document: dict, curves, multiple_bound: int) -> dict:
    records = document["traces"][0]["targets"]
    t0 = QQ(document["parameter"])
    t_ring = PolynomialRing(QQ, "t")
    t = t_ring.gen()
    z_ring = PolynomialRing(QQ, "z")
    z = z_ring.gen()
    quartics = [
        t_ring([QQ(value) for value in row["branch_polynomial_q_coefficients_low_to_high"]])
        for row in records
    ]
    summaries = []
    hits = []
    exact_cross_tests = 0
    for source_index, record in enumerate(records):
        shifted = z_ring(quartics[source_index](z + t0))
        e, d, c, b, a = (QQ(shifted[index]) for index in range(5))
        v0 = QQ(record["cover_rational_witness"]["s"])
        if v0 * v0 != e:
            raise ArithmeticError("pointed quartic constant is not the recorded square")
        curve = EllipticCurve(
            QQ,
            [
                d / v0,
                c - d**2 / (4 * v0**2),
                2 * v0 * b,
                -4 * v0**2 * a,
                a * (d**2 - 4 * v0**2 * c),
            ],
        )
        opposite_x = d**2 / (4 * v0**2) - c
        point = curve.lift_x(opposite_x, all=True)[0]
        if point.has_finite_order():
            raise ArithmeticError("canonical opposite-ordinate point became torsion")
        current = point
        maximum_bits = 0
        modular_cross_survivors = 0
        for multiple in range(2, multiple_bound + 1):
            current += point
            if current.is_zero() or current[1] == 0:
                continue
            x_value, y_value = current[:2]
            s_value = (
                4 * v0**2 * (x_value + c) - d**2
            ) / (2 * v0 * y_value)
            if s_value == 0:
                continue
            quartic_y = (x_value * s_value**2 - d * s_value) / (2 * v0) - v0
            if multiple == 2 and quartic_y**2 != shifted(s_value):
                raise ArithmeticError("pointed-quartic inverse regression failed")
            parameter = s_value + t0
            numerator = int(parameter.numerator())
            denominator = int(parameter.denominator())
            maximum_bits = max(
                maximum_bits,
                abs(numerator).bit_length(),
                denominator.bit_length(),
            )
            split = [source_index + 1]
            for target_index, target_curve in enumerate(curves):
                if target_index == source_index:
                    continue
                if not residue_survives(target_curve, numerator, denominator):
                    continue
                modular_cross_survivors += 1
                exact_cross_tests += 1
                if exact_square_root(target_curve, numerator, denominator) is not None:
                    split.append(target_index + 1)
            if len(split) >= 2:
                hits.append(
                    {
                        "source_target_index_one_based": source_index + 1,
                        "multiple": multiple,
                        "t": rational_text(parameter),
                        "split_target_indices_one_based": sorted(split),
                    }
                )
        summaries.append(
            {
                "source_target_index_one_based": source_index + 1,
                "multiples_tested": [2, multiple_bound],
                "maximum_projective_coordinate_bits": maximum_bits,
                "modular_cross_survivors": modular_cross_survivors,
            }
        )
    return {
        "canonical_point": "opposite ordinate above t0 on the pointed quartic",
        "canonical_point_is_nontorsion_for_every_cover": True,
        "per_cover": summaries,
        "exact_cross_tests": exact_cross_tests,
        "simultaneous_split_hits": hits,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--local-directory", type=Path, default=LOCAL)
    parser.add_argument("--height", type=int, default=10_000)
    parser.add_argument("--multiple-bound", type=int, default=30)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.height < 1 or args.multiple_bound < 2:
        parser.error("search bounds must be positive and --multiple-bound at least 2")

    source = json.loads(args.input.read_text())
    if source.get("status") != "PASS_EXACT_R17_RANK28_GENUS_ONE_BISECTION_PILOT":
        raise SystemExit("the frozen genus-one bisection pilot is missing")
    records = source["traces"][0]["targets"]
    curves = [normalized_integer_quartic(record) for record in records]
    t0 = Fraction(source["parameter"])

    args.local_directory.mkdir(parents=True, exist_ok=True)
    table_path = args.local_directory / "square-residue-masks.bin"
    candidate_path = args.local_directory / "compact-candidates.txt"
    binary_path = args.local_directory / "scan-genus-one-splitting"
    export_residue_tables(table_path, curves)
    compile_command, scanner_summary = compile_and_scan(
        table_path, candidate_path, binary_path, args.height
    )
    compact = exact_compact_results(candidate_path, curves, t0)
    multiples = jacobian_multiple_search(source, curves, args.multiple_bound)
    if compact["new_simultaneous_splits_away_from_t0"] or multiples["simultaneous_split_hits"]:
        status = "PASS_NEW_SIMULTANEOUS_SPLIT_CANDIDATES_REQUIRE_INDEPENDENCE"
        independence = "not run: promote the reported candidates in a separate exact specialization certificate"
    else:
        status = "PASS_EXHAUSTIVE_BOUNDED_NO_NEW_SIMULTANEOUS_SPLIT"
        independence = "vacuous: no new exact simultaneous split survived"

    result = {
        "schema": "elkies-k3.r17-rank28-genus-one-simultaneous-splitting.v1",
        "status": status,
        "input": {
            "path": str(args.input.relative_to(ROOT)),
            "sha256": file_sha256(args.input),
            "frozen_quartic_count": len(curves),
            "original_target_parameter": source["parameter"],
        },
        "compact_projective_scan": {
            "region": f"primitive t=a/b with |a| <= {args.height}, 1 <= b <= {args.height}",
            "height": args.height,
            "residue_primes": list(PRIMES),
            "scanner_summary": scanner_summary,
            **compact,
        },
        "pointed_jacobian_multiple_scan": multiples,
        "specialization_and_quotient_independence": independence,
        "proof_boundary": (
            "The compact t-line result is exhaustive only in the displayed primitive box. "
            "The Jacobian result is exhaustive only in the cyclic subgroups and multiple range "
            "shown. Their union is a bounded negative experiment, not a theorem that the fibre "
            "products have no further rational points."
        ),
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "search_elkies_2026_rank28_simultaneous_splitting.sage"
        ),
        "implementation": {
            "script": str(Path(__file__).relative_to(ROOT)),
            "script_sha256": file_sha256(Path(__file__)),
            "scanner_source": str(CPP.relative_to(ROOT)),
            "scanner_source_sha256": file_sha256(CPP),
            "compile_command": compile_command,
            "residue_table_sha256": file_sha256(table_path),
        },
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        expected = json.loads(args.output.read_text())
        # Paths and hashes are deterministic; compare the full artifact.
        if expected != result:
            raise ArithmeticError("stored simultaneous-splitting artifact changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(
        f"R17RANK28G1SPLIT|height={args.height}|"
        f"modular={compact['modular_candidate_count']}|"
        f"new={len(compact['new_simultaneous_splits_away_from_t0'])}|"
        f"multiples={args.multiple_bound}|status={status}|output={args.output}"
    )


if __name__ == "__main__":
    main()
