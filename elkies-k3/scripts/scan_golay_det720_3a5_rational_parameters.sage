#!/usr/bin/env sage-python
"""Scan small rational values of the G720 marked-family free coordinate.

For every reduced rational ``s6=a/b`` in the declared box whose 7-adic
residue lies in the pinned marked disk, invoke the existing fixed-parameter
Hensel certifier.  Rationally reconstruct all 46 coordinates and retain a row
only after all 55 equations vanish over QQ.  Retained points are also screened
for rational 3-torsion and rational halves of P, Q, or P+Q.

This is a bounded rational-parameter search, not a rational parametrization or
an exclusion of rational points outside the box.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
LIFT_SCRIPT = ROOT / "elkies-k3/scripts/certify_golay_det720_3a5_marked_gf7_lift.sage"
FIBRES = GEN / "elkies-k3-golay-det720-3a5-source-ansatz-mod7-v1.json"
MARKING = GEN / "elkies-k3-golay-det720-3a5-pole0-pairs-mod7-nonsquare-v1.json"
DEFAULT_OUTPUT = GEN / "elkies-k3-golay-det720-3a5-rational-parameter-scan-v1.json"

NAMES = (
    [f"a{index}" for index in range(9)]
    + [f"b{index}" for index in range(13)]
    + [f"p{index}" for index in range(5)]
    + [f"q{index}" for index in range(7)]
    + [f"r{index}" for index in range(5)]
    + [f"s{index}" for index in range(7)]
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_equations():
    coefficient_ring = PolynomialRing(QQ, names=NAMES, order="degrevlex")
    variables = list(coefficient_ring.gens())
    cursor = 0

    def take(length):
        nonlocal cursor
        block = variables[cursor : cursor + length]
        cursor += length
        return block

    a, b = take(9), take(13)
    p, q = take(5), take(7)
    r, s = take(5), take(7)
    function_ring = PolynomialRing(coefficient_ring, "t")
    t = function_ring.gen()

    def polynomial(coefficients):
        return sum(value * t**index for index, value in enumerate(coefficients))

    A, B = polynomial(a), polynomial(b)
    X_P, Y_P = polynomial(p), polynomial(q)
    X_Q, Y_Q = polynomial(r), polynomial(s)
    D = 4 * A**3 + 27 * B**2
    node_P = 2 * A * X_P + 3 * B
    equations = [a[0] + 27]
    equations.extend(D[index] for index in range(6))
    equations.extend(D(t + 1)[index] for index in range(6))
    equations.extend(D[index] for index in range(19, 25))
    equations.extend(
        [
            node_P[0],
            q[0],
            node_P(t + 1)[0],
            Y_P(1),
            node_P[12],
            node_P[11],
            node_P[10],
            q[6],
            q[5],
            q[4],
        ]
    )
    equations.extend((Y_P**2 - X_P**3 - A * X_P - B)[index] for index in range(13))
    equations.extend((Y_Q**2 - X_Q**3 - A * X_Q - B)[index] for index in range(13))
    assert cursor == 46 and len(equations) == 55
    return equations


def candidates(numerator_bound, denominator_bound, required_residue):
    answer = []
    for denominator in range(1, denominator_bound + 1):
        if denominator % 7 == 0:
            continue
        for numerator in range(-numerator_bound, numerator_bound + 1):
            if math.gcd(numerator, denominator) != 1:
                continue
            if numerator * pow(denominator, -1, 7) % 7 != required_residue:
                continue
            answer.append((numerator, denominator))
    return answer


def rational_factor_degrees(polynomial):
    return sorted(
        [int(factor.degree()) for factor, multiplicity in polynomial.factor() for unused in range(multiplicity)]
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--numerator-bound", type=int, default=40)
parser.add_argument("--denominator-bound", type=int, default=40)
parser.add_argument("--lift-precision", type=int, default=32)
parser.add_argument("--workers", type=int, default=8)
parser.add_argument("--marked-model-rank", type=int, default=0)
parser.add_argument("--marked-pair-index", type=int, default=0)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
if arguments.numerator_bound < 1 or arguments.denominator_bound < 1:
    parser.error("rational scan bounds must be positive")
if arguments.lift_precision < 2:
    parser.error("lift precision must be at least two")

output_path = arguments.output.resolve()
marking_payload = json.loads(MARKING.read_text())
marked_models = [row for row in marking_payload["models"] if row["marked_mw2_pairs"]]
if not 0 <= arguments.marked_model_rank < len(marked_models):
    parser.error("--marked-model-rank is outside the marked-model inventory")
marked_model = marked_models[arguments.marked_model_rank]
if not 0 <= arguments.marked_pair_index < len(marked_model["marked_mw2_pairs"]):
    parser.error("--marked-pair-index is outside the selected model")
marked_pair = marked_model["marked_mw2_pairs"][arguments.marked_pair_index]
right_section = marked_model["basis_section_candidates"][1][
    int(marked_pair["right_section_index"])
]
right_y = right_section["Y_coefficients_low_to_high"]
required_residue = int(right_y[6] if len(right_y) > 6 else 0) % 7
parameter_pairs = candidates(
    arguments.numerator_bound, arguments.denominator_bound, required_residue
)
modulus = ZZ(7) ** arguments.lift_precision
sage = shutil.which("sage")
if sage is None:
    raise SystemExit("sage executable not found")
equations = build_equations()


def run_one(task):
    index, (numerator, denominator), temporary_directory = task
    residue = int(ZZ(numerator) * ZZ(denominator).inverse_mod(modulus) % modulus)
    lift_path = Path(temporary_directory) / f"lift-{index:05d}.json"
    command = [
        sage,
        "-python",
        str(LIFT_SCRIPT),
        "--free-parameter-integer",
        str(residue),
        "--lift-precision",
        str(arguments.lift_precision),
        "--output",
        str(lift_path),
        "--marked-model-rank",
        str(arguments.marked_model_rank),
        "--marked-pair-index",
        str(arguments.marked_pair_index),
    ]
    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode:
        return {
            "parameter": f"{numerator}/{denominator}",
            "status": "LIFT_ERROR",
            "stderr_tail": completed.stderr[-500:],
        }
    lift = json.loads(lift_path.read_text())
    finite = lift["finite_precision_lift"]
    try:
        coordinates = [
            ZZ(value).rational_reconstruction(ZZ(finite["modulus"]))
            for value in finite["coordinates_modulus"]
        ]
    except (ArithmeticError, ValueError):
        return {"parameter": f"{numerator}/{denominator}", "status": "NO_FULL_RR"}
    if coordinates[-1] != QQ(numerator) / denominator:
        return {"parameter": f"{numerator}/{denominator}", "status": "NO_FULL_RR"}
    if any(equation(*coordinates) for equation in equations):
        return {"parameter": f"{numerator}/{denominator}", "status": "RR_NOT_EXACT"}

    ring = PolynomialRing(QQ, "t")
    function_field = ring.fraction_field()
    A = ring(coordinates[0:9])
    B = ring(coordinates[9:22])
    X_P, Y_P = ring(coordinates[22:27]), ring(coordinates[27:34])
    X_Q, Y_Q = ring(coordinates[34:39]), ring(coordinates[39:46])
    curve = EllipticCurve(function_field, [0, 0, 0, function_field(A), function_field(B)])
    P = curve(function_field(X_P), function_field(Y_P))
    Q = curve(function_field(X_Q), function_field(Y_Q))
    x_ring = PolynomialRing(function_field, "x")
    x = x_ring.gen()
    division_3 = x_ring(3 * x**4 + 6 * A * x**2 + 12 * B * x - A**2)
    halves = {}
    for name, point in (("P", P), ("Q", Q), ("P+Q", P + Q)):
        x_point = function_field(point[0])
        duplication = (
            x**4
            - 4 * x_point * x**3
            - 2 * A * x**2
            - (4 * A * x_point + 8 * B) * x
            + A**2
            - 4 * B * x_point
        )
        degrees = rational_factor_degrees(duplication)
        halves[name] = {"factor_degrees_over_Q(t)": degrees, "linear_factor": 1 in degrees}
    return {
        "parameter": f"{numerator}/{denominator}",
        "status": "EXACT_QQ_POINT",
        "coordinates": [str(value) for value in coordinates],
        "maximum_numerator_or_denominator": int(
            max(max(abs(value.numerator()), value.denominator()) for value in coordinates)
        ),
        "three_division_factor_degrees_over_Q(t)": rational_factor_degrees(division_3),
        "rational_three_torsion_detected": 1 in rational_factor_degrees(division_3),
        "rational_halves": halves,
    }


with tempfile.TemporaryDirectory(prefix="golay720-rational-scan-") as temporary_directory:
    tasks = [
        (index, pair, temporary_directory) for index, pair in enumerate(parameter_pairs)
    ]
    with ThreadPoolExecutor(max_workers=arguments.workers) as executor:
        rows = list(executor.map(run_one, tasks))

status_counts = {}
for row in rows:
    status_counts[row["status"]] = status_counts.get(row["status"], 0) + 1
errors = [row for row in rows if row["status"] == "LIFT_ERROR"]
if errors:
    raise RuntimeError(f"fixed-parameter lifts failed: {errors[:3]}")
exact_points = [row for row in rows if row["status"] == "EXACT_QQ_POINT"]

payload = {
    "schema": "elkies-k3.golay-det720-3a5-rational-parameter-scan.v1",
    "status": "PASS_BOUNDED_EXACT_RATIONAL_PARAMETER_SCAN",
    "reproduce": (
        "sage -python elkies-k3/scripts/scan_golay_det720_3a5_rational_parameters.sage "
        f"--numerator-bound {arguments.numerator_bound} "
        f"--denominator-bound {arguments.denominator_bound} "
        f"--lift-precision {arguments.lift_precision} --workers {arguments.workers} "
        f"--marked-model-rank {arguments.marked_model_rank} "
        f"--marked-pair-index {arguments.marked_pair_index}"
    ),
    "inputs": {
        relative(LIFT_SCRIPT): digest(LIFT_SCRIPT),
        relative(FIBRES): digest(FIBRES),
        relative(MARKING): digest(MARKING),
    },
    "search": {
        "parameter": "s6",
        "numerator_interval": [-arguments.numerator_bound, arguments.numerator_bound],
        "denominator_interval": [1, arguments.denominator_bound],
        "reduced_fractions_only": True,
        "denominator_prime_to_7": True,
        "marked_model_rank": arguments.marked_model_rank,
        "marked_pair_index": arguments.marked_pair_index,
        "fibre_example_index": int(marked_model["example_index"]),
        "required_residue_mod_7": required_residue,
        "candidate_count": len(parameter_pairs),
        "lift_precision_exponent": arguments.lift_precision,
        "status_counts": status_counts,
    },
    "exact_rational_points": exact_points,
    "proof_boundary": (
        "Every retained row is an exact QQ solution of all 55 marked equations.  "
        "Failure of simultaneous rational reconstruction is not a proof that the "
        "corresponding 7-adic point is irrational, and the bounded s6 box is not a "
        "rational-point or family classification.  Torsion/divisibility screens in "
        "this artifact detect factors over QQ(t), not all constant-field extensions."
    ),
}

encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not output_path.exists() or output_path.read_text() != encoded:
        raise SystemExit(f"stale artifact: {output_path}")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(encoded)

print(
    f"GOLAY720RATSCAN|candidates={len(parameter_pairs)}|"
    f"exact_QQ={len(exact_points)}|status_counts={status_counts}|status=PASS"
)
