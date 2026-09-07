#!/usr/bin/env sage-python
"""Scan small rational values of the NS0031 formal coordinate m9.

For every reduced rational value in the declared box whose 7-adic residue
lies in the certified model-157 disk, fix m9 throughout the exact Hensel lift,
rationally reconstruct all 52 coordinates, and substitute them into all 59
normalized equations over QQ.  Retained points must also pass the exact
``I2+2I8+6I1`` open gates.

This is a bounded rational-parameter search.  A miss is not an obstruction to
rational points on the marked curve, and a hit would still require a separate
Picard-rank and primitive Neron--Severi audit.
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

from sage.all import PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
LIFT_SCRIPT = (
    ROOT
    / "elkies-k3/scripts/"
    "certify_lattice_foundry_ns0031_marked_gf7_hensel.sage"
)
HENSEL = GEN / "elkies-k3-lattice-foundry-ns0031-marked-gf7-hensel-v1.json"
DEFAULT_OUTPUT = (
    GEN / "elkies-k3-lattice-foundry-ns0031-rational-parameter-scan-v1.json"
)

NAMES = (
    [f"a{index}" for index in range(9)]
    + [f"b{index}" for index in range(13)]
    + [f"p{index}" for index in range(5)]
    + [f"q{index}" for index in range(7)]
    + ["c0"]
    + [f"n{index}" for index in range(7)]
    + [f"m{index}" for index in range(10)]
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_system():
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
    c0 = take(1)[0]
    n, m = take(7), take(10)
    polynomial_ring = PolynomialRing(coefficient_ring, "t")
    t = polynomial_ring.gen()

    def polynomial(coefficients):
        return sum(value * t**index for index, value in enumerate(coefficients))

    A, B = polynomial(a), polynomial(b)
    X, Y = polynomial(p), polynomial(q)
    C, N, M = t + c0, polynomial(n), polynomial(m)
    D = 4 * A**3 + 27 * B**2
    equations = [a[0] + 3]
    equations.extend(D[index] for index in range(2))
    equations.extend(D(t + 1)[index] for index in range(8))
    equations.extend(D[index] for index in range(17, 25))
    equations.extend(
        [
            2 * a[0] * p[0] + 3 * b[0],
            q[0],
            2 * a[8] * p[4] + 3 * b[12],
            2 * (a[8] * p[3] + a[7] * p[4]) + 3 * b[11],
            q[6],
            q[5],
            2 * A(1) * N(1) + 3 * B(1) * C(1) ** 2,
            M(1),
        ]
    )
    equations.extend((Y**2 - X**3 - A * X - B)[index] for index in range(13))
    equations.extend(
        (M**2 - N**3 - A * N * C**4 - B * C**6)[index]
        for index in range(19)
    )
    if cursor != 52 or len(equations) != 59:
        raise ArithmeticError("unexpected NS0031 marked-system size")
    return variables, equations


def candidates(numerator_bound, denominator_bound, required_residue):
    answer = []
    for denominator in range(1, denominator_bound + 1):
        if denominator % 7 == 0:
            continue
        for numerator in range(-numerator_bound, numerator_bound + 1):
            if math.gcd(numerator, denominator) != 1:
                continue
            if numerator * pow(denominator, -1, 7) % 7 == required_residue:
                answer.append((numerator, denominator))
    return answer


def exact_surface_gate(coordinates):
    ring = PolynomialRing(QQ, "t")
    t = ring.gen()
    A = ring(coordinates[0:9])
    B = ring(coordinates[9:22])
    c0 = coordinates[34]
    D = 4 * A**3 + 27 * B**2
    shifted = ring(D(t + 1))
    if not (D[2] and shifted[8] and D[16]):
        return None
    if not (A(0) and A(1) and A[8] and 1 + c0):
        return None
    quotient = D // (t**2 * (t - 1) ** 8)
    if D != t**2 * (t - 1) ** 8 * quotient or quotient.degree() != 6:
        return None
    if not quotient.is_squarefree() or quotient.gcd(t * (t - 1)).degree() != 0:
        return None
    return {
        "discriminant_orders": {"zero": 2, "one": 8, "infinity": 8},
        "remaining_discriminant_degree": 6,
        "remaining_factor_degrees_over_Q": sorted(
            int(factor.degree())
            for factor, multiplicity in quotient.factor()
            for unused in range(multiplicity)
        ),
        "remaining_discriminant_squarefree": True,
        "A_units_at_marked_supports": True,
        "pole_denominator_unit_at_one": True,
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--numerator-bound", type=int, default=40)
parser.add_argument("--denominator-bound", type=int, default=40)
parser.add_argument("--lift-precision", type=int, default=40)
parser.add_argument("--workers", type=int, default=8)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
if arguments.numerator_bound < 1 or arguments.denominator_bound < 1:
    parser.error("rational scan bounds must be positive")
if arguments.lift_precision < 2 or arguments.workers < 1:
    parser.error("lift precision must be at least two and workers positive")

sage = shutil.which("sage")
if sage is None:
    raise SystemExit("sage executable not found")
hensel = json.loads(HENSEL.read_text())
if hensel["status"] != (
    "PASS_ONE_DIMENSIONAL_MARKED_TANGENT_AND_EXPLICIT_Z7_"
    "LIFT_TO_REQUESTED_PRECISION"
):
    raise ValueError("missing NS0031 marked Hensel seed")
if hensel["jacobian_certificate"]["omitted_free_variable_names"] != ["m9"]:
    raise ValueError("NS0031 formal coordinate changed")
required_residue = int(hensel["seed"]["coordinates_mod_7"][-1]) % 7
parameter_pairs = candidates(
    arguments.numerator_bound, arguments.denominator_bound, required_residue
)
variables, equations = build_system()
modulus = ZZ(7) ** arguments.lift_precision


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
    ]
    completed = subprocess.run(command, capture_output=True, text=True)
    parameter = f"{numerator}/{denominator}"
    if completed.returncode:
        return {
            "parameter_m9": parameter,
            "status": "LIFT_ERROR",
            "stderr_tail": completed.stderr[-500:],
        }
    lift = json.loads(lift_path.read_text())
    finite = lift["finite_precision_lift"]
    if finite["achieved_precision_exponent"] != arguments.lift_precision:
        return {"parameter_m9": parameter, "status": "FINITE_LIFT_STOPPED"}
    try:
        coordinates = [
            QQ(ZZ(value).rational_reconstruction(ZZ(finite["modulus"])))
            for value in finite["coordinates_modulus"]
        ]
    except (ArithmeticError, ValueError):
        return {"parameter_m9": parameter, "status": "NO_FULL_RR"}
    if coordinates[-1] != QQ(numerator) / denominator:
        return {"parameter_m9": parameter, "status": "NO_FULL_RR"}
    if any(equation(*coordinates) for equation in equations):
        return {"parameter_m9": parameter, "status": "RR_NOT_EXACT"}
    surface_gate = exact_surface_gate(coordinates)
    if surface_gate is None:
        return {"parameter_m9": parameter, "status": "EXACT_QQ_POINT_OPEN_GATE_FAILED"}
    return {
        "parameter_m9": parameter,
        "status": "EXACT_QQ_POINT_REQUIRES_NS_AUDIT",
        "coordinates": [str(value) for value in coordinates],
        "maximum_numerator_or_denominator": int(
            max(
                max(abs(value.numerator()), value.denominator())
                for value in coordinates
            )
        ),
        "surface_gate": surface_gate,
    }


with tempfile.TemporaryDirectory(prefix="ns0031-rational-scan-") as temporary:
    tasks = [
        (index, pair, temporary) for index, pair in enumerate(parameter_pairs)
    ]
    with ThreadPoolExecutor(max_workers=arguments.workers) as executor:
        rows = list(executor.map(run_one, tasks))

status_counts = {}
for row in rows:
    status_counts[row["status"]] = status_counts.get(row["status"], 0) + 1
errors = [row for row in rows if row["status"] == "LIFT_ERROR"]
if errors:
    raise RuntimeError(f"fixed-parameter lifts failed: {errors[:3]}")
exact_points = [
    row for row in rows if row["status"] == "EXACT_QQ_POINT_REQUIRES_NS_AUDIT"
]

payload = {
    "schema": "elkies-k3.lattice-foundry-ns0031-rational-parameter-scan.v1",
    "status": "PASS_BOUNDED_EXACT_RATIONAL_PARAMETER_SCAN",
    "inputs": {
        relative(LIFT_SCRIPT): digest(LIFT_SCRIPT),
        relative(HENSEL): digest(HENSEL),
    },
    "search": {
        "free_parameter": "m9",
        "required_residue_mod_7": required_residue,
        "numerator_bound": arguments.numerator_bound,
        "denominator_bound": arguments.denominator_bound,
        "candidate_count": len(parameter_pairs),
        "lift_precision_exponent": arguments.lift_precision,
        "status_counts": status_counts,
    },
    "exact_rational_points": exact_points,
    "proof_boundary": {
        "proved": (
            "Every reduced rational m9 in the declared box and model-157 "
            "residue disk was tested. Every retained point satisfies all 59 "
            "normalized equations over QQ and the exact I2+2I8+6I1 open gates."
        ),
        "not_proved": (
            "Reconstruction failure does not prove 7-adic irrationality, and "
            "the bounded scan does not classify rational points. Any hit still "
            "requires geometric Picard-rank and primitive NS0031 audits."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/scan_lattice_foundry_ns0031_rational_parameters.sage "
        f"--numerator-bound {arguments.numerator_bound} "
        f"--denominator-bound {arguments.denominator_bound} "
        f"--lift-precision {arguments.lift_precision} "
        f"--workers {arguments.workers}"
    ),
}
serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
output_path = arguments.output.resolve()
if arguments.check:
    if not output_path.exists() or output_path.read_text() != serialized:
        raise SystemExit(f"stale artifact: {output_path}")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)
print(
    f"NS0031RATSCAN|candidates={len(parameter_pairs)}|"
    f"exact_QQ={len(exact_points)}|status_counts={status_counts}|status=PASS"
)
