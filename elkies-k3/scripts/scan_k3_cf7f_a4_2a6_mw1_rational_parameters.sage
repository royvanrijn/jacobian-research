#!/usr/bin/env sage-python
"""Scan small integral m8 parameters in both smooth determinant-714 MW1 disks.

Each candidate fixes the unique nonpivot coordinate of a marked ``GF(7)``
seed, lifts all 47 equations digit by digit, rationally reconstructs all forty
coordinates, and substitutes them in the normalized ``I5+I7+I7`` marked
system over ``QQ``.  A reconstruction is only a candidate for the intended
determinant-714 K3 until torsion, divisibility, Picard-rank, and primitive-
closure gates are checked independently.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
LIFT_SCRIPT = ROOT / "elkies-k3/scripts/certify_k3_04b_a3_a4_a9_marked_gf5_hensel.sage"
FIBRES = GEN / "elkies-k3-k3-cf7f6c91a3a40d32-a4-2a6-mw1-fibre-ansatz-mod7-v1.json"
MARKING = GEN / "elkies-k3-k3-cf7f6c91a3a40d32-a4-2a6-mw1-pole1-marking-mod7-nonsquare-v1.json"
DEFAULT_OUTPUT = GEN / "elkies-k3-k3-cf7f6c91a3a40d32-a4-2a6-mw1-rational-parameter-scan-v1.json"
HENSEL_SCHEMA = "elkies-k3.k3-cf7f-a4-2a6-mw1-marked-gf7-hensel.v1"

NAMES = (
    [f"a{index}" for index in range(9)]
    + [f"b{index}" for index in range(13)]
    + ["c0"]
    + [f"n{index}" for index in range(7)]
    + [f"m{index}" for index in range(10)]
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_equations(twist_integer):
    coefficient_ring = PolynomialRing(QQ, names=NAMES, order="degrevlex")
    variables = list(coefficient_ring.gens())
    a = variables[0:9]
    b = variables[9:22]
    c0 = variables[22]
    n = variables[23:30]
    m = variables[30:40]
    function_ring = PolynomialRing(coefficient_ring, "t")
    t = function_ring.gen()

    def polynomial(coefficients):
        return sum(value * t**index for index, value in enumerate(coefficients))

    A, B = polynomial(a), polynomial(b)
    C, N, M = t + c0, polynomial(n), polynomial(m)
    discriminant = 4 * A**3 + 27 * B**2
    node = 2 * A * N + 3 * B * C**2
    residual = M**2 - N**3 - A * N * C**4 - B * C**6
    equations = [a[0] + 3 * QQ(twist_integer) ** 2]
    equations.extend(discriminant[index] for index in range(5))
    equations.extend(discriminant(t + 1)[index] for index in range(7))
    equations.extend(discriminant[index] for index in range(18, 25))
    equations.extend(node[index] for index in range(1))
    equations.extend(m[index] for index in range(1))
    equations.extend(node(t + 1)[index] for index in range(2))
    equations.extend(M(t + 1)[index] for index in range(2))
    equations.extend(node[index] for index in range(14, 15))
    equations.extend(m[index] for index in range(9, 10))
    equations.extend(residual[index] for index in range(19))
    if len(variables) != 40 or len(equations) != 47:
        raise ArithmeticError("unexpected marked-system size")
    return equations


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--minimum", type=int, default=-40)
parser.add_argument("--maximum", type=int, default=40)
parser.add_argument("--lift-precision", type=int, default=40)
parser.add_argument("--workers", type=int, default=4)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
if arguments.minimum > arguments.maximum:
    parser.error("--minimum must not exceed --maximum")
if arguments.lift_precision < 2 or arguments.workers < 1:
    parser.error("lift precision must be at least two and workers positive")

sage = shutil.which("sage")
if sage is None:
    raise SystemExit("sage executable not found")
marking = json.loads(MARKING.read_text())
if (
    marking["prime"] != 7
    or marking["status"]
    != "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_WITH_MARKED_MW1_SECTION"
    or len(marking["models"]) != 2
):
    raise ValueError("expected two positive determinant-714 marked models")

branches = []
for model_rank, model in enumerate(marking["models"]):
    section = model["pole_one_sections"][0]
    y = list(section["Y_numerator_coefficients_low_to_high"])
    y.extend([0] * (10 - len(y)))
    branches.append(
        {
            "marked_model_rank": model_rank,
            "fibre_example_index": int(model["example_index"]),
            "required_m8_residue_mod_7": int(y[8]) % 7,
        }
    )

tasks = [
    (branch_rank, branch, parameter)
    for branch_rank, branch in enumerate(branches)
    for parameter in range(arguments.minimum, arguments.maximum + 1)
    if parameter % 7 == branch["required_m8_residue_mod_7"]
]
equations = build_equations(int(marking["quadratic_twist"]))


def run_one(task):
    branch_rank, branch, parameter, temporary_directory = task
    lift_path = Path(temporary_directory) / f"lift-b{branch_rank}-m8-{parameter}.json"
    command = [
        sage,
        "-python",
        str(LIFT_SCRIPT),
        "--fibres",
        str(FIBRES),
        "--marking",
        str(MARKING),
        "--schema",
        HENSEL_SCHEMA,
        "--marked-model-rank",
        str(branch["marked_model_rank"]),
        "--free-parameter-integer",
        str(parameter),
        "--lift-precision",
        str(arguments.lift_precision),
        "--output",
        str(lift_path),
    ]
    completed = subprocess.run(command, capture_output=True, text=True)
    base = {
        "branch_rank": branch_rank,
        "fibre_example_index": branch["fibre_example_index"],
        "marked_model_rank": branch["marked_model_rank"],
        "free_parameter_m8": parameter,
    }
    if completed.returncode:
        return {**base, "status": "LIFT_ERROR", "stderr_tail": completed.stderr[-500:]}
    lift = json.loads(lift_path.read_text())
    finite = lift["finite_precision_lift"]
    if finite["achieved_precision_exponent"] != arguments.lift_precision:
        return {**base, "status": "FINITE_LIFT_STOPPED"}
    modulus = ZZ(finite["modulus"])
    try:
        coordinates = [
            QQ(ZZ(value).rational_reconstruction(modulus))
            for value in finite["coordinates_modulus"]
        ]
    except (ArithmeticError, ValueError):
        return {**base, "status": "NO_FULL_RR"}
    if coordinates[38] != parameter:
        return {**base, "status": "NO_FULL_RR"}
    if any(equation(*coordinates) for equation in equations):
        return {**base, "status": "RR_NOT_EXACT"}
    return {
        **base,
        "status": "EXACT_QQ_POINT_REQUIRES_SATURATION_AUDIT",
        "coordinates": [str(value) for value in coordinates],
        "maximum_numerator_or_denominator": int(
            max(max(abs(value.numerator()), value.denominator()) for value in coordinates)
        ),
    }


with tempfile.TemporaryDirectory(prefix="k3cf7f-rational-scan-", dir=GEN) as temporary:
    prepared = [task + (temporary,) for task in tasks]
    with ThreadPoolExecutor(max_workers=arguments.workers) as executor:
        rows = list(executor.map(run_one, prepared))

status_counts = {}
for row in rows:
    status_counts[row["status"]] = status_counts.get(row["status"], 0) + 1
errors = [row for row in rows if row["status"] == "LIFT_ERROR"]
if errors:
    raise RuntimeError(f"fixed-parameter lifts failed: {errors[:3]}")
exact_rows = [row for row in rows if row["status"].startswith("EXACT_QQ_POINT")]

payload = {
    "schema": "elkies-k3.k3-cf7f-a4-2a6-mw1-rational-parameter-scan.v1",
    "status": "PASS_BOUNDED_EXACT_INTEGRAL_PARAMETER_SCAN",
    "inputs": {
        relative(LIFT_SCRIPT): digest(LIFT_SCRIPT),
        relative(FIBRES): digest(FIBRES),
        relative(MARKING): digest(MARKING),
    },
    "search": {
        "free_parameter": "m8",
        "integer_interval": [arguments.minimum, arguments.maximum],
        "branch_count": len(branches),
        "branches": branches,
        "candidate_count": len(tasks),
        "lift_precision_exponent": arguments.lift_precision,
        "status_counts": status_counts,
    },
    "exact_rational_points": exact_rows,
    "proof_boundary": {
        "proved": (
            "Every integer m8 in the declared interval and required residue class was "
            "tested in both smooth marked GF(7) disks. Every retained reconstruction "
            "satisfies all 47 normalized equations over Q."
        ),
        "not_proved": (
            "Reconstruction failure does not prove 7-adic irrationality, and the bounded "
            "integral scan does not classify rational points. Any exact point still "
            "requires torsion, divisibility, Picard-rank, and primitive-closure audits."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/scan_k3_cf7f_a4_2a6_mw1_rational_parameters.sage "
        f"--minimum {arguments.minimum} --maximum {arguments.maximum} "
        f"--lift-precision {arguments.lift_precision} --workers {arguments.workers}"
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
    f"K3CF7FRATSCAN|branches={len(branches)}|candidates={len(tasks)}|"
    f"exact_QQ={len(exact_rows)}|status_counts={status_counts}|status=PASS"
)
