#!/usr/bin/env sage-python
"""Scan small integral parameters in all three smooth GF(7) MW1 disks.

The normalized marked ``A3+A4+A9`` system has one free coordinate, ``m4``.
For each of the three geometrically distinct positive GF(7) models, this
script fixes every integer in the declared interval with the required residue,
invokes the exact digit-by-digit Hensel certifier, rationally reconstructs all
forty coordinates, and substitutes them in all 53 equations over QQ.

The known reconstruction at ``m4=-20`` is additionally matched coefficient by
coefficient against the exact determinant-20 saturation-rejection certificate.
Failure of coefficientwise rational reconstruction is only a bounded negative
search result; it is not evidence that the corresponding 7-adic point is
irrational.
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
FIBRES = GEN / "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-fibre-ansatz-mod7-v1.json"
SQUARE_MARKING = GEN / "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-pole1-marking-mod7-v1.json"
NONSQUARE_MARKING = GEN / "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-pole1-marking-mod7-nonsquare-v1.json"
REJECTION = GEN / "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-source-qq-rejection-v1.json"
DEFAULT_OUTPUT = GEN / "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-rational-parameter-scan-v1.json"

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
    discriminant_core = 4 * A**3 + 27 * B**2
    node_numerator = 2 * A * N + 3 * B * C**2
    section_residual = M**2 - N**3 - A * N * C**4 - B * C**6
    equations = [a[0] + 3 * QQ(twist_integer) ** 2]
    equations.extend(discriminant_core[index] for index in range(4))
    equations.extend(discriminant_core(t + 1)[index] for index in range(5))
    equations.extend(discriminant_core[index] for index in range(15, 25))
    equations.extend(node_numerator[index] for index in range(2))
    equations.extend(m[index] for index in range(2))
    equations.extend(node_numerator[index] for index in range(10, 15))
    equations.extend(m[index] for index in range(5, 10))
    equations.extend(section_residual[index] for index in range(19))
    assert len(variables) == 40 and len(equations) == 53
    return equations


def rejection_coordinates(payload):
    model = payload["weierstrass_model"]
    section = payload["displayed_section"]
    A = [QQ(value) for value in model["A_coefficients_low_to_high"]]
    B = [QQ(value) for value in model["B_coefficients_low_to_high"]]
    C = [QQ(value) for value in section["C_coefficients_low_to_high"]]
    N = [QQ(value) for value in section["X_numerator_coefficients_low_to_high"]]
    M = [QQ(value) for value in section["Y_numerator_coefficients_low_to_high"]]
    return (
        A + [QQ(0)] * (9 - len(A))
        + B + [QQ(0)] * (13 - len(B))
        + [C[0]]
        + N + [QQ(0)] * (7 - len(N))
        + M + [QQ(0)] * (10 - len(M))
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--minimum", type=int, default=-97)
parser.add_argument("--maximum", type=int, default=99)
parser.add_argument("--lift-precision", type=int, default=40)
parser.add_argument("--workers", type=int, default=1)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
if arguments.minimum > arguments.maximum:
    parser.error("--minimum must not exceed --maximum")
if arguments.lift_precision < 2:
    parser.error("--lift-precision must be at least two")
if arguments.workers < 1:
    parser.error("--workers must be positive")

output_path = arguments.output.resolve()
sage = shutil.which("sage")
if sage is None:
    raise SystemExit("sage executable not found")
rejection = json.loads(REJECTION.read_text())
known_rejected_coordinates = rejection_coordinates(rejection)

branches = []
for label, marking_path in (("square", SQUARE_MARKING), ("nonsquare", NONSQUARE_MARKING)):
    marking = json.loads(marking_path.read_text())
    if marking["prime"] != 7:
        raise ValueError("the bounded scan is pinned to GF(7)")
    for model_rank, model in enumerate(marking["models"]):
        section = model["pole_one_sections"][0]
        padded_m = list(section["Y_numerator_coefficients_low_to_high"]) + [0] * 10
        branches.append(
            {
                "branch": label,
                "marking_path": marking_path,
                "marked_model_rank": model_rank,
                "fibre_example_index": int(model["example_index"]),
                "quadratic_twist": int(marking["quadratic_twist"]),
                "required_residue_mod_7": int(padded_m[4]) % 7,
            }
        )
assert len(branches) == 3
equations_by_twist = {
    twist: build_equations(twist) for twist in {row["quadratic_twist"] for row in branches}
}

tasks = []
for branch_rank, branch in enumerate(branches):
    for parameter in range(arguments.minimum, arguments.maximum + 1):
        if parameter % 7 == branch["required_residue_mod_7"]:
            tasks.append((branch_rank, branch, parameter))


def run_one(task):
    branch_rank, branch, parameter, temporary_directory = task
    lift_path = Path(temporary_directory) / f"lift-b{branch_rank}-m4-{parameter}.json"
    command = [
        sage,
        "-python",
        str(LIFT_SCRIPT),
        "--fibres",
        str(FIBRES),
        "--marking",
        str(branch["marking_path"]),
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
        "branch": branch["branch"],
        "fibre_example_index": branch["fibre_example_index"],
        "marked_model_rank": branch["marked_model_rank"],
        "free_parameter_m4": parameter,
    }
    if completed.returncode:
        return {
            **base,
            "status": "LIFT_ERROR",
            "stderr_tail": completed.stderr[-500:],
        }
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
    if coordinates[34] != parameter:
        return {**base, "status": "NO_FULL_RR"}
    equations = equations_by_twist[branch["quadratic_twist"]]
    if any(equation(*coordinates) for equation in equations):
        return {**base, "status": "RR_NOT_EXACT"}

    matches_rejection = coordinates == known_rejected_coordinates
    return {
        **base,
        "status": (
            "EXACT_QQ_POINT_REJECTED_PRIMITIVE_DET20"
            if matches_rejection
            else "EXACT_QQ_POINT_REQUIRES_SATURATION_AUDIT"
        ),
        "coordinates": [str(value) for value in coordinates],
        "maximum_numerator_or_denominator": int(
            max(max(abs(value.numerator()), value.denominator()) for value in coordinates)
        ),
        "matches_exact_saturation_rejection_certificate": matches_rejection,
        "primitive_closure_determinant": 20 if matches_rejection else None,
    }


with tempfile.TemporaryDirectory(
    prefix="k304b-rational-scan-", dir=GEN
) as temporary_directory:
    prepared_tasks = [task + (temporary_directory,) for task in tasks]
    with ThreadPoolExecutor(max_workers=arguments.workers) as executor:
        rows = list(executor.map(run_one, prepared_tasks))

status_counts = {}
for row in rows:
    status_counts[row["status"]] = status_counts.get(row["status"], 0) + 1
errors = [row for row in rows if row["status"] == "LIFT_ERROR"]
if errors:
    raise RuntimeError(f"fixed-parameter lifts failed: {errors[:3]}")
exact_rows = [row for row in rows if row["status"].startswith("EXACT_QQ_POINT")]

payload = {
    "schema": "elkies-k3.k3-04b-a3-a4-a9-rational-parameter-scan.v1",
    "status": "PASS_BOUNDED_EXACT_INTEGRAL_PARAMETER_SCAN",
    "inputs": {
        relative(LIFT_SCRIPT): digest(LIFT_SCRIPT),
        relative(FIBRES): digest(FIBRES),
        relative(SQUARE_MARKING): digest(SQUARE_MARKING),
        relative(NONSQUARE_MARKING): digest(NONSQUARE_MARKING),
        relative(REJECTION): digest(REJECTION),
    },
    "search": {
        "free_parameter": "m4",
        "integer_interval": [arguments.minimum, arguments.maximum],
        "branch_count": len(branches),
        "branches": [
            {
                key: value
                for key, value in row.items()
                if key not in ("marking_path",)
            }
            for row in branches
        ],
        "candidate_count": len(tasks),
        "lift_precision_exponent": arguments.lift_precision,
        "status_counts": status_counts,
    },
    "exact_rational_points": exact_rows,
    "proof_boundary": {
        "proved": (
            "Every integer m4 in the declared interval and required residue class "
            "was tested in each of the three smooth marked GF(7) disks. Every retained "
            "rational reconstruction satisfies all 53 normalized equations over Q. "
            "The determinant-20 label is used only after literal equality with the "
            "independently certified saturation-rejection model."
        ),
        "not_proved": (
            "Coefficientwise reconstruction failure does not prove 7-adic "
            "irrationality, and this integral interval is not a classification of "
            "rational points or rational parametrizations. Any new exact point would "
            "still require torsion, divisibility, and primitive-closure audits."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/scan_k3_04b_a3_a4_a9_rational_parameters.sage "
        f"--minimum {arguments.minimum} --maximum {arguments.maximum} "
        f"--lift-precision {arguments.lift_precision} --workers {arguments.workers}"
    ),
}
serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not output_path.exists() or output_path.read_text() != serialized:
        raise SystemExit(f"stale artifact: {output_path}")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    f"K304BRATSCAN|branches={len(branches)}|candidates={len(tasks)}|"
    f"exact_QQ={len(exact_rows)}|status_counts={status_counts}|status=PASS",
    flush=True,
)
