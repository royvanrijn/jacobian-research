#!/usr/bin/env sage-python
"""Exhaustively scan reduced NS0011 pole-two section systems over GF(p).

Polynomial functions on ``GF(p)^n`` are first reduced by ``x^p=x``.  A
tensor-product evaluation transform then evaluates every point of the affine
space exactly.  The default five systems cover every nonzero-y infinity chart
of the three stored GF(5) fibre models.  Zero-y infinity points remain outside
the scan and are recorded explicitly.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from sage.all import GF, PolynomialRing


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SYSTEM_DIR = ROOT / "artifacts/local/elkies-k3/ns0011-pole2-section-modp"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0011-pole2-sections-mod5.json"
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def function_exponent(exponent: int, prime: int) -> int:
    if exponent == 0:
        return 0
    return 1 + (exponent - 1) % (prime - 1)


def evaluate_everywhere(polynomial, prime: int) -> np.ndarray:
    """Return the exact polynomial-function table on GF(p)^n."""
    variable_count = polynomial.parent().ngens()
    coefficients = np.zeros((prime,) * variable_count, dtype=np.int64)
    for exponents, coefficient in polynomial.dict().items():
        reduced = tuple(function_exponent(int(value), prime) for value in exponents)
        coefficients[reduced] += int(coefficient)
    coefficients %= prime
    transform = np.array(
        [[pow(value, exponent, prime) for exponent in range(prime)] for value in range(prime)],
        dtype=np.int64,
    )
    values = coefficients
    for axis in range(variable_count):
        values = np.tensordot(transform, values, axes=(1, axis)) % prime
        values = np.moveaxis(values, 0, axis)
    return values


def parse_msolve(path: Path):
    lines = path.read_text().splitlines()
    names = lines[0].split(",")
    prime = int(lines[1])
    ring = PolynomialRing(GF(prime), names=names, order="degrevlex")
    equation_text = "\n".join(lines[2:]).replace("^", "**")
    pieces = [piece.strip() for piece in equation_text.split(",") if piece.strip()]
    return prime, ring, [ring(piece) for piece in pieces]


def square_root_coefficients(coefficients: list[int], prime: int):
    """Return one polynomial square root, or None, in odd characteristic."""
    degree = len(coefficients) - 1
    while degree >= 0 and coefficients[degree] % prime == 0:
        degree -= 1
    if degree < 0:
        return [0]
    if degree % 2:
        return None
    root_degree = degree // 2
    leading_roots = [
        value for value in range(prime) if value * value % prime == coefficients[degree] % prime
    ]
    if not leading_roots:
        return None
    root = [0] * (root_degree + 1)
    root[root_degree] = leading_roots[0]
    inverse_leading = pow(2 * root[root_degree], -1, prime)
    for index in range(root_degree - 1, -1, -1):
        target_degree = root_degree + index
        known = 0
        for left in range(root_degree + 1):
            right = target_degree - left
            if not 0 <= right <= root_degree:
                continue
            if (left, right) in ((root_degree, index), (index, root_degree)):
                continue
            known += root[left] * root[right]
        root[index] = (coefficients[target_degree] - known) * inverse_leading % prime
    square = [0] * (2 * root_degree + 1)
    for left, left_value in enumerate(root):
        for right, right_value in enumerate(root):
            square[left + right] = (square[left + right] + left_value * right_value) % prime
    if any(
        square[index] != coefficients[index] % prime
        for index in range(max(len(square), degree + 1))
    ):
        return None
    return root


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--system-dir", type=Path, default=DEFAULT_SYSTEM_DIR)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

system_dir = args.system_dir.resolve()
metadata_paths = sorted(system_dir.glob("example*-elimR09.json"))
zero_y_metadata_paths = sorted(system_dir.glob("example*-zeroY*-elimR00.json"))
if not metadata_paths:
    raise SystemExit(f"no fully reduced system metadata under {system_dir}")

records = []
for metadata_path in metadata_paths:
    metadata = json.loads(metadata_path.read_text())
    if metadata.get("schema") != "elkies-k3.lattice-foundry-ns0011-pole2-section-modp-system.v1":
        continue
    msolve_path = ROOT / metadata["system"]["msolve_input"]
    if digest(msolve_path) != metadata["system"]["msolve_sha256"]:
        raise ArithmeticError(f"stale msolve digest: {msolve_path}")
    prime, ring, equations = parse_msolve(msolve_path)
    if len(equations) != metadata["system"]["equation_count"]:
        raise ArithmeticError("equation accounting mismatch")
    solution_mask = np.ones((prime,) * ring.ngens(), dtype=bool)
    evaluation_order = sorted(
        range(len(equations)), key=lambda index: len(equations[index].monomials())
    )
    surviving_counts = []
    for equation_index in evaluation_order:
        solution_mask &= evaluate_everywhere(equations[equation_index], prime) == 0
        surviving_counts.append(int(solution_mask.sum()))
        if not surviving_counts[-1]:
            break
    solution_tuples = [tuple(int(value) for value in row) for row in np.argwhere(solution_mask)]

    # Materialize and independently verify every displayed section against the
    # original homogeneous Weierstrass identity.
    outer = PolynomialRing(ring, "t")
    t = outer.gen()
    input_payload = json.loads((ROOT / metadata["input"]["artifact"]).read_text())
    example = input_payload["examples"][metadata["input"]["example_index"]]
    field = GF(prime)
    base = PolynomialRing(field, "t")
    A = base(example["A_coefficients_low_to_high"])
    B = base(example["B_coefficients_low_to_high"])
    Z_symbolic = outer(metadata["system"]["Z"])
    X_symbolic = outer(metadata["system"]["X"])
    R_symbolic = outer(metadata["system"]["R"])

    materialized = []
    for values in solution_tuples:
        def specialize(poly):
            return base([field(coefficient(*values)) for coefficient in poly.list()])

        Z = specialize(Z_symbolic)
        X = specialize(X_symbolic)
        R = specialize(R_symbolic)
        Y = t.parent()(t**2 * (t - 1)) * R
        if Y**2 != X**3 + A * X * Z**4 + B * Z**6:
            raise ArithmeticError("materialized section fails the Weierstrass identity")
        materialized.append(
            {
                "variables_in_declared_order": list(values),
                "Z_coefficients_low_to_high": [int(value) for value in Z],
                "X_coefficients_low_to_high": [int(value) for value in X],
                "Y_coefficients_low_to_high": [int(value) for value in Y],
            }
        )

    records.append(
        {
            "system_metadata": relative(metadata_path),
            "system_metadata_sha256": digest(metadata_path),
            "msolve_input": relative(msolve_path),
            "msolve_sha256": digest(msolve_path),
            "example_index": metadata["input"]["example_index"],
            "infinity_branch": metadata["infinity"]["selected_nonzero_y_branch"],
            "infinity_x_y": metadata["infinity"]["selected_x_y"],
            "zero_y_points_excluded": metadata["infinity"][
                "zero_y_points_excluded_from_this_chart"
            ],
            "affine_points_scanned": prime ** ring.ngens(),
            "equation_evaluation_order": evaluation_order,
            "surviving_counts_after_each_equation": surviving_counts,
            "solution_count": len(materialized),
            "solutions": materialized,
        }
    )

# The smooth zero-y point on each nodal I3 cubic is an identity-component
# point, but the leading triangular coefficient vanishes there.  Scan the
# seven base variables directly and test whether the resulting H(t) is a
# polynomial square.  The nodal zero-y point is correctly excluded: after
# resolution it lies on a nonidentity I3 component.
zero_y_records = []
for metadata_path in zero_y_metadata_paths:
    metadata = json.loads(metadata_path.read_text())
    names = metadata["system"]["base_variables_before_R"]
    prime = int(metadata["prime"])
    ring = PolynomialRing(GF(prime), names=names, order="degrevlex")
    outer = PolynomialRing(ring, "t")
    H = outer(metadata["system"]["H"])
    coefficient_tables = [
        evaluate_everywhere(ring(H[index]), prime).reshape(-1)
        for index in range(19)
    ]
    solution_data = []
    for flat_index in range(prime ** len(names)):
        coefficients = [int(table[flat_index]) for table in coefficient_tables]
        root = square_root_coefficients(coefficients, prime)
        if root is None:
            continue
        values = tuple(
            int(value)
            for value in np.unravel_index(flat_index, (prime,) * len(names))
        )
        field = GF(prime)
        base = PolynomialRing(field, "t")
        tb = base.gen()

        def specialize(poly):
            return base([field(coefficient(*values)) for coefficient in poly.list()])

        Z = specialize(outer(metadata["system"]["Z"]))
        X = specialize(outer(metadata["system"]["X"]))
        R = base(root)
        Y = tb**2 * (tb - 1) * R
        input_payload = json.loads((ROOT / metadata["input"]["artifact"]).read_text())
        example = input_payload["examples"][metadata["input"]["example_index"]]
        A = base(example["A_coefficients_low_to_high"])
        B = base(example["B_coefficients_low_to_high"])
        if Y**2 != X**3 + A * X * Z**4 + B * Z**6:
            raise ArithmeticError("zero-y materialized section fails the identity")
        solution_data.append(
            {
                "variables_in_declared_order": list(values),
                "Z_coefficients_low_to_high": [int(value) for value in Z],
                "X_coefficients_low_to_high": [int(value) for value in X],
                "Y_coefficients_low_to_high": [int(value) for value in Y],
            }
        )
    zero_y_records.append(
        {
            "system_metadata": relative(metadata_path),
            "system_metadata_sha256": digest(metadata_path),
            "example_index": metadata["input"]["example_index"],
            "infinity_x_y": metadata["infinity"]["selected_x_y"],
            "affine_points_scanned": prime ** len(names),
            "solution_count": len(solution_data),
            "solutions": solution_data,
        }
    )

declared_branch_counts = {}
for path in metadata_paths:
    record = json.loads(path.read_text())
    declared_branch_counts[record["input"]["example_index"]] = len(
        record["infinity"]["nonzero_y_branches_up_to_sign"]
    )
expected_charts = sum(declared_branch_counts.values())
if len(records) != expected_charts:
    raise ArithmeticError(
        f"found {len(records)} charts but metadata declares {expected_charts}"
    )
total_solutions = sum(record["solution_count"] for record in records + zero_y_records)
output_payload = {
    "schema": "elkies-k3.lattice-foundry-ns0011-pole2-sections-modp-scan.v1",
    "status": (
        "PASS_EXACT_EXHAUSTIVE_NONZERO_Y_MODULAR_SECTION_SCAN_WITH_SOLUTIONS"
        if total_solutions
        else "PASS_EXACT_EXHAUSTIVE_NONZERO_Y_MODULAR_SECTION_SCAN_EMPTY"
    ),
    "prime": records and json.loads((ROOT / records[0]["system_metadata"]).read_text())["prime"],
    "scope": {
        "stored_fibre_models": len({record["example_index"] for record in records}),
        "nonzero_y_infinity_charts_up_to_section_negation": len(records),
        "smooth_zero_y_infinity_charts": len(zero_y_records),
        "affine_points_per_chart": records[0]["affine_points_scanned"],
        "nodal_zero_y_points_excluded_by_I3_identity_component": True,
    },
    "accounting": {
        "total_nonzero_y_chart_solutions": total_solutions,
        "charts_with_solutions": sum(
            bool(record["solution_count"]) for record in records + zero_y_records
        ),
    },
    "charts": records,
    "smooth_zero_y_charts": zero_y_records,
    "proof_boundary": {
        "proved": (
            "The tensor-product calculation exhausts every GF(5) point of every "
            "nonzero-y and smooth zero-y infinity chart exported for the seven "
            "stored fibre models. The nodal zero-y point is excluded by the I3 "
            "identity-component condition. "
            "Each displayed solution is independently checked in the homogeneous "
            "Weierstrass equation."
        ),
        "not_proved": (
            "A characteristic-zero lift, a rational source parameterization, "
            "and the neighbour route are not covered."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/scan_lattice_foundry_ns0011_pole2_sections_modp.sage"
    ),
}
serialized = json.dumps(output_payload, indent=2, sort_keys=True) + "\n"
output_path = args.output.resolve()
if args.check:
    if output_path.read_text() != serialized:
        raise SystemExit("NS0011 pole-two modular section scan is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "FOUNDRYNS0011SECTIONSCAN|"
    f"charts={len(records) + len(zero_y_records)}|"
    f"points={sum(row['affine_points_scanned'] for row in records + zero_y_records)}|"
    f"solutions={total_solutions}|status=PASS",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
