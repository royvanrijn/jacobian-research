#!/usr/bin/env sage-python
"""Probe the pole-zero NS0048 I5+I7+I2+I1* source ansatz modulo p.

Put the I1* fibre at infinity.  Its short-Weierstrass Tate valuations force
``deg(A)<=6`` and ``deg(B)<=9``; exact discriminant order seven at infinity is
checked after interpolating the I5, I7, and I2 nodal jets at ``0,1,lambda``.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

from sage.all import GF, QQ, PolynomialRing, binomial, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-b-v1.json"
POLES = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-rank1-section-poles-v1.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0048-source-ansatz-mod5.json"
SOURCE_ID = "NS0048-S030"
ORDERS = (5, 7, 2, 7)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def local_square_root(unit_coefficients, root0):
    field = root0.parent()
    answer = [field.zero()] * len(unit_coefficients)
    answer[0] = root0
    for degree in range(1, len(answer)):
        known = sum(answer[left] * answer[degree - left] for left in range(1, degree))
        answer[degree] = (unit_coefficients[degree] - known) / (2 * root0)
    return answer


def truncated_product(left, right, precision):
    field = left[0].parent()
    answer = [field.zero()] * precision
    for i, left_value in enumerate(left[:precision]):
        for j, right_value in enumerate(right[: precision - i]):
            answer[i + j] += left_value * right_value
    return answer


def branch_from_table(a_series, square_roots, sign=1):
    field = a_series[0].parent()
    unit = [-value / field(3) for value in a_series]
    root0 = square_roots.get(unit[0])
    if root0 is None or not root0:
        return None
    if sign == -1:
        root0 = -root0
    root = local_square_root(unit, root0)
    return [
        2 * value
        for value in truncated_product(
            truncated_product(root, root, len(root)), root, len(root)
        )
    ]


def taylor_coefficients(coefficients, point, precision):
    field = point.parent()
    return [
        sum(
            coefficients[index]
            * field(binomial(index, jet))
            * point ** (index - jet)
            for index in range(jet, len(coefficients))
        )
        for jet in range(precision)
    ]


def order_at(poly, point):
    if not poly:
        return None
    shifted = poly(poly.parent().gen() + point)
    return min(index for index, value in enumerate(shifted.list()) if value)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=5)
parser.add_argument("--lambda-value", type=int)
parser.add_argument("--max-a-samples-per-lambda", type=int, default=0)
parser.add_argument("--sample-stride", type=int, default=0)
parser.add_argument("--sample-offset", type=int, default=0)
parser.add_argument("--examples", type=int, default=100)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

field = GF(args.prime)
if field.characteristic() in (2, 3):
    raise SystemExit("--prime must differ from 2 and 3")
fixed_lambda = None if args.lambda_value is None else field(args.lambda_value)
if fixed_lambda is not None and fixed_lambda in (0, 1):
    raise SystemExit("--lambda-value must lie outside {0,1}")

source_payload = json.loads(SOURCE.read_text())
source_row = next(
    row
    for row in source_payload["sources"]
    if row["ns_id"] == "NS0048" and row["source_id"] == SOURCE_ID
)
source = source_row["source"]
assert source["root_type"] == "A1+A4+A6+D5"
assert source["mw_height_gram"] == [["37/14"]]
assert source["root_lattice_primitive"] and source["torsion"] == 1
assert source_row["determinant"] == 740
assert [(row["frame_id"], row["mw_rank_for_rho_19"]) for row in source_row["same_ns_high_rank_targets"]] == [("NS0048-F001", 16)]

adapted_gram = matrix(QQ, source["root_adapted_gram"])
root_gram = adapted_gram[:16, :16]
cross = adapted_gram[:16, 16]
seen = set()
component_corrections = []
for start in range(16):
    if start in seen:
        continue
    component = {start}
    stack = [start]
    while stack:
        left = stack.pop()
        for right in range(16):
            if right not in component and root_gram[left, right]:
                component.add(right)
                stack.append(right)
    seen.update(component)
    indices = sorted(component)
    block = root_gram.matrix_from_rows_and_columns(indices, indices)
    block_cross = matrix(QQ, len(indices), 1, [cross[index] for index in indices])
    correction = (block_cross.transpose() * block.inverse() * block_cross)[0, 0]
    component_corrections.append((len(indices), block.det(), correction))
assert sorted(component_corrections) == sorted(
    [(4, 5, QQ(0)), (6, 7, QQ(6) / 7), (5, 4, QQ(0)), (1, 2, QQ(1) / 2)]
)
pole_payload = json.loads(POLES.read_text())
pole = next(
    row
    for row in pole_payload["sources"]
    if row["source_artifact"] == relative(SOURCE)
    and row["source_id"] == SOURCE_ID
)
assert pole["minimum_section_pole_order"] == 0
assert pole["mw_height"] == "37/14"

ring = PolynomialRing(field, "t")
t = ring.gen()
square_roots = {value * value: value for value in field}
per_lambda_total = args.prime**6
if args.sample_stride and math.gcd(args.sample_stride, per_lambda_total) != 1:
    raise SystemExit("--sample-stride must be coprime to prime^6")


def coefficient_digits(index):
    digits = [0] * 6
    for position in range(5, -1, -1):
        digits[position] = index % args.prime
        index //= args.prime
    return tuple(digits)


lambda_values = [fixed_lambda] if fixed_lambda is not None else list(field)
lambda_values = [value for value in lambda_values if value not in (0, 1)]
examples = []
lambda_records = []
accounting = {
    "normalized_A_samples": 0,
    "branch_eligible_with_signs": 0,
    "hermite_compatible_with_signs": 0,
    "exact_prescribed_orders": 0,
    "squarefree_examples_with_signs": 0,
}

for lambda_value in lambda_values:
    rows = []
    for point, precision in (
        (field.zero(), ORDERS[0]),
        (field.one(), ORDERS[1]),
        (lambda_value, ORDERS[2]),
    ):
        for jet in range(precision):
            rows.append(
                [
                    field(binomial(index, jet)) * point ** (index - jet)
                    if index >= jet
                    else field.zero()
                    for index in range(10)
                ]
            )
    hermite = matrix(field, rows)
    if hermite.nrows() != 14 or hermite.rank() != 10:
        raise ArithmeticError("unexpected NS0048 finite-support Hermite rank")
    compatibility = hermite.left_kernel().basis_matrix()
    if compatibility.nrows() != 4:
        raise ArithmeticError("unexpected NS0048 compatibility codimension")

    local = {key: 0 for key in accounting}
    sample_limit = args.max_a_samples_per_lambda or per_lambda_total
    if args.sample_stride:
        coefficient_iterator = (
            (
                (args.sample_offset + sample_index * args.sample_stride) % per_lambda_total,
                coefficient_digits(
                    (args.sample_offset + sample_index * args.sample_stride) % per_lambda_total
                ),
            )
            for sample_index in range(sample_limit)
        )
    else:
        coefficient_iterator = itertools.islice(
            enumerate(itertools.product(range(args.prime), repeat=6)), sample_limit
        )

    for coefficient_index, digits in coefficient_iterator:
        local["normalized_A_samples"] += 1
        coefficients = [field(-3)] + [field(value) for value in digits]
        if not coefficients[6]:
            continue
        A = ring(coefficients)
        series = (
            coefficients[: ORDERS[0]],
            taylor_coefficients(coefficients, field.one(), ORDERS[1]),
            taylor_coefficients(coefficients, lambda_value, ORDERS[2]),
        )
        positive = tuple(branch_from_table(values, square_roots) for values in series)
        if any(branch is None for branch in positive):
            continue
        for signs in itertools.product((1, -1), repeat=2):
            branches = [positive[0]] + [
                [sign * value for value in branch]
                for sign, branch in zip(signs, positive[1:])
            ]
            local["branch_eligible_with_signs"] += 1
            target = vector(field, sum((list(branch) for branch in branches), []))
            if compatibility * target:
                continue
            b_coefficients = list(hermite.solve_right(target))
            local["hermite_compatible_with_signs"] += 1
            if not b_coefficients[9]:
                continue
            B = ring(b_coefficients)
            discriminant_core = 4 * A**3 + 27 * B**2
            orders = (
                order_at(discriminant_core, field.zero()),
                order_at(discriminant_core, field.one()),
                order_at(discriminant_core, lambda_value),
                24 - discriminant_core.degree(),
            )
            if orders != ORDERS:
                continue
            local["exact_prescribed_orders"] += 1
            divisor = t**ORDERS[0] * (t - 1) ** ORDERS[1] * (t - lambda_value) ** ORDERS[2]
            residual, remainder = discriminant_core.quo_rem(divisor)
            if remainder or residual.degree() != 3:
                raise ArithmeticError("unexpected NS0048 residual discriminant")
            if any(residual(point) == 0 for point in (0, 1, lambda_value)):
                continue
            if residual.gcd(residual.derivative()).degree() != 0:
                continue
            local["squarefree_examples_with_signs"] += 1
            if len(examples) < args.examples:
                examples.append(
                    {
                        "lambda": int(lambda_value),
                        "coefficient_index_within_lambda": coefficient_index,
                        "branch_signs_at_one_lambda": list(signs),
                        "A_coefficients_low_to_high": [int(value) for value in coefficients],
                        "B_coefficients_low_to_high": [int(value) for value in b_coefficients],
                        "residual_discriminant_coefficients_low_to_high": [int(value) for value in residual],
                        "geometric_fibre_profile": "I5+I7+I2+I1*+3I1",
                    }
                )
    for key in accounting:
        accounting[key] += local[key]
    exhausted = local["normalized_A_samples"] == per_lambda_total and not args.max_a_samples_per_lambda
    lambda_records.append(
        {
            "lambda": int(lambda_value),
            "samples_consumed": local["normalized_A_samples"],
            "exhausted": exhausted,
            "sample_stride": args.sample_stride or 1,
            "sample_offset": args.sample_offset,
            "accounting": local,
        }
    )

exhausted = all(record["exhausted"] for record in lambda_records)
output = {
    "schema": "elkies-k3.lattice-foundry-ns0048-source-ansatz-modp.v1",
    "status": (
        "PASS_EXACT_EXHAUSTIVE_MODULAR_SOURCE_FIBRE_ANSATZ_WITH_EXAMPLES"
        if exhausted and examples
        else "PASS_EXACT_EXHAUSTIVE_MODULAR_SOURCE_FIBRE_ANSATZ_EMPTY"
        if exhausted
        else "PASS_BOUNDED_MODULAR_SOURCE_FIBRE_ANSATZ"
        if examples
        else "PASS_BOUNDED_NO_MODULAR_SOURCE_FIBRE_ANSATZ"
    ),
    "prime": args.prime,
    "scan": {
        "lambda_slices": lambda_records,
        "exhausted": exhausted,
        "normalized_A_polynomials_per_lambda": per_lambda_total,
    },
    "accounting": accounting | {"stored_examples": len(examples)},
    "ansatz": {
        "short_weierstrass": "y^2=x^3+A(t)x+B(t)",
        "degree_bounds_from_I1star_at_infinity": {"A": 6, "B": 9},
        "normalization": "A(0)=-3; supports 0:I5, 1:I7, lambda:I2, infinity:I1*",
        "finite_nodal_Hermite_conditions": 14,
        "B_coefficient_rank": 10,
        "finite_support_compatibility_equations": 4,
        "exact_infinity_discriminant_order": 7,
        "expected_fibre_stratum_dimension": 2,
        "minimum_section_pole_order_to_impose": 0,
        "section_component_corrections": {
            "0:I5": "0",
            "1:I7": "6/7",
            "lambda:I2": "1/2",
            "infinity:I1*": "0",
        },
    },
    "examples": examples,
    "source": {
        "artifact": relative(SOURCE),
        "artifact_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "source_id": SOURCE_ID,
        "source_gram_sha256": source["gram_sha256"],
        "section_pole_artifact": relative(POLES),
        "section_pole_artifact_sha256": hashlib.sha256(POLES.read_bytes()).hexdigest(),
        "same_ns_high_rank_targets": source_row["same_ns_high_rank_targets"],
    },
    "proof_boundary": {
        "proved": (
            "Every stored model has the exact displayed finite-field fibre profile."
            if examples
            else "The declared finite coefficient slices were checked exactly."
        ),
        "not_proved": (
            "The pole-zero MW section, NS0048 marking, characteristic-zero lift, "
            "rational parameterization, and neighbour corridor remain open."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/probe_lattice_foundry_ns0048_source_ansatz_modp.sage"
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
output_path = args.output.resolve()
if args.check:
    if output_path.read_text() != serialized:
        raise SystemExit("NS0048 source-ansatz artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "FOUNDRYNS0048ANSATZ|"
    f"p={args.prime}|lambdas={len(lambda_records)}|samples={accounting['normalized_A_samples']}|"
    f"compatible={accounting['hermite_compatible_with_signs']}|"
    f"squarefree={accounting['squarefree_examples_with_signs']}|"
    f"exhausted={int(exhausted)}|status={'PASS' if exhausted or examples else 'BOUNDED_NEGATIVE'}",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
