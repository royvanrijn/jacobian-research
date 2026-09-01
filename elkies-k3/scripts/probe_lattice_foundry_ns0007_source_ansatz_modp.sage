#!/usr/bin/env sage-python
"""Probe the NS0007 A1+A3+2A6 fibre ansatz over a finite field.

The selected MW1 source has semistable profile ``I2+I4+2I7+4I1`` and an
exact pole-zero generator.  Supports are normalized to ``0,1,lambda,infinity``.
For every lambda and normalized degree-eight A polynomial, the script solves
the twenty branch jets for the thirteen B coefficients and retains exact
squarefree residual discriminants.

This gate imposes the reducible fibres only.  It does not yet impose the
polynomial MW section, identify the NS0007 marking, or lift to characteristic
zero.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, binomial, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_POLES = (
    ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-rank1-section-poles-v1.json"
)
CANDIDATES = {
    "ns0007": {
        "ns_id": "NS0007",
        "source_id": "NS0007-S025",
        "source_file": "elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-c-v1.json",
        "root_type": "A1+A3+2A6",
        "height": "11/4",
        "orders": (2, 4, 7, 7),
        "profile": "I2+I4+2I7+4I1",
        "corrections": ["1/2", "3/4", "0", "0"],
    },
    "ns0034": {
        "ns_id": "NS0034",
        "source_id": "NS0034-S008",
        "source_file": "elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-a-v1.json",
        "root_type": "A2+A3+A4+A7",
        "height": "19/8",
        "orders": (4, 8, 3, 5),
        "profile": "I4+I8+I3+I5+4I1",
        "corrections": ["3/4", "7/8", "0", "0"],
    },
}


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


def multiplicative_branch(a_series, sign=1):
    field = a_series[0].parent()
    unit = [-value / field(3) for value in a_series]
    if not unit[0] or not unit[0].is_square():
        return None
    root0 = unit[0].sqrt()
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
            coefficients[index] * field(binomial(index, jet)) * point ** (index - jet)
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
parser.add_argument("--candidate", choices=sorted(CANDIDATES), default="ns0007")
parser.add_argument("--source", type=Path)
parser.add_argument("--source-id")
parser.add_argument("--section-poles", type=Path, default=DEFAULT_POLES)
parser.add_argument("--output", type=Path)
parser.add_argument("--prime", type=int, default=5)
parser.add_argument("--examples", type=int, default=20)
parser.add_argument(
    "--max-a-samples-per-lambda",
    type=int,
    default=0,
    help="truncate each lambda slice; zero exhausts all normalized A polynomials",
)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

candidate = CANDIDATES[args.candidate]
source_path = (
    args.source.resolve()
    if args.source
    else ROOT / "artifacts/generated-results" / candidate["source_file"]
)
pole_path = args.section_poles.resolve()
source_id = args.source_id or candidate["source_id"]
output_path = (
    args.output.resolve()
    if args.output
    else ROOT
    / "artifacts/generated-results"
    / f"elkies-k3-lattice-foundry-{args.candidate}-source-ansatz-mod{args.prime}.json"
)
source_payload = json.loads(source_path.read_text())
source_entry = next(
    row
    for row in source_payload["sources"]
    if row["ns_id"] == candidate["ns_id"] and row["source_id"] == source_id
)
source = source_entry["source"]
assert source["root_type"] == candidate["root_type"]
assert source["root_lattice_primitive"] and source["torsion"] == 1
assert source["mw_height_gram"] == [[candidate["height"]]]
pole_payload = json.loads(pole_path.read_text())
pole_row = next(
    row
    for row in pole_payload["sources"]
    if row["source_artifact"] == relative(source_path)
    and row["source_id"] == source_id
)
assert pole_row["minimum_section_pole_order"] == 0

field = GF(args.prime)
if field.characteristic() in (2, 3):
    raise SystemExit("--prime must differ from 2 and 3")
ring = PolynomialRing(field, "t")
t = ring.gen()

examples = []
accounting = {
    "normalized_A_samples": 0,
    "branch_eligible_with_signs": 0,
    "hermite_compatible_with_signs": 0,
    "exact_prescribed_orders": 0,
    "squarefree_examples_with_signs": 0,
}
lambda_records = []
per_lambda_total = args.prime**8
for lambda_value in field:
    if lambda_value in (field.zero(), field.one()):
        continue
    rows = []
    order_zero, order_one, order_lambda, order_infinity = candidate["orders"]
    for point, precision in (
        (field.zero(), order_zero),
        (field.one(), order_one),
        (lambda_value, order_lambda),
    ):
        for jet in range(precision):
            rows.append(
                [
                    field(binomial(index, jet)) * point ** (index - jet)
                    if index >= jet
                    else field.zero()
                    for index in range(13)
                ]
            )
    for jet in range(order_infinity):
        rows.append([field(index == 12 - jet) for index in range(13)])
    hermite = matrix(field, rows)
    if hermite.nrows() != 20 or hermite.ncols() != hermite.rank():
        raise ArithmeticError("unexpected Hermite rank")
    compatibility = hermite.left_kernel().basis_matrix()
    if compatibility.nrows() != 7:
        raise ArithmeticError("unexpected compatibility codimension")

    local = {key: 0 for key in accounting}
    samples = 0
    for digits in itertools.product(range(args.prime), repeat=8):
        samples += 1
        if args.max_a_samples_per_lambda and samples > args.max_a_samples_per_lambda:
            samples -= 1
            break
        local["normalized_A_samples"] += 1
        coefficients = [field(-3)] + [field(value) for value in digits]
        if not coefficients[8]:
            continue
        A = ring(coefficients)
        series = (
            coefficients[:order_zero],
            taylor_coefficients(coefficients, field.one(), order_one),
            taylor_coefficients(coefficients, lambda_value, order_lambda),
            [coefficients[8 - jet] for jet in range(order_infinity)],
        )
        positive = tuple(multiplicative_branch(values, 1) for values in series)
        if any(branch is None for branch in positive):
            continue
        for signs in itertools.product((1, -1), repeat=3):
            branches = (
                positive[0],
                [signs[0] * value for value in positive[1]],
                [signs[1] * value for value in positive[2]],
                [signs[2] * value for value in positive[3]],
            )
            local["branch_eligible_with_signs"] += 1
            target = vector(field, sum((list(branch) for branch in branches), []))
            if compatibility * target:
                continue
            b_coefficients = list(hermite.solve_right(target))
            local["hermite_compatible_with_signs"] += 1
            B = ring(b_coefficients)
            discriminant_core = 4 * A**3 + 27 * B**2
            orders = (
                order_at(discriminant_core, field.zero()),
                order_at(discriminant_core, field.one()),
                order_at(discriminant_core, lambda_value),
                24 - discriminant_core.degree(),
            )
            if orders != candidate["orders"]:
                continue
            local["exact_prescribed_orders"] += 1
            divisor = (
                t**order_zero
                * (t - 1) ** order_one
                * (t - lambda_value) ** order_lambda
            )
            residual, remainder = discriminant_core.quo_rem(divisor)
            if remainder or residual.degree() != 4:
                raise ArithmeticError("unexpected residual discriminant")
            if any(residual(point) == 0 for point in (0, 1, lambda_value)):
                continue
            if residual.gcd(residual.derivative()).degree() != 0:
                continue
            local["squarefree_examples_with_signs"] += 1
            if len(examples) < args.examples:
                examples.append(
                    {
                        "lambda": int(lambda_value),
                        "sample_index_within_lambda": samples,
                        "branch_signs_at_one_lambda_infinity": list(signs),
                        "A_coefficients_low_to_high": [int(value) for value in coefficients],
                        "B_coefficients_low_to_high": [int(value) for value in b_coefficients],
                        "residual_discriminant_coefficients_low_to_high": [
                            int(value) for value in residual
                        ],
                        "geometric_fibre_profile": candidate["profile"],
                    }
                )
    for key in accounting:
        accounting[key] += local[key]
    lambda_records.append(
        {
            "lambda": int(lambda_value),
            "samples_consumed": samples,
            "exhausted": samples == per_lambda_total and not args.max_a_samples_per_lambda,
            "accounting": local,
        }
    )

exhausted = all(record["exhausted"] for record in lambda_records)
output = {
    "schema": f"elkies-k3.lattice-foundry-{args.candidate}-source-ansatz-modp.v1",
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
        "degree_bounds": {"A": 8, "B": 12},
        "normalization": "A(0)=-3; supports at 0,1,lambda,infinity",
        "normalized_reducible_supports": [
            f"0:I{order_zero}",
            f"1:I{order_one}",
            f"lambda:I{order_lambda}",
            f"infinity:I{order_infinity}",
        ],
        "hermite_conditions": 20,
        "B_coefficient_rank": 13,
        "compatibility_equations_on_A": 7,
        "expected_fibre_stratum_dimension": 2,
        "minimum_section_pole_order_to_impose": 0,
        "section_component_corrections": candidate["corrections"],
    },
    "examples": examples,
    "source": {
        "artifact": relative(source_path),
        "artifact_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "source_id": source_id,
        "source_gram_sha256": source["gram_sha256"],
        "section_pole_artifact": relative(pole_path),
        "section_pole_artifact_sha256": hashlib.sha256(pole_path.read_bytes()).hexdigest(),
    },
    "proof_boundary": {
        "proved": (
            "The scan exhausts the displayed normalized finite-field fibre chart "
            "and proves that it contains no squarefree model."
            if exhausted and not examples
            else "Every stored example has the exact displayed finite-field fibre profile."
            if examples
            else "The scan exactly checks the declared bounded coefficient slices."
        ),
        "not_proved": (
            f"The pole-zero MW section, {candidate['ns_id']} marking, rational parameterization, "
            "characteristic-zero lift, and neighbour route are not proved."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/probe_lattice_foundry_ns0007_source_ansatz_modp.sage "
        f"--candidate {args.candidate}"
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if output_path.read_text() != serialized:
        raise SystemExit(f"{candidate['ns_id']} modular source-ansatz artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    f"FOUNDRY{candidate['ns_id']}ANSATZ|"
    f"p={args.prime}|lambdas={len(lambda_records)}|"
    f"samples={accounting['normalized_A_samples']}|"
    f"compatible={accounting['hermite_compatible_with_signs']}|"
    f"squarefree={accounting['squarefree_examples_with_signs']}|"
    f"exhausted={int(exhausted)}|"
    f"status={'PASS' if exhausted or examples else 'BOUNDED_NEGATIVE'}",
    flush=True,
)
