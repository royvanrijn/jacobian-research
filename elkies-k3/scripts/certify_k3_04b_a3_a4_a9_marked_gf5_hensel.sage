#!/usr/bin/env sage-python
"""Certify the tangent and finite Hensel gate for a marked pole-one MW1 seed."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FIBRES = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-fibre-ansatz-mod5-v1.json"
)
DEFAULT_MARKING = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-pole1-marking-mod5-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-marked-gf5-hensel-v1.json"
)
SCHEMA = "elkies-k3.k3-04b-a3-a4-a9-marked-gfp-hensel.v1"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def padded(values, length):
    if len(values) > length:
        raise ValueError("coefficient vector exceeds degree bound")
    return list(values) + [0] * (length - len(values))


NAMES = (
    [f"a{index}" for index in range(9)]
    + [f"b{index}" for index in range(13)]
    + ["c0"]
    + [f"n{index}" for index in range(7)]
    + [f"m{index}" for index in range(10)]
)


def fibre_orders(root_type):
    orders = []
    for term in root_type.split("+"):
        match = re.fullmatch(r"(?:(\d+))?A(\d+)", term)
        if match is None:
            raise ValueError("source is not semistable")
        orders.extend([int(match.group(2)) + 1] * int(match.group(1) or 1))
    return orders


def build_system(base_ring, orders, depths, normalization_scale=1):
    coefficient_ring = PolynomialRing(base_ring, names=NAMES, order="degrevlex")
    variables = list(coefficient_ring.gens())
    cursor = 0

    def take(length):
        nonlocal cursor
        block = variables[cursor : cursor + length]
        cursor += length
        return block

    a = take(9)
    b = take(13)
    c0 = take(1)[0]
    n = take(7)
    m = take(10)
    if cursor != len(variables) or len(variables) != 40:
        raise ArithmeticError("unexpected variable count")

    function_ring = PolynomialRing(coefficient_ring, "t")
    t = function_ring.gen()

    def polynomial(coefficients):
        return sum(value * t**index for index, value in enumerate(coefficients))

    A = polynomial(a)
    B = polynomial(b)
    C = t + c0
    N = polynomial(n)
    M = polynomial(m)
    discriminant_core = 4 * A**3 + 27 * B**2
    node_numerator = 2 * A * N + 3 * B * C**2
    section_residual = M**2 - N**3 - A * N * C**4 - B * C**6
    if len(orders) != 3 or len(depths) != 3:
        raise ValueError("marked system requires three supports")
    blocks = {
        "normalization": [a[0] + 3 * base_ring(normalization_scale)],
        "fibre_at_zero": [discriminant_core[index] for index in range(orders[0])],
        "fibre_at_one": [discriminant_core(t + 1)[index] for index in range(orders[1])],
        "fibre_at_infinity": [
            discriminant_core[index] for index in range(25 - orders[2], 25)
        ],
        "component_marking": (
            [node_numerator[index] for index in range(depths[0])]
            + [m[index] for index in range(depths[0])]
            + [node_numerator(t + 1)[index] for index in range(depths[1])]
            + [M(t + 1)[index] for index in range(depths[1])]
            + [node_numerator[index] for index in range(15 - depths[2], 15)]
            + [m[index] for index in range(10 - depths[2], 10)]
        ),
        "pole_one_section": [section_residual[index] for index in range(19)],
    }
    equations = [equation for block in blocks.values() for equation in block]
    expected_count = 1 + sum(orders) + 2 * sum(depths) + 19
    if len(equations) != expected_count:
        raise ArithmeticError("unexpected equation count")
    return variables, blocks, equations


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--fibres", type=Path, default=DEFAULT_FIBRES)
parser.add_argument("--marking", type=Path, default=DEFAULT_MARKING)
parser.add_argument("--schema", default=SCHEMA)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--lift-precision", type=int, default=8)
parser.add_argument("--marked-model-rank", type=int, default=0)
parser.add_argument("--marked-section-index", type=int, default=0)
parser.add_argument(
    "--free-parameter-integer",
    type=int,
    help="fix the unique nonpivot coordinate to this congruent integer",
)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
if arguments.lift_precision < 1:
    parser.error("lift precision must be positive")

fibres_path = arguments.fibres.resolve()
marking_path = arguments.marking.resolve()
output_path = arguments.output.resolve()
legacy_profile = (
    fibres_path == DEFAULT_FIBRES.resolve()
    and marking_path == DEFAULT_MARKING.resolve()
    and arguments.schema == SCHEMA
)
fibres = json.loads(fibres_path.read_text())
marking = json.loads(marking_path.read_text())
prime = int(fibres.get("prime"))
if prime in (2, 3) or not fibres["scan"]["exhausted"]:
    raise ValueError("Hensel gate requires an exhaustive good-prime fibre artifact")
if marking.get("status") != (
    "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_WITH_MARKED_MW1_SECTION"
):
    raise ValueError("Hensel gate requires a positive marking artifact")
if int(marking.get("prime")) != prime:
    raise ValueError("fibre and marking primes differ")
orders = fibre_orders(marking["source"]["root_type"])
depths = marking["source"].get(
    "component_depths_at_normalized_supports",
    marking["source"].get("component_depths_at_I4_I5_I10"),
)
if len(orders) != 3 or sum(order - 1 for order in orders) != 16:
    raise ValueError("marked seed is not a three-support MW1 source")
if depths is None or len(depths) != 3:
    raise ValueError("marking artifact omits the component depths")
expected_supports = [f"0:I{orders[0]}", f"1:I{orders[1]}", f"infinity:I{orders[2]}"]
if fibres["ansatz"]["normalized_reducible_supports"] != expected_supports:
    raise ValueError("fibre and marking support profiles differ")

if not 0 <= arguments.marked_model_rank < len(marking["models"]):
    parser.error("marked model rank is outside the positive marking inventory")
marked_model = marking["models"][arguments.marked_model_rank]
example_index = int(marked_model["example_index"])
if not 0 <= arguments.marked_section_index < len(marked_model["pole_one_sections"]):
    parser.error("marked section index is outside the selected model")
section = marked_model["pole_one_sections"][arguments.marked_section_index]
fibre = fibres["examples"][example_index]

field = GF(prime)
twist_integer = int(marking["quadratic_twist"])
variables, equation_blocks, equations = build_system(
    field, orders, depths, field(twist_integer) ** 2
)
twist = field(marking["quadratic_twist"])
point_values = (
    [int(twist**2 * field(value)) for value in padded(fibre["A_coefficients_low_to_high"], 9)]
    + [int(twist**3 * field(value)) for value in padded(fibre["B_coefficients_low_to_high"], 13)]
    + [section["C_coefficients_low_to_high"][0]]
    + padded(section["X_numerator_coefficients_low_to_high"], 7)
    + padded(section["Y_numerator_coefficients_low_to_high"], 10)
)
point = [field(value) for value in point_values]
if len(point) != len(variables):
    raise ArithmeticError("marked point has wrong coordinate count")


def evaluate(polynomial_value):
    return field(polynomial_value(*point))


if any(evaluate(equation) for equation in equations):
    raise ArithmeticError("marked point does not solve the integral lift system")
jacobian = matrix(
    field,
    [
        [evaluate(equation.derivative(variable)) for variable in variables]
        for equation in equations
    ],
)
rank = int(jacobian.rank())
pivot_columns = [int(index) for index in jacobian.pivots()]
pivot_rows = [int(index) for index in jacobian.transpose().pivots()]
free_columns = [index for index in range(len(variables)) if index not in pivot_columns]
if arguments.free_parameter_integer is not None:
    if len(free_columns) != 1:
        parser.error("fixed free parameter requires a one-dimensional tangent")
    if arguments.free_parameter_integer % prime != int(point[free_columns[0]]):
        parser.error("fixed free parameter is not congruent to the marked seed")
minor = (
    int(jacobian.matrix_from_rows_and_columns(pivot_rows, pivot_columns).det())
    if len(pivot_columns) == len(pivot_rows) == rank
    else 0
)
block_ranks = {}
row_count = 0
for name, block in equation_blocks.items():
    row_count += len(block)
    block_ranks[name] = int(jacobian[:row_count, :].rank())

integer_variables, unused_blocks, integer_equations = build_system(
    ZZ, orders, depths, twist_integer**2
)
lift_coordinates = [int(value) for value in point]
lift_steps = []
lift_failure = None
for exponent in range(1, arguments.lift_precision):
    modulus = prime**exponent
    values = [ZZ(equation(*lift_coordinates)) for equation in integer_equations]
    if any(value % modulus for value in values):
        raise ArithmeticError("coordinates lost their certified precision")
    right = vector(field, [-(value // modulus) for value in values])
    try:
        if arguments.free_parameter_integer is None:
            correction = jacobian.solve_right(right)
        else:
            free_index = free_columns[0]
            free_digit = field(
                (arguments.free_parameter_integer - lift_coordinates[free_index])
                // modulus
            )
            pivot_matrix = jacobian.matrix_from_rows_and_columns(
                pivot_rows, pivot_columns
            )
            pivot_right = vector(field, [right[index] for index in pivot_rows])
            pivot_right -= vector(
                field,
                [jacobian[index, free_index] * free_digit for index in pivot_rows],
            )
            pivot_correction = pivot_matrix.solve_right(pivot_right)
            entries = [field.zero()] * len(variables)
            entries[free_index] = free_digit
            for index, value in zip(pivot_columns, pivot_correction):
                entries[index] = value
            correction = vector(field, entries)
            if jacobian * correction != right:
                raise ValueError("fixed-parameter correction violates dependent rows")
    except ValueError:
        lift_failure = {
            "failed_lift_from_exponent": exponent,
            "right_hand_side_mod_prime": [int(value) for value in right],
        }
        break
    next_modulus = modulus * prime
    lift_coordinates = [
        int((value + modulus * int(delta)) % next_modulus)
        for value, delta in zip(lift_coordinates, correction)
    ]
    if any(
        ZZ(equation(*lift_coordinates)) % next_modulus
        for equation in integer_equations
    ):
        raise ArithmeticError("linear correction failed the next Hensel digit")
    lift_steps.append(
        {
            "precision_exponent": exponent + 1,
            "modulus": next_modulus,
            "maximum_centered_correction_digit": max(
                min(int(value), prime - int(value)) for value in correction
            ),
        }
    )

achieved = 1 + len(lift_steps)
lift_modulus = prime**achieved
if any(
    ZZ(equation(*lift_coordinates)) % lift_modulus
    for equation in integer_equations
):
    raise ArithmeticError("reported lift has nonzero residuals")
status = (
    f"PASS_ONE_DIMENSIONAL_MARKED_TANGENT_AND_EXPLICIT_Z{prime}_LIFT_TO_REQUESTED_PRECISION"
    if rank == len(variables) - 1 and minor and lift_failure is None
    else (
        "PASS_ONE_DIMENSIONAL_MARKED_TANGENT_BUT_FINITE_LIFT_STOPPED"
        if rank == len(variables) - 1 and minor
        else "PASS_POINT_BUT_EXPECTED_MARKED_TANGENT_DIMENSION_FAILED"
    )
)
payload = {
    "schema": arguments.schema,
    "status": status,
    "prime": prime,
    "inputs": {
        relative(fibres_path): digest(fibres_path),
        relative(marking_path): digest(marking_path),
    },
    "seed": {
        "fibre_example_index": example_index,
        "marked_model_rank": arguments.marked_model_rank,
        "marked_section_index": arguments.marked_section_index,
        "coordinate_names": NAMES,
        "coordinates_mod_prime": [int(value) for value in point],
    },
    "system": {
        "normalization": (
            f"A(0)=-3*{twist_integer}^2; "
            + ",".join(f"I{order}" for order in orders)
            + " at 0,1,infinity; C monic"
        ),
        "variable_count": len(variables),
        "equation_count": len(equations),
        "equation_block_sizes": {
            name: len(block) for name, block in equation_blocks.items()
        },
        "equation_block_cumulative_jacobian_ranks": block_ranks,
        "all_equation_residuals_zero_mod_prime": True,
    },
    "jacobian_certificate": {
        "rank_mod_prime": rank,
        "tangent_dimension": len(variables) - rank,
        "pivot_row_indices": pivot_rows,
        "pivot_column_indices": pivot_columns,
        "pivot_variable_names": [NAMES[index] for index in pivot_columns],
        "omitted_free_variable_names": [
            name for index, name in enumerate(NAMES) if index not in pivot_columns
        ],
        "pivot_minor_determinant_mod_prime": minor,
    },
    "finite_precision_lift": {
        "requested_precision_exponent": arguments.lift_precision,
        "achieved_precision_exponent": achieved,
        "modulus": lift_modulus,
        "coordinates_modulus": lift_coordinates,
        (
            "all_53_residuals_zero_modulus"
            if legacy_profile
            else f"all_{len(equations)}_residuals_zero_modulus"
        ): True,
        "steps": lift_steps,
        "failure": lift_failure,
    },
    "proof_boundary": {
        "proved": (
            f"The displayed GF({prime}) point solves all {len(equations)} normalized integral fibre, "
            "section, and component-jet equations. The exact maximal Jacobian "
            "minor certifies the tangent dimension, and the integer coordinates "
            + (
                "solve every equation modulo the reported power of five."
                if legacy_profile
                else f"solve every equation modulo the reported power of {prime}."
            )
        ),
        "not_proved": (
            "For an overdetermined presentation, a one-dimensional tangent and "
            "finite lift do not alone prove formal smoothness or algebraization. "
            "No Q-rational family or neighbour corridor is constructed."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/certify_k3_04b_a3_a4_a9_marked_gf5_hensel.sage"
        + (
            f" --marked-model-rank {arguments.marked_model_rank}"
            if arguments.marked_model_rank
            else ""
        )
        + (
            f" --marked-section-index {arguments.marked_section_index}"
            if arguments.marked_section_index
            else ""
        )
        + (
            f" --output {relative(output_path)}"
            if output_path != DEFAULT_OUTPUT.resolve()
            else ""
        )
        + (
            f" --fibres {relative(fibres_path)}"
            if fibres_path != DEFAULT_FIBRES.resolve()
            else ""
        )
        + (
            f" --marking {relative(marking_path)}"
            if marking_path != DEFAULT_MARKING.resolve()
            else ""
        )
        + (f" --schema {arguments.schema}" if arguments.schema != SCHEMA else "")
        + (
            f" --lift-precision {arguments.lift_precision}"
            if arguments.lift_precision != 8
            else ""
        )
    ),
}
if arguments.free_parameter_integer is not None:
    payload["finite_precision_lift"]["fixed_free_parameter_integer"] = int(
        arguments.free_parameter_integer
    )
    payload["reproduce"] += (
        f" --free-parameter-integer {arguments.free_parameter_integer}"
    )
serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if output_path.read_text() != serialized:
        raise SystemExit("determinant-500 marked GF(5) Hensel artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)
print(
    "K3MW1HENSEL|"
    f"variables={len(variables)}|equations={len(equations)}|rank={rank}|"
    f"tangent_dimension={len(variables)-rank}|minor={minor}|"
    f"lift_precision={achieved}|status={'PASS' if status.startswith('PASS_ONE') else 'FAILED_GATE'}",
    flush=True,
)
