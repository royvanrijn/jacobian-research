#!/usr/bin/env sage-python
"""Certify the tangent and finite 7-adic lift gates for G720-S0128.

The pinned nonsquare-twist GF(7) seed has fibre profile ``3I6+6I1`` and a
pole-zero MW2 basis with component depths ``(1,1,3)`` and ``(0,0,0)``.  This
script constructs the normalized integral fibre, two-section, and exact
component-jet equations, computes a maximal Jacobian minor modulo seven, and
applies deterministic simultaneous Newton corrections to every equation.

The finite lift is evidence, not a proof of an infinite Z_7 point or a model
over Q.  The overdetermined ideal still requires a formal-smoothness or
localized-dependence certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FIBRES = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-3a5-source-ansatz-mod7-v1.json"
)
DEFAULT_MARKING = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-3a5-pole0-pairs-mod7-nonsquare-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-3a5-marked-gf7-lift-v1.json"
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def display_path(path: Path) -> str:
    path = path.resolve()
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def padded(values, length):
    if len(values) > length:
        raise ValueError("coefficient vector exceeds declared degree bound")
    return list(values) + [0] * (length - len(values))


NAMES = (
    [f"a{index}" for index in range(9)]
    + [f"b{index}" for index in range(13)]
    + [f"p{index}" for index in range(5)]
    + [f"q{index}" for index in range(7)]
    + [f"r{index}" for index in range(5)]
    + [f"s{index}" for index in range(7)]
)


def build_system(base_ring):
    coefficient_ring = PolynomialRing(
        base_ring, names=NAMES, order="degrevlex"
    )
    variables = list(coefficient_ring.gens())
    cursor = 0

    def take(length):
        nonlocal cursor
        block = variables[cursor : cursor + length]
        cursor += length
        return block

    a = take(9)
    b = take(13)
    p = take(5)
    q = take(7)
    r = take(5)
    s = take(7)
    if cursor != len(variables) or len(variables) != 46:
        raise ArithmeticError("unexpected normalized variable count")

    function_ring = PolynomialRing(coefficient_ring, "t")
    t = function_ring.gen()

    def polynomial(coefficients):
        return sum(
            value * t**index for index, value in enumerate(coefficients)
        )

    A = polynomial(a)
    B = polynomial(b)
    X_P = polynomial(p)
    Y_P = polynomial(q)
    X_Q = polynomial(r)
    Y_Q = polynomial(s)
    discriminant_core = 4 * A**3 + 27 * B**2
    node_P = 2 * A * X_P + 3 * B
    blocks = {
        # The seed is the literal quadratic twist d=3 of A(0)=-3, hence -27.
        "normalization": [a[0] + 27],
        "fibre_at_zero": [discriminant_core[index] for index in range(6)],
        "fibre_at_one": [
            discriminant_core(t + 1)[index] for index in range(6)
        ],
        "fibre_at_infinity": [
            discriminant_core[index] for index in range(19, 25)
        ],
        "component_marking_nonidentity_section": [
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
        ],
        "nonidentity_section": [
            (Y_P**2 - X_P**3 - A * X_P - B)[index]
            for index in range(13)
        ],
        "identity_section": [
            (Y_Q**2 - X_Q**3 - A * X_Q - B)[index]
            for index in range(13)
        ],
    }
    equations = [equation for block in blocks.values() for equation in block]
    if len(equations) != 55:
        raise ArithmeticError("unexpected normalized equation count")
    return variables, blocks, equations


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--fibres", type=Path, default=DEFAULT_FIBRES)
parser.add_argument("--marking", type=Path, default=DEFAULT_MARKING)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--lift-precision", type=int, default=8)
parser.add_argument(
    "--marked-model-rank",
    type=int,
    default=0,
    help="zero-based rank among models carrying at least one marked MW2 pair",
)
parser.add_argument(
    "--marked-pair-index",
    type=int,
    default=0,
    help="zero-based marked-pair index inside the selected model",
)
parser.add_argument(
    "--free-parameter-integer",
    type=int,
    help=(
        "fix the unique nonpivot coordinate s6 to this integer throughout "
        "the lift; it must reduce to the marked value modulo seven"
    ),
)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
if arguments.lift_precision < 1:
    parser.error("lift precision must be positive")

fibres_path = arguments.fibres.resolve()
marking_path = arguments.marking.resolve()
output_path = arguments.output.resolve()
fibres = json.loads(fibres_path.read_text())
marking = json.loads(marking_path.read_text())
if fibres.get("prime") != 7 or not fibres["scan"]["exhausted"]:
    raise ValueError("lift gate requires the exhaustive GF(7) fibre artifact")
if (
    marking.get("schema")
    != "elkies-k3.golay-det720-three-support-pole0-pairs-modp.v1"
    or marking.get("prime") != 7
    or marking.get("quadratic_twist") != 3
    or marking.get("source", {}).get("source_id") != "G720-S0128"
):
    raise ValueError("lift gate requires the nonsquare G720-S0128 marking")
if marking.get("status") != (
    "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_WITH_MARKED_MW2_PAIRS"
):
    raise ValueError("G720-S0128 marking artifact has no certified pair")

marked_models = [row for row in marking["models"] if row["marked_mw2_pairs"]]
if not marked_models or int(marked_models[0]["example_index"]) != 54:
    raise ArithmeticError("the pinned first marked model changed")
if not 0 <= arguments.marked_model_rank < len(marked_models):
    parser.error("--marked-model-rank is outside the marked-model inventory")
marked_model = marked_models[arguments.marked_model_rank]
if not 0 <= arguments.marked_pair_index < len(marked_model["marked_mw2_pairs"]):
    parser.error("--marked-pair-index is outside the selected model")
pair = marked_model["marked_mw2_pairs"][arguments.marked_pair_index]
left_index = int(pair["left_section_index"])
right_index = int(pair["right_section_index"])
section_P = marked_model["basis_section_candidates"][0][left_index]
section_Q = marked_model["basis_section_candidates"][1][right_index]
example_index = int(marked_model["example_index"])
fibre = fibres["examples"][example_index]

field = GF(7)
twist = field(3)
A_values = [int(twist**2 * field(value)) for value in fibre["A_coefficients_low_to_high"]]
B_values = [int(twist**3 * field(value)) for value in fibre["B_coefficients_low_to_high"]]
point_values = (
    padded(A_values, 9)
    + padded(B_values, 13)
    + padded(section_P["X_coefficients_low_to_high"], 5)
    + padded(section_P["Y_coefficients_low_to_high"], 7)
    + padded(section_Q["X_coefficients_low_to_high"], 5)
    + padded(section_Q["Y_coefficients_low_to_high"], 7)
)

variables, equation_blocks, equations = build_system(field)
point = [field(value) for value in point_values]
if len(point) != len(variables):
    raise ArithmeticError("marked point has the wrong number of coordinates")


def evaluate(polynomial_value):
    return field(polynomial_value(*point))


if any(evaluate(equation) for equation in equations):
    raise ArithmeticError("marked point does not solve the lift equations")
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
if free_columns != [45] or NAMES[free_columns[0]] != "s6":
    raise ArithmeticError("unexpected marked-system free coordinate")
if (
    arguments.free_parameter_integer is not None
    and arguments.free_parameter_integer % 7 != int(point[free_columns[0]])
):
    parser.error("--free-parameter-integer must be congruent to 3 modulo seven")
minor_determinant = int(
    jacobian.matrix_from_rows_and_columns(
        pivot_rows, pivot_columns
    ).determinant()
)
block_ranks = {}
block_row_count = 0
for name, block in equation_blocks.items():
    block_row_count += len(block)
    block_ranks[name] = int(jacobian[:block_row_count, :].rank())

unused_integer_variables, unused_integer_blocks, integer_equations = build_system(ZZ)
lift_coordinates = [int(value) for value in point]
lift_steps = []
lift_failure = None
for exponent in range(1, arguments.lift_precision):
    modulus = 7**exponent
    values = [ZZ(equation(*lift_coordinates)) for equation in integer_equations]
    if any(value % modulus for value in values):
        raise ArithmeticError("current coordinates lost their certified precision")
    right = vector(field, [-(value // modulus) for value in values])
    try:
        if arguments.free_parameter_integer is None:
            correction = jacobian.solve_right(right)
        else:
            free_index = free_columns[0]
            free_digit = field(
                (
                    arguments.free_parameter_integer
                    - lift_coordinates[free_index]
                )
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
            correction_entries = [field.zero()] * len(variables)
            correction_entries[free_index] = free_digit
            for index, value in zip(pivot_columns, pivot_correction):
                correction_entries[index] = value
            correction = vector(field, correction_entries)
            if jacobian * correction != right:
                raise ValueError("fixed-parameter correction violates dependent rows")
    except ValueError:
        lift_failure = {
            "failed_lift_from_exponent": exponent,
            "right_hand_side_mod_7": [int(value) for value in right],
        }
        break
    next_modulus = 7 * modulus
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
        {"precision_exponent": exponent + 1, "modulus": next_modulus}
    )

achieved_precision = 1 + len(lift_steps)
lift_modulus = 7**achieved_precision
if any(
    ZZ(equation(*lift_coordinates)) % lift_modulus
    for equation in integer_equations
):
    raise ArithmeticError("reported finite lift has a nonzero residual")

# The two section graphs meet in a squarefree degree-two divisor away from the
# discriminant at the seed, so their intersection multiplicity is stable in a
# sufficiently small marked deformation.
seed_ring = PolynomialRing(field, "t")
t = seed_ring.gen()
A_seed = seed_ring(A_values)
B_seed = seed_ring(B_values)
X_P_seed = seed_ring(section_P["X_coefficients_low_to_high"])
Y_P_seed = seed_ring(section_P["Y_coefficients_low_to_high"])
X_Q_seed = seed_ring(section_Q["X_coefficients_low_to_high"])
Y_Q_seed = seed_ring(section_Q["Y_coefficients_low_to_high"])
intersection_divisor = (X_P_seed - X_Q_seed).gcd(Y_P_seed - Y_Q_seed)
discriminant_seed = 4 * A_seed**3 + 27 * B_seed**2
if (
    intersection_divisor.degree() != 2
    or intersection_divisor.gcd(intersection_divisor.derivative()).degree() != 0
    or intersection_divisor.gcd(discriminant_seed).degree() != 0
):
    raise ArithmeticError("pinned section intersection is not etale and smooth")

status = (
    "PASS_ONE_DIMENSIONAL_MARKED_TANGENT_AND_EXPLICIT_Z7_LIFT_TO_REQUESTED_PRECISION"
    if rank == len(variables) - 1
    and minor_determinant
    and lift_failure is None
    else (
        "PASS_ONE_DIMENSIONAL_MARKED_TANGENT_BUT_FINITE_LIFT_STOPPED"
        if rank == len(variables) - 1 and minor_determinant
        else "PASS_POINT_BUT_EXPECTED_MARKED_TANGENT_DIMENSION_FAILED"
    )
)
output = {
    "schema": "elkies-k3.golay-det720-3a5-marked-gf7-lift.v1",
    "status": status,
    "prime": 7,
    "inputs": {
        relative(fibres_path): digest(fibres_path),
        relative(marking_path): digest(marking_path),
    },
    "seed": {
        "fibre_example_index": example_index,
        "left_section_index": left_index,
        "right_section_index": right_index,
        "quadratic_twist": 3,
        "coordinate_names": NAMES,
        "coordinates_mod_7": [int(value) for value in point],
        "section_intersection_divisor_coefficients_low_to_high": [
            int(value) for value in intersection_divisor
        ],
        "section_intersection_etale_and_away_from_singular_fibres": True,
    },
    "system": {
        "normalization": (
            "literal d=3 twist of A(0)=-3, so A(0)=-27; supports at 0,1,infinity"
        ),
        "variable_count": len(variables),
        "equation_count": len(equations),
        "equation_block_sizes": {
            name: len(block) for name, block in equation_blocks.items()
        },
        "equation_block_cumulative_jacobian_ranks": block_ranks,
        "all_equation_residuals_zero_mod_7": True,
    },
    "jacobian_certificate": {
        "rank_mod_7": rank,
        "tangent_dimension": len(variables) - rank,
        "pivot_row_indices": pivot_rows,
        "pivot_column_indices": pivot_columns,
        "pivot_variable_names": [NAMES[index] for index in pivot_columns],
        "omitted_free_variable_names": [
            name for index, name in enumerate(NAMES) if index not in pivot_columns
        ],
        "pivot_minor_determinant_mod_7": minor_determinant,
    },
    "finite_precision_lift": {
        "requested_precision_exponent": arguments.lift_precision,
        "achieved_precision_exponent": achieved_precision,
        "modulus": lift_modulus,
        "coordinates_modulus": lift_coordinates,
        "all_55_residuals_zero_modulus": True,
        "steps": lift_steps,
        "failure": lift_failure,
    },
    "proof_boundary": {
        "proved": (
            "The pinned marked GF(7) point solves all normalized integral fibre, "
            "section, and component-jet equations. The maximal unit Jacobian "
            "minor certifies the reported tangent dimension, all 55 equations "
            "vanish on the displayed lift modulo the reported power of seven, "
            "and the seed section intersections are etale and avoid singular fibres."
        ),
        "not_proved": (
            "Finite precision is not an infinite compatible Z_7 lift. No formal "
            "smoothness of the overdetermined ideal, Q-rational family, rational "
            "parameterization, target multisection spectrum, neighbour corridor, "
            "or specialization rank jump is proved."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/certify_golay_det720_3a5_marked_gf7_lift.sage"
    ),
}
if arguments.free_parameter_integer is not None:
    output["finite_precision_lift"]["fixed_free_parameter_integer"] = int(
        arguments.free_parameter_integer
    )
    output["reproduce"] += (
        f" --free-parameter-integer {arguments.free_parameter_integer} "
        f"--lift-precision {arguments.lift_precision} "
        f"--output {display_path(output_path)}"
    )
if arguments.marked_model_rank or arguments.marked_pair_index:
    output["seed"]["marked_model_rank"] = int(arguments.marked_model_rank)
    output["seed"]["marked_pair_index"] = int(arguments.marked_pair_index)
    output["reproduce"] += (
        f" --marked-model-rank {arguments.marked_model_rank}"
        f" --marked-pair-index {arguments.marked_pair_index}"
    )
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if output_path.read_text() != serialized:
        raise SystemExit("Golay-720 3A5 lift artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "GOLAY7203A5LIFT|"
    f"variables={len(variables)}|equations={len(equations)}|rank={rank}|"
    f"tangent_dimension={len(variables)-rank}|minor={minor_determinant}|"
    f"lift_precision={achieved_precision}|"
    f"status={'PASS' if status.startswith('PASS_ONE') else 'FAILED_GATE'}",
    flush=True,
)
