#!/usr/bin/env sage-python
"""Certify the Jacobian/Hensel gate for the marked NS0031 GF(7) seed.

The normalized short-Weierstrass model has coefficient variables for ``A,B``,
the pole-zero section ``P=(X,Y)``, and the pole-one section
``R=(N/C^2,M/C^3)`` with monic linear ``C``.  Fibre divisibility at
``0,1,infinity``, ``A(0)=-3``, the two section equations, and four node-jet
equations encode the exact component marking.  The Jacobian modulo seven
tests whether the full marked system has the expected one-dimensional tangent
space and records an explicit maximal unit minor.

Tangent dimension one is a necessary Hensel precursor, but for this
overdetermined presentation it is not by itself a formal-lifting theorem.  The
script does not algebraize or descend a locus to Q, parameterize it
rationally, or construct a neighbour corridor.
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
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-source-ansatz-mod7-v1.json"
)
DEFAULT_MARKING = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0031-a1-2a7-marking-mod7-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0031-marked-gf7-hensel-v1.json"
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


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
    + ["c0"]
    + [f"n{index}" for index in range(7)]
    + [f"m{index}" for index in range(10)]
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
    c0 = take(1)[0]
    n = take(7)
    m = take(10)
    if cursor != len(variables) or len(variables) != 52:
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
    C = t + c0
    N_R = polynomial(n)
    M_R = polynomial(m)
    discriminant_core = 4 * A**3 + 27 * B**2
    blocks = {
        "normalization": [a[0] + 3],
        "fibre_at_zero": [discriminant_core[index] for index in range(2)],
        "fibre_at_one": [
            discriminant_core(t + 1)[index] for index in range(8)
        ],
        "fibre_at_infinity": [
            discriminant_core[index] for index in range(17, 25)
        ],
        "component_marking": [
            2 * a[0] * p[0] + 3 * b[0],
            q[0],
            2 * a[8] * p[4] + 3 * b[12],
            2 * (a[8] * p[3] + a[7] * p[4]) + 3 * b[11],
            q[6],
            q[5],
            2 * A(1) * N_R(1) + 3 * B(1) * C(1) ** 2,
            M_R(1),
        ],
        "pole_zero_section": [
            (Y_P**2 - X_P**3 - A * X_P - B)[index]
            for index in range(13)
        ],
        "pole_one_section": [
            (M_R**2 - N_R**3 - A * N_R * C**4 - B * C**6)[index]
            for index in range(19)
        ],
    }
    equations = [equation for block in blocks.values() for equation in block]
    if len(equations) != 59:
        raise ArithmeticError("unexpected normalized equation count")
    return variables, blocks, equations


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--fibres", type=Path, default=DEFAULT_FIBRES)
parser.add_argument("--marking", type=Path, default=DEFAULT_MARKING)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--lift-precision", type=int, default=8)
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
    raise ValueError("Hensel gate requires the exhaustive GF(7) fibre artifact")
if (
    marking.get("schema")
    != "elkies-k3.lattice-foundry-ns0031-a1-2a7-marking-modp.v1"
    or marking.get("prime") != 7
    or marking.get("quadratic_twist_square_class") != "square"
):
    raise ValueError("Hensel gate requires the square-twist NS0031 GF(7) marking")
if marking.get("status") != (
    "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_WITH_MARKED_MW2_PAIRS"
):
    raise ValueError("NS0031 marking artifact has no certified pair")

marked_models = [
    row
    for row in marking["models"]
    if any(section["marked_pairs"] for section in row["pole_one_sections"])
]
if len(marked_models) != 1:
    raise ArithmeticError("expected one normalized model carrying marked pairs")
marked_model = marked_models[0]
example_index = int(marked_model["example_index"])
if example_index != 157:
    raise ArithmeticError("the pinned marked model index changed")
pole_one_index, pole_one = next(
    (index, section)
    for index, section in enumerate(marked_model["pole_one_sections"])
    if section["marked_pairs"]
)
if len(pole_one["marked_pairs"]) != 1:
    raise ArithmeticError("expected one marked pole-zero partner for the chosen sign")
pole_zero_index = int(pole_one["marked_pairs"][0]["pole_zero_index"])
pole_zero = marked_model["pole_zero_sections"][pole_zero_index]
fibre = fibres["examples"][example_index]

field = GF(7)
names = NAMES
variables, equation_blocks, equations = build_system(field)

point_values = (
    padded(fibre["A_coefficients_low_to_high"], 9)
    + padded(fibre["B_coefficients_low_to_high"], 13)
    + padded(pole_zero["X_coefficients_low_to_high"], 5)
    + padded(pole_zero["Y_coefficients_low_to_high"], 7)
    + [pole_one["C_coefficients_low_to_high"][0]]
    + padded(pole_one["X_numerator_coefficients_low_to_high"], 7)
    + padded(pole_one["Y_numerator_coefficients_low_to_high"], 10)
)
point = [field(value) for value in point_values]
if len(point) != len(variables):
    raise ArithmeticError("marked point has the wrong number of coordinates")


def evaluate(polynomial_value):
    return field(polynomial_value(*point))


residuals = [evaluate(equation) for equation in equations]
if any(residuals):
    raise ArithmeticError("the marked GF(7) point does not solve the lift system")

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
minor_determinant = (
    int(
        jacobian.matrix_from_rows_and_columns(
            pivot_rows, pivot_columns
        ).determinant()
    )
    if len(pivot_columns) == len(pivot_rows) == rank
    else 0
)
block_ranks = {}
block_row_count = 0
for name, block in equation_blocks.items():
    block_row_count += len(block)
    block_ranks[name] = int(jacobian[:block_row_count, :].rank())

integer_variables, unused_integer_blocks, integer_equations = build_system(ZZ)
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
        correction = jacobian.solve_right(right)
    except ValueError:
        lift_failure = {
            "failed_lift_from_exponent": exponent,
            "right_hand_side_mod_7": [int(value) for value in right],
        }
        break
    next_modulus = modulus * 7
    lift_coordinates = [
        int((value + modulus * int(delta)) % next_modulus)
        for value, delta in zip(lift_coordinates, correction)
    ]
    next_values = [
        ZZ(equation(*lift_coordinates)) for equation in integer_equations
    ]
    if any(value % next_modulus for value in next_values):
        raise ArithmeticError("linear correction failed the next Hensel digit")
    lift_steps.append(
        {
            "precision_exponent": exponent + 1,
            "modulus": next_modulus,
            "maximum_centered_correction_digit": max(
                min(int(value), 7 - int(value)) for value in correction
            ),
        }
    )

achieved_lift_precision = 1 + len(lift_steps)
lift_modulus = 7**achieved_lift_precision
lift_residuals = [
    int(ZZ(equation(*lift_coordinates)) % lift_modulus)
    for equation in integer_equations
]
if any(lift_residuals):
    raise ArithmeticError("reported finite-precision lift has nonzero residuals")

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
    "schema": "elkies-k3.lattice-foundry-ns0031-marked-gf7-hensel.v1",
    "status": status,
    "prime": 7,
    "inputs": {
        relative(fibres_path): digest(fibres_path),
        relative(marking_path): digest(marking_path),
    },
    "seed": {
        "fibre_example_index": example_index,
        "pole_zero_section_index": pole_zero_index,
        "pole_one_section_index": pole_one_index,
        "coordinate_names": names,
        "coordinates_mod_7": [int(value) for value in point],
    },
    "system": {
        "normalization": "A(0)=-3; supports fixed at 0,1,infinity; C monic",
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
        "pivot_variable_names": [names[index] for index in pivot_columns],
        "omitted_free_variable_names": [
            name for index, name in enumerate(names) if index not in pivot_columns
        ],
        "pivot_minor_determinant_mod_7": minor_determinant,
    },
    "finite_precision_lift": {
        "requested_precision_exponent": arguments.lift_precision,
        "achieved_precision_exponent": achieved_lift_precision,
        "modulus": lift_modulus,
        "coordinates_modulus": lift_coordinates,
        "all_59_residuals_zero_modulus": True,
        "steps": lift_steps,
        "failure": lift_failure,
    },
    "proof_boundary": {
        "proved": (
            "The displayed GF(7) marked point solves all normalized integral "
            "fibre, section, and component-jet equations. The displayed maximal "
            "unit Jacobian minor certifies the reported tangent-space dimension "
            "for this exact overdetermined presentation. The displayed integer "
            "coordinates solve all 59 equations modulo the reported power of seven."
        ),
        "not_proved": (
            "For an overdetermined presentation, tangent dimension one alone does "
            "not prove an infinite compatible Z_7 lift of every equation, and a "
            "finite-precision lift is not such a proof. No Q-rational "
            "or algebraic characteristic-zero point, rational parameterization, "
            "physical neighbour corridor, or specialization rank jump is constructed."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/certify_lattice_foundry_ns0031_marked_gf7_hensel.sage"
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if output_path.read_text() != serialized:
        raise SystemExit("NS0031 GF(7) Hensel certificate is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "FOUNDRYNS0031HENSEL|"
    f"variables={len(variables)}|equations={len(equations)}|rank={rank}|"
    f"tangent_dimension={len(variables)-rank}|minor={minor_determinant}|"
    f"lift_precision={achieved_lift_precision}|"
    f"status={'PASS' if status.startswith('PASS_ONE') else 'FAILED_GATE'}",
    flush=True,
)
