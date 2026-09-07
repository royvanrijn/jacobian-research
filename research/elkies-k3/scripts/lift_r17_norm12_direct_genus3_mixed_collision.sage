#!/usr/bin/env sage-python
"""Hensel-lift one mixed degree-six/degree-eight norm-four cover collision."""

from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
from itertools import product
import json
from pathlib import Path
import sys

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
DEFAULT_MODEL = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)
CHORD_SCRIPT = SCRIPTS / "construct_elkies_2026_bisections.sage"


def load_script(name, path):
    loader = SourceFileLoader(name, str(path))
    spec = spec_from_loader(name, loader)
    if spec is None:
        raise ImportError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    loader.exec_module(module)
    return module


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def relative(path):
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def rational_text(value):
    value = QQ(value)
    return (
        str(value.numerator())
        if value.denominator() == 1
        else f"{value.numerator()}/{value.denominator()}"
    )


def rational_function(record, ring, field):
    numerator = ring(
        [QQ(value) for value in record["numerator_coefficients_low_to_high"]]
    )
    denominator = ring(
        [QQ(value) for value in record["denominator_coefficients_low_to_high"]]
    )
    return field(numerator) / field(denominator)


def primitive_integer_polynomial(polynomial, integer_ring):
    denominator = ZZ(polynomial.denominator())
    result = integer_ring(denominator * polynomial)
    content = ZZ(result.content())
    if content:
        result //= content
    return result


def evaluate(polynomial, values):
    return ZZ(polynomial(*values))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--screen", type=Path, required=True)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument("--degree-six-trace-index", type=int, required=True)
parser.add_argument("--degree-six-survivor-index", type=int, required=True)
parser.add_argument("--degree-eight-trace-index", type=int, required=True)
parser.add_argument("--degree-eight-survivor-index", type=int, required=True)
parser.add_argument("--precision", type=int, default=256)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()

chord = load_script("r17_direct_genus3_mixed_collision_chord", CHORD_SCRIPT)
screen_path = args.screen.resolve()
screen = json.loads(screen_path.read_text())
if screen.get("schema") != (
    "elkies-k3.r17-norm12-direct-genus3-normalization-modp-search.v1"
):
    raise ValueError("unexpected modular screen schema")
prime = ZZ(screen["prime"])
finite_field = GF(prime)


def selected_survivor(trace_index, survivor_index):
    items = [
        survivor
        for survivor in screen["survivors"]
        if int(survivor["trace_index"]) == trace_index
    ]
    if survivor_index < 0 or survivor_index >= len(items):
        raise ValueError("survivor index outside selected trace")
    return items[survivor_index]


degree_six_seed = selected_survivor(
    args.degree_six_trace_index, args.degree_six_survivor_index
)
degree_eight_seed = selected_survivor(
    args.degree_eight_trace_index, args.degree_eight_survivor_index
)
if len(degree_six_seed["branch_coefficients_low_to_high"]) != 7:
    raise ValueError("first seed is not degree six")
if len(degree_six_seed["removed_square_factor_coefficients_low_to_high"]) != 3:
    raise ValueError("first seed does not have a quadratic square factor")
if len(degree_eight_seed["branch_coefficients_low_to_high"]) != 9:
    raise ValueError("second seed is not degree eight")
if len(degree_eight_seed["removed_square_factor_coefficients_low_to_high"]) != 4:
    raise ValueError("second seed does not have a cubic square factor")

model_path = args.model.resolve()
model = json.loads(model_path.read_text())
if model.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
    raise ValueError("expected a certified direct model")
u_ring = PolynomialRing(QQ, "u")
function_field = u_ring.fraction_field()
weierstrass = model["weierstrass_model"]
A = u_ring([QQ(value) for value in weierstrass["A_coefficients_low_to_high"]])
B = u_ring([QQ(value) for value in weierstrass["B_coefficients_low_to_high"]])
curve = EllipticCurve(function_field, [A, B])
basis = [
    curve(
        rational_function(record["X"], u_ring, function_field),
        rational_function(record["Y"], u_ring, function_field),
    )
    for record in model["sections"]["records"]
]

names = (
    "am0", "am1", "am2", "as0", "as1",
    "bm0", "bm1", "bm2", "bs0", "bs1", "bs2",
)
coefficient_ring = PolynomialRing(QQ, names=names)
variables = coefficient_ring.gens()
integer_ring = PolynomialRing(ZZ, names=names)
symbolic_u_ring = PolynomialRing(coefficient_ring, "u")
u = symbolic_u_ring.gen()


def raw_branch(seed, offset):
    trace = sum(
        (
            coefficient * point
            for coefficient, point in zip(seed["basis_coordinates"], basis)
            if coefficient
        ),
        curve(0),
    )
    frame = chord.trace_chord_frame(trace[0], trace[1], u_ring)
    h, Nx, Ny = (frame[key] for key in ("h", "Nx", "Ny"))
    if h.degree() != 0:
        raise ArithmeticError("selected trace is not in the polynomial chart")
    m0, m1, m2 = variables[offset : offset + 3]
    slope = m0 + m1 * u + m2 * u**2
    q = (
        slope**4
        - 6 * slope**2 * symbolic_u_ring(Nx)
        - 8 * slope * symbolic_u_ring(Ny)
        - 3 * symbolic_u_ring(Nx) ** 2
        - 4 * symbolic_u_ring(A)
    )
    return q


left_q = raw_branch(degree_six_seed, 0)
left_s0, left_s1 = variables[3:5]
left_coefficients = [coefficient_ring(left_q[index]) for index in range(9)]
left_c2 = left_coefficients[6]
left_c1 = left_coefficients[5] - 2 * left_s1 * left_c2
left_c0 = (
    left_coefficients[4]
    - 2 * left_s1 * left_c1
    - (left_s1**2 + 2 * left_s0) * left_c2
)
left_square = u**2 + left_s1 * u + left_s0
left_residual = left_q - left_square**2 * (
    left_c2 * u**2 + left_c1 * u + left_c0
)
equations_QQ = [left_coefficients[8], left_coefficients[7]]
equations_QQ.extend(coefficient_ring(left_residual[index]) for index in range(4))

right_q = raw_branch(degree_eight_seed, 5)
right_s0, right_s1, right_s2 = variables[8:11]
right_coefficients = [coefficient_ring(right_q[index]) for index in range(9)]
right_c2 = right_coefficients[8]
right_c1 = right_coefficients[7] - 2 * right_s2 * right_c2
right_c0 = (
    right_coefficients[6]
    - 2 * right_s2 * right_c1
    - (right_s2**2 + 2 * right_s1) * right_c2
)
right_square = u**3 + right_s2 * u**2 + right_s1 * u + right_s0
right_residual = right_q - right_square**2 * (
    right_c2 * u**2 + right_c1 * u + right_c0
)
equations_QQ.extend(coefficient_ring(right_residual[index]) for index in range(6))
equations_QQ.extend(
    (
        left_c1 * right_c2 - right_c1 * left_c2,
        left_c0 * right_c2 - right_c0 * left_c2,
    )
)
equations = [
    primitive_integer_polynomial(equation, integer_ring)
    for equation in equations_QQ
]
jacobian = matrix(
    integer_ring,
    [
        [equation.derivative(variable) for variable in integer_ring.gens()]
        for equation in equations
    ],
)
seed = vector(
    ZZ,
    degree_six_seed["m0_m1_m2"]
    + degree_six_seed["removed_square_factor_coefficients_low_to_high"][:2]
    + degree_eight_seed["m0_m1_m2"]
    + degree_eight_seed["removed_square_factor_coefficients_low_to_high"][:3],
)
if any(evaluate(equation, seed) % prime for equation in equations):
    raise ArithmeticError("mixed collision seed does not solve the reduced system")
jacobian_modp = matrix(
    finite_field,
    [
        [finite_field(evaluate(entry, seed)) for entry in row]
        for row in jacobian.rows()
    ],
)
rank = int(jacobian_modp.rank())
current = vector(ZZ, seed)
checkpoints = []
exact_solution = None
status = "SINGULAR_MODP_COLLISION_SEED"
compatible_count = None
compatible_samples = []
if rank == len(variables):
    for exponent in range(1, args.precision):
        power = prime**exponent
        next_modulus = power * prime
        values = vector(ZZ, [evaluate(equation, current) for equation in equations])
        if any(value % power for value in values):
            raise ArithmeticError("Hensel invariant failed")
        rhs = vector(
            finite_field,
            [finite_field(-(value // power)) for value in values],
        )
        if jacobian_modp.augment(rhs).rank() != rank:
            status = "PASS_HENSEL_OBSTRUCTION_NO_LIFT"
            break
        correction = jacobian_modp.solve_right(rhs)
        current = vector(
            ZZ,
            [
                (current[index] + power * ZZ(correction[index])) % next_modulus
                for index in range(len(variables))
            ],
        )
        achieved = exponent + 1
        if achieved not in {8, 16, 32, 64, 128, args.precision}:
            continue
        reconstructed = []
        try:
            reconstructed = [
                QQ(value.rational_reconstruction(next_modulus)) for value in current
            ]
        except (ArithmeticError, ValueError):
            pass
        exact_zero = bool(reconstructed) and all(
            equation(*reconstructed) == 0 for equation in equations_QQ
        )
        checkpoints.append(
            {
                "exponent": achieved,
                "all_coordinates_reconstructed": bool(reconstructed),
                "exact_system_zero": exact_zero,
            }
        )
        if exact_zero:
            exact_solution = reconstructed
            status = "PASS_EXACT_RATIONAL_COMMON_COVER_ATOM"
            break
    else:
        status = "NO_EXACT_RATIONAL_RECONSTRUCTION_WITHIN_PRECISION"
else:
    values = vector(ZZ, [evaluate(equation, seed) for equation in equations])
    rhs = vector(
        finite_field,
        [finite_field(-(value // prime)) for value in values],
    )
    if jacobian_modp.augment(rhs).rank() != rank:
        status = "PASS_FIRST_ORDER_OBSTRUCTION_NO_PRIME_SQUARED_LIFT"
    else:
        kernel_basis = list(jacobian_modp.right_kernel().basis())
        direction_count = int(prime ** len(kernel_basis))
        if direction_count <= 100000:
            particular = jacobian_modp.solve_right(rhs)
            compatible_count = 0
            compatible_samples = []
            for parameters in product(
                range(int(prime)), repeat=len(kernel_basis)
            ):
                correction = particular + sum(
                    (
                        finite_field(parameter) * basis_vector
                        for parameter, basis_vector in zip(parameters, kernel_basis)
                        if parameter
                    ),
                    vector(finite_field, len(variables)),
                )
                state = vector(
                    ZZ,
                    [
                        int(seed[index]) + int(prime) * ZZ(correction[index])
                        for index in range(len(variables))
                    ],
                )
                values_mod_prime_cubed = vector(
                    ZZ,
                    [evaluate(equation, state) % (prime**3) for equation in equations],
                )
                if any(value % (prime**2) for value in values_mod_prime_cubed):
                    raise ArithmeticError("first-order lift invariant failed")
                second_rhs = vector(
                    finite_field,
                    [
                        finite_field(-(value // (prime**2)))
                        for value in values_mod_prime_cubed
                    ],
                )
                if jacobian_modp.augment(second_rhs).rank() != rank:
                    continue
                compatible_count += 1
                if len(compatible_samples) < 16:
                    compatible_samples.append(list(map(int, state)))
            if compatible_count == 0:
                status = "PASS_SECOND_ORDER_OBSTRUCTION_NO_PRIME_CUBED_LIFT"
        else:
            compatible_count = None
            compatible_samples = []

output_path = args.output if args.output.is_absolute() else ROOT / args.output
output_path.parent.mkdir(parents=True, exist_ok=True)
payload = {
    "schema": "elkies-k3.r17-norm12-direct-genus3-mixed-collision-hensel.v1",
    "status": status,
    "prime": int(prime),
    "precision": args.precision,
    "degree_six_trace_index": args.degree_six_trace_index,
    "degree_six_survivor_index": args.degree_six_survivor_index,
    "degree_eight_trace_index": args.degree_eight_trace_index,
    "degree_eight_survivor_index": args.degree_eight_survivor_index,
    "translation_orbit_masks": [
        int(degree_six_seed["translation_orbit_mask"]),
        int(degree_eight_seed["translation_orbit_mask"]),
    ],
    "seed_variables": list(map(int, seed)),
    "equation_count": len(equations),
    "variable_count": len(variables),
    "jacobian_rank": rank,
    "singular_kernel_dimension": len(variables) - rank,
    "second_order_compatible_direction_count": (
        None if rank == len(variables) else compatible_count
    ),
    "second_order_compatible_state_samples_mod_p2": (
        [] if rank == len(variables) else compatible_samples
    ),
    "rational_reconstruction_checkpoints": checkpoints,
    "exact_variables": (
        None if exact_solution is None else list(map(rational_text, exact_solution))
    ),
    "proof_boundary": (
        "The exact mixed degree-six/degree-eight factorization equations and the "
        "two monic quadratic atom-equality equations are lifted from the displayed "
        "modular seed. A Hensel obstruction excludes this seed from any p-adic, "
        "hence rational, collision. A reconstructed point still requires scalar "
        "squareclass and cover-section verification. Singular compatible seeds "
        "remain unresolved."
    ),
    "inputs": {
        relative(screen_path): digest(screen_path),
        relative(model_path): digest(model_path),
        relative(CHORD_SCRIPT): digest(CHORD_SCRIPT),
    },
}
output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "R17GENUS3MIXEDCOLLISION"
    f"|p={prime}|rank={rank}/{len(variables)}"
    f"|output={relative(output_path)}|status={status}",
    flush=True,
)
