#!/usr/bin/env sage-python
"""Hensel-lift modular norm-four three-node normalization candidates.

For each degree-eight survivor in a modular genus-3 normalization screen, use
the monic cubic square factor to seed the six exact equations

    q(m0,m1,m2;u) = (u^3+s2*u^2+s1*u+s0)^2 * Q2(u).

The coefficients of ``Q2`` are eliminated from degrees 8, 7, and 6.  A seed
whose 6-by-6 Jacobian is nonsingular has a unique p-adic lift.  This script
lifts it digit by digit, attempts rational reconstruction at declared
checkpoints, and accepts only exact QQ solutions whose normalized cover
section satisfies both Weierstrass coefficient identities.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
import json
from itertools import product
from pathlib import Path
import sys

import numpy as np

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
DEFAULT_MODEL = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
DEFAULT_SMOOTH = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisection-collisions-full-v1.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-genus3-normalization-hensel-v1.json"
CHORD_SCRIPT = SCRIPTS / "construct_elkies_2026_bisections.sage"
HASH_SCRIPT = SCRIPTS / "hash_bisection_extensions.py"


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
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def coefficients(polynomial):
    return [rational_text(polynomial[index]) for index in range(polynomial.degree() + 1)]


def rational_function(record, ring, field):
    numerator = ring([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = ring([QQ(value) for value in record["denominator_coefficients_low_to_high"]])
    return field(numerator) / field(denominator)


def primitive_integer_polynomial(polynomial, integer_ring):
    denominator = ZZ(polynomial.denominator())
    result = integer_ring(denominator * polynomial)
    content = ZZ(result.content())
    if content:
        result //= content
    return result


def evaluate_integer_polynomial(polynomial, values):
    return ZZ(polynomial(*values))


def squareclass_decomposition(polynomial, ring):
    factorization = polynomial.factor()
    square_part = ring.one()
    reduced = ring(factorization.unit())
    for factor, exponent in factorization:
        square_part *= factor ** (int(exponent) // 2)
        if int(exponent) % 2:
            reduced *= factor
    if square_part**2 * reduced != polynomial:
        raise ArithmeticError("squareclass decomposition failed")
    return square_part, reduced


def singular_hensel_tree(equations, jacobian_modp, seed, prime, target_exponent, state_cap):
    """Enumerate the complete singular Hensel tree while it remains under cap."""

    field = jacobian_modp.base_ring()
    rank = jacobian_modp.rank()
    kernel_basis = list(jacobian_modp.right_kernel().basis())
    left_kernel_basis = list(jacobian_modp.left_kernel().basis())
    direction_count = int(prime ** len(kernel_basis))
    pivot_columns = jacobian_modp.pivots()
    pivot_rows = jacobian_modp.transpose().pivots()
    pivot_matrix_inverse = jacobian_modp.matrix_from_rows_and_columns(
        pivot_rows, pivot_columns
    ).inverse()
    states = [vector(ZZ, seed)]
    levels = []
    truncated = False
    exact_rational_solutions = []
    for exponent in range(1, target_exponent):
        power = prime**exponent
        next_modulus = power * prime
        compatible = []
        for state in states:
            values = vector(
                ZZ,
                [
                    evaluate_integer_polynomial(equation, state) % next_modulus
                    for equation in equations
                ],
            )
            if any(value % power for value in values):
                raise ArithmeticError("singular Hensel tree invariant failed")
            rhs = vector(
                field,
                [field(-(value // power)) for value in values],
            )
            if any(left_vector.dot_product(rhs) for left_vector in left_kernel_basis):
                continue
            compatible.append((state, rhs))
        output_count = len(compatible) * direction_count
        level = {
            "input_exponent": exponent,
            "input_state_count": len(states),
            "compatible_parent_count": len(compatible),
            "output_exponent": exponent + 1,
            "output_state_count": output_count,
        }
        if len(compatible) <= 1000:
            level["compatible_parent_state_samples"] = [
                list(map(int, state)) for state, unused_rhs in compatible
            ]
        for state, unused_rhs in compatible:
            reconstructed = []
            try:
                reconstructed = [
                    QQ(ZZ(value).rational_reconstruction(ZZ(power)))
                    for value in state
                ]
            except (ArithmeticError, ValueError):
                continue
            if all(equation(*reconstructed) == 0 for equation in equations):
                serialized = list(map(rational_text, reconstructed))
                if serialized not in exact_rational_solutions:
                    exact_rational_solutions.append(serialized)
        if output_count > state_cap:
            level["truncated_at_state_cap"] = True
            levels.append(level)
            truncated = True
            states = []
            break
        level["truncated_at_state_cap"] = False
        next_states = []
        for state, rhs in compatible:
            selected_rhs = vector(field, [rhs[index] for index in pivot_rows])
            pivot_solution = pivot_matrix_inverse * selected_rhs
            particular = vector(field, 6)
            for index, pivot_column in enumerate(pivot_columns):
                particular[pivot_column] = pivot_solution[index]
            for parameters in product(range(int(prime)), repeat=len(kernel_basis)):
                correction = particular + sum(
                    (
                        field(parameter) * basis_vector
                        for parameter, basis_vector in zip(parameters, kernel_basis)
                        if parameter
                    ),
                    vector(field, 6),
                )
                next_states.append(
                    vector(
                        ZZ,
                        [
                            (state[index] + power * ZZ(correction[index]))
                            % next_modulus
                            for index in range(6)
                        ],
                    )
                )
        states = next_states
        levels.append(level)
        if not states:
            break
    return {
        "target_exponent": target_exponent,
        "kernel_dimension": 6 - int(rank),
        "directions_per_compatible_parent": direction_count,
        "state_cap": state_cap,
        "levels": levels,
        "truncated": truncated,
        "final_state_count": len(states),
        "final_state_samples": [list(map(int, state)) for state in states[:16]],
        "exact_rational_solutions": exact_rational_solutions,
    }


def resolved_singular_hensel_branches(
    equations,
    jacobian_integer,
    jacobian_modp,
    seed,
    prime,
    target_exponent,
    state_cap,
):
    """Track only singular states that extend one further p-adic digit.

    For exponent at least two, Taylor expansion modulo p^(k+2) is linear in
    the digit correction because 2*k >= k+2.  Thus each parent's polynomial
    values and Jacobian modulo p^2 are evaluated once; all kernel directions
    are tested using small integer matrix arithmetic.
    """

    field = jacobian_modp.base_ring()
    rank = jacobian_modp.rank()
    kernel_basis = list(jacobian_modp.right_kernel().basis())
    left_kernel_basis = list(jacobian_modp.left_kernel().basis())
    direction_parameters = list(product(range(int(prime)), repeat=len(kernel_basis)))
    direction_count = len(direction_parameters)
    kernel_direction_matrix = np.array(
        [
            [
                int(value)
                for value in sum(
                    (
                        field(parameter) * basis_vector
                        for parameter, basis_vector in zip(parameters, kernel_basis)
                        if parameter
                    ),
                    vector(field, 6),
                )
            ]
            for parameters in direction_parameters
        ],
        dtype=np.int64,
    )
    left_kernel_matrix = np.array(
        [[int(value) for value in basis_vector] for basis_vector in left_kernel_basis],
        dtype=np.int64,
    )
    pivot_columns = jacobian_modp.pivots()
    pivot_rows = jacobian_modp.transpose().pivots()
    pivot_inverse = jacobian_modp.matrix_from_rows_and_columns(
        pivot_rows, pivot_columns
    ).inverse()

    def rhs_for_state(state, exponent):
        power = prime**exponent
        values = [evaluate_integer_polynomial(equation, state) for equation in equations]
        if any(value % power for value in values):
            raise ArithmeticError("resolved singular branch invariant failed")
        return vector(field, [field(-(value // power)) for value in values])

    def is_compatible(rhs):
        return not any(
            left_vector.dot_product(rhs) for left_vector in left_kernel_basis
        )

    def particular_solution(rhs):
        selected_rhs = vector(field, [rhs[index] for index in pivot_rows])
        pivot_solution = pivot_inverse * selected_rhs
        result = vector(field, 6)
        for index, pivot_column in enumerate(pivot_columns):
            result[pivot_column] = pivot_solution[index]
        return result

    states = [vector(ZZ, seed)]
    levels = []
    exact_rational_solutions = []
    truncated = False
    for exponent in range(1, target_exponent):
        power = prime**exponent
        next_modulus = power * prime
        child_test_count = len(states) * direction_count
        if child_test_count > state_cap:
            levels.append(
                {
                    "input_exponent": exponent,
                    "extendable_parent_count": len(states),
                    "child_test_count": child_test_count,
                    "truncated_at_state_cap": True,
                }
            )
            truncated = True
            break
        next_states = []
        for state in states:
            rhs = rhs_for_state(state, exponent)
            if not is_compatible(rhs):
                raise ArithmeticError("retained singular state is no longer compatible")
            particular = particular_solution(rhs)
            if exponent >= 2:
                parent_values = vector(
                    ZZ,
                    [evaluate_integer_polynomial(equation, state) for equation in equations],
                )
                quotient_mod_p2 = vector(
                    ZZ,
                    [(value // power) % (prime**2) for value in parent_values],
                )
                jacobian_mod_p2 = matrix(
                    ZZ,
                    [
                        [
                            evaluate_integer_polynomial(entry, state) % (prime**2)
                            for entry in row
                        ]
                        for row in jacobian_integer.rows()
                    ],
                )
            particular_integer = np.array(
                [int(value) for value in particular], dtype=np.int64
            )
            corrections = (kernel_direction_matrix + particular_integer) % int(prime)
            if exponent == 1:
                selected_corrections = []
                for correction_values in corrections:
                    correction_integer = vector(
                        ZZ, list(map(int, correction_values))
                    )
                    child = vector(
                        ZZ,
                        [
                            (state[index] + power * correction_integer[index])
                            % next_modulus
                            for index in range(6)
                        ],
                    )
                    if is_compatible(rhs_for_state(child, exponent + 1)):
                        selected_corrections.append(correction_values)
            else:
                quotient_array = np.array(
                    list(map(int, quotient_mod_p2)), dtype=np.int64
                )
                jacobian_array = np.array(
                    [[int(value) for value in row] for row in jacobian_mod_p2.rows()],
                    dtype=np.int64,
                )
                next_numerators = quotient_array + corrections @ jacobian_array.transpose()
                if np.any(next_numerators % int(prime)):
                    raise ArithmeticError("Taylor digit correction invariant failed")
                child_quotients = (next_numerators // int(prime)) % int(prime)
                compatibility = np.all(
                    (child_quotients @ left_kernel_matrix.transpose()) % int(prime)
                    == 0,
                    axis=1,
                )
                selected_corrections = corrections[compatibility]
            for correction_values in selected_corrections:
                correction_integer = vector(ZZ, list(map(int, correction_values)))
                next_states.append(
                    vector(
                        ZZ,
                        [
                            (state[index] + power * correction_integer[index])
                            % next_modulus
                            for index in range(6)
                        ],
                    )
                )
        levels.append(
            {
                "input_exponent": exponent,
                "extendable_parent_count": len(states),
                "child_test_count": child_test_count,
                "output_exponent": exponent + 1,
                "extendable_child_count": len(next_states),
                "truncated_at_state_cap": False,
            }
        )
        states = next_states
        for state in states:
            try:
                reconstructed = [
                    QQ(ZZ(value).rational_reconstruction(ZZ(next_modulus)))
                    for value in state
                ]
            except (ArithmeticError, ValueError):
                continue
            if all(equation(*reconstructed) == 0 for equation in equations):
                serialized = list(map(rational_text, reconstructed))
                if serialized not in exact_rational_solutions:
                    exact_rational_solutions.append(serialized)
        if not states:
            break
    return {
        "target_exponent": target_exponent,
        "kernel_dimension": 6 - int(rank),
        "directions_per_parent": direction_count,
        "state_cap": state_cap,
        "levels": levels,
        "truncated": truncated,
        "final_extendable_state_count": len(states),
        "final_extendable_state_samples": [
            list(map(int, state)) for state in states[:32]
        ],
        "exact_rational_solutions": exact_rational_solutions,
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--screen", type=Path, required=True)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument("--smooth-collisions", type=Path, default=DEFAULT_SMOOTH)
parser.add_argument("--precision", type=int, default=256)
parser.add_argument(
    "--singular-second-order",
    action="store_true",
    help="enumerate every first-order kernel direction and test lifting to p^3",
)
parser.add_argument(
    "--singular-tree-exponent",
    type=int,
    default=0,
    help="complete target p-adic exponent for selected singular seeds",
)
parser.add_argument("--only-trace-index", type=int)
parser.add_argument("--only-survivor-index", type=int)
parser.add_argument(
    "--singular-resolved-exponent",
    type=int,
    default=0,
    help="target exponent for one-step-lookahead singular branch tracking",
)
parser.add_argument(
    "--singular-state-cap",
    type=int,
    default=1000000,
    help="maximum p^(6-rank) directions per singular seed",
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if args.precision < 2:
    parser.error("--precision must be at least two")

chord = load_script("r17_direct_genus3_hensel_chord", CHORD_SCRIPT)
hasher = load_script("r17_direct_genus3_hensel_hasher", HASH_SCRIPT)
screen_path = args.screen.resolve()
screen = json.loads(screen_path.read_text())
if screen.get("schema") != "elkies-k3.r17-norm12-direct-genus3-normalization-modp-search.v1":
    raise ValueError("unexpected modular screen schema")
prime = ZZ(screen["prime"])
finite_field = GF(prime)

model_path = args.model.resolve()
model = json.loads(model_path.read_text())
if model.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
    raise ValueError("expected certified direct model")
smooth_path = args.smooth_collisions.resolve()
smooth = json.loads(smooth_path.read_text())
if smooth.get("status") != "PASS_EXTENSION_CANONICALIZATION" or not smooth.get("compact_output"):
    raise ValueError("expected compact smooth collision artifact")
smooth_by_digest = {item["extension_sha256"]: item for item in smooth["extension_manifest"]}

u_ring = PolynomialRing(QQ, "u")
u = u_ring.gen()
function_field = u_ring.fraction_field()
weierstrass = model["weierstrass_model"]
A = u_ring([QQ(value) for value in weierstrass["A_coefficients_low_to_high"]])
B = u_ring([QQ(value) for value in weierstrass["B_coefficients_low_to_high"]])
surface_discriminant = u_ring(-16 * (4 * A**3 + 27 * B**2))
curve = EllipticCurve(function_field, [A, B])
basis = [curve(rational_function(record["X"], u_ring, function_field), rational_function(record["Y"], u_ring, function_field)) for record in model["sections"]["records"]]

variable_names = ("m0", "m1", "m2", "s0", "s1", "s2")
coefficient_ring = PolynomialRing(QQ, names=variable_names)
m0, m1, m2, s0, s1, s2 = coefficient_ring.gens()
integer_ring = PolynomialRing(ZZ, names=variable_names)
symbolic_u_ring = PolynomialRing(coefficient_ring, "u")
symbolic_u = symbolic_u_ring.gen()

survivors_by_trace = {}
for survivor in screen["survivors"]:
    if len(survivor["branch_coefficients_low_to_high"]) != 9:
        continue
    if len(survivor["removed_square_factor_coefficients_low_to_high"]) != 4:
        continue
    survivors_by_trace.setdefault(int(survivor["trace_index"]), []).append(survivor)

lift_records = []
exact_candidates = []
rank_histogram = Counter()
checkpoint_exponents = sorted(set([8, 16, 32, 64, 128, args.precision]))
checkpoint_exponents = [value for value in checkpoint_exponents if value <= args.precision]
for trace_index, survivors in sorted(survivors_by_trace.items()):
    if args.only_trace_index is not None and trace_index != args.only_trace_index:
        continue
    basis_coordinates = survivors[0]["basis_coordinates"]
    trace = sum((coefficient * point for coefficient, point in zip(basis_coordinates, basis) if coefficient), curve(0))
    frame = chord.trace_chord_frame(trace[0], trace[1], u_ring)
    h, Nx, Ny = (frame[key] for key in ("h", "Nx", "Ny"))
    if h.degree() != 0:
        raise ArithmeticError("screened trace is not polynomial")

    M = m0 + m1 * symbolic_u + m2 * symbolic_u**2
    q_symbolic = M**4 - 6 * M**2 * symbolic_u_ring(Nx) - 8 * M * symbolic_u_ring(Ny) - 3 * symbolic_u_ring(Nx)**2 - 4 * symbolic_u_ring(A)
    q_coefficients = [coefficient_ring(q_symbolic[index]) for index in range(9)]
    c2 = q_coefficients[8]
    c1 = q_coefficients[7] - 2 * s2 * c2
    c0 = q_coefficients[6] - 2 * s2 * c1 - (s2**2 + 2 * s1) * c2
    square_factor = symbolic_u**3 + s2 * symbolic_u**2 + s1 * symbolic_u + s0
    residual_quadratic = c2 * symbolic_u**2 + c1 * symbolic_u + c0
    residual = q_symbolic - square_factor**2 * residual_quadratic
    equations_QQ = [coefficient_ring(residual[index]) for index in range(6)]
    equations = [primitive_integer_polynomial(equation, integer_ring) for equation in equations_QQ]
    jacobian = matrix(integer_ring, [[equation.derivative(variable) for variable in integer_ring.gens()] for equation in equations])

    for survivor_index, survivor in enumerate(survivors):
        if (
            args.only_survivor_index is not None
            and survivor_index != args.only_survivor_index
        ):
            continue
        square_coefficients = survivor["removed_square_factor_coefficients_low_to_high"]
        seed = vector(ZZ, survivor["m0_m1_m2"] + square_coefficients[:3])
        if any(evaluate_integer_polynomial(equation, seed) % prime for equation in equations):
            raise ArithmeticError("modular survivor does not satisfy exact reduced system")
        jacobian_modp = matrix(finite_field, [[finite_field(evaluate_integer_polynomial(entry, seed)) for entry in row] for row in jacobian.rows()])
        rank = int(jacobian_modp.rank())
        rank_histogram[str(rank)] += 1
        lift_record = {
            "trace_index": trace_index,
            "survivor_index_within_trace": survivor_index,
            "seed_variables_m0_m1_m2_s0_s1_s2": list(map(int, seed)),
            "jacobian_rank": rank,
            "unique_hensel_lift": rank == 6,
            "rational_reconstruction_checkpoints": [],
        }
        seed_values = vector(
            ZZ,
            [evaluate_integer_polynomial(equation, seed) for equation in equations],
        )
        first_rhs = vector(
            finite_field,
            [finite_field(-(value // prime)) for value in seed_values],
        )
        augmented_rank = int(jacobian_modp.augment(first_rhs).rank())
        lift_record["first_order_augmented_rank"] = augmented_rank
        lift_record["lifts_mod_prime_squared"] = augmented_rank == rank
        lift_record["first_order_lift_dimension"] = (
            6 - rank if augmented_rank == rank else None
        )
        if rank != 6:
            lift_record["status"] = (
                "SINGULAR_MODP_SEED_COMPATIBLE_TO_PRIME_SQUARED"
                if augmented_rank == rank
                else "PASS_FIRST_ORDER_OBSTRUCTION_NO_PRIME_SQUARED_LIFT"
            )
            if args.singular_second_order and augmented_rank == rank:
                kernel_basis = list(jacobian_modp.right_kernel().basis())
                direction_count = int(prime ** len(kernel_basis))
                lift_record["first_order_kernel_direction_count"] = direction_count
                if direction_count > args.singular_state_cap:
                    lift_record["second_order_status"] = "SKIPPED_STATE_CAP"
                else:
                    reduced_augmented = jacobian_modp.augment(first_rhs).rref()
                    pivots = jacobian_modp.pivots()
                    particular = vector(finite_field, 6)
                    for row_index, pivot_column in enumerate(pivots):
                        particular[pivot_column] = reduced_augmented[row_index, 6]
                    compatible_count = 0
                    compatible_samples = []
                    for parameters in product(range(int(prime)), repeat=len(kernel_basis)):
                        lift_digit = particular + sum(
                            (
                                finite_field(parameter) * basis_vector
                                for parameter, basis_vector in zip(parameters, kernel_basis)
                                if parameter
                            ),
                            vector(finite_field, 6),
                        )
                        state = vector(
                            ZZ,
                            [
                                int(seed[index]) + int(prime) * ZZ(lift_digit[index])
                                for index in range(6)
                            ],
                        )
                        values_mod_prime_cubed = vector(
                            ZZ,
                            [
                                evaluate_integer_polynomial(equation, state)
                                % (prime**3)
                                for equation in equations
                            ],
                        )
                        if any(value % (prime**2) for value in values_mod_prime_cubed):
                            raise ArithmeticError("first-order lift enumeration invariant failed")
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
                        if len(compatible_samples) < 8:
                            compatible_samples.append(list(map(int, state)))
                    lift_record["second_order_compatible_direction_count"] = compatible_count
                    lift_record["second_order_compatible_state_samples_mod_p2"] = compatible_samples
                    lift_record["second_order_status"] = (
                        "PASS_COMPATIBLE_LIFTS_TO_PRIME_CUBED"
                        if compatible_count
                        else "PASS_SECOND_ORDER_OBSTRUCTION_NO_PRIME_CUBED_LIFT"
                    )
                    if not compatible_count:
                        lift_record["status"] = "PASS_SECOND_ORDER_OBSTRUCTION_NO_PRIME_CUBED_LIFT"
            if args.singular_tree_exponent > 1 and augmented_rank == rank:
                lift_record["singular_hensel_tree"] = singular_hensel_tree(
                    equations,
                    jacobian_modp,
                    seed,
                    prime,
                    args.singular_tree_exponent,
                    args.singular_state_cap,
                )
                tree = lift_record["singular_hensel_tree"]
                if tree["levels"] and tree["levels"][-1]["compatible_parent_count"] == 0:
                    lift_record["status"] = (
                        "PASS_SINGULAR_HENSEL_TREE_NO_LIFT_TO_TARGET_EXPONENT"
                    )
            if args.singular_resolved_exponent > 1 and augmented_rank == rank:
                lift_record["resolved_singular_hensel_branches"] = (
                    resolved_singular_hensel_branches(
                        equations,
                        jacobian,
                        jacobian_modp,
                        seed,
                        prime,
                        args.singular_resolved_exponent,
                        args.singular_state_cap,
                    )
                )
            lift_records.append(lift_record)
            continue
        inverse_jacobian = jacobian_modp.inverse()
        current = vector(ZZ, seed)
        exact_parameters = None
        for exponent in range(1, args.precision):
            power = prime**exponent
            next_modulus = power * prime
            values = vector(ZZ, [evaluate_integer_polynomial(equation, current) for equation in equations])
            if any(value % power for value in values):
                raise ArithmeticError("Hensel invariant failed")
            rhs = vector(finite_field, [finite_field(-(value // power)) for value in values])
            correction = inverse_jacobian * rhs
            current = vector(ZZ, [(current[index] + power * ZZ(correction[index])) % next_modulus for index in range(6)])
            achieved_exponent = exponent + 1
            if achieved_exponent not in checkpoint_exponents:
                continue
            reconstructed = []
            failed = False
            for value in current:
                try:
                    reconstructed.append(QQ(ZZ(value).rational_reconstruction(ZZ(next_modulus))))
                except (ArithmeticError, ValueError):
                    failed = True
                    break
            exact_zero = False
            if not failed:
                exact_zero = all(equation(*reconstructed) == 0 for equation in equations_QQ)
            lift_record["rational_reconstruction_checkpoints"].append({
                "exponent": achieved_exponent,
                "all_coordinates_reconstructed": not failed,
                "exact_system_zero": bool(exact_zero),
            })
            if exact_zero:
                exact_parameters = reconstructed
                break
        lift_record["final_exponent"] = achieved_exponent
        if exact_parameters is None:
            lift_record["status"] = "NO_EXACT_RATIONAL_RECONSTRUCTION_WITHIN_PRECISION"
            lift_records.append(lift_record)
            continue
        lift_record["status"] = "PASS_EXACT_RATIONAL_FACTORIZATION_PARAMETERS"
        lift_record["exact_variables_m0_m1_m2_s0_s1_s2"] = list(map(rational_text, exact_parameters))
        lift_records.append(lift_record)

        m0_value, m1_value, m2_value = exact_parameters[:3]
        slope = u_ring(m0_value + m1_value * u + m2_value * u**2)
        q = u_ring(slope**4 - 6 * slope**2 * Nx - 8 * slope * Ny - 3 * Nx**2 - 4 * A)
        square_part, q_reduced = squareclass_decomposition(q, u_ring)
        if q_reduced.degree() != 2:
            raise ArithmeticError("exact reconstructed branch is not quadratic")
        sum_x = u_ring(slope**2 - Nx)
        product_x = function_field(((slope * Nx + Ny) ** 2 - B) / Nx)
        if product_x.denominator() != 1:
            raise ArithmeticError("residual x-product is not polynomial")
        product_x = u_ring(product_x)
        if sum_x**2 - 4 * product_x != q:
            raise ArithmeticError("residual discriminant identity failed")
        x0 = u_ring(sum_x / 2)
        x1 = u_ring(square_part / 2)
        intercept = u_ring(-(Ny + slope * Nx))
        y0 = u_ring(slope * x0 + intercept)
        y1 = u_ring(slope * square_part / 2)
        if y0**2 + y1**2 * q_reduced != x0**3 + 3 * x0 * x1**2 * q_reduced + A * x0 + B:
            raise ArithmeticError("constant cover-section identity failed")
        if 2 * y0 * y1 != 3 * x0**2 * x1 + x1**3 * q_reduced + A * x1:
            raise ArithmeticError("linear cover-section identity failed")
        branch = {"numerator_coefficients": coefficients(q_reduced), "denominator_coefficients": ["1"]}
        sympy_u = hasher.sp.Symbol("u")
        extension = hasher.extension_key(branch, sympy_u)
        extension_digest = hasher.key_digest(extension)
        exact_candidates.append({
            "label": f"genus3-hensel-trace-{trace_index:05d}-seed-{survivor_index:03d}",
            "trace_index": trace_index,
            "trace_basis_coordinates": basis_coordinates,
            "trace_translation_orbit_mask": int(survivor["translation_orbit_mask"]),
            "slope_coefficients_m0_m1_m2": list(map(rational_text, exact_parameters[:3])),
            "raw_branch_coefficients_low_to_high": coefficients(q),
            "removed_square_factor_coefficients_low_to_high": coefficients(square_part),
            "branch": branch,
            "extension_squareclass": extension,
            "extension_sha256": extension_digest,
            "smooth_atlas_match": smooth_by_digest.get(extension_digest),
            "branch_coprime_to_surface_discriminant": q_reduced.gcd(surface_discriminant).degree() == 0,
            "lifted_section": {
                "field": "QQ(u,r), r^2=q_reduced(u)",
                "x0_coefficients_low_to_high": coefficients(x0),
                "x1_coefficients_low_to_high": coefficients(x1),
                "y0_coefficients_low_to_high": coefficients(y0),
                "y1_coefficients_low_to_high": coefficients(y1),
                "constant_and_linear_identities_verified": True,
            },
        })

collisions = {}
for candidate in exact_candidates:
    collisions.setdefault(candidate["extension_sha256"], []).append(candidate["label"])
collisions = {key: labels for key, labels in collisions.items() if len(labels) >= 2}
output_path = args.output if args.output.is_absolute() else ROOT / args.output
output_path.parent.mkdir(parents=True, exist_ok=True)
payload = {
    "schema": "elkies-k3.r17-norm12-direct-genus3-normalization-hensel.v1",
    "status": "PASS_EXACT_RATIONAL_NORMALIZATION_CANDIDATES" if exact_candidates else "PASS_BOUNDED_HENSEL_NO_EXACT_RATIONAL_NORMALIZATION",
    "prime": int(prime),
    "precision": args.precision,
    "degree_eight_seed_count": sum(len(items) for items in survivors_by_trace.values()),
    "jacobian_rank_histogram": dict(sorted(rank_histogram.items(), key=lambda item: int(item[0]))),
    "unique_hensel_seed_count": sum(record["unique_hensel_lift"] for record in lift_records),
    "singular_second_order_requested": args.singular_second_order,
    "singular_tree_exponent": args.singular_tree_exponent,
    "singular_resolved_exponent": args.singular_resolved_exponent,
    "singular_second_order_no_lift_count": sum(
        record.get("second_order_compatible_direction_count") == 0
        for record in lift_records
    ),
    "singular_second_order_survivor_count": sum(
        record.get("second_order_compatible_direction_count", 0) > 0
        for record in lift_records
    ),
    "candidate_count": len(exact_candidates),
    "smooth_atlas_match_count": sum(candidate["smooth_atlas_match"] is not None for candidate in exact_candidates),
    "candidate_collision_count": len(collisions),
    "candidate_collisions": collisions,
    "candidates": exact_candidates,
    "lifts": lift_records,
    "proof_boundary": (
        "Every displayed full-rank modular seed is lifted uniquely through the "
        "displayed p-adic exponent and tested for exact rational reconstruction. "
        "Accepted candidates satisfy the factorization and cover-section identities "
        "over QQ. Singular modular seeds and rational points beyond the reconstruction "
        "precision remain unresolved by this artifact."
    ),
    "inputs": {relative(path): digest(path) for path in (screen_path, model_path, smooth_path, CHORD_SCRIPT, HASH_SCRIPT)},
    "reproducing_command": (
        "sage -python elkies-k3/scripts/lift_r17_norm12_direct_genus3_normalizations.sage "
        f"--screen {relative(screen_path)} --precision {args.precision}"
        + (" --singular-second-order" if args.singular_second_order else "")
        + (
            ""
            if not args.singular_tree_exponent
            else f" --singular-tree-exponent {args.singular_tree_exponent}"
        )
        + (
            ""
            if args.only_trace_index is None
            else f" --only-trace-index {args.only_trace_index}"
        )
        + (
            ""
            if args.only_survivor_index is None
            else f" --only-survivor-index {args.only_survivor_index}"
        )
        + (
            ""
            if not args.singular_resolved_exponent
            else f" --singular-resolved-exponent {args.singular_resolved_exponent}"
        )
        + f" --singular-state-cap {args.singular_state_cap}"
    ),
}
output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"R17GENUS3HENSEL|p={prime}|seeds={payload['degree_eight_seed_count']}"
    f"|full_rank={payload['unique_hensel_seed_count']}|precision={args.precision}"
    f"|candidates={len(exact_candidates)}|smooth_matches={payload['smooth_atlas_match_count']}"
    f"|collisions={len(collisions)}|output={relative(output_path)}|status={payload['status']}",
    flush=True,
)
