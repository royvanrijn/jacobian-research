#!/usr/bin/env sage -python
"""Interpolate and reconstruct the exact-gauge long q12 Jacobian."""

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

from sage.all import CRT_list, Matrix, QQ, ZZ


if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--manifest",
    type=Path,
    default=RESULTS / "q80-third-q12-p19-adic-precision1024-sample-manifest.json",
)
parser.add_argument(
    "--transport",
    type=Path,
    default=RESULTS / "q80-third-q12-long-jacobians-exact-quadratic-gauge.json",
)
parser.add_argument(
    "--source",
    type=Path,
    default=RESULTS / "q80-third-q12-exact-pencil-p19-adic-precision1028.json",
)
parser.add_argument(
    "--operands",
    type=Path,
    default=RESULTS / "q80-third-q12-um2-biquadratic-closure-operands-p19-hensel-qq.json",
)
parser.add_argument(
    "--output",
    type=Path,
    default=RESULTS / "q80-third-q12-long-jacobian-p19-adic-reconstructed-qq.json",
)
parser.add_argument("--check", action="store_true")
parser.add_argument(
    "--holdout-prime",
    action="append",
    default=[],
    type=int,
    help="transported prime excluded from CRT/LLL and reserved for literal replay",
)
args = parser.parse_args()
args.manifest = args.manifest.resolve()
args.transport = args.transport.resolve()
args.source = args.source.resolve()
args.operands = args.operands.resolve()
args.output = args.output.resolve()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


manifest = json.loads(args.manifest.read_text())
transport = json.loads(args.transport.read_text())
source = json.loads(args.source.read_text())
operands = json.loads(args.operands.read_text())
if manifest.get("status") != "PASS_EXACT_THIRD_Q12_WEIERSTRASS_P19_ADIC_SAMPLE_BATCH":
    raise ValueError("high-precision sample batch is not certified")
if transport.get("status") != "PASS_EXACT_TRANSPORTED_THIRD_Q12_LONG_JACOBIANS_COMMON_QUADRATIC_GAUGE":
    raise ValueError("finite-prime exact-gauge transport is not certified")
if operands.get("status") != "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_CLOSURE_OPERANDS_P19_HENSEL":
    raise ValueError("exact closure operands are not certified")

prime = 19
digits = int(manifest["specialization"]["digits"])
modulus = prime**digits
omega_square = int(source["quadratic_field"]["omega_square_modulus"]) % modulus
ALL_TRANSPORTED_PRIMES = sorted(int(value) for value in transport["transported_models"])
HELD_OUT_TRANSPORTED_PRIMES = sorted(set(args.holdout_prime))
unknown_holdouts = sorted(set(HELD_OUT_TRANSPORTED_PRIMES) - set(ALL_TRANSPORTED_PRIMES))
if unknown_holdouts:
    raise ValueError(f"holdout primes are absent from the transported models: {unknown_holdouts}")
if prime in HELD_OUT_TRANSPORTED_PRIMES:
    raise ValueError("p=19 supplies the p-adic residue and cannot be a holdout")
RECONSTRUCTION_TRANSPORTED_PRIMES = [
    value for value in ALL_TRANSPORTED_PRIMES if value not in HELD_OUT_TRANSPORTED_PRIMES
]
CRT_AUXILIARY_PRIMES = [
    value for value in RECONSTRUCTION_TRANSPORTED_PRIMES if value != prime
]
reconstruction_modulus = ZZ(modulus)
for local_prime in CRT_AUXILIARY_PRIMES:
    reconstruction_modulus *= local_prime
ZERO = (0, 0)
ONE = (1, 0)


def c(value):
    return int(value[0]) % modulus, int(value[1]) % modulus


def add(left, right):
    return (left[0] + right[0]) % modulus, (left[1] + right[1]) % modulus


def neg(value):
    return -value[0] % modulus, -value[1] % modulus


def sub(left, right):
    return add(left, neg(right))


def mul(left, right):
    return (
        (left[0] * right[0] + omega_square * left[1] * right[1]) % modulus,
        (left[0] * right[1] + left[1] * right[0]) % modulus,
    )


def power(value, exponent):
    result = ONE
    base = value
    while exponent:
        if exponent & 1:
            result = mul(result, base)
        base = mul(base, base)
        exponent >>= 1
    return result


def is_unit(value):
    return (value[0] * value[0] - omega_square * value[1] * value[1]) % prime != 0


def inverse(value):
    norm = (value[0] * value[0] - omega_square * value[1] * value[1]) % modulus
    if norm % prime == 0:
        raise ZeroDivisionError("non-unit quadratic coefficient")
    inverse_norm = pow(norm, -1, modulus)
    return value[0] * inverse_norm % modulus, -value[1] * inverse_norm % modulus


def divide(left, right):
    return mul(left, inverse(right))


samples = []
seen = set()
for record in manifest["samples"]:
    path = ROOT / record["path"]
    if sha256(path) != record["sha256"]:
        raise ArithmeticError(f"sample hash changed: {path}")
    payload = json.loads(path.read_text())
    base = tuple(payload["specialization"]["base_U_coefficients_1_omega"])
    if base in seen:
        raise ArithmeticError(f"duplicate sample base {base}")
    seen.add(base)
    samples.append((path, payload))
samples.sort(key=lambda item: tuple(reversed(item[1]["specialization"]["base_U_coefficients_1_omega"])))
if len(samples) < 20:
    raise ArithmeticError("twenty residue-distinct samples are required")
training = samples[:17]
held_out = samples[17:]


def solve_overdetermined(rows, unknown_count):
    matrix = [[c(value) for value in row] for row in rows]
    rank = 0
    for column in range(unknown_count):
        pivot = next(
            (index for index in range(rank, len(matrix)) if is_unit(matrix[index][column])),
            None,
        )
        if pivot is None:
            raise ArithmeticError(f"interpolation lost a unit pivot at column {column}")
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        scale = inverse(matrix[rank][column])
        matrix[rank] = [mul(value, scale) for value in matrix[rank]]
        for row_index in range(len(matrix)):
            if row_index == rank:
                continue
            factor = matrix[row_index][column]
            if factor == ZERO:
                continue
            matrix[row_index] = [
                sub(matrix[row_index][index], mul(factor, matrix[rank][index]))
                for index in range(unknown_count + 1)
            ]
        rank += 1
    for row in matrix[rank:]:
        if any(value != ZERO for value in row):
            raise ArithmeticError("overdetermined interpolation equations are inconsistent")
    return [matrix[index][-1] for index in range(unknown_count)]


names = ("a1", "a2", "a3", "a4", "a6")
degree_bounds = dict(zip(names, (2, 4, 4, 6, 8)))
interpolated = {}
for coefficient_index, name in enumerate(names):
    degree = degree_bounds[name]
    unknown_count = 2 * degree + 1
    rows = []
    for unused_path, payload in training:
        u_value = c(payload["specialization"]["base_U_coefficients_1_omega"])
        y_value = c(
            payload["weierstrass"]["a1_a2_a3_a4_a6_mod_19_power_1_omega"][
                coefficient_index
            ]
        )
        powers = [power(u_value, exponent) for exponent in range(degree + 1)]
        rows.append(
            powers
            + [neg(mul(y_value, powers[exponent])) for exponent in range(degree)]
            + [mul(y_value, powers[degree])]
        )
    solution = solve_overdetermined(rows, unknown_count)
    numerator = solution[: degree + 1]
    denominator = solution[degree + 1 :] + [ONE]
    for path, payload in held_out:
        u_value = c(payload["specialization"]["base_U_coefficients_1_omega"])
        expected = c(
            payload["weierstrass"]["a1_a2_a3_a4_a6_mod_19_power_1_omega"][
                coefficient_index
            ]
        )
        numerator_value = ZERO
        denominator_value = ZERO
        for exponent in range(degree + 1):
            numerator_value = add(numerator_value, mul(numerator[exponent], power(u_value, exponent)))
            denominator_value = add(
                denominator_value, mul(denominator[exponent], power(u_value, exponent))
            )
        if not is_unit(denominator_value) or divide(numerator_value, denominator_value) != expected:
            raise ArithmeticError(f"held-out replay failed for {name} at {path}")
    interpolated[name] = {
        "degrees_numerator_denominator": [degree, degree],
        "numerator": [list(value) for value in numerator],
        "denominator": [list(value) for value in denominator],
    }


def rational_mod(value, local_prime):
    numerator = int(value.numerator()) % local_prime
    denominator = int(value.denominator()) % local_prime
    if not denominator:
        raise ZeroDivisionError
    return numerator * pow(denominator, -1, local_prime) % local_prime


def validate_candidate_at_finite_primes(candidate, name, selected_primes=None):
    selected_primes = set(selected_primes or ALL_TRANSPORTED_PRIMES)
    for local_prime_text, model in transport["transported_models"].items():
        local_prime = int(local_prime_text)
        if local_prime not in selected_primes:
            continue
        expected = model["weierstrass"][name]
        for exact_key, finite_key in (
            ("numerator", "numerator_coefficients_low_to_high_1_omega"),
            ("denominator", "denominator_coefficients_low_to_high_1_omega"),
        ):
            reduction = [
                [rational_mod(QQ(coordinate), local_prime) for coordinate in pair]
                for pair in candidate[exact_key]
            ]
            if reduction != expected[finite_key]:
                return False
    return True


def reconstruct_projectively(record, name):
    residues = []
    positions = []
    for key in ("numerator", "denominator"):
        for coefficient_index, pair in enumerate(record[key]):
            if key == "denominator" and coefficient_index == len(record[key]) - 1:
                if pair != [1, 0]:
                    raise ArithmeticError(f"{name}: denominator is not normalized")
                continue
            for omega_index, value in enumerate(pair):
                residues.append(ZZ(value) % modulus)
                positions.append((key, coefficient_index, omega_index))
    dimension = len(residues) + 1
    lattice = Matrix(ZZ, dimension, dimension)
    for index in range(dimension - 1):
        lattice[index, index] = modulus
    for index, value in enumerate(residues):
        lattice[dimension - 1, index] = value
    lattice[dimension - 1, dimension - 1] = 1
    reduced = lattice.LLL(delta=0.99)
    diagnostics = []
    for row in sorted(reduced.rows(), key=lambda value: value.dot_product(value)):
        row = list(row)
        if not row[-1] or row[-1] % prime == 0:
            continue
        common = math.gcd(*(abs(int(value)) for value in row))
        if common > 1:
            row = [value // common for value in row]
        if row[-1] < 0:
            row = [-value for value in row]
        candidate = {
            "degrees_numerator_denominator": record["degrees_numerator_denominator"],
            "numerator": [[QQ(0), QQ(0)] for unused in record["numerator"]],
            "denominator": [[QQ(0), QQ(0)] for unused in record["denominator"]],
        }
        candidate["denominator"][-1] = [QQ(1), QQ(0)]
        for index, (key, coefficient_index, omega_index) in enumerate(positions):
            candidate[key][coefficient_index][omega_index] = QQ(row[index]) / QQ(row[-1])
        maximum_bits = max(abs(ZZ(value)).nbits() for value in row)
        diagnostics.append(maximum_bits)
        if validate_candidate_at_finite_primes(candidate, name):
            return candidate, {
                "lattice_dimension": dimension,
                "maximum_primitive_coordinate_bits": maximum_bits,
                "random_lattice_boundary_bits": int(
                    math.ceil(ZZ(modulus).nbits() * (dimension - 1) / dimension)
                ),
                "short_rows_tested": len(diagnostics),
                "validated_primes": sorted(int(value) for value in transport["transported_models"]),
            }
    raise ArithmeticError(
        f"{name}: no projective LLL row validates at all aligned primes; "
        f"short-row bits={diagnostics[:8]}"
    )


def reconstruct_pair_block(modular_pairs, validator, label):
    residues = [
        ZZ(value) % reconstruction_modulus for pair in modular_pairs for value in pair
    ]
    dimension = len(residues) + 1
    lattice = Matrix(ZZ, dimension, dimension)
    for index in range(dimension - 1):
        lattice[index, index] = reconstruction_modulus
    for index, value in enumerate(residues):
        lattice[dimension - 1, index] = value
    lattice[dimension - 1, dimension - 1] = 1
    reduced = lattice.LLL(delta=0.99)
    diagnostics = []
    for row in sorted(reduced.rows(), key=lambda value: value.dot_product(value)):
        row = list(row)
        if not row[-1] or row[-1] % prime == 0:
            continue
        common = math.gcd(*(abs(int(value)) for value in row))
        if common > 1:
            row = [value // common for value in row]
        if row[-1] < 0:
            row = [-value for value in row]
        pairs = [
            [QQ(row[2 * index]) / QQ(row[-1]), QQ(row[2 * index + 1]) / QQ(row[-1])]
            for index in range(len(modular_pairs))
        ]
        maximum_bits = max(abs(ZZ(value)).nbits() for value in row)
        diagnostics.append(maximum_bits)
        if validator(pairs):
            return pairs, {
                "lattice_dimension": dimension,
                "maximum_primitive_coordinate_bits": maximum_bits,
                "random_lattice_boundary_bits": int(
                    math.ceil(reconstruction_modulus.nbits() * (dimension - 1) / dimension)
                ),
                "short_rows_tested": len(diagnostics),
                "validated_primes": sorted(
                    int(value) for value in RECONSTRUCTION_TRANSPORTED_PRIMES
                ),
            }
    raise ArithmeticError(
        f"{label}: no structured projective LLL row validates; "
        f"short-row bits={diagnostics[:8]}"
    )


def reconstruct_scalar_asymmetric(residue, validator, label):
    old_remainder, remainder = ZZ(modulus), ZZ(residue) % modulus
    old_cofactor, cofactor = ZZ(0), ZZ(1)
    diagnostics = []
    while remainder:
        if cofactor and cofactor % prime:
            candidate = QQ(remainder) / QQ(cofactor)
            numerator_bits = abs(ZZ(candidate.numerator())).nbits()
            denominator_bits = ZZ(candidate.denominator()).nbits()
            diagnostics.append([numerator_bits, denominator_bits])
            if validator(candidate):
                return candidate, {
                    "method": "extended-Euclidean asymmetric rational reconstruction",
                    "numerator_bits": numerator_bits,
                    "denominator_bits": denominator_bits,
                    "bit_sum": numerator_bits + denominator_bits,
                    "euclidean_candidates_tested": len(diagnostics),
                    "validated_primes": sorted(
                        int(value) for value in transport["transported_models"]
                    ),
                }
        quotient = old_remainder // remainder
        old_remainder, remainder = remainder, old_remainder - quotient * remainder
        old_cofactor, cofactor = cofactor, old_cofactor - quotient * cofactor
    raise ArithmeticError(
        f"{label}: no asymmetric convergent validates; tail={diagnostics[-8:]}"
    )


def modular_pair_polynomial_multiply(left, right):
    result = [[0, 0] for unused in range(len(left) + len(right) - 1)]
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            product = mul(tuple(left_value), tuple(right_value))
            result[left_index + right_index] = list(
                add(tuple(result[left_index + right_index]), product)
            )
    return result


def modular_pair_polynomial_power(value, exponent):
    result = [[1, 0]]
    base = value
    while exponent:
        if exponent & 1:
            result = modular_pair_polynomial_multiply(result, base)
        base = modular_pair_polynomial_multiply(base, base)
        exponent >>= 1
    return result


denominator_exponents = {"a1": 1, "a2": 2, "a3": 2, "a4": 3, "a6": 4}
modular_H = interpolated["a1"]["denominator"]
for name, exponent in denominator_exponents.items():
    if interpolated[name]["denominator"] != modular_pair_polynomial_power(modular_H, exponent):
        raise ArithmeticError(f"{name}: denominator is not the expected common H power")


def crt_extend_pair(padic_pair, finite_pair_getter):
    moduli = [ZZ(modulus)] + [ZZ(value) for value in CRT_AUXILIARY_PRIMES]
    result = []
    for omega_index, padic_residue in enumerate(padic_pair):
        residues = [ZZ(padic_residue)] + [
            ZZ(finite_pair_getter(local_prime)[omega_index])
            for local_prime in CRT_AUXILIARY_PRIMES
        ]
        result.append(int(CRT_list(residues, moduli)))
    return result


def validate_H(candidate, selected_primes=RECONSTRUCTION_TRANSPORTED_PRIMES):
    selected_primes = set(selected_primes)
    for prime_text, model in transport["transported_models"].items():
        local_prime = int(prime_text)
        if local_prime not in selected_primes:
            continue
        reduction = [
            [rational_mod(value, local_prime) for value in pair] for pair in candidate
        ]
        if reduction != model["weierstrass"]["a1"][
            "denominator_coefficients_low_to_high_1_omega"
        ]:
            return False
    return True


inverse_two = pow(2, -1, modulus)
modular_H_root_constant = list(mul(tuple(modular_H[1]), (inverse_two, 0)))
modular_H_root = [modular_H_root_constant, [1, 0]]
if modular_pair_polynomial_power(modular_H_root, 2) != modular_H:
    raise ArithmeticError("the p-adic common denominator H is not a monic linear square")


def finite_H_root_constant(local_prime):
    h_linear = transport["transported_models"][str(local_prime)]["weierstrass"]["a1"][
        "denominator_coefficients_low_to_high_1_omega"
    ][1]
    inverse_local_two = pow(2, -1, local_prime)
    return [value * inverse_local_two % local_prime for value in h_linear]


def validate_H_root(candidate):
    for local_prime in RECONSTRUCTION_TRANSPORTED_PRIMES:
        try:
            reduction = [rational_mod(value, local_prime) for value in candidate[0]]
        except ZeroDivisionError:
            return False
        if reduction != finite_H_root_constant(local_prime):
            return False
    return True


crt_H_root_constant = crt_extend_pair(modular_H_root_constant, finite_H_root_constant)
reconstructed_H_root, H_root_diagnostic = reconstruct_pair_block(
    [crt_H_root_constant], validate_H_root, "linear square root of common denominator H"
)
exact_H_root = [reconstructed_H_root[0], [QQ(1), QQ(0)]]
exact_H = None
H_reconstruction = {
    "method": "reconstruct monic linear factor and square exactly",
    "identity": "H=(U+r)^2",
    "linear_factor": H_root_diagnostic,
}


def rational_record(record):
    return QQ(ZZ(record["numerator"])) / QQ(ZZ(record["denominator"]))


q1 = rational_record(operands["biquadratic_field"]["q1"])
q2 = rational_record(operands["biquadratic_field"]["q2"])
omega_square_exact = QQ(16) * q1 * q2


def exact_pair_add(left, right):
    return [left[0] + right[0], left[1] + right[1]]


def exact_pair_multiply(left, right):
    return [
        left[0] * right[0] + omega_square_exact * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    ]


def exact_pair_polynomial_multiply(left, right):
    result = [[QQ(0), QQ(0)] for unused in range(len(left) + len(right) - 1)]
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            result[left_index + right_index] = exact_pair_add(
                result[left_index + right_index],
                exact_pair_multiply(left_value, right_value),
            )
    return result


def exact_pair_polynomial_power(value, exponent):
    result = [[QQ(1), QQ(0)]]
    base = value
    while exponent:
        if exponent & 1:
            result = exact_pair_polynomial_multiply(result, base)
        base = exact_pair_polynomial_multiply(base, base)
        exponent >>= 1
    return result


exact_H = exact_pair_polynomial_power(exact_H_root, 2)
if not validate_H(exact_H, RECONSTRUCTION_TRANSPORTED_PRIMES):
    raise ArithmeticError("exact squared common denominator H failed reconstruction-prime replay")
if HELD_OUT_TRANSPORTED_PRIMES and not validate_H(
    exact_H, HELD_OUT_TRANSPORTED_PRIMES
):
    raise ArithmeticError("exact squared common denominator H failed held-out-prime replay")


exact = {}
reconstruction = {"common_denominator_H": H_reconstruction}
for name in names:
    def validate_numerator(candidate, coefficient_name=name):
        for prime_text, model in transport["transported_models"].items():
            local_prime = int(prime_text)
            if local_prime not in RECONSTRUCTION_TRANSPORTED_PRIMES:
                continue
            reduction = [
                [rational_mod(value, local_prime) for value in pair]
                for pair in candidate
            ]
            if reduction != model["weierstrass"][coefficient_name][
                "numerator_coefficients_low_to_high_1_omega"
            ]:
                return False
        return True

    numerator = []
    numerator_reconstruction = {
        "method": "quadratic coefficient pairs reconstructed separately",
        "coefficients": [],
    }
    for coefficient_index, modular_pair in enumerate(interpolated[name]["numerator"]):
        def validate_numerator_coefficient(
            candidate, coefficient_name=name, index=coefficient_index
        ):
            for prime_text, model in transport["transported_models"].items():
                local_prime = int(prime_text)
                if local_prime not in RECONSTRUCTION_TRANSPORTED_PRIMES:
                    continue
                reduction = [rational_mod(value, local_prime) for value in candidate[0]]
                expected = model["weierstrass"][coefficient_name][
                    "numerator_coefficients_low_to_high_1_omega"
                ][index]
                if reduction != expected:
                    return False
            return True

        crt_pair = crt_extend_pair(
            modular_pair,
            lambda local_prime, coefficient_name=name, index=coefficient_index: transport[
                "transported_models"
            ][str(local_prime)]["weierstrass"][coefficient_name][
                "numerator_coefficients_low_to_high_1_omega"
            ][index],
        )
        reconstructed, diagnostic = reconstruct_pair_block(
            [crt_pair],
            validate_numerator_coefficient,
            f"{name} numerator coefficient {coefficient_index}",
        )
        numerator.append(reconstructed[0])
        numerator_reconstruction["coefficients"].append(diagnostic)
    exact[name] = {
        "degrees_numerator_denominator": interpolated[name][
            "degrees_numerator_denominator"
        ],
        "numerator": numerator,
        "denominator": exact_pair_polynomial_power(
            exact_H, denominator_exponents[name]
        ),
    }
    if not validate_candidate_at_finite_primes(
        exact[name], name, RECONSTRUCTION_TRANSPORTED_PRIMES
    ):
        raise ArithmeticError(f"{name}: structured candidate full record failed")
    if HELD_OUT_TRANSPORTED_PRIMES and not validate_candidate_at_finite_primes(
        exact[name], name, HELD_OUT_TRANSPORTED_PRIMES
    ):
        raise ArithmeticError(f"{name}: structured candidate failed held-out-prime replay")
    reconstruction[name] = {
        "denominator_power_of_H": denominator_exponents[name],
        "numerator": numerator_reconstruction,
    }


def encode_exact(candidate):
    return {
        "degrees_numerator_denominator": candidate["degrees_numerator_denominator"],
        "numerator_coefficients_low_to_high_U_1_omega": [
            [str(value) for value in pair] for pair in candidate["numerator"]
        ],
        "denominator_coefficients_low_to_high_U_1_omega": [
            [str(value) for value in pair] for pair in candidate["denominator"]
        ],
    }


output = {
    "schema": "elkies-k3.q80-third-q12-long-jacobian-p19-adic-reconstructed-qq.v1",
    "status": "PASS_CANDIDATE_THIRD_Q12_LONG_JACOBIAN_RECONSTRUCTION_QQ",
    "specialization": {
        "u": "-2",
        "coefficient_field_basis": ["1", "omega"],
        "omega_square": str(omega_square_exact),
    },
    "common_denominator_H_coefficients_low_to_high_U_1_omega": [
        [str(value) for value in pair] for pair in exact_H
    ],
    "weierstrass": {name: encode_exact(exact[name]) for name in names},
    "reconstruction": reconstruction,
    "validation": {
        "p19_digits": digits,
        "training_samples": len(training),
        "held_out_samples": len(held_out),
        "independent_aligned_primes": sorted(
            int(value) for value in transport["transported_models"] if value != "19"
        ),
        "literal_exact_pencil_map_replay": False,
    },
    "inputs": {
        "manifest": {"path": str(args.manifest.relative_to(ROOT)), "sha256": sha256(args.manifest)},
        "transport": {"path": str(args.transport.relative_to(ROOT)), "sha256": sha256(args.transport)},
        "source": {"path": str(args.source.relative_to(ROOT)), "sha256": sha256(args.source)},
        "operands": {
            "path": str(args.operands.relative_to(ROOT)),
            "sha256": sha256(args.operands),
        },
    },
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "unique modular interpolation at the certified long-coefficient degree bounds through 19^digits",
            "three held-out p-adic sample replays",
            "projective LLL candidates reducing to all seven independently transported finite-prime models",
            "the common denominator identities a1:H, a2:H^2, a3:H^2, a4:H^3, a6:H^4",
        ],
        "not_proved": [
            "literal characteristic-zero substitution into the exact genus-one pencil",
            "exact birational maps or a characteristic-zero Jacobian theorem",
            "minimality, discriminant factorization, Kodaira fibres, or transported marking",
        ],
    },
    "reproduce": "sage -python elkies-k3/scripts/reconstruct_q80_third_q12_long_jacobian_p19_adic.sage",
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"reconstructed long-Jacobian artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    f"Q80THIRDQ12LONGRECONSTRUCT|digits={digits}|samples={len(samples)}|"
    "primes=19,61,67,83,89,103,131|"
    "status=PASS_CANDIDATE_THIRD_Q12_LONG_JACOBIAN_RECONSTRUCTION_QQ"
)
