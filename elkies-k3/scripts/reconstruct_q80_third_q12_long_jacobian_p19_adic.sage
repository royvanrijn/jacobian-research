#!/usr/bin/env sage -python
"""Interpolate and reconstruct the exact-gauge long q12 Jacobian."""

import argparse
import hashlib
import json
import math
from pathlib import Path

from sage.all import Matrix, QQ, ZZ


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
    "--output",
    type=Path,
    default=RESULTS / "q80-third-q12-long-jacobian-p19-adic-reconstructed-qq.json",
)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
args.manifest = args.manifest.resolve()
args.transport = args.transport.resolve()
args.source = args.source.resolve()
args.output = args.output.resolve()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


manifest = json.loads(args.manifest.read_text())
transport = json.loads(args.transport.read_text())
source = json.loads(args.source.read_text())
if manifest.get("status") != "PASS_EXACT_THIRD_Q12_WEIERSTRASS_P19_ADIC_SAMPLE_BATCH":
    raise ValueError("high-precision sample batch is not certified")
if transport.get("status") != "PASS_EXACT_TRANSPORTED_THIRD_Q12_LONG_JACOBIANS_COMMON_QUADRATIC_GAUGE":
    raise ValueError("finite-prime exact-gauge transport is not certified")

prime = 19
digits = int(manifest["specialization"]["digits"])
modulus = prime**digits
omega_square = int(source["quadratic_field"]["omega_square_modulus"]) % modulus
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


def validate_candidate_at_finite_primes(candidate, name):
    for local_prime_text, model in transport["transported_models"].items():
        local_prime = int(local_prime_text)
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


exact = {}
reconstruction = {}
for name in names:
    exact[name], reconstruction[name] = reconstruct_projectively(interpolated[name], name)


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
    "specialization": {"u": "-2", "coefficient_field_basis": ["1", "omega"]},
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
