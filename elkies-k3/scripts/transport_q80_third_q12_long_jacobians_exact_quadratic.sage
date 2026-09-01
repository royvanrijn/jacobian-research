#!/usr/bin/env sage -python
"""Transport six local long Jacobians to the exact QQ(ab) pencil gauge."""

import argparse
import hashlib
import json
import re
from pathlib import Path

from sage.all import GF, PolynomialRing


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
PRIMES = (19, 61, 67, 83, 89, 103, 131)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--alignment", type=Path,
    default=RESULTS / "q80-third-q12-um2-exact-quadratic-pencils-p19-legacy-aligned.json",
)
for prime in PRIMES:
    parser.add_argument(
        f"--p{prime}", type=Path,
        default=RESULTS / f"q80-third-q12-p{prime}-jacobian-interpolated.json",
    )
parser.add_argument(
    "--output", type=Path,
    default=RESULTS / "q80-third-q12-long-jacobians-exact-quadratic-gauge.json",
)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
args.alignment = args.alignment.resolve()
args.output = args.output.resolve()
paths = {prime: getattr(args, f"p{prime}").resolve() for prime in PRIMES}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


alignment = json.loads(args.alignment.read_text())
if alignment.get("status") != "PASS_EXACT_QQ_THIRD_Q12_QUADRATIC_PENCIL_DESCENT_AND_LOCAL_ALIGNMENT":
    raise ValueError("exact quadratic pencil alignment is not certified")
models = {prime: json.loads(path.read_text()) for prime, path in paths.items()}
for prime, model in models.items():
    expected_status = (
        "PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_MOD19_QUADRATIC"
        if prime == 19
        else "PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_COMMON_PRODUCER"
    )
    if model.get("status") != expected_status:
        raise ValueError(f"p={prime}: long Jacobian is not certified")


def coordinate_element(field, generator, coordinates):
    return field(coordinates[0]) + field(coordinates[1]) * generator


transported = {}
common_shapes = None
for prime in PRIMES:
    model = models[prime]
    specialization = model["specialization"]
    if specialization["prime"] != prime:
        raise ValueError(f"p={prime}: specialization mismatch")
    modulus_match = re.fullmatch(
        r"r\^2\s*\+\s*(\d+)\*r\s*\+\s*(\d+)", specialization["extension_modulus"]
    )
    if modulus_match is None:
        raise ValueError(f"p={prime}: cannot parse extension modulus")
    linear, constant = map(int, modulus_match.groups())
    prime_field = GF(prime)
    modulus_ring = PolynomialRing(prime_field, "m")
    m = modulus_ring.gen()
    field = GF(prime**2, "r", modulus=m**2 + linear * m + constant)
    r = field.gen()
    u_ring = PolynomialRing(field, "U")
    U = u_ring.gen()
    u_field = u_ring.fraction_field()

    local_alignment = alignment["local_alignments"][str(prime)]
    scale = prime_field(local_alignment["unique_omega_to_local_generator_scale"])
    matrix = [
        [coordinate_element(field, r, value) for value in row]
        for row in local_alignment["base_PGL2_matrix_a_b_c_d"]
    ]
    a, b = matrix[0]
    c, d = matrix[1]
    if a * d - b * c == 0:
        raise ArithmeticError(f"p={prime}: singular base transformation")
    # Alignment convention: V_exact=(a+b*V_local)/(c+d*V_local).
    V_local = u_field(a - c * U) / u_field(d * U - b)

    inverse_two = prime_field(2) ** -1
    inverse_two_scale = (prime_field(2) * scale) ** -1

    def global_coordinates(value):
        coordinates = list(field(value).list())
        coordinates += [prime_field.zero()] * (2 - len(coordinates))
        local_constant, local_r = coordinates[:2]
        # omega=scale*(2*r+linear), so
        # local_constant+local_r*r = constant + omega_coefficient*omega.
        constant_part = local_constant - local_r * linear * inverse_two
        omega_coefficient = local_r * inverse_two_scale
        reconstructed = field(constant_part) + field(omega_coefficient) * field(scale) * (
            2 * r + linear
        )
        if reconstructed != field(value):
            raise ArithmeticError(f"p={prime}: global omega coordinate round trip failed")
        return [int(constant_part), int(omega_coefficient)]

    def polynomial_from_coordinates(values):
        return u_ring([coordinate_element(field, r, value) for value in values])

    def rational_from_record(record):
        numerator = polynomial_from_coordinates(
            record["numerator_coefficients_low_to_high_1_r"]
        )
        denominator = polynomial_from_coordinates(
            record["denominator_coefficients_low_to_high_1_r"]
        )
        return u_field(numerator(V_local) / denominator(V_local))

    def canonical(value):
        value = u_field(value)
        numerator = u_ring(value.numerator())
        denominator = u_ring(value.denominator())
        leading = denominator.leading_coefficient()
        numerator /= leading
        denominator /= leading
        if denominator.leading_coefficient() != 1:
            raise ArithmeticError(f"p={prime}: rational denominator is not monic")
        return u_field(numerator / denominator)

    def record(value):
        value = canonical(value)
        return {
            "numerator_coefficients_low_to_high_1_omega": [
                global_coordinates(coefficient) for coefficient in value.numerator().list()
            ],
            "denominator_coefficients_low_to_high_1_omega": [
                global_coordinates(coefficient) for coefficient in value.denominator().list()
            ],
            "degrees_numerator_denominator": [
                int(value.numerator().degree()), int(value.denominator().degree())
            ],
        }

    a1, a2, a3, a4, a6 = [
        canonical(rational_from_record(model["weierstrass"][name]))
        for name in ("a1", "a2", "a3", "a4", "a6")
    ]
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3**2 + 4 * a6
    b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
    discriminant = canonical(-b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6)
    c4 = canonical(b2**2 - 24 * b4)
    j = canonical(c4**3 / discriminant)
    expected_discriminant = canonical(rational_from_record(model["discriminant"]))
    expected_j = canonical(rational_from_record(model["j"]))
    if discriminant != expected_discriminant or j != expected_j:
        raise ArithmeticError(f"p={prime}: transported Weierstrass invariants do not replay")

    records = {
        name: record(value)
        for name, value in zip(
            ("a1", "a2", "a3", "a4", "a6", "discriminant", "j"),
            (a1, a2, a3, a4, a6, discriminant, j),
        )
    }
    shapes = {
        name: value["degrees_numerator_denominator"]
        for name, value in records.items()
    }
    if common_shapes is None:
        common_shapes = shapes
    elif shapes != common_shapes:
        raise ArithmeticError(f"p={prime}: transported long-model degree shapes disagree")
    transported[str(prime)] = {
        "specialization": specialization,
        "exact_base_coordinate": "U=V_exact",
        "exact_quadratic_basis": "1,omega with omega^2=16*q1*q2",
        "weierstrass": {name: records[name] for name in ("a1", "a2", "a3", "a4", "a6")},
        "discriminant": records["discriminant"],
        "j": records["j"],
        "literal_weierstrass_invariant_replay": True,
    }


output = {
    "schema": "elkies-k3.q80-third-q12-long-jacobians-exact-quadratic-gauge.v1",
    "status": "PASS_EXACT_TRANSPORTED_THIRD_Q12_LONG_JACOBIANS_COMMON_QUADRATIC_GAUGE",
    "specialization": {"u": "-2", "primes": list(PRIMES)},
    "exact_gauge": {
        "base": "U=V_exact from the exact connected-pencil compiler",
        "coefficient_field": "QQ(omega), omega=4*a*b, omega^2=16*q1*q2",
        "long_weierstrass": "unchanged pinned simple-old-infinity Laurent gauge",
    },
    "common_degree_shapes": common_shapes,
    "transported_models": transported,
    "inputs": {
        "alignment": {
            "path": str(args.alignment.relative_to(ROOT)),
            "sha256": sha256(args.alignment),
        },
        "long_models": {
            str(prime): {
                "path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
            }
            for prime, path in paths.items()
        },
    },
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "seven independently computed generic long Jacobians transported to one exact base coordinate",
            "all coefficients expressed in reductions of the one exact quadratic basis (1,omega)",
            "common rational-function degree shapes and literal discriminant/j replay",
        ],
        "not_proved": [
            "transport of either direction of the birational maps",
            "CRT/LLL reconstruction or a characteristic-zero Jacobian equation",
        ],
    },
    "reproduce": (
        "sage -python "
        "elkies-k3/scripts/transport_q80_third_q12_long_jacobians_exact_quadratic.sage"
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"transported long-Jacobian artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    "Q80THIRDQ12LONGTRANSPORT|field=QQ(omega)|base=V_exact|"
    "primes=19,61,67,83,89,103,131|invariants=Delta,j|"
    "status=PASS_EXACT_TRANSPORTED_THIRD_Q12_LONG_JACOBIANS_COMMON_QUADRATIC_GAUGE"
)
