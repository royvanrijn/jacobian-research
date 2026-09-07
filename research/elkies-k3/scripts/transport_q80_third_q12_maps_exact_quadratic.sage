#!/usr/bin/env sage -python
"""Transport six generic q12 birational maps to the exact QQ(ab) gauge."""

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
        default=RESULTS / f"q80-third-q12-p{prime}-birational-maps.json",
    )
parser.add_argument(
    "--output", type=Path,
    default=RESULTS / "q80-third-q12-birational-maps-exact-quadratic-gauge.json",
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
maps_by_prime = {prime: json.loads(path.read_text()) for prime, path in paths.items()}
for prime, maps in maps_by_prime.items():
    expected_status = (
        "PASS_EXACT_GENERIC_THIRD_Q12_BIRATIONAL_MAPS_MOD19_QUADRATIC"
        if prime == 19
        else "PASS_EXACT_GENERIC_THIRD_Q12_BIRATIONAL_MAPS_COMMON_PRODUCER"
    )
    if maps.get("status") != expected_status:
        raise ValueError(f"p={prime}: birational maps are not certified")


def element(field, generator, coordinates):
    return field(coordinates[0]) + field(coordinates[1]) * generator


transported = {}
common_shape = None
for prime in PRIMES:
    maps = maps_by_prime[prime]
    specialization = maps["specialization"]
    match = re.fullmatch(
        r"r\^2\s*\+\s*(\d+)\*r\s*\+\s*(\d+)", specialization["extension_modulus"]
    )
    if match is None or specialization["prime"] != prime:
        raise ValueError(f"p={prime}: map specialization mismatch")
    linear, constant = map(int, match.groups())
    prime_field = GF(prime)
    modulus_ring = PolynomialRing(prime_field, "m")
    m = modulus_ring.gen()
    field = GF(prime**2, "r", modulus=m**2 + linear * m + constant)
    r = field.gen()
    u_ring = PolynomialRing(field, "U")
    U = u_ring.gen()
    u_field = u_ring.fraction_field()
    w_ring = PolynomialRing(u_field, "W")
    W = w_ring.gen()
    w_field = w_ring.fraction_field()

    local_alignment = alignment["local_alignments"][str(prime)]
    scale = prime_field(local_alignment["unique_omega_to_local_generator_scale"])
    matrix = [
        [element(field, r, value) for value in row]
        for row in local_alignment["base_PGL2_matrix_a_b_c_d"]
    ]
    a, b = matrix[0]
    c, d = matrix[1]
    V_local = u_field(a - c * U) / u_field(d * U - b)
    inverse_two = prime_field(2) ** -1
    inverse_two_scale = (prime_field(2) * scale) ** -1

    def global_coordinates(value):
        coordinates = list(field(value).list())
        coordinates += [prime_field.zero()] * (2 - len(coordinates))
        local_constant, local_r = coordinates[:2]
        constant_part = local_constant - local_r * linear * inverse_two
        omega_coefficient = local_r * inverse_two_scale
        if (
            field(constant_part)
            + field(omega_coefficient) * field(scale) * (2 * r + linear)
            != field(value)
        ):
            raise ArithmeticError(f"p={prime}: omega coordinate replay failed")
        return [int(constant_part), int(omega_coefficient)]

    def v_polynomial(values):
        return u_ring([element(field, r, value) for value in values])

    def rational_from_record(record):
        numerator = v_polynomial(record["numerator_coefficients_low_to_high_1_r"])
        denominator = v_polynomial(record["denominator_coefficients_low_to_high_1_r"])
        return u_field(numerator(V_local) / denominator(V_local))

    def canonical_u(value):
        value = u_field(value)
        numerator = u_ring(value.numerator())
        denominator = u_ring(value.denominator())
        leading = denominator.leading_coefficient()
        return u_field((numerator / leading) / (denominator / leading))

    def u_record(value):
        value = canonical_u(value)
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

    def joint_polynomial(values):
        result = w_ring.zero()
        for w_degree, coefficients in enumerate(values):
            polynomial = v_polynomial(coefficients)
            result += u_field(polynomial(V_local)) * W**w_degree
        return result

    def joint_record(source):
        numerator = joint_polynomial(
            source["numerator_coefficients_low_to_high_auxiliary_power_then_V"]
        )
        denominator = joint_polynomial(
            source["denominator_coefficients_low_to_high_auxiliary_power_then_V"]
        )
        value = w_field(numerator / denominator)
        numerator = w_ring(value.numerator())
        denominator = w_ring(value.denominator())
        leading = denominator.leading_coefficient()
        numerator /= leading
        denominator /= leading
        return {
            "numerator_coefficients_low_to_high_W": [
                u_record(coefficient) for coefficient in numerator.list()
            ],
            "denominator_coefficients_low_to_high_W": [
                u_record(coefficient) for coefficient in denominator.list()
            ],
            "degrees_W_numerator_denominator": [
                int(numerator.degree()), int(denominator.degree())
            ],
        }

    forward = {
        coordinate: [joint_record(record) for record in maps["forward_long"][coordinate]]
        for coordinate in ("X", "Y")
    }
    inverse = {}
    for target in ("W", "old_x"):
        source = maps["inverse_long"][target]
        inverse[target] = {
            "weighted_bound": source["weighted_bound"],
            "monomials_X_power_Y_power": source["monomials_X_power_Y_power"],
            "formula": source["formula"],
            "numerator_coefficients": [
                u_record(rational_from_record(record))
                for record in source["numerator_coefficients"]
            ],
            "denominator_coefficients": [
                u_record(rational_from_record(record))
                for record in source["denominator_coefficients"]
            ],
        }

    def shape_u(record):
        return record["degrees_numerator_denominator"]

    shape = {
        "forward": {
            coordinate: [
                {
                    "W": record["degrees_W_numerator_denominator"],
                    "numerator_U": [shape_u(value) for value in record["numerator_coefficients_low_to_high_W"]],
                    "denominator_U": [shape_u(value) for value in record["denominator_coefficients_low_to_high_W"]],
                }
                for record in forward[coordinate]
            ]
            for coordinate in ("X", "Y")
        },
        "inverse": {
            target: {
                "weighted_bound": inverse[target]["weighted_bound"],
                "monomials": inverse[target]["monomials_X_power_Y_power"],
                "numerator_U": [shape_u(value) for value in inverse[target]["numerator_coefficients"]],
                "denominator_U": [shape_u(value) for value in inverse[target]["denominator_coefficients"]],
            }
            for target in ("W", "old_x")
        },
    }
    if common_shape is None:
        common_shape = shape
    elif shape != common_shape:
        raise ArithmeticError(f"p={prime}: transported map shapes disagree")
    transported[str(prime)] = {
        "specialization": specialization,
        "forward_long": forward,
        "inverse_long": inverse,
        "transport_is_field_and_base_substitution_of_certified_generic_identities": True,
    }


output = {
    "schema": "elkies-k3.q80-third-q12-birational-maps-exact-quadratic-gauge.v1",
    "status": "PASS_EXACT_TRANSPORTED_THIRD_Q12_BIRATIONAL_MAPS_COMMON_QUADRATIC_GAUGE",
    "specialization": {"u": "-2", "primes": list(PRIMES)},
    "exact_gauge": {
        "base": "U=V_exact",
        "coefficient_field": "QQ(omega), omega=4*a*b, omega^2=16*q1*q2",
        "long_weierstrass": "pinned simple-old-infinity Laurent gauge",
    },
    "common_shape": common_shape,
    "transported_maps": transported,
    "inputs": {
        "alignment": {
            "path": str(args.alignment.relative_to(ROOT)),
            "sha256": sha256(args.alignment),
        },
        "maps": {
            str(prime): {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
            for prime, path in paths.items()
        },
    },
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "seven certified generic forward and inverse maps transported by exact field/base substitutions",
            "one common exact quadratic coefficient basis and base coordinate",
            "common transported support and degree profiles at all six primes",
        ],
        "not_proved": [
            "CRT/LLL or characteristic-zero reconstruction of the maps",
            "literal characteristic-zero map identities",
        ],
    },
    "reproduce": (
        "sage -python elkies-k3/scripts/transport_q80_third_q12_maps_exact_quadratic.sage"
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"transported map artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    "Q80THIRDQ12MAPTRANSPORT|field=QQ(omega)|base=V_exact|"
    "primes=19,61,67,83,89,103,131|maps=forward,inverse|"
    "status=PASS_EXACT_TRANSPORTED_THIRD_Q12_BIRATIONAL_MAPS_COMMON_QUADRATIC_GAUGE"
)
