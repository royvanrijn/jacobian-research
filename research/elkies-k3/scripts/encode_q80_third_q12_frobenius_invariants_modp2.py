#!/usr/bin/env python3
"""Encode a complete common-prime q12 child in generator-free invariants."""

import argparse
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
for name in ("resolved", "long", "minimal", "maps", "marking", "horizontal", "output"):
    parser.add_argument(f"--{name}", type=Path, required=True)
args = parser.parse_args()
for name in ("resolved", "long", "minimal", "maps", "marking", "horizontal", "output"):
    setattr(args, name, getattr(args, name).resolve())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text())


resolved, long_model, minimal, maps, marking, horizontal = (
    load(getattr(args, name))
    for name in ("resolved", "long", "minimal", "maps", "marking", "horizontal")
)
expected = (
    (resolved, "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_COMMON_PRODUCER"),
    (long_model, "PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_COMMON_PRODUCER"),
    (minimal, "PASS_EXACT_MINIMAL_THIRD_Q12_JACOBIAN_AND_FIBRES_COMMON_PRODUCER"),
    (maps, "PASS_EXACT_GENERIC_THIRD_Q12_BIRATIONAL_MAPS_COMMON_PRODUCER"),
    (marking, "PASS_EXACT_TRANSPORTED_THIRD_Q12_COMPONENT_MARKING_COMMON_PRODUCER"),
    (horizontal, "PASS_EXACT_THIRD_Q12_HORIZONTAL_FROM_COMMON_CLOSURE_PRODUCER"),
)
if any(payload.get("status") != status for payload, status in expected):
    raise ValueError("one or more invariant-encoding inputs are uncertified")
specialization = resolved["specialization"]
if any(payload["specialization"] != specialization for payload in (long_model, minimal, maps, marking)):
    raise ValueError("complete-child specializations disagree")
prime = int(specialization["prime"])
modulus = specialization["extension_modulus"]
match = re.fullmatch(r"r\^2 \+ (\d+)\*r \+ (\d+)", modulus)
if match is None:
    raise ValueError("cannot parse quadratic extension modulus")
linear, constant = map(int, match.groups())
field_discriminant = (linear * linear - 4 * constant) % prime
if pow(field_discriminant, (prime - 1) // 2, prime) != prime - 1:
    raise ArithmeticError("quadratic extension discriminant is not a nonsquare")
inverse_two = pow(2, -1, prime)
inverse_four = pow(4, -1, prime)
encoded_count = 0


def encode(coordinates):
    global encoded_count
    a, b = (int(coordinates[0]) % prime, int(coordinates[1]) % prime)
    trace = (2 * a - linear * b) % prime
    coefficient_discriminant = field_discriminant * b * b % prime
    norm = (trace * trace - coefficient_discriminant) * inverse_four % prime
    reconstructed = ((trace + linear * b) * inverse_two % prime, b)
    if reconstructed != (a, b):
        raise ArithmeticError("invariant coefficient round trip failed")
    conjugate = ((a - linear * b) % prime, (-b) % prime)
    conjugate_trace = (2 * conjugate[0] - linear * conjugate[1]) % prime
    if conjugate_trace != trace:
        raise ArithmeticError("trace is not Frobenius invariant")
    encoded_count += 1
    return {
        "trace": trace,
        "anti_invariant_coefficient": b,
        "coefficient_discriminant": coefficient_discriminant,
        "norm": norm,
    }


def encode_polynomial(values):
    return [encode(value) for value in values]


def encode_rational(record):
    return {
        "numerator_low_to_high": encode_polynomial(record["numerator_coefficients_low_to_high_1_r"]),
        "denominator_low_to_high": encode_polynomial(record["denominator_coefficients_low_to_high_1_r"]),
        "degrees_numerator_denominator": record["degrees_numerator_denominator"],
    }


def encode_joint(record):
    return {
        "numerator_low_to_high_auxiliary_power_then_V": [
            encode_polynomial(value)
            for value in record["numerator_coefficients_low_to_high_auxiliary_power_then_V"]
        ],
        "denominator_low_to_high_auxiliary_power_then_V": [
            encode_polynomial(value)
            for value in record["denominator_coefficients_low_to_high_auxiliary_power_then_V"]
        ],
        "degrees_V_numerator_denominator": record["degrees_V_numerator_denominator"],
        "degrees_W_numerator_denominator": record["degrees_W_numerator_denominator"],
    }


invariant_payload = {
    "resolved_parent_and_pencil": {
        "horizontal": {
            name: encode_polynomial(values)
            for name, values in resolved["horizontal"].items()
        },
        "moving_equation_terms_T_W_old_x": [
            [t_degree, w_degree, x_degree, encode(coefficient)]
            for t_degree, w_degree, x_degree, coefficient
            in resolved["moving_equation"]["terms_T_W_x_coefficient_1_r"]
        ],
    },
    "generic_long_weierstrass": {
        **{
            name: encode_rational(long_model["weierstrass"][name])
            for name in ("a1", "a2", "a3", "a4", "a6")
        },
        "discriminant": encode_rational(long_model["discriminant"]),
        "j": encode_rational(long_model["j"]),
    },
    "minimal_short_weierstrass": {
        "A_low_to_high": encode_polynomial(
            minimal["minimal_short_weierstrass"]["A_coefficients_low_to_high_1_r"]
        ),
        "B_low_to_high": encode_polynomial(
            minimal["minimal_short_weierstrass"]["B_coefficients_low_to_high_1_r"]
        ),
        "discriminant_low_to_high": encode_polynomial(
            minimal["minimal_short_weierstrass"]["discriminant_coefficients_low_to_high_1_r"]
        ),
    },
    "birational_maps": {
        "forward_long": {
            coordinate: [encode_joint(record) for record in maps["forward_long"][coordinate]]
            for coordinate in ("X", "Y")
        },
        "inverse_long": {
            target: {
                "weighted_bound": maps["inverse_long"][target]["weighted_bound"],
                "monomials_X_power_Y_power": maps["inverse_long"][target]["monomials_X_power_Y_power"],
                "numerator_coefficients": [
                    encode_rational(value)
                    for value in maps["inverse_long"][target]["numerator_coefficients"]
                ],
                "denominator_coefficients": [
                    encode_rational(value)
                    for value in maps["inverse_long"][target]["denominator_coefficients"]
                ],
            }
            for target in ("W", "old_x")
        },
    },
}
canonical_json = json.dumps(invariant_payload, sort_keys=True, separators=(",", ":"))
candidate = horizontal["pairwise_producer"]["candidates"][0]
if candidate["relative_sign"] not in (-1, 1):
    raise ArithmeticError("common-producer horizontal sign is not explicit")

output = {
    "schema": "elkies-k3.q80-third-q12-frobenius-invariants-modp2.v2",
    "status": "PASS_EXACT_FROBENIUS_INVARIANT_THIRD_Q12_ENCODING_COMMON_PRODUCER",
    "specialization": specialization,
    "quadratic_extension": {
        "generator_equation": modulus + "=0",
        "frobenius": f"r |-> {-linear % prime}-r",
        "anti_invariant_generator": f"eta=2*r+{linear}",
        "eta_square_quadratic_discriminant": field_discriminant,
        "coefficient_reconstruction": "c=(trace+anti_invariant_coefficient*eta)/2",
        "frobenius_on_encoding": "trace, norm, and coefficient_discriminant fixed; anti coefficient negated",
    },
    "encoded_coefficients": invariant_payload,
    "gauge_ledger": {
        "section_sign": {
            "left_support_index_zero_based": candidate["left_index_zero_based"],
            "right_support_index_zero_based": candidate["right_index_zero_based"],
            "relative_sign": candidate["relative_sign"],
            "action": "P -> -P is independent of coefficient Frobenius",
        },
        "base_PGL2": {
            "current_coordinate": "V=N1/N0 in the RREF-normalized ordered resolved-kernel basis",
            "old_zero_I6_anchor": marking["base_alignment"],
            "canonicalization_status": "kernel columns and RREF pivots fix the displayed base gauge",
        },
        "weierstrass": {
            "long_Laurent_gauge": long_model["gauge"],
            "long_to_minimal": maps["long_to_minimal"],
            "canonicalization_status": "simple-branch Laurent leading terms and monic gauge-pole scaling are fixed",
        },
    },
    "validation": {
        "encoded_coefficient_count": encoded_count,
        "all_round_trips": True,
        "frobenius_trace_discriminant_invariance": True,
        "canonical_invariant_payload_sha256": hashlib.sha256(canonical_json.encode()).hexdigest(),
    },
    "inputs": [
        {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
        for path in (
            args.resolved, args.long, args.minimal, args.maps,
            args.marking, args.horizontal,
        )
    ],
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "lossless trace/discriminant/anti-invariant encoding of every displayed parent, child, and map coefficient",
            "literal coefficient round trip and Frobenius sign law",
            "separate pinned ledgers for section sign, base PGL2, and Weierstrass scaling",
        ],
        "not_proved": [
            "a common characteristic-zero quadratic field from two residues",
            "a characteristic-zero coefficient reconstruction",
        ],
    },
    "reproduce": (
        "python3 elkies-k3/scripts/encode_q80_third_q12_frobenius_invariants_modp2.py "
        f"--resolved {args.resolved} --long {args.long} --minimal {args.minimal} "
        f"--maps {args.maps} --marking {args.marking} --horizontal {args.horizontal} "
        f"--output {args.output}"
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"Q80THIRDQ12COMMONFROBENIUS|prime={prime}|coefficients={encoded_count}|"
    f"Delta={field_discriminant}|gauges=section_sign,base_PGL2,weierstrass_scaling|"
    "status=PASS_EXACT_FROBENIUS_INVARIANT_THIRD_Q12_ENCODING_COMMON_PRODUCER"
)
