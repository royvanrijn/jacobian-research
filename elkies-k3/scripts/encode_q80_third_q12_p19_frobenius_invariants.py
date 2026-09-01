#!/usr/bin/env python3
"""Encode the complete p=19 child and maps in a Frobenius-invariant schema."""

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESOLVED = ROOT / "artifacts/generated-results/q80-third-q12-um2-p19-resolved-pencil.json"
LONG = ROOT / "artifacts/generated-results/q80-third-q12-p19-jacobian-interpolated.json"
MINIMAL = ROOT / "artifacts/generated-results/q80-third-q12-p19-jacobian-minimal.json"
MAPS = ROOT / "artifacts/generated-results/q80-third-q12-p19-birational-maps.json"
MARKING = ROOT / "artifacts/generated-results/q80-third-q12-p19-component-marking.json"
HORIZONTAL = ROOT / "artifacts/generated-results/q80-fixed-u-minus2-p19-po0-rur-third-q12-modp.json"
OUTPUT = ROOT / "artifacts/generated-results/q80-third-q12-p19-frobenius-invariants.json"

PRIME = 19
FIELD_DISCRIMINANT = 18
INVERSE_TWO = pow(2, -1, PRIME)
encoded_count = 0


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


resolved = json.loads(RESOLVED.read_text())
long_model = json.loads(LONG.read_text())
minimal = json.loads(MINIMAL.read_text())
maps = json.loads(MAPS.read_text())
marking = json.loads(MARKING.read_text())
horizontal = json.loads(HORIZONTAL.read_text())
expected = (
    (resolved, "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_MOD19_QUADRATIC"),
    (long_model, "PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_MOD19_QUADRATIC"),
    (minimal, "PASS_EXACT_MINIMAL_THIRD_Q12_JACOBIAN_AND_FIBRES_MOD19_QUADRATIC"),
    (maps, "PASS_EXACT_GENERIC_THIRD_Q12_BIRATIONAL_MAPS_MOD19_QUADRATIC"),
    (marking, "PASS_EXACT_TRANSPORTED_THIRD_Q12_COMPONENT_MARKING_MOD19_QUADRATIC"),
    (horizontal, "PASS_EXACT_MODP2_THIRD_Q12_HORIZONTAL_FROBENIUS_ORBIT"),
)
if any(payload.get("status") != status for payload, status in expected):
    raise ValueError("one or more Frobenius-encoding inputs are not certified")


def encode(coordinates):
    """Encode c=a+b*r as (Tr(c), b, (c-cbar)^2), with eta=2r+12."""
    global encoded_count
    a, b = (int(coordinates[0]) % PRIME, int(coordinates[1]) % PRIME)
    trace = (2 * a + 7 * b) % PRIME
    anti = b
    coefficient_discriminant = FIELD_DISCRIMINANT * anti * anti % PRIME
    norm = (trace * trace - coefficient_discriminant) * pow(4, -1, PRIME) % PRIME

    # c=(trace+anti*eta)/2 and eta=2r+12.
    reconstructed_a = (trace + 12 * anti) * INVERSE_TWO % PRIME
    reconstructed_b = anti
    if (reconstructed_a, reconstructed_b) != (a, b):
        raise ArithmeticError("invariant coefficient round trip failed")

    # Frobenius sends (a,b) to (a+7b,-b), fixing trace and discriminant and
    # negating only the anti-invariant coefficient.
    conjugate_a = (a + 7 * b) % PRIME
    conjugate_b = (-b) % PRIME
    conjugate_trace = (2 * conjugate_a + 7 * conjugate_b) % PRIME
    if conjugate_trace != trace:
        raise ArithmeticError("trace is not Frobenius invariant")
    encoded_count += 1
    return {
        "trace": trace,
        "anti_invariant_coefficient": anti,
        "coefficient_discriminant": coefficient_discriminant,
        "norm": norm,
    }


def encode_polynomial(coordinates):
    return [encode(value) for value in coordinates]


def encode_rational(record):
    return {
        "numerator_low_to_high": encode_polynomial(
            record["numerator_coefficients_low_to_high_1_r"]
        ),
        "denominator_low_to_high": encode_polynomial(
            record["denominator_coefficients_low_to_high_1_r"]
        ),
        "degrees_numerator_denominator": record["degrees_numerator_denominator"],
    }


resolved_encoding = {
    "horizontal": {
        name: encode_polynomial(values)
        for name, values in resolved["horizontal"].items()
    },
    "moving_equation_terms_T_W_old_x": [
        [t_degree, w_degree, x_degree, encode(coefficient)]
        for t_degree, w_degree, x_degree, coefficient in resolved["moving_equation"][
            "terms_T_W_x_coefficient_1_r"
        ]
    ],
}

long_encoding = {
    name: encode_rational(long_model["weierstrass"][name])
    for name in ("a1", "a2", "a3", "a4", "a6")
}
long_encoding["discriminant"] = encode_rational(long_model["discriminant"])
long_encoding["j"] = encode_rational(long_model["j"])

minimal_encoding = {
    "A_low_to_high": encode_polynomial(
        minimal["minimal_short_weierstrass"]["A_coefficients_low_to_high_1_r"]
    ),
    "B_low_to_high": encode_polynomial(
        minimal["minimal_short_weierstrass"]["B_coefficients_low_to_high_1_r"]
    ),
    "discriminant_low_to_high": encode_polynomial(
        minimal["minimal_short_weierstrass"]["discriminant_coefficients_low_to_high_1_r"]
    ),
}


def encode_joint_map_block(record):
    return {
        "numerator_low_to_high_auxiliary_power_then_V": [
            encode_polynomial(polynomial)
            for polynomial in record[
                "numerator_coefficients_low_to_high_auxiliary_power_then_V"
            ]
        ],
        "denominator_low_to_high_auxiliary_power_then_V": [
            encode_polynomial(polynomial)
            for polynomial in record[
                "denominator_coefficients_low_to_high_auxiliary_power_then_V"
            ]
        ],
        "degrees_V_numerator_denominator": record["degrees_V_numerator_denominator"],
        "degrees_W_numerator_denominator": record["degrees_W_numerator_denominator"],
    }


forward_encoding = {
    coordinate: [encode_joint_map_block(record) for record in maps["forward_long"][coordinate]]
    for coordinate in ("X", "Y")
}
inverse_encoding = {}
for target in ("W", "old_x"):
    record = maps["inverse_long"][target]
    inverse_encoding[target] = {
        "weighted_bound": record["weighted_bound"],
        "monomials_X_power_Y_power": record["monomials_X_power_Y_power"],
        "numerator_coefficients": [encode_rational(value) for value in record["numerator_coefficients"]],
        "denominator_coefficients": [encode_rational(value) for value in record["denominator_coefficients"]],
    }

candidate = horizontal["third_q12"]["candidates_up_to_sign"][0]
if candidate["polynomial_section_sign"] not in (-1, 1):
    raise ArithmeticError("horizontal section sign is not explicit")

invariant_payload = {
    "resolved_parent_and_pencil": resolved_encoding,
    "generic_long_weierstrass": long_encoding,
    "minimal_short_weierstrass": minimal_encoding,
    "birational_maps": {
        "forward_long": forward_encoding,
        "inverse_long": inverse_encoding,
    },
}
canonical_invariant_json = json.dumps(invariant_payload, sort_keys=True, separators=(",", ":"))

output = {
    "schema": "elkies-k3.q80-third-q12-frobenius-invariants-modp2.v1",
    "status": "PASS_EXACT_FROBENIUS_INVARIANT_THIRD_Q12_ENCODING_MOD19_QUADRATIC",
    "specialization": {"u": "-2", "prime": PRIME},
    "quadratic_extension": {
        "generator_equation": "r^2+12*r+3=0",
        "frobenius": "r |-> 7-r",
        "anti_invariant_generator": "eta=2*r+12",
        "eta_square_quadratic_discriminant": FIELD_DISCRIMINANT,
        "coefficient_reconstruction": "c=(trace+anti_invariant_coefficient*eta)/2",
        "frobenius_on_encoding": "trace and coefficient_discriminant fixed; anti_invariant_coefficient negated",
    },
    "encoded_coefficients": invariant_payload,
    "gauge_ledger": {
        "section_sign": {
            "chosen_polynomial_section_index_one_based": candidate[
                "polynomial_section_index_one_based"
            ],
            "chosen_sign": candidate["polynomial_section_sign"],
            "action": "horizontal P -> -P fixes x and negates the parent y-coordinate; it is not absorbed into Frobenius",
        },
        "base_PGL2": {
            "current_coordinate": "V=N1/N0 in the ordered resolved-kernel basis",
            "current_I6_anchor": marking["base_alignment"],
            "action": "V -> (a*V+b)/(c*V+d), tracked independently of coefficient conjugation",
            "canonicalization_status": "current kernel gauge retained; no cross-prime PGL2 normalization asserted",
        },
        "weierstrass": {
            "long_Laurent_gauge": long_model["gauge"],
            "long_to_minimal": maps["long_to_minimal"],
            "residual_scaling_action": "(X,Y,A,B) -> (lambda^2*X,lambda^3*Y,lambda^4*A,lambda^6*B)",
            "canonicalization_status": "Laurent leading coefficients and the displayed q(V) scaling are pinned; residual cross-prime scaling remains separate",
        },
    },
    "validation": {
        "encoded_coefficient_count": encoded_count,
        "all_round_trips": True,
        "frobenius_trace_discriminant_invariance": True,
        "canonical_invariant_payload_sha256": hashlib.sha256(
            canonical_invariant_json.encode()
        ).hexdigest(),
    },
    "inputs": [
        {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
        for path in (RESOLVED, LONG, MINIMAL, MAPS, MARKING, HORIZONTAL)
    ],
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "lossless trace/discriminant/anti-invariant encoding of every displayed parent, child, and birational-map coefficient",
            "literal coefficient round trip and Frobenius sign law",
            "separate ledgers for section sign, base PGL2, and Weierstrass scaling",
        ],
        "not_proved": [
            "a canonical PGL2 or Weierstrass alignment at a second prime",
            "a characteristic-zero coefficient reconstruction",
        ],
    },
    "reproduce": "python3 elkies-k3/scripts/encode_q80_third_q12_p19_frobenius_invariants.py",
}
OUTPUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"Q80THIRDQ12FROBENIUS|coefficients={encoded_count}|Delta={FIELD_DISCRIMINANT}|"
    "gauges=section_sign,base_PGL2,weierstrass_scaling|"
    "status=PASS_EXACT_FROBENIUS_INVARIANT_THIRD_Q12_ENCODING_MOD19_QUADRATIC"
)
