#!/usr/bin/env python3
"""Certify full finite-field child alignment between p=19 and p=61."""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
defaults = {
    "producer_alignment": RESULTS / "q80-third-q12-um2-p19-p61-common-producer-alignment.json",
    "p19_resolved": RESULTS / "q80-third-q12-um2-p19-resolved-pencil.json",
    "p19_long": RESULTS / "q80-third-q12-p19-jacobian-interpolated.json",
    "p19_minimal": RESULTS / "q80-third-q12-p19-jacobian-minimal.json",
    "p19_maps": RESULTS / "q80-third-q12-p19-birational-maps.json",
    "p19_marking": RESULTS / "q80-third-q12-p19-component-marking.json",
    "p19_frobenius": RESULTS / "q80-third-q12-p19-frobenius-invariants.json",
    "p61_resolved": RESULTS / "q80-third-q12-um2-p61-resolved-pencil.json",
    "p61_long": RESULTS / "q80-third-q12-p61-jacobian-interpolated.json",
    "p61_minimal": RESULTS / "q80-third-q12-p61-jacobian-minimal.json",
    "p61_maps": RESULTS / "q80-third-q12-p61-birational-maps.json",
    "p61_marking": RESULTS / "q80-third-q12-p61-component-marking.json",
    "p61_frobenius": RESULTS / "q80-third-q12-p61-frobenius-invariants.json",
}
for name, path in defaults.items():
    parser.add_argument("--" + name.replace("_", "-"), type=Path, default=path)
parser.add_argument(
    "--output", type=Path,
    default=RESULTS / "q80-third-q12-um2-p19-p61-full-child-alignment.json",
)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
paths = {name: getattr(args, name).resolve() for name in defaults}
args.output = args.output.resolve()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


data = {name: json.loads(path.read_text()) for name, path in paths.items()}
expected = {
    "producer_alignment": "PASS_EXACT_SECOND_PRIME_COMMON_PRODUCER_ALIGNMENT",
    "p19_resolved": "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_MOD19_QUADRATIC",
    "p19_long": "PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_MOD19_QUADRATIC",
    "p19_minimal": "PASS_EXACT_MINIMAL_THIRD_Q12_JACOBIAN_AND_FIBRES_MOD19_QUADRATIC",
    "p19_maps": "PASS_EXACT_GENERIC_THIRD_Q12_BIRATIONAL_MAPS_MOD19_QUADRATIC",
    "p19_marking": "PASS_EXACT_TRANSPORTED_THIRD_Q12_COMPONENT_MARKING_MOD19_QUADRATIC",
    "p19_frobenius": "PASS_EXACT_FROBENIUS_INVARIANT_THIRD_Q12_ENCODING_MOD19_QUADRATIC",
    "p61_resolved": "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_COMMON_PRODUCER",
    "p61_long": "PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_COMMON_PRODUCER",
    "p61_minimal": "PASS_EXACT_MINIMAL_THIRD_Q12_JACOBIAN_AND_FIBRES_COMMON_PRODUCER",
    "p61_maps": "PASS_EXACT_GENERIC_THIRD_Q12_BIRATIONAL_MAPS_COMMON_PRODUCER",
    "p61_marking": "PASS_EXACT_TRANSPORTED_THIRD_Q12_COMPONENT_MARKING_COMMON_PRODUCER",
    "p61_frobenius": "PASS_EXACT_FROBENIUS_INVARIANT_THIRD_Q12_ENCODING_COMMON_PRODUCER",
}
for name, status in expected.items():
    if data[name].get("status") != status:
        raise ValueError(f"uncertified full-child input: {name}")


def coefficient_shape(value):
    if isinstance(value, dict):
        if set(value) == {"trace", "anti_invariant_coefficient", "coefficient_discriminant", "norm"}:
            return "quadratic_coefficient_slot"
        return {key: coefficient_shape(item) for key, item in sorted(value.items())}
    if isinstance(value, list):
        return [coefficient_shape(item) for item in value]
    return value


shape19 = coefficient_shape(data["p19_frobenius"]["encoded_coefficients"])
shape61 = coefficient_shape(data["p61_frobenius"]["encoded_coefficients"])
if shape19 != shape61:
    raise ArithmeticError("Frobenius-invariant coefficient slots do not align")
counts = {
    "19": data["p19_frobenius"]["validation"]["encoded_coefficient_count"],
    "61": data["p61_frobenius"]["validation"]["encoded_coefficient_count"],
}
if counts != {"19": 1947, "61": 1947}:
    raise ArithmeticError("unexpected encoded coefficient counts")


def long_degrees(payload):
    return {
        name: payload["weierstrass"][name]["degrees_numerator_denominator"]
        for name in ("a1", "a2", "a3", "a4", "a6")
    } | {
        "Delta": payload["discriminant"]["degrees_numerator_denominator"],
        "j": payload["j"]["degrees_numerator_denominator"],
    }


if long_degrees(data["p19_long"]) != long_degrees(data["p61_long"]):
    raise ArithmeticError("generic long-model degree profiles disagree")


def map_signature(payload):
    return {
        "forward": {
            coordinate: [
                [record["degrees_V_numerator_denominator"], record["degrees_W_numerator_denominator"]]
                for record in payload["forward_long"][coordinate]
            ]
            for coordinate in ("X", "Y")
        },
        "inverse_bounds": {
            target: payload["inverse_long"][target]["weighted_bound"]
            for target in ("W", "old_x")
        },
        "validation": {
            key: payload["validation"][key]
            for key in (
                "generic_resolved_cubic_irreducible",
                "generic_forward_weierstrass_identity",
                "generic_inverse_W_identity",
                "generic_inverse_old_x_identity",
            )
        },
    }


if map_signature(data["p19_maps"]) != map_signature(data["p61_maps"]):
    raise ArithmeticError("generic birational-map signatures disagree")

for prime in (19, 61):
    minimal = data[f"p{prime}_minimal"]
    if minimal["minimal_short_weierstrass"]["degrees_A_B"] != [8, 12]:
        raise ArithmeticError(f"p={prime}: wrong minimal K3 degrees")
    if minimal["fibres"]["configuration"] != "I6+I4+3I2+8I1":
        raise ArithmeticError(f"p={prime}: wrong fibre configuration")
    if minimal["fibres"]["root_lattice"] != "A5+A3+3A1":
        raise ArithmeticError(f"p={prime}: wrong root marking")
    marking = data[f"p{prime}_marking"]
    if marking["transported_root_graph"]["type"] != "A5+A3+3A1":
        raise ArithmeticError(f"p={prime}: wrong transported graph")
    if marking["zero_orientation"]["selected_new_zero"] != "R5":
        raise ArithmeticError(f"p={prime}: wrong selected zero")
    if not marking["base_alignment"]["match"]:
        raise ArithmeticError(f"p={prime}: old-zero/I6 base mismatch")

labels19 = data["p19_resolved"]["saturated_ambient"]["column_labels_generator_degree"]
labels61 = data["p61_resolved"]["saturated_ambient"]["column_labels_generator_degree"]
if labels19 != labels61:
    raise ArithmeticError("ordered resolved-kernel columns disagree")
for prime in (19, 61):
    kernel = data[f"p{prime}_resolved"]["resolved_gates"]["kernel"]
    if kernel[0][:2] != [[1, 0], [0, 0]] or kernel[1][:2] != [[0, 0], [1, 0]]:
        raise ArithmeticError(f"p={prime}: resolved kernel is not in the common RREF gauge")

common_shape_json = json.dumps(shape19, sort_keys=True, separators=(",", ":"))
output = {
    "schema": "elkies-k3.q80-third-q12-full-child-prime-alignment.v1",
    "status": "PASS_EXACT_FULL_SECOND_PRIME_CHILD_ALIGNMENT",
    "specialization": {"u": "-2", "primes": [19, 61]},
    "common_child": {
        "source": "D7+D5/MW5",
        "neighbor": "q12",
        "target": "A5+A3+3A1/MW6",
        "minimal_fibres": "I6+I4+3I2+8I1",
        "selected_zero": "R5",
        "maps": "explicit generic birational maps in both directions",
    },
    "coefficient_alignment": {
        "ordered_frobenius_invariant_slots": 1947,
        "slot_shape_sha256": hashlib.sha256(common_shape_json.encode()).hexdigest(),
        "generator_free_values_for_crt": ["trace", "norm", "coefficient_discriminant"],
        "local_only_value": "anti_invariant_coefficient",
        "quadratic_generator_policy": (
            "no cross-prime generator identification is made; anti-invariant signs and scales "
            "are recovered only after a characteristic-zero quadratic field is reconstructed"
        ),
    },
    "gauge_alignment": {
        "section_sign": "selected by the common all-pairs producer and kept separate from Frobenius",
        "base_PGL2": (
            "V=N1/N0 in the common ordered RREF kernel gauge; the old zero lands on the I6 factor at both primes"
        ),
        "weierstrass_scaling": (
            "simple-branch Laurent leading terms fix the long gauge; monic d and q=d^2 fix the minimal scaling"
        ),
    },
    "degree_signatures": {
        "long_weierstrass": long_degrees(data["p19_long"]),
        "birational_maps": map_signature(data["p19_maps"]),
    },
    "inputs": {
        name: {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
        for name, path in paths.items()
    },
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "a complete mapped and minimal p=61 child independently aligned with the complete p=19 child",
            "the same A5+A3+3A1 component and R5 zero marking at both primes",
            "a common ordered generator-free coefficient schema suitable for CRT collection",
            "separate canonical ledgers for section sign, base PGL2, and Weierstrass scaling",
        ],
        "not_proved": [
            "a characteristic-zero quadratic coefficient field",
            "a characteristic-zero child equation or Mordell--Weil rank",
        ],
    },
    "reproduce": "python3 elkies-k3/scripts/align_q80_third_q12_full_children_primes.py",
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"full-child alignment artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    "Q80THIRDQ12FULLALIGN|u=-2|primes=19,61|child=A5+A3+3A1/MW6|"
    "maps=both|slots=1947|status=PASS_EXACT_FULL_SECOND_PRIME_CHILD_ALIGNMENT"
)
