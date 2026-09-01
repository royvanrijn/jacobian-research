#!/usr/bin/env python3
"""Certify prime-independent alignment of the p=19 and p=61 q12 pencils.

This aligns the producer orbit and resolved divisor, not finite-field
coefficient gauges.  The closure RUR may split differently at the two primes.
"""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"


def defaults(prime):
    common = "common-" if prime == 19 else ""
    return {
        "horizontal": RESULTS / f"q80-third-q12-um2-p{prime}-common-producer-horizontal.json",
        "pencil": RESULTS / f"q80-third-q12-um2-p{prime}-{common}resolved-pencil.json",
        "genus": RESULTS / f"q80-third-q12-um2-p{prime}-{common}resolved-genus.json",
    }


parser = argparse.ArgumentParser(description=__doc__)
for prime in (19, 61):
    for kind, path in defaults(prime).items():
        parser.add_argument(f"--p{prime}-{kind}", type=Path, default=path)
parser.add_argument(
    "--output", type=Path,
    default=RESULTS / "q80-third-q12-um2-p19-p61-common-producer-alignment.json",
)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


paths = {
    p: {kind: getattr(args, f"p{p}_{kind}") for kind in defaults(p)}
    for p in (19, 61)
}
data = {
    p: {kind: json.loads(path.read_text()) for kind, path in group.items()}
    for p, group in paths.items()
}


def signature(prime):
    horizontal, pencil, genus = (
        data[prime][kind] for kind in ("horizontal", "pencil", "genus")
    )
    expected_status = {
        "horizontal": "PASS_EXACT_THIRD_Q12_HORIZONTAL_FROM_COMMON_CLOSURE_PRODUCER",
        "pencil": "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_COMMON_PRODUCER",
        "genus": "PASS_EXACT_THIRD_Q12_GENUS_ONE_COMMON_PRODUCER",
    }
    for kind, payload in data[prime].items():
        if payload.get("status") != expected_status[kind]:
            raise ValueError(f"p={prime}: uncertified {kind}")
        specialization = payload["specialization"]
        if specialization["prime"] != prime or specialization["u"] != "-2":
            raise ValueError(f"p={prime}: specialization mismatch")

    pairwise = horizontal["pairwise_producer"]
    ambient = pencil["saturated_ambient"]
    gates = pencil["resolved_gates"]
    return {
        "u": "-2",
        "closure": {
            "quotient_degree_with_multiplicity": horizontal["rur"]["quotient_degree_with_multiplicity"],
            "squarefree_support_degree": horizontal["rur"]["squarefree_support_degree"],
            "decoded_geometric_support_points": len(horizontal["decoded_sections"]),
        },
        "horizontal_orbit": {
            "acceptance_profile": pairwise["acceptance_profile"],
            "pairs_with_relative_sign_tested": pairwise["pairs_with_relative_sign_tested"],
            "unsigned_target_hits": pairwise["unsigned_target_hits"],
            "frobenius_orbits_up_to_sign": pairwise["frobenius_orbits_up_to_sign"],
        },
        "resolved_divisor": {
            "smith_degrees": ambient["smith_degrees"],
            "ambient_dimension": ambient["dimension"],
            "generator_weights": ambient["generator_weights"],
            "column_labels_generator_degree": ambient["column_labels_generator_degree"],
            "D7_gate_rank": gates["D7"]["rank"],
            "combined_gate_rank": gates["combined_rank"],
            "kernel_dimension": gates["kernel_dimension"],
            "moving_degrees_T_W_x": pencil["moving_equation"]["degrees_T_W_x"],
        },
        "genus_one": {
            "lattice": genus["lattice"],
            "linear_system": genus["linear_system"],
        },
    }


signatures = {str(p): signature(p) for p in (19, 61)}
if signatures["19"] != signatures["61"]:
    raise ArithmeticError("the common-producer signatures do not align")
common = signatures["19"]
observed = (
    tuple(common["closure"].values()),
    tuple(common["horizontal_orbit"][key] for key in (
        "pairs_with_relative_sign_tested", "unsigned_target_hits",
        "frobenius_orbits_up_to_sign")),
    (common["resolved_divisor"]["smith_degrees"],
     common["resolved_divisor"]["ambient_dimension"],
     common["resolved_divisor"]["D7_gate_rank"],
     common["resolved_divisor"]["combined_gate_rank"],
     common["resolved_divisor"]["kernel_dimension"],
     common["resolved_divisor"]["moving_degrees_T_W_x"]),
)
expected = ((16, 12, 12), (156, 2, 1), ([0, 0, 6], 7, 4, 5, 2, [2, 9, 3]))
if observed != expected:
    raise ArithmeticError(("unexpected aligned signature", observed))

output = {
    "schema": "elkies-k3.q80-third-q12-common-producer-prime-alignment.v1",
    "status": "PASS_EXACT_SECOND_PRIME_COMMON_PRODUCER_ALIGNMENT",
    "inputs": {
        str(p): {
            kind: {"path": str(path.resolve().relative_to(ROOT)), "sha256": sha256(path)}
            for kind, path in group.items()
        } for p, group in paths.items()
    },
    "aligned_primes": [19, 61],
    "canonical_selection": {
        "specialization": "u=-2",
        "horizontal_equivalence": "section sign and Frobenius",
        "selection_rule": (
            "unique orbit among all pairwise differences of polynomial-closure support "
            "with P.O=2, height 8, and identity components at D7 and D5"
        ),
    },
    "common_signature": common,
    "prime_local_factorization": {
        str(p): {
            "extension_modulus": data[p]["pencil"]["specialization"]["extension_modulus"],
            "decoded_factor_degrees": sorted(
                row["factor_degree"] for row in data[p]["horizontal"]["decoded_sections"]
            ),
        } for p in (19, 61)
    },
    "claim_boundary": {
        "proved": [
            "the same exact producer schema canonically selects one horizontal orbit at p=19 and p=61",
            "both selected horizontals produce the same complete connected divisor signature",
            "both resolved pencils are primitive irreducible genus-one pencils",
        ],
        "not_proved": [
            "a coefficient-wise identification before fixing the quadratic generator and gauges",
            "a p=61 minimal Weierstrass child or birational map to it",
            "a characteristic-zero horizontal, pencil, or Jacobian equation",
        ],
    },
    "reproduce": "python3 elkies-k3/scripts/align_q80_third_q12_common_producer_primes.py",
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"alignment artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    "Q80THIRDQ12ALIGN|u=-2|primes=19,61|horizontal_orbits=1,1|"
    "pencil=primitive_genus_one|status=PASS_EXACT_SECOND_PRIME_COMMON_PRODUCER_ALIGNMENT"
)
