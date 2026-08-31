#!/usr/bin/env sage
"""Regression test for the aligned Q80 endpoint CRT consumer."""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
RECONSTRUCTOR = ROOT / "elkies-k3/scripts/reconstruct_q80_alternate_rootless_crt_qq.sage"
SCHEMA = "elkies-k3.q80-alternate-rootless-aligned-modp.v1"
STATUS = "PASS_ALIGNED_Q80_ALTERNATE_ROOTLESS_MODP"
ROUTE = "Q80-q4-q4-q12-q12-alternate-q4-q6"

ring = PolynomialRing(QQ, "u")
u = ring.gen()
A = u**8 + u + 1
B = u**12 + u**2 + 1
delta = -16 * (4 * A**3 + 27 * B**2)
assert delta.degree() == 24 and delta.is_squarefree()
primes = (101, 103, 107, 109, 113)
normalization = {
    "base": "u",
    "base_marking": ["u=0", "u=1", "u=infinity"],
    "weierstrass_gauge": "short model with the retained compiler scaling",
    "regression_only": True,
}


def reduce_coefficient(coefficient, prime):
    coefficient = QQ(coefficient)
    return int(
        coefficient.numerator()
        * pow(int(coefficient.denominator()), -1, prime)
        % prime
    )


with tempfile.TemporaryDirectory(prefix="q80-rootless-crt-regression-") as directory:
    temporary = Path(directory)
    inputs = []
    for prime in primes:
        payload = {
            "schema": SCHEMA,
            "status": STATUS,
            "route_id": ROUTE,
            "prime": prime,
            "canonical_marking_id": "regression-fixed-base-and-short-gauge-v1",
            "marking_transport": {
                "route_map_chain_id": "regression-q80-complete-map-chain-v1",
                "parent_child_map_digests_modp": {
                    "q80_to_first_q4": f"regression-{prime}-01",
                    "first_q4_to_second_q4": f"regression-{prime}-02",
                    "second_q4_to_third_q12": f"regression-{prime}-03",
                    "third_q12_to_fourth_q12": f"regression-{prime}-04",
                    "fourth_q12_to_alternate_q4": f"regression-{prime}-05",
                    "alternate_q4_to_rootless_q6": f"regression-{prime}-06",
                },
            },
            "normalization": normalization,
            "model": {
                "A_coefficients_low_to_high": [
                    reduce_coefficient(value, prime) for value in A.list()
                ],
                "B_coefficients_low_to_high": [
                    reduce_coefficient(value, prime) for value in B.list()
                ],
                "A_degree": 8,
                "B_degree": 12,
                "root_rank": 0,
            },
        }
        path = temporary / f"endpoint-mod-{prime}.json"
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        inputs.append(path)

    output = temporary / "endpoint-qq.json"
    command = [
        sys.executable,
        str(RECONSTRUCTOR),
        *map(str, inputs),
        "--heldout-count",
        "2",
        "--output",
        str(output),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode:
        raise RuntimeError(completed.stdout + completed.stderr)
    result = json.loads(output.read_text())
    assert result["status"] == "PASS_Q80_ALTERNATE_ROOTLESS_MW17_CRT_HELDOUT"
    family = result["rootless_family"]
    assert ring(family["A_coefficients_low_to_high"]) == A
    assert ring(family["B_coefficients_low_to_high"]) == B
    assert family["classification"]["rootless"]
    assert family["classification"]["root_rank"] == 0
    assert family["classification"]["euler_number"] == 24
    assert family["classification"]["MW_rank_at_Picard_19"] == 17
    assert result["reconstruction"]["training_primes"] == list(primes[:3])
    assert result["reconstruction"]["heldout_primes"] == list(primes[3:])
    assert result["reconstruction"]["all_training_coefficients_replayed"]
    assert result["reconstruction"]["all_heldout_coefficients_replayed"]

    j_record = family["reduced_j"]
    numerator = ring(j_record["primitive_integer_numerator_coefficients_low_to_high"])
    denominator = ring(j_record["primitive_integer_denominator_coefficients_low_to_high"])
    expected_numerator = (-48 * A) ** 3
    expected_denominator = delta
    assert numerator * expected_denominator == denominator * expected_numerator

    checked = subprocess.run(
        command + ["--check"], text=True, capture_output=True, check=False
    )
    if checked.returncode:
        raise RuntimeError(checked.stdout + checked.stderr)

    # A held-out coefficient change must invalidate the reconstruction.
    last_path = inputs[-1]
    last_payload = json.loads(last_path.read_text())
    original_constant = last_payload["model"]["A_coefficients_low_to_high"][0]
    last_payload["model"]["A_coefficients_low_to_high"][0] = (
        original_constant + 1
    ) % primes[-1]
    last_path.write_text(json.dumps(last_payload, indent=2, sort_keys=True) + "\n")
    tampered = subprocess.run(command, text=True, capture_output=True, check=False)
    assert tampered.returncode != 0
    last_payload["model"]["A_coefficients_low_to_high"][0] = original_constant
    last_path.write_text(json.dumps(last_payload, indent=2, sort_keys=True) + "\n")

    # A different prime-specific marking must be rejected before CRT.
    last_payload["canonical_marking_id"] = "deliberately-misaligned-regression"
    last_path.write_text(json.dumps(last_payload, indent=2, sort_keys=True) + "\n")
    misaligned = subprocess.run(command, text=True, capture_output=True, check=False)
    assert misaligned.returncode != 0
    assert "canonical marking identifiers disagree" in (
        misaligned.stdout + misaligned.stderr
    )
    last_payload["canonical_marking_id"] = "regression-fixed-base-and-short-gauge-v1"
    last_payload["marking_transport"]["parent_child_map_digests_modp"].pop(
        "third_q12_to_fourth_q12"
    )
    last_path.write_text(json.dumps(last_payload, indent=2, sort_keys=True) + "\n")
    missing_map_slot = subprocess.run(
        command, text=True, capture_output=True, check=False
    )
    assert missing_map_slot.returncode != 0
    assert "parent/child map slots disagree" in (
        missing_map_slot.stdout + missing_map_slot.stderr
    )

print(
    "Q80ROOTLESSCRTREGRESSION|training=3|heldout=2|"
    "fibres=24I1|root_rank=0|MW=17|j_replayed=1|"
    "tampered_heldout_rejected=1|misaligned_marking_rejected=1|"
    "missing_map_slot_rejected=1|status=PASS",
    flush=True,
)
