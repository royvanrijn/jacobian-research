#!/usr/bin/env sage-python
"""Independently compute n=1,2 Frobenius moments for two singleton twists.

This is the fibrewise quadratic-character calculation, independent of the
toric controlled-reduction backend and its raw-output parser.
"""

from hashlib import sha256
import json
from pathlib import Path
import runpy

import numpy as np
from sage.all import GF, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
BISECTIONS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-cheapest-1024-v1.json"
)
MODEL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)
PRODUCT_AUDITOR = (
    ROOT / "elkies-k3/scripts/audit_r17_norm12_11952_product_twist_finite_field_bounds.sage"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-singleton-twist-finite-field-bound-audit-v1.json"
)
LABELS = ("alternate-orbit-0fda0", "alternate-orbit-1037d")
PRIMES = (131, 137, 151, 157)


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


functions = runpy.run_path(str(PRODUCT_AUDITOR), run_name="r17_product_trace_library")
legendre_table = functions["legendre_table"]
first_nonsquare = functions["first_nonsquare"]
fp_power_sums = functions["fp_power_sums"]
fp2_power_sums = functions["fp2_power_sums"]

bisections = json.loads(BISECTIONS.read_text())
model = json.loads(MODEL.read_text())
if bisections.get("schema") != "elkies-k3.bisection-extension-input.v1":
    raise ValueError("unexpected bisection schema")
if model.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
    raise ValueError("unexpected direct model")
by_label = {record["label"]: record for record in bisections["bisections"]}
weierstrass = model["weierstrass_model"]

records = {label: [] for label in LABELS}
for prime in PRIMES:
    field = GF(prime)
    ring = PolynomialRing(field, "t")
    A = ring([field(QQ(value)) for value in weierstrass["A_coefficients_low_to_high"]])
    B = ring([field(QQ(value)) for value in weierstrass["B_coefficients_low_to_high"]])
    Delta = -field(16) * (field(4) * A**3 + field(27) * B**2)
    twists = [
        ring(
            [
                field(QQ(value))
                for value in by_label[label]["branch"]["numerator_coefficients"]
            ]
        )
        for label in LABELS
    ]
    if (
        A.degree() != 8
        or B.degree() != 12
        or Delta.degree() != 24
        or not Delta.is_squarefree()
        or any(
            twist.degree() != 2
            or not twist.is_squarefree()
            or twist.gcd(Delta).degree()
            for twist in twists
        )
    ):
        raise ArithmeticError(f"bad singleton reduction at {prime}")

    a_coefficients = [int(A[index]) for index in range(9)]
    b_coefficients = [int(B[index]) for index in range(13)]
    # The shared product-twist routines expect the degree-four infinity slot.
    # Padding with zeros leaves every finite evaluation unchanged and suppresses
    # the infinity term, which is added explicitly below with the degree-two
    # leading coefficient.
    padded_twists = [
        [int(twist[index]) for index in range(3)] + [0, 0] for twist in twists
    ]
    first_sums, _local = fp_power_sums(
        a_coefficients, b_coefficients, padded_twists, prime
    )
    second_sums, nonsquare, fp2_local_sha = fp2_power_sums(
        a_coefficients, b_coefficients, padded_twists, prime
    )

    character = legendre_table(prime)
    xs = np.arange(prime, dtype=np.int64)
    infinity_rhs = (
        xs * xs % prime * xs
        + a_coefficients[8] * xs
        + b_coefficients[12]
    ) % prime
    infinity_fibre_sum_fp = int(character[infinity_rhs].sum())

    coordinates = np.arange(prime, dtype=np.int64)
    xa = np.repeat(coordinates, prime)
    xb = np.tile(coordinates, prime)
    x2a = (xa * xa + nonsquare * xb * xb) % prime
    x2b = 2 * xa * xb % prime
    x3a = (x2a * xa + nonsquare * x2b * xb) % prime
    x3b = (x2a * xb + x2b * xa) % prime
    rhs_a = (x3a + a_coefficients[8] * xa + b_coefficients[12]) % prime
    rhs_b = (x3b + a_coefficients[8] * xb) % prime
    infinity_norm = (rhs_a * rhs_a - nonsquare * rhs_b * rhs_b) % prime
    infinity_fibre_sum_fp2 = int(character[infinity_norm].sum())

    for index, (label, twist) in enumerate(zip(LABELS, twists)):
        first_sums[index] += infinity_fibre_sum_fp * int(character[int(twist[2])])
        # Every nonzero F_p element is a square in F_{p^2}.
        second_sums[index] += infinity_fibre_sum_fp2
        records[label].append(
            {
                "prime": prime,
                "good_reduction": True,
                "elliptic_L_frobenius_power_sums_n1_n2": [
                    int(first_sums[index]),
                    int(second_sums[index]),
                ],
                "fp2_nonsquare": nonsquare,
                "fp2_finite_local_character_sums_sha256": fp2_local_sha,
            }
        )

payload = {
    "schema": "elkies-k3.r17-singleton-twist-finite-field-bound-audit.v1",
    "status": "PASS_INDEPENDENT_FIBREWISE_N1_N2_MOMENTS",
    "labels": list(LABELS),
    "primes": list(PRIMES),
    "targets": [
        {"label": label, "reductions": records[label]} for label in LABELS
    ],
    "method": (
        "Exact vectorized quadratic-character sums over P1(F_p) and P1(F_p^2); "
        "independent of ToricControlledReduction"
    ),
    "inputs": {
        str(path.relative_to(ROOT)): digest(path)
        for path in (Path(__file__).resolve(), BISECTIONS, MODEL, PRODUCT_AUDITOR)
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "R17SINGLETONFF|labels=0fda0,1037d|primes=131,137,151,157|moments=2|"
    f"output={OUTPUT}|status={payload['status']}",
    flush=True,
)
