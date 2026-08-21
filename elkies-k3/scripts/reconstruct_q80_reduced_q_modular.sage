#!/usr/bin/env sage
"""Recover the degree-ten Q(u) quotient by modular gcd and CRT.

The exact bridge/subresultant composition has numerator and denominator of
degree 142 with a degree-132 common factor.  Direct characteristic-zero gcds
are needlessly expensive because the input coefficients have about 52,000
decimal digits.  This script computes the reduced quotient modulo deterministic
61-bit primes, rationally reconstructs its 22 coefficients, and accepts a
candidate only after an exact coefficient-wise cross-product identity.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


parser = argparse.ArgumentParser()
parser.add_argument(
    "--input",
    default=(
        "artifacts/generated-results/"
        "q80-cm24-slope-8-87-qq-reduced-composition.json"
    ),
)
parser.add_argument(
    "--output",
    default=(
        "artifacts/generated-results/"
        "q80-cm24-slope-8-87-qq-reduced-Q.json"
    ),
)
parser.add_argument("--maximum-primes", type=int, default=500)
parser.add_argument("--reconstruct-every", type=int, default=50)
arguments = parser.parse_args()


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


input_path = Path(arguments.input)
payload = json.loads(input_path.read_text())
if payload.get("schema") != "q80-cm24-qq-reduced-composition-v1":
    raise ValueError("unexpected reduced-composition schema")

polynomial_ring = PolynomialRing(QQ, "u")
u = polynomial_ring.gen()
unreduced_numerator = polynomial_ring(payload["Q_unreduced_numerator"])
unreduced_denominator = polynomial_ring(payload["Q_unreduced_denominator"])
if (unreduced_numerator.degree(), unreduced_denominator.degree()) != (142, 142):
    raise ValueError("expected the degree-142 subresultant composition")

common_denominator = lcm(
    coefficient.denominator()
    for polynomial in (unreduced_numerator, unreduced_denominator)
    for coefficient in polynomial.coefficients()
)
integral_numerator = polynomial_ring(common_denominator*unreduced_numerator)
integral_denominator = polynomial_ring(common_denominator*unreduced_denominator)
if any(
    coefficient.denominator() != 1
    for polynomial in (integral_numerator, integral_denominator)
    for coefficient in polynomial.coefficients()
):
    raise ArithmeticError("failed to clear the composition denominators")

expected_reduced_degree = 10
expected_gcd_degree = 132
prime_seed = ZZ(2)**60 + 1234567
primes = []
candidate_prime = next_prime(prime_seed)
while len(primes) < arguments.maximum_primes:
    if (
        integral_numerator.leading_coefficient() % candidate_prime
        and integral_denominator.leading_coefficient() % candidate_prime
    ):
        primes.append(ZZ(candidate_prime))
    candidate_prime = next_prime(candidate_prime+2)

# Reduce each huge integer only once modulo the product of all candidate
# primes.  Subsequent per-prime reductions act on numbers bounded by that
# much smaller product.
full_modulus = prod(primes)
numerator_mod_full = [
    ZZ(integral_numerator[index]) % full_modulus
    for index in range(143)
]
denominator_mod_full = [
    ZZ(integral_denominator[index]) % full_modulus
    for index in range(143)
]


def exact_cross_product_holds(reduced_numerator, reduced_denominator):
    maximum_degree = 142 + expected_reduced_degree
    for degree in range(maximum_degree+1):
        left = QQ.zero()
        right = QQ.zero()
        for quotient_degree in range(expected_reduced_degree+1):
            input_degree = degree-quotient_degree
            if 0 <= input_degree <= 142:
                left += (
                    integral_numerator[input_degree]
                    * reduced_denominator[quotient_degree]
                )
                right += (
                    integral_denominator[input_degree]
                    * reduced_numerator[quotient_degree]
                )
        if left != right:
            return False
    return True


crt_modulus = ZZ.one()
numerator_residues = [ZZ.zero()]*(expected_reduced_degree+1)
denominator_residues = [ZZ.zero()]*(expected_reduced_degree+1)
accepted_primes = []
bad_primes = []
reconstructed_pair = None

for prime in primes:
    finite_field = GF(prime)
    finite_ring = PolynomialRing(finite_field, "u")
    finite_numerator = finite_ring(
        [finite_field(value % prime) for value in numerator_mod_full]
    )
    finite_denominator = finite_ring(
        [finite_field(value % prime) for value in denominator_mod_full]
    )
    finite_gcd = finite_numerator.gcd(finite_denominator)
    if finite_gcd.degree() != expected_gcd_degree:
        bad_primes.append({"prime": int(prime), "gcd_degree": int(finite_gcd.degree())})
        continue
    reduced_numerator_mod_prime = finite_numerator//finite_gcd
    reduced_denominator_mod_prime = finite_denominator//finite_gcd
    if (
        reduced_numerator_mod_prime.degree() > expected_reduced_degree
        or reduced_denominator_mod_prime.degree() != expected_reduced_degree
    ):
        bad_primes.append(
            {
                "prime": int(prime),
                "gcd_degree": int(finite_gcd.degree()),
                "quotient_degrees": [
                    int(reduced_numerator_mod_prime.degree()),
                    int(reduced_denominator_mod_prime.degree()),
                ],
            }
        )
        continue
    normalization = ~reduced_denominator_mod_prime.leading_coefficient()
    reduced_numerator_mod_prime *= normalization
    reduced_denominator_mod_prime *= normalization

    inverse_old_modulus = inverse_mod(crt_modulus % prime, prime)
    for residues, polynomial in (
        (numerator_residues, reduced_numerator_mod_prime),
        (denominator_residues, reduced_denominator_mod_prime),
    ):
        for index in range(expected_reduced_degree+1):
            value = ZZ(polynomial[index])
            correction = (
                (value-residues[index]) % prime
                * inverse_old_modulus
            ) % prime
            residues[index] += crt_modulus*correction
    crt_modulus *= prime
    accepted_primes.append(prime)

    accepted_count = len(accepted_primes)
    if (
        accepted_count % arguments.reconstruct_every
        and accepted_count != arguments.maximum_primes
    ):
        continue
    try:
        rational_numerator = [
            residue.rational_reconstruction(crt_modulus)
            for residue in numerator_residues
        ]
        rational_denominator = [
            residue.rational_reconstruction(crt_modulus)
            for residue in denominator_residues
        ]
    except ArithmeticError:
        print(
            "Q80MODQ|stage=reconstruct|"
            f"accepted_primes={accepted_count}|modulus_bits={crt_modulus.nbits()}|"
            "status=insufficient_modulus",
            flush=True,
        )
        continue
    candidate_numerator = polynomial_ring(rational_numerator)
    candidate_denominator = polynomial_ring(rational_denominator)
    if exact_cross_product_holds(candidate_numerator, candidate_denominator):
        reconstructed_pair = candidate_numerator, candidate_denominator
        print(
            "Q80MODQ|stage=reconstruct|"
            f"accepted_primes={accepted_count}|modulus_bits={crt_modulus.nbits()}|"
            "status=PASS_EXACT_CROSS_PRODUCT",
            flush=True,
        )
        break
    print(
        "Q80MODQ|stage=reconstruct|"
        f"accepted_primes={accepted_count}|modulus_bits={crt_modulus.nbits()}|"
        "status=candidate_failed_exact_identity",
        flush=True,
    )

if reconstructed_pair is None:
    raise ArithmeticError(
        "modular quotient did not reconstruct within the configured prime bound"
    )

reduced_numerator, reduced_denominator = reconstructed_pair
coefficient_denominator = lcm(
    coefficient.denominator()
    for polynomial in reconstructed_pair
    for coefficient in polynomial.coefficients()
)
reduced_numerator = polynomial_ring(coefficient_denominator*reduced_numerator)
reduced_denominator = polynomial_ring(coefficient_denominator*reduced_denominator)
joint_content = gcd(
    tuple(
        ZZ(coefficient)
        for polynomial in (reduced_numerator, reduced_denominator)
        for coefficient in polynomial.coefficients()
    )
)
reduced_numerator = reduced_numerator//joint_content
reduced_denominator = reduced_denominator//joint_content
if reduced_denominator.leading_coefficient() < 0:
    reduced_numerator = -reduced_numerator
    reduced_denominator = -reduced_denominator
if not exact_cross_product_holds(reduced_numerator, reduced_denominator):
    raise ArithmeticError("primitive quotient failed the exact identity")

output_payload = {
    "schema": "q80-cm24-qq-reduced-Q-v1",
    "scope": "exact_degree_ten_quotient_of_degree_142_subresultant_composition",
    "status": "exact_cross_product_identity",
    "Q_numerator": str(reduced_numerator),
    "Q_denominator": str(reduced_denominator),
    "degrees": [int(reduced_numerator.degree()), int(reduced_denominator.degree())],
    "cancelled_degree": int(expected_gcd_degree),
    "accepted_primes": [int(prime) for prime in accepted_primes],
    "bad_primes": bad_primes,
    "crt_modulus_bits": int(crt_modulus.nbits()),
    "exact_check": "A_142*D_10 == B_142*N_10 coefficientwise over QQ",
    "input": {"path": str(input_path), "sha256": sha256(input_path)},
}
output_path = Path(arguments.output)
output_path.parent.mkdir(parents=True, exist_ok=True)
output_path.write_text(json.dumps(output_payload, indent=2, sort_keys=True)+"\n")

print(
    "Q80MODQ|degrees="
    f"{reduced_numerator.degree()},{reduced_denominator.degree()}|"
    f"cancelled_degree={expected_gcd_degree}|"
    f"accepted_primes={len(accepted_primes)}|output={output_path}|"
    "status=PASS_EXACT_REDUCED_Q",
    flush=True,
)
