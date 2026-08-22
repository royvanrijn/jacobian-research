#!/usr/bin/env sage -python
"""CRT-reconstruct the exact q8-child q-frame principal-part correction.

The correction is the unique R of degree below deg(Nx) satisfying

    R*h*Dy = Ny (mod Nx).

Raw rational extended Euclid is prohibitively expensive for the pinned
18--27K-bit section coefficients.  This script computes R modulo deterministic
good 61-bit primes, CRT-lifts its 96 coefficients, rationally reconstructs
them, validates representative accepted and withheld residues, and then
proves the exact congruence.  The final exact congruence, rather than the
modular spot checks, is the certificate.  It writes ``complete=false`` if its
requested prime budget is insufficient; such output is not a certificate.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, next_prime, prod


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-pole-normalization-crt.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def monic_power_root(value, exponent):
    root = value.parent().one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        root *= factor.monic() ** (multiplicity // exponent)
    return root.monic()


def residue(numerator, denominator, prime):
    denominator %= prime
    if not denominator:
        raise ZeroDivisionError("coefficient denominator vanishes modulo prime")
    return numerator % prime * denominator.inverse_mod(prime) % prime


def reduction_data(source):
    return tuple((ZZ(value.numerator()), ZZ(value.denominator())) for value in source.list())


def reduce_polynomial(source, finite_ring, prime):
    return finite_ring([
        finite_ring.base_ring()(residue(numerator, denominator, prime))
        for numerator, denominator in source
    ])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--maximum-primes", type=int, default=2000)
parser.add_argument("--minimum-primes", type=int, default=200)
parser.add_argument("--reconstruct-every", type=int, default=100)
parser.add_argument("--withheld-primes", type=int, default=3)
parser.add_argument(
    "--accepted-validation-primes",
    type=int,
    default=5,
    help="number of evenly distributed incorporated primes checked before the exact proof",
)
parser.add_argument(
    "--prime-bits",
    type=int,
    default=61,
    help="bit size of deterministic good primes (29 through 61)",
)
parser.add_argument(
    "--checkpoint",
    type=Path,
    help="write resumable CRT state after this batch (not a mathematical certificate)",
)
parser.add_argument(
    "--resume-from",
    type=Path,
    help="continue a checkpoint written by this exact script and input pair",
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if (
    args.maximum_primes < 1
    or args.minimum_primes < 1
    or args.reconstruct_every < 1
    or args.accepted_validation_primes < 1
):
    raise ValueError("prime counts must be positive")
if not 29 <= args.prime_bits <= 61:
    raise ValueError("prime bit size must lie between 29 and 61")

child = json.loads(CHILD.read_text())
marking = json.loads(MARKING.read_text())
ring = PolynomialRing(QQ, "T")
section = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
nx = polynomial(ring, section["x_numerator_coefficients_low_to_high"])
dx = polynomial(ring, section["x_denominator_coefficients_low_to_high"])
ny = polynomial(ring, section["y_numerator_coefficients_low_to_high"])
dy = polynomial(ring, section["y_denominator_coefficients_low_to_high"])
h = monic_power_root(dx, 2)
assert nx.degree() == 96 and nx.gcd(ny) in QQ and nx.gcd(h*dy) in QQ
nx_data, ny_data, dy_data, h_data = map(reduction_data, (nx, ny, dy, h))


def modular_correction(prime):
    finite = GF(prime)
    finite_ring = PolynomialRing(finite, "T")
    nx_p = reduce_polynomial(nx_data, finite_ring, prime)
    ny_p = reduce_polynomial(ny_data, finite_ring, prime)
    dy_p = reduce_polynomial(dy_data, finite_ring, prime)
    h_p = reduce_polynomial(h_data, finite_ring, prime)
    if nx_p.degree() != 96 or nx_p.gcd(ny_p).degree() or nx_p.gcd(h_p*dy_p).degree():
        raise ArithmeticError("degree or coprimality drop")
    normalizer = (ny_p*(h_p*dy_p).inverse_mod(nx_p)).mod(nx_p)
    if (normalizer*h_p*dy_p-ny_p) % nx_p:
        raise ArithmeticError("modular congruence failure")
    if normalizer.degree() != 95:
        raise ArithmeticError("normalizer leading-degree drop")
    return normalizer


def representative_primes(primes, count):
    """Select deterministic, evenly spread checks from incorporated primes.

    Rechecking every residue costs a fresh modular extended Euclidean
    computation per prime.  It is redundant once the final rational
    congruence is proved exactly, and makes the first true CRT candidate
    unnecessarily slow.  The end points are always included when possible.
    """

    if count >= len(primes):
        return tuple(primes)
    if count == 1:
        return (primes[-1],)
    indices = {
        index * (len(primes) - 1) // (count - 1)
        for index in range(count)
    }
    return tuple(primes[index] for index in sorted(indices))


residues = [ZZ.zero() for _ in range(96)]
modulus = ZZ.one()
accepted = []
rejected = []
candidate = None
prime = next_prime(ZZ(2)**(args.prime_bits-1)+7654321)
if args.resume_from:
    checkpoint = json.loads(args.resume_from.read_text())
    if checkpoint.get("schema") != "elkies-k3.h92-q6-child-q8-q-pole-crt-checkpoint.v1":
        raise ValueError("unrecognized CRT checkpoint schema")
    expected_inputs = {"child": digest(CHILD), "marking": digest(MARKING)}
    if checkpoint.get("inputs") != expected_inputs:
        raise ValueError("CRT checkpoint has different pinned inputs")
    if checkpoint.get("prime_bits") != args.prime_bits:
        raise ValueError("CRT checkpoint has different prime bit size")
    residues = [ZZ(value) for value in checkpoint["residues"]]
    modulus = ZZ(checkpoint["modulus"])
    accepted = [ZZ(value) for value in checkpoint["accepted_primes"]]
    rejected = checkpoint["rejected_primes"]
    prime = ZZ(checkpoint["next_prime"])
    if len(residues) != 96 or modulus.nbits() < len(accepted)*(args.prime_bits-1):
        raise ValueError("malformed CRT checkpoint")
while len(accepted) < args.maximum_primes:
    try:
        correction = modular_correction(prime)
    except (ArithmeticError, ZeroDivisionError) as error:
        rejected.append({"prime": int(prime), "reason": str(error)})
        prime = next_prime(prime+2)
        continue
    inverse_old_modulus = (modulus % prime).inverse_mod(prime)
    for index in range(96):
        value = ZZ(correction[index])
        adjustment = ((value-residues[index]) % prime*inverse_old_modulus) % prime
        residues[index] += modulus*adjustment
    modulus *= prime
    accepted.append(ZZ(prime))
    count = len(accepted)
    if count >= args.minimum_primes and count % args.reconstruct_every == 0:
        try:
            values = [QQ(value.rational_reconstruction(modulus)) for value in residues]
        except (ArithmeticError, ValueError):
            print(
                "H92Q6CHILDQCRT|accepted={}|modulus_bits={}|status=INSUFFICIENT_MODULUS".format(
                    count, modulus.nbits()
                ), flush=True,
            )
        else:
            possible = ring(values)
            if possible.degree() != 95:
                raise ArithmeticError("CRT candidate loses expected leading degree")
            if all(
                reduce_polynomial(reduction_data(possible), PolynomialRing(GF(p), "T"), p)
                == modular_correction(p)
                for p in representative_primes(accepted, args.accepted_validation_primes)
            ):
                candidate = possible
                print(
                    "H92Q6CHILDQCRT|accepted={}|modulus_bits={}|status=CANDIDATE".format(
                        count, modulus.nbits()
                    ), flush=True,
                )
                break
    prime = next_prime(prime+2)

if args.checkpoint:
    checkpoint_payload = {
        "schema": "elkies-k3.h92-q6-child-q8-q-pole-crt-checkpoint.v1",
        "inputs": {"child": digest(CHILD), "marking": digest(MARKING)},
        "prime_bits": args.prime_bits,
        "residues": [str(value) for value in residues],
        "modulus": str(modulus),
        "accepted_primes": [str(value) for value in accepted],
        "rejected_primes": rejected,
        "next_prime": str(prime),
    }
    args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    args.checkpoint.write_text(json.dumps(checkpoint_payload, indent=2, sort_keys=True)+"\n")

withheld = []
while len(withheld) < args.withheld_primes:
    try:
        correction = modular_correction(prime)
    except (ArithmeticError, ZeroDivisionError):
        prime = next_prime(prime+2)
        continue
    withheld.append((ZZ(prime), correction))
    prime = next_prime(prime+2)

complete = False
exact_check = False
largest_bits = None
if candidate is not None:
    for validation_prime, expected in withheld:
        finite_ring = PolynomialRing(GF(validation_prime), "T")
        if reduce_polynomial(reduction_data(candidate), finite_ring, validation_prime) != expected:
            raise ArithmeticError("candidate fails withheld modular validation")
    exact_check = ((candidate*h*dy-ny) % nx == 0)
    if not exact_check:
        raise ArithmeticError("candidate fails exact congruence")
    complete = True
    largest_bits = max(
        max(abs(value.numerator()).nbits(), value.denominator().nbits())
        for value in candidate.coefficients()
    )

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-q-pole-normalization-crt.v1",
    "complete": complete,
    "status": (
        "PASS_EXACT_CRT_PRINCIPAL_PART_NORMALIZATION" if complete
        else "INCOMPLETE_CRT_MODULAR_DATA"
    ),
    "inputs": {"child": digest(CHILD), "marking": digest(MARKING)},
    "modular_reconstruction": {
        "accepted_primes": [str(value) for value in accepted],
        "prime_bits": args.prime_bits,
        "rejected_primes": rejected,
        "withheld_primes": [str(value) for value, _ in withheld],
        "accepted_validation_primes": int(args.accepted_validation_primes),
        "modulus_bits": int(modulus.nbits()),
        "expected_degree_R": int(95),
        "largest_reconstructed_numden_bits": largest_bits,
    },
    "exact_check": exact_check,
    "normalizer": (
        {"R_coefficients_low_to_high": [str(value) for value in candidate.list()]}
        if complete else None
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q6CHILDQCRT|accepted={}|rejected={}|modulus_bits={}|complete={}|"
    "exact_check={}|status={}".format(
        len(accepted), len(rejected), modulus.nbits(), int(complete), int(exact_check),
        payload["status"],
    ), flush=True,
)
