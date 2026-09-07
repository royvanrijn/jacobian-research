#!/usr/bin/env python3
"""Evaluate candidate cubic-field squareclasses in a stored BNF-free map.

The relation collector emits global generators in ascending power-basis
coordinates.  This Sage/PARI bridge evaluates those generators at exactly the
local and auxiliary coordinates recorded in a signature-map artifact, then
writes a manifest directly consumable by ``residual_selmer_quotient.py``.

It deliberately performs no class-group calculation and makes no assertion
that a supplied candidate satisfies Selmer local conditions.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from sage.all import AA, PolynomialRing, QQ, pari


PROTOCOL = "BNFFREESIG"


def rational_text(value: object) -> Fraction:
    if isinstance(value, (int, str)):
        return Fraction(value)
    raise ValueError("power-basis coefficients must be integer or rational strings")


def qpari(value: Fraction):
    return pari(value.numerator) / pari(value.denominator)


def bit_mask(bits: list[int]) -> str:
    return f"0x{sum(int(bit) << index for index, bit in enumerate(bits)):x}"


def power_element(theta, coefficients: list[Fraction]):
    if len(coefficients) != 3:
        raise ValueError("a cubic-field generator needs exactly three coefficients")
    return sum(
        (qpari(coefficient) * theta**index for index, coefficient in enumerate(coefficients)),
        pari(0),
    )


def defining_polynomial(coefficients: list[Fraction]) -> str:
    if len(coefficients) != 4 or coefficients[-1] != 1:
        raise ValueError("expected a monic cubic defining polynomial")
    terms = []
    for degree, coefficient in enumerate(coefficients):
        if not coefficient:
            continue
        if degree == 0:
            terms.append(f"({coefficient})")
        elif degree == 1:
            terms.append(f"({coefficient})*z")
        else:
            terms.append(f"({coefficient})*z^{degree}")
    return "+".join(terms)


def exact_real_roots(coefficients: list[Fraction]):
    """Return increasing real roots as exact algebraic reals.

    A fixed-precision PARI root is not safe here: relation products can have
    large cancelling coefficients, and a false sign creates a fake residual
    signature direction.  Algebraic-real comparisons refine isolating
    intervals as needed.
    """
    ring = PolynomialRing(QQ, "z")
    z = ring.gen()
    polynomial = sum(QQ(value.numerator) / value.denominator * z**degree
                     for degree, value in enumerate(coefficients))
    return list(polynomial.roots(AA, multiplicities=False))


def real_sign_bit(coefficients: list[Fraction], root) -> int:
    value = sum(
        QQ(coefficient.numerator) / coefficient.denominator * root**index
        for index, coefficient in enumerate(coefficients)
    )
    if value == 0:
        raise ArithmeticError("candidate vanishes at a real cubic embedding")
    return 1 if value < 0 else 0


def local_places(nf, rational_prime: int):
    result = []
    for prime_ideal in pari.idealprimedec(nf, rational_prime):
        uniformizer = pari.nfbasistoalg(nf, pari.idealappr(nf, prime_ideal))
        if int(pari.idealval(nf, uniformizer, prime_ideal)) != 1:
            raise ArithmeticError(f"bad uniformizer above {rational_prime}")
        result.append((prime_ideal, uniformizer, pari.nfmodprinit(nf, prime_ideal)))
    return result


def bit_at_finite_place(nf, alpha, place, kind: str) -> int:
    prime_ideal, uniformizer, mod_prime = place
    valuation = int(pari.idealval(nf, alpha, prime_ideal))
    if kind.endswith("valuation_parity"):
        return valuation & 1
    if kind.endswith("unit_squareclass"):
        unit = alpha / uniformizer**valuation
        residue = pari.nfmodpr(nf, unit, mod_prime)
        return 0 if bool(pari.issquare(residue)) else 1
    raise ValueError(f"unexpected finite coordinate kind {kind}")


def two_adic_coordinates(nf, two_primes, alpha, basis):
    for mask in range(1 << len(basis)):
        candidate = alpha
        bits = []
        for index, generator in enumerate(basis):
            bit = (mask >> index) & 1
            bits.append(bit)
            if bit:
                candidate /= generator
        if all(bool(pari.nfislocalpower(nf, prime, candidate, 2)) for prime in two_primes):
            return bits
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--signature-map", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument(
        "--class-quotient-audit",
        type=Path,
        help="attach a matching BNF-free S-class quotient audit to the output manifest",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    signature = json.loads(args.signature_map.read_text())
    candidates_input = json.loads(args.candidates.read_text())
    if signature.get("schema") != "elliptic-curves.bnf-free-signature-map.v1":
        raise ValueError("unexpected signature-map schema")
    if isinstance(candidates_input, dict):
        if candidates_input.get("schema") not in {
            "elliptic-curves.bnf-free-squareclass-candidates.v1",
            "elliptic-curves.bnf-free-norm-filtered-squareclass-candidates.v1",
        }:
            raise ValueError("unexpected candidate-object schema")
        candidates = candidates_input.get("candidates")
    else:
        candidates = candidates_input
    if not isinstance(candidates, list):
        raise ValueError("candidate input must be a JSON list")

    polynomial_coefficients = [
        rational_text(value) for value in signature["defining_polynomial_ascending"]
    ]
    class_quotient_certification = signature["class_quotient_certification"]
    if args.class_quotient_audit:
        class_audit = json.loads(args.class_quotient_audit.read_text())
        if class_audit.get("protocol") != "BNFFREECLASS-v1":
            raise ValueError("unexpected BNF-free class-quotient audit")
        if class_audit.get("defining_polynomial_ascending") != signature.get(
            "defining_polynomial_ascending"
        ):
            raise ValueError("class-quotient audit is for a different cubic field")
        if not class_audit.get("principal_relations_verified"):
            raise ValueError("class-quotient audit did not verify its principal relations")
        class_quotient_certification = class_audit[
            "class_quotient_certification"
        ]
    polynomial = defining_polynomial(polynomial_coefficients)
    nf = pari.nfinit(pari(polynomial))
    theta = pari(f"Mod(z,{polynomial})")

    local_coordinates = signature["local_coordinates"]
    fingerprint_coordinates = signature["fingerprint_coordinates"]
    finite_cache = {}

    def place_for(coordinate):
        rational_prime = int(coordinate["rational_prime"])
        if rational_prime not in finite_cache:
            finite_cache[rational_prime] = local_places(nf, rational_prime)
        place = finite_cache[rational_prime][int(coordinate["place_index"])]
        if str(place[0]) != coordinate["prime_ideal"]:
            raise ArithmeticError("signature map no longer has the recorded prime-ideal order")
        return place

    two_coordinates = [
        coordinate
        for coordinate in local_coordinates
        if coordinate["kind"]
        in {"two_adic_product_basis", "two_adic_product_basis_extension"}
    ]
    two_primes = list(pari.idealprimedec(nf, 2))
    two_basis = [
        power_element(
            theta,
            [rational_text(value) for value in coordinate["generator_power_basis"].strip("()").split(",")],
        )
        for coordinate in two_coordinates
    ]
    roots = exact_real_roots(polynomial_coefficients)

    original_two_dimension = len(two_basis)
    dynamic_two_coordinates = []
    evaluated = []
    for candidate in candidates:
        if not isinstance(candidate, dict):
            raise ValueError("every candidate must be a JSON object")
        coefficients = [rational_text(value) for value in candidate["generator_coefficients"]]
        alpha = power_element(theta, coefficients)
        two_bits = two_adic_coordinates(nf, two_primes, alpha, two_basis)
        if two_bits is None:
            # The stored known-MW image does not span this genuinely new local
            # 2-adic direction.  Extend the coordinate system at its end;
            # every stored known image has zero in the new basis coordinate.
            for prior in evaluated:
                prior["dynamic_two_bits"].append(0)
            two_basis.append(alpha)
            dynamic_two_coordinates.append(
                {
                    "kind": "two_adic_product_basis_extension",
                    "basis_index": len(two_basis) - 1,
                    "generator_power_basis": "("
                    + ",".join(str(value) for value in coefficients)
                    + ")",
                    "two_adic_primes": [str(pr) for pr in two_primes],
                }
            )
            two_bits = [0] * len(two_basis)
            two_bits[-1] = 1
        two_index = 0
        local_bits = []
        for coordinate in local_coordinates:
            kind = coordinate["kind"]
            if kind in {
                "two_adic_product_basis",
                "two_adic_product_basis_extension",
            }:
                local_bits.append(two_bits[two_index])
                two_index += 1
            elif kind == "real_sign":
                local_bits.append(
                    real_sign_bit(coefficients, roots[int(coordinate["embedding_index"])])
                )
            else:
                local_bits.append(bit_at_finite_place(nf, alpha, place_for(coordinate), kind))
        fingerprint_bits = [
            bit_at_finite_place(nf, alpha, place_for(coordinate), coordinate["kind"])
            for coordinate in fingerprint_coordinates
        ]
        evaluated.append(
            {
                "label": str(candidate["label"]),
                "generator": "(" + ",".join(str(value) for value in coefficients) + ")",
                "base_local_bits": local_bits,
                "dynamic_two_bits": two_bits[original_two_dimension:],
                "fingerprint_bits": fingerprint_bits,
            }
        )

    images = [
        {
            "label": candidate["label"],
            "generator": candidate["generator"],
            "local": bit_mask(candidate["base_local_bits"] + candidate["dynamic_two_bits"]),
            "fingerprint": bit_mask(candidate["fingerprint_bits"]),
        }
        for candidate in evaluated
    ]
    output_local_coordinates = local_coordinates + dynamic_two_coordinates

    output = {
        "local_dimension": len(output_local_coordinates),
        "fingerprint_dimension": signature["fingerprint_dimension"],
        "local_coordinates": output_local_coordinates,
        "known_mw_images": signature["known_mw_images"],
        "candidate_images": images,
        "class_quotient_certification": class_quotient_certification,
        "status": "images_computed_not_a_selmer_local_solubility_certificate",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|stage=complete|candidates={len(images)}"
        f"|local_dimension={len(output_local_coordinates)}"
        f"|fingerprint_dimension={signature['fingerprint_dimension']}"
        f"|status=IMAGES_ONLY_NOT_SELMER",
        flush=True,
    )


if __name__ == "__main__":
    main()
