#!/usr/bin/env sage
"""Build an exact BNF-free Kummer signature for a q12o5867 fibre.

For a generalized integral model

    y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6

put ``X=4*x`` and ``Z=4*(2*y+a1*x+a3)``.  Then

    Z^2 = X^3 + b2*X^2 + 8*b4*X + 16*b6,

where ``b2=a1^2+4*a2``, ``b4=a1*a3+2*a4`` and
``b6=a3^2+4*a6``.  Thus the no-rational-2-torsion Kummer element of a point
is exactly ``X-theta`` in the monic cubic field.

This program records the local images of the certified seventeen points at
2, selected odd places, and the real places.  It then greedily adds auxiliary
good-prime fingerprints until their combined image has rank 17.  Auxiliary
coordinates are witnesses, not Selmer conditions; no class-group or Selmer
upper bound is claimed.
"""

from __future__ import annotations

from research_runtime.pari_context import prepared_prime_ideals, prepared_factor

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time
from typing import Sequence


REPOSITORY = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = REPOSITORY / "elliptic-curves"
CAS_ROOT = ELLIPTIC_ROOT / "cas"
sys.path.insert(0, str(CAS_ROOT))

from elliptic_candidate_record import is_on_weierstrass_curve  # noqa: E402
from run_fermigier_rank20_auxiliary_fingerprints import (  # noqa: E402
    append_columns,
    f2_mask,
    f2_rank,
    prime_local_rows,
    qpari,
    two_adic_coords,
)


PROTOCOL = "Q12BNFSIG"
INPUT_STATUS = "PASS_EXACT_Q12O5867_SPECIALIZED_GENERIC_RANK17_LOWER_BOUND"
KNOWN_RANK = 17


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def monic_cubic_coefficients(
    model: Sequence[int | Fraction | str],
) -> tuple[int, int, int, int]:
    """Return ascending coefficients of the integral monic descent cubic."""

    if len(model) != 5:
        raise ValueError("a generalized Weierstrass model needs five coefficients")
    values = tuple(Fraction(value) for value in model)
    if any(value.denominator != 1 for value in values):
        raise ValueError("the model must be integral")
    a1, a2, a3, a4, a6 = (value.numerator for value in values)
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    return 16 * b6, 8 * b4, b2, 1


def point_on_monic_cubic(
    model: Sequence[int | Fraction | str],
    point: Sequence[int | Fraction | str],
) -> tuple[Fraction, Fraction]:
    """Transport one affine generalized-Weierstrass point to ``Z^2=f(X)``."""

    if len(model) != 5 or len(point) != 2:
        raise ValueError("unexpected model or point length")
    a1, _a2, a3, _a4, _a6 = (Fraction(value) for value in model)
    x_coordinate, y_coordinate = (Fraction(value) for value in point)
    return (
        4 * x_coordinate,
        4 * (2 * y_coordinate + a1 * x_coordinate + a3),
    )


def evaluate_cubic(
    coefficients: Sequence[int | Fraction], value: Fraction
) -> Fraction:
    return sum(Fraction(coefficient) * value**index for index, coefficient in enumerate(coefficients))


def parse_primes(text: str) -> tuple[int, ...]:
    try:
        values = tuple(sorted({int(item) for item in text.split(",") if item}))
    except ValueError as error:
        raise argparse.ArgumentTypeError("primes must be comma-separated integers") from error
    if any(value < 3 for value in values):
        raise argparse.ArgumentTypeError("odd local primes must be at least 3")
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--specialization", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--prime-bound", type=int, default=5000)
    parser.add_argument(
        "--odd-local-primes",
        type=parse_primes,
        default=(),
        help="known odd Selmer-relevant primes to include as local coordinates",
    )
    args = parser.parse_args()
    if args.prime_bound < 3:
        parser.error("--prime-bound must be at least 3")
    sys.set_int_max_str_digits(0)

    from sage.all import AA, GF, PolynomialRing, QQ, ZZ, pari, prime_range

    artifact = json.loads(args.specialization.read_text())
    if artifact.get("status") != INPUT_STATUS:
        raise ValueError("input is not an exact q12o5867 rank-17 specialization")
    minimal = artifact["global_minimal_specialization"]
    model = tuple(Fraction(value) for value in minimal["model"])
    if any(value.denominator != 1 for value in model):
        raise ValueError("global minimal model is not integral")
    points = tuple(
        (Fraction(record[0]), Fraction(record[1])) for record in minimal["points"]
    )
    if len(points) != KNOWN_RANK or len(set(points)) != KNOWN_RANK:
        raise ValueError("expected seventeen distinct certified points")
    if any(not is_on_weierstrass_curve(model, point) for point in points):
        raise ArithmeticError("a certified point misses the input minimal model")

    coefficients = monic_cubic_coefficients(model)
    transformed = tuple(point_on_monic_cubic(model, point) for point in points)
    if any(z_coordinate**2 != evaluate_cubic(coefficients, x_coordinate) for x_coordinate, z_coordinate in transformed):
        raise ArithmeticError("the exact monic-cubic point transport failed")

    ring = PolynomialRing(QQ, "z")
    z = ring.gen()
    polynomial = sum(ZZ(value) * z**index for index, value in enumerate(coefficients))
    if not polynomial.is_irreducible():
        raise ValueError("the specialization has rational 2-torsion; cubic-field descent is inapplicable")
    discriminant = ZZ(polynomial.discriminant())
    polynomial_text = "+".join(
        f"({value})*t^{index}" for index, value in enumerate(coefficients) if value
    )
    nf_started = time.monotonic()
    from research_runtime.pari_context import prepared_nf
    nf = prepared_nf(pari(polynomial_text))
    theta = pari(f"Mod(t,{polynomial_text})")
    nf_seconds = time.monotonic() - nf_started
    x_coordinates = [point[0] for point in transformed]
    alphas = [qpari(pari, x_coordinate) - theta for x_coordinate in x_coordinates]

    rows = [[] for _ in alphas]
    local_coordinates = []
    for rational_prime in args.odd_local_primes:
        if discriminant % rational_prime == 0 or rational_prime == 2:
            extra, places = prime_local_rows(pari, nf, alphas, rational_prime)
            rows = append_columns(rows, extra)
            for place_index, prime_ideal in enumerate(places):
                local_coordinates.extend(
                    (
                        {
                            "kind": "odd_valuation_parity",
                            "rational_prime": rational_prime,
                            "prime_ideal": prime_ideal,
                            "place_index": place_index,
                        },
                        {
                            "kind": "odd_unit_squareclass",
                            "rational_prime": rational_prime,
                            "prime_ideal": prime_ideal,
                            "place_index": place_index,
                        },
                    )
                )
        else:
            raise ValueError(
                f"declared odd local prime {rational_prime} does not divide the cubic discriminant"
            )

    two_primes = list(prepared_prime_ideals(nf, 2))
    two_basis, two_basis_origins, two_rows = two_adic_coords(
        pari, nf, two_primes, alphas
    )
    rows = append_columns(rows, two_rows)
    local_coordinates.extend(
        {
            "kind": "two_adic_product_basis",
            "basis_index": index,
            "generator_power_basis": f"({x_coordinates[two_basis_origins[index]]},-1,0)",
            "two_adic_primes": [str(prime) for prime in two_primes],
        }
        for index in range(len(two_basis))
    )

    exact_roots = list(polynomial.roots(AA, multiplicities=False))
    real_rows = [
        [1 if QQ(x_coordinate) - root < 0 else 0 for root in exact_roots]
        for x_coordinate in x_coordinates
    ]
    rows = append_columns(rows, real_rows)
    local_coordinates.extend(
        {
            "kind": "real_sign",
            "embedding_index": index,
            "root_order": "increasing_real_root",
        }
        for index in range(len(exact_roots))
    )
    local_rows = [list(row) for row in rows]
    baseline_rank = f2_rank(rows)

    selected = []
    fingerprint_coordinates = []
    current = rows
    excluded = {2, *args.odd_local_primes}
    for rational_prime in prime_range(3, args.prime_bound + 1):
        rational_prime = int(rational_prime)
        if rational_prime in excluded or discriminant % rational_prime == 0:
            continue
        extra, places = prime_local_rows(pari, nf, alphas, rational_prime)
        trial = append_columns(current, extra)
        old_rank = f2_rank(current)
        new_rank = f2_rank(trial)
        if new_rank > old_rank:
            selected.append((rational_prime, len(places), new_rank - old_rank))
            current = trial
            for place_index, prime_ideal in enumerate(places):
                fingerprint_coordinates.extend(
                    (
                        {
                            "kind": "auxiliary_valuation_parity",
                            "rational_prime": rational_prime,
                            "prime_ideal": prime_ideal,
                            "place_index": place_index,
                        },
                        {
                            "kind": "auxiliary_unit_squareclass",
                            "rational_prime": rational_prime,
                            "prime_ideal": prime_ideal,
                            "place_index": place_index,
                        },
                    )
                )
        if f2_rank(current) == KNOWN_RANK:
            break
    final_rank = f2_rank(current)
    if final_rank != KNOWN_RANK:
        raise RuntimeError(
            f"fingerprints through {args.prime_bound} have rank {final_rank}, not 17"
        )

    fingerprint_width = len(fingerprint_coordinates)
    output = {
        "schema": "elliptic-curves.bnf-free-signature-map.v1",
        "status": "known_kummer_image_only_not_a_selmer_bound",
        "source": {
            "kind": "q12o5867_exact_specialization",
            "path": str(args.specialization.resolve()),
            "sha256": sha256_file(args.specialization),
            "parameter": artifact["parameter"],
            "certified_rank_lower_bound": KNOWN_RANK,
        },
        "generalized_minimal_model": [str(value) for value in model],
        "monic_cubic_point_map": {
            "X": "4*x",
            "Z": "4*(2*y+a1*x+a3)",
        },
        "field_generator": "theta",
        "generator_coordinate_order": ["1", "theta", "theta^2"],
        "defining_polynomial_ascending": [str(value) for value in coefficients],
        "defining_polynomial_discriminant": str(discriminant),
        "nfinit_seconds": nf_seconds,
        "local_dimension": len(local_coordinates),
        "fingerprint_dimension": fingerprint_width,
        "local_coordinates": local_coordinates,
        "fingerprint_coordinates": fingerprint_coordinates,
        "known_mw_images": [
            {
                "label": f"P{index + 1}",
                "generator": f"({x_coordinates[index]},-1,0)",
                "generator_coefficients": [str(x_coordinates[index]), "-1", "0"],
                "local": f"0x{f2_mask(local_rows[index]):x}",
                "fingerprint": f"0x{f2_mask(current[index][len(local_coordinates):]):x}",
            }
            for index in range(KNOWN_RANK)
        ],
        "selected_auxiliary_primes": [prime for prime, _, _ in selected],
        "known_mw_local_rank": baseline_rank,
        "known_mw_target_rank": final_rank,
        "class_quotient_certification": {
            "method": "none",
            "remaining_dimension_upper_bound": None,
        },
        "claim_boundary": [
            "The point transport, Kummer generators, local images, and auxiliary fingerprints are exact.",
            "Auxiliary good-prime coordinates distinguish squareclasses but are not Selmer conditions.",
            "The selected odd local-prime list need not contain every bad prime.",
            "No class-group completeness, Selmer upper bound, or rank upper bound is claimed.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|stage=complete|known_rank={final_rank}"
        f"|local_dimension={len(local_coordinates)}"
        f"|fingerprint_dimension={fingerprint_width}"
        f"|two_adic_dimension={len(two_basis)}"
        f"|auxiliary_primes={output['selected_auxiliary_primes']}"
        f"|nfinit_seconds={nf_seconds:.6f}|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
