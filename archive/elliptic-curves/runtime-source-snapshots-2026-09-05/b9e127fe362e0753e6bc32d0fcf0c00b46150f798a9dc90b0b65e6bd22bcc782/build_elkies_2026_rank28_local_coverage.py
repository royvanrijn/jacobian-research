#!/usr/bin/env sage
"""Build the exact bad-place signature fixture for Elkies rank 28.

The published rank-28 fibre is used as a positive control for the local layer
of cubic 2-descent.  This program reconstructs the generic seventeen points
and the certified eleven-point public complement, recomputes their Kummer
images at every bad finite place and at infinity, and checks the generic rows
against the independently pinned bad-place ledger.

The output deliberately keeps the generic seventeen as ``known_mw_images``
for the residual quotient.  The public complement is recorded separately as
a calibration: on this fibre it adds no bad-place local-image dimension at
any block, despite adding eleven globally independent Mordell--Weil
directions.  This proves why these signatures cannot be promoted to a Selmer
upper bound.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = ROOT / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
BAD_PLACE_LEDGER = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_rank28_bad_place_kummer_ledger_v1.json"
)
DEFAULT_SIGNATURE_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_rank28_generic17_local_signature_v1.json"
)
PROTOCOL = "ELKIESR28COVERAGE"
SIGNATURE_SCHEMA = "elliptic-curves.bnf-free-signature-map.v1"

sys.path[:0] = [str(ELLIPTIC_ROOT), str(CAS)]

from build_elkies_2026_rank28_bad_place_ledger import (  # noqa: E402
    DISCRIMINANT_FACTORIZATION,
)
from build_elkies_2026_rank28_relative_descent_magma import (  # noqa: E402
    load_relative_input,
)
from build_q12o5867_bnf_free_signature import (  # noqa: E402
    evaluate_cubic,
    monic_cubic_coefficients,
    point_on_monic_cubic,
)
from run_fermigier_rank20_auxiliary_fingerprints import (  # noqa: E402
    f2_rank,
    prime_local_rows,
    qpari,
    two_adic_coords,
)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def mask(row: list[int]) -> str:
    return hex(sum((int(value) & 1) << index for index, value in enumerate(row)))


def append_rows(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    if len(left) != len(right):
        raise ArithmeticError("local row blocks have different point counts")
    return [left[index] + right[index] for index in range(len(left))]


def block_index(ledger: dict) -> dict[tuple[str, int | None], dict]:
    answer = {}
    for block in ledger["completed_blocks"]:
        kind = block.get("kind")
        if kind == "factorization_proof":
            continue
        prime = int(block["rational_prime"]) if "rational_prime" in block else None
        key = (str(kind), prime)
        if key in answer:
            raise ValueError(f"duplicate bad-place block {key}")
        answer[key] = block
    return answer


def build_signature() -> dict:
    from sage.all import AA, PolynomialRing, QQ, ZZ, pari

    ledger = json.loads(BAD_PLACE_LEDGER.read_text())
    if not (
        ledger.get("status")
        == "COMPLETE_ALL_BAD_PLACE_KUMMER_IMAGES_NOT_A_SELMER_BOUND"
        and ledger.get("parameter") == "-9529/5471"
        and ledger.get("factorization_product_verified") is True
        and ledger.get("factor_primality_proof_completed") is True
        and ledger.get("all_bad_place_blocks_completed") is True
        and ledger.get("completed_local_block_count") == 13
    ):
        raise ValueError("the exact rank-28 bad-place ledger is stale")

    source = load_relative_input()
    points = source.generic_points + source.public_complement
    if len(source.generic_points) != 17 or len(source.public_complement) != 11:
        raise ArithmeticError("the rank-28 control no longer has a 17+11 decomposition")
    model = tuple(Fraction(value) for value in source.model)
    if [str(value) for value in model] != ledger["specialization"]["global_minimal_model"]:
        raise ArithmeticError("the rank-28 sources use different minimal models")

    coefficients = monic_cubic_coefficients(model)
    if [str(value) for value in coefficients] != ledger[
        "descent_cubic_coefficients_ascending"
    ]:
        raise ArithmeticError("the rank-28 descent cubic changed")
    transformed = tuple(point_on_monic_cubic(model, point) for point in points)
    if any(z * z != evaluate_cubic(coefficients, x) for x, z in transformed):
        raise ArithmeticError("a rank-28 point failed monic-cubic transport")

    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    polynomial = sum(ZZ(value) * x**index for index, value in enumerate(coefficients))
    factor_primes = [prime for prime, _ in DISCRIMINANT_FACTORIZATION]
    pari.addprimes(factor_primes)
    nf = pari.nfinit([pari(polynomial), factor_primes])
    if list(pari.nfcertify(nf)):
        raise ArithmeticError("factor-supplied number field failed certification")
    theta = pari(f"Mod(x,{polynomial})")
    alphas = [qpari(pari, x_coordinate) - theta for x_coordinate, _ in transformed]

    pinned = block_index(ledger)
    all_rows = [[] for _ in points]
    local_coordinates = []
    block_summaries = []
    odd_primes = [
        prime for prime, _ in DISCRIMINANT_FACTORIZATION if prime != 2
    ]
    for rational_prime in odd_primes:
        rows, prime_ideals = prime_local_rows(
            pari, nf, alphas, rational_prime
        )
        expected = pinned[("odd_bad_prime", rational_prime)]
        if rows[:17] != expected["rows_for_generic_P1_through_P17"]:
            raise ArithmeticError(f"generic local rows changed at p={rational_prime}")
        if [str(prime) for prime in prime_ideals] != expected["prime_ideals"]:
            raise ArithmeticError(f"prime-ideal order changed at p={rational_prime}")
        all_rows = append_rows(all_rows, rows)
        for place_index, prime_ideal in enumerate(prime_ideals):
            local_coordinates.extend(
                (
                    {
                        "kind": "odd_valuation_parity",
                        "rational_prime": rational_prime,
                        "place_index": place_index,
                        "prime_ideal": str(prime_ideal),
                    },
                    {
                        "kind": "odd_unit_squareclass",
                        "rational_prime": rational_prime,
                        "place_index": place_index,
                        "prime_ideal": str(prime_ideal),
                    },
                )
            )
        generic_rank = f2_rank(rows[:17])
        combined_rank = f2_rank(rows)
        block_summaries.append(
            {
                "kind": "odd_bad_prime",
                "rational_prime": str(rational_prime),
                "coordinate_dimension": len(rows[0]),
                "generic_17_image_rank": generic_rank,
                "generic_plus_public_11_image_rank": combined_rank,
                "public_complement_incremental_rank": combined_rank - generic_rank,
            }
        )

    two_primes = list(pari.idealprimedec(nf, 2))
    two_basis, two_origins, two_rows = two_adic_coords(
        pari, nf, two_primes, alphas
    )
    expected_two = pinned[("two_adic", 2)]
    if two_rows[:17] != expected_two["rows_for_generic_P1_through_P17"]:
        raise ArithmeticError("generic two-adic rows changed")
    if [str(prime) for prime in two_primes] != expected_two["prime_ideals"]:
        raise ArithmeticError("two-adic prime-ideal order changed")
    all_rows = append_rows(all_rows, two_rows)
    local_coordinates.extend(
        {
            "kind": "two_adic_product_basis",
            "basis_index": index,
            "basis_origin_one_based": int(two_origins[index]) + 1,
            "generator": str(two_basis[index]),
            "two_adic_primes": [str(prime) for prime in two_primes],
        }
        for index in range(len(two_basis))
    )
    generic_two_rank = f2_rank(two_rows[:17])
    combined_two_rank = f2_rank(two_rows)
    block_summaries.append(
        {
            "kind": "two_adic",
            "rational_prime": "2",
            "coordinate_dimension": len(two_rows[0]),
            "generic_17_image_rank": generic_two_rank,
            "generic_plus_public_11_image_rank": combined_two_rank,
            "public_complement_incremental_rank": (
                combined_two_rank - generic_two_rank
            ),
            "coverage_classification": (
                "UNRESOLVED_TWO_ADIC_LOCAL_KUMMER_IMAGE_COVERAGE"
            ),
        }
    )

    roots = list(polynomial.roots(AA, multiplicities=False))
    real_rows = [
        [1 if QQ(x_coordinate) - root < 0 else 0 for root in roots]
        for x_coordinate, _ in transformed
    ]
    expected_real = pinned[("real_places", None)]
    if real_rows[:17] != expected_real["rows_for_generic_P1_through_P17"]:
        raise ArithmeticError("generic real Kummer rows changed")
    all_rows = append_rows(all_rows, real_rows)
    local_coordinates.extend(
        {
            "kind": "real_sign",
            "embedding_index": index,
            "root_order": "increasing_real_root",
        }
        for index in range(len(roots))
    )
    generic_real_rank = f2_rank(real_rows[:17])
    combined_real_rank = f2_rank(real_rows)
    block_summaries.append(
        {
            "kind": "real_places",
            "coordinate_dimension": len(real_rows[0]),
            "generic_17_image_rank": generic_real_rank,
            "generic_plus_public_11_image_rank": combined_real_rank,
            "public_complement_incremental_rank": (
                combined_real_rank - generic_real_rank
            ),
        }
    )

    if f2_rank(all_rows[:17]) != ledger["combined_known_kummer_rank"]:
        raise ArithmeticError("combined generic local-image rank changed")
    combined_rank = f2_rank(all_rows)
    if combined_rank != ledger["combined_known_kummer_rank"]:
        raise ArithmeticError("public complement unexpectedly changes local-image rank")
    if any(row["public_complement_incremental_rank"] for row in block_summaries):
        raise ArithmeticError("public complement changes a bad-place block rank")

    local_dimension = len(local_coordinates)
    if any(len(row) != local_dimension for row in all_rows):
        raise ArithmeticError("local coordinate bookkeeping lost alignment")
    x_coordinates = [point[0] for point in transformed]
    output = {
        "schema": SIGNATURE_SCHEMA,
        "status": "KNOWN_RANK28_BAD_PLACE_IMAGES_NOT_A_SELMER_BOUND",
        "source": {
            "kind": "elkies_2026_rank28_positive_control",
            "parameter": "-9529/5471",
            "bad_place_ledger": {
                "path": str(BAD_PLACE_LEDGER.resolve()),
                "sha256": file_sha256(BAD_PLACE_LEDGER),
            },
            "positive_controls_sha256": source.controls_sha256,
            "generic_point_sequence_sha256": source.generic_point_sequence_sha256,
            "combined_point_sequence_sha256": source.combined_point_sequence_sha256,
        },
        "generalized_minimal_model": [str(value) for value in model],
        "monic_cubic_point_map": {
            "X": "4*x",
            "Z": "4*(2*y+a1*x+a3)",
        },
        "field_generator": "theta",
        "generator_coordinate_order": ["1", "theta", "theta^2"],
        "defining_polynomial_ascending": [str(value) for value in coefficients],
        "defining_polynomial_discriminant": ledger["descent_cubic_discriminant"],
        "local_dimension": local_dimension,
        "fingerprint_dimension": 0,
        "local_coordinates": local_coordinates,
        "fingerprint_coordinates": [],
        "known_mw_images": [
            {
                "label": f"P{index + 1}",
                "generator": f"({x_coordinates[index]},-1,0)",
                "generator_coefficients": [str(x_coordinates[index]), "-1", "0"],
                "local": mask(all_rows[index]),
                "fingerprint": "0x0",
            }
            for index in range(17)
        ],
        "public_positive_control_complement_images": [
            {
                "label": f"Q{index + 1}",
                "generator": f"({x_coordinates[index + 17]},-1,0)",
                "generator_coefficients": [
                    str(x_coordinates[index + 17]),
                    "-1",
                    "0",
                ],
                "local": mask(all_rows[index + 17]),
                "fingerprint": "0x0",
            }
            for index in range(11)
        ],
        "block_rank_comparison": block_summaries,
        "known_mw_local_rank": f2_rank(all_rows[:17]),
        "combined_rank28_local_rank": combined_rank,
        "public_complement_bad_place_incremental_rank": 0,
        "selected_auxiliary_primes": [],
        "class_quotient_certification": {
            "method": "none",
            "remaining_dimension_upper_bound": None,
        },
        "claim_boundary": [
            "Every local row is recomputed exactly and the generic seventeen are checked against the pinned bad-place ledger.",
            "The eleven globally independent public directions add zero bad-place local-image rank; local signatures do not measure the Mordell--Weil quotient gain.",
            "The known local-image span is not the ambient local Kummer image unless a separate coverage audit proves equality at that place.",
            "No class-group completeness, local-solubility completion, Selmer upper bound, or rank upper bound is claimed.",
        ],
    }
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_SIGNATURE_OUTPUT)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    output = build_signature()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if args.overwrite else "x"
    with args.output.open(mode) as handle:
        json.dump(output, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(
        f"{PROTOCOL}|stage=complete|local_dimension={output['local_dimension']}"
        f"|generic_rank={output['known_mw_local_rank']}"
        f"|combined_rank={output['combined_rank28_local_rank']}"
        f"|public_increment={output['public_complement_bad_place_incremental_rank']}"
        f"|output={args.output}"
    )


if __name__ == "__main__":
    main()
