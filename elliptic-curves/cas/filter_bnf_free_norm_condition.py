#!/usr/bin/env sage
"""Filter BNF-free cubic squareclass candidates by exact descent conditions.

For a monic short model ``y^2 = f(x)`` and ``K = Q(theta)``, with
``f(theta)=0``, the Kummer representative is ``x-theta`` and satisfies
``Norm(x-theta)=y^2``.  Thus every 2-descent class must have rational-square
norm.  The script also removes representatives that are exact global squares
in the cubic field.  With ``--generate-norm-kernel``, it first takes the GF(2)
kernel of rational norm parity over the declared Selmer primes, so it finds
norm-square *products* of relation representatives rather than testing only
each closed relation individually.  Both tests are BNF-free and must be
applied before treating a K(S,2) representative as a possible nonzero Selmer
class.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from sage.all import NumberField, PolynomialRing, QQ, ZZ


PROTOCOL = "BNFFREENORM"


def rational(value: str) -> QQ:
    value = Fraction(value)
    return QQ(value.numerator) / QQ(value.denominator)


def cubic_field(coefficients):
    coefficients = [rational(value) for value in coefficients]
    if len(coefficients) != 4 or coefficients[-1] != 1:
        raise ValueError("the candidate field must have a monic cubic polynomial")
    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    return NumberField(sum(value * x**index for index, value in enumerate(coefficients)), "theta")


def is_rational_square(value) -> bool:
    value = QQ(value)
    return value >= 0 and value.numerator().is_square() and value.denominator().is_square()


def norm_parity_mask(norm, selmer_primes: list[int]) -> int:
    """Squareclass of a rational norm, restricted to S and its sign."""
    norm = QQ(norm)
    numerator = abs(ZZ(norm.numerator()))
    denominator = ZZ(norm.denominator())
    mask = 0
    for index, prime in enumerate(selmer_primes):
        if (numerator.valuation(prime) - denominator.valuation(prime)) & 1:
            mask |= 1 << index
    if norm < 0:
        mask |= 1 << len(selmer_primes)
    return mask


def kernel_dependencies(masks: list[int]) -> list[int]:
    """Independent dependencies among packed GF(2) rows, with provenance."""
    pivots: dict[int, tuple[int, int]] = {}
    dependencies = []
    for index, mask in enumerate(masks):
        row = mask
        provenance = 1 << index
        while row:
            pivot = row.bit_length() - 1
            previous = pivots.get(pivot)
            if previous is None:
                pivots[pivot] = (row, provenance)
                break
            row ^= previous[0]
            provenance ^= previous[1]
        if row == 0:
            dependencies.append(provenance)
    return dependencies


def selmer_primes_from_ledger(path: Path, coefficients) -> list[int]:
    ledger = json.loads(path.read_text())
    if ledger.get("schema") != "elliptic-curves.bnf-free-principal-relation-ledger.v1":
        raise ValueError("expected a BNF-free principal-relation ledger")
    if ledger.get("defining_polynomial_ascending") != [str(value) for value in coefficients]:
        raise ValueError("relation ledger belongs to a different cubic field")
    primes = sorted({int(value) for value in ledger.get("selmer_rational_primes", ())})
    if not primes:
        raise ValueError("relation ledger has no declared Selmer rational primes")
    return primes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument(
        "--generate-norm-kernel",
        action="store_true",
        help="emit a basis of norm-square products of all supplied candidates",
    )
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--relation-ledger", type=Path)
    source.add_argument("--selmer-rational-primes")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    record = json.loads(args.candidates.read_text())
    if not isinstance(record, dict) or record.get("schema") != "elliptic-curves.bnf-free-squareclass-candidates.v1":
        raise ValueError("expected a BNF-free squareclass-candidate object")
    field = cubic_field(record["field_polynomial_ascending"])
    theta = field.gen()
    coefficients = [rational(value) for value in record["field_polynomial_ascending"]]
    source_candidates = []
    for candidate in record.get("candidates", []):
        if not isinstance(candidate, dict):
            raise ValueError("candidate list must contain objects")
        coordinates = [rational(value) for value in candidate["generator_coefficients"]]
        if len(coordinates) != 3:
            raise ValueError("a cubic candidate needs three power-basis coordinates")
        alpha = sum(value * theta**index for index, value in enumerate(coordinates))
        norm = QQ(alpha.norm())
        source_candidates.append((candidate, alpha, norm))

    if args.generate_norm_kernel:
        if args.relation_ledger:
            selmer_primes = selmer_primes_from_ledger(args.relation_ledger, coefficients)
        elif args.selmer_rational_primes:
            try:
                selmer_primes = sorted({int(value.strip()) for value in args.selmer_rational_primes.split(",") if value.strip()})
            except ValueError as exc:
                raise ValueError("selmer-rational-primes must be comma-separated integers") from exc
            if not selmer_primes or any(value < 2 or not ZZ(value).is_prime(proof=False) for value in selmer_primes):
                raise ValueError("selmer-rational-primes must contain rational primes")
        else:
            raise ValueError("--generate-norm-kernel requires --relation-ledger or --selmer-rational-primes")
        dependencies = kernel_dependencies(
            [norm_parity_mask(norm, selmer_primes) for _, _, norm in source_candidates]
        )
        candidates_for_filter = []
        for dependency_index, dependency in enumerate(dependencies):
            alpha = field(1)
            source_indices = []
            for index, (_, candidate_alpha, _) in enumerate(source_candidates):
                if (dependency >> index) & 1:
                    alpha *= candidate_alpha
                    source_indices.append(index)
            norm = QQ(alpha.norm())
            if not is_rational_square(norm):
                raise ArithmeticError("computed norm-kernel dependency does not have square norm")
            alpha_coefficients = list(alpha.polynomial().list())
            alpha_coefficients.extend([QQ(0)] * (3 - len(alpha_coefficients)))
            candidates_for_filter.append(
                (
                    {
                        "label": f"norm-kernel-{dependency_index}",
                        "generator_coefficients": [str(value) for value in alpha_coefficients],
                        "source_candidate_indices": source_indices,
                        "source_candidate_labels": [
                            str(source_candidates[index][0]["label"])
                            for index in source_indices
                        ],
                    },
                    alpha,
                    norm,
                )
            )
    else:
        selmer_primes = None
        dependencies = None
        candidates_for_filter = source_candidates

    survivors = []
    norm_rejected = []
    globally_square_rejected = []
    for candidate, alpha, norm in candidates_for_filter:
        if is_rational_square(norm):
            try:
                square_root = alpha.sqrt()
            except ValueError:
                square_root = None
            if square_root is not None and square_root**2 == alpha:
                globally_square_rejected.append(
                    {
                        "label": str(candidate["label"]),
                        "norm": str(norm),
                        "square_root_coefficients": [
                            str(value) for value in square_root.polynomial().list()
                        ],
                    }
                )
            else:
                survivor = dict(candidate)
                survivor["norm"] = str(norm)
                survivors.append(survivor)
        else:
            norm_rejected.append({"label": str(candidate["label"]), "norm": str(norm)})

    output = {
        "schema": "elliptic-curves.bnf-free-norm-filtered-squareclass-candidates.v1",
        "status": "nontrivial_norm_square_candidates_not_local_selmer_certificate",
        "field_polynomial_ascending": record["field_polynomial_ascending"],
        "candidates": survivors,
        "norm_rejected": norm_rejected,
        "globally_square_rejected": globally_square_rejected,
        "source_candidate_count": len(source_candidates),
    }
    if args.generate_norm_kernel:
        output.update(
            {
                "norm_kernel_mode": True,
                "norm_kernel_selmer_rational_primes": selmer_primes,
                "norm_kernel_dimension": len(dependencies),
                "norm_kernel_source_count": len(source_candidates),
            }
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|stage=complete|candidates={len(survivors)}"
        f"|norm_rejected={len(norm_rejected)}"
        f"|globally_square_rejected={len(globally_square_rejected)}"
        f"|norm_kernel_dimension={len(dependencies) if dependencies is not None else 'none'}"
        "|status=GLOBAL_NECESSARY_CONDITIONS_ONLY",
        flush=True,
    )


if __name__ == "__main__":
    main()
