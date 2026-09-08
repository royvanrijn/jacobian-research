#!/usr/bin/env sage
"""Compare exact BNF-free S-class factor-base envelopes for the four targets.

This is a planning calculation, not a class-group or Selmer computation.  For
each fixed no-rational-2-torsion curve it constructs the monic integral cubic
for the 2-division field, computes its exact field discriminant and signature,
and materializes every prime ideal required by Bach's ERH generation bound,
together with the elliptic Selmer primes.  The resulting size orders the
targets for the custom relation collector without using a heuristic score.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from sage.all import NumberField, PolynomialRing, QQ, ZZ, factor, prime_range


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "conductor_first_near_miss_descent_targets_v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "conductor_first_s_class_envelopes_v1.json"
)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def generation_bounds(field) -> tuple[int, int]:
    # Import the shared interval-arithmetic implementation so the bound used
    # here is identical to the one enforced by the relation-ledger auditor.
    from audit_bnf_free_s_class_quotient import bounds

    minkowski, bach = bounds(field)
    return int(minkowski), int(bach)


def invariants(model: list[str]) -> tuple[int, int, int, int]:
    a1, a2, a3, a4, a6 = (ZZ(value) for value in model)
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    delta = -b2 * b2 * b8 - 8 * b4**3 - 27 * b6 * b6 + 9 * b2 * b4 * b6
    if delta == 0:
        raise ValueError("singular target model")
    return int(b2), int(b4), int(b6), int(delta)


def prime_signature(prime) -> tuple[str, int, int, int]:
    return (
        str(prime.pari_hnf()),
        int(prime.norm()),
        int(prime.residue_class_degree()),
        int(prime.smallest_integer()),
    )


def target_record(target: dict[str, Any]) -> dict[str, Any]:
    b2, b4, b6, curve_delta = invariants(target["global_minimal_model"])
    # If z=4x, then 4x^3+b2*x^2+2*b4*x+b6=0 becomes the
    # monic integral cubic z^3+b2*z^2+8*b4*z+16*b6=0.
    coefficients = [16 * b6, 8 * b4, b2, 1]
    ring = PolynomialRing(QQ, "z")
    z = ring.gen()
    polynomial = sum(QQ(value) * z**index for index, value in enumerate(coefficients))
    if not polynomial.is_irreducible():
        raise ValueError(f"{target['id']} has reducible 2-division cubic")
    field = NumberField(polynomial, "theta")
    minkowski_bound, bach_bound = generation_bounds(field)

    selmer_rational_primes = sorted(
        {ZZ(2)} | {ZZ(prime) for prime, _ in factor(abs(ZZ(curve_delta)))}
    )
    bach_rational_primes = tuple(prime_range(2, bach_bound + 1))
    factor_base = {
        prime_signature(prime)
        for rational_prime in bach_rational_primes
        for prime in field.primes_above(rational_prime)
    }
    factor_base_without_large_s = len(factor_base)
    s_prime_ideals = {
        prime_signature(prime)
        for rational_prime in selmer_rational_primes
        for prime in field.primes_above(rational_prime)
    }
    factor_base.update(s_prime_ideals)
    real_places, complex_places = field.signature()
    return {
        "id": target["id"],
        "certified_known_rank": int(target["certified_known_rank"]),
        "rank21_directions_needed": int(target["rank21_directions_needed"]),
        "two_division_cubic_coefficients_ascending": [str(value) for value in coefficients],
        "field_discriminant": str(field.discriminant()),
        "field_discriminant_bits": int(abs(ZZ(field.discriminant())).nbits()),
        "field_signature": [int(real_places), int(complex_places)],
        "unconditional_minkowski_generation_bound": str(minkowski_bound),
        "bach_erh_generation_bound": bach_bound,
        "bach_rational_prime_count": len(bach_rational_primes),
        "bach_prime_ideal_count": factor_base_without_large_s,
        "selmer_rational_primes": [str(value) for value in selmer_rational_primes],
        "selmer_rational_prime_count": len(selmer_rational_primes),
        "selmer_prime_ideal_count": len(s_prime_ideals),
        "large_selmer_prime_ideals_added_to_bach_base": (
            len(factor_base) - factor_base_without_large_s
        ),
        "bach_plus_selmer_prime_ideal_count": len(factor_base),
        "scope": (
            "exact factor-base envelope; Bach generation is conditional on "
            "ERH/GRH and no principal-relation completeness is asserted"
        ),
    }


def build(manifest_path: Path) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text())
    records = [target_record(target) for target in manifest["targets"]]
    priority = sorted(
        records,
        key=lambda record: (
            record["bach_plus_selmer_prime_ideal_count"],
            record["id"],
        ),
    )
    payload: dict[str, Any] = {
        "schema": "elliptic-curves.conductor-first-s-class-envelopes.v1",
        "classification": "EXACT_PLANNING_COMPUTATION_WITH_ERH_CONDITIONAL_BOUND",
        "source_manifest": str(manifest_path.relative_to(ROOT)),
        "source_manifest_sha256": file_sha256(manifest_path),
        "generation_method": {
            "unconditional": "Minkowski prime-ideal generation bound",
            "conditional": "Bach 12*(log |Delta_K|)^2 under ERH/GRH",
            "factor_base": (
                "all prime ideals over rational primes through the Bach bound, "
                "union all prime ideals over 2*Delta(E)"
            ),
        },
        "targets": records,
        "collector_priority": [record["id"] for record in priority],
        "interpretation": (
            "This orders the BNF-free relation collectors by an exact materialized "
            "factor-base size. It is not a class-group, Selmer, or rank bound."
        ),
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elliptic-curves/cas/build_conductor_first_s_class_envelopes.py --check"
        ),
    }
    payload["result_sha256"] = stable_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build(args.manifest.resolve())
    if args.check:
        existing = json.loads(args.output.read_text())
        if existing != payload:
            raise SystemExit(f"stale artifact: {args.output}")
        print(
            "CFNMSCLASS|stage=check|status=PASS"
            f"|priority={','.join(payload['collector_priority'])}"
            f"|sha256={payload['result_sha256']}"
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "CFNMSCLASS|stage=build|status=PASS"
        f"|priority={','.join(payload['collector_priority'])}"
        f"|sha256={payload['result_sha256']}"
    )


if __name__ == "__main__":
    main()
