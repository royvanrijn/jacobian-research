#!/usr/bin/env sage-python
"""Audit the first two Frobenius moments for one product twist at one prime.

status: EXACT_AUXILIARY_AUDIT
claim: independent fibrewise n=1,2 moments for a declared extra good prime
inputs: direct alternate-Q80 model, V4 shortlist, and exact rank-one screen
outputs: one compact audit accepted by certify_r17_product_toric_frobenius.sage

The main seventeen-target audit is deliberately frozen at the first two common
good primes 131 and 137.  This helper reuses its exact point-counting routines
without changing that pinned artifact, so a theorem-directed third-prime
calculation can retain the verifier's independent-moment gate.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import importlib.util
import json
from pathlib import Path
import shlex
import sys

from sage.all import ZZ


ROOT = Path(__file__).resolve().parents[2]
BASE_AUDITOR = Path(__file__).with_name(
    "audit_r17_norm12_11952_product_twist_finite_field_bounds.sage"
)
DEFAULT_OUTPUT_DIR = (
    ROOT / "artifacts/generated-results/elkies-k3-r17-product-extra-prime-audits"
)
SCHEMA = "elkies-k3.r17-product-twist-extra-prime-audit.v1"


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def display_path(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def load_auditor():
    specification = importlib.util.spec_from_loader(
        "r17_product_base_auditor",
        SourceFileLoader("r17_product_base_auditor", str(BASE_AUDITOR)),
    )
    if specification is None or specification.loader is None:
        raise ImportError(BASE_AUDITOR)
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


def build_payload(pair_key: str, prime: int):
    auditor = load_auditor()
    model, targets = auditor.load_inputs()
    target = next((row for row in targets if row["pair_key"] == pair_key), None)
    if target is None:
        raise ValueError("--pair-key is not one of the seventeen exact-rank-one targets")

    a_values = model["A_coefficients_low_to_high"]
    b_values = model["B_coefficients_low_to_high"]
    ring, coefficient_a, coefficient_b, discriminant_core = auditor.reduce_model(
        a_values, b_values, prime
    )
    twist = auditor.reduce_twist(
        ring,
        target["product_quartic_coefficients_low_to_high"],
        discriminant_core,
        prime,
    )
    a_coefficients = [int(coefficient_a[index]) for index in range(9)]
    b_coefficients = [int(coefficient_b[index]) for index in range(13)]
    twist_coefficients = [[int(twist[index]) for index in range(5)]]

    first_sums, fp_local_sums = auditor.fp_power_sums(
        a_coefficients, b_coefficients, twist_coefficients, prime
    )
    second_sums, nonsquare, fp2_local_sha = auditor.fp2_power_sums(
        a_coefficients, b_coefficients, twist_coefficients, prime
    )
    first_sum = ZZ(first_sums[0])
    second_sum = ZZ(second_sums[0])
    coefficient_two_numerator = first_sum**2 - second_sum
    if coefficient_two_numerator % 2:
        raise ArithmeticError("Newton coefficient is not integral")

    reduction = {
        "prime": prime,
        "good_reduction": True,
        "twist_factor_degrees": auditor.factor_degrees(twist),
        "twist_discriminant_gcd_degree": int(twist.gcd(discriminant_core).degree()),
        "fibre_configuration_geometric": "4I0*+24I1",
        "elliptic_L_frobenius_power_sums_n1_n2": [
            int(first_sum),
            int(second_sum),
        ],
        "elliptic_L_coefficients_through_T2": [
            1,
            int(-first_sum),
            int(coefficient_two_numerator // 2),
        ],
        "moments_computed": 2,
        "moments_needed_with_functional_equation": 14,
        "rank_zero_decided": False,
        "status": "UNKNOWN_PARTIAL_FROBENIUS_DATA",
    }
    return {
        "schema": SCHEMA,
        "status": "PASS_EXACT_EXTRA_PRIME_N1_N2_AUDIT",
        "pair_key": pair_key,
        "prime": prime,
        "targets": [
            {
                "shortlist_rank": int(target["shortlist_rank"]),
                "pair_key": pair_key,
                "labels": target["labels"],
                "product_quartic_coefficients_low_to_high": target[
                    "product_quartic_coefficients_low_to_high"
                ],
                "reductions": [reduction],
            }
        ],
        "reduction_summary": {
            "fp2_nonsquare": int(nonsquare),
            "original_discriminant_factor_degrees": auditor.factor_degrees(
                discriminant_core
            ),
            "original_nodal_fibre_count_geometric": 24,
            "fp_local_character_sums_sha256": sha256(
                ",".join(map(str, fp_local_sums)).encode("ascii")
            ).hexdigest(),
            "fp2_local_character_sums_sha256": fp2_local_sha,
        },
        "method": {
            "local_trace_identity": (
                "a_(p^n)(E^d_u)=chi_(p^n)(d(u))*a_(p^n)(E_u)"
            ),
            "independence": (
                "direct fibrewise character sums, independent of toric controlled reduction"
            ),
        },
        "inputs": {
            display_path(path): digest(path)
            for path in (
                Path(__file__).resolve(),
                BASE_AUDITOR,
                auditor.DIRECT,
                auditor.SHORTLIST,
                auditor.RANK_SCREEN,
            )
        },
        "proof_boundary": (
            "The good-reduction checks and first two Frobenius power sums are exact. "
            "This audit alone supplies neither the complete degree-28 polynomial nor "
            "a Mordell--Weil rank bound."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair-key", required=True)
    parser.add_argument("--prime", required=True, type=int)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.prime < 5:
        parser.error("--prime must be an odd prime at least five")
    if args.output is None:
        tag = args.pair_key.replace(":", "--")
        args.output = DEFAULT_OUTPUT_DIR / f"{tag}-p{args.prime}-v1.json"

    payload = build_payload(args.pair_key, int(args.prime))
    payload["reproducing_command"] = shlex.join(sys.argv)
    if args.check:
        stored = json.loads(args.output.read_text())
        stored.pop("reproducing_command", None)
        payload.pop("reproducing_command", None)
        if stored != payload:
            raise ArithmeticError("stored extra-prime audit does not replay")
        print(
            f"R17PRODUCTEXTRACHECK|pair={args.pair_key}|p={args.prime}"
            f"|moments=2/14|status=PASS",
            flush=True,
        )
        return

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"R17PRODUCTEXTRA|pair={args.pair_key}|p={args.prime}"
        f"|moments=2/14|output={display_path(args.output)}|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
