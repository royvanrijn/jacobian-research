#!/usr/bin/env python3
"""Certify the numerical hypotheses of the NS0024 rational-marking no-go.

The geometric descent argument and the theorem on rational points of the
Fricke quotient are written in the canonical proof note.  This checker pins
the exact NS0024 Inose frame and verifies that N=475 satisfies Momose's
published criterion via p=19.  It does not reprove Momose's theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
SOURCE = GENERATED / "elkies-k3-ns0024-new-rootless-source-route-v1.json"
DIRECT = GENERATED / "elkies-k3-ns0024-direct-qq-inose-obstruction-v1.json"
OUTPUT = GENERATED / "elkies-k3-ns0024-qq-marking-obstruction-v1.json"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_payload():
    source = json.loads(SOURCE.read_text())
    direct = json.loads(DIRECT.read_text())
    equation_source = source["equation_source"]

    assert equation_source["source_id"] == "NS0024-S001"
    assert equation_source["frame_type"] == "2E8/MW1"
    assert equation_source["mw_height"] == 950
    assert direct["status"] == "PASS_DIRECT_QQ_INOSE_SOURCE_OBSTRUCTION"
    assert direct["required_cyclic_isogeny_degree"] == 475

    level = equation_source["mw_height"] // 2
    witness_prime = 19
    published_finiteness_exceptions_below_300 = (37, 151, 199, 227, 277)
    assert level == 475 == 5**2 * witness_prime
    assert level > 1 and level % witness_prime == 0
    assert witness_prime >= 17
    assert witness_prime not in published_finiteness_exceptions_below_300

    return {
        "schema": "elkies-k3.ns0024-qq-marking-obstruction.v1",
        "status": "PASS_NS0024_QQ_RATIONAL_MARKING_OBSTRUCTION",
        "ns_id": "NS0024",
        "ns_rank": 19,
        "ns_determinant": 950,
        "inose_frame": {
            "source_id": "NS0024-S001",
            "root_type": "2E8",
            "mw_rank": 1,
            "mw_height": 950,
            "isogeny_degree": level,
            "fricke_modular_curve": "X0+(475)",
        },
        "momose_criterion": {
            "level_is_composite": True,
            "prime_divisor": witness_prime,
            "prime_is_at_least_17": True,
            "j0_minus_rational_points_finite": True,
            "published_exceptional_primes_below_300": list(
                published_finiteness_exceptions_below_300
            ),
            "conclusion": (
                "Every QQ-rational point of X0+(475) is a cusp or a CM point."
            ),
        },
        "rank_gate": {
            "geometric_mw_rank": 1,
            "cm_inose_mw_rank_at_least": 2,
            "required_modular_point_type": "noncuspidal and non-CM",
            "allowed_by_momose": False,
        },
        "arithmetic_conclusion": (
            "No characteristic-zero K3 over QQ with geometric NS=NS0024 can "
            "have all nineteen NS divisor classes rational. In particular a "
            "rootless NS0024 fibration cannot have a saturated rational MW17 "
            "basis over QQ(t)."
        ),
        "external_theorem": {
            "name": "Momose's composite-level theorem for X0+(N)(QQ)",
            "reference": (
                "F. Momose, Rational points on the modular curves X_0^+(N), "
                "J. Math. Soc. Japan 39 (1987), Theorem 0.1"
            ),
            "url": "https://doi.org/10.2969/jmsj/03920269",
        },
        "input_hashes": {
            relative(SOURCE): digest(SOURCE),
            relative(DIRECT): digest(DIRECT),
        },
        "proof_boundary": {
            "proved": (
                "Using the cited Momose theorem and the standard Inose/Fricke "
                "moduli correspondence, a full QQ-rational NS0024 marking is "
                "impossible; the exact N=475,p=19 hypotheses and the rank gate "
                "are pinned here."
            ),
            "not_proved": (
                "Nonexistence of geometric NS0024 surfaces, models over number "
                "fields, or QQ models with a proper Galois-invariant sublattice. "
                "The checker does not independently reprove Momose's theorem."
            ),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    output = args.output if args.output.is_absolute() else ROOT / args.output
    rendered = json.dumps(build_payload(), indent=2, sort_keys=True) + "\n"
    if args.check:
        assert output.read_text() == rendered
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered)
    print("PASS ns0024 QQ rational-marking obstruction")


if __name__ == "__main__":
    main()
