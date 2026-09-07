#!/usr/bin/env python3
"""Reproduce one exact CRT/lattice fixture in the Elkies--Klagsbrun family.

This is a regression fixture, not a target hit.  Python's exact rational
arithmetic verifies the congruences, valuations, curve equations, and lattice
identity.  PARI/GP supplies the minimal model, conductor, local reduction data,
and a 2-descent rank interval.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys
from typing import Any

from crt_lattice import crt_pair, gauss_reduce, short_rational_representatives
from ek_k3 import EKK3Family, rational_to_string, valuation
from pari_bridge import minimal_curve_data, pari_version


Q = Fraction
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/elliptic-curves/elliptic_ek_k3_crt_fixture.json"
)
CONDUCTOR_TARGET_LOG = "182.72"


def point_record(point: tuple[Fraction, Fraction]) -> dict[str, str]:
    return {
        "x": rational_to_string(point[0]),
        "y": rational_to_string(point[1]),
    }


def kodaira_symbol(code: int) -> str:
    """Decode the multiplicative PARI Kodaira codes used by this fixture."""

    if code >= 5:
        return f"I_{code - 4}"
    raise ValueError(f"fixture expected multiplicative reduction, got code {code}")


def build_artifact() -> dict[str, Any]:
    family = EKK3Family(u=Q(2, 5), m=Q(2))
    requested = (
        (11, 2, 5, 5),
        (17, 2, 109, 8),
        (19, 2, 102, 7),
    )

    constraint_records: list[dict[str, Any]] = []
    residue, modulus = 0, 1
    for prime, exponent, expected_root, expected_factor in requested:
        roots = family.power_roots(prime, exponent)
        matches = tuple(
            root
            for root in roots
            if root.residue == expected_root and root.factor_index == expected_factor
        )
        if len(matches) != 1:
            raise AssertionError(
                f"missing expected p={prime} root {expected_root} of B_{expected_factor}"
            )
        root = matches[0]
        if root.split_multiplicative:
            raise AssertionError("this fixture expects nonsplit multiplicative roots")
        residue, modulus = crt_pair(residue, modulus, root.residue, root.modulus)
        constraint_records.append(
            {
                "prime": prime,
                "exponent": exponent,
                "modulus": root.modulus,
                "residue": root.residue,
                "B_factor_index": root.factor_index,
                "predicted_reduction": "nonsplit multiplicative",
            }
        )

    if (residue, modulus) != (7_814_669, 12_623_809):
        raise AssertionError("CRT regression changed")
    reduced_basis = gauss_reduce((modulus, 0), (residue, 1))
    if reduced_basis != ((-1468, 21), (-187, 8602)):
        raise AssertionError("Gauss-reduced lattice basis regression changed")
    representatives = short_rational_representatives(
        residue, modulus, coefficient_radius=8, limit=1
    )
    if len(representatives) != 1:
        raise AssertionError("rational reconstruction returned no representative")
    representative = representatives[0]
    if (representative.numerator, representative.denominator) != (-1468, 21):
        raise AssertionError("short rational representative regression changed")

    parameter = Q(representative.numerator, representative.denominator)
    if (parameter.numerator - residue * parameter.denominator) % modulus != 0:
        raise AssertionError("rational representative does not satisfy the CRT class")
    if not family.is_nonsingular(parameter):
        raise AssertionError("the reconstructed specialization is singular")

    invariants = family.invariants(parameter)
    factors = family.b_factors(parameter)
    exact_local: dict[str, Any] = {}
    for prime, _, _, factor_index in requested:
        factor_valuation = valuation(factors[factor_index - 1], prime)
        discriminant_valuation = valuation(invariants["discriminant"], prime)
        c4_valuation = valuation(invariants["c4"], prime)
        if (factor_valuation, discriminant_valuation, c4_valuation) != (2, 4, 0):
            raise AssertionError(f"unexpected exact valuations at p={prime}")
        # local_data receives t mod p, so divide by the rational denominator.
        residue_mod_prime = (
            parameter.numerator * pow(parameter.denominator, -1, prime) % prime
        )
        local = family.local_data(residue_mod_prime, prime)
        if local.good_reduction or local.split_multiplicative is not False:
            raise AssertionError(f"unexpected local classification at p={prime}")
        exact_local[str(prime)] = {
            "B_factor_index": factor_index,
            "B_factor_valuation": factor_valuation,
            "c4_valuation": c4_valuation,
            "discriminant_valuation": discriminant_valuation,
            "split_multiplicative": local.split_multiplicative,
        }

    points = family.known_points(parameter)
    if len(points) != 9:
        raise AssertionError("expected nine published sections")
    pari = minimal_curve_data(
        family.coefficients(parameter),
        timeout=60.0,
        rank_effort=0,
        known_points=points,
        local_primes=tuple(item[0] for item in requested),
    )
    if pari["conductor"] != int(
        "35290917445780946083484251246019551569185126876962439304281102406"
    ):
        raise AssertionError("conductor regression changed")
    if pari["pari_ellrank"]["lower_bound"] != 10:
        raise AssertionError("PARI failed to recover the recorded rank lower bound")
    if pari["pari_ellrank"]["upper_bound"] != 10:
        raise AssertionError("PARI rank upper bound regression changed")
    if not pari["supplied_points"]["all_on_curve"]:
        raise AssertionError("PARI rejected a published section")
    for prime, _, _, _ in requested:
        reduction = pari["local_reduction"][str(prime)]
        if (
            reduction["conductor_exponent"] != 1
            or reduction["kodaira_code"] != 8
            or reduction["minimal_c4_valuation"] != 0
            or reduction["minimal_discriminant_valuation"] != 4
            or reduction["ellap"] != -1
        ):
            raise AssertionError(f"PARI local reduction regression changed at p={prime}")
        reduction["kodaira_symbol"] = kodaira_symbol(reduction["kodaira_code"])
        reduction["multiplicative_splitness_from_ellap"] = "nonsplit"

    log_conductor = pari["log_conductor"]
    below_conductor_threshold = Q(log_conductor) < Q(CONDUCTOR_TARGET_LOG)
    target_hit = (
        pari["pari_ellrank"]["lower_bound"] >= 21 and below_conductor_threshold
    ) or pari["pari_ellrank"]["lower_bound"] >= 30
    if target_hit:
        raise AssertionError("the fixed regression fixture must not be labeled a target hit")

    return {
        "schema_version": 1,
        "artifact_kind": "elliptic_curve_computational_regression_fixture",
        "status": "verified_regression_fixture_not_target_hit",
        "claim_scope": {
            "exact": (
                "CRT, rational reconstruction, rational curve equations, and "
                "p-adic valuations were checked with exact integer/rational arithmetic"
            ),
            "computer_algebra": (
                "minimal model, conductor, local reduction, and ellrank bounds are "
                "reported by PARI/GP and are not an independent proof certificate"
            ),
        },
        "reproduction": {
            "command": (
                "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
                "elliptic-curves/cas/verify_ek_k3_fixture.py"
            ),
            "python": sys.version.split()[0],
            "pari_gp": pari_version(),
        },
        "family": {
            "name": "Elkies--Klagsbrun rank-nine K3 family",
            "u": "2/5",
            "m": "2",
            "parameter_t": rational_to_string(parameter),
        },
        "local_constraints": constraint_records,
        "crt": {
            "residue": residue,
            "modulus": modulus,
            "congruence": "a - residue*b == 0 (mod modulus)",
            "gauss_reduced_basis": [list(vector) for vector in reduced_basis],
            "selected_vector": [representative.numerator, representative.denominator],
            "projective_height": representative.height,
        },
        "exact_local_checks": exact_local,
        "curve": {
            "input_coefficients_a1_a2_a3_a4_a6": [
                rational_to_string(value) for value in family.coefficients(parameter)
            ],
            "known_sections": [point_record(point) for point in points],
        },
        "pari": pari,
        "target_evaluation": {
            "criterion": "rank >= 21 and log(N) < 182.72, or rank >= 30",
            "log_conductor_below_182_72": below_conductor_threshold,
            "certified_rank_lower_bound_at_least_21": False,
            "certified_rank_lower_bound_at_least_30": False,
            "target_hit": False,
            "note": (
                "The conductor is below the requested threshold, but PARI reports "
                "rank bounds [10,10], far below either rank target."
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    artifact = build_artifact()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    rank = artifact["pari"]["pari_ellrank"]
    print(f"wrote {args.output}")
    print(
        f"t={artifact['family']['parameter_t']} "
        f"log(N)={artifact['pari']['log_conductor']} "
        f"PARI rank bounds=[{rank['lower_bound']},{rank['upper_bound']}]"
    )
    print("target_hit=false")


if __name__ == "__main__":
    main()
