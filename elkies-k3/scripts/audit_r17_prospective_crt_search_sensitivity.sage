#!/usr/bin/env sage-python
"""Post-experiment positive-control sensitivity audit for the frozen CRT search.

This diagnostic applies the already-frozen v2 completed-square point-search
call to historical record fibres 356 and 385.  It does not alter the frozen
candidate lists, protocol, or analysis contrasts.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import runpy
from typing import Any

from sage.all import PolynomialRing, QQ, ZZ, pari


ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-search-protocol-v2.json"
FINGERPRINTS = ROOT / "artifacts/generated-results/elkies-k3-r17-residual-selmer-fingerprints-v1.json"
LOCAL_IMPLEMENTATION = ROOT / "elkies-k3/scripts/audit_r17_prospective_crt_local_stability.sage"
DIRECT_IMPLEMENTATION = ROOT / "elkies-k3/scripts/run_r17_prospective_crt_direct_point_search.sage"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-search-sensitivity-v1.json"

SCHEMA = "elkies-k3.r17-prospective-crt-search-sensitivity.v1"
EXPECTED_PROTOCOL_HASH = "63d6b9e83f52bc7208b9057298e05941dfcedc85d53f5681186c953498947d4b"
X_HEIGHT = 10_000
CERTIFICATE_PRIME_BOUND = 1000


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def build():
    protocol = json.loads(PROTOCOL.read_text())
    if protocol.get("protocol_definition_sha256") != EXPECTED_PROTOCOL_HASH:
        raise ArithmeticError("the reviewed amended search protocol changed")
    amended = protocol["amended_uniform_bounded_search"]
    if (
        amended["engine"] != "PARI hyperellratpoints"
        or amended["x_numerator_denominator_height"] != X_HEIGHT
        or amended["finite_quotient_certificate_prime_bound"] != CERTIFICATE_PRIME_BOUND
    ):
        raise ArithmeticError("the positive-control diagnostic does not match protocol v2")

    fingerprints = json.loads(FINGERPRINTS.read_text())
    fingerprint_by_curve = {
        int(row["curve_id"]): row for row in fingerprints["fingerprints"]
    }
    if any(
        fingerprint_by_curve[anchor]["certified_known_residual_dimension"] != 12
        for anchor in (356, 385)
    ):
        raise ArithmeticError("the historical +12 positive-control labels changed")

    local = runpy.run_path(str(LOCAL_IMPLEMENTATION))
    direct = runpy.run_path(str(DIRECT_IMPLEMENTATION))
    family = local["Family"]()
    point_key = direct["point_key"]
    python_point = direct["python_point"]
    certificate_record = direct["certificate_record"]
    find_mod2_reduction_certificate = direct["find_mod2_reduction_certificate"]
    combined_mod2_rank = direct["combined_mod2_rank"]
    find_two_torsion_certificate_prime = direct["find_two_torsion_certificate_prime"]

    controls = []
    for anchor in (356, 385):
        parameter = family.target_parameters[anchor]
        curve, known = family.specialize(parameter)
        search_curve = curve.local_data(2).minimal_model()
        isomorphisms = curve.isomorphisms(search_curve)
        if not isomorphisms:
            raise ArithmeticError("no exact isomorphism to the fixed p=2-minimal model")
        to_search = isomorphisms[0]
        from_search = ~to_search
        ainvs = [ZZ(value) for value in search_curve.a_invariants()]
        if any(QQ(value).denominator() != 1 for value in search_curve.a_invariants()):
            raise ArithmeticError("the fixed search model is not integral")
        a1, a2, a3, a4, a6 = ainvs
        polynomial_ring = PolynomialRing(QQ, "x")
        x_variable = polynomial_ring.gen()
        completed_square = (
            4 * x_variable**3
            + (a1**2 + 4 * a2) * x_variable**2
            + (2 * a1 * a3 + 4 * a4) * x_variable
            + (a3**2 + 4 * a6)
        )
        raw_points = list(pari(completed_square).hyperellratpoints(X_HEIGHT))
        candidates = {}
        for raw_point in raw_points:
            x_coordinate = QQ(raw_point[0])
            completed_y = QQ(raw_point[1])
            if completed_y**2 != completed_square(x_coordinate):
                raise ArithmeticError("PARI returned a point off the completed-square cubic")
            y_coordinate = (completed_y - a1 * x_coordinate - a3) / 2
            search_point = search_curve(x_coordinate, y_coordinate)
            point = from_search(search_point)
            if point not in curve or to_search(point) != search_point:
                raise ArithmeticError("a returned point failed exact model transport")
            candidates.setdefault(point_key(point), point)

        coefficients = [
            Fraction(0),
            Fraction(0),
            Fraction(0),
            Fraction(str(curve.a4())),
            Fraction(str(curve.a6())),
        ]
        known_python = [python_point(point) for point in known]
        selected = []
        certificates = []
        uncertified = []
        torsion_prime = None
        for candidate in candidates.values():
            trial = known_python + [python_point(point) for point in selected + [candidate]]
            signatures = find_mod2_reduction_certificate(
                coefficients, trial, prime_bound=CERTIFICATE_PRIME_BOUND
            )
            rank = combined_mod2_rank(signatures, len(trial))
            if rank == len(trial):
                if torsion_prime is None:
                    torsion_prime = find_two_torsion_certificate_prime(
                        coefficients, prime_bound=200
                    )
                    if torsion_prime is None:
                        raise ArithmeticError("no rational-2-torsion exclusion certificate")
                selected.append(candidate)
                certificates.append(certificate_record(signatures, len(trial)))
            else:
                uncertified.append(
                    {
                        "point": [str(candidate[0]), str(candidate[1])],
                        "achieved_rank": rank,
                        "column_count": len(trial),
                    }
                )
        controls.append(
            {
                "curve_id": anchor,
                "parameter": str(parameter),
                "historically_certified_residual_MW_dimension": 12,
                "search_model_ainvs": [str(value) for value in ainvs],
                "search_model_sha256": canonical_hash([str(value) for value in ainvs]),
                "completed_square_coefficients_low_to_high": [
                    str(value) for value in completed_square.list()
                ],
                "raw_hyperellratpoints_count": len(raw_points),
                "distinct_exact_candidate_count": len(candidates),
                "certified_escape_count": len(selected),
                "certified_points": [
                    {
                        "point": [str(point[0]), str(point[1])],
                        "finite_quotient_certificate": certificate,
                    }
                    for point, certificate in zip(selected, certificates)
                ],
                "uncertified_candidates": uncertified,
                "rational_two_torsion_exclusion_prime": torsion_prime,
                "sensitivity_status": (
                    "REDETECTED_CERTIFIED_ESCAPE"
                    if selected
                    else "NO_ESCAPE_REDETECTED_AT_FROZEN_BOUND"
                ),
            }
        )

    body = {
        "schema": SCHEMA,
        "status": (
            "PASS_POSITIVE_CONTROL_SENSITIVITY"
            if all(row["certified_escape_count"] for row in controls)
            else "FAILED_TO_REDETECT_BOTH_KNOWN_PLUS12_POSITIVE_CONTROLS"
        ),
        "search_protocol_sha256": EXPECTED_PROTOCOL_HASH,
        "post_experiment_diagnostic_only": True,
        "changes_frozen_candidates_protocol_or_contrasts": False,
        "limits": {
            "engine": "PARI hyperellratpoints",
            "x_numerator_denominator_height": X_HEIGHT,
            "same_completed_square_model_as_frozen_protocol": True,
        },
        "controls": controls,
        "interpretation": (
            "Failure to rediscover a point on either known +12 fibre shows that the all-zero "
            "prospective ledger is detector-limited at this bound. It does not rescue an "
            "enrichment claim and does not alter the frozen experiment."
        ),
        "claim_boundary": [
            "This post-experiment positive-control audit cannot redefine the frozen protocol or cohorts.",
            "A bounded failure to rediscover a known point says nothing about exact rank.",
            "The historical residual dimension 12 is an exact known-MW lower bound, not a full Selmer upper bound.",
        ],
    }
    return {
        **body,
        "sensitivity_definition_sha256": canonical_hash(body),
        "inputs": {
            relative(PROTOCOL): digest(PROTOCOL),
            relative(FINGERPRINTS): digest(FINGERPRINTS),
        },
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "command": "sage -python elkies-k3/scripts/audit_r17_prospective_crt_search_sensitivity.sage",
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    document = build()
    serialized = json.dumps(document, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != serialized:
            raise ArithmeticError("stored positive-control sensitivity audit differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        "R17CRTSEARCHSENSITIVITY"
        f"|356={document['controls'][0]['certified_escape_count']}"
        f"|385={document['controls'][1]['certified_escape_count']}"
        f"|status={document['status']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
