#!/usr/bin/env python3
"""Verify the bounded t0-open directional Fitting reconstructions."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts" / "generated-results"
OUTPUT = ARTIFACTS / "two_pair_sic_bidegree33_t0_fitting_degree_scout.json"
VARIABLES = ("s1", "s2", "t1", "t2", "u")
EXPECTED = {
    "s1": {
        "c0_constant": (145, 59),
        "c0_s3": (145, 59),
        "c1_constant": (155, 66),
        "c1_s3": (155, 66),
    },
    "s2": {
        "c0_constant": (94, 50),
        "c0_s3": (92, 50),
        "c1_constant": (99, 54),
        "c1_s3": (97, 54),
    },
    "t1": {
        "c0_constant": (124, 47),
        "c0_s3": (123, 47),
        "c1_constant": (129, 50),
        "c1_s3": (128, 50),
    },
    "t2": {
        "c0_constant": (83, 41),
        "c0_s3": (81, 41),
        "c1_constant": (85, 42),
        "c1_s3": (84, 42),
    },
    "u": {
        "c0_constant": (177, 91),
        "c0_s3": (176, 91),
        "c1_constant": (185, 96),
        "c1_s3": (184, 96),
    },
}


def artifact_path(variable: str, prime: int, base: int) -> Path:
    suffix = "" if base == 1 else "_base2"
    return ARTIFACTS / (
        f"two_pair_sic_bidegree33_t0_fitting_{variable}"
        f"{suffix}_mod{prime}.json"
    )


def main() -> None:
    records = []
    for base, primes in ((1, (1019, 2039)), (2, (1019,))):
        for prime in primes:
            for variable in VARIABLES:
                path = artifact_path(variable, prime, base)
                payload = json.loads(path.read_text(encoding="utf-8"))
                assert payload["prime"] == prime
                assert payload["varied_base_variable"] == variable
                assert payload["sample_count"] == 450
                reconstruction = payload["rational_reconstruction"]
                observed = {
                    coordinate: (
                        datum["numerator_degree"],
                        datum["denominator_degree"],
                    )
                    for coordinate, datum in reconstruction.items()
                }
                assert observed == EXPECTED[variable]
                assert all(
                    datum["held_out_count"] == 50
                    and datum["held_out_verified"]
                    for datum in reconstruction.values()
                )
                denominator = payload["denominator_model"]
                assert denominator["verified"]
                assert denominator["models"]["c0"]["formula"] == (
                    "a2^41*Q^3*J^3"
                )
                assert denominator["models"]["c1"]["formula"] == (
                    "a2^42*Q^4*J^4"
                )
                assert all(
                    all(model["coordinate_matches"].values())
                    for model in denominator["models"].values()
                )
                records.append(
                    {
                        "base": base,
                        "prime": prime,
                        "variable": variable,
                        "degree_profile": observed,
                        "candidate_interval": payload[
                            "candidate_interval"
                        ],
                        "artifact": str(path.relative_to(ROOT)),
                        "sha256": hashlib.sha256(
                            path.read_bytes()
                        ).hexdigest(),
                    }
                )
    summary = {
        "format": "two-pair-sic-bidegree33-t0-fitting-degree-scout-v1",
        "status": (
            "exact bounded modular directional reconstructions; "
            "not a multivariate Fitting certificate"
        ),
        "coefficient_definitions": {
            "c0": "det(M_mu6)",
            "c1": "coefficient of z in det(M_mu6+z*M_mu7)",
        },
        "residue_basis": ["1", "s3"],
        "directional_degree_profiles": EXPECTED,
        "universal_denominator_model_observed": {
            "c0": "a2^41*Q^3*J^3",
            "c1": "a2^42*Q^4*J^4",
            "a2": "coefficient of s3^2 in mu_3",
        },
        "bases": 2,
        "primes": [1019, 2039],
        "line_reconstructions": len(records),
        "paired_samples_per_line": 450,
        "training_samples_per_line": 400,
        "held_out_samples_per_line": 50,
        "records": records,
        "interpretation_limit": (
            "Matching directional degrees and denominator factors give "
            "stable interpolation bounds. They do not reconstruct the "
            "five-variable numerators or prove a common-root exclusion."
        ),
    }
    OUTPUT.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("PASS 15 bounded directional Fitting reconstructions")
    print("PASS 6750 paired-root samples and 750 held-out checks")
    print("PASS denominator models match at two bases and two primes")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
