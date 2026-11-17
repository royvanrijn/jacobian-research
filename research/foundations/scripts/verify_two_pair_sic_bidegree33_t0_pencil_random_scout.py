#!/usr/bin/env python3
"""Verify the bounded random determinant-pencil scout on the t0-open."""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from explore_two_pair_sic_bidegree33_full_anchor import (
    PARAMETERS,
    moment_terms,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts" / "generated-results"
OUTPUT = ARTIFACTS / "two_pair_sic_bidegree33_t0_pencil_random_scout.json"
EXPECTED = {
    (43, 401): (1070, 0),
    (43, 402): (1102, 0),
    (43, 403): (1051, 0),
    (43, 404): (1004, 1),
    (43, 405): (1127, 2),
    (43, 406): (1110, 0),
    (43, 407): (1126, 0),
    (43, 408): (1027, 0),
    (43, 409): (1109, 0),
    (43, 410): (1073, 0),
    (43, 411): (1042, 1),
    (43, 412): (1073, 1),
    (43, 413): (1101, 0),
    (43, 414): (1057, 1),
    (43, 415): (1070, 0),
    (43, 416): (1079, 0),
    (43, 417): (1052, 0),
    (43, 418): (1069, 1),
    (43, 419): (1058, 0),
    (43, 420): (1033, 1),
    (43, 421): (1079, 2),
    (43, 422): (1025, 0),
    (43, 423): (1084, 0),
    (43, 424): (1113, 1),
    (47, 101): (1074, 0),
    (47, 102): (1026, 1),
    (47, 103): (1098, 0),
    (47, 104): (1026, 0),
    (47, 105): (1024, 1),
    (47, 106): (1128, 0),
    (47, 107): (1112, 0),
    (47, 108): (1118, 1),
    (47, 109): (1156, 0),
    (47, 110): (1064, 0),
    (47, 111): (1092, 2),
    (47, 112): (1024, 1),
    (59, 201): (1007, 0),
    (59, 202): (989, 1),
    (59, 203): (989, 1),
    (59, 204): (1031, 1),
    (71, 301): (941, 0),
    (71, 302): (932, 0),
    (71, 303): (951, 0),
    (71, 304): (959, 0),
}
EXPECTED_STRATA = {
    ("Q", 501): (1044, {"4": 15, "5": 885}, 1),
    ("J", 502): (954, {"3": 1, "4": 22, "5": 877}, 1),
    ("K", 503): (997, {"6": 900}, 0),
    ("H", 504): (1055, {"6": 900}, 0),
}


def artifact_path(prime: int, seed: int) -> Path:
    return ARTIFACTS / (
        "two_pair_sic_bidegree33_t0_pencil_random_"
        f"p{prime}_seed{seed}.json"
    )


def stratum_artifact_path(stratum: str, seed: int) -> Path:
    return ARTIFACTS / (
        "two_pair_sic_bidegree33_t0_stratum_"
        f"{stratum}_random_p43_seed{seed}.json"
    )


@lru_cache(maxsize=None)
def cached_moment_terms(
    order: int,
    prime: int,
) -> dict[tuple[int, ...], int]:
    return moment_terms(order, prime)


def evaluate_moment(
    order: int,
    point: dict[str, int],
    prime: int,
) -> int:
    result = 0
    for exponents, coefficient in cached_moment_terms(order, prime).items():
        term = coefficient
        for variable, exponent in zip(PARAMETERS, exponents, strict=True):
            term = term * pow(point[variable], exponent, prime) % prime
        result = (result + term) % prime
    return result


def main() -> None:
    records = []
    candidate_records = []
    pivot_signatures = Counter()
    for (prime, seed), (attempted_count, candidate_count) in EXPECTED.items():
        path = artifact_path(prime, seed)
        payload = json.loads(path.read_text(encoding="utf-8"))
        assert payload["format"] == (
            "two-pair-sic-bidegree33-t0-fitting-samples-v2"
        )
        assert payload["sampling_mode"] == "random"
        assert payload["prime"] == prime
        assert payload["random_seed"] == seed
        assert payload["sample_count"] == 450
        assert payload["attempted_base_count"] == attempted_count
        assert payload["common_pencil_zero_sample_count"] == candidate_count
        assert "pairs" not in payload
        candidates = payload["common_pencil_zero_samples"]
        replays = payload["candidate_replays"]
        assert len(candidates) == len(replays) == candidate_count
        for candidate, replay in zip(candidates, replays, strict=True):
            assert candidate["pencil_coefficients"] == [0] * 7
            pivot = candidate["pivot_signature"]
            assert pivot["rank_B7"] == 5
            assert pivot["rank_B8"] == 6
            assert pivot["pivot_rows_zero_based"] == [0, 1, 2, 3, 4]
            assert pivot["mu8_completion_column_zero_based"] == 0
            pivot_key = tuple(pivot["pivot_columns_B7_zero_based"])
            assert pivot_key in {
                (0, 1, 2, 3, 4),
                (0, 1, 2, 3, 6),
            }
            pivot_signatures[pivot_key] += 1
            assert replay["sample_index"] == candidate["sample_index"]
            assert replay["dimension"] == 0
            assert replay["vector_space_dimension"] == 1
            assert replay["standard_basis_size"] == 2
            assert len(replay["standard_basis"]) == 2
            higher = replay["higher_moment_evaluation"]
            assert higher["evaluated_orders"] == list(range(2, 9))
            assert higher["first_nonzero_order"] == 8
            assert not higher["all_corrected_moments_zero"]
            point = higher["full_parameter_point"]
            assert point["t0"] == 1
            assert point["s3"] == candidate["s3"]
            assert point["s0"] * candidate["base"]["u"] % prime == 1
            for variable in ("s1", "s2", "t1", "t2"):
                assert point[variable] == candidate["base"][variable]
            independently_evaluated = {
                str(order): evaluate_moment(order, point, prime)
                for order in range(2, 9)
            }
            assert independently_evaluated == higher["moment_values"]
            assert all(
                independently_evaluated[str(order)] == 0
                for order in range(2, 8)
            )
            assert independently_evaluated["8"] != 0
            candidate_records.append(
                {
                    "prime": prime,
                    "seed": seed,
                    "base": candidate["base"],
                    "s3": candidate["s3"],
                    "standard_basis": replay["standard_basis"],
                    "full_parameter_point": point,
                    "mu8": independently_evaluated["8"],
                    "pivot_signature": pivot,
                }
            )
        records.append(
            {
                "prime": prime,
                "seed": seed,
                "accepted_paired_bases": payload["sample_count"],
                "attempted_bases": attempted_count,
                "common_pencil_zero_points": candidate_count,
                "artifact": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    assert pivot_signatures == {
        (0, 1, 2, 3, 4): 19,
        (0, 1, 2, 3, 6): 1,
    }
    stratum_records = []
    stratum_candidate_records = []
    for (
        stratum,
        seed,
    ), (
        attempted_count,
        expected_lengths,
        expected_candidates,
    ) in EXPECTED_STRATA.items():
        path = stratum_artifact_path(stratum, seed)
        payload = json.loads(path.read_text(encoding="utf-8"))
        assert payload["format"] == (
            "two-pair-sic-bidegree33-t0-specialized-direct-scout-v1"
        )
        assert payload["prime"] == 43
        assert payload["stratum"] == stratum
        assert payload["random_seed"] == seed
        assert payload["attempted_base_count"] == attempted_count
        assert payload["accepted_paired_base_count"] == 450
        assert payload["evaluated_mu3_root_count"] == 900
        assert payload["mu4_mu5_length_distribution"] == expected_lengths
        assert payload["common_through_mu7_point_count"] == (
            expected_candidates
        )
        candidates = payload["common_through_mu7_points"]
        assert len(candidates) == expected_candidates
        for candidate in candidates:
            assert candidate["through_mu7_dimension"] == 0
            assert candidate["through_mu7_length"] == 1
            assert candidate["through_mu7_standard_basis_size"] == 2
            higher = candidate["higher_moment_evaluation"]
            assert higher["first_nonzero_order"] == 8
            assert not higher["all_corrected_moments_zero"]
            point = higher["full_parameter_point"]
            independently_evaluated = {
                str(order): evaluate_moment(order, point, 43)
                for order in range(2, 9)
            }
            assert independently_evaluated == higher["moment_values"]
            assert all(
                independently_evaluated[str(order)] == 0
                for order in range(2, 8)
            )
            assert independently_evaluated["8"] != 0
            stratum_candidate_records.append(
                {
                    "stratum": stratum,
                    "seed": seed,
                    "base": candidate["base"],
                    "s3": candidate["s3"],
                    "standard_basis": candidate[
                        "through_mu7_standard_basis"
                    ],
                    "full_parameter_point": point,
                    "mu8": independently_evaluated["8"],
                }
            )
        stratum_records.append(
            {
                "stratum": stratum,
                "prime": 43,
                "seed": seed,
                "attempted_bases": attempted_count,
                "accepted_paired_bases": 450,
                "mu4_mu5_length_distribution": expected_lengths,
                "common_through_mu7_points": expected_candidates,
                "artifact": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    accepted_pairs = sum(record["accepted_paired_bases"] for record in records)
    attempted_bases = sum(record["attempted_bases"] for record in records)
    summary = {
        "format": "two-pair-sic-bidegree33-t0-pencil-random-scout-v2",
        "status": (
            "exact bounded modular random search; not a characteristic-zero "
            "common-root exclusion"
        ),
        "primes": sorted({prime for prime, _ in EXPECTED}),
        "shard_count": len(records),
        "attempted_base_count": attempted_bases,
        "accepted_paired_base_count": accepted_pairs,
        "evaluated_mu3_root_count": 2 * accepted_pairs,
        "common_mu4_through_mu7_point_count": len(candidate_records),
        "points_excluded_by_mu8": len(candidate_records),
        "survivors_through_mu8": 0,
        "pivot_signature_counts": {
            "M6_columns_1_through_5": 19,
            "M6_columns_1_through_4_then_M7_column_1": 1,
        },
        "determinant_pencil": (
            "det(M_mu6+z*M_mu7), with all seven coefficients tested"
        ),
        "records": records,
        "candidate_records": candidate_records,
        "specialized_strata": {
            "prime": 43,
            "attempted_base_count": sum(
                record["attempted_bases"] for record in stratum_records
            ),
            "accepted_paired_base_count": 450 * len(stratum_records),
            "evaluated_mu3_root_count": 900 * len(stratum_records),
            "common_through_mu7_point_count": len(
                stratum_candidate_records
            ),
            "points_excluded_by_mu8": len(stratum_candidate_records),
            "records": stratum_records,
            "candidate_records": stratum_candidate_records,
        },
        "interpretation_limit": (
            "The twenty generic and two specialized direct finite-field "
            "common roots through mu7 are all excluded by mu8. Random "
            "sampling neither reconstructs the global Fitting ideal nor "
            "proves that no characteristic-zero component survives."
        ),
    }
    OUTPUT.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"PASS {len(records)} deterministic bounded random shards")
    print(
        f"PASS {attempted_bases} attempted bases, "
        f"{accepted_pairs} paired bases, {2 * accepted_pairs} mu3 roots"
    )
    print(
        f"PASS {len(candidate_records)} direct common roots through mu7; "
        "all excluded by mu8"
    )
    print(
        "PASS two generic pivot signatures, with the leading M6 chart "
        "covering 19 of 20 points"
    )
    print(
        f"PASS {len(stratum_records)} specialized divisor scouts and "
        f"{len(stratum_candidate_records)} direct candidates, both "
        "excluded by mu8"
    )
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
