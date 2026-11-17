#!/usr/bin/env python3
"""Verify the bounded t0-open stratum, border, and rank continuation."""

from __future__ import annotations

from collections import defaultdict
from functools import reduce
import hashlib
from itertools import product
import json
from pathlib import Path
import re

from verify_two_pair_sic_bidegree33_t0_pencil_random_scout import (
    evaluate_moment,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts" / "generated-results"
OUTPUT = (
    ARTIFACTS
    / "two_pair_sic_bidegree33_t0_strata_rank_continuation.json"
)
PRIME = 43
BASE_VARIABLES = ("s1", "s2", "t1", "t2", "u")

DIRECT_SEEDS = {
    "Q": (501, 811, 821, 831, 841, 842, 843),
    "J": (502, 812, 822, 832, 851, 852, 853),
    "K": (503,),
    "H": (504,),
    "KH": (601,),
    "QJH": (602, 813, 823, 833),
    "JH": (603, 814, 824, 834),
    "JK": (604, 815, 825, 835),
    "a2": (701,),
    "discriminant": (702,),
}
EXPECTED_DIRECT = {
    "Q": (6300, 132, 4),
    "J": (6300, 136, 6),
    "K": (900, 0, 0),
    "H": (900, 0, 0),
    "KH": (900, 0, 0),
    "QJH": (3600, 171, 5),
    "JH": (3600, 151, 0),
    "JK": (3600, 77, 1),
    "a2": (450, 0, 0),
    "discriminant": (450, 0, 1),
}
EXPECTED_LEADING = {
    "Q": (5, ("s6^2", "s5^3", "s6*s5^2"), 578, 36, 132, 132),
    "J": (5, ("s6^2", "s5^3", "s6*s5^2"), 1245, 36, 136, 136),
    "K": (6, ("s6*s5", "s6^3", "s5^4"), 16, 22, 0, 0),
    "H": (6, ("s6^2", "s5^3"), 15, 36, 0, 0),
    "KH": (6, ("s5^2", "s6^3"), 1, 4, 0, 0),
    "QJH": (
        4,
        ("s6*s5", "s6^2", "s5^3"),
        2997,
        59,
        171,
        263,
    ),
    "JH": (
        4,
        ("s6*s5", "s6^2", "s5^3"),
        1399,
        57,
        151,
        226,
    ),
    "JK": (5, ("s6*s5", "s5^3", "s6^3"), 166, 29, 77, 77),
    "a2": (6, ("s6^2", "s6*s5^2", "s5^4"), 187, 33, 0, 0),
}
RANK_SEED_SAMPLES = {
    901: 225,
    902: 450,
    903: 225,
    904: 450,
    905: 225,
    906: 225,
    907: 225,
    908: 225,
    909: 225,
    910: 225,
    911: 225,
    912: 225,
}
EXPECTED_RESULTANTS = {
    "Q": {
        "resultant": (5429, 76),
        "residual_factors": ((20, 195, 2),),
        "pivot": (256, 33, 525, 38),
    },
    "J": {
        "resultant": (16222, 80),
        "residual_factors": ((24, 612, 2),),
        "pivot": (608, 35, 1127, 40),
    },
    "JK": {
        "resultant": (972, 66),
        "residual_factors": ((24, 161, 2),),
        "pivot": (78, 28, 136, 33),
    },
    "QJH": {
        "resultant": (11558, 136),
        "residual_factors": (
            (20, 101, 1),
            (24, 161, 1),
            (24, 163, 3),
        ),
        "pivot": (862, 58, 1122, 63),
    },
    "JH": {
        "resultant": (4865, 128),
        "residual_factors": (
            (20, 54, 3),
            (24, 143, 1),
            (24, 143, 1),
        ),
        "pivot": (413, 56, 549, 61),
    },
}


def direct_path(stratum: str, seed: int) -> Path:
    return ARTIFACTS / (
        "two_pair_sic_bidegree33_t0_stratum_"
        f"{stratum}_random_p43_seed{seed}.json"
    )


def leading_path(stratum: str) -> Path:
    return ARTIFACTS / (
        f"two_pair_sic_bidegree33_t0_stratum_{stratum}_leading_mod43.json"
    )


def resultant_path(stratum: str) -> Path:
    return ARTIFACTS / (
        "two_pair_sic_bidegree33_t0_stratum_"
        f"{stratum}_border_resultant_mod43.json"
    )


def polynomial_terms(
    serialized: str,
    variables: tuple[str, ...],
) -> dict[tuple[int, ...], int]:
    while serialized.startswith("(") and serialized.endswith(")"):
        serialized = serialized[1:-1]
    result: dict[tuple[int, ...], int] = defaultdict(int)
    variable_index = {
        variable: index for index, variable in enumerate(variables)
    }
    for term in re.findall(r"[+-]?[^+-]+", serialized):
        sign = -1 if term.startswith("-") else 1
        body = term[1:] if term[:1] in "+-" else term
        coefficient = sign
        exponents = [0] * len(variables)
        for factor in body.split("*"):
            if factor.isdigit():
                coefficient = coefficient * int(factor) % PRIME
                continue
            match = re.fullmatch(r"([a-z]\w*)(?:\^(\d+))?", factor)
            assert match is not None, factor
            variable = match.group(1)
            assert variable in variable_index, (variable, variables)
            exponents[variable_index[variable]] += int(match.group(2) or 1)
        monomial = tuple(exponents)
        result[monomial] = (result[monomial] + coefficient) % PRIME
    return {
        monomial: coefficient
        for monomial, coefficient in result.items()
        if coefficient
    }


def evaluate_polynomial(
    polynomial: dict[tuple[int, ...], int],
    variables: tuple[str, ...],
    point: dict[str, int],
) -> int:
    return sum(
        coefficient
        * reduce(
            lambda value, datum: (
                value * pow(point[datum[0]] % PRIME, datum[1], PRIME)
            )
            % PRIME,
            zip(variables, monomial, strict=True),
            1,
        )
        for monomial, coefficient in polynomial.items()
    ) % PRIME


def associates_mod(
    left: dict[tuple[int, ...], int],
    right: dict[tuple[int, ...], int],
) -> bool:
    if set(left) != set(right) or not left:
        return False
    first = next(iter(left))
    scalar = right[first] * pow(left[first], -1, PRIME) % PRIME
    return all(
        right[monomial] == scalar * coefficient % PRIME
        for monomial, coefficient in left.items()
    )


def serialized_profile(
    serialized: str,
    variables: tuple[str, ...],
) -> tuple[int, int]:
    while serialized.startswith("(") and serialized.endswith(")"):
        serialized = serialized[1:-1]
    terms = re.findall(r"[+-]?[^+-]+", serialized)
    degrees = []
    for term in terms:
        degree = 0
        for variable, exponent in re.findall(
            r"([a-z]\w*)(?:\^(\d+))?",
            term,
        ):
            if variable in variables:
                degree += int(exponent or 1)
        degrees.append(degree)
    return len(terms), max(degrees)


def coefficient_point(
    stratum: str,
    record: dict[str, object],
) -> dict[str, int]:
    base = record["base"]
    assert isinstance(base, dict)
    point = {key: int(value) % PRIME for key, value in base.items()}
    point["s3"] = int(record["s3"]) % PRIME
    point["ell"] = (
        point["s1"] * point["u"] - point["t1"]
    ) % PRIME
    if stratum == "H":
        x_value = (
            point["s1"] * point["s1"] * point["u"] - point["s2"]
        ) % PRIME
        q_value = (3 * x_value - 13 * point["u"]) % PRIME
        numerator = (
            q_value + 155 * point["u"] * pow(33, -1, PRIME)
        ) % PRIME
        point["r"] = (
            0
            if point["ell"] == 0
            else numerator * pow(point["ell"], -1, PRIME) % PRIME
        )
    return point


def matrix_rank_mod(matrix: list[list[int]]) -> int:
    rows = [[entry % PRIME for entry in row] for row in matrix]
    rank = 0
    for column in range(len(rows[0]) if rows else 0):
        pivot = next(
            (
                row
                for row in range(rank, len(rows))
                if rows[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, PRIME)
        rows[rank] = [entry * inverse % PRIME for entry in rows[rank]]
        for row in range(len(rows)):
            if row == rank or not rows[row][column]:
                continue
            scalar = rows[row][column]
            rows[row] = [
                (left - scalar * right) % PRIME
                for left, right in zip(rows[row], rows[rank], strict=True)
            ]
        rank += 1
        if rank == len(rows):
            break
    return rank


def monomial_exponents(variable_count: int, degree: int) -> list[tuple[int, ...]]:
    return [
        exponents
        for exponents in product(range(degree + 1), repeat=variable_count)
        if sum(exponents) <= degree
    ]


def verify_candidate(candidate: dict[str, object]) -> None:
    assert candidate["through_mu7_dimension"] == 0
    assert candidate["through_mu7_length"] == 1
    higher = candidate["higher_moment_evaluation"]
    assert isinstance(higher, dict)
    assert higher["first_nonzero_order"] == 8
    assert not higher["all_corrected_moments_zero"]
    point = higher["full_parameter_point"]
    assert isinstance(point, dict)
    independently_evaluated = {
        str(order): evaluate_moment(order, point, PRIME)
        for order in range(2, 9)
    }
    assert independently_evaluated == higher["moment_values"]
    assert all(
        independently_evaluated[str(order)] == 0
        for order in range(2, 8)
    )
    assert independently_evaluated["8"] != 0


def main() -> None:
    direct_records = []
    direct_samples: dict[str, list[dict[str, object]]] = defaultdict(list)
    candidate_count = 0
    for stratum, seeds in DIRECT_SEEDS.items():
        total_roots = 0
        total_drops = 0
        total_candidates = 0
        for seed in seeds:
            path = direct_path(stratum, seed)
            payload = json.loads(path.read_text(encoding="utf-8"))
            assert payload["prime"] == PRIME
            assert payload["stratum"] == stratum
            assert payload["random_seed"] == seed
            total_roots += payload["evaluated_mu3_root_count"]
            total_drops += payload["length_drop_sample_count"]
            total_candidates += payload["common_through_mu7_point_count"]
            if stratum != "discriminant":
                samples = payload["direct_samples"]
                assert len(samples) == payload["evaluated_mu3_root_count"]
                direct_samples[stratum].extend(samples)
            for candidate in payload["common_through_mu7_points"]:
                verify_candidate(candidate)
            direct_records.append(
                {
                    "stratum": stratum,
                    "seed": seed,
                    "artifact": str(path.relative_to(ROOT)),
                    "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                }
            )
        assert (total_roots, total_drops, total_candidates) == (
            EXPECTED_DIRECT[stratum]
        )
        candidate_count += total_candidates
    assert candidate_count == 17

    leading_records = []
    for stratum, expected in EXPECTED_LEADING.items():
        (
            quotient_length,
            leading_monomials,
            expected_terms,
            expected_degree,
            expected_drop_zeros,
            expected_total_zeros,
        ) = expected
        path = leading_path(stratum)
        payload = json.loads(path.read_text(encoding="utf-8"))
        assert payload["prime"] == PRIME
        assert payload["stratum"] == stratum
        assert payload["quotient_length"] == quotient_length
        assert tuple(payload["leading_monomials"]) == leading_monomials
        variables = tuple(payload["coefficient_variables"])
        polynomial = polynomial_terms(
            payload["leading_coefficient_lcm"],
            variables,
        )
        assert len(polynomial) == expected_terms
        assert max(map(sum, polynomial)) == expected_degree
        assert len(payload["leading_coefficient_factors"]) == 1
        samples = direct_samples[stratum]
        generic_length = max(
            int(sample["mu4_mu5_length"]) for sample in samples
        )
        drop_zeros = 0
        total_zeros = 0
        missed_drops = 0
        for sample in samples:
            point = coefficient_point(stratum, sample)
            zero = evaluate_polynomial(polynomial, variables, point) == 0
            drop = int(sample["mu4_mu5_length"]) < generic_length
            total_zeros += int(zero)
            drop_zeros += int(zero and drop)
            missed_drops += int(drop and not zero)
        assert missed_drops == 0
        assert (drop_zeros, total_zeros) == (
            expected_drop_zeros,
            expected_total_zeros,
        )
        leading_records.append(
            {
                "stratum": stratum,
                "quotient_length": quotient_length,
                "border_term_count": expected_terms,
                "border_total_degree": expected_degree,
                "sampled_drop_points_on_border": drop_zeros,
                "sampled_border_points": total_zeros,
                "artifact": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )

    degree_four_monomials = monomial_exponents(len(BASE_VARIABLES), 4)
    assert len(degree_four_monomials) == 126
    quartic_records = []
    for stratum in ("Q", "J"):
        points = [
            coefficient_point(stratum, sample)
            for sample in direct_samples[stratum]
            if int(sample["mu4_mu5_length"]) < 5
        ]
        evaluation_matrix = [
            [
                reduce(
                    lambda value, datum: (
                        value
                        * pow(point[datum[0]] % PRIME, datum[1], PRIME)
                    )
                    % PRIME,
                    zip(BASE_VARIABLES, exponents, strict=True),
                    1,
                )
                for exponents in degree_four_monomials
            ]
            for point in points
        ]
        rank = matrix_rank_mod(evaluation_matrix)
        assert rank == 120
        quartic_records.append(
            {
                "stratum": stratum,
                "point_count": len(points),
                "monomial_count": 126,
                "rank": rank,
                "nullity": 126 - rank,
                "interpretation": (
                    "the six relations are the ambient cubic component "
                    "equation and its five linear multiples"
                ),
            }
        )

    rank_records = []
    selected_complement_count = 0
    rank_at_most_four_count = 0
    rank_candidate_count = 0
    for seed, sample_count in RANK_SEED_SAMPLES.items():
        path = ARTIFACTS / (
            "two_pair_sic_bidegree33_t0_rank_complement_random_"
            f"p43_seed{seed}.json"
        )
        payload = json.loads(path.read_text(encoding="utf-8"))
        assert payload["prime"] == PRIME
        assert payload["random_seed"] == seed
        assert payload["sample_count"] == sample_count
        selected_complement_count += (
            payload["selected_pivot_complement_sample_count"]
        )
        rank_at_most_four_count += payload["rank_at_most_four_sample_count"]
        rank_candidate_count += payload["common_pencil_zero_sample_count"]
        assert not payload["rank_at_most_four_samples"]
        for sample in payload["selected_pivot_complement_samples"]:
            assert sample["ranks"]["M6_M7"] == 6
            assert sample["ranks"]["M6_M7_M8"] == 6
        for replay in payload["candidate_replays"]:
            higher = replay["higher_moment_evaluation"]
            assert higher["first_nonzero_order"] == 8
            assert not higher["all_corrected_moments_zero"]
        rank_records.append(
            {
                "seed": seed,
                "accepted_paired_bases": sample_count,
                "selected_pivot_complement_points": payload[
                    "selected_pivot_complement_sample_count"
                ],
                "rank_at_most_four_points": payload[
                    "rank_at_most_four_sample_count"
                ],
                "common_pencil_zero_points": payload[
                    "common_pencil_zero_sample_count"
                ],
                "artifact": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    assert sum(RANK_SEED_SAMPLES.values()) == 3150
    assert selected_complement_count == 4
    assert rank_at_most_four_count == 0
    assert rank_candidate_count == 4

    resultant_records = []
    for stratum, expected in EXPECTED_RESULTANTS.items():
        path = resultant_path(stratum)
        payload = json.loads(path.read_text(encoding="utf-8"))
        assert payload["prime"] == PRIME
        assert payload["stratum"] == stratum
        assert (
            payload["resultant_term_count"],
            payload["resultant_total_degree"],
        ) == expected["resultant"]
        residual_factors = sorted(
            (
                factor["total_degree"],
                factor["term_count"],
                factor["multiplicity"],
            )
            for factor in payload["factors"]
            if factor["total_degree"] >= 10
        )
        assert residual_factors == sorted(expected["residual_factors"])
        for factor in payload["factors"]:
            if factor["total_degree"] >= 10:
                assert factor["gcd_A"] == factor["gcd_B"] == "1"
        pivot = payload["linear_subresultant"]
        assert (
            pivot["A_term_count"],
            pivot["A_total_degree"],
            pivot["B_term_count"],
            pivot["B_total_degree"],
        ) == expected["pivot"]
        resultant_records.append(
            {
                "stratum": stratum,
                "resultant_term_count": payload["resultant_term_count"],
                "resultant_total_degree": payload["resultant_total_degree"],
                "residual_factors": residual_factors,
                "dense_linear_pivot_on_every_residual_factor": True,
                "artifact": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )

    exact_leading_path = (
        ARTIFACTS
        / "two_pair_sic_bidegree33_t0_stratum_Q_leading_exact.json"
    )
    exact_leading = json.loads(
        exact_leading_path.read_text(encoding="utf-8")
    )
    assert exact_leading["prime"] == 0
    assert exact_leading["quotient_length"] == 5
    exact_leading_variables = tuple(exact_leading["coefficient_variables"])
    assert serialized_profile(
        exact_leading["leading_coefficient_lcm"],
        exact_leading_variables,
    ) == (588, 36)
    assert len(exact_leading["leading_coefficient_factors"]) == 1
    modular_leading = json.loads(
        leading_path("Q").read_text(encoding="utf-8")
    )
    assert associates_mod(
        polynomial_terms(
            exact_leading["leading_coefficient_lcm"],
            exact_leading_variables,
        ),
        polynomial_terms(
            modular_leading["leading_coefficient_lcm"],
            tuple(modular_leading["coefficient_variables"]),
        ),
    )

    exact_resultant_path = (
        ARTIFACTS
        / "two_pair_sic_bidegree33_t0_stratum_Q_border_resultant_exact.json"
    )
    exact_resultant = json.loads(
        exact_resultant_path.read_text(encoding="utf-8")
    )
    assert exact_resultant["prime"] == 0
    assert (
        exact_resultant["resultant_term_count"],
        exact_resultant["resultant_total_degree"],
    ) == (5563, 76)
    exact_residual_factors = [
        factor
        for factor in exact_resultant["factors"]
        if factor["total_degree"] >= 10
    ]
    assert [
        (
            factor["total_degree"],
            factor["term_count"],
            factor["multiplicity"],
        )
        for factor in exact_residual_factors
    ] == [(20, 200, 2)]
    assert all(
        factor["gcd_A"] == factor["gcd_B"] == "1"
        for factor in exact_residual_factors
    )
    exact_pivot = exact_resultant["linear_subresultant"]
    assert (
        exact_pivot["A_term_count"],
        exact_pivot["A_total_degree"],
        exact_pivot["B_term_count"],
        exact_pivot["B_total_degree"],
    ) == (262, 33, 535, 38)
    modular_resultant = json.loads(
        resultant_path("Q").read_text(encoding="utf-8")
    )
    assert associates_mod(
        polynomial_terms(
            exact_resultant["resultant"],
            tuple(exact_resultant["base_variables"]),
        ),
        polynomial_terms(
            modular_resultant["resultant"],
            tuple(modular_resultant["base_variables"]),
        ),
    )
    modular_residual = next(
        factor
        for factor in modular_resultant["factors"]
        if factor["total_degree"] == 20
    )
    assert associates_mod(
        polynomial_terms(
            exact_residual_factors[0]["factor"],
            tuple(exact_resultant["base_variables"]),
        ),
        polynomial_terms(
            modular_residual["factor"],
            tuple(modular_resultant["base_variables"]),
        ),
    )
    exact_q_record = {
        "leading_border_term_count": 588,
        "leading_border_total_degree": 36,
        "projected_resultant_term_count": 5563,
        "projected_resultant_total_degree": 76,
        "residual_factor": {
            "total_degree": 20,
            "term_count": 200,
            "multiplicity": 2,
        },
        "dense_linear_pivot": True,
        "good_reduction_matches_mod43": True,
        "leading_artifact": str(exact_leading_path.relative_to(ROOT)),
        "leading_sha256": hashlib.sha256(
            exact_leading_path.read_bytes()
        ).hexdigest(),
        "resultant_artifact": str(exact_resultant_path.relative_to(ROOT)),
        "resultant_sha256": hashlib.sha256(
            exact_resultant_path.read_bytes()
        ).hexdigest(),
    }

    summary = {
        "format": (
            "two-pair-sic-bidegree33-t0-strata-rank-continuation-v1"
        ),
        "status": (
            "exact characteristic-zero Q-border calculation plus bounded "
            "finite-field strata and sampled incidence checks; not a "
            "global common-root exclusion"
        ),
        "prime": PRIME,
        "direct_stratum_artifact_count": len(direct_records),
        "direct_candidate_count": candidate_count,
        "direct_candidates_excluded_by_mu8": candidate_count,
        "direct_records": direct_records,
        "leading_border_records": leading_records,
        "quartic_point_cloud_records": quartic_records,
        "rank_complement": {
            "shard_count": len(rank_records),
            "accepted_paired_base_count": 3150,
            "evaluated_mu3_root_count": 6300,
            "selected_two_pivot_complement_point_count": (
                selected_complement_count
            ),
            "rank_at_most_four_point_count": rank_at_most_four_count,
            "common_pencil_zero_point_count": rank_candidate_count,
            "common_pencil_points_excluded_by_mu8": rank_candidate_count,
            "records": rank_records,
        },
        "projected_border_resultants": resultant_records,
        "characteristic_zero_Q_border": exact_q_record,
        "interpretation_limit": (
            "Only the Q border and its projected residual factor have "
            "been promoted to characteristic zero. No residual component "
            "is excluded through later moments. The next exact step is "
            "arithmetic in the degree-five and degree-six residual "
            "extensions, followed by mu6--mu8 Fitting tests."
        ),
    }
    OUTPUT.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("PASS 31 specialized direct-stratum artifacts")
    print("PASS 17 direct common roots through mu7, all excluded by mu8")
    print("PASS nine exact modular leading-border calculations")
    print("PASS Q and J have no extra degree-at-most-four drop relation")
    print("PASS 6300 rank-complement roots and no rank-at-most-four point")
    print("PASS five factored border resultants with dense linear pivots")
    print("PASS characteristic-zero Q border and degree-20 residual factor")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
