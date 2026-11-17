#!/usr/bin/env python3
"""Central exact audit of the degree-six/nine marked-root rank ladder.

The support formulas are combinatorial identities.  The kernel and rank laws
are asserted only for the four explicitly reduced rows; the output labels
them as a finite verified pattern rather than an all-degree theorem.
"""

from __future__ import annotations

import argparse
from math import ceil
import json
from pathlib import Path

from sympy.polys.domains import GF, QQ

from verify_degree_six_relative_quantization_obstruction import (
    S2_SUPPORT as S2_6,
    S4_SUPPORT as S4_6,
    T2_SUPPORT as T2_6,
    T4_SUPPORT as T4_6,
    rank_record as rank_6,
)
from verify_degree_seven_relative_quantization_obstruction import (
    S2_SUPPORT as S2_7,
    S4_SUPPORT as S4_7,
    T2_SUPPORT as T2_7,
    T4_SUPPORT as T4_7,
    rank_record as rank_7,
)
from verify_degree_eight_relative_quantization_obstruction import (
    S2_SUPPORT as S2_8,
    S4_SUPPORT as S4_8,
    T2_SUPPORT as T2_8,
    T4_SUPPORT as T4_8,
    rank_record as rank_8,
)
from verify_degree_nine_relative_quantization_obstruction import (
    S2_SUPPORT as S2_9,
    S4_SUPPORT as S4_9,
    S6_SUPPORT as S6_9,
    T2_SUPPORT as T2_9,
    T4_SUPPORT as T4_9,
    T6_SUPPORT as T6_9,
    rank_record as rank_9,
)


SUPPORTS = {
    6: ((S2_6, T2_6), (S4_6, T4_6)),
    7: ((S2_7, T2_7), (S4_7, T4_7)),
    8: ((S2_8, T2_8), (S4_8, T4_8)),
    9: ((S2_9, T2_9), (S4_9, T4_9), (S6_9, T6_9)),
}
RANKERS = {6: rank_6, 7: rank_7, 8: rank_8, 9: rank_9}


def support_counts(degree: int, correction_half_order: int):
    r = correction_half_order
    left = degree - 2 * r
    right = degree - 1 - 2 * r
    s_count = ((left + 1) * (9 * left + 8 * r - 4) - 2 * ceil(left / 2)) // 4
    t_count = ((right + 1) * (9 * right + 8 * r) - 2 * ceil(right / 2)) // 4
    return s_count, t_count


def normalized_record(record):
    correction_rank = record["h5_correction_rank"]
    return {
        "h3_rank": record["h3_rank"],
        "h3_kernel_dimension": record["h3_kernel_dimension"],
        "h5_correction_rank": correction_rank,
        "h5_correction_kernel_dimension": record.get(
            "h5_correction_kernel_dimension"
        ),
        "h5_strong_span_rank": record["h5_strong_span_rank"],
        "h5_augmented_rank": record["h5_augmented_rank"],
        "h5_output_dimension": record["h5_output_dimension"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    rows = {}
    for degree in range(6, 10):
        support_record = []
        for r, pair in enumerate(SUPPORTS[degree], start=1):
            actual = tuple(map(len, pair))
            predicted = support_counts(degree, r)
            assert actual == predicted
            support_record.append({
                "correction_order": 2 * r,
                "S": actual[0],
                "T": actual[1],
            })
        record = normalized_record(RANKERS[degree](QQ, QQ.one, QQ.zero))
        correction_dimension = support_record[1]["S"] + support_record[1]["T"]
        correction_kernel = correction_dimension - record["h5_correction_rank"]
        if record["h5_correction_kernel_dimension"] is None:
            record["h5_correction_kernel_dimension"] = correction_kernel
        assert record["h3_kernel_dimension"] == 2 * degree - 6
        assert correction_kernel == 2 * degree - 12
        assert record["h5_strong_span_rank"] - record["h5_correction_rank"] == 7
        assert record["h5_augmented_rank"] == record["h5_strong_span_rank"] + 1
        rows[str(degree)] = {"supports": support_record, "ranks": record}

    f17 = GF(17)
    f19 = GF(19)
    bad_reductions = {
        "degree_6_p17": normalized_record(rank_6(f17, f17.one, f17.zero)),
        "degree_8_p17": normalized_record(rank_8(f17, f17.one, f17.zero)),
        "degree_9_p19": normalized_record(rank_9(f19, f19.one, f19.zero)),
    }
    assert bad_reductions["degree_6_p17"]["h5_strong_span_rank"] == 40
    assert bad_reductions["degree_8_p17"]["h3_rank"] == 167
    assert bad_reductions["degree_9_p19"]["h3_rank"] == 226

    certificate = {
        "scope": (
            "exact four-row rank audit and combinatorial support formulas; "
            "no all-degree rank theorem"
        ),
        "verified_degrees": [6, 7, 8, 9],
        "rows": rows,
        "verified_row_patterns": {
            "dim_kernel_d3": "2*n-6",
            "dim_kernel_D5": "2*n-12",
            "rank_M5_minus_rank_D5": 7,
            "generic_constant_rank_jump": 1,
        },
        "bad_reductions": bad_reductions,
        "bad_prime_conclusion": (
            "the audited exceptions do not follow a uniform 2*n+1 rule"
        ),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print("PASS: exact degree-six-through-nine rank ladder and support formulas")
    print("SCOPE: four verified rows; no all-degree rank theorem")


if __name__ == "__main__":
    main()
