#!/usr/bin/env python3
"""Audit pinned ICARM rank-30/31 artifacts without recomputing arithmetic."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
EXPECTED_HASHES = {
    "elliptic-curves/cas/verify_icarm_curve273_rank30.py":
        "8076eb41d755472f623df920c80c3ab94e4599305cfeab38e2337050ba9a0c70",
    "elliptic-curves/cas/icarm_curve273.py":
        "1752a091340a561ac73ca1f07d81fbbe207d45aaf7620e07d0332060c27b01ba",
    "elliptic-curves/cas/mod2_reduction_independence.py":
        "085f9f0893a0c48244fd4b32203b27ea57d8f8733a6da3835ee0808748e017c8",
    "elliptic-curves/cas/pari_bridge.py":
        "6feaf108994f05709bb3737f14fffc744ad039f8975e15bd1f59df92fa93e219",
    "elliptic-curves/cas/search_extra_points.py":
        "8719a8868f60a055ecfc2c68e95b253386a0f5e856c2686749bd8d8e1f83775e",
    "elliptic-curves/scripts/verify_icarm_curve273_rank30_sage.py":
        "51c5724ccea8868ec713c800485d84c85bcdb1fcef85c7a9db45103931879e37",
    "elliptic-curves/cas/check_icarm_curve302_rank31_pinned.py":
        "056dbc90c1490d6457c6573f0367e0fd573da6a2fbafc95e9017f85a7f587e5a",
    "elliptic-curves/cas/verify_icarm_curve302_rank31.py":
        "f06304f1991992323f19c2c695873afc99d8fb697a5085d2cfa9aa3c523bc0cd",
    "elliptic-curves/cas/icarm_curve302.py":
        "7b3c5f92f92278b7f0823114ef8e9967a9e8c867f1435238a2521938380ec73f",
    "elliptic-curves/cas/verify_icarm_curve398_rank30.py":
        "91004711de569b0ab7f57c7850bfaa7e29c21e3942e02a9ef85dd7967e260a31",
    "elliptic-curves/cas/icarm_curve398.py":
        "6a2769de8822bf652668f9749e68229256596a97625b509388394f2f6439156f",
    "elliptic-curves/data/icarm_curve398_known_a1_mod179.json":
        "663071aadfdb0325f623dd04e5f05a9515f61f765c0e0495d8888826e66acee9",
    "artifacts/generated-results/elliptic-curves/icarm_curve273_rank30_v1.json":
        "e2a7a322fbd4703af4239f497749a69a68f9d5149aa8a1f696b39ab3941a3284",
    "artifacts/generated-results/elliptic-curves/icarm_curve302_rank31_v1.json.gz":
        "fc50b4b9ec5fe1dd1fe31aa299f13d8bc3476d43f3ed98e2ade5a4fc8972aa04",
    "artifacts/generated-results/elliptic-curves/icarm_curve398_rank30_and_construction_v1.json":
        "1fd4f23ff2167321be0e3a7bf12b693f0a9ebe26d1e2125ce131da30ad05bf60",
}
EXPECTED_302_JSON_HASH = (
    "3be0d6fe82c58e0f9284df5d9340332944a1d906508ea986d4abe00357036991"
)


def sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def flattened_rows(signatures: list[dict[str, object]]) -> list[list[int]]:
    return [
        row
        for signature in signatures
        for row in signature["matrix_rows"]
    ]


def assert_matrix_shape(rows: list[list[int]], height: int, width: int) -> None:
    assert len(rows) == height
    assert all(len(row) == width for row in rows)
    assert all(value in (0, 1) for row in rows for value in row)


def main() -> None:
    for relative_path, expected_hash in EXPECTED_HASHES.items():
        actual_hash = sha256((ROOT / relative_path).read_bytes())
        assert actual_hash == expected_hash, (
            relative_path,
            actual_hash,
            expected_hash,
        )

    curve273 = json.loads(
        (
            ROOT
            / "artifacts/generated-results/elliptic-curves/"
            "icarm_curve273_rank30_v1.json"
        ).read_text(encoding="utf-8")
    )
    assert curve273["artifact_kind"] == "exact_elliptic_curve_rank_lower_bound"
    assert curve273["claim"] == "rank E(Q) >= 30"
    assert curve273["claim_status"] == (
        "exact unconditional lower bound; no exact-rank claim"
    )
    assert len(curve273["points"]) == 30
    assert curve273["point_membership_checks"] == 30
    certificate273 = curve273["independence_certificate"]
    assert certificate273["method"] == (
        "finite good-reduction quotients E(F_p)/2E(F_p)"
    )
    assert certificate273["no_rational_2_torsion_witness_prime"] == 23
    assert certificate273["combined_binary_rank"] == 30
    rows273 = flattened_rows(certificate273["rows"])
    assert_matrix_shape(rows273, 31, 30)
    assert curve273["curve"]["torsion_order"] == 1
    assert curve273["curve"]["minimal_model_same"] is True
    assert curve273["height_diagnostic"]["used_for_rank_claim"] is False
    assert curve273["generation"]["checker_sha256"] == EXPECTED_HASHES[
        "elliptic-curves/cas/verify_icarm_curve273_rank30.py"
    ]
    assert curve273["generation"]["model_data_sha256"] == EXPECTED_HASHES[
        "elliptic-curves/cas/icarm_curve273.py"
    ]

    compressed302 = (
        ROOT
        / "artifacts/generated-results/elliptic-curves/"
        "icarm_curve302_rank31_v1.json.gz"
    ).read_bytes()
    rendered302 = gzip.decompress(compressed302)
    assert sha256(rendered302) == EXPECTED_302_JSON_HASH
    curve302 = json.loads(rendered302)
    assert curve302["artifact_kind"] == "exact_elliptic_curve_rank_lower_bound"
    assert curve302["claim"] == "rank E(Q) >= 31"
    assert curve302["claim_status"] == (
        "exact unconditional lower bound; no unconditional exact-rank claim"
    )
    assert curve302["point_count"] == 31
    assert curve302["exact_membership_checks_passed"] == 31
    certificate302 = curve302["independence_certificate"]
    assert certificate302["combined_binary_rank"] == 31
    assert certificate302["matrix_row_count"] == 32
    rows302 = flattened_rows(certificate302["signatures"])
    assert_matrix_shape(rows302, 32, 31)
    cross_check = curve302["quadratic_character_cross_check"]
    assert cross_check["combined_binary_rank"] == 31
    assert_matrix_shape(cross_check["matrix_rows"], 32, 31)
    assert curve302["torsion_certificate"]["gcd"] == 1
    assert curve302["curve"]["torsion_subgroup"] == "trivial"
    assert curve302["curve"]["global_minimal_model"] is True
    assert curve302["public_conditional_statement"]["used_for_rank_lower_bound"] is False
    assert curve302["height_diagnostic"]["used_for_rank_claim"] is False
    assert curve302["construction_status"][
        "identified_as_H3_or_R17_specialization"
    ] is False
    assert curve302["generation"]["checker_sha256"] == EXPECTED_HASHES[
        "elliptic-curves/cas/verify_icarm_curve302_rank31.py"
    ]
    assert curve302["generation"]["model_data_sha256"] == EXPECTED_HASHES[
        "elliptic-curves/cas/icarm_curve302.py"
    ]

    curve398 = json.loads(
        (
            ROOT
            / "artifacts/generated-results/elliptic-curves/"
            "icarm_curve398_rank30_and_construction_v1.json"
        ).read_text(encoding="utf-8")
    )
    assert curve398["artifact_kind"] == (
        "exact_elliptic_curve_rank_lower_bound_and_construction_dissection"
    )
    assert curve398["claim"] == "rank E(Q) >= 30"
    assert curve398["claim_status"] == (
        "exact unconditional lower bound; no unconditional exact-rank claim"
    )
    assert curve398["point_membership_checks"] == 30
    assert len(curve398["points"]) == 30
    certificate398 = curve398["independence_certificate"]
    assert certificate398["combined_binary_rank"] == 30
    rows398 = flattened_rows(certificate398["rows"])
    assert_matrix_shape(rows398, 31, 30)
    assert certificate398["full_torsion_witness_group_orders"] == {
        "11": 18,
        "23": 31,
    }
    assert curve398["curve"]["torsion_order"] == 1
    assert curve398["curve"]["global_minimal_model"] is True
    assert curve398["curve"]["root_number"] == 1
    assert curve398["rational_isogeny_certificate"] == {
        "conclusion": (
            "The Q-isogeny class contains only curve 398's Q-isomorphism "
            "class; there is no nontrivial rational isogeny to another curve."
        ),
        "isomorphism_class_count": 1,
        "method": "PARI ellisomat over Q",
        "minimal_isogeny_degree_matrix": "Mat(1)",
    }
    assert len(curve398["local_reduction"]) == 18
    assert all(row["conductor_exponent"] == 1 for row in curve398["local_reduction"])
    assert curve398["height_diagnostic"]["used_for_rank_claim"] is False
    assert curve398["construction_provenance"][
        "exact_fibration_parameter_and_section_map"
    ] == "NOT_PUBLIC_UNKNOWN"
    a1_test = curve398["construction_provenance"]["known_equation_explicit_a1_test"]
    assert a1_test["prime"] == 179
    assert a1_test["finite_roots"] == []
    assert a1_test["root_at_infinity"] is False
    assert curve398["generation"]["checker_sha256"] == EXPECTED_HASHES[
        "elliptic-curves/cas/verify_icarm_curve398_rank30.py"
    ]
    assert curve398["generation"]["model_data_sha256"] == EXPECTED_HASHES[
        "elliptic-curves/cas/icarm_curve398.py"
    ]
    assert curve398["generation"]["known_a1_reduction_sha256"] == EXPECTED_HASHES[
        "elliptic-curves/data/icarm_curve398_known_a1_mod179.json"
    ]

    print(
        "PASS: pinned ECR30/ECR31 and curve-398 artifacts and their exact source provenance "
        "match; stored point, matrix-shape, torsion, and claim-boundary fields "
        "are internally consistent"
    )
    print(
        "SCOPE: unconditional rank lower bounds 30 and 31 only; no exact-rank, "
        "BSD/GRH, height-regulator, or K3-family inference is accepted"
    )
    print(
        "NO RECOMPUTATION: no curve arithmetic, finite-group enumeration, "
        "matrix-rank calculation, PARI, Sage, or artifact rewrite"
    )


if __name__ == "__main__":
    main()
