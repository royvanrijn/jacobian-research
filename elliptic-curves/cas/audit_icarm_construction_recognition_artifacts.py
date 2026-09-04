#!/usr/bin/env python3
"""Audit committed ICARM construction-recognition ledgers without searching."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
EXPECTED_HASHES = {
    "elliptic-curves/cas/analyze_icarm_construction_fingerprints.py":
        "fdaa2de2aa8c312e131dd55760855bc1d701a851d6b2d506d04a6c93359b6b42",
    "elliptic-curves/scripts/discover_record_families.py":
        "dbb4ae2580c45d4a25800225e11301c51158c97253b3c3631054a15c60769101",
    "elliptic-curves/ecsearch/family_discovery.py":
        "d4766c7f1da8d3e14d09ab8173b1ddce181d6623cf09c0b767952ed319f2c34a",
    "elliptic-curves/ecsearch/conductor_engineering.py":
        "24f0da76e99728ccb36b662c97f79993de6727cd5ad58fb7a45fd1e4af12b48f",
    "elliptic-curves/cas/icarm_curve245_mestre.py":
        "0d8b6823ba1cb80de4251adb3ff944bc3465bff58ea2bf82d2e0691a0194f225",
    "elliptic-curves/cas/mestre_root_tuples.py":
        "83dd78f2ae1c652c5a9f33c46f2475c180371076fa149fc955cb653480496d75",
    "elliptic-curves/data/family-discovery/icarm_273_282_302.json":
        "ecf80809eb27c386ac92341c5e0d97b9a95d20a69bfe6fb0e82372a9385cca6b",
    "elliptic-curves/families/fermigier_mestre_rank12.json":
        "705351f0ee70c774e00521ddd37fe838f0320688293a14aeb4aca9a5008f369b",
    "archive/elliptic-curves/artifacts/generated-results/"
    "elliptic_mestre_root_tuple_scale_max200_census.json":
        "7270769007f9c130fce8b1813164373de9c6a5eb1c6d86cfe71b8c96fada161b",
    "archive/elliptic-curves/artifacts/generated-results/"
    "elliptic_mestre_root_tuple_scale_max300_census.json":
        "c5a68905977f059182efc1233e7301c039b2164b45bdfef1f8fd106b13d263ea",
    "artifacts/generated-results/elliptic-curves/"
    "icarm_7fff_zip_public_source_281_282_285_286.json":
        "b722a795491fa96755506b53c8261194ad1bcff28d253f5638a643dc11c12ac4",
    "artifacts/generated-results/elliptic-curves/"
    "icarm_7fff_zip_independence_analysis_v1.json":
        "cc5af39c55cc68ab291b3c318b171ae8cbe987f47a1957a94d60098cd8833af8",
    "artifacts/generated-results/elliptic-curves/"
    "icarm_construction_fingerprints_v1.json":
        "66d04fe6fb95eba5c1049ea8337cdfa68bdb10aec0d99aa55b259971912bebfc",
    "artifacts/generated-results/elliptic-curves/"
    "icarm_construction_fingerprints_v2.json":
        "ced63ed67c61bb23484039237259127ffd0864426ae41429cd005e6989bfdc4a",
    "artifacts/generated-results/elliptic-curves/"
    "icarm_273_282_302_family_discovery_v1.json":
        "5d19571dc74f9e8c270bcaaa943f19e3876bb2a13ffc33c37d479a3206fd8770",
}


def load(relative_path: str) -> dict[str, object]:
    return json.loads((ROOT / relative_path).read_text(encoding="utf-8"))


def exact_matches(payload: dict[str, object]) -> dict[int, list[dict[str, object]]]:
    return {
        row["curve_id"]: row["exact_j_matches"]
        for row in payload["six_root_mestre_recognition"]["targets"]
    }


def main() -> None:
    for relative_path, expected_hash in EXPECTED_HASHES.items():
        actual_hash = hashlib.sha256((ROOT / relative_path).read_bytes()).hexdigest()
        assert actual_hash == expected_hash, (
            relative_path,
            actual_hash,
            expected_hash,
        )

    v1 = load(
        "artifacts/generated-results/elliptic-curves/"
        "icarm_construction_fingerprints_v1.json"
    )
    v2 = load(
        "artifacts/generated-results/elliptic-curves/"
        "icarm_construction_fingerprints_v2.json"
    )
    assert v1["schema"] == "elliptic-curves.icarm-construction-fingerprints.v1"
    assert v2["schema"] == "elliptic-curves.icarm-construction-fingerprints.v2"
    assert sorted(map(int, v1["targets"])) == [273, 281, 282, 285, 286]
    assert sorted(map(int, v2["targets"])) == [273, 281, 282, 285, 286, 302]
    for payload in (v1, v2):
        recognition = payload["six_root_mestre_recognition"]
        assert recognition["census_family_count"] == 2330
        assert recognition["diameter_at_most_300_family_count"] == 2329
        assert "No-match results exclude only this fixed-root Mestre census" in (
            recognition["boundary"]
        )
        matches = exact_matches(payload)
        assert len(matches[282]) == 1
        assert matches[282][0]["roots"] == [0, 29, 658, 722, 981, 1036]
        assert matches[282][0]["parameter_T"] == "11671/21"
        assert all(not rows for curve_id, rows in matches.items() if curve_id != 282)
        assert "necessary but not sufficient" in payload["repository_model_scan"][
            "interpretation"
        ]
        assert "does not exclude an isogenous quotient" in payload[
            "forced_torsion_exclusion"
        ]["boundary"]

    discovery = load(
        "artifacts/generated-results/elliptic-curves/"
        "icarm_273_282_302_family_discovery_v1.json"
    )
    assert discovery["schema"] == "elliptic-curves.generated-family-discovery.v1"
    assert discovery["claim_level"] == "complete bounded exact computation"
    assert discovery["specification_sha256"] == EXPECTED_HASHES[
        "elliptic-curves/data/family-discovery/icarm_273_282_302.json"
    ]
    construction = discovery["construction_space"]
    assert construction["polynomial_weierstrass_family_count"] == 1
    assert construction["six_root_mestre_family_count"] == 2333
    assert construction["total_family_count"] == 2334
    assert construction["six_root_generation_summary"] == {
        "duplicate_family_emission_count": 6,
        "emitted_family_count_before_deduplication": 2339,
        "generator_degenerate_parameter_pairs": 15,
        "generator_parameter_pairs_tested": 25,
    }
    rows = {row["target"]: row for row in discovery["targets"]}
    assert set(rows) == {"ICARM curve 273", "ICARM curve 282", "ICARM curve 302"}
    assert [rows[label]["exact_factorization_survivor_count"] for label in rows] == [
        113,
        114,
        146,
    ]
    assert not rows["ICARM curve 273"]["q_isomorphism_matches"]
    assert not rows["ICARM curve 302"]["q_isomorphism_matches"]
    matches282 = rows["ICARM curve 282"]["q_isomorphism_matches"]
    assert len(matches282) == 2
    assert {row["parameter"] for row in matches282} == {"11671/42", "11671/21"}
    assert {row["q_isomorphism_invariant_scale"] for row in matches282} == {
        "882",
        "147",
    }
    assert "declared bounded generators" in discovery["interpretation"]["boundary"]

    fingerprint_source = (
        ROOT / "elliptic-curves/cas/analyze_icarm_construction_fingerprints.py"
    ).read_text(encoding="utf-8")
    assert 'GENERATED.rglob("*.json")' in fingerprint_source

    print(
        "PASS: committed fixed-root and generated-space ICARM recognition "
        "ledgers and their declared inputs match exact hashes and bounded scopes"
    )
    print(
        "SCOPE: modular no-root exits are one-sided exact exclusions and every "
        "survivor is retained for exact factorization; no-match results apply "
        "only to the declared 2,330/2,334 family spaces"
    )
    print(
        "CORRECTION: the historical repository-model diagnostic scanned only "
        "uncompressed artifacts/generated-results/**/*.json, not gzip files, "
        "the archive, or every repository model"
    )
    print("NO RECOMPUTATION: no family generation, modular sieve, or factorization")


if __name__ == "__main__":
    main()
