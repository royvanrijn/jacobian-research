#!/usr/bin/env python3
"""Verify the pinned low-degree census ledgers.

The default mode performs the historical full deterministic replay, including
the Singular and SymPy coefficient-ideal calculations.  The
``--audit-existing-only`` mode is deliberately solver-free: it checks the
manifest, exact support/orbit labels, stage-to-stage routing, and the internal
arithmetic of the committed records without recomputing any Groebner basis.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from jcsearch.low_degree_pipeline import (  # noqa: E402
    MANIFEST_FILENAME,
    STAGE_FILENAMES,
    build_low_degree_census,
)
from jcsearch.low_degree_census import (  # noqa: E402
    Support,
    degree_profiles_below,
    determinant3,
    profile_rank_function,
    sha256_json,
)


ARTIFACT_ROOT = ROOT / "artifacts/generated-results"
PINNED_PARAMETERS = {
    "max_degree": 7,
    "target_profile": [7, 6, 4],
    "max_nonlinear_support": 6,
    "primes": [11, 13, 17],
}
STAGE_SCHEMAS = (
    "global-low-degree-census.profiles.v1",
    "global-low-degree-census.supports.v1",
    "global-low-degree-census.buckets.v1",
    "global-low-degree-census.valuations.v1",
    "global-low-degree-census.smt.v1",
    "global-low-degree-census.modular.v1",
    "global-low-degree-census.exact.v1",
    "global-low-degree-census.boundary.v1",
)


def _unique_ids(rows: list[dict[str, object]], label: str) -> tuple[str, ...]:
    identifiers = tuple(str(row["support_id"]) for row in rows)
    assert len(identifiers) == len(set(identifiers)), f"duplicate {label} support id"
    return identifiers


def _support_from_json(payload: object) -> Support:
    assert isinstance(payload, list) and len(payload) == 3
    rows = tuple(
        tuple(tuple(int(value) for value in exponent) for exponent in row)
        for row in payload
    )
    return Support(rows)  # type: ignore[arg-type]


def _load_pinned_artifacts(
    artifact_root: Path,
) -> tuple[dict[str, object], dict[str, dict[str, object]]]:
    manifest_path = artifact_root / MANIFEST_FILENAME
    manifest = json.loads(manifest_path.read_text())
    assert manifest["schema"] == "global-low-degree-census.manifest.v1"
    assert manifest["generator"] == "scripts/compile_global_low_degree_census.py"
    assert manifest["parameters"] == PINNED_PARAMETERS
    assert tuple(manifest["stage_sha256"]) == STAGE_FILENAMES
    assert "cardinality-unbounded support census remains open" in manifest["claim_boundary"]

    stages: dict[str, dict[str, object]] = {}
    for filename, schema in zip(STAGE_FILENAMES, STAGE_SCHEMAS):
        path = artifact_root / filename
        content = path.read_bytes()
        actual_hash = hashlib.sha256(content).hexdigest()
        assert actual_hash == manifest["stage_sha256"][filename], filename
        payload = json.loads(content)
        assert payload["schema"] == schema, filename
        assert payload["parameters"] == {
            "dimension": 3,
            "maximum_coordinate_degree": 7,
            "target_invariant_degree_profile": [7, 6, 4],
            "maximum_nonlinear_support": 6,
            "collision_axis": 1,
            "normalization": ["F(0)=0", "JF(0)=I", "F(e1)=0"],
        }, filename
        stages[filename] = payload
    return manifest, stages


def audit_existing_only(artifact_root: Path) -> None:
    """Audit committed records and their routing without invoking a solver."""

    _manifest, stages = _load_pinned_artifacts(artifact_root)
    profiles = stages[STAGE_FILENAMES[0]]
    supports = stages[STAGE_FILENAMES[1]]
    buckets = stages[STAGE_FILENAMES[2]]
    valuations = stages[STAGE_FILENAMES[3]]
    signs = stages[STAGE_FILENAMES[4]]
    modular = stages[STAGE_FILENAMES[5]]
    exact = stages[STAGE_FILENAMES[6]]
    boundary = stages[STAGE_FILENAMES[7]]

    expected_profiles = degree_profiles_below((7, 6, 4))
    assert profiles["profile_count"] == len(expected_profiles) == 74
    assert profiles["profiles"] == [list(profile) for profile in expected_profiles]
    assert profiles["profile_rank_gates"] == [
        {
            "profile": list(profile),
            "rank_above_degree": {
                str(threshold): profile_rank_function(profile, threshold)
                for threshold in range(8)
            },
        }
        for profile in expected_profiles
    ]

    orbit_rows: list[dict[str, object]] = []
    representative_by_id: dict[str, Support] = {}
    labelled_counts: dict[str, int] = {}
    orbit_counts: dict[str, int] = {}
    for size in range(1, 7):
        rows = supports["orbits"][str(size)]
        assert isinstance(rows, list)
        size_labelled_count = 0
        for row in rows:
            support = _support_from_json(row["support"])
            assert support.nonlinear_size == size
            assert support == support.canonical_under_collision_stabilizer()
            assert row["support_id"] == support.identifier
            swapped = support.swapped_23()
            expected_members = sorted({support.identifier, swapped.identifier})
            assert row["member_ids"] == expected_members
            assert row["orbit_size"] == len(expected_members)
            assert row["second_member_rule"] == (
                None
                if len(expected_members) == 1
                else "simultaneous x2<->x3, F2<->F3"
            )
            assert support.identifier not in representative_by_id
            representative_by_id[support.identifier] = support
            size_labelled_count += len(expected_members)
            orbit_rows.append(row)
        labelled_counts[str(size)] = size_labelled_count
        orbit_counts[str(size)] = len(rows)

    representative_ids = _unique_ids(orbit_rows, "support-orbit")
    assert supports["determinant_balanced_supports_by_size"] == labelled_counts == {
        "1": 0,
        "2": 0,
        "3": 0,
        "4": 30,
        "5": 85,
        "6": 1694,
    }
    assert supports["determinant_balanced_orbits_by_size"] == orbit_counts == {
        "1": 0,
        "2": 0,
        "3": 0,
        "4": 15,
        "5": 47,
        "6": 851,
    }
    assert len(representative_ids) == 913
    assert supports["full_boolean_support_space"]["optional_support_atoms"] == 348

    bucket_rows = buckets["representatives"]
    assert isinstance(bucket_rows, list)
    assert _unique_ids(bucket_rows, "bucket") == representative_ids
    assert buckets["representative_count"] == len(bucket_rows) == 913
    for row in bucket_rows:
        support = representative_by_id[str(row["support_id"])]
        stored_buckets = row["buckets"]
        assert isinstance(stored_buckets, list) and stored_buckets
        assert row["bucket_sha256"] == sha256_json(stored_buckets)
        assert row["bucket_count"] == len(stored_buckets)
        exponents = [tuple(bucket["exponent"]) for bucket in stored_buckets]
        assert exponents == sorted(exponents, key=lambda exponent: (sum(exponent), exponent))
        assert len(exponents) == len(set(exponents))
        nonconstant = [bucket for bucket in stored_buckets if bucket["exponent"] != [0, 0, 0]]
        histogram = Counter(str(len(bucket["terms"])) for bucket in nonconstant)
        assert row["nonconstant_bucket_count"] == len(nonconstant)
        assert row["singleton_bucket_count"] == sum(
            len(bucket["terms"]) == 1 for bucket in nonconstant
        ) == 0
        assert row["contribution_count_histogram"] == dict(
            sorted(histogram.items(), key=lambda item: int(item[0]))
        )
        for bucket in stored_buckets:
            bucket_exponent = tuple(bucket["exponent"])
            assert bucket["terms"]
            for term in bucket["terms"]:
                alpha, beta, gamma = (
                    tuple(term[name]) for name in ("alpha", "beta", "gamma")
                )
                assert alpha in support.full_rows[0]
                assert beta in support.full_rows[1]
                assert gamma in support.full_rows[2]
                assert term["multiplier"] == determinant3(alpha, beta, gamma) != 0
                assert tuple(
                    alpha[index] + beta[index] + gamma[index] - 1
                    for index in range(3)
                ) == bucket_exponent

    valuation_rows = valuations["representatives"]
    assert isinstance(valuation_rows, list)
    assert _unique_ids(valuation_rows, "valuation") == representative_ids
    valuation_histogram = Counter()
    for row in valuation_rows:
        assert row["valuation_class_count"] == len(row["classes"])
        assert row["classes_sha256"] == sha256_json(row["classes"])
        valuation_histogram[str(len(row["classes"]))] += 1
    assert valuations["representative_count_by_valuation_class_count"] == dict(
        sorted(valuation_histogram.items(), key=lambda item: int(item[0]))
    )
    assert valuations["supports_without_candidate_valuation"] == [
        row["support_id"] for row in valuation_rows if not row["classes"]
    ]

    sign_rows = signs["representatives"]
    assert isinstance(sign_rows, list)
    assert _unique_ids(sign_rows, "sign-SMT") == representative_ids
    assert signs["status_counts"] == dict(
        sorted(Counter(row["status"] for row in sign_rows).items())
    )
    assert all(
        row["scope"] == "necessary over ordered coefficient fields; not a complex-field exclusion"
        for row in sign_rows
    )

    assert modular["primes"] == [11, 13, 17]
    assert "routing evidence only" in modular["logical_status"]
    for prime in modular["primes"]:
        rows = modular["representatives_by_prime"][str(prime)]
        assert _unique_ids(rows, f"mod-{prime}") == representative_ids
        assert all(row["field"] == f"F_{prime}" for row in rows)
        summary = {
            "unit_ideal": sum(bool(row["unit_ideal"]) for row in rows),
            "isolated_or_zero_dimensional": sum(
                not row["unit_ideal"] and row["dimension"] == 0 for row in rows
            ),
            "positive_dimensional": sum(
                not row["unit_ideal"]
                and row["dimension"] is not None
                and row["dimension"] > 0
                for row in rows
            ),
        }
        assert modular["summary"][str(prime)] == summary
        assert summary["unit_ideal"] == 913

    exact_rows = exact["singular_results"]
    assert isinstance(exact_rows, list)
    assert _unique_ids(exact_rows, "exact-QQ") == representative_ids
    assert all(row["field"] == "QQ" and row["unit_ideal"] for row in exact_rows)
    assert exact["representative_count"] == len(exact_rows) == 913
    assert exact["unit_ideal_count"] == 913
    assert exact["surviving_support_ids"] == []
    assert exact["dense_quadratic_collision_ideal"]["unit_ideal_over_QQ"] is True
    assert exact["dense_quadratic_collision_ideal"]["basis"] == ["1"]
    assert exact["sparse_frontier_theorem"] == {
        "lower_bound_nonlinear_support": 7,
        "attainment_proved": False,
        "scope": (
            "Every normalized characteristic-zero Keller collision of raw coordinate "
            "degree at most 7 has at least 7 nonlinear monomial occurrences. Hence "
            "the same lower bound holds for every invariant degree profile below "
            "(7, 6, 4). No support of size 7 is constructed or asserted to exist."
        ),
    }
    assert exact["completely_eliminated_profiles"] == [
        [1, 1, 1],
        [2, 1, 1],
        [2, 2, 1],
        [2, 2, 2],
    ]

    assert boundary["coefficient_boundary"]["downward_support_boundary_complete"] is True
    assert boundary["supports_requiring_component_boundary_audit"] == []
    assert boundary["projective_charts_requiring_audit_after_exact_lifting"] == 0
    assert boundary["survives_projective_boundary_analysis"] is False
    assert "support-at-most-6 stratum" in boundary["scope_warning"]

    print(
        "PASS committed global low-degree census audit: manifest and 8 stages; "
        "913 support ids route unchanged through valuation, sign, modular, and QQ records; "
        "no solver replay"
    )


def full_replay(artifact_root: Path) -> None:
    expected = build_low_degree_census(
        progress=lambda message: print(message, flush=True)
    )
    manifest_path = artifact_root / MANIFEST_FILENAME
    manifest = json.loads(manifest_path.read_text())
    assert manifest["schema"] == "global-low-degree-census.manifest.v1"
    assert tuple(manifest["stage_sha256"]) == STAGE_FILENAMES

    for filename in STAGE_FILENAMES:
        path = artifact_root / filename
        content = path.read_bytes()
        actual_hash = hashlib.sha256(content).hexdigest()
        assert actual_hash == manifest["stage_sha256"][filename], filename
        assert json.loads(content) == expected[filename], filename

    profiles = expected[STAGE_FILENAMES[0]]
    supports = expected[STAGE_FILENAMES[1]]
    exact = expected[STAGE_FILENAMES[6]]
    boundary = expected[STAGE_FILENAMES[7]]
    assert profiles["profile_count"] == 74
    assert len(profiles["profile_rank_gates"]) == 74
    assert profiles["profile_rank_gates"][-1] == {
        "profile": [7, 6, 3],
        "rank_above_degree": {
            "0": 3,
            "1": 3,
            "2": 3,
            "3": 2,
            "4": 2,
            "5": 2,
            "6": 1,
            "7": 0,
        },
    }
    assert supports["determinant_balanced_supports_by_size"] == {
        "1": 0,
        "2": 0,
        "3": 0,
        "4": 30,
        "5": 85,
        "6": 1694,
    }
    assert supports["determinant_balanced_orbits_by_size"] == {
        "1": 0,
        "2": 0,
        "3": 0,
        "4": 15,
        "5": 47,
        "6": 851,
    }
    buckets = expected[STAGE_FILENAMES[2]]
    assert all(row["buckets"] for row in buckets["representatives"])
    assert all(
        row["bucket_sha256"] == sha256_json(row["buckets"])
        for row in buckets["representatives"]
    )
    assert exact["representative_count"] == 913
    assert exact["unit_ideal_count"] == 913
    assert exact["surviving_support_ids"] == []
    assert exact["sparse_frontier_theorem"]["lower_bound_nonlinear_support"] == 7
    assert exact["sparse_frontier_theorem"]["attainment_proved"] is False
    assert exact["completely_eliminated_profiles"] == [
        [1, 1, 1],
        [2, 1, 1],
        [2, 2, 1],
        [2, 2, 2],
    ]
    assert not boundary["survives_projective_boundary_analysis"]
    print(
        "PASS global low-degree census: 74 profiles; sparse supports "
        "30/85/1694; 913 exact orbit ideals are units; support lower bound seven"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifacts-dir", type=Path, default=ARTIFACT_ROOT)
    parser.add_argument(
        "--audit-existing-only",
        action="store_true",
        help="validate committed records and routing without Singular, SymPy, or Z3",
    )
    args = parser.parse_args()
    if args.audit_existing_only:
        audit_existing_only(args.artifacts_dir)
    else:
        full_replay(args.artifacts_dir)


if __name__ == "__main__":
    main()
