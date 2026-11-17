#!/usr/bin/env python3
"""Audit the committed bounded HC4 finite-field search ledgers.

This is a provenance and coverage check only.  It does not enumerate
potentials, expand a Hessian determinant, invoke Singular, or promote a
finite-field search to a characteristic-zero statement.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "artifacts" / "generated-results"
ARTIFACT_HASHES = {
    "hc4_finite_field_sparse_search.json": (
        "2658309db3ecad6df24a3ab5c80d9b0b69499b30c64ce6582653d32d688fd77c"
    ),
    "hc4_finite_field_dense_support_search.json": (
        "be0f1f6fc353d418e34abc8a73538f316bff54056ef16790a3830279e4de9409"
    ),
    "hc4_finite_field_axis_support_search.json": (
        "a90543c537000f74e69f8ffa5586296678aa90f460a3884d8a6d48867420b495"
    ),
    "hc4_finite_field_cone_bridge_search.json": (
        "4f6ebf4bc13d4fb8c9e6dc87f61142f006fe7c82a187ab5a7bbe21559783794a"
    ),
    "hc4_finite_field_oblique_cone_bridge_search.json": (
        "fe78703eeeff9988acb574cc76fa09538df95d39eaf1c8c629e8bcf7100feadf"
    ),
}
SOURCE_HASHES = {
    "scripts/search_hc4_finite_field_potentials.py": (
        "f3b76a0f32fcca96f75abce9c4fdc12641178593a7bc23a8cb20d5c8a1b727ee"
    ),
    "scripts/search_hc4_finite_field_dense_supports.py": (
        "783b4c836185012b149e8c3a9b69dfed94139c04de6f176bce11a140959f0dd6"
    ),
    "scripts/search_hc4_oblique_cone_bridges.py": (
        "d88f601439befe293ab6c8e9ff9d616f61e8c90750707251c07b6abe05db4361"
    ),
}
DEGREES = [5, 6, 7, 8]
PRIMES = [11, 13]
SUPPORT_SIZES = [6, 8, 10, 12]
NORMALIZATION = {
    "collision": "grad(psi)(0) = grad(psi)(1,0,0,0)",
    "constant_hessian_determinant": 1,
    "quadratic_part": "x0*x1 + x2*x3",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_artifact(name: str) -> dict[str, Any]:
    path = GENERATED / name
    actual = digest(path)
    assert actual == ARTIFACT_HASHES[name], (
        f"{name} drifted: expected {ARTIFACT_HASHES[name]}, got {actual}"
    )
    return json.loads(path.read_text(encoding="utf-8"))


def verify_internal_hash(payload: dict[str, Any], field: str) -> None:
    copy = json.loads(json.dumps(payload))
    recorded = copy.pop(field)
    if field == "deterministic_content_sha256":
        for record in copy["records"]:
            record.pop("solver_seconds", None)
    canonical = json.dumps(copy, sort_keys=True, separators=(",", ":"))
    assert hashlib.sha256(canonical.encode()).hexdigest() == recorded


def verify_sparse(payload: dict[str, Any]) -> int:
    assert payload["status"] == "bounded finite-field experiment; not a proof"
    assert payload["normalization"] == NORMALIZATION
    assert payload["supports_occurring_at_multiple_primes"] == []
    verify_internal_hash(payload, "content_sha256_before_hash_field")

    expected_direction_counts = {5: 107, 6: 191, 7: 311, 8: 476}
    records = payload["records"]
    assert {(row["prime"], row["degree_bound"]) for row in records} == set(
        itertools.product(PRIMES, DEGREES)
    )
    total = 0
    for row in records:
        prime = row["prime"]
        degree = row["degree_bound"]
        directions = expected_direction_counts[degree]
        expected_count = directions * (prime - 1) + (
            directions * (directions - 1) // 2
        ) * (prime - 1) ** 2
        assert row["direction_count"] == directions
        assert row["support_bound"] == 2
        assert len(row["deterministic_points"]) == 6
        assert row["potentials_evaluated"] == expected_count
        assert 0 <= row["evaluation_survivors"] <= expected_count
        assert row["exact_candidate_count"] == 0
        assert row["exact_candidates"] == []
        total += expected_count
    assert total == 45_181_194
    return total


def verify_sampled(
    payload: dict[str, Any], strategies: list[str], trials: int
) -> tuple[int, int]:
    assert "bounded" in payload["status"] and "not a proof" in payload["status"]
    assert payload["degrees"] == DEGREES
    assert payload["primes"] == PRIMES
    assert payload["support_sizes"] == SUPPORT_SIZES
    assert payload["strategies"] == strategies
    assert payload["trials"] == trials
    assert payload["normalization"] == NORMALIZATION
    verify_internal_hash(payload, "deterministic_content_sha256")

    expected_keys = set(
        itertools.product(DEGREES, PRIMES, SUPPORT_SIZES, strategies, range(trials))
    )
    rows = payload["records"]
    actual_keys = {
        (
            row["degree_bound"],
            row["prime"],
            row["support_size"],
            row["strategy"],
            row["trial"],
        )
        for row in rows
    }
    assert actual_keys == expected_keys
    assert len(rows) == len(expected_keys)
    assert payload["status_counts"] == {"unit": len(rows)}
    for row in rows:
        support = row["support_indices"]
        assert support == sorted(support)
        assert len(support) == len(set(support)) == row["support_size"]
        assert row["prime"] > row["degree_bound"]
        assert row["solver_status"] == "unit"
    return len(rows) // len(PRIMES), len(rows)


def verify_oblique(payload: dict[str, Any]) -> tuple[int, int]:
    assert payload["status"] == "bounded oblique-cone finite-field experiment; not a proof"
    assert payload["degrees"] == DEGREES
    assert payload["primes"] == PRIMES
    assert payload["support_sizes"] == SUPPORT_SIZES
    assert payload["slopes"] == ["-1", "1", "2"]
    assert payload["trials"] == 3
    assert {
        key: payload["normalization"][key] for key in NORMALIZATION
    } == NORMALIZATION
    assert payload["normalization"]["top_cone"] == "u=x2+slope*x3"
    verify_internal_hash(payload, "deterministic_content_sha256")

    expected_keys = set(
        itertools.product(
            DEGREES, payload["slopes"], PRIMES, SUPPORT_SIZES, range(3)
        )
    )
    rows = payload["records"]
    actual_keys = {
        (
            row["degree_bound"],
            row["slope"],
            row["prime"],
            row["support_size"],
            row["trial"],
        )
        for row in rows
    }
    assert actual_keys == expected_keys
    assert len(rows) == len(expected_keys)
    assert payload["status_counts"] == {"unit": len(rows)}
    for row in rows:
        support = row["support_indices"]
        assert support == sorted(support)
        assert len(support) == len(set(support)) == row["support_size"]
        assert row["top_direction_count"] == row["support_size"] // 2
        assert row["bridge_direction_count"] == (
            row["support_size"] - row["top_direction_count"]
        )
        assert row["solver_status"] == "unit"
    return len(rows) // len(PRIMES), len(rows)


def main() -> None:
    for relative, expected in SOURCE_HASHES.items():
        actual = digest(ROOT / relative)
        assert actual == expected, (
            f"{relative} drifted: expected {expected}, got {actual}"
        )

    sparse_total = verify_sparse(
        load_artifact("hc4_finite_field_sparse_search.json")
    )
    family_total = 0
    ideal_total = 0
    for name, strategies, trials in (
        (
            "hc4_finite_field_dense_support_search.json",
            ["uniform", "homogeneous", "mixed"],
            2,
        ),
        ("hc4_finite_field_axis_support_search.json", ["axis"], 2),
        (
            "hc4_finite_field_cone_bridge_search.json",
            ["cone2", "cone3"],
            4,
        ),
    ):
        families, ideals = verify_sampled(load_artifact(name), strategies, trials)
        family_total += families
        ideal_total += ideals
    families, ideals = verify_oblique(
        load_artifact("hc4_finite_field_oblique_cone_bridge_search.json")
    )
    family_total += families
    ideal_total += ideals
    assert family_total == 400
    assert ideal_total == 800

    print(
        "HC4_FINITE_FIELD_COMMITTED_AUDIT_PASS "
        f"sparse_choices={sparse_total} sampled_families={family_total} "
        f"unit_ideals={ideal_total}; bounded GF(11)/GF(13) evidence only; "
        "no enumeration or Singular replay"
    )


if __name__ == "__main__":
    main()
