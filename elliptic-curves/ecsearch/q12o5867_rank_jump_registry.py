#!/usr/bin/env python3
"""Exact-only registry support for q12o5867 specialization rank jumps.

Registry admission requires an exact finite-quotient certificate for at least
eighteen points.  Search scores, numerical height tests, raw point hits, and
bounded non-results are deliberately outside the schema.
"""

from __future__ import annotations

from fractions import Fraction
import gzip
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any, Iterable, Sequence


ELLIPTIC_ROOT = Path(__file__).resolve().parents[1]
CAS = ELLIPTIC_ROOT / "cas"
if str(CAS) not in sys.path:
    sys.path.insert(0, str(CAS))

from elliptic_candidate_record import weierstrass_invariants  # noqa: E402


Q = Fraction
REGISTRY_SCHEMA = "elliptic-curves.q12o5867-certified-rank-jump-registry.v1"


def canonical_json_sha256(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return sha256(encoded).hexdigest()


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_global_minimal_q_isomorphism_key(
    model: Sequence[Fraction | int | str],
) -> str:
    """Key a global minimal Q-isomorphism class by exact ``(c4,c6)``.

    Global minimal Q-isomorphic models differ by an admissible change with
    scale ``u=+/-1``.  Their exact ``c4`` and ``c6`` therefore agree, including
    at the exceptional ``j=0,1728`` values.  Twists do not share both values.
    """

    if len(model) != 5:
        raise ValueError("a global minimal Weierstrass model has five coefficients")
    coefficients = tuple(Q(value) for value in model)
    if any(value.denominator != 1 for value in coefficients):
        raise ValueError("the canonical key requires an integral global minimal model")
    invariants = weierstrass_invariants(coefficients)
    payload = {
        "normalization": "global-minimal-Q-model-exact-c4-c6-v1",
        "c4": str(invariants["c4"]),
        "c6": str(invariants["c6"]),
    }
    return "qmin-c4c6-sha256:" + canonical_json_sha256(payload)


def empty_registry() -> dict[str, Any]:
    return {
        "schema": REGISTRY_SCHEMA,
        "status": "EMPTY_NO_CERTIFIED_Q12O5867_RANK_JUMP",
        "entry_count": 0,
        "entries": [],
        "admission_rule": {
            "minimum_exact_certified_rank_lower_bound": 18,
            "minimum_exact_quotient_gain_beyond_generic_rank17": 1,
            "deduplication_key": "global-minimal exact (c4,c6) Q-isomorphism key",
            "repository_exclusion": (
                "reject a key already present in repository status/generated manifests"
            ),
        },
        "claim_boundary": [
            "Scores, numerical points, and uncertified point searches are never entries.",
            "An entry is a certified lower bound, not an exact-rank claim.",
        ],
    }


def validate_registry(registry: dict[str, Any]) -> None:
    if registry.get("schema") != REGISTRY_SCHEMA:
        raise ValueError("unexpected q12o5867 rank-jump registry schema")
    entries = registry.get("entries")
    if not isinstance(entries, list):
        raise ValueError("registry entries must be a list")
    if registry.get("entry_count") != len(entries):
        raise ValueError("registry entry_count disagrees with entries")
    keys = []
    for entry in entries:
        if int(entry.get("exact_certified_rank_lower_bound", 0)) < 18:
            raise ValueError("a registry entry lacks exact certified rank at least 18")
        if int(entry.get("exact_quotient_gain_beyond_generic_rank17", 0)) < 1:
            raise ValueError("a registry entry lacks exact quotient escape")
        key = entry.get("canonical_global_minimal_q_isomorphism_key")
        expected = canonical_global_minimal_q_isomorphism_key(
            entry["global_minimal_model"]
        )
        if key != expected:
            raise ValueError("a registry entry has a noncanonical Q-isomorphism key")
        keys.append(key)
    if len(keys) != len(set(keys)):
        raise ValueError("the registry contains a duplicate Q-isomorphism key")
    expected_status = (
        "EMPTY_NO_CERTIFIED_Q12O5867_RANK_JUMP"
        if not entries
        else "PASS_EXACT_CERTIFIED_Q12O5867_RANK_JUMP_REGISTRY"
    )
    if registry.get("status") != expected_status:
        raise ValueError("registry status disagrees with entry inventory")


def exact_admission_gate(
    *, certified_independent: bool, point_count: int, quotient_gain: int
) -> None:
    if not certified_independent:
        raise ValueError("registry admission requires an exact independence certificate")
    if point_count < 18:
        raise ValueError("registry admission requires at least 18 certified points")
    if quotient_gain < 1 or point_count < 17 + quotient_gain:
        raise ValueError("registry admission requires exact gain beyond generic rank 17")


def merge_exact_entry(registry: dict[str, Any], entry: dict[str, Any]) -> dict[str, Any]:
    """Insert or merge one already validated exact entry by canonical key."""

    validate_registry(registry)
    exact_admission_gate(
        certified_independent=bool(entry.get("certified_independent")),
        point_count=int(entry.get("exact_certified_rank_lower_bound", 0)),
        quotient_gain=int(entry.get("exact_quotient_gain_beyond_generic_rank17", 0)),
    )
    expected_key = canonical_global_minimal_q_isomorphism_key(
        entry["global_minimal_model"]
    )
    if entry.get("canonical_global_minimal_q_isomorphism_key") != expected_key:
        raise ValueError("new entry does not carry its canonical Q-isomorphism key")
    answer = json.loads(json.dumps(registry))
    old = next(
        (
            item
            for item in answer["entries"]
            if item["canonical_global_minimal_q_isomorphism_key"] == expected_key
        ),
        None,
    )
    if old is None:
        answer["entries"].append(entry)
    else:
        for field in ("parameters", "provenance"):
            merged = {
                canonical_json_sha256(value): value
                for value in (*old.get(field, []), *entry.get(field, []))
            }
            old[field] = [merged[key] for key in sorted(merged)]
        if int(entry["exact_certified_rank_lower_bound"]) > int(
            old["exact_certified_rank_lower_bound"]
        ):
            preserved_parameters = old.get("parameters", [])
            preserved_provenance = old.get("provenance", [])
            old.clear()
            old.update(entry)
            old["parameters"] = preserved_parameters
            old["provenance"] = preserved_provenance
    answer["entries"].sort(
        key=lambda item: item["canonical_global_minimal_q_isomorphism_key"]
    )
    answer["entry_count"] = len(answer["entries"])
    answer["status"] = "PASS_EXACT_CERTIFIED_Q12O5867_RANK_JUMP_REGISTRY"
    validate_registry(answer)
    return answer


def _walk_models(value: Any, parent_key: str = "") -> Iterable[Sequence[Any]]:
    if isinstance(value, dict):
        for key, child in value.items():
            lower = str(key).lower()
            if (
                isinstance(child, list)
                and len(child) == 5
                and ("model" in lower or "ainv" in lower)
            ):
                yield child
            yield from _walk_models(child, lower)
    elif isinstance(value, list):
        for child in value:
            yield from _walk_models(child, parent_key)


def documented_q_isomorphism_keys(
    paths: Iterable[Path], *, excluded_paths: Iterable[Path] = ()
) -> dict[str, list[str]]:
    """Conservatively inventory exact models in status/generated manifests."""

    excluded = {path.resolve() for path in excluded_paths}
    answer: dict[str, list[str]] = {}
    for path in paths:
        path = path.resolve()
        if path in excluded or not path.is_file():
            continue
        try:
            if path.suffix == ".gz":
                with gzip.open(path, "rt") as handle:
                    value = json.load(handle)
            elif path.suffix == ".json":
                value = json.loads(path.read_text())
            else:
                continue
        except (OSError, UnicodeDecodeError, json.JSONDecodeError):
            continue
        for model in _walk_models(value):
            try:
                key = canonical_global_minimal_q_isomorphism_key(model)
            except (ValueError, TypeError, ZeroDivisionError):
                continue
            answer.setdefault(key, []).append(str(path))
    return answer

