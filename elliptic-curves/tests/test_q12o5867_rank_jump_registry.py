#!/usr/bin/env python3
"""Tests for the exact-only q12o5867 rank-jump registry."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = ROOT / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
sys.path.insert(0, str(ELLIPTIC_ROOT))
sys.path.insert(0, str(CAS))

from ecsearch.q12o5867_rank_jump_registry import (  # noqa: E402
    canonical_global_minimal_q_isomorphism_key,
    documented_q_isomorphism_keys,
    empty_registry,
    exact_admission_gate,
    merge_exact_entry,
    validate_registry,
)
from elliptic_candidate_record import (  # noqa: E402
    WeierstrassChange,
    change_weierstrass_model,
)


class Q12O5867RankJumpRegistryTests(unittest.TestCase):
    def test_empty_registry_and_exact_admission_gate(self) -> None:
        registry = empty_registry()
        validate_registry(registry)
        with self.assertRaises(ValueError):
            exact_admission_gate(
                certified_independent=False, point_count=18, quotient_gain=1
            )
        with self.assertRaises(ValueError):
            exact_admission_gate(
                certified_independent=True, point_count=17, quotient_gain=0
            )
        exact_admission_gate(
            certified_independent=True, point_count=18, quotient_gain=1
        )

    def test_key_is_invariant_under_unit_scale_minimal_change(self) -> None:
        source = tuple(Fraction(value) for value in (1, -1, 1, -10, 20))
        target = change_weierstrass_model(
            source, WeierstrassChange(Fraction(-1), Fraction(3), Fraction(2), Fraction(5))
        )
        self.assertTrue(all(value.denominator == 1 for value in target))
        self.assertEqual(
            canonical_global_minimal_q_isomorphism_key(source),
            canonical_global_minimal_q_isomorphism_key(target),
        )

    def test_merge_deduplicates_by_q_isomorphism_key(self) -> None:
        model = [1, -1, 1, -10, 20]
        key = canonical_global_minimal_q_isomorphism_key(model)

        def entry(parameter: str) -> dict[str, object]:
            return {
                "canonical_global_minimal_q_isomorphism_key": key,
                "global_minimal_model": model,
                "parameters": [{"affine_value": parameter}],
                "provenance": [{"reproducing_command": f"probe {parameter}"}],
                "certified_independent": True,
                "exact_certified_rank_lower_bound": 18,
                "exact_quotient_gain_beyond_generic_rank17": 1,
            }

        registry = merge_exact_entry(empty_registry(), entry("1/2"))
        registry = merge_exact_entry(registry, entry("2/3"))
        self.assertEqual(registry["entry_count"], 1)
        self.assertEqual(len(registry["entries"][0]["parameters"]), 2)
        validate_registry(registry)

    def test_documented_model_inventory(self) -> None:
        model = [1, -1, 1, -10, 20]
        key = canonical_global_minimal_q_isomorphism_key(model)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            path.write_text(json.dumps({"global_minimal_model": model}))
            inventory = documented_q_isomorphism_keys((path,))
            excluded = documented_q_isomorphism_keys(
                (path,), excluded_paths=(path,)
            )
        self.assertEqual(inventory[key], [str(path.resolve())])
        self.assertEqual(excluded, {})


if __name__ == "__main__":
    unittest.main()
