from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "elliptic-curves" / "ecsearch"))

from rank_jump_benchmark import (  # noqa: E402
    GROUP_FIELDS,
    SCHEMA,
    evaluate_manifest,
    validate_manifest,
)


FEATURE_NAMES = (
    "family_residual_s0_cumulative_B200",
    "family_residual_s5_cumulative_B200",
    "family_residual_s0_window_101_200",
    "family_residual_s5_window_101_200",
    "conductor_scaled_s0_B200",
    "root_number",
    "predicted_local_conductor",
    "quotient_escape",
)


def candidate(
    identifier: str, family: str, jump: int, signal: float, twist: str
) -> dict[str, object]:
    discovery = {
        name: signal if "root_number" not in name else (-1 if signal > 0 else 1)
        for name in FEATURE_NAMES
    }
    held = {name: value * 0.9 for name, value in discovery.items()}
    return {
        "id": identifier,
        "family": family,
        "root_shape": f"shape-{family}",
        "parametrization_component": f"component-{family}",
        "quadratic_twist_class": twist,
        "exceptional_quotient_rank_lower_bound": jump,
        "discovery_features": discovery,
        "held_forward_features": held,
    }


def manifest() -> dict[str, object]:
    # Synthetic values exercise only the evaluator's leakage logic.  They are
    # not elliptic-curve evidence or a rank-labelled research corpus.
    return {
        "schema": SCHEMA,
        "discovery_primes_disjoint_from_held_forward": True,
        "feature_names": list(FEATURE_NAMES),
        "candidates": [
            candidate("a-plus", "a", 10, 5.0, "twist-a"),
            candidate("a-control", "a", 0, -2.0, "twist-a-control"),
            candidate("b-plus", "b", 8, 4.0, "twist-b"),
            candidate("b-control", "b", 0, -1.0, "twist-b-control"),
            candidate("c-plus", "c", 6, 3.0, "twist-c"),
            candidate("c-control", "c", 0, -3.0, "twist-c-control"),
            candidate("d-plus", "d", 10, 6.0, "twist-d"),
            candidate("d-control", "d", 0, -4.0, "twist-d-control"),
        ],
    }


class RankJumpBenchmarkTests(unittest.TestCase):
    def test_manifest_requires_declared_residual_and_auxiliary_features(self) -> None:
        self.assertEqual(validate_manifest(manifest()), FEATURE_NAMES)
        incomplete = manifest()
        incomplete["feature_names"] = ["family_residual_s0_cumulative_B200"]
        with self.assertRaisesRegex(ValueError, "family_residual_s5"):
            validate_manifest(incomplete)

    def test_all_structural_protocols_exclude_the_held_out_group(self) -> None:
        result = evaluate_manifest(manifest(), (1, 2))
        self.assertEqual([item["group_field"] for item in result["protocols"]], list(GROUP_FIELDS))
        for protocol in result["protocols"]:
            field = protocol["group_field"]
            for fold in protocol["folds"]:
                self.assertTrue(
                    fold["leakage_check"]["training_excludes_entire_held_out_group"]
                )
                self.assertFalse(
                    fold["leakage_check"]["held_forward_features_used_for_fitting"]
                )
                held = fold["held_out_group"]
                training_ids = set(fold["training_candidate_ids"])
                for row in manifest()["candidates"]:
                    if row["id"] in training_ids:
                        self.assertNotEqual(row[field], held)

    def test_metric_is_top_k_enrichment_not_accuracy(self) -> None:
        result = evaluate_manifest(manifest(), (1,))
        family_fold = result["protocols"][0]["folds"][0]
        evaluation = family_fold["targets"]["6"]
        self.assertEqual(evaluation["status"], "evaluated")
        top = evaluation["rankings"]["held_forward_features"]["top_k_enrichment"][0]
        self.assertEqual(top["requested_k"], 1)
        self.assertIn("enrichment_over_random", top)
        self.assertNotIn("accuracy", top)


if __name__ == "__main__":
    unittest.main()
