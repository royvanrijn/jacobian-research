from __future__ import annotations

from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
SCRIPT = CAS / "search_nagao_rank21_anchor_union.py"
ARTIFACT = ROOT / "artifacts/generated-results/elliptic_nagao_rank21_anchor_union.json"
sys.path.insert(0, str(CAS))
SPEC = importlib.util.spec_from_file_location("search_nagao_rank21_anchor_union", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class AnchorUnionTests(unittest.TestCase):
    def test_anchor_parameters(self) -> None:
        self.assertEqual(
            tuple(anchor.parameter for anchor in MODULE.ANCHORS),
            (
                Fraction(14721, 188),
                Fraction(1393, 108),
                Fraction(1649, 6),
                Fraction(6629, 174),
            ),
        )

    def test_trace_fingerprints_replay(self) -> None:
        tables = MODULE.build_residue_tables(200)
        for anchor in MODULE.ANCHORS:
            self.assertEqual(
                MODULE.learn_trace_fingerprint(anchor, tables),
                MODULE.EXPECTED_TRACE_FINGERPRINTS[anchor.label],
            )

    def test_prior_population_manifest_is_pinned(self) -> None:
        excluded, audit = MODULE.prior_parameter_exclusions()
        self.assertTrue(all(abs(anchor.parameter) in excluded for anchor in MODULE.ANCHORS))
        self.assertEqual(len(audit["pinned_artifacts"]), 6)

    def test_generated_artifact(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("generated artifact not present")
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(data["target"]["rank_at_least"], 21)
        self.assertEqual(data["target"]["strict_log_conductor_upper_bound"], "182.72")
        self.assertEqual(
            len(data["exact_b2000_population"]),
            min(data["generation"]["retained_after_stratified_cap"], MODULE.B200_KEEP),
        )
        self.assertEqual(
            set(map(int, data["root_ball_unions"])),
            set(MODULE.ROOT_PRIMES),
        )
        for checkpoint in data["exact_checkpoints_stable_numerical_rank_at_least_18"]:
            self.assertGreaterEqual(checkpoint["stable_numerical_rank"], 18)
            self.assertEqual(checkpoint["exact_rank_certificate"]["status"], "certified")


if __name__ == "__main__":
    unittest.main()
