from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "elliptic-curves/cas/verify_icarm_curve394_rank21.py"
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "icarm_curve394_rank21_v1.json"
)


def load_checker():
    spec = importlib.util.spec_from_file_location("icarm394", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class IcarmCurve394Rank21Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.checker = load_checker()
        cls.certificate = cls.checker.build_certificate()

    def test_exact_rank_lower_bound_and_conductor(self) -> None:
        record = self.certificate
        self.assertEqual(record["rank_lower_bound"]["certified_rank"], 21)
        self.assertTrue(
            record["compact_specialization"]["exact_public_minimal_model_equality"]
        )
        self.assertEqual(record["public_point_replay"]["count"], 21)
        self.assertEqual(record["curve"]["root_number"], -1)
        self.assertEqual(len(record["curve"]["local_reductions"]), 12)
        self.assertTrue(record["curve"]["exact_log_bound"]["strict_target_proved_exactly"])

    def test_pinned_artifact_matches_checker(self) -> None:
        expected = json.dumps(self.certificate, indent=2, sort_keys=True) + "\n"
        self.assertEqual(ARTIFACT.read_text(), expected)


if __name__ == "__main__":
    unittest.main()
