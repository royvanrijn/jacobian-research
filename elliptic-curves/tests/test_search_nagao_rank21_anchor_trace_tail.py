from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
SCRIPT = CAS / "search_nagao_rank21_anchor_trace_tail.py"
ARTIFACT = ROOT / "artifacts/generated-results/elliptic_nagao_rank21_anchor_trace_tail.json"
sys.path.insert(0, str(CAS))
SPEC = importlib.util.spec_from_file_location("search_nagao_rank21_anchor_trace_tail", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class AnchorTraceTailTests(unittest.TestCase):
    def test_parent_point_population_is_pinned(self) -> None:
        excluded, audit = MODULE.parent_point_exclusions()
        self.assertEqual(len(excluded), 28)
        self.assertEqual(audit["sha256"], MODULE.BROAD_ARTIFACT_SHA256)

    def test_generated_artifact(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("generated artifact not present")
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(data["target"]["rank_at_least"], 21)
        self.assertEqual(
            data["generation"]["exactly_deduplicated_proxy_survivors"],
            len(data["exact_b2000_population"]),
        )
        parent = set(data["parent_broad_search"]["constructor_parameters"])
        searched = set(data["point_population_selection"]["constructor_parameters"])
        self.assertFalse(parent & searched)
        for checkpoint in data["exact_checkpoints_stable_numerical_rank_at_least_18"]:
            self.assertEqual(checkpoint["exact_rank_certificate"]["status"], "certified")


if __name__ == "__main__":
    unittest.main()
