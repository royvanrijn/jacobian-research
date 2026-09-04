from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
TRUTH = ART / "icarm_curve398_hidden_a1_mw16_v1.json"
INPUT = ROOT / "elliptic-curves/data/icarm_curve398_mw16_blind_input_v1.json"
BLIND = ART / "curve398_mw16_adaptive_half_lattice_blind_v1.json"
VERIFIED = ART / "curve398_mw16_adaptive_half_lattice_verification_v1.json"
SCREEN = ROOT / "artifacts/generated-results/elkies-k3-curve398-11952-norm8-a1-modular-screen-v1.json"
SURVIVORS = ROOT / "artifacts/generated-results/elkies-k3-curve398-11952-norm8-a1-exact-survivors-v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


class Curve398MW16AdaptiveHalfLatticeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.truth = json.loads(TRUTH.read_text())
        cls.input = json.loads(INPUT.read_text())
        cls.blind = json.loads(BLIND.read_text())
        cls.verified = json.loads(VERIFIED.read_text())
        cls.screen = json.loads(SCREEN.read_text())
        cls.survivors = json.loads(SURVIVORS.read_text())

    def test_complete_discovery_screen(self) -> None:
        self.assertEqual(self.screen["status"], "PASS_COMPLETE_MODULAR_SCREEN_TWO_SURVIVORS")
        self.assertEqual(self.screen["search"]["priority_table_class_count"], 63917)
        self.assertEqual(self.screen["search"]["excluded_count"], 63915)
        self.assertEqual(self.screen["search"]["survivor_priority_ranks"], [16875, 63669])
        self.assertEqual(self.survivors["status"], "PASS_EXACT_RATIONAL_PARAMETER_CANDIDATES")
        self.assertEqual(self.survivors["rational_candidate_count"], 2)
        self.assertTrue(
            all(
                record["specializations"][0]["isomorphic_to_curve398_over_Q"]
                for record in self.survivors["records"]
            )
        )

    def test_hidden_fibration_and_generic_subgroup(self) -> None:
        self.assertEqual(
            self.truth["status"],
            "PASS_EXACT_HIDDEN_A1_MW16_FIBRATION_PARAMETER_AND_PUBLIC_SUBGROUP",
        )
        self.assertEqual(self.truth["fibration"]["fibre_configuration"], "I2 at infinity + 22 I1")
        self.assertEqual(self.truth["fibration"]["priority_rank"], 16875)
        self.assertEqual(self.truth["generic_mw16"]["rank"], 16)
        self.assertEqual(self.truth["generic_mw16"]["height_gram_determinant"], "474")
        self.assertTrue(self.truth["generic_mw16"]["saturated"])
        self.assertTrue(self.truth["parameter_recovery"]["isomorphic_to_curve398_over_Q"])
        self.assertEqual(self.truth["parameter_recovery"]["factor_degrees_with_multiplicity"], [[1, 1], [23, 1]])

    def test_blind_input_is_redacted(self) -> None:
        self.assertEqual(self.input["status"], "PASS_REDACTED_GENERIC_MW16_INPUT")
        self.assertEqual(self.input["generic_rank"], 16)
        self.assertEqual(len(self.input["generic_points"]), 16)
        self.assertEqual(self.input["redaction"]["held_out_point_count"], 14)
        self.assertFalse(self.input["redaction"]["contains_public_embedding_coordinates"])
        self.assertNotIn("coordinates_in_public_rank30_points", INPUT.read_text())
        self.assertNotIn("heldout_complement", INPUT.read_text())

    def test_blind_transition_and_complete_searches(self) -> None:
        self.assertEqual(self.blind["status"], "STOPPED_AT_DECLARED_LIFT_LIMIT")
        self.assertEqual(self.blind["initial_search"]["basis_rank_after"], 21)
        self.assertEqual(len(self.blind["iterations"]), 1)
        self.assertEqual(self.blind["iterations"][0]["basis_rank_after"], 30)
        records = self.blind["initial_search"]["cover_records"] + self.blind["iterations"][0]["cover_records"]
        self.assertEqual(len(records), 384)
        self.assertTrue(all(row["search"]["status"] == "bounded_search_complete" for row in records))
        self.assertFalse(self.blind["blindness_boundary"]["public_rank30_fixture_loaded"])

    def test_post_search_exact_rediscovery(self) -> None:
        self.assertEqual(
            self.verified["status"],
            "PASS_EXACT_CROSS_FIBRATION_RANK30_REDISCOVERY",
        )
        self.assertEqual(self.verified["blind_transition"]["initial_independent_directions"], 5)
        self.assertEqual(self.verified["blind_transition"]["adaptive_independent_directions"], 9)
        self.assertEqual(self.verified["heldout_complement"]["count"], 14)
        self.assertTrue(self.verified["heldout_complement"]["all_rediscovered"])
        self.assertTrue(self.verified["subgroup_equality"]["blind_group_equals_public_rank30_subgroup"])

    def test_hash_chain(self) -> None:
        inputs = self.verified["inputs"]
        for path in (TRUTH, INPUT, BLIND):
            self.assertEqual(inputs[str(path.relative_to(ROOT))], digest(path))


if __name__ == "__main__":
    unittest.main()
