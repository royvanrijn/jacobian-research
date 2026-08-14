from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from certify_nagao_rank17_frontier import (  # noqa: E402
    CANDIDATE_PARAMETERS,
    exact_log_conductor_certificate,
    input_subset,
    load_candidates,
)


Q = Fraction


class NagaoRank17FrontierCertificateTests(unittest.TestCase):
    def setUp(self) -> None:
        generated = ROOT / "artifacts/generated-results"
        self.input_paths = (
            generated / "elliptic_nagao_rank13_rank_gain_search.json",
            generated / "elliptic_nagao_rank13_rank_gain_mutations.json",
        )

    def test_candidate_loader_pins_four_frontier_parameters(self) -> None:
        records = load_candidates(self.input_paths)
        self.assertEqual(tuple(records), CANDIDATE_PARAMETERS)
        self.assertEqual(tuple(input_subset(records[item]) and item for item in records), CANDIDATE_PARAMETERS)

    def test_each_input_has_17_exact_coordinate_records(self) -> None:
        records = load_candidates(self.input_paths)
        for parameter, record in records.items():
            points = input_subset(record)
            self.assertEqual(len(points), 17)
            self.assertEqual(Q(record["parameter_u"]), parameter)
            self.assertTrue(record["conductor_probe"]["below_strict_log_conductor_target"])

    def test_exact_log_bound_uses_only_rational_inequalities(self) -> None:
        certificate = exact_log_conductor_certificate(10**66 - 1)
        self.assertEqual(certificate["decimal_digit_count"], 66)
        self.assertEqual(certificate["deduced_log_conductor_upper_bound"], "7623/50")
        self.assertTrue(certificate["strict_target_proved_exactly"])
        with self.assertRaises(AssertionError):
            exact_log_conductor_certificate(10**80)

    def test_generated_artifact_contains_exact_certificates_when_present(self) -> None:
        path = ROOT / "artifacts/generated-results/elliptic_nagao_rank17_frontier_certificate.json"
        if not path.exists():
            self.skipTest("generated certificate has not been reproduced")
        data = json.loads(path.read_text())
        script = CAS / "certify_nagao_rank17_frontier.py"
        engine = CAS / "mod2_reduction_independence.py"
        self.assertEqual(data["script_sha256"], hashlib.sha256(script.read_bytes()).hexdigest())
        self.assertEqual(
            data["certificate_engine_sha256"],
            hashlib.sha256(engine.read_bytes()).hexdigest(),
        )
        for input_record, input_path in zip(data["inputs"], self.input_paths):
            self.assertEqual(
                input_record["sha256"], hashlib.sha256(input_path.read_bytes()).hexdigest()
            )
        self.assertEqual(len(data["certificates"]), 4)
        for record in data["certificates"]:
            certificate = record["finite_reduction_certificate"]
            self.assertEqual(certificate["combined_exact_rank_over_F2"], 17)
            self.assertEqual(certificate["certified_algebraic_rank_lower_bound"], 17)
            self.assertEqual(len(record["saturated_basis"]), 17)
            self.assertTrue(
                record["direct_conductor_replay"]["conductor_matches_discovery_artifact"]
            )
            self.assertTrue(record["exact_log_conductor_bound"]["strict_target_proved_exactly"])


if __name__ == "__main__":
    unittest.main()
