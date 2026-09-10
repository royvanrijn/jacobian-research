import importlib.util
from pathlib import Path
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "rank-jump" / "prospective_constructed_strict_seed_gate.py"
SPEC = importlib.util.spec_from_file_location("constructed_strict_seed_gate", SCRIPT)
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)


def signal(**overrides):
    value = {
        "schema": gate.SCHEMA,
        "provenance": "PROSPECTIVE",
        "historical_point_inputs": 0,
        "candidate_id": "curve302:prospective-principal-dependency-001",
        "generic_strict_dimension": 0,
        "principal_dependency": {"status": "CERTIFIED", "certificate": "principal.json"},
        "constructed_strict_class": {"status": "CERTIFIED_NONZERO", "certificate": "strict.json"},
        "cover_solubility": {"status": "UNKNOWN"},
        "new_rational_point": {"status": "UNKNOWN"},
        "independent_quotient_direction": {"status": "UNKNOWN"},
    }
    value.update(overrides)
    return value


class ConstructedStrictSeedGateTests(unittest.TestCase):
    def test_curve302_first_nonzero_strict_class_is_a_cover_seed(self):
        result = gate.evaluate(signal())
        self.assertTrue(result["strict_nonzero_is_outside_generic_strict_subgroup"])
        self.assertEqual(result["next_action"], "CONSTRUCT_AND_AUDIT_COVER")

    def test_every_gate_is_required_before_foundry_admission(self):
        row = signal(
            cover_solubility={"status": "CERTIFIED_SOLUBLE", "certificate": "cover.json"},
            new_rational_point={"status": "CERTIFIED_NEW", "certificate": "point.json"},
            independent_quotient_direction={"status": "CERTIFIED_INDEPENDENT", "certificate": "quotient.json"},
        )
        self.assertEqual(gate.evaluate(row)["next_action"], "ADMIT_CERTIFIED_QUOTIENT_DIRECTION_TO_FOUNDRY")

    def test_no_independence_without_a_rational_point(self):
        row = signal(independent_quotient_direction={"status": "CERTIFIED_INDEPENDENT", "certificate": "quotient.json"})
        with self.assertRaisesRegex(ValueError, "requires a certified new rational point"):
            gate.evaluate(row)

    def test_oracle_provenance_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "forbidden"):
            gate.evaluate(signal(provenance="RETROSPECTIVE"))
