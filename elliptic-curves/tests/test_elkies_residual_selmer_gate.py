from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "cas/elkies_residual_selmer_gate.py"
SPEC = importlib.util.spec_from_file_location("elkies_residual_selmer_gate", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
GATE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = GATE
SPEC.loader.exec_module(GATE)


class ElkiesResidualSelmerGateTests(unittest.TestCase):
    def test_pari_ellrank_fields_recover_total_selmer_dimension(self) -> None:
        self.assertEqual(
            GATE.pari_ellrank_total_two_selmer_dimension(
                rank_lower=28,
                rank_upper=30,
                cassels_pairing_rank=2,
                two_torsion_dimension=0,
            ),
            32,
        )

    def test_pari_ellrank_cassels_rank_is_not_arbitrary_sha_dimension(self) -> None:
        with self.assertRaisesRegex(GATE.ResidualSelmerGateError, "odd Cassels"):
            GATE.pari_ellrank_total_two_selmer_dimension(
                rank_lower=28,
                rank_upper=31,
                cassels_pairing_rank=1,
                two_torsion_dimension=0,
            )

    def test_pari_ellrank_reversed_interval_fails_closed(self) -> None:
        with self.assertRaisesRegex(GATE.ResidualSelmerGateError, "reversed"):
            GATE.pari_ellrank_total_two_selmer_dimension(
                rank_lower=29,
                rank_upper=28,
                cassels_pairing_rank=0,
                two_torsion_dimension=0,
            )

    def test_dimension_below_fifteen_is_exact_rejection(self) -> None:
        record = GATE.gate_record(total_two_selmer_dimension=31)
        self.assertEqual(record["residual_two_selmer_quotient_dimension"], 14)
        self.assertEqual(record["status"], GATE.REJECT_STATUS)
        self.assertFalse(record["expensive_search_authorized"])

    def test_dimension_fifteen_passes(self) -> None:
        record = GATE.gate_record(total_two_selmer_dimension=32)
        self.assertEqual(record["residual_two_selmer_quotient_dimension"], 15)
        self.assertEqual(record["status"], GATE.PASS_STATUS)

    def test_signature_or_incomplete_backend_cannot_authorize_search(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "gate.json"
            path.write_text(
                json.dumps(
                    {
                        "schema": GATE.SCHEMA,
                        "status": GATE.PASS_STATUS,
                        "parameter": "-9529/5471",
                        "global_minimal_model": ["1", "-1", "1", "-2", "3"],
                        "descent_backend": {
                            "unconditional": True,
                            "class_group_completeness_completed": False,
                            "all_local_solubility_conditions_completed": True,
                        },
                        "gate": {
                            "residual_two_selmer_quotient_dimension": 15,
                            "required_residual_dimension": 15,
                            "expensive_search_authorized": True,
                        },
                    }
                )
            )
            with self.assertRaisesRegex(GATE.ResidualSelmerGateError, "global/local"):
                GATE.require_expensive_search_gate(path)

    def test_gate_is_bound_to_exact_minimal_curve(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "gate.json"
            path.write_text(
                json.dumps(
                    {
                        "schema": GATE.SCHEMA,
                        "status": GATE.PASS_STATUS,
                        "parameter": "-9529/5471",
                        "global_minimal_model": ["1", "-1", "1", "-2", "3"],
                        "descent_backend": {
                            "unconditional": True,
                            "class_group_completeness_completed": True,
                            "all_local_solubility_conditions_completed": True,
                        },
                        "gate": {
                            "residual_two_selmer_quotient_dimension": 15,
                            "required_residual_dimension": 15,
                            "expensive_search_authorized": True,
                        },
                    }
                )
            )
            document = GATE.require_expensive_search_gate(
                path, expected_model=(1, -1, 1, -2, 3)
            )
            self.assertEqual(document["status"], GATE.PASS_STATUS)
            with self.assertRaisesRegex(GATE.ResidualSelmerGateError, "different minimal"):
                GATE.require_expensive_search_gate(
                    path, expected_model=(1, -1, 1, -2, 4)
                )

    def test_specialization_adapter_binds_parameter_and_model(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "gate.json"
            path.write_text(
                json.dumps(
                    {
                        "schema": GATE.SCHEMA,
                        "status": GATE.PASS_STATUS,
                        "parameter": "-9529/5471",
                        "global_minimal_model": ["1", "-1", "1", "-2", "3"],
                        "descent_backend": {
                            "unconditional": True,
                            "class_group_completeness_completed": True,
                            "all_local_solubility_conditions_completed": True,
                        },
                        "gate": {
                            "residual_two_selmer_quotient_dimension": 15,
                            "required_residual_dimension": 15,
                            "expensive_search_authorized": True,
                        },
                    }
                )
            )
            specialization = {
                "parameter": {"affine_value": "-9529/5471"},
                "global_minimal_specialization": {
                    "model": ["1", "-1", "1", "-2", "3"]
                },
            }
            document = GATE.require_gate_for_specialization(path, specialization)
            self.assertEqual(document["status"], GATE.PASS_STATUS)
            specialization["parameter"]["affine_value"] = "-2/377"
            with self.assertRaisesRegex(GATE.ResidualSelmerGateError, "different fibre"):
                GATE.require_gate_for_specialization(path, specialization)


if __name__ == "__main__":
    unittest.main()
