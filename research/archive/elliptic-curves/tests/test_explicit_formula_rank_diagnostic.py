#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import hashlib
import json
import sys
import unittest


CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))

from explicit_formula_rank_diagnostic import (  # noqa: E402
    DELTA,
    E20_MODEL,
    PRIME_LIMIT,
    T6793_MODEL,
    calibrated_bound,
    good_prime_power_sums,
    gp_program,
)


class ExplicitFormulaDiagnosticTests(unittest.TestCase):
    def test_power_sum_recurrence(self) -> None:
        self.assertEqual(good_prime_power_sums(-3, 11, 5), (-3, -13, 72, -73, -573))
        self.assertEqual(good_prime_power_sums(1, 5, 0), ())
        with self.assertRaises(ValueError):
            good_prime_power_sums(0, 1, 2)

    def test_gp_program_has_declared_support_and_models(self) -> None:
        program = gp_program((("E20", E20_MODEL), ("T6793", T6793_MODEL)))
        self.assertIn(f"D={DELTA};", program)
        self.assertIn(f"LIM={PRIME_LIMIT};", program)
        self.assertIn("p^k", program)
        self.assertIn("ROW|E20|", program)
        self.assertIn("ROW|T6793|", program)

    def test_pinned_relative_bound(self) -> None:
        records = {
            "E20": {
                "log_conductor": "170.08766484224888113543147761107035486432595000237",
                "prime_sum": "-53.527198246590105182372104818094548896473885897316",
            },
            "T6793": {
                "log_conductor": "158.57264848930317785280111545077626825294311412594",
                "prime_sum": "-51.437008629739551796598916661179959027140591528834",
            },
        }
        result = calibrated_bound(records)
        self.assertTrue(result["total_difference"].startswith("-1.2489998956987"))
        self.assertLess(float(result["calibrated_upper_value"]), 20.452)

    def test_generated_artifact_matches_script(self) -> None:
        root = Path(__file__).resolve().parents[2]
        script = root / "elliptic-curves/cas/explicit_formula_rank_diagnostic.py"
        artifact_path = (
            root
            / "archive/elliptic-curves/artifacts/generated-results/elliptic_nagao_rank21_t6793_explicit_formula.json"
        )
        artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
        self.assertEqual(
            artifact["script_sha256"], hashlib.sha256(script.read_bytes()).hexdigest()
        )
        self.assertTrue(artifact["comparison"]["calibrated_upper_value_less_than_21"])
        self.assertIn("no algebraic-rank upper bound", artifact["interpretation"]["unconditional"])


if __name__ == "__main__":
    unittest.main()
