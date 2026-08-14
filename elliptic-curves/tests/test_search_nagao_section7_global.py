from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from search_nagao_section7_global import (  # noqa: E402
    expected_primitive_pair_count,
    exact_predeclared_seeds,
    parse_scanner_output,
    run_global_scanner,
)


SCANNER = ROOT / "elliptic-curves/cas/scan_nagao_section7_global.cpp"
SCRIPT = ROOT / "elliptic-curves/cas/search_nagao_section7_global.py"
ARTIFACT = ROOT / "artifacts/generated-results/elliptic_nagao_section7_global.json"


def synthetic_table(prime: int, multiplier: int):
    return tuple(
        SimpleNamespace(contribution=multiplier * residue / 100.0)
        for residue in range(prime + 1)
    )


class NagaoSection7GlobalTests(unittest.TestCase):
    def test_scanner_exhausts_primitive_box_and_does_not_select_on_validation(self):
        training = {5: synthetic_table(5, 1), 7: synthetic_table(7, -1)}
        validation_a = {11: synthetic_table(11, 1)}
        validation_b = {11: synthetic_table(11, -100)}
        first, first_audit = run_global_scanner(
            source=SCANNER,
            a_max=31,
            b_max=9,
            global_keep=12,
            per_denominator_keep=2,
            training=training,
            validation=validation_a,
            compile_timeout=20,
            scan_timeout=20,
        )
        second, second_audit = run_global_scanner(
            source=SCANNER,
            a_max=31,
            b_max=9,
            global_keep=12,
            per_denominator_keep=2,
            training=training,
            validation=validation_b,
            compile_timeout=20,
            scan_timeout=20,
        )
        self.assertEqual(
            first_audit["primitive_pairs_enumerated"],
            expected_primitive_pair_count(31, 9),
        )
        self.assertEqual(
            second_audit["primitive_pairs_enumerated"],
            first_audit["primitive_pairs_enumerated"],
        )
        self.assertEqual(
            {candidate.parameter for candidate in first},
            {candidate.parameter for candidate in second},
        )
        self.assertTrue(
            any(
                left.validation_scaled != right.validation_scaled
                for left, right in zip(
                    sorted(first, key=lambda item: item.parameter),
                    sorted(second, key=lambda item: item.parameter),
                )
            )
        )

    def test_scanner_parser_rejects_nonprimitive_rows(self):
        output = "SUMMARY\t1\t1\t1\t1\nROW\t2\t2\t1\t1\n"
        with self.assertRaisesRegex(AssertionError, "nonprimitive"):
            parse_scanner_output(output, a_max=2, b_max=2)

    def test_all_twenty_one_generic_sections_are_exactly_predeclared(self):
        quartic, jacobian, coefficients = exact_predeclared_seeds(Fraction(19, 7))
        self.assertEqual(len(quartic), 21)
        self.assertEqual(len(jacobian), 21)
        self.assertEqual(len(coefficients), 5)

    def test_generated_artifact_is_pinned_when_present(self):
        if not ARTIFACT.exists():
            self.skipTest("the global section-7 artifact is absent")
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            data["reproducibility"]["script_sha256"],
            hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            data["reproducibility"]["scanner_sha256"],
            hashlib.sha256(SCANNER.read_bytes()).hexdigest(),
        )
        self.assertEqual(data["population"]["primitive_pairs_enumerated"], 18_244_819)
        self.assertEqual(data["leakage_free_scoring"]["training_prime_band"], [5, 199])
        self.assertEqual(
            data["leakage_free_scoring"]["heldout_validation_prime_band"],
            [211, 397],
        )
        self.assertEqual(data["leakage_free_scoring"]["exact_final_prime_band"], [401, 1999])
        self.assertTrue(data["point_triage"]["all_21_generic_sections_decontaminated"])
        self.assertEqual(data["target"]["certified_hits"], [])
        self.assertEqual(data["final_frontier"][0]["constructor_parameter_T"], "599/2")
        self.assertEqual(data["final_frontier"][0]["stable_numerical_rank"], 17)
        checkpoint = data["exact_checkpoints_stable_numerical_rank_at_least_17"]
        self.assertEqual(set(checkpoint), {"section7-global-426-1", "section7-global-599-2"})
        self.assertTrue(
            all(record["status"] == "certified" for record in checkpoint.values())
        )
        self.assertTrue(
            all(
                record["certified_algebraic_rank_lower_bound"] == 17
                and record["combined_exact_rank_over_F2"] == 17
                for record in checkpoint.values()
            )
        )
        self.assertEqual(
            data["final_frontier"][0]["conductor"]["conductor"],
            756855624206125019617445466192515119979543579877522270,
        )


if __name__ == "__main__":
    unittest.main()
