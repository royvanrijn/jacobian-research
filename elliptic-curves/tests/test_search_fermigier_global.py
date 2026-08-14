from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from search_fermigier_global import (  # noqa: E402
    ModularTable,
    build_modular_tables,
    expected_primitive_pair_count,
    parse_scanner_output,
    polynomial_mod,
    run_scanner,
)


SCANNER = CAS / "scan_fermigier_global.cpp"
SCRIPT = CAS / "search_fermigier_global.py"
ARTIFACT = ROOT / "artifacts/generated-results/elliptic_fermigier_global.json"


def synthetic_table(prime: int, rank_multiplier: int, power_residue: int) -> ModularTable:
    modulus = prime * prime
    rank = tuple(rank_multiplier * residue for residue in range(prime)) + (0,)
    power = tuple(
        rank_multiplier if residue == power_residue else 0
        for residue in range(modulus)
    ) + (0,)
    return ModularTable(prime, rank, power)


class FermigierGlobalTests(unittest.TestCase):
    def test_scanner_exhausts_box_and_never_selects_on_held_band(self):
        discovery = (synthetic_table(5, 7, 6), synthetic_table(7, -3, 8))
        held_a = (synthetic_table(11, 1, 12),)
        held_b = (synthetic_table(11, -1000, 13),)
        first, first_audit = run_scanner(
            source=SCANNER,
            a_max=31,
            b_max=9,
            frontier_keep=12,
            per_denominator_keep=2,
            discovery=discovery,
            held=held_a,
            compile_timeout=20,
            scan_timeout=20,
        )
        second, second_audit = run_scanner(
            source=SCANNER,
            a_max=31,
            b_max=9,
            frontier_keep=12,
            per_denominator_keep=2,
            discovery=discovery,
            held=held_b,
            compile_timeout=20,
            scan_timeout=20,
        )
        self.assertEqual(
            first_audit["primitive_pairs_enumerated"],
            expected_primitive_pair_count(31, 9),
        )
        self.assertEqual(second_audit["primitive_pairs_enumerated"], 185)
        self.assertEqual(
            {candidate.parameter for candidate in first},
            {candidate.parameter for candidate in second},
        )
        by_parameter_a = {candidate.parameter: candidate for candidate in first}
        by_parameter_b = {candidate.parameter: candidate for candidate in second}
        self.assertTrue(
            any(
                by_parameter_a[parameter].held_rank_scaled
                != by_parameter_b[parameter].held_rank_scaled
                for parameter in by_parameter_a
            )
        )

    def test_parser_rejects_nonprimitive_row(self):
        output = "SUMMARY\t1\t1\t1\t1\t1\t1\nROW\t2\t2\t1\t2\t3\t4\n"
        with self.assertRaisesRegex(AssertionError, "not primitive"):
            parse_scanner_output(output, a_max=2, b_max=2)

    def test_power_table_marks_exact_p_squared_congruences(self):
        table = build_modular_tables(7, 7)[0]
        marked = [
            residue
            for residue, weight in enumerate(table.power_weights[:-1])
            if weight > 0
        ]
        self.assertEqual(marked, [0, 7, 14, 21, 28, 35, 42])
        self.assertTrue(all(polynomial_mod(residue, 49) == 0 for residue in marked))
        self.assertTrue(
            all(
                table.power_weights[residue] == 0
                for residue in range(49)
                if polynomial_mod(residue, 49) != 0
            )
        )

    def test_generated_artifact_is_pinned(self):
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            data["reproducibility"]["script_sha256"],
            hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            data["reproducibility"]["scanner_sha256"],
            hashlib.sha256(SCANNER.read_bytes()).hexdigest(),
        )
        self.assertEqual(data["population"]["a_max"], 100_000)
        self.assertEqual(data["population"]["b_max"], 1_000)
        self.assertEqual(data["population"]["primitive_pairs_enumerated"], 60_815_684)
        self.assertEqual(
            data["population"]["primitive_pairs_enumerated"],
            data["population"]["expected_primitive_pair_count_reference"],
        )
        self.assertEqual(data["population"]["retained_union_count"], 16_133)
        self.assertEqual(data["population"]["novel_nonsingular_discovery_union"], 16_110)
        self.assertEqual(data["leakage_free_features"]["discovery_prime_band"], [401, 499])
        self.assertEqual(data["leakage_free_features"]["held_prime_band"], [503, 599])
        self.assertTrue(data["leakage_free_features"]["bands_disjoint"])
        self.assertFalse(data["leakage_free_features"]["held_information_used_by_scanner_retention"])
        self.assertEqual(data["conductor_summary"]["completed"], 38)
        self.assertEqual(data["conductor_summary"]["below_strict_target"], 22)
        self.assertEqual(data["frontier_summary"]["best_stable_numerical_rank"], 13)
        self.assertEqual(
            [stage["input_count"] for stage in data["point_triage"]["stages"]],
            [12, 4, 2],
        )
        self.assertEqual(data["target"]["certified_hits"], [])
        self.assertEqual(data["point_triage"]["exact_certificates"], {})
        self.assertTrue(data["point_triage"]["all_calls_no_retry"])
        self.assertNotIn(
            Fraction(39508, 39),
            {
                Fraction(row["parameter_t"])
                for row in data["selection"]["held_validation_candidates"]
            },
        )


if __name__ == "__main__":
    unittest.main()
