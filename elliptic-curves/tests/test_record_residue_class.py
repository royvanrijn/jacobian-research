from __future__ import annotations

from fractions import Fraction
from math import gcd
from pathlib import Path
import sys
import unittest

CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from fermigier_mestre import NORMALIZED_RECORD_PARAMETER  # noqa: E402
from fermigier_mestre import DISCRIMINANT_FACTOR_COEFFICIENTS  # noqa: E402
from multiple_root_lifting import (  # noqa: E402
    affine_variable_coefficients,
    fixed_divisor_valuation,
)
from search_record_residue_class import (  # noqa: E402
    CRT_MODULUS,
    CRT_RESIDUE,
    RESIDUE_CONSTRAINTS,
    build_score_tables,
    exact_candidate_record,
    primitive_rationals_in_class,
    root_multiplicity_mod,
    score_rational,
)


class RecordResidueClassTests(unittest.TestCase):
    def test_crt_class_and_multiple_root_multiplicities(self) -> None:
        self.assertEqual(CRT_MODULUS, 7 * 17 * 37)
        for prime, residue in RESIDUE_CONSTRAINTS:
            self.assertEqual(CRT_RESIDUE % prime, residue)
        self.assertEqual(
            {
                prime: root_multiplicity_mod(prime, residue)
                for prime, residue in RESIDUE_CONSTRAINTS
            },
            {7: 16, 17: 4, 37: 3},
        )
        self.assertEqual(
            {
                prime: fixed_divisor_valuation(
                    affine_variable_coefficients(
                        DISCRIMINANT_FACTOR_COEFFICIENTS,
                        residue,
                        prime,
                    ),
                    prime,
                )
                for prime, residue in RESIDUE_CONSTRAINTS
            },
            {7: 18, 17: 4, 37: 3},
        )

    def test_enumerator_matches_a_direct_small_box(self) -> None:
        height = 200
        expected = []
        for denominator in range(1, height + 1):
            for numerator in range(-height, height + 1):
                if gcd(numerator, denominator) != 1:
                    continue
                if (numerator - CRT_RESIDUE * denominator) % CRT_MODULUS == 0:
                    expected.append((numerator, denominator))
        self.assertEqual(list(primitive_rationals_in_class(height)), expected)
        with self.assertRaises(ValueError):
            tuple(primitive_rationals_in_class(0))

    def test_score_excludes_denominator_and_bad_primes(self) -> None:
        tables = build_score_tables(37, "fermigier-good")
        denominator_prime = score_rational(1, 5, tables, include_traces=True)
        self.assertEqual(denominator_prime["skipped_denominator_primes"], 1)
        self.assertNotIn(
            5,
            [term["prime"] for term in denominator_prime["traces"]],
        )

        benchmark = score_rational(
            NORMALIZED_RECORD_PARAMETER.numerator,
            NORMALIZED_RECORD_PARAMETER.denominator,
            tables,
            include_traces=True,
        )
        used_primes = {term["prime"] for term in benchmark["traces"]}
        self.assertTrue({7, 17, 37}.isdisjoint(used_primes))
        self.assertGreaterEqual(benchmark["skipped_bad_primes"], 3)

    def test_benchmark_has_exact_constraint_valuations(self) -> None:
        self.assertEqual(
            (
                NORMALIZED_RECORD_PARAMETER.numerator
                - CRT_RESIDUE * NORMALIZED_RECORD_PARAMETER.denominator
            )
            % CRT_MODULUS,
            0,
        )
        tables = build_score_tables(37, "fermigier-good")
        record = exact_candidate_record(
            NORMALIZED_RECORD_PARAMETER.numerator,
            NORMALIZED_RECORD_PARAMETER.denominator,
            tables,
            {7: 18, 17: 4, 37: 3},
        )
        self.assertEqual(Fraction(record["t"]), NORMALIZED_RECORD_PARAMETER)
        self.assertTrue(record["is_published_benchmark"])
        self.assertEqual(record["h_valuations"], {"7": 18, "17": 4, "37": 3})


if __name__ == "__main__":
    unittest.main()
