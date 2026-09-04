"""Small adversarial arithmetic checks added by the external audit."""

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
from types import SimpleNamespace

PROGRAM = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(PROGRAM), str(PROGRAM / 'cas')]

from alternate_quartic_covers import short_add
from ecsearch.rank_certification import add_rational_points
from verify_record_prime_factors import prove_primes
from elkies_residual_selmer_gate import monotone_sieve_gate_record, ResidualSelmerGateError


class AuditRegressions(unittest.TestCase):
    def test_legacy_search_gate_rejects_nonfinite_limits(self):
        for value in (float('inf'), float('nan'), 0, True):
            with self.assertRaises(ResidualSelmerGateError):
                monotone_sieve_gate_record(
                    stages=[{'residual_upper_bound': None,
                             'proof_status': 'NO_FINITE_UPPER_BOUND_YET'}],
                    search_limits={'wall_seconds': value})

    def test_integer_point_addition_stays_rational(self):
        model = (0, 0, 0, -36, 0)
        expected = (Q(-144, 25), Q(-504, 125))
        for add in (short_add, add_rational_points):
            self.assertEqual(add(model, (-3, 9), (12, 36)), expected)

    def test_prime_proof_rejects_composite_and_partial_backend_output(self):
        for output in ('[2, 15, 4]\n15|0\n', '[2, 15, 4]\n'):
            with patch('verify_record_prime_factors.subprocess.run', return_value=
                       SimpleNamespace(stdout=output, stderr='')):
                with self.assertRaises(ArithmeticError):
                    prove_primes([15])

    def test_prime_proof_checks_each_factor_and_rejects_unbounded_timeout(self):
        with patch('verify_record_prime_factors.subprocess.run', return_value=
                   SimpleNamespace(stdout='[2, 15, 4]\n2|1\n3|1\n', stderr='')):
            self.assertEqual(prove_primes([3, 2])[1], ['2', '3'])
        for limit in (float('inf'), float('nan'), 0):
            with self.assertRaises(ValueError):
                prove_primes([2], timeout=limit)


if __name__ == '__main__':
    unittest.main()
