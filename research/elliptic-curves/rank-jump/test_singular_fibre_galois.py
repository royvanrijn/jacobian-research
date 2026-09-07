"""Regression controls for the prime-cycle and reduced-pole argument."""
from fractions import Fraction as F
import unittest
import singular_fibre_galois as s


class SingularFibreGaloisTests(unittest.TestCase):
    def test_actual_frobenius_power(self):
        self.assertEqual(s.power_cycle_type([1, 6, 17], 6), [1]*7+[17])
        rows, proved = s.primitive_gate(24, 17)
        self.assertTrue(proved)
        self.assertEqual([r['block_size'] for r in rows], [2, 3, 4, 6, 8, 12])

    def test_imprimitive_control_is_not_excluded(self):
        # D4 on four vertices is transitive and contains transpositions.
        # A 2-cycle does not force primitive action: opposite-pair blocks survive.
        rows, proved = s.primitive_gate(4, 2)
        self.assertFalse(proved)
        self.assertEqual(rows, [{'block_size': 2, 'number_of_blocks': 2,
                                'excluded_by_prime_cycle': False}])

    def test_decomposable_reduced_pole_control(self):
        # 1/(u^4-4u^2+1) = 1/(h(u)^2-3), h=u^2-2.
        h = list(map(F, [-2, 0, 1]))
        denominator = s.multiply(h, h)
        denominator[0] -= 3
        self.assertEqual(denominator, [1, 0, -4, 0, 1])
        derivative = [F(i)*denominator[i] for i in range(1, len(denominator))]
        self.assertEqual(s.polynomial_gcd(denominator, derivative), [1])
        # Squarefree poles alone do not exclude a proper composition.
        self.assertFalse(s.primitive_gate(4, 2)[1])

    def test_cancellation_is_detected(self):
        common = list(map(F, [2, 1]))
        a = s.multiply(common, list(map(F, [3, 0, 1])))
        b = s.multiply(common, list(map(F, [1, 1])))
        self.assertEqual(s.polynomial_gcd(a, b), common)


if __name__ == '__main__':
    unittest.main()
