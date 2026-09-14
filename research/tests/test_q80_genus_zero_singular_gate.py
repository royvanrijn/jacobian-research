"""Regression gates for discriminant boundaries that could invalidate exclusion."""
import importlib.util
from pathlib import Path
import unittest

SCRIPT = Path(__file__).resolve().parents[1]/'elkies-k3/scripts/verify_q80_genus_zero_singular_gate.py'
SPEC = importlib.util.spec_from_file_location('singular_gate',SCRIPT)
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)


class DiscriminantBoundaryTests(unittest.TestCase):
    def test_simple_and_repeated_old_base_infinity(self):
        # t^3+t has one simple binary root at infinity; t^2+1 has two.
        self.assertNotEqual(V.quartic_discriminant([0,1,0,1,0],521),0)
        self.assertEqual(V.quartic_discriminant([1,0,1,0,0],521),0)

    def test_even_roots_are_retained(self):
        self.assertEqual(V.root_multiplicity([1,519,1],1,521),2)
        self.assertEqual(V.root_multiplicity([520,3,518,1],1,521),3)

    def test_zero_polynomial_is_not_exclusion(self):
        with self.assertRaises(ValueError):
            V.check_discriminant({'discriminant_coefficients':[0]*25},{},521)

    def test_parameter_infinity_cannot_be_omitted(self):
        # A constant affine quartic family has Q(t;u,v)=v^4*(t^4+1).
        frame = {'branch_matrix':[[1,0,0,0,0],[0]*5,[0]*5,[0]*5,[1,0,0,0,0]]}
        row = {'discriminant_coefficients':[256]+[0]*24,
               'finite_roots_with_multiplicity':[], 'infinity_multiplicity':24}
        self.assertFalse(V.check_discriminant(row,frame,521))
        row['infinity_multiplicity'] = 0
        with self.assertRaises(ValueError):
            V.check_discriminant(row,frame,521)


if __name__ == '__main__':
    unittest.main()
