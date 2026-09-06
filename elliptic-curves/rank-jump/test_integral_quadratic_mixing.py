import unittest
import integral_quadratic_mixing as m
import verify_integral_quadratic_mixing as v


class IntegralMixingTests(unittest.TestCase):
    def test_sign_and_regular_lattices_have_different_index(self):
        split, regular = m.module_counts(1, 1, 0), m.module_counts(1, 1, 1)
        self.assertEqual(split['integral_eigenspace_index'], 1)
        self.assertEqual(regular['integral_eigenspace_index'], 2)
        self.assertEqual(split['norm_index_plus'], 2)
        self.assertEqual(regular['norm_index_plus'], 1)

    def test_impossible_mixing_rejected(self):
        with self.assertRaises(AssertionError): m.module_counts(1, 3, 2)

    def test_standard_labels_matter(self):
        G0, G1 = m.span([4]), m.span([1, 2, 12])
        self.assertEqual(G0 & G1, {0})
        self.assertEqual(G0 & m.span([1, 2, 4]), {0, 4})

    def test_quadratic_point_is_not_an_eigenpoint(self):
        T = (v.cq(0, 1), v.cq(1, 2))
        Tc = (v.cq(0, -1), v.cq(1, -2))
        self.assertNotEqual(T, Tc)
        self.assertNotEqual(T, (Tc[0], v.neg(Tc[1])))
        self.assertEqual(v.point_add(T, (T[0], v.neg(T[1]))), None)

    def test_shared_Kummer_identity_by_integer_convolution(self):
        def product(a, b):
            out = [0]*(len(a)+len(b)-1)
            for i, x in enumerate(a):
                for j, y in enumerate(b): out[i+j] += x*y
            return out
        lhs = product([4, -1], [1, 1])+[0, 0]
        rhs = product([2, 0, 1], [2, 0, 1])
        self.assertEqual([a-b for a, b in zip(lhs, rhs)], product([0, -1], [-3, 5, 0, 1]))


if __name__ == '__main__': unittest.main()
