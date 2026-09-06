import unittest
import quadratic_norm_blocks as n


class QuadraticNormBlockTests(unittest.TestCase):
    def test_full_Galois_action_removes_rotations(self):
        cycle, swap = (0, 1, 1, 1), (0, 1, 1, 0)
        self.assertEqual(len(n.centralizer([cycle])), 3)
        self.assertEqual(n.centralizer([cycle, swap]), [[1, 0, 0, 1]])

    def test_split_two_torsion_has_all_six_identifications(self):
        self.assertEqual(len(n.centralizer([])), 6)

    def test_reducible_cubic_does_not_pass_S3_gate(self):
        with self.assertRaises(AssertionError):
            n.cubic_certificate([0, -1, 0, 1], 31)

    def test_cyclic_control_is_not_mislabelled_S3(self):
        cert = n.cubic_certificate([-1, -14, -11, 1], 31)
        self.assertEqual(cert['group'], 'C3')
        self.assertEqual(int(cert['discriminant']), 163**2)


if __name__ == '__main__':
    unittest.main()
