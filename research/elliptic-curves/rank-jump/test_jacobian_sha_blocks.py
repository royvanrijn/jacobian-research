import unittest
from fractions import Fraction as F
from jacobian_sha_blocks import geometry, nonisogeny_witness, symplectic_basis
import local_collision as lc
import retrospective as r


class JacobianShaBlockTests(unittest.TestCase):
    def test_all_alternating_four_dimensional_forms(self):
        for mask in range(64):
            B = [[0]*4 for _ in range(4)]
            k = 0
            for i in range(4):
                for j in range(i+1, 4):
                    B[i][j] = B[j][i] = (mask >> k) & 1
                    k += 1
            pairs, radical = symplectic_basis(B)
            chosen = [v for pair in pairs for v in pair]
            self.assertEqual(r.rank(chosen + radical), 4)
            for v in chosen:
                self.assertTrue(any(lc.pairing(v, w, B) for w in chosen))
            for v in radical:
                self.assertFalse(any(lc.pairing(v, 1 << j, B) for j in range(4)))

    def test_smooth_synthetic_geometry_and_degenerations(self):
        self.assertTrue(geometry(F(-2), F(2), 2)["both_quotient_equations_verified"])
        with self.assertRaises(AssertionError):
            geometry(F(-2), F(2), 0)
        with self.assertRaises(AssertionError):
            geometry(F(-2), F(1), 1)  # D=0

    def test_equal_counts_are_not_an_isogeny_certificate(self):
        result = nonisogeny_witness(F(-2), F(2), 0, [11, 19, 23])
        self.assertEqual(result["status"], "UNKNOWN")
        self.assertEqual(result["extension_splitting_status"], "UNKNOWN")


if __name__ == "__main__":
    unittest.main()
