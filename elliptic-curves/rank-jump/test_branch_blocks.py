import unittest
from fractions import Fraction as F
import branch_blocks as b


class BranchGeometryTests(unittest.TestCase):
    def test_disjoint_quadratics(self):
        g = b.character_data([[1,0,1], [2,0,1], [3,0,1]])
        self.assertEqual((g["geometric_character_rank"], g["branch_points"], g["genus"]), (3,6,5))

    def test_shared_branch_gives_genus_one_triple(self):
        g = b.character_data([[0,-1,1], [0,-2,1], [0,-3,1]])
        self.assertEqual((g["geometric_character_rank"], g["branch_points"], g["genus"]), (3,4,1))

    def test_dependent_triangle_and_infinity(self):
        g = b.character_data([[0,-1,1], [0,-2,1], [2,-3,1]])
        self.assertEqual((g["geometric_character_rank"], g["branch_points"], g["genus"]), (2,3,0))
        h = b.character_data([[0,1], [-1,1]])
        self.assertTrue(h["infinity_branched"])
        self.assertEqual((h["branch_points"], h["genus"]), (3,0))

    def test_geometric_constants_ignored_but_not_rational_solubility(self):
        # t^2+1 and 2(t^2+1) cannot both be nonzero rational squares at
        # rational t. The production analysis separately verifies roots.
        g = b.character_data([[1,0,1], [2,0,2]])
        self.assertEqual(g["geometric_character_rank"], 1)

    def test_specialization_kernel(self):
        coords = [[F(1),F(0)], [F(0),F(1)], [F(2),F(2)], [F(0),F(0)]]
        result = b.specialization_blocks(coords)
        self.assertEqual(result["kernel_dimension_modulo_generic"], 2)
        self.assertEqual(result["kernel_basis_integer_coefficients"], [[-2,-2,1,0], [0,0,0,1]])

    def test_reject_repeated_branch(self):
        with self.assertRaises(AssertionError):
            b.character_data([[1,2,1]])


if __name__ == "__main__":
    unittest.main()
