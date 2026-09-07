import unittest
from itertools import permutations
from torsion_difference import (add, compose, inject, section,
                                canonical_subset, baer_coordinates)


class TorsionDifferenceTests(unittest.TestCase):
    def test_subset_section_is_linear(self):
        for a in range(4):
            for b in range(4):
                self.assertEqual(canonical_subset(section(a) ^ section(b)), section(a ^ b))
                self.assertEqual(canonical_subset(inject(a) ^ inject(b)), inject(a ^ b))

    def test_baer_coordinates_ignore_diagonal_two_torsion(self):
        for a in range(4):
            for b in range(4):
                x = (a, b)
                y = ((a+2) % 4, b)
                for t in range(4):
                    xt = tuple((x[i]+2*((t >> i) & 1)) % 4 for i in range(2))
                    yt = tuple((y[i]+2*((t >> i) & 1)) % 4 for i in range(2))
                    self.assertEqual(baer_coordinates(x, y), baer_coordinates(xt, yt))

    def test_discriminant_scalar_cocycle_is_a_coboundary(self):
        H = (2, 3)
        for perm in permutations(range(3)):
            g = (perm[0]+1, perm[1]+1)
            sign = sum(perm[i] > perm[j] for i in range(3) for j in range(i+1, 3)) % 2
            self.assertEqual(add(compose(g, H), compose(H, g)), g if sign else (0, 0))


if __name__ == "__main__":
    unittest.main()
