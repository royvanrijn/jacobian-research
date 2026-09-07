from __future__ import annotations

import unittest

from audit_mestre_diameter235_displayed_lattice import displayed_basis
from search_mestre_root_tuple_scale import point_on_short_curve


class Diameter235DisplayedLatticeTest(unittest.TestCase):
    def test_exact_seed_basis_membership(self) -> None:
        coefficients, points = displayed_basis()
        self.assertEqual(len(points), 11)
        self.assertTrue(all(point_on_short_curve(coefficients, point) for point in points))


if __name__ == "__main__":
    unittest.main()
