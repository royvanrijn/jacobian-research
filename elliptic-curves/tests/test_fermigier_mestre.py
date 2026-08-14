from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import shutil
import sys
import unittest

CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from fermigier_mestre import (  # noqa: E402
    FermigierMestreFamily,
    NORMALIZED_RECORD_PARAMETER,
    ROOTS,
)
from pari_bridge import minimal_curve_data  # noqa: E402


class FermigierMestreTests(unittest.TestCase):
    FAMILY = FermigierMestreFamily

    def test_mestre_polynomial_identity(self) -> None:
        for t, x in ((Q(1), Q(2)), (Q(7, 3), Q(-4, 5)), (Q(5), Q(100))):
            product = Q(1)
            for root in ROOTS:
                product *= (x - t - root) * (x + t - root)
            g = self.FAMILY.square_approximant(t, x)
            quartic = self.FAMILY.quartic_value(t, x)
            self.assertEqual(product, g**2 - (50616 * t) ** 2 * quartic)

    def test_binary_invariant_discriminant_identity(self) -> None:
        for t in (Q(1), Q(2), Q(7, 3), NORMALIZED_RECORD_PARAMETER):
            invariant_i, invariant_j = self.FAMILY.binary_invariants(t)
            self.assertEqual(
                4 * invariant_i**3 - invariant_j**2,
                432 * self.FAMILY.discriminant_factor(t),
            )

    def test_thirteen_visible_points_satisfy_quartic(self) -> None:
        for t in (Q(1), Q(7, 3)):
            points = self.FAMILY.known_quartic_points(t)
            self.assertEqual(len(points), 13)
            for x, y in points:
                self.assertEqual(y**2, self.FAMILY.quartic_value(t, x))
        with self.assertRaises(ValueError):
            self.FAMILY.visible_quartic_points(Q(0))

    def test_hensel_roots_and_valuations(self) -> None:
        expected = {89: {1111, 6810}, 131: {510, 16651}}
        for prime, residues in expected.items():
            roots = self.FAMILY.power_roots(prime, 2, split_only=True)
            self.assertEqual({root.residue for root in roots}, residues)
            for root in roots:
                self.assertEqual(root.modulus, prime**2)
                self.assertTrue(root.split_multiplicative)
                self.assertEqual(
                    self.FAMILY.verify_power_constraint(
                        root.residue, 1, prime, 2
                    ),
                    2,
                )

    def test_bad_and_good_local_data_are_distinct(self) -> None:
        bad = self.FAMILY.local_data(43, 89)
        self.assertFalse(bad.good_reduction)
        self.assertEqual((bad.point_count, bad.trace), (89, 1))
        self.assertTrue(bad.split_multiplicative)

        good = self.FAMILY.local_data(0, 89)
        self.assertTrue(good.good_reduction)
        self.assertEqual((good.point_count, good.trace), (105, -15))
        self.assertIsNone(good.split_multiplicative)

    @unittest.skipUnless(shutil.which("gp"), "PARI/GP is optional")
    def test_published_record_minimal_model(self) -> None:
        data = minimal_curve_data(
            self.FAMILY.coefficients(NORMALIZED_RECORD_PARAMETER),
            timeout=10,
        )
        self.assertEqual(
            data["minimal_model"],
            (
                1,
                0,
                1,
                -940299517776391362903023121165864,
                10707363070719743033425295515449274534651125011362,
            ),
        )
        self.assertEqual(
            data["conductor"],
            22720638514787473197194583889675055980109503436060704437972911338086049759883790,
        )
        self.assertTrue(data["log_conductor"].startswith("182.724910950637428796"))


if __name__ == "__main__":
    unittest.main()
