from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from fermigier_mestre import FermigierMestreFamily  # noqa: E402


class FermigierCovariantTests(unittest.TestCase):
    FAMILY = FermigierMestreFamily

    def test_covariants_have_pinned_exact_values(self) -> None:
        self.assertEqual(
            self.FAMILY.quartic_covariants_at(Q(1), Q(2)),
            (
                Q(2139124873018494148237496979, 4),
                Q(-12454576363441822574258795056406969749005, 2),
            ),
        )

    def test_pinned_quartic_point_maps_to_jacobian(self) -> None:
        quartic_point = self.FAMILY.visible_quartic_points(Q(1))[0]
        self.assertEqual(quartic_point, (Q(-1), Q(135510147)))
        self.assertEqual(
            self.FAMILY.quartic_point_to_jacobian(Q(1), quartic_point),
            (
                Q(684548705792558025841219, 627988096849),
                Q(
                    -139306783692409064523024544135831690,
                    497653563264667993,
                ),
            ),
        )

    def test_all_known_images_satisfy_the_exact_jacobian_equation(self) -> None:
        for parameter in (Q(1), Q(7, 3)):
            points = self.FAMILY.known_jacobian_points(parameter)
            self.assertEqual(len(points), 13)
            self.assertEqual(len(set(points)), 13)
            _, _, _, coefficient_a, coefficient_b = self.FAMILY.coefficients(
                parameter
            )
            for x, y in points:
                self.assertEqual(y**2, x**3 + coefficient_a * x + coefficient_b)

    def test_changing_quartic_sign_negates_jacobian_image(self) -> None:
        parameter = Q(1)
        x, z = self.FAMILY.visible_quartic_points(parameter)[0]
        positive = self.FAMILY.quartic_point_to_jacobian(parameter, (x, z))
        negative = self.FAMILY.quartic_point_to_jacobian(parameter, (x, -z))
        self.assertEqual(positive[0], negative[0])
        self.assertEqual(positive[1], -negative[1])

    def test_affine_map_rejects_invalid_source_points(self) -> None:
        with self.assertRaisesRegex(ValueError, "requires z != 0"):
            self.FAMILY.quartic_point_to_jacobian(Q(1), (Q(0), Q(0)))
        with self.assertRaisesRegex(ValueError, "not on the quartic"):
            self.FAMILY.quartic_point_to_jacobian(Q(1), (Q(0), Q(1)))


if __name__ == "__main__":
    unittest.main()
