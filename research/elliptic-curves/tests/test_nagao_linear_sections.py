from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))

from nagao_1994 import (  # noqa: E402
    RANK13_CONSTRUCTION,
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    rank13_base_parameter,
    rank13_extra_point,
    short_jacobian_coefficients,
)
from nagao_linear_sections import (  # noqa: E402
    COMPANION_JACOBIAN_RELATIONS,
    LINEAR_COMPANION_SECTIONS,
    companion_abscissae,
    omitted_companion_sections,
    point_on_short_curve,
)


Q = Fraction


def add_short(
    coefficients: tuple[Fraction, ...],
    left: tuple[Fraction, Fraction] | None,
    right: tuple[Fraction, Fraction] | None,
) -> tuple[Fraction, Fraction] | None:
    if left is None:
        return right
    if right is None:
        return left
    _, _, _, coefficient_a, _ = coefficients
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if y1 == -y2:
            return None
        slope = (3 * x1**2 + coefficient_a) / (2 * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope**2 - x1 - x2
    return x3, -y1 + slope * (x1 - x3)


class NagaoLinearSectionTests(unittest.TestCase):
    def test_all_six_formulas_hold_exactly(self) -> None:
        for parameter in (Q(1), Q(42), Q(3631, 14), Q(-17, 3)):
            points = tuple(section.point(parameter) for section in LINEAR_COMPANION_SECTIONS)
            self.assertEqual(len({point[0] for point in points}), 6)
            coefficients = short_jacobian_coefficients(RANK13_CONSTRUCTION, parameter)
            for section in LINEAR_COMPANION_SECTIONS:
                self.assertTrue(point_on_short_curve(coefficients, section.jacobian_point(parameter)))

    def test_first_formula_is_nagao_printed_extra_section(self) -> None:
        for parameter in (Q(2), Q(3631, 14)):
            self.assertEqual(LINEAR_COMPANION_SECTIONS[0].point(parameter), rank13_extra_point(parameter))

    def test_u42_companions_explain_five_old_nonvisible_abscissae(self) -> None:
        parameter = rank13_base_parameter(Q(42))
        expected_omitted = {
            Q(6211, 210),
            Q(1829, 10),
            Q(-355, 6),
            Q(47189, 70),
            Q(-39983, 210),
        }
        self.assertEqual(
            {section.point(parameter)[0] for section in omitted_companion_sections()},
            expected_omitted,
        )
        self.assertEqual(len(companion_abscissae(parameter)), 6)

    def test_five_pinned_dependencies_hold_exactly(self) -> None:
        # A second specialization, distinct from the symbolic verifier's scope,
        # keeps the unit test fast while checking the public formulas end to end.
        parameter = rank13_base_parameter(Q(42))
        coefficients = short_jacobian_coefficients(RANK13_CONSTRUCTION, parameter)
        visible = tuple(
            quartic_point_to_short_jacobian(RANK13_CONSTRUCTION, parameter, point)
            for point in primitive_visible_points(RANK13_CONSTRUCTION, parameter)
        )
        companions = {
            section.label: section.jacobian_point(parameter)
            for section in LINEAR_COMPANION_SECTIONS
        }
        basis = visible[:11] + (companions["plus-1/15"],)
        for label, relation in COMPANION_JACOBIAN_RELATIONS.items():
            total = None
            for point, coefficient in zip(basis, relation):
                if coefficient:
                    summand = point if coefficient == 1 else (point[0], -point[1])
                    total = add_short(coefficients, total, summand)
            self.assertEqual(total, companions[label])


if __name__ == "__main__":
    unittest.main()
