from copy import deepcopy
from fractions import Fraction as Q
from math import gcd, isqrt
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"cas"))
from search_observability import point_visibility, masked_control, multiply, transform


def fixture(m=(1, 0, 0, 1), centre=(1, 2), den=1, k=0, H=5):
    """An independent small direct enumerator, not the production search code."""
    m = tuple(map(Q, m))
    xq, yq = map(Q, centre)
    raw = (-3*xq*xq+28, -8*yq, -6*xq, Q(0), Q(1))
    raw_matrix = multiply((Q(den), Q(k, den), Q(0), Q(1)), m)
    final = tuple(f/(den*den) for f in transform(raw, raw_matrix))
    assert all(f.denominator == 1 for f in final)
    record = {"input": {"curve": ["0", "0", "0", "-7", "10"]},
        "short_model": ["0", "0", "0", "-7", "10"], "short_model_x_shift": "0",
        "base_point": dict(zip(("x", "y"), map(str, centre))),
        "pointed_chart": {"curve_coordinate_scale": "1", "point_denominator_root": str(den),
            "shift_mod_denominator_squared": str(k), "unimodular_horizontal_matrix": ["1", "0", "0", "1"]},
        "horizontal_matrix": list(map(str, m)), "ordinate_scale": "1", "coefficients": list(map(str, final)),
        "height_bound": H, "denominator_start": 1, "denominator_end": H, "completed_denominator": H,
        "infinity_checked": True, "primitive_square_hits": [], "finite_curve_points": []}
    seen = set()
    coordinates = [(n, d) for d in range(1, H+1) for n in range(-H, H+1) if gcd(n, d) == 1]+[(1, 0)]
    for n, d in coordinates:
        f = sum(c*n**i*d**(4-i) for i, c in enumerate(final))
        if f < 0 or isqrt(f.numerator)**2 != f:
            continue
        r = isqrt(f.numerator)
        record["primitive_square_hits"].append([str(n), str(d), str(r)])
        a, b, c, e = raw_matrix
        lower = c*n+e*d
        if not lower:
            continue
        slope = (a*n+b*d)/lower
        for signed in {r, -r}:
            w = Q(den)*signed/lower**2
            x = (slope*slope-xq+w)/2
            y = slope*(x-xq)-yq
            assert y*y == x**3-7*x+10
            if (x, y) != (xq, yq):
                seen.add((x, y))
    record["finite_curve_points"] = [dict(zip(("x", "y"), map(str, p))) for p in sorted(seen)]
    return record


class ObservabilityTests(unittest.TestCase):
    def test_raw_known_point_and_complement(self):
        rec = fixture()
        for p in ((2, 2), (13, 46)):
            result = point_visibility(rec, p)
            self.assertEqual(result["coordinate"], ["4", "1"])
            self.assertEqual(result["minimum_affine_height"], 4)
            self.assertEqual(result["status"], "VISIBLE_AND_RECORDED")

    def test_tangent_removable_singularity(self):
        result = point_visibility(fixture(), (1, -2))
        self.assertEqual(result["coordinate"], ["1", "1"])
        self.assertEqual(result["status"], "VISIBLE_AND_RECORDED")

    def test_rational_horizontal_matrix_and_ordinate_scaling(self):
        rec = fixture()
        rec["horizontal_matrix"] = ["1/2", "0", "0", "1"]
        rec["ordinate_scale"] = "4"
        rec["coefficients"] = ["400", "-128", "-24", "0", "1"]
        del rec["height_bound"]
        self.assertEqual(point_visibility(rec, (2, 2))["coordinate"], ["8", "1"])

    def test_curve_scaling_is_undone_in_raw_slope(self):
        rec = fixture()
        rec["input"]["curve"] = rec["short_model"] = ["0", "0", "0", "-7/16", "5/32"]
        rec["base_point"] = {"x": "1/4", "y": "1/4"}
        rec["pointed_chart"]["curve_coordinate_scale"] = "2"
        rec["finite_curve_points"] = [
            {"x": str(Q(p["x"])/4), "y": str(Q(p["y"])/8)} for p in rec["finite_curve_points"]]
        self.assertEqual(point_visibility(rec, (Q(1, 2), Q(1, 4)))["status"], "VISIBLE_AND_RECORDED")

    def test_known_endpoints_are_not_failures(self):
        for p in (None, (1, 2)):
            self.assertEqual(point_visibility(fixture(), p)["status"], "KNOWN_POINTED_ENDPOINT")

    def test_coordinate_infinity_is_checked_outside_affine_box(self):
        rec = fixture((4, 1, 1, 0), H=1)
        result = point_visibility(rec, (2, 2))
        self.assertEqual(result["coordinate"], ["1", "0"])
        self.assertEqual(result["status"], "VISIBLE_AND_RECORDED")
        self.assertIsNone(result["minimum_affine_height"])

    def test_height_is_coordinate_dependent(self):
        self.assertEqual(point_visibility(fixture(H=3), (2, 2))["status"], "OUTSIDE_BOX")
        self.assertEqual(point_visibility(fixture((4, 0, 0, 1), H=1), (2, 2))["status"], "VISIBLE_AND_RECORDED")

    def test_inverse_transport_with_denominators(self):
        rec = fixture((1, 1, 0, 1), centre=(Q(9, 4), Q(19, 8)), den=2, k=1)
        self.assertTrue(rec["finite_curve_points"])
        for p in rec["finite_curve_points"]:
            self.assertEqual(point_visibility(rec, p)["status"], "VISIBLE_AND_RECORDED")

    def test_incomplete_coverage_is_not_a_miss(self):
        rec = fixture()
        rec["completed_denominator"] = 0
        rec["primitive_square_hits"] = []
        rec["finite_curve_points"] = []
        self.assertEqual(point_visibility(rec, (2, 2))["status"], "UNSEARCHED_INTERVAL")

    def test_visible_missed_point_is_detected_without_resieving(self):
        rec = fixture()
        rec["finite_curve_points"] = []
        result = point_visibility(rec, (2, 2))
        self.assertEqual(result["status"], "VISIBLE_NOT_RECORDED")
        self.assertTrue(result["recorded_square_hit"])

    def test_chart_only_no_transcript(self):
        rec = fixture()
        del rec["height_bound"]
        self.assertEqual(point_visibility(rec, (2, 2))["status"], "OBSERVABLE_WITHOUT_TRANSCRIPT")

    def test_tampering(self):
        for key, value in (("ordinate_scale", "2"), ("horizontal_matrix", [1, 0, 0, 2]),
                           ("short_model_x_shift", "1"), ("coefficients", [1, 2, 3, 4, 5])):
            rec = fixture()
            rec[key] = value
            with self.assertRaises(ArithmeticError):
                point_visibility(rec, (2, 2))
        with self.assertRaises(ValueError):
            point_visibility(fixture(), (2, 3))
        with self.assertRaises(ValueError):
            point_visibility(fixture(), (2.0, 2))

    def test_all_retained_hits_under_several_coordinates(self):
        for m in ((1, 0, 0, 1), (1, 1, 0, 1), (0, 1, 1, 0), (2, -1, 1, 1), (-1, 0, 0, 1)):
            rec = fixture(m, H=8)
            for p in rec["finite_curve_points"]:
                self.assertEqual(point_visibility(rec, p)["status"], "VISIBLE_AND_RECORDED")

    def test_general_weierstrass_transport(self):
        rec = fixture()
        rec["input"]["curve"] = ["2", "2", "4", "-8", "0"]
        rec["short_model_x_shift"] = "1"
        rec["base_point"] = {"x": "0", "y": "0"}
        rec["finite_curve_points"] = [
            {"x": str(Q(p["x"])-1), "y": str(Q(p["y"])-Q(p["x"])-1)}
            for p in rec["finite_curve_points"]]
        self.assertEqual(point_visibility(rec, (1, -1))["status"], "VISIBLE_AND_RECORDED")

    def test_two_torsion_centre(self):
        rec = fixture()
        rec["input"]["curve"] = rec["short_model"] = ["0", "0", "0", "-1", "0"]
        rec["base_point"] = {"x": "0", "y": "0"}
        rec["coefficients"] = ["4", "0", "0", "0", "1"]
        del rec["height_bound"]
        self.assertEqual(point_visibility(rec, (1, 0))["coordinate"], ["0", "1"])

    def test_masking_separates_input_from_oracle(self):
        pts = [(1, 2), (2, 2), (3, 4)]
        gram = [[4, 1, 2], [1, 6, 1], [2, 1, 8]]
        search, oracle = masked_control([-7, 10], pts, gram, [1])
        self.assertEqual(search["points"], [{"x": "1", "y": "2"}, {"x": "3", "y": "4"}])
        self.assertEqual(search["metric_gram"], [["4", "2"], ["2", "8"]])
        self.assertEqual(oracle["withheld_points"], [{"x": "2", "y": "2"}])
        self.assertEqual(set(search), {"curve", "points", "metric_gram"})
        self.assertIn("NOT_NEW_RANK", oracle["endpoint"])
        for bad in ([], [0, 1, 2], [3], [1, 1], [True]):
            with self.assertRaises(ValueError):
                masked_control([-7, 10], pts, gram, bad)


if __name__ == "__main__":
    unittest.main()
