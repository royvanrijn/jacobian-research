"""Run with sage -python -m unittest discover -s elliptic-curves/tests -p test_mw16_model_size.py."""
from pathlib import Path
from importlib.machinery import SourceFileLoader
import sys
import unittest

from sage.all import EllipticCurve, QQ

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "cas"))
import mw16_model_size as m


class ModelSizeTest(unittest.TestCase):
    def setUp(self):
        self.E = EllipticCurve(QQ, [-1, 1])
        self.P = self.E(0, 1)
        self.Q = 2*self.P

    def test_nonintegral_ordinate_scale(self):
        f = m.R([9, 0, 0, 0, 9])
        g, scale = m.normalize(f)
        self.assertEqual(scale, QQ(1)/3)
        self.assertEqual(g, m.z**4+1)

    def test_nonsquare_content_retains_only_the_nonsquare_factor(self):
        k = QQ(101)**100
        f = 2*k*k*(m.z**4+1)
        g, scale = m.normalize(f)
        self.assertEqual(g, 2*(m.z**4+1))
        self.assertEqual(scale, 1/k)

    def test_general_weierstrass_transport(self):
        F = self.E.change_weierstrass_model([QQ(3)/2, 7, 2, -3])
        phi, points = m.transport(self.E, F, [self.P, self.Q, 3*self.P])
        u, r, s, t = map(QQ, m.map_record(phi)["u_r_s_t"])
        for old, new in zip([self.P, self.Q, 3*self.P], points):
            self.assertEqual(old[0], u*u*new[0]+r)
            self.assertEqual(old[1], u**3*new[1]+s*u*u*new[0]+t)

    def test_full_projective_point_transport(self):
        f = m.pointed(self.E, self.Q)
        matrix = list(map(QQ, [2, 3, 1, 2]))
        g, scale = m.normalize(m.binary_transform(f, matrix))
        self.assertEqual(m.quartic_j(g), self.E.j_invariant())
        t = (self.P[1]+self.Q[1])/(self.P[0]-self.Q[0])
        v = 2*self.P[0]+self.Q[0]-t*t
        a, b, c, d = matrix
        z = (b-d*t)/(c*t-a)
        w = scale*(c*z+d)**2*v
        record = {"matrix_a_b_c_d": list(map(str, matrix)), "ordinate_scale": str(scale),
                  "integral_coefficients_ascending": [str(g[i]) for i in range(5)]}
        self.assertEqual(m.quartic_point_to_source(record, z, w, self.E, self.Q,
                         self.E.isomorphism_to(self.E)), self.P)
        with self.assertRaises(ArithmeticError):
            m.quartic_point_to_source(record, z, w+1, self.E, self.Q, self.E.isomorphism_to(self.E))

    def test_x_translation_cancels_from_pointed_quartic(self):
        r = QQ(17)/3
        F = EllipticCurve(QQ, [0, 3*r, 0, 3*r*r-1, r**3-r+1])
        q = F(self.Q[0]-r, self.Q[1])
        self.assertEqual(m.pointed(self.E, self.Q), m.pointed(F, q))

    def test_size_selection_is_exact_and_outcome_free(self):
        record = m.select_chart(self.E, self.Q, [self.P, 3*self.P, 4*self.P])
        self.assertEqual((record["maximum_bits"], record["total_bits"], record["name"]),
            min((t["maximum_bits"], t["total_bits"], t["name"]) for t in record["trials"]))
        self.assertLessEqual(m.integral_quartic_bit_lower_bound(self.E.j_invariant()), record["maximum_bits"])
        with self.assertRaises(ArithmeticError):
            m.binary_transform(m.pointed(self.E, self.Q), [1, 2, 2, 4])

    def test_direct_and_reduced_backends_return_exact_points(self):
        runner = SourceFileLoader("mw16_test_runner", str(
            Path(__file__).resolve().parents[1] / "cas/prepare_mw16_short_models.sage")).load_module()
        raw = self.E.change_weierstrass_model([QQ(1)/2, 0, 0, 0])
        phi = raw.isomorphism_to(self.E)
        candidate = {"raw_short_model": list(map(str, raw.ainvs()))}
        row = {"selected_short_model": list(map(str, self.E.ainvs()))}
        chart = {"base_point": m.point_record(self.Q),
                 "raw_base_point": m.point_record((~phi)(self.Q)),
                 "selected": m.select_chart(self.E, self.Q, [self.P, 3*self.P, 4*self.P])}
        for mode in ("raw_direct", "short_direct", "selected_direct", "selected_reduced"):
            result = runner.search_chart(candidate, row, chart, mode, 100, 2)
            self.assertEqual(result["status"], "bounded_search_complete")
            self.assertGreater(len(result["points_transported_to_raw"]), 0)
            for point in result["points_transported_to_raw"]:
                self.assertTrue(m.read_point(raw, point) in raw)


if __name__ == "__main__":
    unittest.main()
