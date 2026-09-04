"""Small exact regressions; run with sage -python -m unittest discover."""
from pathlib import Path
import runpy
import unittest

from sage.all import EllipticCurve, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
M = runpy.run_path(str(ROOT / "elliptic-curves/cas/certify_exceptional_soluble_selmer_panel.sage"))


class SolublePanelTests(unittest.TestCase):
    def setUp(self):
        self.model = list(map(QQ, [0, 0, 0, -1, 1]))
        self.curve = EllipticCurve(QQ, self.model)
        self.point = self.curve(0, 1)
        self.ring = PolynomialRing(QQ, names=("u", "v", "w", "z"))

    def test_cover_witness_and_cubic_norm_identity(self):
        f, points = M["cubic_data"](self.model, [[self.point[0], self.point[1]]])
        cover = M["build_cover"](f, points[0], "P", self.ring)
        self.assertTrue(cover["witness_verified"])
        for equation in cover["primitive_quadrics"]:
            self.assertEqual(self.ring(equation)(1, 0, 0, 1), 0)
        # Independent numeric check of the full covering norm identity, away
        # from the rational witness and without imposing either quadric.
        from sage.all import Matrix
        alpha = list(map(QQ, cover["alpha_coefficients"]))
        mul = M["multiply_mod_cubic"]
        basis = ([1, 0, 0], [0, 1, 0], [0, 0, 1])
        norm = lambda v: Matrix(QQ, [mul(v, e, f.list()) for e in basis]).det()
        for beta in ([2, 3, 5], [-1, 2, 0], [QQ(1)/3, 0, 2]):
            product = mul(alpha, mul(beta, beta, f.list()), f.list())
            self.assertEqual(norm(product), points[0][1]**2 * norm(beta)**2)

    def test_bad_point_and_bad_cover_witness_rejected(self):
        with self.assertRaises((TypeError, ValueError, ArithmeticError)):
            M["cubic_data"](self.model, [[QQ(0), QQ(2)]])
        f, points = M["cubic_data"](self.model, [[self.point[0], self.point[1]]])
        with self.assertRaises(ArithmeticError):
            M["verify_rational_cover_witness"]([points[0][0], -1, 0], f.list(), [1, 0, 0, 2], self.ring)

    def test_doubling_and_negation_do_not_create_kummer_directions(self):
        p, q = self.point, 2*self.point
        f, points = M["cubic_data"](self.model, [[p[0], p[1]], [q[0], q[1]], [p[0], -p[1]]])
        rows, blocks = M["signatures"](f, points, 1, 2, 101)
        self.assertTrue(blocks)
        self.assertEqual(M["f2_rank"](rows), 1)
        self.assertTrue(all(bit == 0 for bit in rows[1]))
        self.assertEqual(rows[0], rows[2])

    def test_rational_two_torsion_hypothesis_rejected(self):
        with self.assertRaises(ArithmeticError):
            M["cubic_data"](list(map(QQ, [0, 0, 0, -1, 0])), [])

    def test_universal_pointed_quartic_map_and_discriminant(self):
        ring = PolynomialRing(QQ, names=("a", "b", "A", "t", "w"))
        a, b, A, t, w = ring.gens()
        B = b*b-a**3-A*a
        f = t**4-6*a*t*t-8*b*t-3*a*a-4*A
        x = (t*t-a+w)/2
        y = t*(t*t-3*a+w)/2-b
        self.assertEqual(ring.ideal(w*w-f).reduce(y*y-x**3-A*x-B), 0)
        poly_ring = PolynomialRing(ring, "T")
        T = poly_ring.gen()
        quartic = T**4-6*a*T*T-8*b*T-3*a*a-4*A
        self.assertEqual(quartic.discriminant(), 256*(-16)*(4*A**3+27*B**2))


if __name__ == "__main__":
    unittest.main()
