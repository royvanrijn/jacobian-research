#!/usr/bin/env sage
"""Exact mapping regressions; run with sage -python -m unittest."""
from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"elliptic-curves/cas"))
from sage.all import EllipticCurve, QQ
import run_fixed_field_point_realization as m
from fixed_cubic_field_curve_family import fixed_field_cubic_coefficients, field_multiply


class FixedFieldPointRealizationTests(unittest.TestCase):
    def test_pair_policy_keeps_all_classes(self):
        basis = [{"mask":1},{"mask":6},{"mask":10}]
        self.assertEqual([v[0] for v in m.class_masks(basis,2)],[1,6,10,7,11,12])
        with self.assertRaises(ArithmeticError):
            m.class_masks([{"mask":1},{"mask":1}],2)

    def test_change_with_scaled_ordinate_and_shift(self):
        # Old Y^2=x^4+1, Y=2y+1 gives y^2+y=x^4/4.
        old = [["1","0","0","0","1"],["0"]*3]
        new = [["0","0","0","0","1/4"],["1","0","0"]]
        change = {"e":"2","matrix":[["1","0"],["0","1"]],"H":["1","0","0"]}
        m.verify_change(old,new,change)
        bad = deepcopy(change)
        bad["H"][0] = "2"
        with self.assertRaises(AssertionError):
            m.verify_change(old,new,bad)

    def test_inverse_map_retains_parameter_infinity(self):
        change = {"e":"1","matrix":[["0","-1"],["1","0"]],"H":["0"]*3}
        self.assertEqual(m.inverse_change([0,1,1],change),(-1,0,1))
        self.assertEqual(m.inverse_change([1,0,1],change),(0,1,1))

    def test_curve_and_actual_kummer_identity(self):
        A,B = QQ(-7),QQ(31)
        E = EllipticCurve(QQ,[A,B])
        # (5,11) lies on y^2=x^3-7x+31.
        result = m.verify_point((5,-1,0),11,(1,0,0),1,A,B,0,E)
        self.assertEqual(result["raw_point"],["5","11"])
        with self.assertRaises(AssertionError):
            m.verify_point((5,-1,0),11,(2,0,0),1,A,B,0,E)

    def test_translation_recovers_original_class_exactly(self):
        A,B = QQ(-7),QQ(11)
        c0,c1,c2,_ = fixed_field_cubic_coefficients(A,B,-1)
        E = EllipticCurve(QQ,[0,QQ(c2),0,QQ(c1),QQ(c0)])
        Q = E(A+1,A-B+1)
        T,P = 2*Q,3*Q
        eta = (A+1,QQ(-1),QQ(1))
        slope = (T[1]-Q[1])/(T[0]-Q[0])
        line = (Q[1]-slope*Q[0],slope,-slope)
        # delta(2Q)=(theta-1/2)^2 and the chord norm identity.
        g = m.divide(line,field_multiply(eta,(-QQ(1)/2,1,0),A,B),A,B)
        if (A-B+1)*m.norm(g,A,B) != P[1]:
            g = tuple(-v for v in g)
        point = m.verify_point(eta,A-B+1,g,1,A,B,-1,E)
        self.assertEqual(point["raw_point"],list(map(str,P[:2])))
        recovered = m.translate_back(point,(1,0,0),QQ(1),eta,A,B,E)
        self.assertEqual(recovered["raw_point"],list(map(str,T[:2])))

    def test_universal_point_has_separating_valuation(self):
        data,run,A,B,E = m.context(m.SOURCE,-1)
        cert = m.universal_point_certificate(data,A,B,E)
        self.assertEqual(cert["prime_ideal_residue_degrees"],[1,2])
        self.assertEqual(cert["point_kummer_valuations"],[0,1])
        self.assertTrue(all(all(v%2==0 for v in row) for row in cert["anchor_kummer_valuations"]))

    def test_pinned_positive_controls_and_corrupt_point(self):
        path = ROOT/"artifacts/generated-results/elliptic-curves/fixed_field_point_realization_positive_controls_v1.json"
        controls = json.loads(path.read_text())
        self.assertEqual(controls["source_sha256"],m.digest(m.SOURCE))
        data,run,A,B,E = m.context(m.SOURCE,0)
        self.assertEqual([r["mask"] for r in controls["cases"]],[1,6])
        for row in controls["cases"]:
            self.assertGreater(len(row["points"]),0)
            m.replay_row(row,data,A,B,0,E)
        self.assertGreater(controls["cases"][0]["rational_parameter_infinity_point_count"],0)
        bad = deepcopy(controls["cases"][0])
        bad["points"][0]["raw_point"][0] = "0"
        with self.assertRaises(AssertionError):
            m.replay_row(bad,data,A,B,0,E)


if __name__ == "__main__":
    unittest.main()
