"""Small exact regressions for the population experiment, run with Sage."""
from importlib.machinery import SourceFileLoader
from pathlib import Path
import unittest

try:
    from sage.all import EllipticCurve, QQ, ZZ
except ImportError:
    EllipticCurve = None

ROOT = Path(__file__).resolve().parents[2]


@unittest.skipIf(EllipticCurve is None, 'requires sage -python')
class FibreHeightPopulationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = SourceFileLoader('height_population_tests',str(
            ROOT/'elliptic-curves/cas/fibre_height_population.sage')).load_module()
        cls.protocol = cls.module.read(cls.module.PROTOCOL)
        cls.families = cls.module.families(cls.protocol)

    def test_new_coordinate_changes_population_and_preserves_exact_fibre_identity(self):
        m = self.module
        for family in self.families:
            scale, shift = family['coordinate']
            self.assertNotEqual(scale,0)
            for parameter in (QQ(-1),QQ(1),QQ(1)/2):
                mapped = shift+scale*parameter
                self.assertEqual((mapped-shift)/scale,parameter)
                old, new = m.exact_fibre(family,parameter),m.exact_fibre(family,mapped)
                self.assertNotEqual(old.j_invariant(),new.j_invariant())
                self.assertLess(m.logheight(new.j_invariant()),0.6*m.logheight(old.j_invariant()))
                hints=[v for q in (scale,shift,mapped) for v in (q.numerator(),q.denominator())]
                normalized=m.remove_obvious_scale(new,hints)
                self.assertTrue(new.is_isomorphic(normalized))
                self.assertEqual(new.j_invariant(),normalized.j_invariant())

    def test_scale_removal_handles_square_mixed_with_nonsquare_content(self):
        m=self.module
        # Deliberately mix an enormous coordinate scale with a nonsquare
        # twist factor, the case that defeated whole-content root tests.
        large=ZZ(2)**127-1
        u=large*QQ(7)/11
        E=EllipticCurve(QQ,[-25*u**4,125*u**6])
        S=m.remove_obvious_scale(E,[large,7,11])
        self.assertTrue(S.is_isomorphic(E))
        self.assertLess(max(abs(a.numerator()).nbits() for a in S.ainvs()),20)

    def test_same_j_does_not_collapse_distinct_rational_twists(self):
        E=EllipticCurve(QQ,[-1,1])
        twist=E.quadratic_twist(2)
        self.assertEqual(E.j_invariant(),twist.j_invariant())
        self.assertFalse(E.is_isomorphic(twist))
        self.assertTrue(E.is_isomorphic(EllipticCurve(QQ,[-16,64])))

    def test_prospective_score_ignores_discovery_fields(self):
        m=self.module
        rows=[]
        for coordinate in ('identity','improved'):
            for i in range(1,13):
                rows.append(dict(id=f'{coordinate}-{i}',coordinate=coordinate,box_pair=[i,1],
                    j_log2_height=i,normalized_weierstrass_log2_size=20-i,sample_cost_seconds=i/100,
                    nagao=dict(block_score_units_1e12=[i,i,i],block_good_bad_counts=[[4,0]]*3)))
        protocol={**self.protocol,'arm_size':2}
        before=m.select_arms(rows,protocol)
        for row in rows:
            row['certified_new_directions']=1000
            row['finite_curve_points']=[{'x':'0','y':'1'}]
        self.assertEqual(m.select_arms(rows,protocol),before)
        self.assertTrue(set(before['nagao']) <= {r['id'] for r in rows})

    def test_independence_certificates_survive_json_round_trip(self):
        import json
        result=self.module.independence(['0','0','0','-2','1'],[{'x':'0','y':'1'}],100)
        self.assertEqual(json.loads(json.dumps(result)),result)

    def test_logarithmic_height_does_not_overflow_at_prospective_sizes(self):
        import math
        self.assertEqual(self.module.logheight(QQ(2)**20000),20000)
        self.assertTrue(math.isfinite(self.module.logheight(QQ(2)**20000/3)))

    def test_nonfinite_measurements_fail_closed_before_freezing(self):
        import tempfile
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'bad.json'
            with self.assertRaises(ValueError):
                self.module.write(output,{'height':float('inf')})
            self.assertFalse(output.exists())

    def test_frozen_families_have_no_cross_family_curve_overlap(self):
        if not self.module.POPULATION.exists():
            self.skipTest('population preparation has not completed')
        p=self.module.read(self.module.POPULATION)
        left={r['j_invariant'] for r in p['families']['mw16']['rows']}
        right={r['j_invariant'] for r in p['families']['mw18']['rows']}
        self.assertFalse(left & right)

    def test_retained_mw16_centres_match_fixed_generic_gram(self):
        from sage.all import matrix,vector
        template=self.module.read(self.module.TEMPLATE)
        p=next(p for p in template['presentations'] if p['presentation_id']==self.protocol['mw16_presentation'])
        gram=matrix(QQ,p['generic_height_gram'])
        representatives=self.protocol['centres']['mw16']
        self.assertEqual(len(representatives),12)
        self.assertEqual({vector(QQ,r)*gram*vector(QQ,r) for r in representatives},{QQ(23)/2})
        masks={sum((int(c)%2)<<i for i,c in enumerate(r)) for r in representatives}
        self.assertEqual(len(masks),12)


if __name__=='__main__':
    unittest.main()
