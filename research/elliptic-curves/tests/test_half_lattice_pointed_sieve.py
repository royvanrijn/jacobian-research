from fractions import Fraction as Q
from math import gcd, isqrt
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "cas"))
import half_lattice_pointed_sieve as sieve
from alternate_quartic_covers import alternate_cover, short_add


class PointedSieveTests(unittest.TestCase):
    model = (Q(0), Q(0), Q(0), Q(-1), Q(1))

    def brute(self, f, h, first=1, last=None):
        out = set()
        for d in range(first, (h if last is None else last)+1):
            for n in range(-h,h+1):
                if gcd(n,d) != 1:
                    continue
                value = sum(c*n**i*d**(4-i) for i,c in enumerate(f))
                if value >= 0 and isqrt(value)**2 == value:
                    out.add((n,d,isqrt(value)))
        return out

    def test_modular_sieve_matches_exhaustive_including_bad_primes_and_zero(self):
        # Quartics with zero, negative values, prime-divisible content and
        # square leading coefficient exercise all projective residue cases.
        for f in ((4,-8,0,0,1), (0,1,-2,3,0), (-7,0,3,0,1),
                  (9,0,18,0,9), (105,-210,0,0,105)):
            record, points = sieve.search_box(f,71,2)
            self.assertEqual(record["status"],"bounded_search_complete")
            self.assertEqual(set(points),self.brute(f,71))
            self.assertEqual(record["integer_pairs_covered"],71*143)

    def test_denominator_shards_are_disjoint_and_complete(self):
        f = (4,-8,0,0,1)
        _, one = sieve.search_box(f,35,2,1,17)
        _, two = sieve.search_box(f,35,2,18,35)
        self.assertFalse(set(one)&set(two))
        self.assertEqual(set(one)|set(two),self.brute(f,35))

    def test_timeout_never_claims_complete_box(self):
        record, points = sieve.search_box((4,-8,0,0,1),10000,1e-9)
        self.assertEqual(record["status"],"bounded_search_timeout")
        self.assertEqual(record["completed_denominator"],0)
        self.assertEqual(points,())

    def test_exact_denominator_transform_and_maps_under_huge_scaling(self):
        point = (Q(0),Q(1))
        for _ in range(5):
            point = short_add(self.model,point,(Q(0),Q(1)))
            for scale in (Q(1),Q(7,13),Q(2**1200,3**400)):
                model = (0,0,0,self.model[3]/scale**4,self.model[4]/scale**6)
                base = (point[0]/scale**2,point[1]/scale**3)
                chart = sieve.make_chart(model,base)
                integral_scale, integral = sieve.integral_short_scale(tuple(map(Q,model)))
                self.assertEqual(sieve.invariants(chart.coefficients),(-48*integral[3],-1728*integral[4]))
                a,b,c,d = chart.matrix
                for s in (Q(-4,7),Q(0),Q(2,3)):
                    if c*s+d == 0:
                        continue
                    t = (chart.denominator*(a*s+b)/(c*s+d)+Q(chart.shift,chart.denominator))/integral_scale
                    expected = alternate_cover(model,base).value(t)
                    expected *= (c*s+d)**4*integral_scale**4/chart.denominator**2
                    actual = sum(value*s**i for i,value in enumerate(chart.coefficients))
                    self.assertEqual(actual,expected)
                for n,den,r in self.brute(chart.coefficients,12):
                    for sign in {r,-r}:
                        answer = chart.map_point(n,den,sign)
                        if answer is not None:
                            self.assertTrue(sieve.point_on_short_curve(model,answer))

    def test_points_at_transformed_infinity_and_raw_pole(self):
        outcome = sieve.run_quartic_search(mask=1,representative=[1],
            short_model=self.model,generic_points=[(Q(0),Q(1))],
            height_bound=12,timeout_seconds=2,stack_bytes=0)
        self.assertTrue(outcome.record["infinity_checked"])
        self.assertIn((Q(1),Q(1)),outcome.curve_points)
        self.assertIn((Q(3),Q(5)),outcome.curve_points)
        # A deliberately exchanged horizontal basis makes transformed
        # infinity map to raw slope zero, a finite, nontrivial point.
        chart = sieve.make_chart(self.model,(Q(0),Q(1)))
        exchanged = sieve.PointedChart(chart.model,chart.base_point,chart.curve_scale,
            chart.denominator,chart.shift,(0,1,1,0),(1,0,0,-8,4))
        self.assertEqual(exchanged.map_point(1,0,2),(Q(1),Q(-1)))
        self.assertIsNone(exchanged.map_point(0,1,1))

    def test_input_rejections(self):
        with self.assertRaises(ValueError):
            sieve.make_chart(self.model,(Q(2),Q(1)))
        with self.assertRaises(ValueError):
            sieve.integral_short_scale((Q(0),)*5)
        with self.assertRaises(ValueError):
            sieve.search_box((1,0,0,0,1),0,2)

    def test_gmp_multiscalar_matches_independent_rational_group_law(self):
        p = (Q(0),Q(1))
        points = [p,short_add(self.model,p,p),(Q(3),Q(5))]
        for coefficients in ([0,0,0],[1,-2,3],[-7,11,-18],[1,0,0]):
            expected = sieve.linear_combination_python(self.model,points,coefficients)
            self.assertEqual(sieve.linear_combination(self.model,points,coefficients),expected)
        # The same exact relation on a model with thousands of coefficient bits.
        scale = Q(2**900,3**301)
        model = (0,0,0,self.model[3]*scale**4,self.model[4]*scale**6)
        large = [(x*scale**2,y*scale**3) for x,y in points]
        self.assertEqual(sieve.linear_combination(model,large,[1,-2,3]),
                         sieve.linear_combination_python(model,large,[1,-2,3]))
        # Reduced denominators may conceal a fourth power in A while B
        # retains the common scaling valuation.
        u, integral = sieve.integral_short_scale((Q(0),Q(0),Q(0),Q(-1,101**3),Q(1,101**6)))
        self.assertEqual(u,Q(101))
        self.assertTrue(all(v.denominator == 1 for v in integral))

    def test_checkpoint_resume_budgets_and_corruption(self):
        kwargs = dict(mask=1,representative=[1],short_model=self.model,
                      generic_points=[(Q(0),Q(1))],height_bound=12,
                      timeout_seconds=2,stack_bytes=0)
        with tempfile.TemporaryDirectory() as directory:
            backend = sieve.CheckpointedBackend(directory)
            original = backend.run_quartic_search(**kwargs)
            with patch.object(sieve,"run_quartic_search",side_effect=RuntimeError("fresh search")):
                self.assertEqual(backend.run_quartic_search(**kwargs),original)
                with self.assertRaises(RuntimeError):
                    backend.run_quartic_search(**{**kwargs,"height_bound":13})
            checkpoint = next(Path(directory).glob("*.json"))
            text = checkpoint.read_text().replace('"infinity_checked": true','"infinity_checked": false')
            checkpoint.write_text(text)
            with self.assertRaises(ArithmeticError):
                backend.run_quartic_search(**kwargs)


if __name__ == "__main__":
    unittest.main()
