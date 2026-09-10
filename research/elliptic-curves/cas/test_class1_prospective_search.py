import unittest
from run_class1_prospective_search import allowance,selection,primes,control_window,integral_projection
from fractions import Fraction

class PolicyTests(unittest.TestCase):
    def test_rational_parent_integral_marking(self):
        def rf(values):return {'numerator_coefficients_low_to_high':values,'denominator_coefficients_low_to_high':['1','1']}
        parent={'A_coefficients_low_to_high':['1/3','-2'], 'B_coefficients_low_to_high':['2/27','1/9'],
                'sections':[{'X':rf(['1/5','1']),'Y':rf(['-1/7','2'])}],'generic_height_gram':[[4]]}
        marked=integral_projection(parent);scale=marked['input_coordinate_scaling']['scale']
        self.assertEqual(scale,3)
        for field,power in [('A_coefficients_low_to_high',4),('B_coefficients_low_to_high',6)]:
            self.assertEqual([Fraction(c)/scale**power for c in marked[field]],list(map(Fraction,parent[field])))
        for coordinate,power in [('X',2),('Y',3)]:
            actual=marked['sections'][0][coordinate];original=parent['sections'][0][coordinate]
            self.assertEqual(actual['denominator_coefficients_low_to_high'],original['denominator_coefficients_low_to_high'])
            self.assertEqual([Fraction(c)/scale**power for c in actual['numerator_coefficients_low_to_high']],list(map(Fraction,original['numerator_coefficients_low_to_high'])))
        self.assertEqual(marked['generic_height_gram'],parent['generic_height_gram'])
        self.assertIs(integral_projection(marked),marked)
    def test_thresholds(self):
        self.assertEqual([allowance(r) for r in (17,19,20,22,23,25,27,31)],[24,24,256,256,1024,2048,8192,8192])
    def test_controls_ignore_score(self):
        rows=[dict(index=i,control=i%8==7,score_units=i,model_bits=1) for i in range(64)]
        ordered=selection(rows)
        self.assertEqual([r['index'] for r in ordered[7::8]],list(range(7,64,8)))
        self.assertEqual(len({r['index'] for r in ordered}),64)
        for r in rows:
            if r['control']:r['score_units']=10**20-r['index']
        self.assertEqual([r['index'] for r in selection(rows)[7::8]],list(range(7,64,8)))
    def test_ranked_stream_is_sorted(self):
        rows=[dict(index=i,control=i%8==7,score_units=i,model_bits=1) for i in range(64)]
        scores=[r['score_units'] for r in selection(rows) if not r['control']]
        self.assertEqual(scores,sorted(scores,reverse=True))
    def test_prime_roster(self):self.assertEqual(primes(20),[5,7,11,13,17,19])
    def test_signed_controls_are_exact_fraction_and_not_sign_locked(self):
        c=control_window(65536,65536)
        self.assertEqual(len(c),8192)
        negative=sum(i%2 for i in c)
        self.assertTrue(3500<negative<4700)
        rows=[dict(index=i,control=i in c,control_order=c.get(i),score_units=i,model_bits=1) for i in range(65536,131072)]
        ordered=selection(rows)
        self.assertEqual([r['index'] for r in ordered[7::8]],list(c))

if __name__=='__main__':unittest.main()
