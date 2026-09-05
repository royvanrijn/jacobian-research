from fractions import Fraction as Q
from math import gcd,isqrt
from pathlib import Path
import json
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
import mw16_sensitivity_backend as s
import verify_mw16_sensitivity as verify
from freeze_mw16_sensitivity_setting import freeze
from build_mw16_sensitivity_gate import build as build_gate
from freeze_mw16_adaptive_policy import freeze as freeze_adaptive
from alternate_quartic_covers import alternate_cover,short_add


class SensitivityTests(unittest.TestCase):
    model=(Q(0),Q(0),Q(0),Q(-1),Q(1))

    def test_rational_pgl2_and_ordinate_map_against_original_cover(self):
        p=(Q(0),Q(1))
        for _ in range(3):
            p=short_add(self.model,p,(Q(0),Q(1)))
        for scaling in (Q(1),Q(7,13)):
            model=(0,0,0,-1/scaling**4,1/scaling**6)
            point=(p[0]/scaling**2,p[1]/scaling**3)
            chart=s.base.make_chart(model,point)
            for spec in ('gauss','red','metric:1/16','metric:16','raw',
                         'red:1,2,0','red:2,1,0','red:2,2,1'):
                matrix,f,scale=s.horizontal(chart,spec)
                a,b,c,d=s.multiply(chart.matrix,matrix)
                den,k,u=chart.denominator,chart.shift,chart.curve_scale
                raw_matrix=(Q(den*den*a+k*c,den)/u,Q(den*den*b+k*d,den)/u,c,d)
                expected=s.base.binary_transform(alternate_cover(model,point).coefficients,raw_matrix)
                self.assertEqual(tuple(x*u**4/den**2*scale**2 for x in expected),f)
                aa,bb,cc,dd=matrix
                for n in range(-9,10):
                    for q in range(1,10):
                        if gcd(n,q)!=1: continue
                        value=sum(v*n**i*q**(4-i) for i,v in enumerate(f))
                        if value<0 or isqrt(value)**2!=value: continue
                        for root in {isqrt(value),-isqrt(value)}:
                            actual=chart.map_point(aa*n+bb*q,cc*n+dd*q,Q(root)/scale)
                            upper,lower=raw_matrix[0]*n+raw_matrix[1]*q,c*n+d*q
                            if lower==0:
                                self.assertIsNone(actual)
                            else:
                                raw=(upper/lower,den*Q(root)/scale/(u**2*lower**2))
                                self.assertEqual(actual,alternate_cover(model,point).cover_point_to_curve(raw))

    def test_square_content_with_nonsquare_factor(self):
        k=37**5*101**3
        values=tuple(Q(2*k*k*c,13**2) for c in (1,2,3,4,5))
        f,scale=s.normalize(values)
        self.assertEqual(f,(2,4,6,8,10))
        self.assertEqual(tuple(v*scale**2 for v in values),f)

    def test_complete_sieve_variant_and_checkpoint_tamper(self):
        args=dict(model=self.model,points=[(Q(0),Q(1))],representative=[1],mask=1,
                  specification='red:2,2,1',height=30,seconds=2)
        with tempfile.TemporaryDirectory() as directory:
            record=s.checkpoint(directory,**args)
            f=tuple(map(int,record['coefficients']))
            brute=set()
            for n in range(-30,31):
                for d in range(1,31):
                    if gcd(n,d)!=1: continue
                    value=sum(v*n**i*d**(4-i) for i,v in enumerate(f))
                    if value>=0 and isqrt(value)**2==value: brute.add((n,d,isqrt(value)))
            if f[4]>=0 and isqrt(f[4])**2==f[4]: brute.add((1,0,isqrt(f[4])))
            self.assertEqual(set(tuple(map(int,p)) for p in record['primitive_square_hits']),brute)
            with patch.object(s,'run',side_effect=RuntimeError('not cached')):
                self.assertEqual(s.checkpoint(directory,**args),record)
                with self.assertRaises(RuntimeError): s.checkpoint(directory,**{**args,'height':31})
            path=next(Path(directory).glob('*.json')); saved=json.loads(path.read_text())
            saved['record']['coefficients'][0]='0'; path.write_text(json.dumps(saved))
            with self.assertRaises(ArithmeticError): s.checkpoint(directory,**args)

    def test_incomplete_checkpoint_is_replayed(self):
        args=dict(model=self.model,points=[(Q(0),Q(1))],representative=[1],mask=1,
                  specification='gauss',height=10000,seconds=1e-9)
        with tempfile.TemporaryDirectory() as directory:
            self.assertEqual(s.checkpoint(directory,**args)['status'],'bounded_search_timeout')
            with patch.object(s,'run',side_effect=RuntimeError('must replay')):
                with self.assertRaises(RuntimeError): s.checkpoint(directory,**args)

    def test_certificate_rejects_changed_polynomial_hit_and_coverage(self):
        record=s.run(model=self.model,points=[(Q(0),Q(1))],representative=[1],mask=1,
                     specification='red:2,2,1',height=30,seconds=2)
        verify.check_chart(self.model,[(Q(0),Q(1))],record)
        for key,value in (('coefficients',['0']*5),('finite_curve_points',[]),
                          ('completed_denominator',29),('representative',[2])):
            with self.assertRaises(ArithmeticError):
                verify.check_chart(self.model,[(Q(0),Q(1))],{**record,key:value})

    def test_selection_prefers_certified_rank_even_when_slower(self):
        rows=[]
        for curve,high,low in zip((398,400,401,542,548),(5,5,11,10,8),(4,4,10,10,8)):
            settings=[]
            for key,gain,seconds in (('sensitive',high,1000),('fast',low,0.001)):
                settings.append({'key':key,'height':100000,'centre':'generic',
                    'specification':key,'exact_quotient_rank_recovered':gain,'wall_seconds':seconds,
                    'classification':{'status':'PASS_BASIS_EQUALS_DISCOVERED_GROUP'},'current_basis':[]})
            rows.append({'curve_id':curve,'settings':settings})
        payload={'status':'COMPLETE','declared_budget':{'mode':'initial'},'results':rows}
        chosen=freeze(payload)
        self.assertEqual(chosen['total_control_quotient_rank'],39)
        self.assertEqual(chosen['selection_ranking'][0]['key'],'sensitive')
        with self.assertRaises(ArithmeticError): freeze({**payload,'results':rows[:-1]})
        with self.assertRaises(ArithmeticError): freeze({**payload,'status':'SEARCHING'})

    def test_prospective_gate_fails_closed_on_recovery_and_replay(self):
        rows=[{'curve_id':curve,'exact_quotient_rank_recovered':gain,'settings':[]}
              for curve,gain in zip((398,400,401,542,548),(14,12,10,10,8))]
        campaign={'status':'COMPLETE','results':rows}
        replay={'status':'PASS_EXACT_REPLAY'}
        self.assertEqual(build_gate(campaign,campaign,replay)['total_control_quotient_rank'],54)
        weaker={**campaign,'results':[{**r,'exact_quotient_rank_recovered':r['exact_quotient_rank_recovered']-1} for r in rows]}
        with self.assertRaises(ArithmeticError): build_gate(campaign,weaker,replay)
        with self.assertRaises(ArithmeticError): build_gate(campaign,campaign,{'status':'UNKNOWN'})

    def test_adaptive_trigger_uses_only_the_initial_recovered_rank(self):
        initial={'status':'COMPLETE','results':[{'parent_id':str(i),'exact_quotient_rank_recovered':g}
            for i,g in enumerate((0,5,7,8,11))]}
        trials={'status':'COMPLETE','declared_budget':{},'claim_boundary':[],
            'results':[{'parent_id':str(i),'exact_quotient_rank_recovered':12,'status':'COMPLETE'} for i in (1,2)]}
        result=freeze_adaptive(initial,trials,8)
        self.assertEqual(result['declared_budget']['active_adaptive_parent_ids'],['1','2'])
        self.assertEqual(result['total_control_quotient_rank'],43)
        with self.assertRaises(ArithmeticError): freeze_adaptive(initial,{**trials,'results':trials['results'][:1]},8)


if __name__=='__main__': unittest.main()
