"""Real policy/files/scoring; CAS-only integration is mocked, never a Sage claim."""
from fractions import Fraction as F
import gzip
import importlib
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
import broad_rank_policy as policy
import broad_rank_runtime as io
import run_broad_rank_search as runner


def parent():
    return {'family':'toy','generic_rank_lower_bound':17,
        'A_coefficients_low_to_high':['1','2','0','1','0','0','0','0','3'],
        'B_coefficients_low_to_high':['5','1','7','0','0','0','0','0','0','0','0','0','2'],
        'sections':[], 'generic_height_gram':[[str(int(i==j)*4) for j in range(17)] for i in range(17)]}


def state(rank=17, calls=98, **kwargs):
    return {'status':'CERTIFIED','rank':rank,'calls':calls,'stale_calls':0,'round':0,
            'rescue_used':False,'history':[],**kwargs}


class PolicyTests(unittest.TestCase):
    def test_addresses_prefix(self):
        self.assertEqual([(a,b) for i,a,b in policy.addresses(0,8)],
                         [(1,1),(-1,1),(1,2),(-1,2),(2,1),(-2,1),(1,3),(-1,3)])

    def test_address_random_access(self):
        allrows=list(policy.addresses(0,10000))
        for offset in (0,1,2,25,1357,9990):
            self.assertEqual(list(policy.addresses(offset,10)),allrows[offset:offset+10])

    def test_address_binary_oracle(self):
        for i,a,b in policy.addresses(131073,2048):
            x=y=1
            for bit in bin(i//2+1)[3:]:
                if bit=='0': y+=x
                else:x+=y
            self.assertEqual((a,b),(x if i%2==0 else -x,y))

    def test_unique_reduced(self):
        import math
        rows=list(policy.addresses(123,20000))
        self.assertEqual(len({F(a,b) for i,a,b in rows}),len(rows))
        self.assertTrue(all(b>0 and math.gcd(a,b)==1 for i,a,b in rows))

    def test_bad_intervals(self):
        for offset,count in ((-1,2),(True,2),(1,0),(2**40,1)):
            with self.assertRaises(ValueError):list(policy.addresses(offset,count))

    def test_control_is_score_independent(self):
        rows=[{'index':i,'score_units':i,'model_bits':20} for i in range(1000)]
        a=policy.select(rows,20,10)
        b=policy.select([{**r,'score_units':-r['score_units']} for r in rows[::-1]],20,10)
        self.assertEqual([r['index'] for r in a if r['arm']=='control'],
                         [r['index'] for r in b if r['arm']=='control'])
        self.assertTrue(any(r['index']%2 for r in a if r['arm']=='control'))
        self.assertTrue(any(r['index']%2==0 for r in a if r['arm']=='control'))

    def test_selection_disjoint_invariant(self):
        rows=[{'index':i,'score_units':i%3,'model_bits':20} for i in range(101)]
        self.assertEqual(policy.select(rows,32,8),policy.select(rows[::-1],32,8))
        chosen=policy.select(rows,32,8)
        self.assertEqual(len({r['index'] for r in chosen}),40)
        self.assertEqual(sum(r['arm']=='control' for r in chosen),8)

    def test_empty_arms(self):
        rows=[{'index':i,'score_units':i,'model_bits':20} for i in range(10)]
        self.assertEqual(len(policy.select(rows,3,0)),3)
        self.assertEqual(len(policy.select(rows,0,3)),3)
        with self.assertRaises(ValueError):policy.select(rows,10,1)
        with self.assertRaises(ValueError):policy.select(rows+rows[:1],1,1)

    def test_low_rank_seed_gets_amplification(self):
        for rank in (18,19):self.assertIsNotNone(policy.next_stage(state(rank=rank),False,96))

    def test_rescue_is_predeclared_once(self):
        self.assertIsNone(policy.next_stage(state(),False,96))
        self.assertTrue(policy.next_stage(state(),True,96)['revival'])
        self.assertGreaterEqual(policy.next_stage(state(),True,96)['bank_index'],4)
        self.assertIsNone(policy.next_stage(state(rescue_used=True),True,96))

    def test_caps_and_censoring(self):
        self.assertIsNone(policy.next_stage(state(rank=32),True,96))
        self.assertIsNone(policy.next_stage(state(rank=27,calls=8192),False,96))
        self.assertIsNone(policy.next_stage(state(rank=20,stale_calls=200),False,96))
        for status in ('UNKNOWN_TIMEOUT','CENSORED_WITH_CERTIFIED_LOWER_BOUND'):
            self.assertIsNone(policy.next_stage(state(rank=25,status=status),True,96))
        self.assertIsNone(policy.next_stage(state(rank=22,round=96),True,96))

    def test_budget_tail(self):
        self.assertEqual(policy.next_stage(state(rank=27,calls=8190),False,96)['allowance'],2)

    def test_prime_roster(self):
        self.assertEqual(policy.primes(19),[5,7,11,13,17,19])
        for n in (2,5000,True):
            with self.assertRaises(ValueError):policy.primes(n)

    def test_rational_scaling(self):
        p=parent();p['A_coefficients_low_to_high']=['1/2','2/3'];p['B_coefficients_low_to_high']=['3/5']
        A,B,d=policy.integral_score_model(p)
        self.assertEqual(policy.polynomial(A,F(2,7)),policy.polynomial(p['A_coefficients_low_to_high'],F(2,7))*d**4)
        self.assertEqual(B[0],F(3,5)*d**6)

    def test_j_scaling_and_twist(self):
        a=[0,0,0,-1,1];b=[0,0,0,-16,64];twist=[0,0,0,-4,8]
        self.assertEqual(policy.j_invariant(a),policy.j_invariant(b))
        self.assertEqual(policy.j_invariant(a),policy.j_invariant(twist))
        self.assertIsNone(policy.j_invariant([0]*5))

    def test_homogeneous_model(self):
        p=parent();t=F(-2,3);model=policy.model_at(p,t)
        self.assertEqual(F(model[3]),policy.homogeneous(list(map(F,p['A_coefficients_low_to_high'])),-2,3,8))
        self.assertEqual(F(model[4]),policy.homogeneous(list(map(F,p['B_coefficients_low_to_high'])),-2,3,12))


class ScoreTests(unittest.TestCase):
    def test_all_point_counts_independent(self):
        p=parent();tables=policy.score_tables(p,43);A,B,d=policy.integral_score_model(p)
        for tab in tables['tables']:
            prime=tab['prime'];k=tab['constant_p_scaling']
            aa=[c//prime**(4*k) for c in A];bb=[c//prime**(6*k) for c in B]
            for r,ap in enumerate(tab['traces']):
                av=policy.homogeneous(aa,r,1,8)%prime if r<prime else (aa+[0]*9)[8]%prime
                bv=policy.homogeneous(bb,r,1,12)%prime if r<prime else (bb+[0]*13)[12]%prime
                points=1+sum((y*y-x*x*x-av*x-bv)%prime==0 for x in range(prime) for y in range(prime))
                self.assertEqual(prime+1-points,ap)

    def test_feature_scalar_score(self):
        p=parent();tabs=policy.score_tables(p,43)
        rows=list(policy.feature_rows(p,tabs,121,44,chunk_size=7))
        for row in rows:
            t=F(row['parameter']);total=0;bad=[]
            for tab in tabs['tables']:
                prime=tab['prime'];r=t.numerator*pow(t.denominator,-1,prime)%prime if t.denominator%prime else prime
                total+=tab['score_units'][r]
                if not tab['smooth'][r]:bad.append(prime)
            self.assertEqual(row['score_units'],total);self.assertEqual(row['local_singular_primes'],bad)

    def test_chunk_invariance(self):
        p=parent();tabs=policy.score_tables(p,19)
        self.assertEqual(list(policy.feature_rows(p,tabs,13,103,7)),list(policy.feature_rows(p,tabs,13,103,64)))

    def test_exact_singular_is_not_local_singular(self):
        p=parent();p['A_coefficients_low_to_high']=['0'];p['B_coefficients_low_to_high']=['0','1']
        tabs=policy.score_tables(p,19)
        row=list(policy.feature_rows(p,tabs,0,1))[0]
        self.assertTrue(row['nonsingular'])
        self.assertFalse(tabs['tables'][0]['smooth'][0])


class IntegrityTests(unittest.TestCase):
    def test_immutable(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x';io.write(p,{'a':1});io.write(p,{'a':1})
            with self.assertRaises(ValueError):io.write(p,{'a':2})

    def test_freeze_tamper(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);root=p/'rt';root.mkdir();(root/'x.py').write_text('x=1')
            io.write(p/'plan.json',{'root':str(root)})
            io.write(p/'manifest.json',{'plan_sha256':io.sha(p/'plan.json'),'files':{'x.py':io.sha(root/'x.py')},'executables':{}})
            io.guard(p);(root/'x.py').write_text('x=2')
            with self.assertRaises(ValueError):io.guard(p)

    def test_receipt_tamper(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);io.write(p/'dispatch.json',{});io.write(p/'result.json',{'status':'PASS'})
            io.write(p/'seal.json',{'completed':True,'dispatch_sha256':io.sha(p/'dispatch.json'),
                                  'files':{'result.json':io.sha(p/'result.json')}})
            self.assertIsNotNone(io.completed(p))
            (p/'result.json').write_text('{}')
            with self.assertRaises(ValueError):io.completed(p)

    def test_failure_never_passes_on_stray_result(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);io.write(p/'dispatch.json',{});io.write(p/'result.json',{'status':'PASS'})
            io.write(p/'seal.json',{'completed':False,'dispatch_sha256':io.sha(p/'dispatch.json'),'files':{}})
            self.assertIsNone(io.completed(p)[1])

    def test_projection_drops_oracles(self):
        p=parent();p['sections']=[{'X':{},'Y':{}}]*17;p['exceptional_points']=['SECRET']
        q=io.normalize_parent(p,p['generic_height_gram'],'test')
        self.assertNotIn('exceptional_points',q)
        self.assertEqual(q['family'],'test')

    def test_short_model_projection(self):
        z={'numerator':['0'],'denominator':['1']}
        p={'a_invariants':[z,z,z,{'numerator':['2'],'denominator':['3']},z],
           'basis_weierstrass_coordinates':[[z,z]]*17}
        q=io.normalize_parent(p,parent()['generic_height_gram'],'test')
        self.assertEqual(q['A_coefficients_low_to_high'],['2/3'])
        p['a_invariants'][0]={'numerator':['1'],'denominator':['1']}
        with self.assertRaises(ValueError):io.normalize_parent(p,parent()['generic_height_gram'],'test')


class ControllerTests(unittest.TestCase):
    def fixture(self,folder):
        root=folder/'runtime/research';p=parent();p['family']='x1092-class1'
        io.write(root/'broad-inputs/parents/toy.json',p)
        io.write(root/'broad-inputs/exclusions.json',{'models':[]})
        plan={'root':str(root),'sage':'/mock/sage','parents':[{'id':'x1092-class1','backend':'generic',
              'path':'broad-inputs/parents/toy.json','sha256':io.sha(root/'broad-inputs/parents/toy.json'),'ranked':3,'controls':2}],
              'offset':1234,'window':32,'prime_bound':13,'workers':2,'max_rounds':0,'min_free_gib':0,
              'job_wall_seconds':3,'rss_bytes':100000,'height':125000,'consecutive_failure_limit':4}
        io.write(folder/'plan.json',plan)
        io.write(folder/'manifest.json',{'plan_sha256':io.sha(folder/'plan.json'),
            'files':{str(x.relative_to(root)):io.sha(x) for x in root.rglob('*') if x.is_file()},'executables':{}})
        calls=[]
        def execute(folder,job,command,seconds):
            if (job/'seal.json').exists():return io.completed(job)
            calls.append(command);job.mkdir(parents=True,exist_ok=True)
            io.write(job/'dispatch.json',{'command':command});io.write(job/'dispatch-intent.json',{})
            if '_score' in command:runner.score_task(folder,job,'x1092-class1')
            elif '_preflight' in command:io.write(job/'result.json',{'status':'PASS_INTERFACE_PREFLIGHT'})
            else:
                req=io.read(job/'request.json');model=policy.model_at(p,req['parameter'])
                pts=[[str(i),str(i+1)] for i in range(17)] # Fake CAS only, not a point/rank proof.
                io.write(job/'packet.json',{'curve':model,'points':pts,'rank_lower_bound':17})
                digest=io.sha(job/'packet.json')
                io.write(job/'verified.json',{'status':'PASS_TWO_FINITE_IMPLEMENTATIONS','packet_sha256':digest})
                io.write(job/'result.json',{'status':'PASS_CERTIFIED_PARENT_EVALUATION','family':'x1092-class1',
                    'parameter':req['parameter'],'request_sha256':io.sha(job/'request.json'),'calls':req['allowance'],
                    'rank_lower_bound':17,'packet_sha256':digest,'segments':[],'gain_timeline':[],'exposure_complete':True})
            files={str(x.relative_to(job)):io.sha(x) for x in job.rglob('*') if x.is_file()}
            io.write(job/'seal.json',{'completed':True,'dispatch_sha256':io.sha(job/'dispatch.json'),
                'files':files,'cpu_seconds':1.0,'outcome':'completed','returncode':0})
            return io.completed(job)
        return root,calls,execute

    def modules(self):
        return {'certify_compact_r17_candidates':types.SimpleNamespace(isomorphic=lambda a,b:a==b)}

    def test_complete_resume_uses_no_new_jobs(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);root,calls,execute=self.fixture(folder)
            with patch.object(runner,'execute',execute),patch.dict(sys.modules,self.modules()):
                runner.run(folder);count=len(calls);runner.run(folder)
            self.assertEqual(count,len(calls));self.assertEqual(count,7)
            self.assertTrue((folder/'COMPLETE.json').exists())
            report=io.read(folder/'REPORT.json')
            self.assertEqual(report['groups']['x1092-class1/control']['finished'],2)
            self.assertEqual(report['groups']['x1092-class1/ranked']['point_calls'],3*198)
            with gzip.open(root/'broad-pools/x1092-class1/features.jsonl.gz','rt') as stream:
                self.assertEqual(len(stream.readlines()),32)

    def test_stop_before_dispatch(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);root,calls,execute=self.fixture(folder);(folder/'STOP').touch()
            with patch.object(runner,'execute',execute):runner.run(folder)
            self.assertEqual(calls,[])

    def test_cache_tampering_rejected_on_resume(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);root,calls,execute=self.fixture(folder)
            with patch.object(runner,'execute',execute),patch.dict(sys.modules,self.modules()):runner.run(folder)
            packet=next(root.glob('broad-cases/*/batch-000/packet.json'));packet.write_text('{}')
            with patch.object(runner,'execute',execute),patch.dict(sys.modules,self.modules()):
                with self.assertRaises(ValueError):runner.run(folder)

    def test_unknown_generic_has_no_fake_rank(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);root,calls,execute=self.fixture(folder)
            def unknown(folder,job,command,seconds):
                result=execute(folder,job,command,seconds)
                if 'parent_foundry_worker.py' in ' '.join(command):
                    (job/'result.json').unlink();io.write(job/'result.json',{'status':'UNRESOLVED_GENERIC_SPECIALIZATION'})
                    seal=io.read(job/'seal.json');seal['files']['result.json']=io.sha(job/'result.json')
                    io.write(job/'seal.json',seal,False);return io.completed(job)
                return result
            with patch.object(runner,'execute',unknown),patch.dict(sys.modules,self.modules()):runner.run(folder)
            states=[io.read(p) for p in root.glob('broad-cases/*/state.json')]
            self.assertTrue(all(r['rank'] is None for r in states))
            self.assertTrue(all(r['status'].startswith('UNKNOWN') for r in states))

    def test_request_adapters(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);root,calls,execute=self.fixture(folder);plan=io.read(folder/'plan.json')
            plan['parents'].append({'id':'074d9','path':'unused','sha256':'abc'})
            row={'backend':'native','parent_id':'074d9','id':'test','family':'074d9','parameter':'1','model':['0','0','0','1','2']}
            req,cmd=runner.request_for(plan,row,policy.next_stage(None,False,96),None,root/'job')
            self.assertEqual(req['allowance'],100);self.assertEqual(req['candidate']['rank'],17)
            self.assertTrue(cmd[-2].endswith('high_rank_foundry_job.py'))
            old=state(rank=18,packet='old.json',packet_sha256='123',head=None)
            req,cmd=runner.request_for(plan,row,policy.next_stage(old,False,96),old,root/'job2')
            self.assertEqual(req['candidate']['packet_sha256'],'123')
            self.assertEqual(req['kind'],'continuation')



class CertificateAdapterTests(unittest.TestCase):
    def packet(self,job,*,native=False,rank=18):
        row={'backend':'native' if native else 'generic','family':'toy','parameter':'2','model':['0','0','0','1','2']}
        io.write(job/'request.json',{})
        io.write(job/'packet.json',{'curve':row['model'],'points':[[str(i),'1'] for i in range(rank)],'rank_lower_bound':rank})
        digest=io.sha(job/'packet.json')
        io.write(job/('packet-verified.json' if native else 'verified.json'),
                 {'status':'PASS_TWO_FINITE_IMPLEMENTATIONS','packet_sha256':digest})
        result={'status':'PASS_CERTIFIED_SEARCH' if native else 'PASS_CERTIFIED_PARENT_EVALUATION',
            'family':'toy','parameter':'2','request_sha256':io.sha(job/'request.json'),
            'packet_sha256':digest,'rank_lower_bound':rank,'calls':100,
            'gain_timeline':[{'before':17,'after':18,'call':5}] if rank==18 else [],
            'point_timeouts':0,'map_timeouts':0,'unresolved_cloud':False,
            'segments':[],'exposure_complete':True}
        return row,result,policy.next_stage(None,False,96)

    def test_native_and_generic_gain(self):
        for native in (False,True):
            with tempfile.TemporaryDirectory() as d:
                p=Path(d);row,result,spec=self.packet(p,native=native)
                state=io.absorb(None,result,row,spec,p,p)
                self.assertEqual(state['rank'],18);self.assertEqual(state['stale_calls'],95)

    def test_censored_lower_bound_is_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);row,result,spec=self.packet(p,native=True);result['point_timeouts']=1
            state=io.absorb(None,result,row,spec,p,p)
            self.assertEqual(state['rank'],18)
            self.assertEqual(state['status'],'CENSORED_WITH_CERTIFIED_LOWER_BOUND')

    def test_unaccounted_gain_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);row,result,spec=self.packet(p);result['gain_timeline']=[]
            with self.assertRaises(ValueError):io.absorb(None,result,row,spec,p,p)

    def test_wrong_curve_and_overbudget_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);row,result,spec=self.packet(p)
            with self.assertRaises(ValueError):io.absorb(None,{**result,'calls':999},row,spec,p,p)
            with self.assertRaises(ValueError):io.absorb(None,result,{**row,'parameter':'3'},spec,p,p)

    def test_derived_state_tamper(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);fixture=ControllerTests();root,calls,execute=fixture.fixture(p)
            with patch.object(runner,'execute',execute),patch.dict(sys.modules,fixture.modules()):runner.run(p)
            path=next(root.glob('broad-cases/*/batch-000/broad-state.json'));path.write_text('{}')
            with patch.object(runner,'execute',execute),patch.dict(sys.modules,fixture.modules()):
                with self.assertRaises(ValueError):runner.run(p)


if __name__=='__main__':unittest.main()
