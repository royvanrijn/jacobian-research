"""Cheap tests: exact adapters, orbit enumeration, and fail-closed replay gates.

These tests do NOT run the Sage/PARI research workload or establish new ranks.
"""
import copy
import gzip
import itertools
import json
import sys
import tempfile
import unittest
from fractions import Fraction as F
from pathlib import Path

CAS = Path(__file__).resolve().parents[1]/'cas'
sys.path.insert(0,str(CAS))
import v3_transfer_contract as c
from v3_transfer_orbits import shell_bank


class AdapterTests(unittest.TestCase):
    def test_literal_export_and_no_execution(self):
        text = """E = EllipticCurve(QQ, [0,0,0,-2,1])
points = [E([QQ(x),QQ(y)]) for x,y in [['0','1']]]
raise RuntimeError('must not execute')
"""
        model, points = c.extract_sage_literals(text)
        self.assertEqual(points,((F(0),F(1)),))
        self.assertTrue(c.on_curve(model,points[0]))

    def test_point_coordinates_export(self):
        model, points = c.extract_sage_literals("E=EllipticCurve(QQ,['0','0','0','-2','1'])\npoint_coordinates=[['0','1']]\n")
        self.assertEqual(len(model),5)
        self.assertEqual(len(points),1)

    def test_multiple_curves_rejected(self):
        with self.assertRaises(ValueError):
            c.extract_sage_literals("E=EllipticCurve(QQ,[0,0,0,-2,1])\nF=EllipticCurve(QQ,[0,0,0,-2,1])\npoint_coordinates=[['0','1']]\n")

    def test_nonliteral_input_rejected(self):
        with self.assertRaises((ValueError,TypeError)):
            c.extract_sage_literals("E=EllipticCurve(QQ,load_model())\npoint_coordinates=[['0','1']]\n")

    def test_off_curve_rejected(self):
        with self.assertRaises(ValueError):
            c.extract_sage_literals("E=EllipticCurve(QQ,[0,0,0,-2,1])\npoint_coordinates=[['0','2']]\n")

    def test_exact_transport(self):
        model=(0,0,0,-2,1)
        pts, proof=c.transport_to_short(model,model,[(0,1)])
        self.assertEqual(pts,((F(0),F(1)),));self.assertEqual(proof['scale'],'6')

    def test_general_model_transport(self):
        # (0,1) lies on y²+xy+y=x³+x²-2x+2.
        model=(1,1,1,-2,2)
        a1,a2,a3,a4,a6=map(F,model)
        b2=a1*a1+4*a2;b4=2*a4+a1*a3;b6=a3*a3+4*a6
        target=(0,0,0,-27*(b2*b2-24*b4),-54*(-b2**3+36*b2*b4-216*b6))
        pts,_=c.transport_to_short(model,target,[(0,1)])
        self.assertTrue(c.on_curve(target,pts[0]))

    def test_twist_not_merely_j_match(self):
        with self.assertRaises(ValueError):
            c.transport_to_short((0,0,0,-2,1),(0,0,0,-8,8),[(0,1)])

    def test_nonshort_rejected(self):
        with self.assertRaises(ValueError):
            c.transport_to_short((0,0,0,-2,1),(1,0,0,-2,1),[(0,1)])

    def test_squared_rational(self):
        self.assertEqual(c.sqrt_rational(F(4,9)),F(2,3))
        with self.assertRaises(ValueError):c.sqrt_rational(F(2))

    def test_immutable_publish(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'data.json';c.write_new(p,{'a':1})
            with self.assertRaises(FileExistsError):c.write_new(p,{'a':2})
            self.assertEqual(c.read_json(p),{'a':1})
            self.assertEqual(len(list(Path(td).iterdir())),1)

    def test_gzip_fallback(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'x.json'
            Path(str(p)+'.gz').write_bytes(gzip.compress(b'{"a": 2}'))
            self.assertEqual(c.read_json(p),{'a':2})

    def test_path_escape(self):
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaises(ValueError):c.within(Path(td),'../escape')

    def test_roster_scope(self):
        self.assertEqual([r['initial_rank'] for r in c.CASE_SPECS],[17,27,27,27])
        self.assertEqual(len({r['id'] for r in c.CASE_SPECS}),4)
        self.assertNotIn('new-20260906-188',[r.get('inventory_id') for r in c.CASE_SPECS])


class OrbitTests(unittest.TestCase):
    def test_exact_small_bank(self):
        r=shell_bank([[4,0],[0,4]],[[1,0],[0,1]])
        self.assertEqual(r['norm_counts'],{0:1,4:4,8:4})
        self.assertEqual(r['rows'],[{'mask':3,'norm':8,'word':[1,-1]}])

    def test_shear_transports_actual_cosets(self):
        # U=[[1,1],[0,1]], H=U*(4I)*U^t.
        a=shell_bank([[4,0],[0,4]],[[1,0],[0,1]])
        b=shell_bank([[8,4],[4,4]],[[1,1],[0,1]])
        self.assertEqual(a['rows'],b['rows'])
        self.assertEqual(a['norm_counts'],b['norm_counts'])

    def test_against_independent_bruteforce(self):
        for G in ([[4,1],[1,4]],[[4,1,0],[1,4,1],[0,1,4]],[[6,2],[2,4]]):
            n=len(G);best={};counts={}
            for w in itertools.product(range(-4,5),repeat=n):
                q=sum(w[i]*G[i][j]*w[j] for i in range(n) for j in range(n))
                if q>10:continue
                counts[q]=counts.get(q,0)+1
                mask=sum((x%2)<<j for j,x in enumerate(w))
                if mask:best[mask]=min(best.get(mask,100),q)
            r=shell_bank(G,[[int(i==j) for j in range(n)] for i in range(n)])
            self.assertEqual(r['norm_counts'],counts)
            self.assertEqual(r['all_represented_parity_minima'],{str(k):v for k,v in best.items()})

    def test_node_limit_no_partial_answer(self):
        with self.assertRaisesRegex(RuntimeError,'NODE_LIMIT'):
            shell_bank([[4,0],[0,4]],[[1,0],[0,1]],node_limit=1)

    def test_indefinite_rejected(self):
        with self.assertRaises(ValueError):shell_bank([[0,1],[1,0]],[[1,0],[0,1]])

    def test_wrong_parent_word_rejected(self):
        G=[[4*int(i==j) for j in range(17)] for i in range(17)]
        row={'mask':3,'norm':8,'word':[1,1]+[0]*15}
        c.check_orbit_rows(G,[row])
        wrong=copy.deepcopy(row);wrong['word'][0]=2
        with self.assertRaises(ValueError):c.check_orbit_rows(G,[wrong])
        with self.assertRaises(ValueError):c.check_orbit_rows(G,[row,row])


class StopTests(unittest.TestCase):
    def status(self,**kw):
        d=dict(rank=27,target=32,charts=10,max_charts=20,executed=10,scheduled=10,censored=0,gain=False)
        return c.classify_stop(**(d|kw))
    def test_no_gain_is_finite(self):self.assertEqual(self.status(),'COMPLETE_FINITE_NO_GAIN')
    def test_budget_is_not_absence(self):self.assertEqual(self.status(charts=20),'CHART_BUDGET_EXHAUSTED')
    def test_censored_is_not_negative(self):self.assertEqual(self.status(censored=1),'CENSORED_SEARCH')
    def test_gain_rebuilds(self):self.assertEqual(self.status(rank=28,gain=True),'REBUILD_AFTER_CERTIFIED_GAIN')
    def test_target_is_lower_bound(self):self.assertEqual(self.status(rank=32,gain=True),'TARGET_LOWER_BOUND_REACHED')
    def test_unfinished_epoch(self):self.assertEqual(self.status(executed=9),'INCOMPLETE_EPOCH')


class GateTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.d=self.root/'artifacts/local/elliptic-curves/adaptive-visibility-cascade-v3'
        self.d.mkdir(parents=True)
        self.cas=self.root/'elliptic-curves/cas';self.cas.mkdir(parents=True)
        for name in ('check_visibility_cascade_v3.sage','check_visibility_metric_v3.sage','engine.py'):
            (self.cas/name).write_text('# test fixture only\n')
        (self.d/'generic.json').write_text('{}')
        policy={'sources':{'elliptic-curves/cas/engine.py':c.sha(self.cas/'engine.py')},
                'inputs':{str((self.d/'generic.json').relative_to(self.root)):c.sha(self.d/'generic.json')}}
        self.save(self.d/'protocol.json',policy);ph=c.sha(self.d/'protocol.json')
        stages=[];rs=[];ms=[]
        for i,(before,after) in enumerate(((17,24),(24,31))):
            wd=self.d/f'replay-M17/epoch-{i:02d}';wd.mkdir(parents=True)
            self.save(wd/'selection.json',{'basis_size':before});self.save(wd/'mod2-000.json',{'rank_lower_bound':after})
            stage={'epoch':i,'before':before,'after':after,'charts':1,'audit':'mod2-000.json',
                   'audit_sha256':c.sha(wd/'mod2-000.json'),'full_cosets_scored':32}
            self.save(wd/'stage.json',stage);stages.append(stage)
            rs.append({'epoch':i,'before':before,'after':after,'charts_replayed':1,'full_cosets_replayed':32,
                       'independent_modl_ranks':{'3':after,'5':after},
                       'checkpoint_hashes':{str(x.relative_to(self.root)):c.sha(x) for x in wd.iterdir()}})
            ms.append({'rank':before,'selection':str((wd/'selection.json').relative_to(self.root)),
                       'sha256':c.sha(wd/'selection.json')})
        self.t=self.d/'replay-M17/terminal.json';self.r=self.d/'final-replay.json';self.m=self.d/'metric-replay-M17.json'
        self.save(self.t,{'status':'COMPLETE_BOUNDED_CALIBRATION','final_rank_lower_bound':31,'charts':2,'stages':stages,'protocol_sha256':ph})
        self.save(self.r,{'schema':'visibility-cascade-v3-replay','rank_lower_bound':31,'protocol_sha256':ph,
                         'checker_sha256':c.sha(self.cas/'check_visibility_cascade_v3.sage'),'stages':rs,'guard_rejected_unredacted_parent':True})
        self.save(self.m,{'protocol_sha256':ph,'checker_sha256':c.sha(self.cas/'check_visibility_metric_v3.sage'),'stages':ms})
    def tearDown(self):self.tmp.cleanup()
    def save(self,p,obj):p.write_text(json.dumps(obj))
    def gate(self):return c.validate_gate(self.root,self.d,self.r,self.m)
    def mutate(self,p,fn):
        d=c.read_json(p);fn(d);self.save(p,d)
    def test_complete_gate(self):self.assertEqual(self.gate()['rank_lower_bound'],31)
    def test_progress_file_not_gate(self):
        self.mutate(self.r,lambda d:d.pop('schema'))
        with self.assertRaisesRegex(ValueError,'final full replay'):self.gate()
    def test_partial_replay_not_gate(self):
        self.mutate(self.r,lambda d:d.update(rank_lower_bound=29))
        with self.assertRaises(ValueError):self.gate()
    def test_missing_odd_proof(self):
        self.mutate(self.r,lambda d:d['stages'][1]['independent_modl_ranks'].pop('5'))
        with self.assertRaisesRegex(ValueError,'odd-prime'):self.gate()
    def test_changed_checkpoint(self):
        (self.d/'replay-M17/epoch-00/selection.json').write_text('{"basis_size":99}')
        with self.assertRaises(ValueError):self.gate()
    def test_new_unbound_epoch_file(self):
        (self.d/'replay-M17/epoch-00/late.json').write_text('{}')
        with self.assertRaisesRegex(ValueError,'full epoch'):self.gate()
    def test_changed_source(self):
        (self.cas/'engine.py').write_text('changed')
        with self.assertRaisesRegex(ValueError,'Changed bound'):self.gate()
    def test_missing_metric_epoch(self):
        self.mutate(self.m,lambda d:d['stages'].pop())
        with self.assertRaisesRegex(ValueError,'coverage'):self.gate()
    def test_wrong_start_rank(self):
        self.mutate(self.r,lambda d:d['stages'][0].update(before=30))
        with self.assertRaisesRegex(ValueError,'rooted at M17'):self.gate()


if __name__=='__main__':unittest.main()
