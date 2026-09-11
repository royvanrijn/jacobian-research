import json, os, sys, tempfile, unittest
from fractions import Fraction as F
from pathlib import Path

CAS=Path(__file__).resolve().parents[1]/'cas'
sys.path.insert(0,str(CAS))
import curve302_chart_replay_adapter as a

NAMES=tuple([f'recovered-local-{i:02d}' for i in range(1,5)] +
            [f'recovered-strict-{i:02d}' for i in range(1,4)] +
            [f'residual-strict-{i:02d}' for i in range(1,8)])

def p(t): return (F(t*t),F(t*t*t))
def pj(P): return {'x':str(P[0]),'y':str(P[1])}

class FakeRec:
    def __init__(self,ainv,m17,exceptional,**kw):
        self.map={a.point_text(P):tuple(1 if j==i else 0 for j in range(14)) for i,P in enumerate(exceptional)}
    def recognize(self,P):
        w=self.map.get(a.point_text(P))
        return None if w is None else a.Recognition(P,tuple([0]*17+list(w)),w,128)

class Tests(unittest.TestCase):
    def test_epoch_from_parent_path(self):
        self.assertEqual(a.context_epoch(Path('seed/replay-M17/epoch-07/cloud-003.json')),7)
        self.assertIsNone(a.context_epoch(Path('seed/no-epoch/cloud-003.json')))

    def test_point_parser(self):
        self.assertEqual(a.parse_point({'point':['4','8']}),(F(4),F(8)))
        self.assertEqual(a.parse_point('(9 : 27 : 1)'),(F(9),F(27)))
        self.assertIsNone(a.parse_point('(9 : 27 : 2)'))

    def test_infer_model(self):
        inv=a.infer_ainvariants([p(i) for i in range(1,18)])
        self.assertEqual(inv,(F(0),)*5)
        self.assertTrue(all(a.on_curve(p(i),inv) for i in range(1,30)))

    def make_raw(self,root, *, cumulative=False, direct_search_point=False):
        common=[p(i) for i in range(1,18)]
        extras=[]
        for si,seed in enumerate(NAMES):
            extra=p(100+si); extras.append(extra)
            d=Path(root)/seed/'replay-M17'/'epoch-00'; d.mkdir(parents=True)
            (d/'state.json').write_text(json.dumps({'schema':a.SCHEMA_MW,'basis':[pj(x) for x in common+[extra]]}))
            chart={'centre':[si,0],'index':0,'mapping':{'kind':'synthetic'},'search':{'complete':True}}
            if direct_search_point: chart['search']['found_point']=pj(extra)
            (d/'cloud-000.json').write_text(json.dumps({'charts':[chart]}))
            if cumulative:
                chart2={'centre':[si,1],'index':1,'mapping':{'kind':'synthetic'},'search':{'complete':True}}
                (d/'cloud-001.json').write_text(json.dumps({'charts':[chart,chart2]}))
            (d/'hit-index-0.json').write_text(json.dumps({'schema':a.SCHEMA_POINT,'index':0,'recorded_point':pj(extra),'complete':True}))
        return common,extras

    def test_derive_seeded_basis(self):
        with tempfile.TemporaryDirectory() as t:
            common,extras=self.make_raw(t)
            b=a.derive_seeded_basis(Path(t),NAMES)
            self.assertEqual(len(b['m17']),17); self.assertEqual(len(b['exceptional']),14)
            self.assertEqual(tuple(a.point_text(x) for x in b['exceptional']),tuple(a.point_text(x) for x in extras))

    def test_cumulative_snapshot_uses_largest(self):
        with tempfile.TemporaryDirectory() as t:
            self.make_raw(t,cumulative=True)
            charts=a.load_stage_charts(Path(t),NAMES,{s:{0} for s in NAMES})
            self.assertTrue(all(len(charts[(s,0)])==2 for s in NAMES))

    def test_join_recorded_hits(self):
        with tempfile.TemporaryDirectory() as t:
            _,extras=self.make_raw(t)
            charts=a.load_stage_charts(Path(t),NAMES,{s:{0} for s in NAMES})
            joined,unresolved=a.join_hits_to_charts(Path(t),NAMES,charts)
            self.assertFalse(unresolved)
            for i,s in enumerate(NAMES):
                self.assertEqual(a.point_text(joined[f'{s}:0:0'][0]['point']),a.point_text(extras[i]))

    def test_join_direct_search_points(self):
        with tempfile.TemporaryDirectory() as t:
            _,extras=self.make_raw(t,direct_search_point=True)
            charts=a.load_stage_charts(Path(t),NAMES,{s:{0} for s in NAMES})
            joined,unresolved=a.join_hits_to_charts(Path(t),NAMES,charts)
            self.assertFalse(unresolved)
            # direct-search + recorded file are deduplicated by exact point
            self.assertTrue(all(len(joined[f'{s}:0:0'])==1 for s in NAMES))

    def test_build_ledger_fake_exact_recognizer(self):
        with tempfile.TemporaryDirectory() as t:
            self.make_raw(t)
            led=a.build_replay_ledger(Path(t),NAMES,{s:{0} for s in NAMES},recognizer_factory=FakeRec)
            self.assertEqual(led['schema'],a.LEDGER_SCHEMA)
            self.assertEqual(led['adapter']['unique_recorded_points'],14)
            self.assertEqual(led['adapter']['recognized_exposures'],14)
            for i,run in enumerate(led['runs']):
                w=run['stages'][0]['charts'][0]['exposures'][0]['quotient_word']
                self.assertEqual(w,[1 if j==i else 0 for j in range(14)])

    def test_schema_audit(self):
        with tempfile.TemporaryDirectory() as t:
            self.make_raw(t)
            out=a.replay_schema_audit(Path(t),NAMES,{s:{0} for s in NAMES})
            self.assertEqual(out['status'],'PASS_REPLAY_SCHEMA_ADAPTER')
            self.assertEqual(out['stages'],14); self.assertEqual(out['charts'],14)
            self.assertEqual(out['basis']['m17'],17); self.assertEqual(out['basis']['exceptional'],14)

    def test_one_pass_index_and_ledger(self):
        with tempfile.TemporaryDirectory() as t:
            self.make_raw(t)
            expected={s:{0} for s in NAMES}
            idx=a.scan_replay_index(Path(t),NAMES,expected)
            self.assertEqual(len(idx['best_charts']),14)
            self.assertEqual(len(idx['hits']),14)
            audit=a.replay_schema_audit_index(idx,NAMES,expected)
            self.assertEqual(audit['status'],'PASS_REPLAY_SCHEMA_ADAPTER')
            led=a.build_replay_ledger_index(idx,NAMES,expected,recognizer_factory=FakeRec)
            self.assertEqual(led['adapter']['recognized_exposures'],14)

    def test_completion_negative_wins(self):
        self.assertFalse(a.completion_evidence({'complete':True,'timeout':True}))
        self.assertTrue(a.completion_evidence({'returncode':0,'status':'done'}))
        self.assertIsNone(a.completion_evidence({'foo':1}))

if __name__=='__main__': unittest.main()
