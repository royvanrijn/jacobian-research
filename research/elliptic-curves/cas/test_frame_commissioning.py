import unittest
import tempfile
from pathlib import Path
from run_frame_commissioning import panel_rows, decision, prepare, run, write, read

class CommissioningTests(unittest.TestCase):
    def test_fixed_prefix_caps(self):
        rows=[dict(index=i,control=i%8==7,control_order=i,score_units=i,model_bits=1) for i in range(1024)]
        panel=panel_rows(rows,64,16)
        self.assertEqual(len(panel),80)
        self.assertEqual(sum(r['control'] for r in panel),16)
        self.assertEqual([r['index'] for r in panel if r['control']],list(range(7,128,8)))
        self.assertEqual([r['index'] for r in panel if not r['control']],sorted([r['index'] for r in rows if not r['control']],reverse=True)[:64])

    def test_failed_panel_never_releases_next_frame(self):
        self.assertEqual(decision([dict(status='UNKNOWN_TIMEOUT',rank_lower_bound=17)]),'UNKNOWN_INCOMPLETE_COMMISSIONING')
        self.assertEqual(decision([dict(status='COMPLETE_BOUNDED',rank_lower_bound=19)]),'NO_EVIDENCE_CURRENT_SEARCH_PRODUCTIVE')
        self.assertEqual(decision([dict(status='COMPLETE_BOUNDED',rank_lower_bound=20)]),'CERTIFIED_GE20_ESCALATED_REVIEW_BEFORE_NEXT_FRAME')

    def test_controller_caps_and_preserves_full_scores(self):
        # Synthetic worker tests scheduling only; it is not a rank certificate.
        for fake_rank in (17,20):
            with self.subTest(rank=fake_rank), tempfile.TemporaryDirectory() as tmp:
                folder=Path(tmp);rt=folder/'runtime';window=rt/'ordinary-search/window-000'
                write(folder/'plan.json',{'root':str(rt),'workers':2})
                rows=[dict(index=i,control=i%8==7,control_order=i,score_units=i,model_bits=1) for i in range(32)]
                write(window/'scores.json',{'rows':rows})
                write(window/'queue.json',{'indices':[r['index'] for r in panel_rows(rows,28,4)]})
                script=rt/'elliptic-curves/cas/run_class1_prospective_search.py'
                script.parent.mkdir(parents=True)
                script.write_text('''from pathlib import Path
import json
def guard(folder): pass
def fibre(folder,window,row):
    case=window/'cases'/str(row['index']);case.mkdir(parents=True)
    result={'index':row['index'],'control':row['control'],'status':'COMPLETE_BOUNDED',
            'rank_lower_bound':RANK,'calls':CALLS,'search_cpu_seconds':1}
    (case/'terminal.json').write_text(json.dumps(result))
    return result
'''.replace('RANK',str(fake_rank)).replace('CALLS','256' if fake_rank==20 else '24'))
                prepare(folder,4,2)
                run(folder)
                result=read(folder/'commissioning-result.json')
                self.assertEqual(len(list((window/'cases').iterdir())),6)
                self.assertEqual((result['ranked'],result['controls']),(4,2))
                self.assertEqual(result['next_frame_released'],fake_rank==17)
                import gzip,json
                frozen=json.loads(gzip.decompress((folder/'matched-baseline/scores.json.gz').read_bytes()))
                self.assertEqual(frozen['rows'],rows)

if __name__=='__main__':unittest.main()
