import sys
from pathlib import Path
import unittest

CAS=Path(__file__).resolve().parents[1]/'cas'
if str(CAS) not in sys.path: sys.path.insert(0,str(CAS))
from curve302_positive_exposure_tables import span_intersection_rank,saturation_contains,stage_positive_words,extract_all

class PositiveEvidenceTablesTest(unittest.TestCase):
    def test_intersection_rank(self):
        e1=(1,0,0); e2=(0,1,0); e3=(0,0,1)
        self.assertEqual(span_intersection_rank((e1,),(e1,e2),3),1)
        self.assertEqual(span_intersection_rank((e3,),(e1,e2),3),0)
        self.assertTrue(saturation_contains(((1,1,0),(1,-1,0)),(e1,e2),3))

    def test_saturation_enables_integral_axis_via_mixed_word(self):
        # Prefix <(1,1)> plus exposed (1,-1) has saturated span Z^2, hence e1/e2.
        self.assertTrue(saturation_contains(((1,1),(1,-1)),((1,0),),2))
        self.assertTrue(saturation_contains(((1,1),(1,-1)),((0,1),),2))

    def test_positive_word_index_counts_only_rationally_new(self):
        prefix=((1,0,0),)
        charts=({'chart_id':'c0','list_order':0,'order':0,'words':((1,0,0),(0,1,0))},
                {'chart_id':'c1','list_order':1,'order':1,'words':((0,1,0),(0,0,1))})
        got=stage_positive_words(charts,prefix,3)
        self.assertNotIn((1,0,0),got)
        self.assertEqual(got[(0,1,0)]['chart_count'],2)
        self.assertEqual(got[(0,0,1)]['chart_count'],1)
        self.assertEqual(got[(0,1,0)]['first_list_order'],0)

    def test_extract_all_positive_only_synthetic(self):
        n=14
        e=[tuple(int(i==j) for i in range(n)) for j in range(n)]
        data={
          "runs":[{"seed":"recovered-local-01","seed_index":0,
                   "stages":[{"epoch":0,"gains":[{"word":e[1],"primitive":e[1]}],"terminal_no_gain":False}],
                   "events":[{"word":e[1],"primitive":e[1],"epoch":0}],"final_dimension":2,"charts":2}],
          "cores":{1:(e[0],),2:(e[0],e[1]),5:tuple(e[:5]),6:tuple(e[:6])}
        }
        ledger={"schema":"curve302-chart-exposure-ledger.v1","direction_ids":[f"d{i}" for i in range(n)],
                "runs":[{"seed":"recovered-local-01","stages":[{"epoch":0,"charts":[
                  {"chart_id":"c0","order":0,"exposures":[{"word":list(e[2])}]},
                  {"chart_id":"c1","order":1,"exposures":[{"word":list(e[1])}]}]}]}]}
        out=extract_all(data,ledger,n)
        self.assertEqual(out["coexposure_choices"]["summary"]["gain_stages"],1)
        self.assertEqual(out["coexposure_choices"]["summary"]["stages_with_positive_alternatives"],1)
        self.assertEqual(out["coexposure_choices"]["summary"]["stages_alternative_precedes_all_actual"],1)
        self.assertEqual(out["saturation_impact"]["summary"]["positive_candidate_stage_pairs"],2)

    def test_initial_seed_containment_is_not_delayed(self):
        n=14
        e=[tuple(int(i==j) for i in range(n)) for j in range(n)]
        data={"runs":[{"seed":"recovered-local-01","seed_index":0,
                       "stages":[{"epoch":0,"gains":[{"word":e[1],"primitive":e[1]}]}],
                       "events":[{"word":e[1],"primitive":e[1],"epoch":0}],"final_dimension":2,"charts":1}],
              "cores":{1:(e[0],),5:tuple(e[:5]),6:tuple(e[:6])}}
        ledger={"schema":"curve302-chart-exposure-ledger.v1","direction_ids":[f"d{i}" for i in range(n)],
                "runs":[{"seed":"recovered-local-01","stages":[{"epoch":0,"charts":[{"chart_id":"c","order":0,"exposures":[{"word":list(e[1])}]}]}]}]}
        out=extract_all(data,ledger,n)["lead_times"]
        l1=next(r for r in out["runs"] if r["target"]=="L1")
        self.assertEqual(l1["first_integral_containment"]["stage_index"],-1)
        self.assertIsNone(l1["first_direct_positive"])

if __name__=='__main__': unittest.main()
