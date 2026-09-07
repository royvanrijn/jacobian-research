"""Positive witness replay must not restart CVP or the point census."""
import copy
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
import run_mw_search as command
from research_runtime.finite_reduction import ReductionCache
from research_runtime.store import FactStore
from research_runtime.cvp import VoronoiIterator

class LazySearchTests(unittest.TestCase):
    def request(self):
        return json.loads((Path(__file__).resolve().parents[1]/'data/runtime_mw_search_control.json').read_text())

    def test_portable_replay_checks_points_without_search_or_cvp(self):
        with TemporaryDirectory() as directory:
            root=Path(directory);first=ReductionCache(FactStore(root/'cache'))
            with patch.object(command,'reduction_cache',return_value=first),patch('pointed_quartic_search.reduction_cache',return_value=first):
                discovered=command.search_request(self.request(),root/'discovery.json')
            fresh=ReductionCache(FactStore(root/'empty'))
            with patch.object(command,'reduction_cache',return_value=fresh),patch('pointed_quartic_search.search_box',side_effect=AssertionError('re-sieve forbidden')),patch.object(VoronoiIterator,'next_holes',side_effect=AssertionError('CVP regeneration forbidden')):
                replay=command.search_request(self.request(),root/'replay.json',retained=discovered)
                self.assertEqual(replay['final_state'],discovered['final_state'])
                tampered=copy.deepcopy(discovered)
                tampered['charts'][0]['search']['coefficients'][0]=str(int(tampered['charts'][0]['search']['coefficients'][0])+1)
                with self.assertRaises(ArithmeticError):command.search_request(self.request(),root/'bad.json',retained=tampered)

    def test_budget_checkpoint_can_request_next_unseen_holes(self):
        with TemporaryDirectory() as directory:
            root=Path(directory);cache=ReductionCache(FactStore(root/'cache'))
            request=self.request();request['cvp_node_budget']=1
            with patch.object(command,'reduction_cache',return_value=cache),patch('pointed_quartic_search.reduction_cache',return_value=cache):
                partial=command.search_request(request,root/'partial.json')
                self.assertEqual(partial['status'],'CVP_BUDGET_REACHED')
                self.assertEqual(partial['charts'],[])
                request['cvp_node_budget']=1000;request['cvp_checkpoint']=partial['cvp_checkpoint']
                complete=command.search_request(request,root/'complete.json')
                self.assertEqual(len(complete['charts']),3)
                request['cvp_checkpoint']=complete['cvp_checkpoint']
                exhausted=command.search_request(request,root/'exhausted.json')
                self.assertEqual(exhausted['charts'],[])
                self.assertFalse(exhausted['full_bnf_requested'])

    def test_manifest_retains_incremental_states_and_replays_without_census(self):
        import run_pointed_quartic_search as jobs
        from research_runtime.search_state import reduction_cache
        with TemporaryDirectory() as directory:
            root=Path(directory)
            manifest={'schema':jobs.SCHEMA,'jobs':[{'id':str(i),'curve':[0,0,0,-7,10],
                'subgroup':[[1,2],[2,2]],'centre':{'coefficients':vector},'coordinate_policy':'metric:16',
                'height':20,'seconds':2} for i,vector in enumerate(([1,0],[0,1]))]}
            result=jobs.execute(manifest,root/'result.json',root/'charts')
            self.assertEqual(len(result['mw_states']),1)
            with patch('pointed_quartic_search.search_box',side_effect=AssertionError('census forbidden')):
                self.assertEqual(jobs.verify_jobs(result)['chart_count'],2)
                bad=copy.deepcopy(result);bad['results'][0]['state_admission']['certified_rank_gain']+=1
                with self.assertRaises(ArithmeticError):jobs.verify_jobs(bad)

if __name__=='__main__':unittest.main()
