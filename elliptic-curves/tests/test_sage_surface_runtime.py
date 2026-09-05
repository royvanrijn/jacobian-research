"""Portable toric output replay and regulator reuse without a Frobenius census."""
from copy import deepcopy
from pathlib import Path
import json
import sys
from tempfile import TemporaryDirectory
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from research_runtime.regulator import Surface,pre_search_gate
from research_runtime.surface_repository import SurfaceProofRepository
from research_runtime.toric_surface import replay_toric
from research_runtime.store import FactStore


class SurfaceReplayTests(unittest.TestCase):
    def setUp(self):
        self.request=json.loads((Path(__file__).resolve().parents[1]/'data/runtime_surface_control.json').read_text())
        self.surface=Surface(**self.request['surface'])

    def test_portable_exact_two_prime_regulator(self):
        with TemporaryDirectory() as directory:
            repository=SurfaceProofRepository(FactStore(directory))
            for proof in self.request['reductions']:
                repository.retain(self.surface,proof,verify=replay_toric)
            reductions=repository.replay(self.surface)
            gate=pre_search_gate(self.surface,reductions,candidate_rank=1)
            self.assertEqual(gate['two_prime_regulator_test'],'INCOMPATIBLE')
            self.assertEqual(gate['comparisons'][0]['ratio'],'141/77')
            self.assertEqual(gate['arithmetic_rank_upper'],0)
            self.assertEqual(len(repository.packets(self.surface)),2)

    def test_wrong_model_polynomial_moments_and_driver_rejected(self):
        for location in ['surface_key','coefficients','moments','raw_output','provenance']:
            proof=deepcopy(self.request['reductions'][0])
            if location in ['surface_key','raw_output']:proof[location]='changed'
            elif location=='provenance':proof[location]['nondegenerate_driver_completed']=False
            else:proof[location][0]='0'
            with self.assertRaises((ArithmeticError,ValueError)):
                replay_toric(self.surface,proof)


if __name__=='__main__':unittest.main()
