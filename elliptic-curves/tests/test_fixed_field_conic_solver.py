"""Positive witnesses and fail-closed replay for the bounded conic experiment."""
import copy
import gzip
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
import solve_fixed_field_conic as solver
from sage.all import matrix


class ConicSolverReplay(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.evidence=json.loads(gzip.decompress(solver.EVIDENCE.read_bytes()))
        cls.K,cls.a,cls.b,cls.C,cls.T=solver.context()

    def test_retained_candidates_and_separate_controls(self):
        s=solver.replay(self.evidence)
        self.assertEqual(s['exact_target_candidates'],25920)
        self.assertEqual(s['target_reconstruction_cells'],320)
        self.assertEqual(s['target_conic_points'],0)
        self.assertEqual(s['seeded_reconstruction_control_points'],1)
        self.assertEqual(s['small_field_norm_control_points'],1)
        self.assertEqual(s['genuine_higher_covers'],0)
        self.assertEqual(set(s['point_or_sha'].values()),{'UNKNOWN'})

    def test_zero_and_fabricated_global_points_rejected(self):
        for p in ([0,0,0],[1,0,0]):
            with self.assertRaises(AssertionError):
                solver.witness(self.K,self.a,self.b,p,self.C,self.T)

    def test_positive_point_map_exactly_checked(self):
        # Same large labelled a, explicitly manufactured b, exact norm point.
        D=matrix.diagonal(self.K,[1,-self.a,self.a-4])
        solver.witness(self.K,self.a,4-self.a,[2,1,1],D,matrix.identity(self.K,3))
        D[0,0]+=1
        with self.assertRaises(AssertionError):
            solver.witness(self.K,self.a,4-self.a,[2,1,1],D,matrix.identity(self.K,3))

    def test_wrong_local_lift_rejected(self):
        e=copy.deepcopy(self.evidence)
        cell=next(z for z in e['runs']['lift']['result']['events'] if z['phase']=='reconstruction_cell')
        cell['local_point']=[['0']*3,['0']*3,['1','0','0']]
        with self.assertRaises(AssertionError):
            solver.replay(e)

    def test_wrong_reconstruction_basis_rejected(self):
        e=copy.deepcopy(self.evidence)
        cell=next(z for z in e['runs']['lift']['result']['events'] if z['phase']=='reconstruction_cell')
        cell['reconstruction_basis'][0]=copy.deepcopy(cell['reconstruction_basis'][1])
        with self.assertRaises(AssertionError):
            solver.replay(e)

    def test_changed_field_label_rejected(self):
        e=copy.deepcopy(self.evidence)
        e['runs']['lift']['result']['field_polynomial'][0]='1'
        with self.assertRaises(AssertionError):
            solver.replay(e)

    def test_timeout_cannot_be_declared_a_point(self):
        e=copy.deepcopy(self.evidence)
        e['runs']['norm']['result']['status']='EXACT_POINT'
        with self.assertRaises(AssertionError):
            solver.replay(e)

    def test_control_cannot_be_promoted_to_target(self):
        e=copy.deepcopy(self.evidence)
        e['runs']['lift-control']['result']['target_mask']=1047173
        with self.assertRaises(AssertionError):
            solver.replay(e)

    def test_long_search_transcripts_and_boundaries(self):
        evidence=json.loads(gzip.decompress(solver.LONG_EVIDENCE.read_bytes()))
        result=solver.long_replay(evidence)
        self.assertEqual(result['deep_reconstruction_cells'],131072)
        self.assertEqual(result['deep_exact_candidates'],10616832)
        self.assertEqual(result['exhaustive_w1_residues'],22973)
        self.assertEqual(result['exhaustive_exact_candidates'],1860813)
        self.assertEqual(result['target_conic_points'],0)
        self.assertEqual(result['genuine_higher_covers'],0)

    def test_long_transcript_tampering_is_rejected(self):
        evidence=json.loads(gzip.decompress(solver.LONG_EVIDENCE.read_bytes()))
        cell=evidence['runs']['deep_2700']['result']['events'][2]
        self.assertEqual(cell['phase'],'deep_reconstruction_cell')
        cell['candidate_sha256']='0'*64
        with self.assertRaises(AssertionError):
            solver.long_replay(evidence)

    def test_long_timeout_cannot_be_reported_as_a_solution(self):
        evidence=json.loads(gzip.decompress(solver.LONG_EVIDENCE.read_bytes()))
        evidence['runs']['norm_2700']['result']['status']='EXACT_POINT'
        with self.assertRaises(AssertionError):
            solver.long_replay(evidence)


if __name__=='__main__':
    unittest.main()
