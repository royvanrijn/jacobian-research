"""Synthetic regressions. These do NOT reproduce the user's new real-data run."""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

CAS = Path(__file__).resolve().parents[1]/"cas"
sys.path.insert(0, str(CAS))
import run_curve302_closure_followup as lab


def synthetic_costs(n):
    rows = []
    for mask in range(1 << n):
        values = []
        for j in range(n):
            if mask >> j & 1:
                values.append(None)
            else:
                # An empty seed needs threshold 12; supplied seeds permit
                # a cheaper chain. All costs monotone under inclusion.
                needed = 1 << ((j-1) % n)
                values.append(2+j if mask & needed else 12+j)
        rows.append(values)
    return rows


def fixture(folder, n=14):
    folder.mkdir(parents=True, exist_ok=True)
    names = lab.NAMES if n == 14 else tuple(f"test-{i}" for i in range(n))
    costs = synthetic_costs(n)
    source = {"schema": "elliptic-curves.curve302-exceptional-subgroup-landscape.v1", "status": "PASS_RETROSPECTIVE_FINITE_SUBGROUP_LANDSCAPE",
              "synthetic_fixture": True,
              "fixed_basis": {"direction_ids": names, "generic_rank": 17, "rounded_metric_scale": 1000000},
              "subset_states": [{"state_mask": i, "retained_numerators": v} for i,v in enumerate(costs)]}
    strict = tuple(1 << i for i in range(4, n))
    form = [[str((i+1 if i == j else 0) + lab.F(1, 5)) for j in range(n)] for i in range(n)]
    relations = {"schema": "curve302-quotient-near-relations.v1", "status": "PASS_QUOTIENT_RELATION_ANALYSIS", "direction_ids": names,
                 "schur_quotient": form, "rounded_metric_scale": 1000000, "strict_quotient_rref_bitmasks": strict}
    runs = []
    for seed in range(n):
        acquired = [i for i in range(n) if i != seed]
        if seed == 0:
            acquired = acquired[:-2]
        masks = [1 << seed]
        stages = []
        for epoch, j in enumerate(acquired):
            word = [0]*n
            word[j] = 1
            word[seed] = 1  # mixed vectors, not simply coordinate-axis paths
            mask = lab.parity(word)
            masks.append(mask)
            entry = {"denominator": 1, "quotient_word": word, "integral": True,
                     "primitive_quotient_word": list(lab.primitive(word)), "mod2_bitmask": mask,
                     "mod2_strict": lab.in_span(mask, strict, n)}
            stages.append({"epoch": epoch, "before": 18+epoch, "after": 19+epoch, "new": [entry],
                           "quotient_mod2_rank": epoch+2, "quotient_subspace_rref": str(lab.rref(masks,n))})
        runs.append({"seed": names[seed], "final_rank": 18+len(acquired), "charts": 100+seed,
                     "stop_reason": "SYNTHETIC_FIXTURE", "stages": stages})
    trajectories = {"schema": "curve302-actual-v3-closure-trajectories.v1", "status": "PASS_ALL_14_SEEDED_TRAJECTORIES_RECONCILED",
                    "direction_ids": names, "runs": runs}
    lab.atomic(folder/"landscape.json", source)
    laws = {"schema": "curve302-closure-laws.v1", "status": "PASS_EXHAUSTIVE_16384_CLOSURE_ANALYSIS", "direction_ids": names,
            "input_sha256": lab.digest(folder/"landscape.json")}
    for name, obj in (("closure-laws", laws), ("quotient-relations", relations), ("trajectories", trajectories)):
        lab.atomic(folder/f"{name}.json", obj)
    lab.atomic(folder/"REPORT.json", {"status": "PASS_THREE_CLOSURE_EXPERIMENTS", "outputs": {
        name: lab.digest(folder/f"{name}.json") for name in ("closure-laws", "quotient-relations", "trajectories")}})
    return folder


def span_set(rows):
    out = {0}
    for v in rows:
        out |= {x^v for x in list(out)}
    return out


class BinaryTests(unittest.TestCase):
    def test_reference_rref(self):
        self.assertEqual(lab.rref([3,5,6],3), (5,6))

    def test_rref_permutations(self):
        for rows in ([3,5,6], [12,7,8], [1,1,0]):
            self.assertEqual(lab.rref(rows,4), lab.rref(list(reversed(rows)),4))

    def test_exhaustive_intersections(self):
        rng = random.Random(19)
        for _ in range(100):
            a, b = [rng.randrange(16) for _ in range(3)], [rng.randrange(16) for _ in range(2)]
            meet = lab.intersection(a,b,4)
            self.assertEqual(span_set(meet), span_set(a)&span_set(b))
            for v in meet:
                word = lab.span_word(v,a,4)
                result = 0
                for i in word: result ^= a[i]
                self.assertEqual(result,v)

    def test_orthogonality(self):
        for row in lab.orthogonal([3,5],4):
            self.assertTrue(all((row & v).bit_count()%2 == 0 for v in [3,5]))

    def test_invalid_mask(self):
        with self.assertRaises(ValueError): lab.rref([16],4)
        with self.assertRaises(ValueError): lab.rref([True],4)

    def test_primitive_sign_and_gcd(self):
        self.assertEqual(lab.primitive([-6,0,9]), (2,0,-3))
        self.assertEqual(lab.primitive([0,0,0]), (0,0,0))


class PersistenceTests(unittest.TestCase):
    def test_all_cutoffs_against_independent_reachability(self):
        costs = synthetic_costs(4)
        for start in range(16):
            arrivals, full, _ = lab.minimax_arrivals(costs,start)
            for threshold in range(17):
                reachable = {start}
                for mask in range(16):
                    if mask not in reachable: continue
                    for j in range(4):
                        if not mask >> j & 1 and costs[mask][j] <= threshold:
                            reachable.add(mask | (1 << j))
                union = 0
                for mask in reachable: union |= mask
                self.assertEqual(union, lab.closure(costs,start,threshold))
                self.assertEqual(union, sum(1 << j for j,a in enumerate(arrivals) if a <= threshold))

    def test_seeded_advantage_nonvacuous(self):
        costs = synthetic_costs(3)
        data = {"names": ("a","b","c"), "costs": costs}
        result = lab.persistence(data,{})
        self.assertGreater(result["strictly_helpful_singletons"],0)
        self.assertGreater(result["empty_full_threshold"], lab.minimax_arrivals(costs,1)[1])

    def test_zero_costs(self):
        costs = [[None if m >> j & 1 else 0 for j in range(3)] for m in range(8)]
        self.assertEqual(lab.minimax_arrivals(costs,0)[0], [0,0,0])
        self.assertEqual(lab.closure(costs,0,0),7)

    def test_random_monotone_costs(self):
        rng = random.Random(33)
        for _ in range(15):
            n = 4
            fresh = [[None if m >> j & 1 else rng.randrange(30) for j in range(n)] for m in range(1<<n)]
            for mask in range(1<<n):
                for j in range(n):
                    if mask >> j & 1: continue
                    for i in lab.bits(mask):
                        fresh[mask][j] = min(fresh[mask][j],fresh[mask^(1<<i)][j])
            for start in (0,1,2,4,8):
                arrival,_,_ = lab.minimax_arrivals(fresh,start)
                for t in range(30):
                    self.assertEqual(lab.closure(fresh,start,t), sum(1<<j for j,a in enumerate(arrival) if a<=t))


class GeometryTests(unittest.TestCase):
    def test_exact_ldl(self):
        self.assertEqual(lab.rational_spd([["2","1"],["1","2"]])[0][0],2)
        for matrix in ([[1,2],[2,1]], [[1,1],[1,1]], [[2,1],[0,2]]):
            with self.assertRaises(ValueError): lab.rational_spd(matrix)

    def test_projection_against_rational_schur(self):
        import numpy as np
        q = [[4,1,1],[1,3,0],[1,0,2]]
        result = lab.projection_scores(q,[[1,0,0]],[[0,1,0],[0,0,1],[0,1,1]])
        for v,norm in zip([[0,1,0],[0,0,1],[0,1,1]], result['residual_norm']):
            base = sum(lab.F(v[i]*q[i][j]*v[j]) for i in range(3) for j in range(3))
            cross = sum(v[i]*q[i][0] for i in range(3))
            exact = base - lab.F(cross*cross,q[0][0])
            self.assertAlmostEqual(norm,float(exact),places=10)

    def test_last_dimension_directional_tie(self):
        values = lab.projection_scores([[1,0,0],[0,1,0],[0,0,1]], [[1,0,0],[0,1,0]], [[0,0,1],[1,0,1],[-1,2,3]])
        self.assertEqual(lab.rank_interval(values['axis_fractional_unlock'],0,False), (1,3))

    def test_tie_intervals(self):
        self.assertEqual(lab.rank_interval([1,1,2],0,True), (1,2))
        self.assertEqual(lab.rank_interval([1,1,2],2,False), (1,1))
        self.assertEqual(lab.all_percentiles([1,1,1],True).tolist(),[.5,.5,.5])

    def test_projection_rejects_dependent(self):
        with self.assertRaises(ValueError): lab.projection_scores([[1,0],[0,1]],[[1,0]],[[2,0]])

    def test_vocab_independent_and_primitive(self):
        values = lab.fixed_vocabulary(4)
        self.assertEqual(len(values),len(set(values)))
        self.assertTrue(all(lab.primitive(v)==v for v in values))
        self.assertNotIn((19,7,0,0),values)


class BundleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.folder = fixture(Path(cls.tmp.name)/'source')
        cls.data = lab.validate_results(cls.folder)
    @classmethod
    def tearDownClass(cls): cls.tmp.cleanup()

    def test_source_shape(self):
        self.assertEqual(sum(len(r['events']) for r in self.data['runs']),180)
        self.assertEqual(CounterForTests(r['final_rank'] for r in self.data['runs']),{31:13,29:1})

    def test_mixed_vectors_create_strict_combinations(self):
        result = lab.filtration(self.data,{})
        first = result['runs'][0]
        self.assertEqual(first['endpoint'],{'quotient_dimension':12,'strict_dimension':8,'local_image_dimension':4})
        stages = first['states']
        # Mixed vectors e0+ej become strict after combining with the seed e0.
        self.assertTrue(any(s.get('strict_increment') == 1 and not s['acquired_vector_itself_strict'] for s in stages[1:]))
        for s in stages:
            self.assertEqual(s['quotient_dimension'],s['strict_dimension']+s['local_image_dimension'])
        final = next(r for r in result['subspace_comparisons'] if r['quotient_dimension']==14)
        self.assertTrue(final['terminal_convergence_is_tautological'])

    def test_tampered_output_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            import shutil
            p = Path(td)/'bundle'; shutil.copytree(self.folder,p)
            obj = lab.read(p/'trajectories.json'); obj['runs'][0]['final_rank'] += 1
            lab.atomic(p/'trajectories.json',obj)
            with self.assertRaisesRegex(ValueError,'hash mismatch'): lab.validate_results(p)

    def test_missing_bundle_not_synthetic(self):
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaisesRegex(ValueError,'No complete'): lab.find_results(Path(td),Path(td))

    def test_duplicate_keys_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)/'bad.json'; p.write_text('{"x":1,"x":2}')
            with self.assertRaises(ValueError): lab.read(p)

    def test_nan_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'bad.json'; p.write_text('{"x":NaN}')
            with self.assertRaises(ValueError): lab.read(p)

    def test_failclosed_nonintegral(self):
        with tempfile.TemporaryDirectory() as td:
            import shutil
            p=Path(td)/'bundle'; shutil.copytree(self.folder,p)
            obj=lab.read(p/'trajectories.json'); obj['runs'][0]['stages'][0]['new'][0]['denominator']=2
            lab.atomic(p/'trajectories.json',obj)
            rep=lab.read(p/'REPORT.json'); rep['outputs']['trajectories']=lab.digest(p/'trajectories.json');lab.atomic(p/'REPORT.json',rep)
            with self.assertRaisesRegex(ValueError,'nonintegral'):lab.validate_results(p)

    def test_production_size_vocabulary(self):
        self.assertEqual(len(lab.fixed_vocabulary(14)),18760)

    def test_lock_conflict(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)
            with lab.lock(p):
                with self.assertRaisesRegex(RuntimeError,'another controller'):
                    with lab.lock(p): pass

    def test_timeout_preserved(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td); phase=p/'phase'; phase.mkdir()
            row=lab.supervise([sys.executable,'-c','import time;time.sleep(5)'],p,phase,.15,1024**3)
            self.assertEqual(row['outcome'],'TIMEOUT_CENSORED')
            self.assertTrue((phase/'receipt.json').is_file())

    def test_no_output_is_not_sealed(self):
        with tempfile.TemporaryDirectory() as td:
            self.assertFalse(lab.sealed_phase(Path(td),'persistence'))


def CounterForTests(values):
    from collections import Counter
    return dict(Counter(values))


class PredictionTests(unittest.TestCase):
    def test_oov_is_not_injected_and_all_events_count(self):
        # Small standalone geometry probe avoids pretending this is actual EC data.
        data={'names':('a','b','c'), 'form':[[3,1,0],[1,2,0],[0,0,1]], 'runs':[
          {'seed':'a','seed_index':0,'events':[{'word':(19,7,1),'primitive':(19,7,1),'epoch':0}]},
          {'seed':'b','seed_index':1,'events':[{'word':(1,0,1),'primitive':(1,0,1),'epoch':0}]}]}
        policy={'profile':'smoke','random_draws':20,'random_seed':8}
        out=lab.next_moves(data,policy)
        for arm in out['arms'].values():
            self.assertEqual(arm['summary']['events'],2)
            row=arm['events'][0]
            self.assertFalse(row['covered'])
            self.assertEqual(row['scores']['residual_norm']['percentile_loss'],1)
            self.assertFalse(row['scores']['residual_norm']['conservative_top5'])
        self.assertEqual(out,lab.next_moves(data,policy))


class SafetyTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.source = fixture(self.root/'source')
        self.folder = self.root/'run'
        self.policy = {'phase_seconds':60,'memory_mib':1024,'random_draws':20,'random_seed':3020911,'profile':'smoke'}
    def tearDown(self): self.tmp.cleanup()
    def prepared(self):
        lab.prepare(self.folder,self.root,self.source,self.source/'landscape.json',self.policy)

    def test_prepare_and_guard(self):
        self.prepared()
        self.assertEqual(lab.guard(self.folder)['source_acquisitions'],180)
        self.assertEqual(lab.guard(self.folder)['dataset_kind'],'SYNTHETIC_TEST_FIXTURE')

    def test_frozen_inputs_isolate_external_changes(self):
        self.prepared()
        (self.source/'trajectories.json').write_text('{}')
        lab.guard(self.folder)

    def test_changed_snapshot_rejected(self):
        self.prepared()
        (self.folder/'inputs/trajectories.json').write_text('{}')
        with self.assertRaisesRegex(ValueError,'frozen input changed'):lab.guard(self.folder)

    def test_changed_plan_rejected(self):
        self.prepared()
        plan=lab.read(self.folder/'plan.json'); plan['policy']['random_seed']+=1
        lab.atomic(self.folder/'plan.json',plan)
        with self.assertRaisesRegex(ValueError,'plan changed'):lab.guard(self.folder)

    def test_changed_source_rejected(self):
        self.prepared()
        with (self.folder/'source/runner.py').open('a') as f:f.write('\n# mutation\n')
        with self.assertRaisesRegex(ValueError,'frozen source changed'):lab.guard(self.folder)

    def test_interrupted_phase_never_restarted(self):
        self.prepared()
        (self.folder/'phases/persistence').mkdir(parents=True)
        with mock.patch.object(lab,'supervise') as launch:
            with self.assertRaisesRegex(ValueError,'unreceipted'):lab.execute(self.folder)
            launch.assert_not_called()

    def test_stop_prevents_launch(self):
        self.prepared()
        (self.folder/'STOP').touch()
        with mock.patch.object(lab,'supervise') as launch:
            with self.assertRaisesRegex(ValueError,'STOP'):lab.execute(self.folder)
            launch.assert_not_called()

    def test_existing_folder_not_overwritten(self):
        self.folder.mkdir()
        with self.assertRaisesRegex(ValueError,'exists'):
            lab.prepare(self.folder,self.root,self.source,self.source/'landscape.json',self.policy)

    def test_failed_stage_not_sealed(self):
        self.prepared()
        with mock.patch.object(lab,'supervise',return_value={'outcome':'FAILED','returncode':2}):
            with self.assertRaisesRegex(RuntimeError,'FAILED'):lab.execute(self.folder)
        self.assertFalse((self.folder/'phases/persistence/seal.json').exists())
        self.assertEqual(lab.read(self.folder/'REPORT.json')['status'],'UNKNOWN_FAILED_OR_CENSORED')


if __name__=='__main__': unittest.main()
