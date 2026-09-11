"""Regressions and an end-to-end *synthetic*, not Curve302-data, smoke."""
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

CAS = Path(__file__).resolve().parents[1]/'cas'
sys.path.insert(0, str(CAS))
import curve302_short_core_controls as c
import run_curve302_short_core_controls as r
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form


def fixture(parent):
    """14 runs, 180 actual acquisitions, 194 prefixes, 196 complete directions."""
    parent = Path(parent)
    short, structure = parent/'short', parent/'structure'
    short.mkdir(); structure.mkdir()
    n = 14
    axes = [tuple(int(i == j) for i in range(n)) for j in range(n)]
    runs, normalized = [], []
    for seed in range(n):
        order = [j for j in range(n) if j != seed]
        if seed == 0:
            order = order[:11]
        events, stages = [], []
        for step, j in enumerate(order):
            word = list(axes[j]); word[seed] = (-1 if step % 3 == 0 else 1) if step % 2 else 0
            word = tuple(word)
            p = c.primitive(word)
            events.append({'word': word, 'primitive': p, 'epoch': step})
            stages.append({'epoch': step, 'new':[{'integral': True, 'denominator':1,
                                                'quotient_word': word, 'primitive_quotient_word':p}]})
        runs.append({'seed':r.NAMES[seed], 'final_rank':18+len(order), 'stages':stages})
        normalized.append({'seed':r.NAMES[seed], 'seed_index':seed, 'events':events, 'final_dimension':len(order)+1, 'order':order})
    relation = {'status':'PASS_QUOTIENT_RELATION_ANALYSIS', 'direction_ids':r.NAMES,
                'schur_quotient':[['1' if i == j else '0' for j in range(n)] for i in range(n)]}
    r.atomic(structure/'quotient-relations.json', relation)
    r.atomic(structure/'trajectories.json', {'status':'PASS_ALL_14_SEEDED_TRAJECTORIES_RECONCILED', 'direction_ids':r.NAMES, 'runs':runs})
    r.atomic(structure/'REPORT.json', {'status':'PASS_THREE_CLOSURE_EXPERIMENTS', 'outputs':{
        s:r.sha(structure/f'{s}.json') for s in ('quotient-relations','trajectories')}})
    vectors = [(F(1),v) for v in axes]
    for i,j in itertools.combinations(range(n),2):
        for sign in (-1,1):
            v = list(axes[i]); v[j]=sign
            vectors.append((F(2),tuple(v)))
    vectors.sort()
    with (short/'primitive-directions.tsv').open('w') as out:
        out.write('norm\tvector\n')
        for norm,v in vectors:
            out.write(f'{norm}\t'+','.join(map(str,v))+'\n')
    r.atomic(short/'enumeration.json', {'schema':'curve302-short-vector-enumeration.v1',
                                      'status':'PASS_COMPLETE_EXACT_ENUMERATION','bound':'2','direction_count':len(vectors),
                                      'enumeration_sha256':r.sha(short/'primitive-directions.tsv')})
    cores=[]
    for d in range(1,15):
        participating=[x for x in normalized if x['final_dimension']>=d]
        common=set(range(n))
        for x in participating:
            common &= {x['seed_index'], *x['order'][:d-1]}
        basis=[axes[j] for j in sorted(common)]
        cores.append({'quotient_dimension':d,'run_count':len(participating),'rank':len(basis),'basis':basis})
    r.atomic(short/'filtration.json', {'status':'PASS_INTRINSIC_AND_OBSERVED_FILTRATIONS','observed_common_integral_cores':cores})
    basins=[]
    for core in cores:
        d=core['quotient_dimension']
        if d < 2 or not core['rank']:
            continue
        for run in normalized:
            if run['final_dimension']<d:
                continue
            # Independent fixture oracle: these triangular acquisitions have
            # exactly coordinate-axis spans, so set algebra gives basin counts.
            prefix={run['seed_index'], *run['order'][:d-2]}
            targets={next(i for i,x in enumerate(v) if x) for v in core['basis']}
            missing=targets-prefix
            deficit=len(missing)
            actual_norm=sum(x*x for x in run['events'][d-2]['primitive'])
            eligible=[(norm,v) for norm,v in vectors if any(x and i not in prefix for i,x in enumerate(v))]
            def hit(v):
                if deficit==0: return True
                return all(not x or i in prefix or i in missing for i,x in enumerate(v))
            bound=[(norm,v) for norm,v in eligible if norm<=actual_norm]
            basins.append({'landmark_dimension':d,'seed':run['seed'],'core_rank':len(core['basis']),
                           'deficit_before':deficit,'actual_next_norm':str(actual_norm),
                           'eligible_at_actual_norm':len(bound),'hits_at_actual_norm':sum(hit(v) for _,v in bound),
                           'eligible_at_complete_bound':len(eligible),'hits_at_complete_bound':sum(hit(v) for _,v in eligible)})
    r.atomic(short/'ranks-basins.json', {'status':'PASS_COMPLETE_RANK_AND_BASIN_CENSUS', 'basins':basins})
    r.atomic(short/'plan.json', {'source':str(structure), 'source_hashes':{name:r.sha(structure/name) for name in r.STRUCTURE_FILES}})
    r.atomic(short/'REPORT.json', {'status':'PASS_THREE_SHORT_VECTOR_CORE_EXPERIMENTS',
                                   'outputs':{s:r.sha(short/f'{s}.json') for s in ('enumeration','filtration','ranks-basins')},
                                   'enumeration_sha256':r.sha(short/'primitive-directions.tsv')})
    return short,structure


def small_vocab(parent, rows, wanted=(), n=2, bound=None):
    path=Path(parent)/'vectors.tsv'
    rows=sorted(rows)
    path.write_text('norm\tvector\n'+''.join(str(a)+'\t'+','.join(map(str,v))+'\n' for a,v in rows))
    return c.Vocabulary.load(path,n,len(rows),bound or rows[-1][0],wanted)


class ExactMathTests(unittest.TestCase):
    def test_reject_silent_integer_truncation(self):
        for value in (True,False,1.0,1.2,'3/2'):
            with self.assertRaises(c.InvalidEvidence): c.integer(value)
        self.assertEqual(c.integer('10000000000000000000000000'),10**25)

    def test_primitive_sign(self):
        self.assertEqual(c.primitive((0,-6,9)),(0,2,-3))

    def test_nontrivial_index_despite_primitive_generators(self):
        out=c.smith_packet(((1,1),(1,-1)),2)
        self.assertEqual(out['index_in_ambient'],'2')
        self.assertFalse(out['generates_ambient_integrally'])
        self.assertEqual([w['smallest_positive_multiple'] for w in out['axis_multiple_witnesses']],['2','2'])

    def test_actual_integral_generation_witnesses(self):
        rows=((2,1),(3,1),(1,0))
        out=c.smith_packet(rows,2)
        self.assertEqual(out['index_in_ambient'],'1')
        for w in out['axis_multiple_witnesses']:
            self.assertEqual([sum(a*b[j] for a,b in zip(w['coefficients'],rows)) for j in range(2)],
                             [int(j==w['axis']) for j in range(2)])

    def test_partial_rank_index_not_ambient(self):
        out=c.smith_packet(((2,0,0),(0,3,0)),3)
        self.assertEqual(out['index_in_saturation'],'6')
        self.assertIsNone(out['index_in_ambient'])

    def test_primitive_quotient_rule_rejects_hidden_halving(self):
        state=c.QuotientMap(2).extend((1,1))
        self.assertEqual(state.admission((1,-1))[0],'NONPRIMITIVE')
        with self.assertRaisesRegex(c.InvalidEvidence,'NONPRIMITIVE'): state.extend((1,-1))
        self.assertEqual(state.admission((0,1))[0],'PRIMITIVE')

    def test_common_basis_order_and_unimodular_presentation_do_not_matter(self):
        plane=c.QuotientMap(3).extend((1,0,0)).extend((0,1,0))
        c.verify_common_core([plane,plane],((0,1,0),(1,1,0)),3)

    def test_no_parity_only_containment(self):
        state=c.QuotientMap(2).extend((1,2))
        self.assertFalse(state.contains(((1,0),)))

    def test_primitive_rule_matches_independent_smith_exhaustively(self):
        prefix=((1,2,1),(0,1,1))
        state=c.QuotientMap(3)
        for v in prefix: state=state.extend(v)
        for v in itertools.product(range(-2,3),repeat=3):
            mat=Matrix([*prefix,v])
            if mat.rank()==2: expected='DEPENDENT'
            else:
                d=smith_normal_form(mat,domain=ZZ)
                expected='PRIMITIVE' if abs(math.prod(int(d[i,i]) for i in range(3)))==1 else 'NONPRIMITIVE'
            self.assertEqual(state.admission(v)[0],expected)

    def test_quotient_map_updates_preserve_prefix_exactly(self):
        rows=((1,1,0,1),(0,1,1,0),(1,0,0,0))
        state=c.QuotientMap(4)
        for v in rows:
            state=state.extend(v)
        self.assertTrue(state.contains(rows))
        self.assertEqual(state.rank,3)
        for v in itertools.product(range(-1,2),repeat=4):
            self.assertEqual(state.contains((v,)),Matrix([*rows,v]).rank()==3)

    def test_rational_rank_int_matches_sympy(self):
        rng=random.Random(47)
        for _ in range(20):
            rows=[tuple(rng.randrange(-5,6) for _ in range(5)) for _ in range(rng.randrange(1,8))]
            self.assertEqual(c.rational_rank(rows,5),Matrix(rows).rank())

    def test_common_core_verifies_integral_not_just_span(self):
        states=[c.QuotientMap(3).extend((1,0,0)).extend(v) for v in ((0,1,0),(0,0,1))]
        c.verify_common_core(states,((1,0,0),),3)
        with self.assertRaisesRegex(c.InvalidEvidence,'saturated'):
            c.verify_common_core(states,((2,0,0),),3)


class VocabularyTests(unittest.TestCase):
    def test_ties_never_split_bands(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows=[(F(1),(1,i)) for i in range(14)]
            v=small_vocab(tmp,rows)
            self.assertEqual(v.band(0)['index_start'],0)
            self.assertEqual(v.band(0)['index_stop'],14)
            self.assertTrue(v.band(0)['tie_expanded'])

    def test_rank_decade_boundary(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows=[(F(i+1),(1,i)) for i in range(120)]
            v=small_vocab(tmp,rows)
            self.assertEqual((v.band(9)['index_start'],v.band(9)['index_stop']),(0,10))
            self.assertEqual((v.band(10)['index_start'],v.band(10)['index_stop']),(10,100))
            self.assertEqual((v.band(100)['index_start'],v.band(100)['index_stop']),(100,120))

    def test_29th_shell_separate(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows=[(F(1),(1,i)) for i in range(35)]
            v=small_vocab(tmp,rows)
            out=c.first_direction_index(v)
            self.assertEqual(out['literal_first_29']['rows'],29)
            self.assertEqual(out['whole_29th_norm_shell']['rows'],35)
            self.assertTrue(out['tie_boundary'])

    def test_corrupt_norm_witness_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'v.tsv';p.write_text('norm\tvector\n2\t1,0\n')
            with self.assertRaisesRegex(c.InvalidEvidence,'norm disagrees'):
                c.Vocabulary.load(p,2,1,F(2),(),((F(1),F(0)),(F(0),F(1))))

    def test_duplicate_enumeration_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(c.InvalidEvidence,'duplicated'):
                small_vocab(tmp,[(F(1),(1,0)),(F(1),(1,0))])

    def test_storage_overflow_is_not_silent(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(c.InvalidEvidence,'int64'):
                small_vocab(tmp,[(F(1),(1,2**70))])

    def test_missing_acquisition_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(c.InvalidEvidence,'missing'):
                small_vocab(tmp,[(F(1),(1,0))],wanted=((0,1),))


class BasinAndSamplingTests(unittest.TestCase):
    def basin(self,seed,deficit,num,den):
        return {'landmark_dimension':6,'seed':seed,'deficit_before':deficit,
                'eligible_at_actual_norm':den,'hits_at_actual_norm':num,
                'eligible_at_complete_bound':den*2,'hits_at_complete_bound':num*2}

    def test_exact_integer_text_basin_fields_normalize_without_truncation(self):
        row=self.basin('a',1,1,4)
        row={key:str(value) if isinstance(value,int) else value for key,value in row.items()}
        self.assertEqual(c.corrected_basins([row])['overall']['nontrivial_cases'],1)

    def test_automatic_basins_excluded_from_both_denominators(self):
        out=c.corrected_basins([self.basin('a',0,100,100),self.basin('b',1,1,4),self.basin('c',1,2,4)])
        stats=out['overall'];self.assertEqual(stats['nontrivial_cases'],2)
        self.assertEqual(stats['bounds']['actual_norm']['pooled_fraction'],'3/8')
        self.assertEqual(stats['bounds']['actual_norm']['hits'],3)
        self.assertEqual(stats['bounds']['actual_norm']['eligible'],8)

    def test_all_automatic_basins_report_no_estimate(self):
        out=c.corrected_basins([self.basin('a',0,10,10)])
        self.assertIsNone(out['overall']['bounds']['actual_norm']['pooled_fraction'])

    def test_pooled_vs_equal_case_mean(self):
        out=c.corrected_basins([self.basin('a',1,1,2),self.basin('b',1,1,10)])
        b=out['overall']['bounds']['actual_norm']
        self.assertEqual(b['pooled_fraction'],'1/6');self.assertEqual(b['mean_case_fraction'],'3/10')

    def test_duplicate_basin_rejected(self):
        with self.assertRaisesRegex(c.InvalidEvidence,'duplicate'):
            c.corrected_basins([self.basin('a',0,1,1)]*2)

    def test_reference_does_not_drop_censored_panels(self):
        out=c.boolean_reference([True,False,None,None])
        self.assertEqual(out['empirical_fraction_lower'],'1/4')
        self.assertEqual(out['empirical_fraction_upper'],'3/4')
        self.assertEqual(out['panels'],4)

    def test_no_admissible_band_is_censored_not_widened(self):
        with tempfile.TemporaryDirectory() as tmp:
            vocab=small_vocab(tmp,[(F(1),(1,0)),(F(2),(0,1))])
            state=c.QuotientMap(2).extend((1,0))
            out=c.sample_extension(vocab,state,{'index_start':0,'index_stop':1},random.Random(2),max_draws=4,rescue_after=2)
            self.assertEqual(out['status'],'CENSORED_EMPTY_ADMISSIBLE_BAND')
            self.assertIsNone(out['index'])

    def test_large_unresolved_band_uses_unknown_not_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            vocab=small_vocab(tmp,[(F(1),(1,0))])
            state=c.QuotientMap(2).extend((1,0))
            out=c.sample_extension(vocab,state,{'index_start':0,'index_stop':1},random.Random(2),max_draws=4,rescue_after=2,rescue_scan_limit=0)
            self.assertEqual(out['status'],'CENSORED_DRAW_CAP')

    def test_only_primitive_extensions_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            vocab=small_vocab(tmp,[(F(1),(1,1)),(F(2),(1,-1)),(F(3),(0,1))])
            state=c.QuotientMap(2).extend((1,1))
            for seed in range(20):
                out=c.sample_extension(vocab,state,{'index_start':0,'index_stop':3},random.Random(seed))
                self.assertEqual(vocab.vector(out['index']),(0,1))

    def test_random_streams_deterministic_and_run_local(self):
        a=c.seeded_rng('x',0,'a');b=c.seeded_rng('x',0,'a');d=c.seeded_rng('x',1,'a')
        self.assertEqual([a.randrange(1000) for _ in range(10)],[b.randrange(1000) for _ in range(10)])
        self.assertNotEqual(c.seeded_rng('x',0,'a').getstate(),d.getstate())

    def test_core_evaluation_does_not_affect_sampling(self):
        with tempfile.TemporaryDirectory() as tmp:
            v=small_vocab(tmp,[(F(1),(1,0)),(F(1),(0,1)),(F(2),(1,1))])
            runs=[{'seed':'a','seed_index':0,'final_dimension':2,'events':[{'band':{'index_start':0,'index_stop':3}}]}]
            policy={'max_draws':100,'rescue_after':5,'rescue_scan_limit':10}
            cores=[{'quotient_dimension':1,'run_count':1,'basis':((1,0),)}]
            a=c.simulate_panel(v,runs,cores,2,(), 's',0,policy)
            b=c.simulate_panel(v,runs,[],2,(), 's',0,policy)
            self.assertEqual(a['runs'],b['runs'])

    def test_late_cohort_is_not_fourteen(self):
        runs=[{'seed':'a','final_dimension':2},{'seed':'b','final_dimension':1}]
        a=c.QuotientMap(2).extend((1,0));full=a.extend((0,1))
        maps={'a':{1:a,2:full},'b':{1:a}}
        metrics=c.panel_metrics(maps,runs,[{'quotient_dimension':2,'basis':((1,0),(0,1))}],2,())
        self.assertEqual(metrics['cores'][0]['participants'],1)
        self.assertTrue(metrics['cores'][0]['trivial_target'])


class ControllerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory(prefix='curve302-control-tests-')
        cls.base=Path(cls.temp.name)
        cls.short,cls.structure=fixture(cls.base)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_source_fixture_complete(self):
        data=r.validate_bundle(self.short,self.structure)
        maps,count=c.track_observed(data['runs'],14,data['cores'])
        self.assertEqual(count,194)
        self.assertEqual(sum(len(x['events']) for x in data['runs']),180)

    def test_source_hash_corruption_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            import shutil
            dest=Path(temp)/'short';shutil.copytree(self.short,dest)
            with (dest/'filtration.json').open('a') as out:out.write(' ')
            with self.assertRaisesRegex(c.InvalidEvidence,'hash mismatch'):
                r.validate_bundle(dest,self.structure)

    def test_duplicate_seed_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            short,structure=fixture(temp)
            traj=r.read(structure/'trajectories.json');traj['runs'][1]['seed']=traj['runs'][0]['seed']
            r.atomic(structure/'trajectories.json',traj)
            report=r.read(structure/'REPORT.json');report['outputs']['trajectories']=r.sha(structure/'trajectories.json');r.atomic(structure/'REPORT.json',report)
            plan=r.read(short/'plan.json');plan['source_hashes']={p:r.sha(structure/p) for p in r.STRUCTURE_FILES};r.atomic(short/'plan.json',plan)
            with self.assertRaisesRegex(c.InvalidEvidence,'distinct seeds'):
                r.validate_bundle(short,structure)

    def test_lock_rejects_second_controller(self):
        with tempfile.TemporaryDirectory() as temp:
            with r.locked(temp):
                with self.assertRaisesRegex(c.InvalidEvidence,'another controller'):
                    with r.locked(temp): pass

    def test_timeout_supervisor_reports_unknown(self):
        with tempfile.TemporaryDirectory() as temp:
            receipt=r.supervise([sys.executable,'-c','import time; time.sleep(30)'],Path(temp)/'log',.1,1024**3)
            self.assertEqual(receipt['outcome'],'UNKNOWN_TIMEOUT')

    def test_missing_output_folder_is_not_prepared(self):
        with self.assertRaises(FileNotFoundError):r.guard(self.base/'not-there')

    def test_snapshot_corruption_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp)/'out';r.prepare(out,self.short,panels=1)
            with (out/'inputs/short/filtration.json').open('a') as f:f.write(' ')
            with self.assertRaisesRegex(c.InvalidEvidence,'snapshot changed'):r.guard(out)

    def test_interrupted_stage_not_reused(self):
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp)/'out';r.prepare(out,self.short,panels=1)
            r.atomic(out/'phases/index/started.json',{'stage':'index'})
            with self.assertRaisesRegex(c.InvalidEvidence,'unsealed'):
                r.run_stage(out,'index',r.guard(out))

    def test_full_synthetic_run_and_deterministic_check(self):
        out=self.base/'integration-output'
        cmd=[sys.executable,str(r.SELF)]
        args=['run','--source',str(self.short),'--folder',str(out),'--panels','2','--stage-seconds','120']
        executed=subprocess.run(cmd+args,capture_output=True,text=True,timeout=160)
        if executed.returncode:
            logs='\n'.join(p.read_text() for p in out.glob('phases/*/worker.log'))
            self.fail(executed.stdout+executed.stderr+logs)
        checked=subprocess.run(cmd+['check','--folder',str(out)],capture_output=True,text=True,timeout=160)
        self.assertEqual(checked.returncode,0,checked.stdout+checked.stderr)
        report=r.read(out/'REPORT.json');self.assertTrue(report['status'].startswith('PASS_THREE_STATIC_CORE_CONTROLS'))
        self.assertEqual(r.read(out/'index.json')['literal_first_29']['index_in_ambient'],'1')
        summary=r.read(out/'sampling.json')['summary']
        self.assertEqual(summary['panels'],2)
        self.assertEqual(len(r.read(out/'bands.json')['runs']),14)
        hashes={p:r.sha(out/p) for p in ('index.json','basins.json','sampling.json','SUMMARY.md')}
        resumed=subprocess.run(cmd+['resume','--folder',str(out)],capture_output=True,text=True,timeout=30)
        self.assertEqual(resumed.returncode,0,resumed.stdout+resumed.stderr)
        self.assertEqual(hashes,{p:r.sha(out/p) for p in hashes})
        # Original sources can move after prepare: the private snapshots suffice.
        self.assertTrue((out/'code'/r.SELF.name).exists())


if __name__=='__main__':unittest.main()
