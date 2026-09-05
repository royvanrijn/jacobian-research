#!/usr/bin/env sage-python
"""Finite compact-R17 selection and point search, with immutable per-chart witnesses.

This is a new, explicitly family-informed experiment, not an amendment to any
historical frozen cohort. Scores and numerical heights only schedule work.
"""
import argparse
from dataclasses import asdict
from fractions import Fraction as F
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from math import gcd, log
from pathlib import Path
import sys
import time

from sage.all import EllipticCurve, QQ, prime_range
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT/'elliptic-curves/cas'
sys.path[:0] = [str(CAS), str(ROOT/'elliptic-curves')]
sys.set_int_max_str_digits(0)
from ecsearch.q12o5867_specialization import load_q12o5867_data, evaluate_projective_specialization
from research_runtime.store import checkpoint
from research_runtime.search_state import raw_state, reduction_cache
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from pointed_quartic_search import PointedQuarticSearch, sources

MODEL = ROOT/'elkies-k3/data/fibrations/elkies_2026_published_r17_model.json'
SECTIONS = ROOT/'elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json'
HOLES = ROOT/'artifacts/generated-results/elliptic-curves/half_lattice_fake_descent_r17_matrix_blind_v1.json'
DIRECTORY = ROOT/'artifacts/local/elliptic-curves/compact-r17-prospective-v1'


def read(path): return json.loads(path.read_text())
def hashed(path): return sha256(path.read_bytes()).hexdigest()
def identity(data): return sha256(json.dumps(data, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def provenance():
    return {**sources(), **{str(p.relative_to(ROOT)): hashed(p) for p in
        (Path(__file__), MODEL, SECTIONS, HOLES, CAS/'half_lattice_fake_descent_replay.sage',
         ROOT/'elliptic-curves/ecsearch/q12o5867_specialization.py',CAS/'research_runtime/memory_store.py')}}


def trace_table(model, p):
    """Exact character sums on all p+1 residues; singular rows are censored.

    p <= 4093 here, so all int64 products have ample headroom. The infinity
    model uses the homogeneous weights 8 and 12 of this K3 family.
    """
    if not 3 < p <= 4093: raise ValueError('table prime outside declared safe range')
    t = np.arange(p, dtype=np.int64)
    def evaluate(key):
        cs = [int(F(c).numerator * pow(F(c).denominator, -1, p) % p) for c in model[key]]
        v = np.zeros(p, dtype=np.int64)
        for c in reversed(cs): v = (v*t+c)%p
        return np.append(v, cs[-1])
    aa, bb = evaluate('A_coefficients_low_to_high'), evaluate('B_coefficients_low_to_high')
    chi = np.full(p, -1, dtype=np.int64); chi[(t*t)%p] = 1; chi[0] = 0
    cube = (t*t%p)*t%p
    traces = np.array([-int(chi[(cube+a*t+b)%p].sum()) for a,b in zip(aa,bb)], dtype=np.int64)
    good = (4*((aa*aa%p)*aa%p)+27*(bb*bb%p))%p != 0
    return traces, good


def add_scores(rows, model, primes, directory, label):
    n = np.array([r['numerator'] for r in rows],dtype=np.int64)
    d = np.array([r['denominator'] for r in rows],dtype=np.int64)
    scores = np.zeros(len(rows)); good_counts = np.zeros(len(rows),dtype=np.int64)
    for k,p in enumerate(primes):
        p=int(p); path=directory/'trace-tables'/f'{p}.json'
        key={'model_sha256':hashed(MODEL),'prime':p,'algorithm':'homogeneous-character-sum-v1'}
        if path.exists():
            table=read(path)
            if table['input']!=key: raise ArithmeticError('trace-table binding changed')
            traces=np.array(table['traces']); good=np.array(table['good'],dtype=bool)
        else:
            traces,good=trace_table(model,p)
            checkpoint(path,{'input':key,'traces':traces.tolist(),'good':good.tolist()})
        inv=np.array([0]+[pow(x,-1,p) for x in range(1,p)],dtype=np.int64)
        ix=(n%p)*inv[d%p]%p; ix[d%p==0]=p
        # The classical finite Nagao/Mestre score; no rank inference.
        scores+=np.where(good[ix],(2-traces[ix])/(p+1-traces[ix])*log(p),0)
        good_counts+=good[ix]
        if k%50==0: print(f'SCORE {label} prime={p} rows={len(rows)}',flush=True)
    for i,r in enumerate(rows): r[label]={'score':float(scores[i]),'good_primes':int(good_counts[i])}


def prepare(directory):
    path=directory/'population.json'
    if path.exists(): raise FileExistsError('population is immutable')
    protocol={'schema':'elliptic-curves.compact-r17-prospective-protocol.v1',
        'height':256,'first_prime_bound':997,'first_survivors':128,'second_prime_bound':4093,
        'finalists':16,'score':'sum_good (2-a_p)*log(p)/(p+1-a_p)',
        'control_parameters':['3/8','-2/377'],
        'charts_per_initial_fibre':43,'chart_height':100000,'seconds_per_chart':4,
        'metric':'specialized numerical canonical height Gram rounded at 1e6; scheduling only',
        'centres':'retained 43 generic deepest parities; specialized CVP representatives',
        'coordinate_policy':{'kind':'metric','weight':'16'},
        'novelty':'post-selection exact j and rational isomorphism comparison; no universal novelty claim',
        'certificate':'exact points, no-rational-2-torsion witness, independent finite quotient columns',
        'sources':provenance(), 'claim_boundary':'All misses are bounded. Rank >=28 is a near-record target; >=32 is the record target.'}
    protocol_path=directory/'protocol.json'
    if protocol_path.exists() and read(protocol_path)!=protocol: raise ArithmeticError('protocol changed')
    checkpoint(protocol_path,protocol)
    h=protocol['height']
    rows=[{'numerator':n,'denominator':d,'parameter':str(F(n,d))}
          for d in range(1,h+1) for n in range(-h,h+1) if gcd(n,d)==1]
    model=read(MODEL)
    add_scores(rows,model,list(prime_range(5,998)),directory,'prefix')
    key=lambda r:(-r['prefix']['score'],max(abs(r['numerator']),r['denominator']),r['numerator'],r['denominator'])
    rows.sort(key=key)
    checkpoint(directory/'prefix-population.json',{'protocol_hash':identity(protocol),'rows':rows})
    selected=rows[:128]
    add_scores(selected,model,list(prime_range(998,4094)),directory,'extension')
    for r in selected: r['combined_score']=r['prefix']['score']+r['extension']['score']
    selected.sort(key=lambda r:(-r['combined_score'],max(abs(r['numerator']),r['denominator']),r['numerator'],r['denominator']))
    checkpoint(path,{'protocol_hash':identity(protocol),'candidate_count':len(rows),
        'retained_candidates':selected,'finalists':selected[:16],
        'public_points_or_record_equations_used_for_selection':False})
    print('FROZEN',[(r['parameter'],round(r['combined_score'],3)) for r in selected[:16]],flush=True)


def run_fibre(directory, label, parameter):
    protocol=read(directory/'protocol.json')
    if protocol['sources']!=provenance():
        amendment=read(directory/'implementation-amendment.json')
        if amendment['original_protocol_hash']!=identity(protocol) or amendment['sources']!=provenance():
            raise ArithmeticError('unrecorded implementation change')
    folder=directory/label; output=folder/'result.json'
    if output.exists() and read(output).get('status')=='COMPLETE':
        print('RETAINED',label,read(output)['rank_lower_bound'],flush=True);return
    family=load_q12o5867_data(MODEL,SECTIONS); t=F(parameter)
    spec=evaluate_projective_specialization(family,t.numerator,t.denominator)
    # An integral short model is enough: no conductor factorization gate.
    curve=EllipticCurve(QQ,list(spec.model)); points=[curve(list(p)) for p in spec.points]
    minimal2=curve.local_data(2).minimal_model()
    transform=curve.isomorphism_to(minimal2)
    a1,a2,a3,a4,a6=minimal2.a_invariants()
    c2,c4,c6=a2+a1*a1/4,a4+a1*a3/2,a6+a3*a3/4
    model=tuple(F(str(v)) for v in (0,0,0,c4-c2*c2/3,c6-c2*c4/3+2*c2**3/27))
    generic=[]
    for P in points:
        x,y=transform(P).xy(); generic.append((F(str(x+c2/3)),F(str(y+(a1*x+a3)/2))))
    cache=ReductionCache(MemoryFactStore())
    initial=raw_state(model,generic,cache=cache,prime_bound=1000)
    if initial.rank!=17: raise ArithmeticError('generic subgroup not certified rank 17')
    engine=SourceFileLoader('compact_r17_height_engine',str(CAS/'half_lattice_fake_descent_replay.sage')).load_module()
    gram=engine.canonical_height_gram(model,generic)
    rounded=[[int((c*1000000).to_integral_value()) for c in row] for row in gram]
    oracle=engine.CosetOracle(rounded)
    masks=sorted(c['mask'] for c in read(HOLES)['fibres'][0]['cover_records'])
    if len(masks)!=43 or len(set(masks))!=43: raise ArithmeticError('generic parity list changed')
    centres=[]
    for mask in masks:
        norm,rep,error=oracle.solve(mask); centres.append({'mask':mask,'representative':list(rep),'metric_norm':norm})
    centres.sort(key=lambda r:(-r['metric_norm'],r['mask']))
    result={'schema':'elliptic-curves.compact-r17-prospective-fibre.v1','label':label,'parameter':parameter,
        'protocol_hash':identity(protocol),'curve':list(map(str,model)),
        'generic_points':[list(map(str,p)) for p in generic], 'initial_state':initial.record(),
        'metric_gram':[[str(c) for c in row] for row in gram], 'centres':centres,
        'charts':[], 'status':'RUNNING','rank_lower_bound':17,'sources':provenance(),
        'claim_boundary':'Mod-2 dependence is inconclusive; ambiguous hits remain available. Numerical heights are not rank proofs.'}
    state=initial
    if output.exists():
        from research_runtime.mw_state import MWState
        from research_runtime.supervisor import preserve_previous
        saved=read(output)
        if any(saved[k]!=result[k] for k in ('curve','generic_points','centres','parameter','protocol_hash','initial_state')):
            raise ArithmeticError('partial worker input or centre plan changed')
        if 'final_state' in saved:
            cache.store.import_snapshot(saved['arithmetic_facts'])
            state=MWState.from_record(saved['final_state'],cache=cache)
            result=saved
            result.setdefault('implementation_history',[]).append(provenance())
        preserve_previous(output)
    checkpoint(output,result)
    # This is a proof-search allowance, never a dependence test. Retain every
    # ambiguous point. Large point lists otherwise create millions of unused
    # per-point/per-prime cache files while proving no additional directions.
    admission_prime_bound=251
    primes=tuple(map(int,prime_range(3,admission_prime_bound+1)))
    for i in range(len(result['charts']),len(centres)):
        centre=centres[i]
        # Frozen initial centres remain in the original 17-dimensional basis.
        rep=centre['representative']+[0]*(state.rank-17)
        search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=protocol['coordinate_policy'])
        outcome=search.search(protocol['chart_height'],protocol['seconds_per_chart'],
            checkpoint_dir=folder/'charts'/state.key)
        for point in outcome.curve_points:state=state.adjoin(point,cache=cache,extra_primes=primes)
        result['charts'].append({'centre':centre,'search':outcome.record,'rank_lower_bound':state.rank,
            'admission_prime_bound':admission_prime_bound})
        result.update(final_state=state.record(),arithmetic_facts=cache.store.snapshot(),rank_lower_bound=state.rank)
        checkpoint(output,result)
        print(f'FIBRE {label} t={parameter} chart={i+1}/43 points={len(outcome.curve_points)} rank>={state.rank}',flush=True)
        if state.rank>=32: break
    result['status']='COMPLETE' if len(result['charts'])==43 else 'RANK32_TARGET_REACHED'
    result['coverage_complete']=all(c['search']['status']=='bounded_search_complete' for c in result['charts'])
    checkpoint(output,result)


def verify(path):
    from research_runtime.mw_state import MWState
    result=read(path); cache=ReductionCache(MemoryFactStore()); cache.store.import_snapshot(result['arithmetic_facts'])
    initial=MWState.from_record(result['initial_state'],cache=cache)
    state=initial
    for row in result['charts']:
        primes=tuple(map(int,prime_range(3,row.get('admission_prime_bound',1000)+1)))
        rep=row['centre']['representative']+[0]*(state.rank-17)
        search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy={'kind':'metric','weight':'16'})
        outcome=search.verify_record(row['search'])
        for point in outcome.curve_points: state=state.adjoin(point,cache=cache,extra_primes=primes)
        if state.rank!=row['rank_lower_bound']: raise ArithmeticError('rank admission changed')
    if state.record()!=result['final_state'] or state.rank!=result['rank_lower_bound']:
        raise ArithmeticError('final subgroup certificate changed')
    print('VERIFIED',result['label'],'rank >=',state.rank)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['prepare','run','control','verify'])
    p.add_argument('--directory',type=Path,default=DIRECTORY);p.add_argument('--index',type=int,default=0)
    p.add_argument('--result',type=Path)
    a=p.parse_args()
    if a.stage=='prepare':prepare(a.directory)
    elif a.stage=='verify':verify(a.result)
    elif a.stage=='control':
        protocol=read(a.directory/'protocol.json');run_fibre(a.directory,f'control-{a.index}',protocol['control_parameters'][a.index])
    else:
        pop=read(a.directory/'population.json');protocol=read(a.directory/'protocol.json')
        if identity(protocol)!=pop['protocol_hash']:raise ArithmeticError('population binding changed')
        run_fibre(a.directory,f'candidate-{a.index:02}',pop['finalists'][a.index]['parameter'])


if __name__=='__main__':main()
