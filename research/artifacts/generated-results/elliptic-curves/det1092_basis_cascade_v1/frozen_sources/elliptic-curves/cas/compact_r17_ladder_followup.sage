#!/usr/bin/env sage-python
"""Broader 301-class follow-up, starting only with search-recovered directions.

The generic-class/quotient-word rule follows the earlier blind ladder. This
pilot uses the current universal sieve at four seconds per chart, and must
not inherit the old implementation's coverage or success claims.
"""
import argparse
from importlib.machinery import SourceFileLoader
from pathlib import Path
from fractions import Fraction as F
import sys

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path[:0]=[str(CAS),str(ROOT/'elliptic-curves')]
m=SourceFileLoader('ladder_common',str(CAS/'compact_r17_prospective.sage')).load_module()
engine=SourceFileLoader('ladder_geometry',str(CAS/'run_curve385_iterated_half_lattice_search.sage')).load_module()
generic=SourceFileLoader('ladder_generic',str(CAS/'half_lattice_fake_descent_replay.sage')).load_module()
import certify_compact_r17_candidates as certificate
from research_runtime.store import checkpoint
from research_runtime.memory_store import MemoryFactStore
from research_runtime.finite_reduction import ReductionCache
from research_runtime.search_state import raw_state
from research_runtime.mw_state import MWState
from pointed_quartic_search import PointedQuarticSearch
from mod2_reduction_independence import _is_prime
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/compact-r17-ladder-followup-v1'


def sources():
    return {**m.provenance(),**{str(p.relative_to(ROOT)):m.hashed(p) for p in
        (Path(__file__),CAS/'run_curve385_iterated_half_lattice_search.sage',CAS/'certify_compact_r17_candidates.py')}}


def prepare(directory):
    if (directory/'protocol.json').exists():raise FileExistsError('frozen follow-up already exists')
    initial_paths=[('known_rank28_control','compact-r17-h16384-v1/control-rank28/result.json'),
        ('prospective','compact-r17-prospective-v1/candidate-01/result.json'),
        ('prospective','compact-r17-wide-v1/candidate-03/result.json'),
        ('prospective','compact-r17-top64-v1/candidate-33/result.json')]
    rows=[]
    for role,name in initial_paths:
        path=ROOT/'artifacts/local/elliptic-curves'/name;d=m.read(path)
        model=tuple(map(F,d['curve']));points=[tuple(map(F,p)) for p in d['final_state']['state']['reductions']['points']]
        rank=certificate.checked_rank(model,points);certificate.family_check(d['parameter'],model,points)
        if d['status']!='COMPLETE' or len(points)<24:raise ArithmeticError('follow-up input incomplete or below gate')
        rows.append({'role':role,'parameter':d['parameter'],'curve':d['curve'],'points':[list(map(str,p)) for p in points],
            'rank_certificate':rank,'input_path':str(path.relative_to(ROOT)),'input_sha256':m.hashed(path)})
    # The retained generic census has exactly 43 norm-12 maxima. Add the
    # first 258 norm-10 masks in mask order, without resieving all 2^17 classes.
    masks=sorted(c['mask'] for c in m.read(m.HOLES)['fibres'][0]['cover_records'])
    oracle=generic.CosetOracle(generic.GENERIC_GRAM);top=[];trials=[]
    for mask in masks:
        norm,rep,error=oracle.solve(mask)
        if norm!=12:raise ArithmeticError('published generic maximum-mask mismatch')
        top.append({'mask':mask,'generic_norm':norm,'generic_representative':list(rep)})
    for mask in range(1<<17):
        if mask in masks:continue
        norm,rep,error=oracle.solve(mask);trials.append([mask,norm])
        if norm==10:top.append({'mask':mask,'generic_norm':norm,'generic_representative':list(rep)})
        if len(top)==301:break
    if len(top)!=301:raise ArithmeticError('generic follow-up pool incomplete')
    protocol={'schema':'elliptic-curves.compact-r17-ladder-followup.v1','sources':sources(),'rows':rows,
        'generic_centres':top,'generic_next_shell_trials':trials,
        'pool_rule':'first 43 retained generic norm-12 masks, then first 258 norm-10 masks; cyclic nonzero quotient words ordered by Hamming weight and integer value',
        'specialized_ranking':'decreasing specialized numerical CVP norm at 1e6; ties by generic mask and quotient word',
        'charts':301,'chart_height':100000,'seconds_per_chart':4,'admission_prime_bound':251,
        'coordinate_policy':{'kind':'metric','weight':'16'},'workers':4,'worker_wall_limit_seconds':1500,
        'worker_rss_limit_bytes':1610612736,'target_rank_lower_bound':28,
        'claim_boundary':'Control recovery is not novelty. Bounded chart coverage is not inherited from the earlier ladder. Every new lower bound needs exact independent replay.'}
    checkpoint(directory/'protocol.json',protocol)
    checkpoint(directory/'population.json',{'protocol_hash':m.identity(protocol),'finalists':rows})
    print('FOLLOWUP FROZEN',[(r['role'],r['parameter'],len(r['points'])) for r in rows],flush=True)


def run(directory,index):
    protocol=m.read(directory/'protocol.json')
    if protocol['sources']!=sources():raise ArithmeticError('follow-up sources changed')
    row=protocol['rows'][index]
    if m.hashed(ROOT/row['input_path'])!=row['input_sha256']:raise ArithmeticError('follow-up input changed')
    model=tuple(map(F,row['curve']));points=tuple(tuple(map(F,p)) for p in row['points']);dimension=len(points)
    certificate.checked_rank(model,points);certificate.family_check(row['parameter'],model,points)
    cache=ReductionCache(MemoryFactStore());state=raw_state(model,points,cache=cache,prime_bound=1000)
    if state.basis!=tuple(tuple(map(str,p)) for p in points):raise ArithmeticError('starting basis changed')
    folder=directory/f'candidate-{index:02}';output=folder/'result.json'
    if output.exists():
        result=m.read(output)
        if result['protocol_hash']!=m.identity(protocol) or result['initial_state']!=state.record():raise ArithmeticError('resume input changed')
        if result['status']!='RUNNING':print('RETAINED',row['parameter'],result['rank_lower_bound'],flush=True);return
        cache.store.import_snapshot(result['arithmetic_facts']);state=MWState.from_record(result['final_state'],cache=cache)
    else:
        gram,asymmetry=engine.canonical_height_gram(model,points)
        oracle=engine.CosetOracle(engine.rounded_gram(gram,1000000))
        words=sorted(range(1,1<<(dimension-17)),key=lambda w:(w.bit_count(),w));pool=[]
        for i,c in enumerate(protocol['generic_centres']):
            word=words[i%len(words)];mask=c['mask']|(word<<17)
            norm,rep,error=oracle.solve([(mask>>j)&1 for j in range(dimension)])
            pool.append({'mask':c['mask'],'quotient_word':word,'parity':mask,'representative':list(rep),
                'metric_norm':norm,'cvp_error':error})
        pool.sort(key=lambda c:(-c['metric_norm'],c['mask'],c['quotient_word']))
        result={'schema':'elliptic-curves.compact-r17-ladder-followup-result.v1','parameter':row['parameter'],'role':row['role'],
            'curve':row['curve'],'generic_points':row['points'][:17],'protocol_hash':m.identity(protocol),'initial_state':state.record(),
            'initial_dimension':dimension,'rank_lower_bound':dimension,'charts':[],'status':'RUNNING',
            'rounds':[{'round':0,'rank_before':dimension,'state_key':state.key,'metric_gram':[[str(v) for v in r] for r in gram],
                'maximum_asymmetry':str(asymmetry),'pool':pool,'centres':pool,'completed':0}],
            'final_state':state.record(),'arithmetic_facts':cache.store.snapshot()}
        checkpoint(output,result)
    plan=result['rounds'][0];primes=tuple(p for p in range(3,252) if _is_prime(p))
    for i in range(plan['completed'],len(plan['centres'])):
        c=plan['centres'][i];rep=c['representative']+[0]*(state.rank-dimension)
        search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=protocol['coordinate_policy'])
        outcome=search.search(protocol['chart_height'],protocol['seconds_per_chart'],checkpoint_dir=folder/'charts'/state.key)
        for point in outcome.curve_points:state=state.adjoin(point,cache=cache,extra_primes=primes)
        result['charts'].append({'round':0,'chart_index':i,'centre':c,'search':outcome.record,
            'rank_lower_bound':state.rank,'admission_prime_bound':251})
        plan['completed']=i+1;result.update(final_state=state.record(),arithmetic_facts=cache.store.snapshot(),rank_lower_bound=state.rank)
        checkpoint(output,result)
        print('LADDER',row['role'],row['parameter'],'chart',i+1,'points',len(outcome.curve_points),'rank',state.rank,flush=True)
        if state.rank>=28:result['status']='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY';checkpoint(output,result);return
    result['status']='COMPLETE_DECLARED_PILOT';checkpoint(output,result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run']);p.add_argument('--directory',type=Path,default=DIRECTORY);p.add_argument('--index',type=int,default=0);a=p.parse_args()
    if a.stage=='prepare':prepare(a.directory)
    else:run(a.directory,a.index)
