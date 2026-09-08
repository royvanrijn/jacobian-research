#!/usr/bin/env sage-python
"""Frozen, two-round sparse quotient-centre pilot on certified fresh R17 curves.

Only certified prospective points enter the enlarged lattice. Numerical height
and floating CVP schedule charts; exact finite quotients alone admit new rank.
"""
import argparse
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
from itertools import combinations
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
CAS=ROOT/'elliptic-curves/cas'
sys.path[:0]=[str(CAS),str(ROOT/'elliptic-curves')]
from sage.all import prime_range
from research_runtime.store import checkpoint
from research_runtime.memory_store import MemoryFactStore
from research_runtime.finite_reduction import ReductionCache
from research_runtime.search_state import raw_state
from research_runtime.mw_state import MWState
from pointed_quartic_search import PointedQuarticSearch
import certify_compact_r17_candidates as certificate

m=SourceFileLoader('compact_adaptive_common',str(CAS/'compact_r17_prospective.sage')).load_module()
engine=SourceFileLoader('compact_adaptive_geometry',str(CAS/'run_curve385_iterated_half_lattice_search.sage')).load_module()
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/compact-r17-adaptive-v2'
INPUTS=[ROOT/'artifacts/generated-results/elliptic-curves'/f for f in
        ('compact_r17_new_curves_v1.json','compact_r17_wide_new_curves_v1.json')]


def provenance():
    return {**m.provenance(),**{str(p.relative_to(ROOT)):m.hashed(p) for p in
        (Path(__file__),CAS/'run_curve385_iterated_half_lattice_search.sage',CAS/'certify_compact_r17_candidates.py')}}


def prepare(directory):
    if (directory/'protocol.json').exists():raise FileExistsError('use the frozen protocol')
    rows=[]
    for path in INPUTS:
        certificate.check(path)
        for row in m.read(path)['curves']:
            if len(row['points'])==24 and not row['icarm_snapshot_isomorphism_matches']:
                rows.append({'parameter':row['parameter'],'curve':row['curve'],'points':row['points'],
                    'certificate':str(path.relative_to(ROOT)),'certificate_sha256':m.hashed(path)})
    if [r['parameter'] for r in rows]!=['33/119','-695/97']:raise ArithmeticError('pilot cohort changed')
    protocol={'schema':'elliptic-curves.compact-r17-adaptive-protocol.v1','sources':provenance(),
        'rows':rows,'rounds':2,'retained_charts_per_round':128,'chart_height':100000,'seconds_per_chart':3,
        'admission_prime_bound':251,'coordinate_policy':{'kind':'metric','weight':'16'},
        'pool':'43 retained generic deep masks, and zero, crossed with all weight-one/two words in the current certified quotient basis',
        'ranking':'decreasing norm of floating CVP representative for specialized numerical height rounded at 1e6; ties by integer parity mask',
        'round_refresh':'After 128 charts rebuild the enlarged geometry only if rank increased; otherwise stop this bounded pilot.',
        'budget':'At most two workers, 1100 wall seconds and 1.5 GiB each, including geometry and all charts.',
        'target_rank_lower_bound':28,'claim_boundary':'Prospective lower bounds only; neither CVP optimality nor search completeness is a proof gate.'}
    checkpoint(directory/'protocol.json',protocol)
    checkpoint(directory/'population.json',{'protocol_hash':m.identity(protocol),'finalists':rows})
    print('FROZEN',[(r['parameter'],len(r['points'])) for r in rows],flush=True)


def run(directory,index):
    protocol=m.read(directory/'protocol.json')
    if protocol['sources']!=provenance():raise ArithmeticError('adaptive sources changed')
    row=protocol['rows'][index]
    if m.hashed(ROOT/row['certificate'])!=row['certificate_sha256']:raise ArithmeticError('input certificate changed')
    folder=directory/f'candidate-{index:02}';output=folder/'result.json'
    cache=ReductionCache(MemoryFactStore())
    model=tuple(map(F,row['curve']));points=tuple(tuple(map(F,p)) for p in row['points'])
    certificate.checked_rank(model,points)
    certificate.family_check(row['parameter'],model,points)
    state=raw_state(model,points,cache=cache,prime_bound=1000)
    if state.basis!=tuple(tuple(map(str,p)) for p in points) or state.rank!=24:raise ArithmeticError('certified starting basis changed')
    result={'schema':'elliptic-curves.compact-r17-adaptive-result.v1','parameter':row['parameter'],
        'curve':row['curve'],'generic_points':row['points'][:17],'protocol_hash':m.identity(protocol),
        'initial_state':state.record(),'rank_lower_bound':state.rank,'charts':[],'rounds':[],
        'status':'RUNNING','final_state':state.record(),'arithmetic_facts':cache.store.snapshot()}
    if output.exists():
        saved=m.read(output)
        if saved['protocol_hash']!=result['protocol_hash'] or saved['initial_state']!=result['initial_state']:
            raise ArithmeticError('resume input changed')
        if saved['status']!='RUNNING':print('RETAINED',row['parameter'],saved['rank_lower_bound'],flush=True);return
        cache.store.import_snapshot(saved['arithmetic_facts'])
        state=MWState.from_record(saved['final_state'],cache=cache);result=saved
    checkpoint(output,result)
    masks=[0]+sorted(c['mask'] for c in m.read(m.HOLES)['fibres'][0]['cover_records'])
    primes=tuple(map(int,prime_range(3,protocol['admission_prime_bound']+1)))
    for round_index in range(protocol['rounds']):
        if round_index>=len(result['rounds']):
            if round_index and state.rank<=result['rounds'][-1]['rank_before']:break
            dimension=state.rank
            basis=tuple(tuple(map(F,p)) for p in state.basis)
            if basis[:17]!=points[:17]:raise ArithmeticError('generic prefix changed')
            gram,asymmetry=engine.canonical_height_gram(model,basis)
            rounded=engine.rounded_gram(gram,1000000);oracle=engine.CosetOracle(rounded)
            pool=[]
            for weight in (1,2):
                for bits in combinations(range(17,dimension),weight):
                    quotient=sum(1<<bit for bit in bits)
                    for mask in masks:
                        parity=mask|quotient
                        norm,rep,error=oracle.solve([(parity>>j)&1 for j in range(dimension)])
                        if any((rep[j]&1)!=((parity>>j)&1) for j in range(dimension)):
                            raise ArithmeticError('CVP changed the requested parity')
                        pool.append({'parity':parity,'representative':list(rep),'metric_norm':norm,'cvp_error':error})
            pool.sort(key=lambda c:(-c['metric_norm'],c['parity']))
            plan={'round':round_index,'rank_before':dimension,'state_key':state.key,
                'metric_gram':[[str(c) for c in line] for line in gram],
                'maximum_asymmetry':str(asymmetry),'pool':pool,
                'centres':pool[:protocol['retained_charts_per_round']],'completed':0}
            result['rounds'].append(plan);checkpoint(output,result)
            print('GEOMETRY',row['parameter'],'round',round_index,'rank',dimension,'pool',len(pool),flush=True)
        plan=result['rounds'][round_index]
        for chart_index in range(plan['completed'],len(plan['centres'])):
            centre=plan['centres'][chart_index]
            rep=centre['representative']+[0]*(state.rank-plan['rank_before'])
            search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=protocol['coordinate_policy'])
            outcome=search.search(protocol['chart_height'],protocol['seconds_per_chart'],checkpoint_dir=folder/'charts'/state.key)
            for point in outcome.curve_points:state=state.adjoin(point,cache=cache,extra_primes=primes)
            result['charts'].append({'round':round_index,'chart_index':chart_index,'centre':centre,
                'search':outcome.record,'rank_lower_bound':state.rank,'admission_prime_bound':protocol['admission_prime_bound']})
            plan['completed']=chart_index+1
            result.update(final_state=state.record(),arithmetic_facts=cache.store.snapshot(),rank_lower_bound=state.rank)
            checkpoint(output,result)
            print('ADAPTIVE',row['parameter'],'round',round_index,'chart',chart_index+1,'points',len(outcome.curve_points),'rank',state.rank,flush=True)
            if state.rank>=protocol['target_rank_lower_bound']:
                result['status']='TARGET_REACHED_PENDING_INDEPENDENT_EXPORT';checkpoint(output,result);return
    result['status']='COMPLETE_DECLARED_PILOT';checkpoint(output,result)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run'])
    p.add_argument('--directory',type=Path,default=DIRECTORY);p.add_argument('--index',type=int,default=0);a=p.parse_args()
    if a.stage=='prepare':prepare(a.directory)
    else:run(a.directory,a.index)


if __name__=='__main__':main()
