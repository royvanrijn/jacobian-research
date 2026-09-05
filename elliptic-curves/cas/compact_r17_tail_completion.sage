#!/usr/bin/env sage-python
"""Finish unvisited denominator tails for four highest-score zero-gain fibres.

The exact chart and coordinate map must agree before a retained prefix can
justify starting later. Completed boxes still give no rational rank upper bound.
"""
import argparse
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path[:0]=[str(CAS),str(ROOT/'elliptic-curves')]
m=SourceFileLoader('tail_common',str(CAS/'compact_r17_prospective.sage')).load_module()
import certify_compact_r17_candidates as certificate
from research_runtime.store import checkpoint
from research_runtime.memory_store import MemoryFactStore
from research_runtime.finite_reduction import ReductionCache
from research_runtime.search_state import raw_state
from pointed_quartic_search import PointedQuarticSearch
from mod2_reduction_independence import _is_prime
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/compact-r17-tail-completion-v1'


def prepare(directory):
    if (directory/'protocol.json').exists():raise FileExistsError('tail cohort already frozen')
    folder=ROOT/'artifacts/local/elliptic-curves/compact-r17-h16384-v1';pop=m.read(folder/'population.json');ledger=m.read(folder/'ledger.json')
    eligible=[]
    for row in ledger['rows']:
        if row.get('rank_lower_bound')==17 and row['status'] in ('COMPLETE','REUSED_COMPLETE_MEASUREMENT'):
            eligible.append({**row,'score_units':pop['finalists'][row['index']]['score_units']})
    eligible.sort(key=lambda r:(-r['score_units'],r['parameter']));rows=[]
    for row in eligible[:4]:
        path=ROOT/row['path'];d=m.read(path)
        if d['status']!='COMPLETE' or len(d['charts'])!=43:raise ArithmeticError('incomplete initial pass')
        rows.append({'parameter':row['parameter'],'input_path':row['path'],'input_sha256':m.hashed(path),
            'score_units':row['score_units'],'old_completed_denominators':[r['search']['completed_denominator'] for r in d['charts']]})
    if len(rows)!=4:raise ArithmeticError('tail cohort incomplete')
    protocol={'schema':'elliptic-curves.compact-r17-tail-completion.v1','sources':m.provenance(),
        'runner_sha256':m.hashed(Path(__file__)),'rows':rows,'selection':'highest full-prime score among completed zero-gain uncatalogued large-cohort fibres',
        'height':100000,'seconds_per_remaining_tail':4,'charts_per_fibre':43,'workers':4,
        'worker_wall_limit_seconds':300,'worker_rss_limit_bytes':1610612736,'admission_prime_bound':251,
        'coordinate_policy':{'kind':'metric','weight':'16'},
        'claim_boundary':'Only coverage and exact point/rank lower bounds. A completed height box is not a rank exclusion.'}
    checkpoint(directory/'protocol.json',protocol);checkpoint(directory/'population.json',{'protocol_hash':m.identity(protocol),'finalists':rows})
    print('TAIL COHORT FROZEN',[(r['parameter'],sum(r['old_completed_denominators'])/4300000) for r in rows],flush=True)


def run(directory,index):
    protocol=m.read(directory/'protocol.json');row=protocol['rows'][index]
    if protocol['sources']!=m.provenance() or protocol['runner_sha256']!=m.hashed(Path(__file__)):raise ArithmeticError('tail source changed')
    path=ROOT/row['input_path']
    if m.hashed(path)!=row['input_sha256']:raise ArithmeticError('retained prefix changed')
    old=m.read(path);model=tuple(map(F,old['curve']));points=[tuple(map(F,p)) for p in old['final_state']['state']['reductions']['points']]
    certificate.checked_rank(model,points);certificate.family_check(row['parameter'],model,points)
    folder=directory/f'candidate-{index:02}';output=folder/'result.json'
    if output.exists():raise FileExistsError('retain this tail attempt; create an explicit continuation if needed')
    cache=ReductionCache(MemoryFactStore());state=raw_state(model,points,cache=cache,prime_bound=1000)
    result={'schema':'elliptic-curves.compact-r17-tail-completion-result.v1','parameter':row['parameter'],'curve':old['curve'],
        'protocol_hash':m.identity(protocol),'initial_state':state.record(),'rank_lower_bound':state.rank,'charts':[],'status':'RUNNING'}
    checkpoint(output,result);primes=tuple(p for p in range(3,252) if _is_prime(p))
    for i,r in enumerate(old['charts']):
        prefix=r['search'];done=prefix['completed_denominator'];last=prefix['denominator_end']
        if prefix['denominator_start']!=1 or last!=100000 or prefix['height_bound']!=100000:raise ArithmeticError('unexpected initial box')
        rep=r['centre']['representative']+[0]*(state.rank-17)
        search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=protocol['coordinate_policy'])
        current=search.chart_record()
        for key in ('base_point','coefficients','horizontal_matrix','ordinate_scale','pointed_chart','short_model','short_model_x_shift'):
            if current[key]!=prefix[key]:raise ArithmeticError('tail chart transport differs: '+key)
        record={'chart_index':i,'centre':r['centre'],'retained_prefix_end':done,'union_completed_denominator':done}
        if done<last:
            outcome=search.search(100000,4,denominator_start=done+1,denominator_end=last,checkpoint_dir=folder/'charts'/state.key)
            for point in outcome.curve_points:state=state.adjoin(point,cache=cache,extra_primes=primes)
            record.update(search=outcome.record,union_completed_denominator=outcome.record['completed_denominator'])
        record['rank_lower_bound']=state.rank;result['charts'].append(record)
        result.update(final_state=state.record(),arithmetic_facts=cache.store.snapshot(),rank_lower_bound=state.rank)
        checkpoint(output,result)
        print('TAIL',row['parameter'],'chart',i+1,'through',record['union_completed_denominator'],'rank',state.rank,flush=True)
    result['status']='COMPLETE_DECLARED_TAIL_PASS';result['full_boxes_completed']=sum(r['union_completed_denominator']==100000 for r in result['charts']);checkpoint(output,result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run']);p.add_argument('--directory',type=Path,default=DIRECTORY);p.add_argument('--index',type=int,default=0);a=p.parse_args()
    if a.stage=='prepare':prepare(a.directory)
    else:run(a.directory,a.index)
