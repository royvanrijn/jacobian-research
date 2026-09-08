#!/usr/bin/env sage-python
"""A frozen 301-chart follow-up on the first independently checked fresh rank-25 atlas curve."""
import argparse
from importlib.machinery import SourceFileLoader
from pathlib import Path
import sys
from sage.all import prime_range

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path[:0]=[str(CAS),str(ROOT/'elliptic-curves')]
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from research_runtime.store import checkpoint,digest
from research_runtime.memory_store import MemoryFactStore
from research_runtime.finite_reduction import ReductionCache
from research_runtime.mw_state import MWState
from pointed_quartic_search import PointedQuarticSearch,sources as pointed_sources
engine=SourceFileLoader('atlas_adaptive_geometry',str(CAS/'half_lattice_fake_descent_replay.sage')).load_module()
geometry=SourceFileLoader('atlas_adaptive_general_geometry',str(CAS/'run_curve385_iterated_half_lattice_search.sage')).load_module()
PARENT=ROOT/'artifacts/local/elliptic-curves/compact-six-r17-h1024-v2'
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/compact-atlas-first25-followup-v2'


def sources():
    paths=(Path(__file__).resolve(),Path(spec.__file__).resolve(),spec.ATLAS,Path(cert.__file__).resolve(),
           CAS/'half_lattice_fake_descent_replay.sage',CAS/'run_curve385_iterated_half_lattice_search.sage',CAS/'research_runtime/memory_store.py')
    return {**pointed_sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}}


def prepare(directory):
    if (directory/'protocol.json').exists():raise FileExistsError('follow-up already frozen')
    path=PARENT/'07ca9/candidate-01/result.json';result=cert.read(path);check=cert.read(path.parent/'independent-rank-check.json')
    if result['status']!='COMPLETE_DECLARED_PILOT' or result['rank_lower_bound']!=25 or check['source_sha256']!=cert.hashed(path):raise ArithmeticError('rank-25 starting point changed')
    family=next(r for r in cert.read(spec.ATLAS)['families'] if r['family']==result['family'])
    model=tuple(map(cert.F,result['curve']));points=tuple(tuple(map(cert.F,p)) for p in result['final_state']['state']['reductions']['points'])
    proof=check['rank_certificate'];cert.checked_rank(model,points,[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
    spec.family_check(family,result['parameter'],model,points)
    masks=cert.read(PARENT/'protocol.json')['generic_masks'][result['family']]
    oracle=engine.CosetOracle(family['generic_height_gram']);pool=[];trials=[]
    for mask in masks:
        norm,rep,error=oracle.solve(mask)
        if norm!=12:raise ArithmeticError('deep generic mask changed')
        pool.append({'mask':mask,'generic_norm':norm})
    for mask in range(1<<17):
        if mask in masks:continue
        norm,rep,error=oracle.solve(mask);trials.append([mask,norm])
        if norm==10:pool.append({'mask':mask,'generic_norm':norm})
        if len(pool)==301:break
    if len(pool)!=301:raise ArithmeticError('generic next-shell pool incomplete')
    checkpoint(directory/'protocol.json',{'schema':'elliptic-curves.compact-atlas-first25-followup.v1','sources':sources(),
        'input_path':str(path.relative_to(ROOT)),'input_sha256':cert.hashed(path),'independent_check_sha256':cert.hashed(path.parent/'independent-rank-check.json'),
        'family':result['family'],'parameter':result['parameter'],'initial_rank_lower_bound':25,
        'selection':'first completed prospective atlas candidate with independently replayed rank at least 25; choose this sole candidate before its adaptive outcomes',
        'generic_pool':pool,'next_shell_trials':trials,'quotient_rule':'cyclic nonzero words ordered by Hamming weight then integer, crossed with the frozen 301 generic classes; rank by specialized numerical CVP norm',
        'charts':301,'chart_height':100000,'seconds_per_chart':4,'admission_prime_bound':251,
        'coordinate_policy':{'kind':'metric','weight':'16'},'workers':1,'worker_wall_seconds':1500,'worker_rss_bytes':1610612736,
        'target_rank_lower_bound':28,'scope':'One additional prospective follow-up; original nineteen-fibre initial batch continues under its own limits. No public exceptional points. Bounded coverage is not an upper rank bound.'})
    print('FROZEN FIRST25 FOLLOWUP',result['family'],result['parameter'],flush=True)


def run(directory):
    protocol=cert.read(directory/'protocol.json');path=ROOT/protocol['input_path'];source=cert.read(path)
    if protocol['sources']!=sources() or cert.hashed(path)!=protocol['input_sha256']:raise ArithmeticError('follow-up source changed')
    folder=directory/'candidate-00';output=folder/'result.json'
    model=tuple(map(cert.F,source['curve']));cache=ReductionCache(MemoryFactStore());cache.store.import_snapshot(source['arithmetic_facts'])
    initial=MWState.from_record(source['final_state'],cache=cache);dimension=initial.rank
    if dimension!=25:raise ArithmeticError('initial rank changed')
    if output.exists():
        result=cert.read(output)
        if result['protocol_hash']!=digest(protocol) or result['initial_state']!=initial.record():raise ArithmeticError('resume mismatch')
        if result['status']!='RUNNING':print('RETAINED FOLLOWUP',result['status'],flush=True);return
        cache.store.import_snapshot(result['arithmetic_facts']);state=MWState.from_record(result['final_state'],cache=cache)
    else:
        points=tuple(tuple(map(cert.F,p)) for p in initial.basis);gram,asymmetry=geometry.canonical_height_gram(model,points)
        oracle=geometry.CosetOracle(geometry.rounded_gram(gram,1000000))
        words=sorted(range(1,1<<(dimension-17)),key=lambda w:(w.bit_count(),w));centres=[]
        for i,g in enumerate(protocol['generic_pool']):
            word=words[i%len(words)];mask=g['mask']|(word<<17);residue=[(mask>>j)&1 for j in range(dimension)]
            norm,rep,error=oracle.solve(residue)
            if len(rep)!=dimension or any((rep[j]-residue[j])%2 for j in range(dimension)):raise ArithmeticError('adaptive centre dimension or parity mismatch')
            centres.append({'generic_mask':g['mask'],'quotient_word':word,'parity':mask,'representative':list(rep),'metric_norm':norm})
        centres.sort(key=lambda r:(-r['metric_norm'],r['generic_mask'],r['quotient_word']));state=initial
        result={'schema':'elliptic-curves.compact-atlas-adaptive-result.v1','family':protocol['family'],'parameter':protocol['parameter'],
            'curve':source['curve'],'protocol_hash':digest(protocol),'initial_dimension':dimension,'generic_points':source['generic_points'],
            'initial_state':initial.record(),'metric_gram':[[str(x) for x in row] for row in gram],
            'centres':centres,'charts':[],'status':'RUNNING','rank_lower_bound':dimension,'final_state':state.record(),'arithmetic_facts':cache.store.snapshot()}
        checkpoint(output,result)
    primes=tuple(map(int,prime_range(3,252)))
    for i in range(len(result['charts']),301):
        c=result['centres'][i];rep=c['representative']+[0]*(state.rank-dimension)
        search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=protocol['coordinate_policy'])
        outcome=search.search(100000,4,checkpoint_dir=folder/'charts'/state.key)
        for point in outcome.curve_points:state=state.adjoin(point,cache=cache,extra_primes=primes)
        result['charts'].append({'centre':c,'search':outcome.record,'rank_lower_bound':state.rank,'admission_prime_bound':251})
        result.update(final_state=state.record(),arithmetic_facts=cache.store.snapshot(),rank_lower_bound=state.rank);checkpoint(output,result)
        print('ATLAS FOLLOWUP',i+1,'rank',state.rank,flush=True)
        if state.rank>=28:result['status']='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY';checkpoint(output,result);return
    result['status']='COMPLETE_DECLARED_PILOT';checkpoint(output,result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run']);p.add_argument('--directory',type=Path,default=DIRECTORY);a=p.parse_args()
    prepare(a.directory) if a.stage=='prepare' else run(a.directory)
