#!/usr/bin/env sage-python
"""Three separately frozen 301-chart follow-ups on the other certified rank-23/24 atlas curves."""
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
CERTIFICATE=ROOT/'artifacts/generated-results/elliptic-curves/compact_atlas_new_curves_v1.json'
PARENT=ROOT/'artifacts/local/elliptic-curves/compact-six-r17-h1024-v2'
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/compact-atlas-other-gains-followup-v1'


def sources():
    paths=(Path(__file__).resolve(),Path(spec.__file__).resolve(),spec.ATLAS,CERTIFICATE,Path(cert.__file__).resolve(),
           CAS/'half_lattice_fake_descent_replay.sage',CAS/'run_curve385_iterated_half_lattice_search.sage',CAS/'research_runtime/memory_store.py')
    return {**pointed_sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}}


def prepare_case(directory,path,certified):
    if (directory/'protocol.json').exists():raise FileExistsError('follow-up already frozen')
    result=cert.read(path);check={'source_sha256':certified['discovery_witness']['sha256'],'rank_certificate':certified['rank_certificate']}
    if result['status']!='COMPLETE_DECLARED_PILOT' or not 23<=result['rank_lower_bound']<=24 or check['source_sha256']!=cert.hashed(path):raise ArithmeticError('rank-23/24 starting point changed')
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
    checkpoint(directory/'protocol.json',{'schema':'elliptic-curves.compact-atlas-other-gain-followup.v1','sources':sources(),
        'input_path':str(path.relative_to(ROOT)),'input_sha256':cert.hashed(path),'independent_check_sha256':cert.hashed(CERTIFICATE),
        'family':result['family'],'parameter':result['parameter'],'initial_rank_lower_bound':result['rank_lower_bound'],
        'selection':'all three remaining independently certified compact-atlas candidates of rank at least 23 and below 25, selected together before their adaptive outcomes',
        'generic_pool':pool,'next_shell_trials':trials,'quotient_rule':'cyclic nonzero words ordered by Hamming weight then integer, crossed with the frozen 301 generic classes; rank by specialized numerical CVP norm',
        'charts':301,'chart_height':100000,'seconds_per_chart':4,'admission_prime_bound':251,
        'coordinate_policy':{'kind':'metric','weight':'16'},'workers':1,'worker_wall_seconds':1500,'worker_rss_bytes':1610612736,
        'target_rank_lower_bound':28,'scope':'One of three additional prospective follow-ups, each 301 charts and 1500 seconds/1.5 GiB. At most three run concurrently, alongside the separately frozen first rank-25 follow-up. No public exceptional points. Bounded coverage is not an upper rank bound.'})
    print('FROZEN OTHER-GAIN FOLLOWUP',result['family'],result['parameter'],flush=True)


def prepare(directory):
    if (directory/'cohort.json').exists():raise FileExistsError('cohort already frozen')
    data=cert.read(CERTIFICATE)
    for name,h in data['sources'].items():
        if cert.hashed(ROOT/name)!=h:raise ArithmeticError('point certificate dependency changed')
    selected=sorted((r for r in data['curves'] if 23<=r['rank_certificate']['rank_lower_bound']<25),
                    key=lambda r:(-r['rank_certificate']['rank_lower_bound'],r['family'],r['parameter']))
    if len(selected)!=3:raise ArithmeticError('declared three-candidate cohort changed')
    cases=[]
    for i,r in enumerate(selected):
        case=directory/f'case-{i:02}'
        prepare_case(case,ROOT/r['discovery_witness']['path'],r)
        cases.append({'index':i,'family':r['family'],'parameter':r['parameter'],
            'rank_lower_bound':r['rank_certificate']['rank_lower_bound'],'protocol_sha256':cert.hashed(case/'protocol.json')})
    checkpoint(directory/'cohort.json',{'schema':'elliptic-curves.compact-atlas-other-gains-cohort.v1',
        'sources':sources(),'certificate_sha256':cert.hashed(CERTIFICATE),'cases':cases,
        'maximum_concurrent_workers':3,'charts_per_worker':301,'worker_wall_seconds':1500,'worker_rss_bytes':1610612736,
        'selection':'All three newly certified atlas curves with initial lower bound at least 23 and below 25. The first rank-25 curve has its own earlier follow-up. No adaptive outcome selected this cohort.'})


def run(directory):
    protocol=cert.read(directory/'protocol.json');path=ROOT/protocol['input_path'];source=cert.read(path)
    if protocol['sources']!=sources() or cert.hashed(path)!=protocol['input_sha256']:raise ArithmeticError('follow-up source changed')
    folder=directory/'candidate-00';output=folder/'result.json'
    model=tuple(map(cert.F,source['curve']));cache=ReductionCache(MemoryFactStore());cache.store.import_snapshot(source['arithmetic_facts'])
    initial=MWState.from_record(source['final_state'],cache=cache);dimension=initial.rank
    if dimension!=protocol['initial_rank_lower_bound']:raise ArithmeticError('initial rank changed')
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
