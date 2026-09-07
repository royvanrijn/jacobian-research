#!/usr/bin/env sage-python
"""Fixed48-chart extension of the generic classes omitted by the43-chart cap."""
import argparse
from importlib.machinery import SourceFileLoader
from pathlib import Path
import sys
from sage.all import prime_range
ROOT=Path(__file__).resolve().parents[2]; CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from research_runtime.store import checkpoint,digest
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.cached_observation_state import CachedObservationMWState as MWState
from research_runtime.pointed_orbit_compression import compress
from pointed_quartic_search import PointedQuarticSearch,sources as pointed_sources
geometry=SourceFileLoader('omitted_r17_geometry',str(CAS/'prospective_half_lattice_v2.sage')).load_module()
PARENT=ROOT/'artifacts/local/elliptic-curves/compact-six-r17-h4096-v1'
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/r17-omitted-generic-classes-v1'


def sources():
    paths=[Path(__file__).resolve(),Path(spec.__file__).resolve(),spec.ATLAS,Path(cert.__file__).resolve(),
           CAS/'prospective_half_lattice_v2.sage',CAS/'research_runtime/cached_observation_state.py',
           CAS/'research_runtime/quotient_only_reduction.py',CAS/'research_runtime/memory_store.py',
           CAS/'research_runtime/pointed_orbit_compression.py',CAS/'alternate_quartic_covers.py']
    return {**pointed_sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}}


def prepare(directory):
    if (directory/'protocol.json').exists():raise FileExistsError('preserve omitted-class protocol')
    ledger=cert.read(PARENT/'point-ledger.json');replay=cert.read(PARENT/'verification/ledger.json')
    if ledger['status']!='COMPLETE_FIXED_BATCH_ATTEMPTS' or len(ledger['rows'])!=24:raise ArithmeticError('parent batch unfinished')
    if len(replay['rows'])!=24 or any(r['status']!='PASS' for r in replay['rows']):raise ArithmeticError('parent replay incomplete')
    parent_protocol=cert.read(PARENT/'point-protocol.json');families=['11952','08f72'];roster=[];classes={};censuses={}
    for family in families:
        path=PARENT/family/'generic-census.json';census=cert.read(path)
        if cert.hashed(path)!=parent_protocol['generic_census_hashes'][family] or census['status']!='COMPLETE_DECLARED_CENSUS':raise ArithmeticError('generic census changed')
        used=set(parent_protocol['generic_masks'][family]);chosen=[r['mask'] for r in census['records'] if cert.F(r['norm'])==12 and r['mask'] not in used]
        if len(used)!=43 or len(chosen)!=6:raise ArithmeticError('omitted-class gate changed')
        classes[family]=chosen;censuses[family]={'path':str(path.relative_to(ROOT)),'sha256':cert.hashed(path)}
        for index in range(4):
            entry=next(r for r in ledger['rows'] if r['family']==family and r['index']==index)
            p=ROOT/entry['result_path'];data=cert.read(p)
            if cert.hashed(p)!=entry['result_sha256'] or len(data['charts'])!=43 or data['status']!='COMPLETE_DECLARED_PILOT':raise ArithmeticError('parent43 input changed')
            roster.append({'family':family,'index':index,'parameter':data['parameter'],'input_path':str(p.relative_to(ROOT)),'input_sha256':cert.hashed(p)})
    checkpoint(directory/'protocol.json',{'schema':'elliptic-curves.r17-omitted-generic-classes.v1','sources':sources(),
        'parent_ledger_sha256':cert.hashed(PARENT/'point-ledger.json'),'parent_replay_sha256':cert.hashed(PARENT/'verification/ledger.json'),
        'roster':roster,'classes':classes,'censuses':censuses,'charts_per_candidate':6,'height':100000,'seconds_per_chart':4,
        'admission_prime_bound':997,'coordinate_policy':{'kind':'metric','weight':'16'},'worker_wall_seconds':90,'worker_rss_bytes':1073741824,'maximum_workers':2,
        'metric':'Reuse each parent recorded384-bit numerical generic17 height Gram; round its exact decimal entries at scale1e6, choose exact-parity CVP representatives, order by computed specialized norm then mask. No new public data or known-record points.',
        'admission':'Keep all raw points. Admit one representative of each pointed involution orbit, retaining the exact integral centre word and partner relation. Audit the entire new cloud independently through997.',
        'mathematical_gate':'Two fresh generic censuses have49 computed-norm12 representatives while the fixed parent cap searched43. Test all six omitted classes on all four already fixed addresses in each family, independent of measured rank. This tests conditional point visibility, not candidate incidence or global Selmer-cover solubility.',
        'claim_boundary':'Exactly48 further bounded charts. Computed generic norms are not certified minima. No catalogue access or replacement in the worker. No upper rank, exhaustive height-box claim or novelty assertion follows from a finite miss.'})
    print('FROZEN OMITTED R17 CLASSES',classes,flush=True)


def run(directory,family,index):
    protocol=cert.read(directory/'protocol.json')
    if protocol['sources']!=sources():raise ArithmeticError('frozen worker source changed')
    entry=next(r for r in protocol['roster'] if r['family']==family and r['index']==index)
    path=ROOT/entry['input_path'];source=cert.read(path)
    if cert.hashed(path)!=entry['input_sha256']:raise ArithmeticError('parent point input changed')
    output=directory/family/f'candidate-{index:02}'/'result.json';folder=output.parent
    if output.exists():raise FileExistsError('preserve omitted-class attempt')
    model=tuple(map(cert.F,source['curve']));generic=tuple(tuple(map(cert.F,p)) for p in source['generic_points'])
    f=next(f for f in cert.read(spec.ATLAS)['families'] if f['family']==family)
    spec.family_check(f,entry['parameter'],model,generic)
    cache=ReductionCache(MemoryFactStore());cache.store.import_snapshot(source['arithmetic_facts'])
    state=MWState.from_record(source['final_state'],cache=cache)
    if tuple(tuple(map(cert.F,p)) for p in state.basis[:17])!=generic:raise ArithmeticError('generic prefix changed')
    gram=tuple(tuple(round(cert.F(x)*1000000) for x in row) for row in source['metric_gram'])
    oracle=geometry.CosetOracle(gram);centres=[]
    for mask in protocol['classes'][family]:
        norm,rep,error=oracle.solve(tuple((mask>>j)&1 for j in range(17)))
        if len(rep)!=17 or any((rep[j]-(mask>>j))%2 for j in range(17)):raise ArithmeticError('parity representative failed')
        centres.append({'mask':mask,'representative':list(rep),'metric_norm':norm,'cvp_error':error})
    centres.sort(key=lambda r:(-r['metric_norm'],r['mask']))
    result={'schema':'elliptic-curves.r17-omitted-generic-result.v1','protocol_hash':digest(protocol),'family':family,'parameter':entry['parameter'],
        'curve':source['curve'],'generic_points':source['generic_points'],'family_to_curve_scale_u':source['family_to_curve_scale_u'],
        'parent_input':entry,'initial_state':state.record(),'initial_dimension':state.rank,'metric_gram':source['metric_gram'],
        'centres':centres,'charts':[],'status':'RUNNING','rank_lower_bound':state.rank,'final_state':state.record(),'arithmetic_facts':cache.store.snapshot()}
    checkpoint(output,result);primes=tuple(map(int,prime_range(3,998)))
    for centre in centres:
        rep=centre['representative']+[0]*(state.rank-17)
        search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=protocol['coordinate_policy'])
        outcome=search.search(100000,4,checkpoint_dir=folder/'charts'/state.key)
        compression=compress(model,state.basis,rep,outcome.curve_points)
        for i in compression['kept_indices']:
            state=state.adjoin(outcome.curve_points[i],cache=cache,extra_primes=primes)
            if not isinstance(state,MWState):state=MWState.from_record(state.record(),cache=cache)
        result['charts'].append({'centre':centre,'search':outcome.record,'admission_compression':compression,'rank_lower_bound':state.rank,'admission_prime_bound':997})
        result.update(rank_lower_bound=state.rank,final_state=state.record(),arithmetic_facts=cache.store.snapshot());checkpoint(output,result)
        print('OMITTED R17',family,index,len(result['charts']),'rank',state.rank,flush=True)
    result['status']='COMPLETE_DECLARED_PILOT';checkpoint(output,result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run']);p.add_argument('--directory',type=Path,default=DIRECTORY);p.add_argument('--family');p.add_argument('--index',type=int);a=p.parse_args()
    prepare(a.directory) if a.stage=='prepare' else run(a.directory,a.family,a.index)
