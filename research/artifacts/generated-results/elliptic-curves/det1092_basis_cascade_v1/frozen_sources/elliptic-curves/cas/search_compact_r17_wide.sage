#!/usr/bin/env sage-python
"""Fixed balanced H4096 point batch on six compact R17 families."""
import argparse
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
from pathlib import Path
import shutil
import sys
from sage.all import EllipticCurve, QQ, prime_range
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import compact_atlas_specialization as spec
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.search_state import raw_state
from research_runtime.cached_observation_state import CachedObservationMWState as MWState
from pointed_quartic_search import PointedQuarticSearch,sources as pointed_sources
geometry=SourceFileLoader('wide_r17_geometry',str(CAS/'prospective_half_lattice_v2.sage')).load_module()
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/compact-six-r17-h4096-v1'
CENSUS=ROOT/'artifacts/local/elliptic-curves/compact-r17-fresh-generic-census-v1'

def sources():
    paths=[Path(__file__).resolve(),Path(spec.__file__).resolve(),spec.ATLAS,Path(cert.__file__).resolve(),
           CAS/'research_runtime/quotient_only_reduction.py',CAS/'research_runtime/memory_store.py',
           CAS/'research_runtime/search_state.py',CAS/'prospective_half_lattice_v2.sage',CAS/'research_runtime/cached_observation_state.py']
    return {**pointed_sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}}

def prepare(directory):
    if (directory/'point-protocol.json').exists():raise FileExistsError('preserve point protocol')
    selection=cert.read(directory/'selection-ledger.json');replays=cert.read(CENSUS/'verification/ledger.json')
    if selection['status']!='COMPLETE_DECLARED_ATTEMPTS' or len(selection['rows'])!=6 or any(r['status']!='PASS' for r in selection['rows']):raise ArithmeticError('selection unfinished')
    if replays['status']!='COMPLETE_DECLARED_REPLAY_ATTEMPTS' or len(replays['rows'])!=6 or any(r['status']!='PASS' for r in replays['rows']):raise ArithmeticError('generic census replay unfinished')
    families=cert.read(directory/'protocol.json')['families'];masks={};hashes={};replay_hashes={}
    for family in families:
        census=cert.read(CENSUS/family/'generic-census.json')
        proof=next(r for r in replays['rows'] if r['family']==family);replay=cert.read(ROOT/proof['output'])
        if census['status']!='COMPLETE_DECLARED_CENSUS' or len(census['records'])!=131072 or replay['census_sha256']!=cert.hashed(CENSUS/family/'generic-census.json'):raise ArithmeticError('fresh generic proof differs')
        destination=directory/family/'generic-census.json'
        if destination.exists():raise FileExistsError('generic census copy already exists')
        shutil.copyfile(CENSUS/family/'generic-census.json',destination)
        hashes[family]=cert.hashed(destination);replay_hashes[family]=cert.hashed(ROOT/proof['output'])
        masks[family]=[r['mask'] for r in census['selected']]
    checkpoint(directory/'point-protocol.json',{'schema':'elliptic-curves.compact-six-r17-wide-point.v1','sources':sources(),
        'selection_protocol_sha256':cert.hashed(directory/'protocol.json'),'families':families,
        'population_hashes':{f:cert.hashed(directory/f/'population.json') for f in families},
        'generic_census_hashes':hashes,'generic_census_replay_hashes':replay_hashes,'generic_masks':masks,
        'charts_per_fibre':43,'chart_height':100000,'seconds_per_chart':4,'admission_prime_bound':251,
        'post_batch_admission_prime_bound':1000,'height_precision_bits':384,
        'coordinate_policy':{'kind':'metric','weight':'16'},'point_worker_wall_seconds':300,'point_worker_rss_bytes':1610612736,
        'maximum_point_workers':4,'fixed_addresses':24,'stop_rank_per_worker':32,
        'scope':'Exactly four fixed addresses per family,43 charts each, with no catalogue read in selection or execution. Every address receives a new attempt; repeated prior addresses are disclosed after the batch, with no refill. Use exact membership and quotient-only caches,384-bit numerical heights, exact retained-chart/rank replay and complete-cloud admission through997. Bounded failures prove no upper rank.'})
    print('FROZEN SIX-FAMILY R17 H4096 POINT BATCH',flush=True)

def bindings(directory):
    p=cert.read(directory/'point-protocol.json')
    if p['sources']!=sources() or p['selection_protocol_sha256']!=cert.hashed(directory/'protocol.json'):raise ArithmeticError('point sources changed')
    return p

def run(directory,family,index):
    protocol=bindings(directory)
    for name,key in [('population.json','population_hashes'),('generic-census.json','generic_census_hashes')]:
        if cert.hashed(directory/family/name)!=protocol[key][family]:raise ArithmeticError('frozen family input changed')
    t=cert.read(directory/family/'population.json')['finalists'][index]['parameter'];f=next(r for r in cert.read(spec.ATLAS)['families'] if r['family']==family)
    folder=directory/family/f'candidate-{index:02}';output=folder/'result.json'
    if output.exists() and cert.read(output)['status']!='RUNNING':raise FileExistsError('terminal point checkpoint exists')
    original,generic=spec.specialize(f,t)
    curve=EllipticCurve(QQ,[QQ(str(q)) for q in original]);pts=[curve([QQ(str(x)),QQ(str(y))]) for x,y in generic]
    minimal=curve.local_data(2).minimal_model();iso=curve.isomorphism_to(minimal)
    a1,a2,a3,a4,a6=minimal.a_invariants();c2,c4,c6=a2+a1*a1/4,a4+a1*a3/2,a6+a3*a3/4
    model=tuple(F(str(v)) for v in (0,0,0,c4-c2*c2/3,c6-c2*c4/3+2*c2**3/27));points=[]
    for P in pts:
        x,y=iso(P).xy();points.append((F(str(x+c2/3)),F(str(y+(a1*x+a3)/2))))
    # Both endpoints are short models, so the exact composite transport is a scale.
    u=F(str(iso.u)) if hasattr(iso,'u') else F(str(iso.tuple()[0]))
    if tuple((x/u**2,y/u**3) for x,y in generic)!=tuple(points):raise ArithmeticError('generic point transport failed')
    if model[3]!=original[3]/u**4 or model[4]!=original[4]/u**6:raise ArithmeticError('model transport failed')
    cache=QuotientOnlyReductionCache(MemoryFactStore());initial=raw_state(model,points,cache=cache,prime_bound=1000)
    initial=MWState.from_record(initial.record(),cache=cache)
    if output.exists():
        result=cert.read(output)
        if result['protocol_hash']!=digest(protocol) or result['initial_state']!=initial.record():raise ArithmeticError('resume mismatch')
        cache.store.import_snapshot(result['arithmetic_facts']);state=MWState.from_record(result['final_state'],cache=cache)
    else:
        state=initial;result={'schema':'elliptic-curves.compact-r17-wide-point-search.v1','family':family,'parameter':t,'protocol_hash':digest(protocol),
            'curve':list(map(str,model)),'generic_points':[list(map(str,p)) for p in points],'family_to_curve_scale_u':str(u),
            'initial_state':initial.record(),'final_state':initial.record(),'arithmetic_facts':cache.store.snapshot(),
            'centres':[],'charts':[],'status':'RUNNING','rank_lower_bound':initial.rank}
        checkpoint(output,result)
    if initial.rank!=17:
        result['status']='INCOMPLETE_GENERIC_MOD2_CERTIFICATE';checkpoint(output,result)
        print('R17 GENERIC PROOF INCOMPLETE',family,t,initial.rank,'of17; all inputs retained',flush=True);return
    if not result['centres']:
        gram,asymmetry=geometry.canonical_height_gram(model,points);rounded=geometry.rounded_gram(gram,1000000);oracle=geometry.CosetOracle(rounded);centres=[]
        for mask in protocol['generic_masks'][family]:
            norm,rep,error=oracle.solve(tuple((mask>>j)&1 for j in range(17)))
            if len(rep)!=17 or any((rep[j]-(mask>>j))%2 for j in range(17)):raise ArithmeticError('specialized parity mismatch')
            centres.append({'mask':mask,'representative':list(rep),'metric_norm':norm,'cvp_error':error})
        centres.sort(key=lambda r:(-r['metric_norm'],r['mask']))
        result.update(metric_gram=[[str(x) for x in row] for row in gram],maximum_gram_asymmetry=str(asymmetry),centres=centres);checkpoint(output,result)
    primes=tuple(map(int,prime_range(3,252)))
    for i in range(len(result['charts']),43):
        centre=result['centres'][i];rep=centre['representative']+[0]*(state.rank-17)
        search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=protocol['coordinate_policy'])
        outcome=search.search(100000,4,checkpoint_dir=folder/'charts'/state.key)
        for point in outcome.curve_points:state=state.adjoin(point,cache=cache,extra_primes=primes)
        result['charts'].append({'centre':centre,'search':outcome.record,'rank_lower_bound':state.rank,'admission_prime_bound':251})
        result.update(final_state=state.record(),arithmetic_facts=cache.store.snapshot(),rank_lower_bound=state.rank);checkpoint(output,result)
        print('PROSPECTIVE R17 H4096',family,t,'chart',i+1,'rank',state.rank,flush=True)
        if state.rank>=32:result['status']='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY';checkpoint(output,result);return
    result['status']='COMPLETE_DECLARED_PILOT';checkpoint(output,result)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run']);p.add_argument('--directory',type=Path,default=DIRECTORY);p.add_argument('--family');p.add_argument('--index',type=int,default=0);a=p.parse_args()
    prepare(a.directory) if a.stage=='prepare' else run(a.directory,a.family,a.index)
