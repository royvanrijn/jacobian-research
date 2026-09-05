#!/usr/bin/env sage-python
"""Fixed next12 addresses from the successful compact-MW16 family05 population."""
import argparse
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
from pathlib import Path
import shutil
import sys
from sage.all import EllipticCurve, QQ, prime_range
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import compact_mw16_specialization as spec
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.search_state import raw_state
from research_runtime.mw_state import MWState
from pointed_quartic_search import PointedQuarticSearch,sources as pointed_sources
geometry=SourceFileLoader('next12_geometry',str(CAS/'prospective_half_lattice_v2.sage')).load_module()
PARENT=ROOT/'artifacts/local/elliptic-curves/prospective-mw16-h1024-v1'
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/prospective-mw16-f5-next12-v1'
FAMILY='a1-fibration-05'

def sources():
    paths=[Path(__file__).resolve(),Path(spec.__file__).resolve(),spec.ATLAS,Path(cert.__file__).resolve(),
           CAS/'research_runtime/quotient_only_reduction.py',CAS/'research_runtime/memory_store.py',
           CAS/'research_runtime/search_state.py',CAS/'prospective_half_lattice_v2.sage']
    return {**pointed_sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}}

def prepare(directory):
    if (directory/'point-protocol.json').exists():raise FileExistsError('point protocol already frozen')
    ledger=cert.read(PARENT/'point-ledger.json')
    if ledger['status']!='COMPLETE_FIXED_BATCH_ATTEMPTS':raise ArithmeticError('parent batch unfinished')
    measured=[r for r in ledger['rows'] if r['family']==FAMILY]
    if len(measured)!=4 or sorted(r['rank_lower_bound'] for r in measured)!=[24,24,24,25]:raise ArithmeticError('observed family-selection gate differs')
    population=cert.read(PARENT/FAMILY/'population.json');census=cert.read(PARENT/FAMILY/'generic-census.json')
    if census['status']!='COMPLETE_DECLARED_CENSUS' or len(census['records'])!=65536:raise ArithmeticError('generic census incomplete')
    selected=population['retained_candidates'][4:16]
    if len(selected)!=12 or population['retained_candidates'][:4]!=population['finalists']:raise ArithmeticError('parent ranking differs')
    folder=directory/FAMILY;folder.mkdir(parents=True)
    shutil.copyfile(PARENT/FAMILY/'generic-census.json',folder/'generic-census.json')
    checkpoint(folder/'population.json',{'family':FAMILY,'finalists':selected,'parent_population_sha256':cert.hashed(PARENT/FAMILY/'population.json'),'parent_rank_positions':list(range(5,17))})
    checkpoint(directory/'protocol.json',{'schema':'elliptic-curves.prospective-mw16-next12-selection.v1',
        'parent_population_sha256':cert.hashed(PARENT/FAMILY/'population.json'),'parent_ledger_sha256':cert.hashed(PARENT/'point-ledger.json'),
        'family':FAMILY,'fixed_addresses':selected,'population_count':population['candidate_count'],
        'gate':'The prior four fixed family05 addresses all certified lower bounds24–25; select ranks5–16 from its already frozen full-prime height1024 population.',
        'claim_boundary':'An adaptive choice of family based on observed point yield. No new parameter enumeration, refitting of scores, catalogue prefilter, or replacement addresses. Every fixed address receives an attempt.'})
    checkpoint(directory/'point-protocol.json',{'schema':'elliptic-curves.prospective-mw16-next12-point.v1','sources':sources(),
        'selection_protocol_sha256':cert.hashed(directory/'protocol.json'),'families':[FAMILY],
        'population_hashes':{FAMILY:cert.hashed(folder/'population.json')},'generic_census_hashes':{FAMILY:cert.hashed(folder/'generic-census.json')},
        'generic_masks':{FAMILY:[r['mask'] for r in census['selected']]},
        'charts_per_fibre':43,'chart_height':100000,'seconds_per_chart':4,'admission_prime_bound':251,
        'post_batch_admission_prime_bound':1000,'height_precision_bits':384,
        'coordinate_policy':{'kind':'metric','weight':'16'},'point_worker_wall_seconds':300,'point_worker_rss_bytes':1610612736,
        'maximum_point_workers':4,'fixed_addresses':12,'stop_rank_per_worker':32,
        'scope':'Exactly12 additional fixed addresses,43 charts each, no catalogue read in selection or execution. Use the validated quotient-only cache and explicit384-bit PARI heights; exact point/rank replay and complete-cloud admission through997 follow terminal batch attempts. Bounded misses are not upper bounds.'})
    print('FROZEN MW16 FAMILY05 NEXT12',[r['parameter'] for r in selected],flush=True)

def bindings(directory):
    p=cert.read(directory/'point-protocol.json')
    if p['sources']!=sources() or p['selection_protocol_sha256']!=cert.hashed(directory/'protocol.json'):raise ArithmeticError('point sources changed')
    return p

def run(directory,family,index):
    protocol=bindings(directory)
    for name,key in [('population.json','population_hashes'),('generic-census.json','generic_census_hashes')]:
        if cert.hashed(directory/family/name)!=protocol[key][family]:raise ArithmeticError('frozen family input changed')
    t=cert.read(directory/family/'population.json')['finalists'][index]['parameter'];f=next(r for r in cert.read(spec.ATLAS)['families'] if r['fibration_id']==family)
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
    if output.exists():
        result=cert.read(output)
        if result['protocol_hash']!=digest(protocol) or result['initial_state']!=initial.record():raise ArithmeticError('resume mismatch')
        cache.store.import_snapshot(result['arithmetic_facts']);state=MWState.from_record(result['final_state'],cache=cache)
    else:
        state=initial;result={'schema':'elliptic-curves.prospective-mw16-point-search.v1','family':family,'parameter':t,'protocol_hash':digest(protocol),
            'curve':list(map(str,model)),'generic_points':[list(map(str,p)) for p in points],'family_to_curve_scale_u':str(u),
            'initial_state':initial.record(),'final_state':initial.record(),'arithmetic_facts':cache.store.snapshot(),
            'centres':[],'charts':[],'status':'RUNNING','rank_lower_bound':initial.rank}
        checkpoint(output,result)
    if initial.rank!=16:
        result['status']='INCOMPLETE_GENERIC_MOD2_CERTIFICATE';checkpoint(output,result)
        print('MW16 GENERIC PROOF INCOMPLETE',family,t,initial.rank,'of16; all inputs retained',flush=True);return
    if not result['centres']:
        gram,asymmetry=geometry.canonical_height_gram(model,points);rounded=geometry.rounded_gram(gram,1000000);oracle=geometry.CosetOracle(rounded);centres=[]
        for mask in protocol['generic_masks'][family]:
            norm,rep,error=oracle.solve(tuple((mask>>j)&1 for j in range(16)))
            if len(rep)!=16 or any((rep[j]-(mask>>j))%2 for j in range(16)):raise ArithmeticError('specialized parity mismatch')
            centres.append({'mask':mask,'representative':list(rep),'metric_norm':norm,'cvp_error':error})
        centres.sort(key=lambda r:(-r['metric_norm'],r['mask']))
        result.update(metric_gram=[[str(x) for x in row] for row in gram],maximum_gram_asymmetry=str(asymmetry),centres=centres);checkpoint(output,result)
    primes=tuple(map(int,prime_range(3,252)))
    for i in range(len(result['charts']),43):
        centre=result['centres'][i];rep=centre['representative']+[0]*(state.rank-16)
        search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=protocol['coordinate_policy'])
        outcome=search.search(100000,4,checkpoint_dir=folder/'charts'/state.key)
        for point in outcome.curve_points:state=state.adjoin(point,cache=cache,extra_primes=primes)
        result['charts'].append({'centre':centre,'search':outcome.record,'rank_lower_bound':state.rank,'admission_prime_bound':251})
        result.update(final_state=state.record(),arithmetic_facts=cache.store.snapshot(),rank_lower_bound=state.rank);checkpoint(output,result)
        print('PROSPECTIVE MW16',family,t,'chart',i+1,'rank',state.rank,flush=True)
        if state.rank>=32:result['status']='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY';checkpoint(output,result);return
    result['status']='COMPLETE_DECLARED_PILOT';checkpoint(output,result)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run']);p.add_argument('--directory',type=Path,default=DIRECTORY);p.add_argument('--family',default=FAMILY);p.add_argument('--index',type=int,default=0);a=p.parse_args()
    prepare(a.directory) if a.stage=='prepare' else run(a.directory,a.family,a.index)
