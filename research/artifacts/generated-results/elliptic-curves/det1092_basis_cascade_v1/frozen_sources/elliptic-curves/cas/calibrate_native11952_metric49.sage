#!/usr/bin/env sage-python
"""Generic-point-only native11952 calibration of the current49-chart metric policy."""
import argparse
from importlib.machinery import SourceFileLoader
from pathlib import Path
import sys
from sage.all import prime_range
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.cached_observation_state import CachedObservationMWState as MWState
from research_runtime.search_state import raw_state
from research_runtime.pointed_orbit_compression import compress
from pointed_quartic_search import PointedQuarticSearch,sources as pointed_sources
geometry=SourceFileLoader('native11952_metric_geometry',str(CAS/'prospective_half_lattice_v2.sage')).load_module()
D=ROOT/'artifacts/local/elliptic-curves/native11952-metric49-control-v1'
BLIND=ROOT/'elliptic-curves/data/half_lattice_rank29_control_inputs_v1.json'
CENSUS=ROOT/'artifacts/local/elliptic-curves/compact-six-r17-h4096-v1/11952/generic-census.json'
RADIUS=ROOT/'artifacts/generated-results/elliptic-curves/r17_exact_parity_radius_11952_v1.json'


def sources():
    return {**pointed_sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in (
        Path(__file__).resolve(),Path(cert.__file__).resolve(),CAS/'prospective_half_lattice_v2.sage',
        CAS/'research_runtime/memory_store.py',CAS/'research_runtime/quotient_only_reduction.py',
        CAS/'research_runtime/cached_observation_state.py',CAS/'research_runtime/pointed_orbit_compression.py',CAS/'alternate_quartic_covers.py')}}


def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve native-control protocol')
    boundary=cert.read(BLIND)
    if boundary['boundary']['output_contains_exceptional_point_coordinates']:raise ArithmeticError('exceptional points in input')
    record=next(r for r in boundary['cases'] if r['label']=='curve12-2024-rank29')
    census=cert.read(CENSUS);radius=cert.read(RADIUS)
    if record['generic_height_gram']!=census['gram'] or radius['status']!='PASS' or radius['classes_with_exact_minimum12']!=49:raise ArithmeticError('native generic geometry gate failed')
    initial=[r['mask'] for r in census['selected']];all_masks=[r['mask'] for r in census['records'] if cert.F(r['norm'])==12]
    extra=[m for m in all_masks if m not in initial]
    if len(initial)!=43 or len(extra)!=6:raise ArithmeticError('fixed43-plus6 roster changed')
    original=tuple(map(cert.F,record['short_model']));model=(cert.F(0),cert.F(0),cert.F(0),original[3]/6**4,original[4]/6**6)
    points=[(cert.F(x)/6**2,cert.F(y)/6**3) for x,y in record['generic_points']]
    if len(points)!=17 or any(not cert.is_on_weierstrass_curve(model,p) for p in points):raise ArithmeticError('generic-only coordinate transport failed')
    checkpoint(D/'redacted-input.json',{'curve':list(map(str,model)),'generic_points':[list(map(str,p)) for p in points],
        'generic_height_gram':record['generic_height_gram'],'boundary_input_sha256':cert.hashed(BLIND),
        'boundary_to_curve_scale_u':'6','boundary_short_model':record['short_model'],'boundary_generic_points':record['generic_points'],
        'claim_boundary':'Only one curve and its17 reconstructed generic points. No exceptional point, displayed-rank field or target parameter is included.'})
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.native11952-current-metric-calibration.v1','sources':sources(),
        'input_sha256':cert.hashed(D/'redacted-input.json'),'census_sha256':cert.hashed(CENSUS),'radius_proof_sha256':cert.hashed(RADIUS),
        'initial43_masks':initial,'additional6_masks':extra,'height':100000,'seconds_per_chart':4,'admission_prime_bound':997,
        'coordinate_policy':{'kind':'metric','weight':'16'},'worker_wall_seconds':400,'worker_rss_bytes':1610612736,'maximum_workers':1,
        'gate':'A recorded native11952 generic-point-only control is available, and its displayed generic Gram exactly matches the newly proved49-class maximum stratum. Test the current metric/GMP implementation with the original43 classes first and six additional classes last. This is a calibration of conditional visibility, separate from prospective candidate selection.',
        'selection':'Compute all representatives before searching from the17-point384-bit numerical height Gram; order each phase by computed norm then mask. All49 charts are attempted without a target-rank stop.',
        'claim_boundary':'A known-curve generic-point-only control, not a new curve or rank-record result. All public exceptional points remain outside the search. Incomplete boxes or low recovery prove no upper rank. Comparison with older implementations is diagnostic, not a controlled causal claim about the class cap alone.'})
    print('FROZEN NATIVE11952 METRIC43 PLUS6 CONTROL',flush=True)


def run():
    protocol=cert.read(D/'protocol.json');p=D/'redacted-input.json'
    if protocol['sources']!=sources() or cert.hashed(p)!=protocol['input_sha256']:raise ArithmeticError('native control binding changed')
    source=cert.read(p);model=tuple(map(cert.F,source['curve']));points=tuple(tuple(map(cert.F,p)) for p in source['generic_points'])
    output=D/'candidate-00/result.json';folder=output.parent
    if output.exists():raise FileExistsError('preserve control attempt')
    cache=ReductionCache(MemoryFactStore());state=raw_state(model,points,cache=cache,prime_bound=1000);state=MWState.from_record(state.record(),cache=cache)
    if state.rank!=17:raise ArithmeticError('native generic independence failed')
    gram,asymmetry=geometry.canonical_height_gram(model,points);oracle=geometry.CosetOracle(geometry.rounded_gram(gram,1000000));centres=[]
    for label,key in [('initial43','initial43_masks'),('additional6','additional6_masks')]:
        block=[]
        for mask in protocol[key]:
            norm,rep,error=oracle.solve(tuple((mask>>j)&1 for j in range(17)))
            block.append({'phase':label,'mask':mask,'representative':list(rep),'metric_norm':norm,'cvp_error':error})
        centres.extend(sorted(block,key=lambda r:(-r['metric_norm'],r['mask'])))
    data={'schema':'elliptic-curves.native11952-metric-control-result.v1','protocol_hash':digest(protocol),'family':'native11952-control','parameter':'generic-only-fixture',
        'curve':source['curve'],'generic_points':source['generic_points'],'initial_state':state.record(),'initial_dimension':17,
        'metric_gram':[[str(x) for x in row] for row in gram],'centres':centres,'charts':[],'status':'RUNNING','rank_lower_bound':17,'final_state':state.record(),'arithmetic_facts':cache.store.snapshot()}
    checkpoint(output,data);primes=tuple(map(int,prime_range(3,998)))
    for c in centres:
        rep=c['representative']+[0]*(state.rank-17)
        search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=protocol['coordinate_policy'])
        outcome=search.search(100000,4,checkpoint_dir=folder/'charts'/state.key);compression=compress(model,state.basis,rep,outcome.curve_points)
        for j in compression['kept_indices']:
            state=state.adjoin(outcome.curve_points[j],cache=cache,extra_primes=primes)
            if not isinstance(state,MWState):state=MWState.from_record(state.record(),cache=cache)
        data['charts'].append({'centre':c,'search':outcome.record,'admission_compression':compression,'rank_lower_bound':state.rank,'admission_prime_bound':997})
        data.update(rank_lower_bound=state.rank,final_state=state.record(),arithmetic_facts=cache.store.snapshot());checkpoint(output,data)
        print('NATIVE11952 CURRENT CONTROL',len(data['charts']),c['phase'],'rank',state.rank,flush=True)
    data['status']='COMPLETE_DECLARED_PILOT';checkpoint(output,data)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run']);a=p.parse_args();prepare() if a.stage=='prepare' else run()
