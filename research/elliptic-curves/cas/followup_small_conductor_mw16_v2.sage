#!/usr/bin/env sage-python
"""Frozen301-chart attempt to raise the new small-conductor curve from rank22 to23."""
import argparse
from importlib.machinery import SourceFileLoader
from pathlib import Path
import sys
from sage.all import prime_range

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path[:0]=[str(CAS),str(ROOT/'elliptic-curves')]
import certify_compact_r17_candidates as cert
import compact_mw16_specialization as spec
from research_runtime.store import checkpoint,digest
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.mw_state import MWState
from pointed_quartic_search import PointedQuarticSearch,sources as pointed_sources
geometry=SourceFileLoader('atlas_adaptive_general_geometry',str(CAS/'prospective_half_lattice_v2.sage')).load_module()
PARENT=ROOT/'artifacts/local/elliptic-curves/prospective-mw16-f5-next12-v1'
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/prospective-mw16-small-conductor-followup-v2'


def sources():
    paths=(Path(__file__).resolve(),Path(spec.__file__).resolve(),spec.ATLAS,Path(cert.__file__).resolve(),
           CAS/'prospective_half_lattice_v2.sage',CAS/'research_runtime/memory_store.py',CAS/'research_runtime/quotient_only_reduction.py')
    return {**pointed_sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}}


def prepare(directory):
    if (directory/'protocol.json').exists():raise FileExistsError('follow-up already frozen')
    path=PARENT/'a1-fibration-05/candidate-00/result.json';result=cert.read(path)
    certificate=ROOT/'artifacts/generated-results/elliptic-curves/prospective_mw16_next12_results_v1.json'
    record=next(r for r in cert.read(certificate)['curves'] if r['family']==result['family'] and r['parameter']==result['parameter'])
    if result['status']!='COMPLETE_DECLARED_PILOT' or result['rank_lower_bound']!=22 or record['discovery_witness']['sha256']!=cert.hashed(path):raise ArithmeticError('rank22 input changed')
    if record['icarm_matches'] or record['previous_matches']:raise ArithmeticError('not the new prospective rank22 curve')
    model=tuple(map(cert.F,result['curve']));points=tuple(tuple(map(cert.F,p)) for p in result['final_state']['state']['reductions']['points'])
    proof=record['rank_certificate'];cert.checked_rank(model,points,[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
    if result['curve']!=record['curve'] or [list(map(str,p)) for p in points]!=record['points']:raise ArithmeticError('certified basis changed')
    census_path=PARENT/result['family']/'generic-census.json';census=cert.read(census_path)
    if census['status']!='COMPLETE_DECLARED_CENSUS':raise ArithmeticError('generic census incomplete')
    selected=sorted(census['records'][1:],key=lambda r:(-cert.F(r['norm']),r['mask']))[:301]
    if selected[:43]!=census['selected']:raise ArithmeticError('initial class prefix changed')
    pool=[{'mask':r['mask'],'generic_norm':r['norm']} for r in selected]
    checkpoint(directory/'protocol.json',{'schema':'elliptic-curves.prospective-mw16-rank22-followup.v1','sources':sources(),
        'input_path':str(path.relative_to(ROOT)),'input_sha256':cert.hashed(path),'certificate_sha256':cert.hashed(certificate),
        'generic_census_path':str(census_path.relative_to(ROOT)),'generic_census_sha256':cert.hashed(census_path),
        'family':result['family'],'parameter':result['parameter'],'initial_rank_lower_bound':22,'generic_dimension':16,
        'selection':'The new3/17 rank22 curve has proven exact conductor below all but two listed rank>=22 curves; one further independent point would beat every listed rank>=23 conductor. Use only its certified discovered points.',
        'generic_pool':pool,'quotient_rule':'Top301 computed generic parity norms, mask tie order; cyclic nonzero6-bit words ordered by Hamming weight then integer, crossed with those generic classes; order by specialized numerical CVP norm.',
        'charts':301,'chart_height':100000,'seconds_per_chart':4,'admission_prime_bound':251,
        'coordinate_policy':{'kind':'metric','weight':'16'},'workers':1,'worker_wall_seconds':1800,'worker_rss_bytes':1610612736,
        'target_rank_lower_bound':23,'scope':'One additional adaptive follow-up after the entire initial batch and its exact point-certificate replay. No public exceptional points. Bounded coverage is not an upper rank bound.'})
    print('FROZEN MW16 RANK22 FOLLOWUP',result['family'],result['parameter'],flush=True)


def run(directory):
    protocol=cert.read(directory/'protocol.json');path=ROOT/protocol['input_path'];source=cert.read(path)
    if protocol['sources']!=sources() or cert.hashed(path)!=protocol['input_sha256']:raise ArithmeticError('follow-up source changed')
    folder=directory/'candidate-00';output=folder/'result.json'
    model=tuple(map(cert.F,source['curve']));cache=ReductionCache(MemoryFactStore());cache.store.import_snapshot(source['arithmetic_facts'])
    initial=MWState.from_record(source['final_state'],cache=cache);dimension=initial.rank
    if dimension!=22:raise ArithmeticError('initial rank changed')
    if output.exists():
        result=cert.read(output)
        if result['protocol_hash']!=digest(protocol) or result['initial_state']!=initial.record():raise ArithmeticError('resume mismatch')
        if result['status']!='RUNNING':print('RETAINED FOLLOWUP',result['status'],flush=True);return
        cache.store.import_snapshot(result['arithmetic_facts']);state=MWState.from_record(result['final_state'],cache=cache)
    else:
        points=tuple(tuple(map(cert.F,p)) for p in initial.basis);gram,asymmetry=geometry.canonical_height_gram(model,points)
        oracle=geometry.CosetOracle(geometry.rounded_gram(gram,1000000))
        words=sorted(range(1,1<<(dimension-16)),key=lambda w:(w.bit_count(),w));centres=[]
        for i,g in enumerate(protocol['generic_pool']):
            word=words[i%len(words)];mask=g['mask']|(word<<16);residue=[(mask>>j)&1 for j in range(dimension)]
            norm,rep,error=oracle.solve(residue)
            if len(rep)!=dimension or any((rep[j]-residue[j])%2 for j in range(dimension)):raise ArithmeticError('adaptive centre dimension or parity mismatch')
            centres.append({'generic_mask':g['mask'],'quotient_word':word,'parity':mask,'representative':list(rep),'metric_norm':norm})
        centres.sort(key=lambda r:(-r['metric_norm'],r['generic_mask'],r['quotient_word']));state=initial
        result={'schema':'elliptic-curves.prospective-mw16-adaptive-result.v1','family':protocol['family'],'parameter':protocol['parameter'],
            'curve':source['curve'],'protocol_hash':digest(protocol),'initial_dimension':dimension,'generic_points':source['generic_points'],'family_to_curve_scale_u':source['family_to_curve_scale_u'],
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
        print('MW16 FOLLOWUP',i+1,'rank',state.rank,flush=True)
        if state.rank>=23:result['status']='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY';checkpoint(output,result);return
    result['status']='COMPLETE_DECLARED_PILOT';checkpoint(output,result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run']);p.add_argument('--directory',type=Path,default=DIRECTORY);a=p.parse_args()
    prepare(a.directory) if a.stage=='prepare' else run(a.directory)
