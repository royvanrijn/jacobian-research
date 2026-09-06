#!/usr/bin/env python3
"""Finite rank certificates for already constructed marked-carrier fibres."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from research_runtime.store import checkpoint,digest
from research_runtime.search_state import raw_state
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from memory_rank_certificate import checked_rank
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/native-rank3-carrier-subgroups-v1'
INPUT=ART/'native_rank3_carrier_images_v1.json';REPLAY=ART/'native_rank3_carrier_images_replay_v1.json';COMPARISON=ART/'native_rank3_carrier_images_comparison_v1.json';OUT=ART/'native_rank3_carrier_subgroups_v1.json'
def sources():
    paths=[Path(__file__).resolve(),INPUT,REPLAY,COMPARISON,Path(cert.__file__),Path(spec.__file__),spec.ATLAS,ROOT/'elliptic-curves/cas/memory_rank_certificate.py',ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',ROOT/'elliptic-curves/cas/research_runtime/search_state.py',ROOT/'elliptic-curves/cas/research_runtime/finite_reduction.py',ROOT/'elliptic-curves/cas/research_runtime/memory_store.py',ROOT/'elliptic-curves/cas/research_runtime/store.py']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve subgroup proof protocol')
    replay=cert.read(REPLAY);comparison=cert.read(COMPARISON);image=cert.read(INPUT)
    if replay['status']!='PASS' or comparison['status']!='PASS' or len(image['rows'])!=12:raise ArithmeticError('all twelve exact constructions required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.native-rank3-carrier-subgroups.v1','sources':sources(),'rows':[{'word':r['word'],'compact_parameter':r['compact_parameter']} for r in image['rows']],'prime_bound':997,'seconds_per_stage':180,'rss_bytes':2147483648,'scope':'Certify the span of17 exactly specialized generic points and the two already constructed native points on every fixed image. This is finite arithmetic on existing witnesses; the earlier400-bit gate governed prospective point-search cost. Check every point identity and exact finite quotient proof, preserving any lower bound below19. No additional point search, parameter, auxiliary group enumeration, saturation, original exact rank or rank-density claim.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('subgroup sources changed')
    return p

def one(row,p):
    f=next(f for f in cert.read(spec.ATLAS)['families'] if f['family']=='08234');model,generic=spec.specialize(f,row['compact_parameter']);extra=tuple(tuple(map(cert.F,P)) for P in row['supplied_points'])
    if list(map(str,model))!=row['curve'] or len(generic)!=17 or len(extra)!=2 or any(not cert.is_on_weierstrass_curve(model,P) for P in extra):raise ArithmeticError('complete inherited and native point witnesses required')
    state=raw_state(model,tuple(generic)+extra,cache=ReductionCache(MemoryFactStore()),prime_bound=p['prime_bound']);proof=checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime)
    return {'word':row['word'],'family':'08234','parameter':row['compact_parameter'],'curve':row['curve'],'generic_points':[[str(v) for v in P] for P in generic],'supplied_points':row['supplied_points'],'independent_points':[[str(v) for v in P] for P in state.basis],'rank_lower_bound':state.rank,'rank_certificate':proof,'model_coefficient_bits':row['model_coefficient_bits']}
def build():
    p=protocol()
    if OUT.exists() or (D/'checkpoint.json').exists():raise FileExistsError('preserve subgroup run')
    d={'status':'RUNNING','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'rows':[]};checkpoint(D/'checkpoint.json',d)
    for r in cert.read(INPUT)['rows']:
        d['rows'].append(one(r,p));checkpoint(D/'checkpoint.json',d);print('CONSTRUCTED CARRIER SUBGROUP',r['word'],'rank >=',d['rows'][-1]['rank_lower_bound'],flush=True)
    d.update(schema='elliptic-curves.native-rank3-carrier-subgroups-result.v1',status='PASS',claim_boundary=p['scope']);checkpoint(OUT,d)
def check():
    p=protocol();d=cert.read(OUT)
    if d['status']!='PASS' or d['sources']!=sources() or d['protocol_sha256']!=cert.hashed(D/'protocol.json') or len(d['rows'])!=12:raise ArithmeticError('complete subgroup proof binding differs')
    expected=[one(r,p) for r in cert.read(INPUT)['rows']]
    if json.loads(json.dumps(expected))!=d['rows']:raise ArithmeticError('finite subgroup replay differs')
    print('ALL12 CONSTRUCTED CARRIER SUBGROUP PROOFS REPLAY',flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['prepare','build','check']);v=a.parse_args();globals()[v.stage]()
