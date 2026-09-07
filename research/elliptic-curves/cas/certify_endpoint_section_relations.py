#!/usr/bin/env python3
"""Exact rational span of specialized generic sections, with bounded proposals."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
from memory_rank_certificate import checked_rank
from half_lattice_pointed_sieve import linear_combination_python
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/endpoint-section-relations-v1';INPUT=ART/'compact_atlas_endpoints_v2.json';OUT=ART/'endpoint_section_relations_v1.json'

def sources():
    names=['certify_endpoint_section_relations.py','propose_endpoint_section_relations.sage','prospective_half_lattice_v2.sage','half_lattice_pointed_sieve.py','memory_rank_certificate.py','certify_compact_r17_candidates.py']
    return {**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in names},str(INPUT.relative_to(ROOT)):cert.hashed(INPUT)}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve exact endpoint relation protocol')
    rows=[{'family':r['family'],'endpoint':r['endpoint'],'rank':r['rank_lower_bound']} for r in cert.read(INPUT)['rows'] if r['status']=='CERTIFIED_SPECIALIZED_SUBGROUP']
    if len(rows)!=21:raise ArithmeticError('all21 nonsingular endpoints required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.endpoint-section-relations.v1','sources':sources(),'rows':rows,'maximum_denominator':64,'maximum_relation_coefficient':64,'precision_bits':384,'wall_seconds':600,'rss_bytes':1610612736,'gate':'The full specialized section sets give the same11to17 lower bounds modulo2,3,5, below generic16/17 for most endpoints. Finite tests alone do not establish dependence. Use numerical heights only to propose small rational coordinates, then verify exact integer group relations against each independently certified subset.','scope':'A verified relation dQ=sum c_i P_i places Q in the rational span of the certified subset. If all generic sections have such relations, its rational span has exactly the certified dimension; this is not the whole elliptic curve rank or a saturation claim. Failed rational approximations remain UNKNOWN. All21 endpoints and coefficient limits are fixed; no searches or new points are generated.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen exact relation sources differ')
    return p

def expected():
    p=protocol();proposal=cert.read(D/'proposals.json');rows=[]
    if proposal['status']!='COMPLETE_DECLARED_PROPOSALS' or len(proposal['rows'])!=21:raise ArithmeticError('complete finite proposals required')
    original=[r for r in cert.read(INPUT)['rows'] if r['status']=='CERTIFIED_SPECIALIZED_SUBGROUP']
    for r,q in zip(original,proposal['rows']):
        if (q['family'],q['endpoint'])!=(r['family'],r['endpoint']) or len(q['relations'])!=len(r['generic_points']):raise ArithmeticError('all section relation roster differs')
        model=tuple(map(cert.F,r['curve']));basis=[tuple(map(cert.F,P)) for P in r['points']];proof=r['rank_certificate']
        checked_rank(model,basis,[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
        verified=[]
        for j,rel in enumerate(q['relations']):
            if rel['section_index']!=j:raise ArithmeticError('section order differs')
            if rel['status']=='UNKNOWN':verified.append(rel);continue
            d=rel['denominator'];word=rel['coefficients']
            if type(d)!=int or not 1<=d<=p['maximum_denominator'] or len(word)!=len(basis) or any(type(c)!=int or abs(c)>p['maximum_relation_coefficient'] for c in word):raise ArithmeticError('bounded integral relation differs')
            Q=tuple(map(cert.F,r['generic_points'][j]))
            if linear_combination_python(model,[*basis,Q],[*word,-d]) is not None:raise ArithmeticError('proposed relation is not exact')
            verified.append({'section_index':j,'status':'EXACT_INTEGER_GROUP_RELATION','denominator':d,'coefficients':word})
        complete=all(s['status']=='EXACT_INTEGER_GROUP_RELATION' for s in verified)
        rows.append({'family':r['family'],'endpoint':r['endpoint'],'curve':r['curve'],'independent_points':r['points'],'generic_points':r['generic_points'],'rank_certificate':proof,'certified_lower_bound':len(basis),'generic_section_span_rank':len(basis) if complete else 'UNKNOWN','relations':verified})
    return {'schema':'elliptic-curves.endpoint-section-relations-result.v1','status':'PASS','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'proposal_sha256':cert.hashed(D/'proposals.json'),'rows':rows,'completed_span_ranks':sum(r['generic_section_span_rank']!='UNKNOWN' for r in rows),'claim_boundary':p['scope']}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','build','check']);a=p.parse_args()
    if a.stage=='prepare':prepare()
    else:
        d=expected()
        if a.stage=='check':
            if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('exact endpoint section span proof differs')
        else:
            if OUT.exists():raise FileExistsError('preserve endpoint relation proof')
            checkpoint(OUT,d)
        print('EXACT ENDPOINT SECTION SPANS',d['completed_span_ranks'],'of21',flush=True)
