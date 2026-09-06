#!/usr/bin/env python3
"""Public-data reproduction and finite-quotient union with the earlier local27 subgroup."""
import argparse,json
from dataclasses import asdict
from pathlib import Path
import certify_compact_r17_candidates as cert
import memory_rank_certificate as memory
from mod2_reduction_independence import _primes_up_to
from audit_recorded_point_mod2_rank_v3 import insert
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/inventory200-current-catalogue-v1'
OUT=ART/'inventory188_public28_reproduction_v1.json'
INDEX=ART/'new_high_rank_curve_index_v20.json'
DB=D/'database.json'


def inputs():
    replay=cert.read(ART/'inventory200_current_catalogue_sage_replay_v1.json')
    if replay['status']!='PASS':raise ArithmeticError('independent catalogue match required')
    row=next(r for r in replay['comparisons'] if r['id']=='new-20260906-188')
    if row['q_isomorphism_matches']!=[619]:raise ArithmeticError('fixed public match differs')
    intake=cert.read(D/'metadata.json')
    if cert.hashed(DB)!=intake['sha256']:raise ArithmeticError('public snapshot differs')
    old=next(r for r in cert.read(INDEX)['curves'] if r['id']==row['id'])
    public=next(r for r in cert.read(DB)['curves'] if r['id']==619)
    if old['rank_lower_bound']!=27 or len(old['points'])!=27 or len(public['points'])!=28:
        raise ArithmeticError('fixed local27 and public28 inputs required')
    model=tuple(map(cert.F,old['curve']));source=tuple(map(cert.F,public['ainvs']))
    a=cert.weierstrass_invariants(source);b=cert.weierstrass_invariants(model)
    u=cert.square_root(a['c6']*b['c4']/(a['c4']*b['c6']))
    if u is None or a['c4']!=u**4*b['c4'] or a['c6']!=u**6*b['c6'] or any(model[:3]):
        raise ArithmeticError('exact rational short-model transport required')
    points=[]
    for raw in public['points']:
        x,y=map(cert.F,raw)
        if not cert.is_on_weierstrass_curve(source,(x,y)):raise ArithmeticError('public point off source equation')
        points.append(((x+a['b2']/12)/u**2,(y+(source[0]*x+source[2])/2)/u**3))
    if any(not cert.is_on_weierstrass_curve(model,p) for p in points):raise ArithmeticError('transported point off local equation')
    return old,public,model,points,u


def sources():
    files=[Path(__file__).resolve(),Path(cert.__file__),Path(memory.__file__),INDEX,DB,D/'metadata.json',
           ART/'inventory200_current_catalogue_sage_replay_v1.json',
           ROOT/'elliptic-curves/cas/audit_recorded_point_mod2_rank_v3.py',
           ROOT/'elliptic-curves/cas/research_runtime/finite_reduction.py',
           ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',
           ROOT/'elliptic-curves/cas/elliptic_candidate_record.py']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in files}


def build():
    if OUT.exists():raise FileExistsError('preserve public reproduction')
    old,public,model,transport,u=inputs();prior=[tuple(map(cert.F,p)) for p in old['points']]
    oldproof=old['rank_certificate'];tp=oldproof['no_rational_2_torsion_prime']
    if json.dumps(memory.checked_rank(model,prior,[r['prime'] for r in oldproof['signatures']],tp),sort_keys=True)!=json.dumps(oldproof,sort_keys=True):
        raise ArithmeticError('old27 certificate differs')
    points=prior+transport;cache=ReductionCache(MemoryFactStore());pivots={};signatures=[];processed=[]
    for prime in _primes_up_to(997):
        if prime==2:continue
        try:sig=cache.signature(model,points,prime)
        except ValueError:continue
        before=len(pivots)
        for row in sig.rows:insert(pivots,row)
        if len(pivots)>before:signatures.append(asdict(sig))
        processed.append(prime)
    indices=sorted(pivots);selected=[points[i] for i in indices];primes=[s['prime'] for s in signatures]
    union_proof=memory.checked_rank(model,selected,primes,tp)
    public_proof=memory.checked_rank(model,transport,primes,tp)
    result={'schema':'elliptic-curves.inventory188-public28.v1','status':'PASS','sources':sources(),
            'local_id':old['id'],'public_id':619,'public_url':'https://elliptic-rank.icarm.cloud/curve/619',
            'public_created_at':public['created_at'],'public_submitter':public['submitter'],
            'family':old['family'],'parameter':old['parameter'],'curve':old['curve'],
            'public_curve':public['ainvs'],'public_points':public['points'],'transport_scale_u':str(u),
            'transported_public_points':[list(map(str,p)) for p in transport],
            'public_rank_certificate':public_proof,'old_rank_lower_bound':27,
            'points':[list(map(str,p)) for p in points],'prime_bound':997,'processed_primes':processed,
            'union_signatures':signatures,'independent_column_indices':indices,
            'independent_points':[list(map(str,p)) for p in selected],
            'rank_certificate':union_proof,'rank_lower_bound':len(indices),
            'claim_boundary':'Exact public-point reproduction on the earlier locally found curve188. The first27 union points are its previous certified subgroup; the last28 come from ICARM619 and are explicitly public oracle data. Fixed good-prime finite quotients through997 certify the displayed union lower bound and independently the public28 basis. This is not a new curve discovery, a prospective recovery, an exact rank or an upper bound on the union. No candidate scoring, validation-prime selection, point search or frozen campaign change.'}
    checkpoint(OUT,result)


def check():
    d=cert.read(OUT);old,public,model,transport,u=inputs()
    if d['sources']!=sources() or d['transport_scale_u']!=str(u) or d['public_points']!=public['points']:
        raise ArithmeticError('public reproduction binding differs')
    expected=old['points']+[list(map(str,p)) for p in transport]
    if d['points']!=expected or d['transported_public_points']!=expected[27:]:raise ArithmeticError('point provenance differs')
    points=[tuple(map(cert.F,p)) for p in d['points']];cache=ReductionCache(MemoryFactStore());pivots={}
    for s in d['union_signatures']:
        actual=asdict(cache.signature(model,points,s['prime']))
        if json.dumps(actual,sort_keys=True)!=json.dumps(s,sort_keys=True):raise ArithmeticError('finite signature differs')
        for row in actual['rows']:insert(pivots,row)
    if sorted(pivots)!=d['independent_column_indices']:raise ArithmeticError('independent union columns differ')
    chosen=[points[i] for i in sorted(pivots)]
    if [list(map(str,p)) for p in chosen]!=d['independent_points']:raise ArithmeticError('exported basis differs')
    for pts,key in [(chosen,'rank_certificate'),(transport,'public_rank_certificate')]:
        proof=d[key];actual=memory.checked_rank(model,pts,[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
        if json.dumps(actual,sort_keys=True)!=json.dumps(proof,sort_keys=True):raise ArithmeticError('rank proof replay differs')
    if d['rank_lower_bound']!=len(chosen):raise ArithmeticError('claimed lower bound differs')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args()
    check() if a.check else build()
    print('PUBLIC28 REPRODUCTION PASS; UNION LOWER BOUND',cert.read(OUT)['rank_lower_bound'])
