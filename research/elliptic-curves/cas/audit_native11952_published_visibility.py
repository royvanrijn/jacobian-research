#!/usr/bin/env python3
"""Separate retrospective published-point visibility; never a prospective selector."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
from memory_rank_certificate import checked_rank
from search_observability import point_visibility
from audit_recorded_point_mod2_rank_v3 import signature,insert
from mod2_reduction_independence import _primes_up_to
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/native11952-published-visibility-v1'
PUBLIC=ART/'elliptic_elkies_klagsbrun_rank29_certificate.json';CONTROL=ROOT/'artifacts/local/elliptic-curves/native11952-pari49-control-v1/candidate-00/result.json';CLOUD=ART/'native11952_pari49_recorded_mod2_v1.json';OUT=ART/'native11952_published_visibility_v1.json'
THRESHOLDS=(100000,1000000)
def bindings():
    paths=[Path(__file__).resolve(),PUBLIC,CONTROL,CLOUD]
    paths.extend(ROOT/'elliptic-curves/cas'/s for s in ['search_observability.py','memory_rank_certificate.py','audit_recorded_point_mod2_rank_v3.py','mod2_reduction_independence.py','research_runtime/finite_reduction.py','research_runtime/memory_store.py','elliptic_candidate_record.py'])
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}
def finite(model,points,torsion):
    cache=ReductionCache(MemoryFactStore());pivots={};primes=[]
    for p in _primes_up_to(997):
        if p==2:continue
        try:s=signature(cache,model,points,p)
        except ValueError:continue
        before=len(pivots)
        for row in s.rows:insert(pivots,row)
        if len(pivots)>before:primes.append(p)
        if len(pivots)==len(points):break
    ids=sorted(pivots);chosen=[points[i] for i in ids]
    return {'rank_lower_bound':len(ids),'independent_column_indices':ids,'points':[list(map(str,p)) for p in chosen],'rank_certificate':checked_rank(model,chosen,primes,torsion)}
def expected():
    protocol=cert.read(D/'protocol.json')
    if protocol['bindings']!=bindings() or protocol['height_thresholds']!=list(THRESHOLDS):raise ArithmeticError('frozen retrospective protocol differs')
    public,control,cloud=map(cert.read,(PUBLIC,CONTROL,CLOUD));model=tuple(map(cert.F,control['curve']))
    if len(control['charts'])!=49 or cloud['input_sha256']!=cert.hashed(CONTROL) or cloud['rank_lower_bound']!=27:raise ArithmeticError('completed blind27 control required')
    source=tuple(map(cert.F,public['general_weierstrass_coefficients']));inv=cert.weierstrass_invariants(source)
    if model!=(0,0,0,-inv['c4']/48,-inv['c6']/864):raise ArithmeticError('published/native short model differs')
    a1,a2,a3,a4,a6=source;points=[]
    for p in public['published_points']:
        x,y=cert.F(p['x']),cert.F(p['y']);points.append((x+inv['b2']/12,y+(a1*x+a3)/2))
    if len(points)!=29:raise ArithmeticError('published roster differs')
    old=public['finite_reduction_certificate'];public_proof=checked_rank(model,points,old['certificate_primes'],old['two_torsion_certificate_prime'])
    baseline=[tuple(map(cert.F,p)) for p in cloud['independent_points']];bp=cloud['rank_certificate']
    if json.dumps(checked_rank(model,baseline,[s['prime'] for s in bp['signatures']],bp['no_rational_2_torsion_prime']),sort_keys=True)!=json.dumps(bp,sort_keys=True):raise ArithmeticError('blind27 proof differs')
    rows=[]
    for i,(x,y) in enumerate(points):
        observations=[]
        for j,chart in enumerate(control['charts']):
            record=chart['search']
            if record['status']!='bounded_search_complete' or record['height_bound']!=100000 or record['completed_denominator']!=100000 or not record['infinity_checked']:raise ArithmeticError('incomplete baseline box')
            for sign in (1,-1):
                v=point_visibility(record,(x,sign*y))
                if v['status'] in ('VISIBLE_NOT_RECORDED','UNSEARCHED_INTERVAL'):raise ArithmeticError('baseline visibility discrepancy')
                observations.append({'chart':j,'sign':sign,**v})
        heights=[v['minimum_affine_height'] for v in observations if v.get('minimum_affine_height') is not None]
        endpoint=any(v['status']=='KNOWN_POINTED_ENDPOINT' or v.get('at_parameter_infinity',False) for v in observations)
        rows.append({'published_index':i,'point':list(map(str,(x,y))),'minimum_height_over49_and_signs':0 if endpoint else min(heights),'endpoint_or_infinity':endpoint,'observations':observations})
        print('PUBLIC VISIBILITY',i,rows[-1]['minimum_height_over49_and_signs'],flush=True)
    tiers=[]
    for h in THRESHOLDS:
        ids=[r['published_index'] for r in rows if r['minimum_height_over49_and_signs']<=h]
        tiers.append({'height':h,'visible_published_indices':ids,'blind27_plus_visible':finite(model,baseline+[points[i] for i in ids],bp['no_rational_2_torsion_prime'])})
    return {'schema':'elliptic-curves.native11952-published-visibility.v1','status':'PASS','bindings':bindings(),'protocol_sha256':cert.hashed(D/'protocol.json'),'curve':control['curve'],'published_rank29_proof':public_proof,'rows':rows,'height_tiers':tiers,'all_published_plus_blind27':finite(model,baseline+points,bp['no_rational_2_torsion_prime']),'claim_boundary':'Separate retrospective oracle audit after the fixed generic17-only control completed. Exact coordinates and independently certified lower bounds of indicated unions. No search at larger heights occurred here, no oracle point entered prospective selection, and a point outside these boxes does not exclude other representatives of its Mordell-Weil coset or any unseen point. Not a new curve or new record.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();data=expected()
    if a.check:
        if cert.read(OUT)!=data:raise ArithmeticError('published visibility replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve retrospective certificate')
        checkpoint(OUT,data)
    print('PUBLIC VISIBILITY PASS',[(r['height'],len(r['visible_published_indices']),r['blind27_plus_visible']['rank_lower_bound']) for r in data['height_tiers']],flush=True)
