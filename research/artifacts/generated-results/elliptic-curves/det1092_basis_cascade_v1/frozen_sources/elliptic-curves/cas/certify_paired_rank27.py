#!/usr/bin/env python3
"""Rank27, global minimality and post-batch equation comparison, without Sage."""
import argparse,json
from math import gcd
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/fresh-r17-paired-pari-v1';CLOUD=ART/'fresh_r17_paired_074d9_007_mod2_v1.json';OUT=ART/'new_paired_rank27_proof_v1.json';DATABASE=ROOT/'artifacts/local/elliptic-curves/next12-current-catalogue-v1/database.json';OLD=ART/'compact_r17_wide_results_v1.json'

def sources():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),Path(cert.__file__).resolve(),ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',ROOT/'elliptic-curves/cas/elliptic_candidate_record.py')}

def verify(d):
    if d['sources']!=sources():raise ArithmeticError('rank27 proof source changed')
    model=tuple(map(cert.F,d['minimal_curve']));inv=cert.weierstrass_invariants(model)
    if any(q.denominator!=1 for q in model) or not inv['discriminant'] or gcd(int(inv['c4']),int(inv['c6']))!=16:raise ArithmeticError('integral minimality criterion differs')
    delta=abs(int(inv['discriminant']))
    if delta%256 or (delta//256)%2==0:raise ArithmeticError('minimal discriminant at2 criterion differs')
    short=tuple(map(cert.F,d['discovery_curve']))
    if short!=(0,0,0,-inv['c4']/48,-inv['c6']/864):raise ArithmeticError('short equation differs')
    a1,a2,a3,a4,a6=model;b2=inv['b2'];points=tuple(tuple(map(cert.F,p)) for p in d['points']);transported=[(x+b2/12,y+(a1*x+a3)/2) for x,y in points]
    if any(y*y+a1*x*y+a3*y!=x*x*x+a2*x*x+a4*x+a6 for x,y in points) or [list(map(str,p)) for p in transported]!=d['discovery_points']:raise ArithmeticError('point transport differs')
    proof=d['rank_certificate_on_discovery_model'];actual=cert.checked_rank(short,transported,[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
    if len(points)!=27 or json.dumps(actual,sort_keys=True)!=json.dumps(proof,sort_keys=True):raise ArithmeticError('27-point independence differs')
    if d['icarm_matches'] or d['previous_matches'] or any(cert.isomorphic(model,r['ainvs']) for r in d['catalogue']['equations']) or any(cert.isomorphic(model,r['curve']) for r in d['previous_equations']):raise ArithmeticError('catalogue/prior comparison failed')
    print('PASS:27 independent points; global minimality; no Q-isomorphism match in',len(d['catalogue']['equations']),'catalogue and',len(d['previous_equations']),'earlier equations',flush=True)

def build():
    if OUT.exists():raise FileExistsError('preserve rank27 proof')
    ledger=cert.read(D/'ledger.json');p=cert.read(D/'protocol.json');v=cert.read(D/'074d9-007/verification.json');cloud=cert.read(CLOUD)
    if ledger['status']!='COMPLETE_FIXED_BATCH_ATTEMPTS' or [r['id'] for r in ledger['rows']]!=[r['id'] for r in p['rows']] or any(r['status']=='PENDING' for r in ledger['rows']) or v['status']!='PASS' or v['cloud_sha256']!=cert.hashed(CLOUD):raise ArithmeticError('terminal batch and rank27 replay required')
    if cloud['rank_lower_bound']!=27 or cloud['family']!='074d9' or cloud['parameter']!='2818/1535':raise ArithmeticError('rank27 discovery differs')
    short=tuple(map(cert.F,cloud['curve']));inv=cert.weierstrass_invariants(short);a1,a2,a3=0,1,0;b2=cert.F(4);b4=(b2*b2-inv['c4'])/24;b6=(-b2**3+36*b2*b4-inv['c6'])/216;model=(cert.F(a1),cert.F(a2),cert.F(a3),b4/2,b6/4);points=[(cert.F(x)-b2/12,cert.F(y)) for x,y in cloud['independent_points']]
    old=cert.read(OLD);previous=old['previous_equations']+[{'address':OLD.name+':'+r['family']+':'+r['parameter'],'curve':r['curve']} for r in old['curves']];projection=[{'id':r['id'],'ainvs':r['ainvs']} for r in cert.read(DATABASE)['curves']]
    d={'schema':'elliptic-curves.new-paired-rank27-proof.v1','sources':sources(),'inputs':{str(path.relative_to(ROOT)):cert.hashed(path) for path in (CLOUD,D/'protocol.json',D/'ledger.json',D/'074d9-007/verification.json',OLD,DATABASE)},'family':cloud['family'],'parameter':cloud['parameter'],'minimal_curve':list(map(str,model)),'points':[list(map(str,P)) for P in points],'discovery_curve':cloud['curve'],'discovery_points':cloud['independent_points'],'rank_certificate_on_discovery_model':cloud['rank_certificate'],'transport':'x_short = X + 1/3; y_short = Y','minimality':'The integral model has gcd(c4,c6)=16, excluding nonminimality at every odd prime. Its discriminant has2-adic valuation8<12, excluding nonminimality at2.','catalogue':{'url':'https://elliptic-rank.icarm.cloud/database.json','raw_sha256':cert.hashed(DATABASE),'curve_count':len(projection),'equations':projection},'previous_equations':previous,'icarm_matches':[r['id'] for r in projection if cert.isomorphic(model,r['ainvs'])],'previous_matches':[r['address'] for r in previous if cert.isomorphic(model,r['curve'])],'claim_boundary':'Rank at least27 and a global minimal integral equation. No matching Q-isomorphism class in the pinned586-equation catalogue or299 earlier measured equations. Exact rank, exact conductor and universal novelty remain unknown; no new rank28/32 is asserted.'}
    verify(d);checkpoint(OUT,d)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();verify(cert.read(OUT)) if a.check else build()
