#!/usr/bin/env python3
"""Separate post-terminal comparison of the new60 cohort with the fresh620 snapshot."""
import json,hashlib
from pathlib import Path
from sage.all import QQ,EllipticCurve
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
CAT=ART/'inventory200_current_catalogue_comparison_v1.json'
COHORT=ART/'nearcut60v2_mw16_novelty_v1.json'
OUT=ART/'nearcut60_current_catalogue_sage_v1.json'
if OUT.exists():raise FileExistsError('preserve current-cohort comparison')
catalogue=json.loads(CAT.read_bytes());cohort=json.loads(COHORT.read_bytes())
for data in (catalogue,cohort):
    if data['status']!='PASS':raise ArithmeticError('completed comparison required')
    for n,h in data['sources'].items():
        if hashlib.sha256((ROOT/n).read_bytes()).hexdigest()!=h:raise ArithmeticError('frozen evidence changed')
byj={}
for r in catalogue['equations']:
    E=EllipticCurve(QQ,r['ainvs']);byj.setdefault(E.j_invariant(),[]).append((r['id'],E))
rows=[]
for r in cohort['rows']:
    E=EllipticCurve(QQ,r['curve']);j=E.j_invariant()
    rows.append({'id':r['id'],'rank_lower_bound':r['rank_lower_bound'],
                 'q_isomorphism_matches':[ident for ident,F in byj.get(j,[]) if E.is_isomorphic(F)]})
if len(rows)!=60 or catalogue['catalogue_count']!=620:raise ArithmeticError('fixed60 against620 required')
result={'schema':'elliptic-curves.nearcut60-current-catalogue.v1','status':'PASS',
        'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__).resolve(),CAT,COHORT)},
        'catalogue_count':620,'comparisons':rows,
        'claim_boundary':'Separate Sage Q-isomorphism comparison against the fresh620-equation snapshot after the fixed60 trial and its independent accounting completed. No selection, score, point, validation or budget changes. Prior593-equation proofs are preserved. No universal novelty or rank recertification.'}
with OUT.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
print('NEARCUT60 CURRENT620 CATALOGUE PASS; MATCHES',[r for r in rows if r['q_isomorphism_matches']])
