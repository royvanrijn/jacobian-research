#!/usr/bin/env python3
"""Independent Sage isomorphism replay of the frozen current-catalogue comparison."""
import json,hashlib
from pathlib import Path
from sage.all import QQ,EllipticCurve

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'inventory200_current_catalogue_comparison_v1.json'
OUT=ART/'inventory200_current_catalogue_sage_replay_v1.json'
if OUT.exists():raise FileExistsError('preserve independent Sage replay')
d=json.loads(INPUT.read_bytes())
for n,h in d['sources'].items():
    if hashlib.sha256((ROOT/n).read_bytes()).hexdigest()!=h:raise ArithmeticError('source changed')
byj={}
for r in d['equations']:
    E=EllipticCurve(QQ,r['ainvs']);byj.setdefault(E.j_invariant(),[]).append((r['id'],E))
checked=[]
for r in d['inventory_comparisons']:
    E=EllipticCurve(QQ,r['curve']);j=E.j_invariant();candidates=byj.get(j,[])
    samej=[ident for ident,F in candidates];matches=[ident for ident,F in candidates if E.is_isomorphic(F)]
    if str(j)!=r['j_invariant'] or samej!=r['same_j_ids'] or matches!=r['q_isomorphism_matches']:
        raise ArithmeticError('independent Sage comparison differs')
    checked.append({'id':r['id'],'q_isomorphism_matches':matches})
if len(checked)!=200 or sum(len(v) for v in byj.values())!=620:raise ArithmeticError('incomplete comparison')
result={'schema':'elliptic-curves.inventory200-current-catalogue-sage.v1','status':'PASS',
        'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__).resolve(),INPUT)},
        'inventory_count':200,'catalogue_count':620,'comparisons':checked,
        'claim_boundary':'Independent Sage rational-isomorphism replay for every inventory200 equation against every equal-j class of the620-equation snapshot. Distinct j excludes all other equations. No point search, rank recertification, first-discovery priority or universal novelty.'}
with OUT.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
print('SAGE INVENTORY200 VS620 PASS', [r for r in checked if r['q_isomorphism_matches']])
