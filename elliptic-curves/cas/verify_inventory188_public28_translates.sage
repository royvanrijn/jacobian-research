#!/usr/bin/env python3
"""Independent group additions and inverse-Mobius coordinates for all55 translates."""
import json,hashlib
from pathlib import Path
from sage.all import QQ,EllipticCurve,matrix
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'inventory188_public28_translates_v1.json'
OUT=ART/'inventory188_public28_translates_sage_v1.json'
if OUT.exists():raise FileExistsError('preserve independent translation replay')
d=json.loads(INPUT.read_bytes())
for n,h in d['sources'].items():
    if hashlib.sha256((ROOT/n).read_bytes()).hexdigest()!=h:raise ArithmeticError('bound source changed')
source=json.loads((ART/'full11952_late64_r17_results_v1.json').read_bytes())
r=next(r for r in source['curves'] if r['parameter']=='110314/102227')
rawpath=ROOT/r['discovery_witness']['path']
if hashlib.sha256(rawpath.read_bytes()).hexdigest()!=r['discovery_witness']['sha256']:raise ArithmeticError('original charts changed')
raw=json.loads(rawpath.read_bytes());E=EllipticCurve(QQ,d['curve']);P=E(d['public_point']);basis=[E(p) for p in d['old_basis']]
maps=[]
for c in raw['charts']:
    s=c['search'];saved=s['pointed_chart'];den=QQ(saved['point_denominator_root']);scale=QQ(saved['curve_coordinate_scale']);k=QQ(saved['shift_mod_denominator_squared'])
    M=matrix(QQ,[[den/scale,k/(den*scale)],[0,1]])*matrix(QQ,2,2,saved['unimodular_horizontal_matrix'])*matrix(QQ,2,2,s['horizontal_matrix'])
    maps.append((E(s['base_point']),M))
checked=0;all_outside=True;best_rows=[]
for row in d['rows']:
    i=row['basis_index'];Q=P if i is None else P+row['coefficient']*basis[i]
    if Q!=E(row['point']):raise ArithmeticError('independent group addition differs')
    best=None
    for chart,(C,M) in enumerate(maps):
        for sign in (-1,1):
            R=sign*Q
            if R[0]==C[0]:raise ArithmeticError('unexpected endpoint in extra coset')
            t=(R[1]+C[1])/(R[0]-C[0]);n=M[1,1]*t-M[0,1];v=M[0,0]-M[1,0]*t
            if v==0:raise ArithmeticError('unexpected parameter infinity')
            s=n/v;height=max(abs(s.numerator()),s.denominator());key=(int(height),chart,sign)
            all_outside=all_outside and height>125000;checked+=1
            if best is None or key<best[0]:best=(key,str(s))
    recorded=row['best_observation'];v=recorded['visibility']
    if best[0]!=(v['minimum_affine_height'],recorded['chart'],recorded['ordinate_sign']) or QQ(best[1])!=QQ(v['coordinate'][0])/QQ(v['coordinate'][1]):
        raise ArithmeticError('independent coordinate minimum differs')
    best_rows.append({'word_index':row['word_index'],'minimum_height':best[0][0]})
if checked!=5390 or not all_outside:raise ArithmeticError('claimed complete outside-box result differs')
result={'schema':'elliptic-curves.inventory188-translates-sage.v1','status':'PASS',
        'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__).resolve(),INPUT,rawpath)},
        'translations_checked':55,'coordinates_checked':checked,'all_outside_height125000':True,'minima':best_rows,
        'claim_boundary':'Independent Sage group additions and rational inverse-Mobius coordinates for exactly the55 frozen translations, both signs and49 original charts. No point search, untested-translate exclusion, upper rank or prospective recovery.'}
with OUT.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
print('INDEPENDENT SAGE55 TRANSLATES;5390 OUTSIDE-BOX COORDINATES PASS')
