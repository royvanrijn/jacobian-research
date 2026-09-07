#!/usr/bin/env python3
"""Exact maps and post-freeze public-witness coordinates on the own27 control."""
import json,hashlib,sys
from decimal import Decimal,localcontext
from pathlib import Path
from sage.all import QQ,ZZ,matrix,EllipticCurve,PolynomialRing
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import inventory188_own27_geometry_control as control
from research_runtime.store import digest
ART=control.ART;D=control.D
OUT=ART/'inventory188_own27_geometry_visibility_v1.json'
if OUT.exists():raise FileExistsError('preserve own27 geometry/witness audit')
p=control.protocol();m=json.loads((D/'maps.json').read_bytes());seed=json.loads(control.SEED.read_bytes())
supervision=json.loads((D/'geometry.supervisor.json').read_bytes())
if supervision['outcome']!='completed' or supervision['returncode']!=0 or m['status']!='COMPLETE_DECLARED_MAPS' or m['protocol_hash']!=digest(p):
    raise ArithmeticError('all geometry must finish before opening public witness')
g=matrix(ZZ,m['rounded_gram']);u=matrix(ZZ,m['change_of_basis']);h=matrix(ZZ,m['reduced_gram'])
with localcontext() as ctx:
    ctx.prec=110;rounded=[[int((Decimal(v)*1000000).to_integral_value()) for v in r] for r in m['metric_gram']]
if g!=matrix(ZZ,rounded) or not g.is_symmetric() or not g.is_positive_definite() or abs(u.det())!=1 or h!=u*g*u.transpose():
    raise ArithmeticError('rounded metric/unimodular basis transport differs')
if [r['parity'] for r in m['sample']]!=control.masks(p):raise ArithmeticError('fixed2048 sample differs')
for r in m['sample']:
    w=matrix(ZZ,1,27,r['representative']);v=matrix(ZZ,1,27,r['reduced_representative'])
    if w!=v*u or int((w*g*w.transpose())[0,0])!=r['metric_norm'] or any((int(w[0,j])-(r['parity']>>j))%2 for j in range(27)):
        raise ArithmeticError('exact sampled parity or norm differs')
selected=sorted(m['sample'],key=lambda r:(-r['metric_norm'],r['parity']))[:49]
if m['centres']!=selected or [r['centre'] for r in m['rows']]!=selected:raise ArithmeticError('fixed largest-norm selection differs')
E=EllipticCurve(QQ,seed['curve']);basis=[E(P) for P in seed['points']];R=PolynomialRing(QQ,'z');z=R.gen();checked=[]
for i,r in enumerate(m['rows']):
    C=sum((int(a)*P for a,P in zip(r['centre']['representative'],basis)),E(0));x,y=C.xy()
    raw=[-3*x*x-4*E.a4(),-8*y,-6*x,0,1]
    if raw!=list(map(QQ,r['raw_coefficients'])):raise ArithmeticError('centre/raw quartic differs')
    M=matrix(QQ,2,2,r['matrix']);M1=matrix(QQ,2,2,r['first_matrix']);M2=matrix(QQ,2,2,r['second_matrix'])
    if M!=M1*M2 or M.det()==0:raise ArithmeticError('horizontal map composition differs')
    A=M[0,0]*z+M[0,1];B=M[1,0]*z+M[1,1]
    transformed=sum(raw[k]*A**k*B**(4-k) for k in range(5));quartic=R(list(map(QQ,r['discriminant_quartic'])))
    ratio=QQ(r['square_ratio'])
    if ratio<=0 or not ratio.is_square() or transformed!=ratio*quartic:raise ArithmeticError('exact quartic identity differs')
    checked.append((C,M))
# Only after every selection and geometry identity passes is this oracle read.
public_path=ART/'inventory188_public28_reproduction_v1.json';public_replay_path=ART/'inventory188_public28_sage_replay_v1.json'
public=json.loads(public_path.read_bytes());replay=json.loads(public_replay_path.read_bytes())
for data in (public,replay):
    if data['status']!='PASS':raise ArithmeticError('independent public-point proof required')
    for n,hsh in data['sources'].items():
        if hashlib.sha256((ROOT/n).read_bytes()).hexdigest()!=hsh:raise ArithmeticError('public source changed')
if public['curve']!=seed['curve'] or public['points'][:27]!=seed['points'] or public['independent_column_indices']!=list(range(27))+[53]:
    raise ArithmeticError('same old27 and certified extra representative required')
P=E(public['transported_public_points'][26]);observations=[]
for chart,(C,M) in enumerate(checked):
    for sign in (-1,1):
        Q=sign*P
        if Q[0]==C[0]:raise ArithmeticError('extra point cannot be an old-subgroup endpoint')
        t=(Q[1]+C[1])/(Q[0]-C[0]);n=M[1,1]*t-M[0,1];den=M[0,0]-M[1,0]*t
        if den==0:coordinate=['1','0'];height=1;infinity=True
        else:
            s=n/den;coordinate=[str(s.numerator()),str(s.denominator())];height=int(max(abs(s.numerator()),s.denominator()));infinity=False
        observations.append({'chart':chart,'sign':sign,'coordinate':coordinate,'projective_height':height,'at_infinity':infinity,
                             'inside_height125000_or_infinity':infinity or height<=125000})
best=min(observations,key=lambda r:(r['projective_height'],r['chart'],r['sign']))
paths=[Path(__file__).resolve(),D/'protocol.json',D/'maps.json',D/'geometry.supervisor.json',control.SEED,public_path,public_replay_path]
result={'schema':'elliptic-curves.inventory188-own27-geometry-visibility.v1','status':'PASS_EXACT_GEOMETRY_AND_COORDINATES',
        'sources':{str(q.relative_to(ROOT)):hashlib.sha256(q.read_bytes()).hexdigest() for q in paths},
        'sample_masks_checked':2048,'exact_maps_checked':49,'public_point_index':26,'observations':observations,
        'best_observation':best,'visible_observations':sum(r['inside_height125000_or_infinity'] for r in observations),
        'claim_boundary':'Fixed own27 geometry was constructed from the old local subgroup without public-point inputs. Exact parity/norm transports, quartic identities and all49 maps are checked before the independently certified public witness is opened for98 coordinate observations. Numerical canonical heights and CVP guide selection only. No point enumeration, prospective recovery, rank increase, absence of cheaper translates or optimality theorem.'}
with OUT.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
print('OWN27 GEOMETRY PASS; VISIBLE',result['visible_observations'],'BEST',best)
