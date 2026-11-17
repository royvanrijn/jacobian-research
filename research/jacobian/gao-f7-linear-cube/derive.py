"""Derive the seven-two necessary coefficient system; no exclusion decision."""
from pathlib import Path
import hashlib,json,time
import sympy as S
HERE=Path(__file__).resolve().parent;START=time.monotonic()
t,r,q,z,h,k,Y1,Y2,Y3,Y4=S.symbols('t r q z h k Y1 Y2 Y3 Y4')
b=Y4-3*Y2*t*t+2*Y1*t**3-t**4
a=Y3-2*Y2*t+Y1*t*t-t**3-t*t*b
R=S.expand(a*a+a+t*b+t**3+Y1*t-Y2)
P=S.Poly(S.expand((t-r)**7*(t-q)**2*(t**3+z*t*t+h*t+k)),t)
y1=-P.nth(11)/4;y2=(P.nth(10)-4*y1*y1)/6
assert S.diff(P.nth(9),k)==1
kfit=S.expand(-12*y1*y2-2-P.nth(9).subs(k,0))
P=S.Poly(S.expand(P.as_expr().subs(k,kfit)),t)
y4=(6*y1+9*y2*y2-P.nth(8))/2
y3=(P.nth(6)-14*y1*y2+6*y2*y4-2)/2
res=S.Poly(S.expand(R.subs({Y1:y1,Y2:y2,Y3:y3,Y4:y4})-P.as_expr()),t)
for j in [12,11,10,9,8,6]:assert res.nth(j)==0
rows=[]
for j in [7,5,4,3,2,1,0]:
 poly=S.Poly(res.nth(j),r,q,z,h,domain=S.QQ).clear_denoms()[1].primitive()[1]
 rows.append({'coefficient':j,'polynomial':str(poly.as_expr()),'total_degree':poly.total_degree(),'terms':len(poly.terms())})
old=json.loads((HERE.parent/'gao-f7-three-point/linear_cube_cubic_square-input.json').read_text())
u,v,w=S.symbols('u v w');transport={u:-2*r-q,v:r*r+2*r*q,w:-r*r*q}
for row in rows:
 oldrow=next(x for x in old['equations'] if x['coefficient']==row['coefficient'])
 value=S.expand(S.sympify(oldrow['polynomial']).subs(transport))
 ratio=S.cancel(value/S.sympify(row['polynomial']))
 assert ratio.is_Rational and ratio!=0
f7=S.sympify(next(x['polynomial'] for x in old['equations'] if x['coefficient']==7))
Delta=(z-u)*(z+3*r)+2*(v-h)
assert S.expand(S.diff(f7,w)-384*Delta)==0
result={'variables':['r','q','z','h'],'factorization':'R=(t-r)^7*(t-q)^2*S3, S3 monic cubic','equations':rows,'transport_from_general_family':{str(x):str(y) for x,y in transport.items()},'scope':'Necessary equations only. No squarefree/disjointness assumptions or r!=q restriction. No common-zero exclusion is claimed.'}
(HERE/'seven-two-input.json').write_text(json.dumps(result,indent=2)+'\n')
receipt={'result':'NECESSARY_EQUATIONS_DERIVED','equation_count':7,'variables':4,'source_coefficient_transport':'PASS','general_linear_w_coefficient':'384*((z-u)*(z+3r)+2*(v-h))','seconds':time.monotonic()-START,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'theorem_exclusion':False}
(HERE/'derivation-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt));print([(row['coefficient'],row['total_degree'],row['terms']) for row in rows])
