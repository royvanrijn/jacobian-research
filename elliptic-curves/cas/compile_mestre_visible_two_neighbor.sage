#!/usr/bin/env sage-python
"""Compile one fixed effective pencil O+C0 on u11, no neighbour enumeration."""
import json,hashlib
from pathlib import Path
from sage.all import QQ,PolynomialRing,EllipticCurve
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
BUNDLE=ART/'mestre_468_replay_bundle_v1.json';OUT=ART/'mestre_u11_visible_two_neighbor_v1.json'
row=next(r for r in json.loads(BUNDLE.read_text())['rows'] if r['outer_u']=='11');h=row['generic_heights']
R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();A=R(h['curve_A_coefficients']);B=R(h['curve_B_coefficients']);E=EllipticCurve(K,[A,B]);P=E([K(c) for c in h['covariant_points'][0]])
q=P[0].denominator().sqrt();assert q in R and q.degree()==2
q=R(q);a=R(P[0]*q*q);b=R(P[1]*q**3)
assert b*b==a**3+A*a*q**4+B*q**6 and a.gcd(q)==1
c=(b*a.inverse_mod(q*q))%(q*q)
S=PolynomialRing(R,'z');z=S.gen();mnum=q*q*z-c
N=mnum**4-6*a*mnum*mnum-8*b*mnum-3*a*a-4*A*q**4
coeff=[]
for f in N.list():
 div,rem=f.quo_rem(q**6);assert not rem;coeff.append(div)
quartic=S(coeff);assert max(f.degree() for f in coeff)==4
Z=PolynomialRing(QQ,'z');zz=Z.gen();quartic_T=[Z([f[i] for f in coeff]) for i in range(5)]
e,d,cc,bb,aa=quartic_T;I=12*aa*e-3*bb*d+cc*cc;J=72*aa*cc*e+9*bb*cc*d-27*aa*d*d-27*bb*bb*e-2*cc**3
newA=-27*I;newB=-27*J;delta=-16*(4*newA**3+27*newB**2)
print('Jacobian degrees',newA.degree(),newB.degree(),'Delta',delta.degree(),flush=True)
print('factor discriminant',flush=True);fac=delta.factor();print([(f.degree(),k) for f,k in fac],flush=True)
print('common c4',delta.gcd(newA).degree(),flush=True)
old_delta=-16*(4*A**3+27*B**2);double=next(f for f,k in old_delta.squarefree_decomposition() if k==2)
points=[]
for t in double.roots(QQ,multiplicities=False):
 val=sum(f*t**i for i,f in enumerate(quartic_T));yy=val.sqrt();assert yy in Z and yy*yy==val
 points += [[str(t),str(yy)],[str(t),str(-yy)]]
result={'schema':'elliptic-curves.mestre-visible-two-neighbor.v1','status':'PASS','outer_u':'11','old_A':list(map(str,A.list())),'old_B':list(map(str,B.list())),
'pole_point':list(map(str,P.xy())),'q':list(map(str,q.list())),'a':list(map(str,a.list())),'b':list(map(str,b.list())),'c':list(map(str,c.list())),
'quartic_T_coefficients':[list(map(str,f.list())) for f in quartic_T], 'Jacobian_A':list(map(str,newA.list())),'Jacobian_B':list(map(str,newB.list())),
'discriminant_factors':[{'coefficients':list(map(str,f.list())),'multiplicity':int(k)} for f,k in fac], 'discriminant_unit':str(fac.unit()),'c4_gcd_degree':int(delta.gcd(newA).degree()),
'quartic_sections':points,'map':'z=(q*(y+yP)/(x-xP)+c)/q^2; W=(2*x+xP-m^2)/q, m=q*z-c/q',
'inverse':'m=q*z-c/q; x=(q*W-xP+m^2)/2; y=m*(x-xP)-yP',
'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).resolve(),BUNDLE]},
'scope':'One exact chord substitution and its rational inverse, with rational quartic sections. Generic MW independence, minimal Jacobian fibre geometry and useful search admission require separate checks. No point search or other pencil.'}
assert not OUT.exists();OUT.write_text(json.dumps(result,indent=2)+'\n');print('PASS explicit u11 O+C0 pencil',flush=True)
