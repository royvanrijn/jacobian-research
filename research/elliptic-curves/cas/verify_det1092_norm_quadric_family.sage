#!/usr/bin/env sage-python
"""Universal short-cubic quadric and unchanged nine-control comparison.

No seed data, searches, field/ideal backend or factorization.25s cap.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,PolynomialRing,EllipticCurve,matrix,vector
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_seed_norm_quadric_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
ROSTER=ART/'det1092_rr_generic_point_controls_v2/protocol.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    p=OUT/name
    if p.exists():assert read(p)==data
    else:
        with p.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
save('family-protocol.json',{'classification':'universal equation identity and fixed equation-only control test',
    'rule':'Use only the parent equation and the same nine original parameters. Check the explicit rational quadric point (0,-3,2A,-A) and quadratic projection for the depressed cubic at every address. This tests no exceptional point, ideal class, unit or Selmer condition.',
    'limits':{'seconds':25,'fixed_parameters':9,'new_parameters':0,'point_searches':0,'factorizations':0},
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [PARENT,ROSTER,Path(__file__)]}})
S=PolynomialRing(QQ,['A','B','v1','v2','a','b','r','s','h'])
A,B,v1,v2,aa,bb,r,s,h=S.gens()
Q=matrix(S,[[A,3*B/2,0,0],[3*B/2,-A*A/3,0,0],[0,0,1,QQ(1)/2],[0,0,QQ(1)/2,1]])
X=vector(S,[v1,v2,aa,bb]);p=vector(S,[0,-3,2*A,-A])
assert Q.det()==(-4*A**3-27*B*B)/16 and p*Q*p==0
T=matrix(S,[[0,0,-B],[1,0,-A],[0,1,0]]);Id=matrix.identity(S,3)
V=v1*T+v2*(T*T+2*A/3*Id)
assert (V+aa*Id).det()-(V+bb*Id).det()==(aa-bb)*(X*Q*X)
q=A*r*r+s*s+s*h+h*h;ell=-9*B*r+3*A*s
Phi=vector(S,[-ell*r,-3*q,2*A*q-ell*s,-A*q-ell*h])
assert Phi*Q*Phi==0
assert Phi==(vector(S,[r,0,s,h])*Q*vector(S,[r,0,s,h]))*p-ell*vector(S,[r,0,s,h])
# At p the ratio (V+a)/(V+b) equals 1+(A/B)*theta.
numerator=-3*T*T;denominator=-3*(T*T+A*Id)
assert B*numerator==denominator*(B*Id+A*T)
assert (B*Id+A*T).det()==B**3
R=PolynomialRing(QQ,'t');F=R.fraction_field()
parent=read(PARENT)
def dec(v):return F(R(v['numerator']))/R(v['denominator'])
generic=EllipticCurve(F,[dec(v) for v in parent['a_invariants']])
rows=[]
cases=read(ROSTER)['cases'];assert len(cases)==9
for case in cases:
    tau=QQ(case['parameter']);E=EllipticCurve(QQ,[v(tau) for v in generic.a_invariants()])
    # Depressed 2-division cubic in eta=4x+b2/3; at302 it is theta+b2/3.
    a=-E.c4()/3;b=-2*E.c6()/27;disc=-4*a**3-27*b*b
    assert disc and a and b
    gram=matrix(QQ,[[a,3*b/2,0,0],[3*b/2,-a*a/3,0,0],[0,0,1,QQ(1)/2],[0,0,QQ(1)/2,1]])
    point=vector(QQ,[0,-3,2*a,-a]);assert point*gram*point==0
    assert gram.det()==disc/16
    companion=matrix(QQ,[[0,0,-b],[1,0,-a],[0,1,0]])
    unit=matrix.identity(QQ,3)+(a/b)*companion
    assert unit.det()==1 and (unit-matrix.identity(QQ,3)).det()!=0
    rows.append({'label':case['label'],'parameter':str(tau),'short_cubic_A':str(a),'short_cubic_B':str(b),
        'quadric_point':list(map(str,point)),'cubic_discriminant_sign':int(disc.sign()),
        'cubic_discriminant_square':bool(disc.is_square()),'rational_norm_quadric':True})
report={'status':'PASS_UNIVERSAL_NORM_QUADRIC_AND_NINE_RATIONAL_CONTROLS',
    'classification':'universal elementary identity and generic-only fixed-panel application',
    'formula':'A*v1^2+3*B*v1*v2-A^2*v2^2/3+a^2+a*b+b^2=0',
    'rational_point':'(0,-3,2*A,-A)','discriminant_ratio':'1/16',
    'projection':'q=A*r^2+s^2+s*h+h^2, ell=-9*B*r+3*A*s; Phi=(-ell*r,-3*q,2*A*q-ell*s,-A*q-ell*h)',
    'cases':rows,'boundary':'All nine norm quadrics are rational by a uniform equation-only formula. This does not construct integral norm solutions, units, ideal-class relations, elliptic points or Selmer/Sha classes, and is not a seed discriminator.',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [PARENT,ROSTER,Path(__file__),OUT/'family-protocol.json']}}
save('family-replay.json',report)
print(report['status'],[(r['label'],r['cubic_discriminant_sign']) for r in rows],flush=True)
