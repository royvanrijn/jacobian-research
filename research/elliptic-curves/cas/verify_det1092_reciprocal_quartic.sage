#!/usr/bin/env sage-python
"""Universal pointed-quartic identities, then one generic parent application.

No exceptional point, member calibration, specialization roster, or search
input. This identifies both the curve and its degree-four covering map.
Run under timeout25s; existing deterministic certificate is replayed.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,PolynomialRing,EllipticCurve

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_reciprocal_quartic_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
PENCIL=ART/'det1092_norm8_seed_cover_v2/generic.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,data):
    if p.exists():assert json.loads(p.read_text())==data
    else:
        with p.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
def invariants(f):
    e,d,c,b,a=f.list()
    return 12*a*e-3*b*d+c*c,72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3

OUT.mkdir(exist_ok=True)
save(OUT/'protocol.json',{
    'classification':'verified application of classical binary-quartic descent',
    'scope':'Universal trace-centred quartic and generic norm8 pencil; no specialized point input.',
    'limits':{'seconds':25,'generic_pencils':1,'point_searches':0,'exceptional_point_inputs':0},
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [PARENT,PENCIL,Path(__file__)]}})

# Independent symbolic coefficients, with the trace constrained only to E.
A=PolynomialRing(QQ,['a','c','d']);a,c,d=A.gens()
L=PolynomialRing(A,'m');m=L.gen();b=d*d-c**3-a*c
q=m**4-6*c*m*m-8*d*m-3*c*c-4*a
G=c*m**4+4*d*m**3+(6*c*c+4*a)*m*m+4*c*d*m+c**3+4*b
H=(q.derivative()*G-q*G.derivative())/2;H=L(H)
I,J=invariants(q)
assert I==-48*a and J==-1728*b
assert q.discriminant()==256*(-16*(4*a**3+27*b*b))
assert H*H==4*G**3+4*a*G*q*q+4*b*q**3
assert G.derivative()*q-G*q.derivative()==-2*H
assert G.degree()==4 and H.degree()==6 and H[6]==2*d
V=PolynomialRing(L,'w');w=V.gen()
X=(w-c+m*m)/2;Y=m*(X-c)-d
assert (Y*Y-X**3-a*X-b)%(w*w-q)==0
assert (q.derivative()+4*m*w)/4==2*Y
# The two signs lie on the same chord through -Z and add to Z.
Xbar=(-w-c+m*m)/2;Ybar=m*(Xbar-c)-d
assert X+Xbar==m*m-c and Y-Ybar==m*w
universal={
    'classification':'established theory, independently checked universal identities',
    'quartic':'q(m)=m^4-6*c*m^2-8*d*m-3*c^2-4*a',
    'trace_constraint':'d^2=c^3+a*c+b',
    'invariants':{'I':'-48*a','J':'-1728*b','quartic_discriminant':'256*Delta(E)'},
    'curve_isomorphism':'X=(w-c+m^2)/2; Y=m*(X-c)-d',
    'inverse':'m=(Y+d)/(X-c); w=2*X+c-m^2',
    'origin':'positive infinity maps to O; negative infinity maps to Z=(c,d)',
    'G':str(G),'H':str(H),
    'cover_map':'Phi=(G/q,H/(2*q*w))=Z-2*P',
    'differential':'isomorphism pulls dX/(2Y) to dm/w; Phi pulls it to -2*dm/w',
    'involution':'P maps to Z-P',
    'branch_condition':'w=0 iff 2*P=Z',
    'covering_class':'delta(Z); its image in H^1(K,E) and Sha is zero',
    'boundary':'The covering class is not the Kummer class of the marked point P. No full Selmer group is computed.'}

# Same universal identities applied to the entire generic pencil, transposed.
parent=json.loads(PARENT.read_text());pencil=json.loads(PENCIL.read_text())
R=PolynomialRing(QQ,'t');t=R.gen();K=R.fraction_field()
def value(v):return K(R(v['numerator']))/R(v['denominator'])
E=EllipticCurve(K,[value(v) for v in parent['a_invariants']])
aa,bb=-E.c4()/48,-E.c6()/864
old=[E([value(v) for v in point]) for point in parent['basis_weierstrass_coordinates']]
Z=old[14]-old[15]
h,nx,ny,shift=[R(pencil[key]) for key in ['pole_h','nx','ny','shift']]
cx,cy=K(nx)/(h*h),K(ny)/(h**3)
assert cx==Z[0]+E.b2()/12 and cy==Z[1]+(E.a1()*Z[0]+E.a3())/2
S=PolynomialRing(K,'z');z=S.gen()
f=sum((S(list(map(QQ,row)))*t**i for i,row in enumerate(pencil['quartic_t_coefficients_in_z'])),S.zero())
slope=h*z-shift/h
assert h*h*f==slope**4-6*cx*slope*slope-8*cy*slope-3*cx*cx-4*aa
ii,jj=invariants(f)
assert ii==-48*aa and jj==-1728*bb and f[4]==h*h
assert -27*ii==6**4*aa and -27*jj==6**6*bb
def rec(v):return {'numerator':list(map(str,v.numerator().list())),'denominator':list(map(str,v.denominator().list()))}
result={
    'status':'PASS_UNIVERSAL_COVER_CLASS_AND_GENERIC_RECIPROCAL_IDENTITY',
    'classification':'verified application and new programme obstruction',
    'universal':universal,
    'generic':{'I':rec(ii),'J':rec(jj),'leading_coefficient':rec(f[4]),
               'Jacobian_scaling':'X_J=36*X_E; Y_J=216*Y_E',
               'curve_with_origin':'For fixed t, W^2=F_z(t), origin W/z^2=+h(t), is E_t.',
               'known_centre_word':[0]*14+[1,-1,0]},
    'consequence':'The reciprocal member-selection quartic is not a smaller auxiliary descent problem. Its covering class is the inherited centre class on every smooth fibre in the chart.',
    'checker_sha256':sha(Path(__file__)),
    'boundary':'Does not identify Jac(C_z) with a fixed original fibre. Fixed-z carriers still provide genuine base changes, and independent marked points can exist on these inherited covering classes.'}
save(OUT/'generic.json',result)
print(result['status'])
