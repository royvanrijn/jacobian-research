#!/usr/bin/env sage-python
"""Universal sixth-order identity and three exact normal-form applications.

No geometric black-box or search. Symbolic polynomial identities over Q.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,PolynomialRing
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_euclidean_bisection_survivors_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    text=json.dumps(data,indent=2,sort_keys=True)+'\n';dest=OUT/name
    if dest.exists():assert read(dest)==data
    else:
        with dest.open('x') as stream:stream.write(text)
save('identity-protocol.json',{'classification':'universal polynomial verification and same three constructed conics',
    'limits':{'seconds':25,'new_orbits':0,'new_parameters':0,'point_searches':0,
              'factorizations':0,'symbolic_identity_families':1},
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [PARENT,OUT/'summary.json',Path(__file__)]}})
S=PolynomialRing(QQ,['h','m','b','k','q','W']);h,m,b,k,q,W=S.gens()
A=m*k-3*b*b/4-h*h*q/4
B=(m*m*q+b**3-2*m*b*k+h*h*(k*k-b*q))/4
nx=m*m-h*h*b
ny=-m**3+3*h*h*m*b/2-h**4*k/2
g=(m*b-h*h*k)/2
assert ny*ny==nx**3+A*nx*h**4+B*h**6
assert m*nx+ny==h*h*g
assert m**4-6*nx*m*m-8*ny*m-3*nx*nx-4*A*h**4==h**6*q
X=(b+h*W)/2;Y=-(h*k+m*W)/2
assert -g+m*X+h*Y==0
assert 8*(Y*Y-X**3-A*X-B)==(W*W-q)*(2*m*m-3*b*h*h-h**3*W)
R=PolynomialRing(QQ,'t');K=R.fraction_field();parent=read(PARENT)
def dec(r):return K(R(r['numerator']))/R(r['denominator'])
a1,a2,a3,a4,a6=[dec(r) for r in parent['a_invariants']]
b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6
AA=R(b4/2-b2*b2/48);BB=R(b6/4-b2*b4/24+b2**3/864)
applications=[]
for mask in [61,107,111]:
    rec=read(OUT/('orbit-%d.json'%mask));f0,mm,hh=map(R,rec['line'])
    x0,x1,y0,y1=map(R,rec['elliptic_maps']);bb=2*x0;qq=R(rec['branch_quadratic'])
    kk,remainder=(-2*y0).quo_rem(hh);assert not remainder
    assert kk.degree()<=3 and mm.degree()<=5 and bb.degree()<=4 and qq.degree()==2
    assert AA==mm*kk-3*bb*bb/4-hh*hh*qq/4
    assert BB==(mm*mm*qq+bb**3-2*mm*bb*kk+hh*hh*(kk*kk-bb*qq))/4
    assert -f0==(mm*bb-hh*hh*kk)/2
    applications.append({'mask':mask,'normal_form_degrees':[int(f.degree()) for f in [hh,mm,bb,kk,qq]],
                         'k':list(map(str,kk.list()))})
report={'status':'PASS_UNIVERSAL_EUCLIDEAN_BISECTION_IDENTITY',
    'classification':'new exact normal form and symbolic verification, not a new seed family selection',
    'formulas':{'A':'m*k-3*b^2/4-h^2*q/4',
                'B':'(m^2*q+b^3-2*m*b*k+h^2*(k^2-b*q))/4',
                'X':'(b+h*W)/2','Y':'-(h*k+m*W)/2','cover':'W^2=q'},
    'identities':['trace equation','RR line','sixth-order cancellation','elliptic lift modulo W^2-q'],
    'applications':applications,
    'boundary':'Parent coefficients and17 sections remain fixed. The universal form does not supply a new determinant1092 parameter or resolve rational incidence.',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [OUT/'identity-protocol.json',OUT/'summary.json',Path(__file__)]}}
save('universal-identity.json',report)
print(report['status'],'all three normal forms verified',flush=True)
