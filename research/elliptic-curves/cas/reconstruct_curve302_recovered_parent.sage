"""Reconstruct the normalized MW17 parent from its quartic and four conics.

Frozen exact construction, not a new search. One worker, 120 seconds.
The separate verifier proves rank and saturation and checks Frobenius.
"""
import json
import runpy
import signal
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing, EllipticCurve, matrix, vector
from sage.schemes.elliptic_curves.jacobian import Jacobian

signal.alarm(120)
ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'artifacts/generated-results/elliptic-curves'
data = json.loads((BASE/'curve302_recovered_mw17_construction_inputs_v1.json').read_text())
d = data['quartic']
cs = data['initial_conics']
R0 = PolynomialRing(QQ, 't')
t = R0.gen()
K = R0.fraction_field()
R = PolynomialRing(K, names=('X','Y','Z'))
X,Y,Z = R.gens()


def decode(records, W=0):
    return sum((K(row['coefficient'])*X**row['powers'][0]*Y**row['powers'][1]*Z**row['powers'][2]*W**row['powers'][3] for row in records),R(0))


L = decode(d['residual_line'])
C = decode(d['plane_cubic'])
f = C+t*decode(d['quartic_W_quotient'],t*L)
assert decode(d['quartic'],t*L) == L*f
O = vector(K,[1,1,0])
assert f(*O) == 0
sections, words = [], []
for xy, uv, word in zip(d['plane_points'],d['line_directions'],d['polarization']['line_section_words']):
    x,y,u,v = map(K,xy+uv)
    z,w = 1-t*L(u,v,0),t*L(x,y,1)
    sections.append(vector(K,[x*z+u*w,y*z+v*w,z]))
    words.append(word)
ellz = K(L.monomial_coefficient(Z))
assert L == X-Y+ellz*Z
s = K(cs['s'])
a,b,c = map(K,cs['A'])
dd,e,ff = map(K,cs['B'])
r = 1-t*ellz*s
wa = (r+(a+b)*t)*t
xa = r*r+b*r*t-c*t*t
ya = xa-r*(r+(a+b)*t)
wb = (r+(dd+e)*t)*t
xb = e*r*t-ff*t*t
yb = xb-r*(r+(dd+e)*t)
sections += [vector(K,[xa,ya,s*wa]),vector(K,[xb,yb,s*wb])]
words += [d['polarization']['P_word'],d['polarization']['Q_word']]
assert matrix(ZZ,words).rank() == 15
for conic in data['additional_conics']:
    assert len(conic['solutions']) == 1
    sol = conic['solutions'][0]
    a,b,c,dd = map(K,sol['plane'])
    assert a == 1 and b != -1
    def q(p):
        return sum(K(row['coefficient'])*p[0]**row['powers'][0]*p[1]**row['powers'][1]*p[2]**row['powers'][2] for row in sol['conic'])
    base = vector(K,[ellz-c,b+1,0])
    direction = vector(K,[1+dd*t,0,-t*(b+1)])
    assert q(base) == 0
    p = q(direction)*base-(q(base+direction)-q(direction))*direction
    assert q(p) == 0
    yy,zz,ww = p
    xx = -b*yy-c*zz-dd*ww
    assert ww == t*(xx-yy+ellz*zz)
    sections.append(vector(K,[xx,yy,zz]))
    words.append(conic['word'])
assert len(sections) == 19 and all(f(*P)==0 for P in sections)
W = matrix(ZZ,words)
indices = data['basis_indices']
assert abs(W.matrix_from_rows(indices).det()) == 1

# Non-flex plane-cubic conversion with the prescribed origin. The tangent
# construction and Cremona identity are verified as polynomial identities.
J = Jacobian(f)
V = PolynomialRing(K,'v')
v = V.gen()
dy,dz = [f.derivative(z)(*O) for z in [Y,Z]]
direction = vector(K,[0,1,-dy/dz])
restriction = V(f(*(O+v*direction)))
assert restriction[0] == restriction[1] == 0 and restriction[3]
P2 = O-restriction[2]/restriction[3]*direction
assert f(*P2) == 0 and P2 != O
dy,dz = [f.derivative(z)(*P2) for z in [Y,Z]]
P3 = vector(K,[0,1,-dy/dz])
M = matrix(K,[O,P2,P3]).transpose()
assert M.det()
F2 = R(M.act_on_polynomial(f))
pulled = R(F2(X*X,Y*Z,X*Z))
assert all(e[0]>=2 and e[2]>=1 for e in pulled.dict())
F3 = R({(e[0]-2,e[1],e[2]-1):c for e,c in pulled.dict().items()})
assert F3*X*X*Z == pulled
aa,bb = map(K,[F3.monomial_coefficient(X**3),F3.monomial_coefficient(Y*Y*Z)])
WW = R(F3(-X,Y/bb,aa*bb*Z)/aa)
intermediate = EllipticCurve(WW(X,Y,1))
assert WW == intermediate.defining_polynomial()(X,Y,Z)
MI = M.inverse()
iso = intermediate.isomorphism_to(J)
u,r,s,tau = map(K,[-6,15,-3,-108])
normalized = J.change_weierstrass_model([u,r,s,tau])
E,basis,t0 = runpy.run_path(str(ROOT/'elliptic-curves/cas/load_curve302_recovered_parent.sage'))['load_curve302_recovered_parent']()
assert normalized == E and t0 == 0
for idx, expected in zip(indices,basis):
    v,w,z = MI*sections[idx]
    raw = iso(intermediate(-aa*bb*v/z,aa*bb*bb*v*w/z**2))
    actual = E((raw[0]-r)/u**2,(raw[1]-s*(raw[0]-r)-tau)/u**3)
    assert actual == expected
print('PASS quartic pencil and four conics reconstruct the exact normalized parent and all17 basis sections')
