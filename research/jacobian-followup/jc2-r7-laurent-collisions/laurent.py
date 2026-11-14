"""Canonical polynomially scaled Laurent gap-21 family and collision inputs."""
from pathlib import Path
import json
import sys
import sympy as S0

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE.parent/'jc2-r7-image-family'))
import family as F

z,t,v,w,V,U,S=S0.symbols('z t v w V U S')

def scaled(poly,degree):
    p=S0.Poly(S0.expand(poly.subs(F.laurent).subs(F.t,t/F.x)*F.x**degree),F.x,t)
    assert all(m[0]%3==0 for m,c in p.terms())
    return S0.expand(sum(c*z**(m[0]//3)*t**m[1] for m,c in p.terms()))

P,Q=scaled(F.p,6),scaled(F.q,10)
hh=[S0.Integer(1),v]
for n in range(2,10):hh.append(S0.expand(v*hh[-1]-w*hh[-2]))

def divided(poly):
    p=S0.Poly(poly,t)
    return S0.expand(sum(p.nth(i)*hh[i-1] for i in range(1,p.degree()+1)))

A,B=divided(P),divided(Q)
sub=S0.subresultants(A,B,w)
H=S0.cancel(81*sub[-1]/v**10)
linear=S0.cancel(-12*sub[-2]/v**3)
k=S0.Poly(linear,w).nth(1)
n=S0.Poly(linear,w).nth(0)
J=S0.expand(S0.diff(A,v)*S0.diff(B,w)-S0.diff(A,w)*S0.diff(B,v))
C1=3*z+2
C2=15*z-4
W=25*z*z+5*z+4
mods=json.loads((HERE/'exceptional_moduli.json').read_text())
Fproj=S0.sympify(mods['projection'],locals={'z':z})
Tangent=S0.sympify(mods['tangent'],locals={'z':z})
point=json.loads((HERE/'triple_point_certificate.json').read_text())

def univariate(rows):
    return S0.Add(*(S0.Rational(c)*z**m[0] for m,c in rows))

D=univariate(point['modulus'])
triple=[]
cube=t**3-S*t**2+U*t-V

def primitive(f,variables):
    return S0.Poly(f,*variables).clear_denoms()[1].primitive()[1].as_expr()

for poly in [P,Q]:
    rem=S0.Poly(S0.rem(poly,cube,t),t)
    triple.extend(primitive(rem.nth(i),(V,U,S,z)) for i in [2,1])
Vgraph=S0.solve(triple[0],V)[0]
first=primitive(S0.together(triple[1].subs(V,Vgraph)).as_numer_denom()[0],(U,S,z))
reduced=[first]
for f in triple[2:]:
    raw=primitive(S0.together(f.subs(V,Vgraph)).as_numer_denom()[0],(U,S,z))
    reduced.append(primitive(S0.rem(raw,first,U),(U,S,z)))
