"""Exact unordered collisions for the polynomial gap-21 coefficient family."""
from pathlib import Path
import sys
import sympy as S

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE.parent/'jc2-r7-image-family'))
import family as F

d,t=F.d,F.t
v,w,u,r=S.symbols('v w u r')
p=S.expand(F.p.subs(F.line))
q=S.expand(F.q.subs(F.line))
H=[S.Integer(1),v]
for n in range(2,10):H.append(S.expand(v*H[-1]-w*H[-2]))

def divided(poly):
    f=S.Poly(poly,t)
    return S.expand(sum(f.nth(i)*H[i-1] for i in range(1,f.degree()+1)))

A=divided(p)
B=divided(q)
A1=v*v-w+d
A2=v**3-3*v*w+d*v+1
H9=S.expand(-162*B.subs(w,v*v+d))
H3=2*v**3+9*d*v+6
rem=S.Poly(S.rem(q,t**3+d*t-u,t),t)
r2,r1,r0=[rem.nth(j) for j in [2,1,0]]
upair=-v*(v*v+d)
qpair=S.expand((r0-r2*(v*v+d)).subs(u,upair))
otherU=-1-upair
other_q=S.expand(rem.as_expr().subs({u:otherU,t:r})-qpair)
mixed=S.factor(S.resultant(r**3+d*r-otherU,other_q,r))
mixed_remainder=S.rem(mixed,H9,v)

T=(3761479876608*d**18+15898976317440*d**15+25767272444544*d**12
   +20621286822727*d**9+8479009232730*d**6+1615258153728*d**3+91424858112)
U=2187*d**9+3402*d**6+1512*d**3+196
M=243*d**3+64
C=3*d**3+2

def primitive(f,var):
    return S.Poly(f,d,var).clear_denoms()[1].primitive()[1].as_expr()

GCD_INPUTS={
    'mixed':(M,H9,primitive(mixed_remainder,v),v),
    'sameU':(U,primitive(r2,u),primitive(r1,u),u),
    'tangent':(T,H9,S.diff(H9,v),v),
}
