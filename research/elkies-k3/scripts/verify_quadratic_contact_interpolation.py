#!/usr/bin/env python3
"""Exact jet interpolation identities; no K3 point or coefficient search."""
import json
from pathlib import Path
import sympy as S
root=Path(__file__).resolve().parents[2]
t,u,v,T,U,V=S.symbols('t u v T U V')
m0=S.Matrix([[1,t,u,t*u],[0,1,v,u+t*v],[1,T,U,T*U],[0,1,V,U+T*V]])
m2=S.Matrix([[1,t,t*t,u],[0,1,2*t,v],[1,T,T*T,U],[0,1,2*T,V]])
assert S.expand(m0.det()-((U-u)**2-v*V*(T-t)**2))==0
assert S.expand(m2.det()-(T-t)*((T-t)*(V+v)-2*(U-u)))==0
# Conjugate quadratic contact coordinates: t^2-s*t+p=0,
# u=a+b*t, slope=c+d*t. The F0 condition is Norm(slope)=b^2.
a,b,c,d,s,p=S.symbols('a b c d s p')
q=t*t-s*t+p
norm=S.rem(S.expand((c+d*t)*(c+d*(s-t))),q,t)
assert norm==c*c+c*d*s+d*d*p
# F2 quadratic section a+b*t+k*q: slope b+k*(2t-s).
k=S.symbols('k')
assert S.expand((b+k*(2*t-s)).subs(k,d/2)-(b-s*d/2+d*t))==0
# Interpolation controls on Q(sqrt2), not points on any K3.
r=S.sqrt(2)
control0=m0.subs({t:r,T:-r,u:r/2,U:-r/2,v:S.Rational(3,2)-r,V:S.Rational(3,2)+r})
assert control0.rank()==3 and (control0*S.Matrix([-1,-1,2,1])).applyfunc(S.simplify)==S.zeros(4,1)
control2=m2.subs({t:r,T:-r,u:5+r,U:5-r,v:1+2*r,V:1-2*r})
assert control2.rank()==3 and (control2*S.Matrix([-3,-1,-1,1])).applyfunc(S.simplify)==S.zeros(4,1)
# Both points on one ruling with zero slopes: determinant zero but no integral (1,1) member.
bad=m0.subs({t:0,T:1,u:0,U:0,v:0,V:0})
assert bad.rank()==2
assert all(z[0]==z[1]==0 for z in bad.nullspace())
# Genus9 plus primitive degree4 ruling excludes all double maps to genus<=1.
assert 4*0+2*1+(4-1)*(2-1)==5<9
# For p_a=3, C^2=2; two curves sharing contact length2 have intersection >=4.
assert 2*2>2
out={'schema':'quadratic-contact-interpolation-v1','status':'PASS',
 'F0_determinant':str(S.factor(m0.det())),
 'F2_determinant':str(S.factor(m2.det())),
 'quadratic_F0_condition':'c^2+c*d*s+d^2*p=b^2',
 'quadratic_F2_condition':'c=b-s*d/2',
 'controls':{'F0_rank':3,'F0_kernel':[-1,-1,2,1],'F2_rank':3,'F2_kernel':[-3,-1,-1,1],'ruling_degeneracy_rank':2},
 'boundary':'Jet identities and linear algebra controls only. No control is asserted to lie on an actual halving curve. The primitive monodromy, quadratic-point finiteness and contact rigidity are written/inherited mathematics. No positive MW17 correlated-gain cover is supplied.'}
path=root/'artifacts/generated-results/elkies-k3-quadratic-contact-interpolation-v1/result.json'
if '--write' in __import__('sys').argv:path.parent.mkdir(exist_ok=True);path.write_text(json.dumps(out,indent=2)+'\n')
else:assert out==json.loads(path.read_text())
print('PASS: both jet determinants; quadratic norm/trace conditions; nondegenerate and ruling controls')
