#!/usr/bin/env python3
"""Exact symbolic check of the Laurent-coordinate reduction for JC2 case (8,28), case (2)."""
import sympy as sp

t,z=sp.symbols('t z')
A,B,C,D,E,F,G=[sp.Function(n)(t) for n in 'ABCDEFG']
P=A*z**2+B*z+C
Q=D*z**3+E*z**2+F*z+G
Jtz=sp.expand(sp.diff(P,t)*sp.diff(Q,z)-sp.diff(P,z)*sp.diff(Q,t))
coeff=[sp.expand(Jtz.coeff(z,k)) for k in range(5)]
expected=[
    -B*sp.diff(G,t)+sp.diff(C,t)*F,
    -2*A*sp.diff(G,t)-B*sp.diff(F,t)+2*sp.diff(C,t)*E+sp.diff(B,t)*F,
    -2*A*sp.diff(F,t)-B*sp.diff(E,t)+3*sp.diff(C,t)*D+2*sp.diff(B,t)*E+sp.diff(A,t)*F,
    -2*A*sp.diff(E,t)-B*sp.diff(D,t)+3*sp.diff(B,t)*D+2*sp.diff(A,t)*E,
    -2*A*sp.diff(D,t)+3*sp.diff(A,t)*D,
]
assert all(sp.expand(x-y)==0 for x,y in zip(coeff,expected))
# t=xy^2,z=y^-1 has Jacobian -1; x^2=t^2 z^4.
x,y=sp.symbols('x y')
txy=x*y**2; zxy=y**-1
jac_tz=sp.simplify(sp.diff(txy,x)*sp.diff(zxy,y)-sp.diff(txy,y)*sp.diff(zxy,x))
assert jac_tz==-1
assert sp.simplify(x**2-txy**2*zxy**4)==0
# lattice support counts for the two stated Newton polygons
def inside_convex(pt, verts):
    # boundary-inclusive, vertices counterclockwise or clockwise
    signs=[]
    for (x1,y1),(x2,y2) in zip(verts,verts[1:]+verts[:1]):
        c=(x2-x1)*(pt[1]-y1)-(y2-y1)*(pt[0]-x1)
        if c: signs.append(1 if c>0 else -1)
    return not signs or all(s==signs[0] for s in signs)
Pverts=[(0,0),(1,0),(8,14),(8,16)]
Qverts=[(0,0),(2,1),(12,21),(12,24)]
Ppts=[(i,j) for i in range(9) for j in range(17) if inside_convex((i,j),Pverts)]
Qpts=[(i,j) for i in range(13) for j in range(25) if inside_convex((i,j),Qverts)]
assert len(Ppts)==25, len(Ppts)
assert len(Qpts)==47, len(Qpts)
# every P point is in one of j=2i-2,2i-1,2i bands (plus origin), Q in four bands.
assert all((i,j)==(0,0) or j in (2*i-2,2*i-1,2*i) for i,j in Ppts)
assert all((i,j)==(0,0) or j in (2*i-3,2*i-2,2*i-1,2*i) for i,j in Qpts)
print('PASS: coordinate Jacobian = -1')
print('PASS: x^2 = t^2 z^4')
print('PASS: five coefficient identities verified exactly')
print(f'PASS: lattice supports contain {len(Ppts)} and {len(Qpts)} points')
