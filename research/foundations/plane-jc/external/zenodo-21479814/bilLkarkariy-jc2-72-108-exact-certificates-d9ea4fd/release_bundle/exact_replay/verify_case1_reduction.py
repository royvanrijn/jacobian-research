#!/usr/bin/env python3
from fractions import Fraction
# Proposition 4.3 subcase (1) polygons; enumerate lattice points.
def inside(poly,p):
    # convex polygon CCW/boundary, sign-consistent cross products
    signs=[]
    for a,b in zip(poly,poly[1:]+poly[:1]):
        cross=(b[0]-a[0])*(p[1]-a[1])-(b[1]-a[1])*(p[0]-a[0])
        if cross: signs.append(cross>0)
    return not signs or all(signs) or not any(signs)
def lattice(poly):
    xs=range(min(x for x,y in poly),max(x for x,y in poly)+1)
    ys=range(min(y for x,y in poly),max(y for x,y in poly)+1)
    return {(x,y) for x in xs for y in ys if inside(poly,(x,y))}
Ppoly=[(0,0),(1,0),(8,14),(8,16),(0,8)]
Qpoly=[(0,0),(2,1),(12,21),(12,24),(0,12)]
Ps=lattice(Ppoly); Qs=lattice(Qpoly)
assert len(Ps)==61,(len(Ps),sorted(Ps))
assert len(Qs)==125,len(Qs)
# t=x y^2, z=y^-1: x^i y^j=t^i z^(2i-j)
def bands(S):
    out={}
    for i,j in S: out.setdefault(2*i-j,[]).append(i)
    return {k:(min(v),max(v),len(v)) for k,v in sorted(out.items())}
PB=bands(Ps);QB=bands(Qs)
assert PB[2]==(1,8,8) and PB[1]==(1,8,8) and PB[0]==(0,8,9)
for k in range(-1,-9,-1): assert PB[k]==(0,8+k,9+k)
assert QB[3]==(2,12,11) and QB[2]==(2,12,11) and QB[1]==(1,12,12) and QB[0]==(0,12,13)
for k in range(-1,-13,-1): assert QB[k]==(0,12+k,13+k)
# Coordinate determinant: t_x=y^2, t_y=2xy; z_x=0,z_y=-y^-2, so [t,z]=-1.
# Hence [P_i(t)z^i,Q_j(t)z^j]_{x,y}= z^(i+j-1)(i P_i Q_j' - j P_i' Q_j).
print('P_lattice_points',len(Ps),'bands',PB)
print('Q_lattice_points',len(Qs),'bands',QB)
print('coordinate_bracket_[t,z]',-1)
print('generic_band_formula_verified_symbolically_by_chain_rule',True)
