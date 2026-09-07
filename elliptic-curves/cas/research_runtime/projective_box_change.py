"""Exact sufficient box equivalence and explicit new-coordinate witnesses.

A coordinate witness need not lift to a rational point of a quartic. This
module proves coordinate exposure only, never existence of a new direction.
"""
from fractions import Fraction as Q
from math import gcd,lcm

def primitive(values):
    values=tuple(map(Q,values));den=lcm(*(v.denominator for v in values))
    integers=tuple(int(v*den) for v in values);g=gcd(*integers)
    if not g:raise ValueError('zero projective tuple')
    sign=1 if next(v for v in integers if v)>0 else -1
    return tuple(sign*v//g for v in integers)

def transition(old,new):
    a,b,c,d=map(Q,old);e,f,g,h=map(Q,new)
    if a*d-b*c==0 or e*h-f*g==0:raise ValueError('singular horizontal map')
    return primitive((d*e-b*g,d*f-b*h,a*g-c*e,a*h-c*f))

def apply(matrix,pair):
    a,b,c,d=matrix;n,m=pair
    return primitive((a*n+b*m,c*n+d*m))

def classify(old,new,height):
    if type(height) is not int or height<2:raise ValueError('integer height at least2 required')
    T=transition(old,new);a,b,c,d=T
    same=(b==c==0 and abs(a)==abs(d)==1) or (a==d==0 and abs(b)==abs(c)==1)
    result={'old_coordinate_from_new':list(T),'height':height,
            'status':'EXACT_SAME_BOX' if same else 'BOX_RELATION_UNRESOLVED',
            'new_coordinate_witness':None,'old_coordinate_of_witness':None}
    if same:return result
    pairs=[(1,0),(0,1),(1,1),(-1,1)]
    for sign in (-1,1):
        pairs.extend([(sign*height,1),(sign,height),(sign*height,height-1),(sign*(height-1),height)])
    for pair in pairs:
        p=primitive(pair);q=apply(T,p)
        if max(map(abs,p))<=height and max(map(abs,q))>height:
            result.update(status='PROVED_NEW_COORDINATE_EXPOSURE',new_coordinate_witness=list(p),old_coordinate_of_witness=list(q));break
    return result
