"""Coordinate-free tensor wave and necessary highest-direction obstruction."""
import sympy as S
def wave_tensor(A,b,c,w,z):
    return S.expand(S.diff(A*b*b,w,2)+2*S.diff(A*b*c,w,z)+S.diff(A*c*c,z,2))
def homogeneous_part(f,w,z):
    poly=S.Poly(S.expand(f),w,z)
    if poly.is_zero:return S.Integer(0),-1
    degree=poly.total_degree()
    return S.expand(sum(co*w**i*z**j for (i,j),co in poly.terms() if i+j==degree)),degree
def leading_direction(b,c,w,z):
    bt,bd=homogeneous_part(b,w,z);ct,cd=homogeneous_part(c,w,z)
    D=max(bd,cd)
    if bd<D:bt=S.Integer(0)
    if cd<D:ct=S.Integer(0)
    common=S.gcd(bt,ct)
    if common==0:raise ValueError("nonzero direction required")
    B=S.cancel(bt/common);C=S.cancel(ct/common)
    degree=max(S.Poly(B,w,z).total_degree(),S.Poly(C,w,z).total_degree())
    if degree>=2:
        return dict(rejected=True,reason="primitive leading projective degree>=2",
                    D=D,common=common,primitive_degree=degree,B=B,C=C)
    if degree==0:
        derivative=S.expand(B*S.diff(common,w)+C*S.diff(common,z))
        return dict(rejected=derivative!=0,reason="constant leading direction requires transverse common factor",
                    D=D,common=common,primitive_degree=0,B=B,C=C,directional_derivative=derivative)
    M=S.Matrix([[S.diff(B,w),S.diff(B,z)],[S.diff(C,w),S.diff(C,z)]])
    if S.expand(S.trace(M)**2-4*M.det())==0:
        return dict(rejected=True,reason="repeated projective fixed direction",
                    D=D,common=common,primitive_degree=degree,B=B,C=C)
    derivative=S.expand(B*S.diff(common,w)+C*S.diff(common,z))
    quotient=S.cancel(derivative/common)
    if quotient.has(w,z):
        return dict(rejected=True,reason="leading common factor is not supported on eigenlines",
                    D=D,common=common,primitive_degree=degree,B=B,C=C)
    return dict(rejected=False,reason="passes necessary leading-direction tests; no existence claim",
                D=D,common=common,primitive_degree=degree,B=B,C=C)
def conjugate_affine_pair(b,c,T,shift,w,z):
    q=S.Matrix([w,z]);old=T*q+shift
    vector=S.Matrix([b,c]).subs({w:old[0],z:old[1]},simultaneous=True)
    return (T.inv()*vector).applyfunc(S.expand)
