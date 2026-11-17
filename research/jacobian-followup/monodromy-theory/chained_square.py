"""Canonical chained-square conjugacy and exact trace obstructions."""
import sympy as S
from unmatched_quadratic import canonical_forms
def direction(A,ps):
    return S.Matrix(A)*canonical_forms(ps)["A"]
def feedback(V,ps,eta,lam):
    return S.Matrix(V).jacobian(ps)+S.Matrix(V)*S.Matrix([list(eta)])-lam*ps[0]*S.eye(4)
def shift_matrices(t):
    C=S.eye(4);C[1,0]=t
    B=S.Matrix([[1,0,0,0],[-t,1,0,0],[t*t,-2*t,1,0],[0,0,0,1]])
    return C,B
def rank3_family(ps,s):
    u,x,y,z=ps;lam=S.Integer(s+2)
    V=direction(S.diag(1,lam,lam,lam),ps);F=u**s
    X=S.Matrix([-lam*u*F,*list(F*(V-lam*u*S.Matrix(ps)))])
    P=X[1:,:].jacobian(ps)-S.Matrix(ps)*X[:1,:].jacobian(ps)
    return dict(V=V,F=F,lam=lam,X=X,P=P)
