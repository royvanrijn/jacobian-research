"""Last-coordinate adjugate multiplicity and the mixed cubic family."""
import sympy as S
from nilpotent_normal import cubic_form
def last_adjugate(V,ps):
    V=S.Matrix(V);head=ps[:-1];z=ps[-1]
    v=V[:-1,:];c=S.diff(V[-1],z);R=S.expand(V[-1]-c*z)
    D=c*S.eye(len(head))-v.jacobian(head)
    det=D.det(method="domain-ge");adj=D.adjugate();rowc=S.Matrix([[S.diff(c,q) for q in head]])
    rowR=S.Matrix([[S.diff(R,q) for q in head]])
    W=det*(c*z+R)+((z*rowc+rowR)*adj*v)[0]
    K=c*det+(rowc*adj*v)[0]
    return dict(c=c,R=R,D=D,W=S.expand(W),K=S.expand(K),v=v)
def direction(A,ps):
    return S.Matrix(A)*cubic_form(ps)
