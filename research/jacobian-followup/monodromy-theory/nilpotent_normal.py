"""Linear feedback, initial normal forms and exact cubic controls."""
import sympy as S
def linear_feedback(B,q,J):
    B=S.Matrix(B);q=S.Matrix(q)
    return B+B*q*S.Matrix([[S.diff(J,v)/J for v in q]])
def moments(B,q,J):
    B=S.Matrix(B);q=S.Matrix(q)
    row=S.Matrix([[S.diff(J,v) for v in q]])
    return [(row*B**j*q)[0] for j in range(1,len(q))]
def cubic_form(ps):
    u,x,y,z=ps
    return S.Matrix([u**3,u*u*x,u*u*y+u*x*x,u*u*z+x*y*y])
def cubic_control(ps,s):
    u,x,y,z=ps;lam=S.Integer(s+3);A=S.diag(1,lam,lam,lam)
    V=A*cubic_form(ps);F=u**s;h=lam*u*u
    X=S.Matrix([-F*h,*list(F*(V-h*S.Matrix(ps)))])
    P=X[1:,:].jacobian(ps)-S.Matrix(ps)*X[:1,:].jacobian(ps)
    return dict(A=A,V=V,F=F,h=h,lam=lam,X=X,P=P)
