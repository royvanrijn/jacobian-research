"""Rational degree-one fiber graph and nonlinear critical-factor directions."""
import sympy as S
def graph_shift(C,q,mu,T):
    q=S.Matrix(q);C=S.Matrix(C)
    dT=S.Matrix([[S.diff(T,v) for v in q]])
    shift=(dT*(mu*S.eye(len(q))-C).inv()*q)[0]
    dS=S.Matrix([[S.diff(shift,v) for v in q]])
    defect=dS*C+dT-mu*dS
    return dict(S=S.cancel(shift),dS=dS,dT=dT,defect=defect.applyfunc(S.cancel))
def rank_one_direction(Q,R,M,p,b):
    p=S.Matrix(p);M=S.Matrix(M);b=S.Matrix(b)
    return Q*M*p+(M*b)*R
def explicit_cubic():
    x,y,w,z=S.symbols("x y w z");ps=[x,y,w,z];p=S.Matrix(ps)
    Q=x*x+y*y+w*w;R=x**3;b=S.Matrix([0,0,0,1])
    C=S.Matrix([[2,1,0],[-1,2,0],[0,0,2]]);M=S.diag(C,S.Integer(5))
    V=rank_one_direction(Q,R,M,p,b)
    return dict(ps=ps,p=p,Q=Q,R=R,b=b,C=C,M=M,V=V,mu=S.Integer(5))
