"""Rational fiber decoupling, flags and translation-kernel perturbations."""
import sympy as S
from rational_fiber import graph_shift
def nonlinear_fiber(v,q,mu,T):
    v=S.Matrix(v);q=S.Matrix(q);B=v.jacobian(list(q))
    data=graph_shift(B,q,mu,T);data.update(v=v,B=B)
    return data
def multi_direction(Q,R,M,p,B):
    return S.Matrix(M)*(Q*S.Matrix(p)+S.Matrix(B)*S.Matrix(R))
def split_cubic():
    u,x,y,z=S.symbols("u x y z");ps=[u,x,y,z];p=S.Matrix(ps)
    Q=u*x;B=S.eye(4)[:,2:4];R=S.Matrix([u**3,x**3]);M=S.diag(1,2,3,4)
    V=multi_direction(Q,R,M,p,B)
    return dict(ps=ps,p=p,Q=Q,B=B,R=R,M=M,V=V)
