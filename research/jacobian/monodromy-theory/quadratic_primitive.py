"""Quadratic direction closure: translation and zero-stratum limits."""
import sympy as S
from unmatched_quadratic import canonical_forms
def feedback(V,ps,F,h):
    V=S.Matrix(V)
    return V.jacobian(ps)+V*S.Matrix([[S.diff(F,v)/F for v in ps]])-h*S.eye(len(ps))
def trace_equation(V,ps,F,R):
    V=S.Matrix(V)
    return S.expand(sum(V[i]*S.diff(F,ps[i]) for i in range(len(ps)))+
                    S.trace(V.jacobian(ps))*F-R*F)
def normal_limit(V,ps,F,h,c,eps):
    normal=ps[:c];sub={v:eps*v for v in normal}
    N=feedback(V,ps,F,h).subs(sub,simultaneous=True)
    T=S.diag(*([eps]*c+[1]*(len(ps)-c)))
    gauge=(T.inv()*N*T).applyfunc(S.cancel)
    return dict(gauge=gauge,limit=gauge.applyfunc(lambda v:S.limit(v,eps,0)))
def lower_rank_family(ps,s,name):
    u,x,y,z=ps;lam=S.Integer(s+2);A=S.diag(1,lam,lam,lam)
    V=A*canonical_forms(ps)[name];F=u**s;h=lam*u
    X=S.Matrix([-F*h,*list(F*(V-h*S.Matrix(ps)))])
    P=X[1:,:].jacobian(ps)-S.Matrix(ps)*X[:1,:].jacobian(ps)
    return dict(V=V,F=F,h=h,lam=lam,X=X,P=P,A=A)
