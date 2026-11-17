"""Exact projective rank-one curvature and universal companion curl data."""
import sympy as S
from primitive_multiplier import feedback_form
def homogenize(v,u,xs,d):
    v=S.Matrix(v)
    if d<1 or any(S.Poly(e,*xs).total_degree()>d for e in v):raise ValueError("invalid homogenization degree")
    return S.Matrix([u**d,*[S.expand(u**d*e.subs({x:x/u for x in xs},simultaneous=True)) for e in v]])
def companion():
    u,x,y,z,lam,tau,sig,delta=S.symbols("u x y z lam tau sig delta")
    ps=[u,x,y,z];v=S.Matrix([y,z-sig*x+tau*y+x*x,delta*x])
    V=homogenize(v,u,[x,y,z],2)
    f=feedback_form(V,ps,lam*u)
    D=S.cancel(f["detC"]/(delta*u**8)).subs(u,1)
    eta=f["eta"].subs(u,1)
    curl=S.cancel(D**2*(S.diff(eta[2],x)-S.diff(eta[1],y)))
    return dict(ps=ps,parameters=[lam,tau,sig,delta],v=v,V=V,feedback=f,D=D,eta=eta,curl_polynomial=curl)
def companion_frame(M,c,ell):
    delta=M.det();tau=S.trace(M);sig=(tau*tau-S.trace(M*M))/2
    frame=S.Matrix.vstack(ell,ell*M,delta*ell*M.inv())
    return dict(delta=delta,tau=tau,sigma=S.factor(sig),a=(ell*M*c)[0],frame=frame,
                transformed_M=(frame*M*frame.inv()).applyfunc(S.cancel),
                transformed_c=(frame*c).applyfunc(S.cancel))
