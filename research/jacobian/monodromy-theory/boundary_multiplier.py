"""Weighted boundary matrices for critical-coordinate scalar factors."""
import sympy as S
def weighted_matrix(Q,ell,xs,s):
    Q=S.Matrix(Q);ell=S.Matrix([list(ell)]);xx=S.Matrix(xs)
    return S.BlockMatrix([[S.Matrix([[(s+1)*(ell*xx)[0]]]),ell],
                         [s*Q,Q.jacobian(xs)]]).as_explicit()
def boundary_limit(V,ps,e,lam,G):
    V=S.Matrix(V);u=ps[0];xs=ps[1:]
    eta=S.Matrix([[e/u+S.diff(G,u)/G,*[S.diff(G,x)/G for x in xs]]])
    N=V.jacobian(ps)+V*eta-lam*u*S.eye(len(ps))
    T=S.diag(u,*([1]*len(xs)))
    rescaled=(T.inv()*N*T).applyfunc(S.cancel)
    return dict(N=N,rescaled=rescaled,limit=rescaled.subs(u,0).applyfunc(S.cancel))
def rank3_family(ps,s):
    u,x,y,z=ps;lam=S.Integer(s+2)
    V=S.Matrix([u*u,lam*u*x,lam*u*y+x*x,lam*u*z+x*y])
    F=u**s;h=lam*u
    X=S.Matrix([-F*h,*list(F*(V-h*S.Matrix(ps)))])
    P=X[1:,:].jacobian(ps)-S.Matrix(ps)*X[:1,:].jacobian(ps)
    return dict(V=V,F=F,h=h,X=X,P=P,lam=lam)
