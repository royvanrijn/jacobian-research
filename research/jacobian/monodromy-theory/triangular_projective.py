"""Exact homogeneous triangular projective data and cyclic-coordinate identities."""
import sympy as S
def architecture(ps,degree,weights,curvatures):
    r=len(ps);n=r-1;u=ps[0];xs=ps[1:]
    if degree<1 or len(weights)!=n or len(curvatures)!=n:raise ValueError("invalid architecture")
    qs=[S.sympify(q) for q in curvatures]
    for i,q in enumerate(qs):
        if any(q.has(x) for x in xs[i:]):raise ValueError("curvature must use earlier coordinates")
        if i and S.Poly(q,*xs[:i]).total_degree()>degree:raise ValueError("homogenization degree too small")
    v=S.Matrix([a*x+q for a,x,q in zip(weights,xs,qs)])
    V=S.Matrix([u**degree,*[S.expand(u**degree*e.subs({x:x/u for x in xs},simultaneous=True)) for e in v]])
    B=v.jacobian(xs);W=v-S.Matrix(xs)
    A=V.jacobian(ps).subs(u,1);vv=V.subs(u,1)
    return dict(parameters=ps,degree=degree,weights=weights,Q=qs,v=v,V=V,A=A,chart_V=vv,B=B,W=W)
def resolvent(data,T):
    d=data["degree"];xs=data["parameters"][1:];B=data["B"];W=data["W"]
    high=(S.Matrix(xs)+T*(T*S.eye(len(xs))-B).inv()*W)/(T-d)
    return S.Matrix([1/(T-d),*list(high)]).applyfunc(S.cancel)
def last_coordinate(data,k):
    # Coordinates after k must already be diagonal. k is a zero-based projective index.
    xs=data["parameters"][1:];B=data["B"];W=data["W"];aa=data["weights"];c=aa[k]
    for j in range(k+1,len(xs)):
        if data["Q"][j]!=0:raise ValueError("later projective coordinates must be diagonal")
    B0=B[:k,:k];W0=W[:k,:];Q=data["Q"][k]
    grad=S.Matrix([[S.diff(Q,x) for x in xs[:k]]]) if k else S.zeros(1,0)
    chi0=(c*S.eye(k)-B0).det() if k else S.Integer(1)
    suffix=S.prod(c-a for a in aa[k+1:])
    adjterm=(grad*(c*S.eye(k)-B0).adjugate()*W0)[0] if k else S.Integer(0)
    L=S.expand(c*suffix*((c-1)*chi0*xs[k]+chi0*Q+adjterm))
    out=dict(c=c,chi0=chi0,suffix=suffix,L=L,K=c*(c-1)*chi0*suffix)
    if chi0!=0 and c!=1:
        row=(grad*(c*S.eye(k)-B0).inv()).applyfunc(S.cancel) if k else S.zeros(1,0)
        Sshift=S.cancel((Q+(row*W0)[0])/ (c-1)) if k else S.cancel(Q/(c-1))
        ell=S.Matrix([[*list(row),1,*([0]*(len(xs)-k-1))]])
        out.update(shift=Sshift,q=xs[k]+Sshift,ell=ell)
    return out
def rigidity_residuals(v,xs,q,c):
    B=v.jacobian(xs);W=v-S.Matrix(xs);dq=S.Matrix([[S.diff(q,x) for x in xs]])
    H=S.hessian(q,xs);R=dq*B-c*dq;Z=(dq*W)[0]-(c-1)*q
    DR=R.jacobian(xs)
    return dict(B=B,W=W,H=H,R=R,Z=Z,
     symmetric_identity=H*B-B.T*H+DR-DR.T,
     differentiated_identity=H*W-S.Matrix([S.diff(Z,x) for x in xs])+R.T)
