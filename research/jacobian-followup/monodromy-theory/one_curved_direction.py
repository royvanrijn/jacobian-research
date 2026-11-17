"""Homogeneous one-curved directions and their degree-independent arrow identities."""
import sympy as S

def architecture(ps,degree,weights,c,Q):
    u=ps[0];xs=ps[1:-1];t=ps[-1]
    if len(weights)!=len(xs):raise ValueError("one weight per earlier coordinate required")
    poly=S.Poly(Q,*xs) if xs else None
    if poly is not None and poly.total_degree()>degree:raise ValueError("projective curvature degree too large")
    if degree<1:raise ValueError("positive degree required")
    Qh=S.expand(u**degree*Q.subs({x:x/u for x in xs},simultaneous=True))
    V=S.Matrix([u**degree,*[a*u**(degree-1)*x for a,x in zip(weights,xs)],c*u**(degree-1)*t+Qh])
    A=V.jacobian(ps).subs(u,1);vv=V.subs(u,1)
    return dict(parameters=ps,degree=degree,weights=weights,c=c,Q=Q,Qh=Qh,V=V,A=A,v=vv)

def chart_identities(data):
    ps=data["parameters"];xs=ps[1:-1];t=ps[-1];aa=data["weights"];c=data["c"];Q=data["Q"];d=data["degree"]
    T=S.Symbol("T")
    prod=S.prod(c-a for a in aa)
    K=c*(c-1)*prod
    R=S.expand(c*(Q*prod+sum((a-1)*x*S.diff(Q,x)*S.prod(c-b for j,b in enumerate(aa) if j!=i) for i,(a,x) in enumerate(zip(aa,xs)))))
    chi=(T-d)*(T-c)*S.prod(T-a for a in aa)
    resolvent=S.Matrix([1/(T-d),*[a*x*(T-1)/((T-d)*(T-a)) for a,x in zip(aa,xs)],
     (c*t*(T-1)+T*Q+sum(T*(a-1)*x*S.diff(Q,x)/(T-a) for a,x in zip(aa,xs)))/((T-d)*(T-c))])
    adj=(chi*resolvent).applyfunc(S.cancel)
    determinant=S.prod(a*(a-1)*x for a,x in zip(aa,xs))*S.prod(aa[j]-aa[i] for i in range(len(aa)) for j in range(i+1,len(aa)))*(K*t+R)
    return dict(T=T,K=K,R=R,chi=chi,resolvent=resolvent,adjugate_vector=adj,krylov_determinant=determinant)
