"""Polynomial one-parameter wave solutions; every nonconstant h sweep is finite."""
import sympy as S
def polynomial_chart(h,top,w=None,z=None):
    w=S.Symbol("w") if w is None else w
    z=S.Symbol("z") if z is None else z
    if S.degree(h,w)<1:raise ValueError("h must be nonconstant")
    s=len(top)-1
    # top is a list of pairs (slope,constant), one for each z coefficient.
    pairs=top
    if s<0 or pairs[s]==(0,0):raise ValueError("nonzero top coefficient required")
    fs={s:pairs[s][0]*w+pairs[s][1],s+1:S.Integer(0)}
    for b in range(s-1,-1,-1):
        forcing=-2*(b+1)*S.diff(h*fs[b+1],w)-(b+1)*(b+2)*h*h*fs.get(b+2,0)
        fs[b]=S.integrate(S.integrate(forcing,w),w)+pairs[b][0]*w+pairs[b][1]
    lam=S.expand(sum(fs[b]*z**b for b in range(s+1)))
    rho=h
    x1z=-S.diff(lam,w)-S.diff(rho*lam,z)
    x1w=S.diff(rho*lam,w)+S.diff(rho*rho*lam,z)
    X1=S.expand(S.integrate(x1z,z)+S.integrate(x1w.subs(z,0),w))
    X2=S.expand(S.integrate(w*S.diff(X1,z)+lam,z)+S.integrate((w*S.diff(X1,w)-rho*lam).subs(z,0),w))
    X3=S.expand(S.integrate(z*S.diff(X1,z)+rho*lam,z)+S.integrate((-rho*rho*lam).subs(z,0),w))
    X=S.Matrix([X1,X2,X3]);d=int(S.degree(h,w));j=d+1;eps=int(pairs[s][0]!=0)
    n=s*j+eps;a=n+j-1
    degree=(s+1)*((s+2)*j-1) if eps==0 else (s+2)*((s+1)*j+1)
    exceptional=False;m=None
    if s==0 and eps==1:
        center=-S.sympify(pairs[0][1])/S.sympify(pairs[0][0])
        v=S.Symbol("v")
        shifted=S.expand(h.subs(w,v+center))
        exceptional=S.expand(shifted-shifted.subs(v,-v))==0
        if exceptional:m=d//2
    if exceptional:
        group="SL2(F3) on8 letters" if m==1 else f"2^{2*m+1}:A{2*m+2} on{4*m+4} letters"
    else:group="A"+str(degree)
    return dict(h=h,lambda_polynomial=lam,coefficients=fs,X=X,w=w,z=z,
       d=d,s=s,epsilon=eps,j=j,n=n,weight_gamma=a,
       source_weights=(a,1,j),target_weights=(a,a+1,a+j),
       N=degree,delta=2 if exceptional else 1,geometric_group=group)
