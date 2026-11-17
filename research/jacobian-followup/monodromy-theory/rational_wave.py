"""Exact reconstruction in the complete z-dependent rational wave chart."""
import sympy as S
def rational_chart(r,s,R,w=None,z=None,k=None):
    if not isinstance(r,int) or not isinstance(s,int) or r<2 or s<0:
        raise ValueError("r>=2 and s>=0 must be integers")
    w=S.Symbol("w") if w is None else w
    z=S.Symbol("z") if z is None else z
    k=S.Symbol("kappa") if k is None else k
    minimal=r*(r-1)+2*k*(r-1)*(s+1)+k*k*(s+1)*(s+2)
    def red(e):
        num,den=S.fraction(S.cancel(e))
        return S.rem(S.rem(num,minimal,k)*S.invert(S.rem(den,minimal,k),minimal,k),minimal,k).expand()
    def C(a,b):return a*(a-1)+2*k*(a-1)*(b+1)+k*k*(b+1)*(b+2)
    def laurent_terms(e):
        out={}
        for term in S.Add.make_args(S.expand(e)):
            a=int(term.as_powers_dict().get(w,0))
            co=S.cancel(term/w**a)
            if co.has(w):raise ValueError("R must be Laurent polynomial in w")
            out[a]=red(out.get(a,0)+co)
        return {a:co for a,co in out.items() if co!=0}
    terms=laurent_terms(R)
    if terms and min(terms)<0 and r+(s+2)*min(terms)<2:
        raise ValueError("pole exceeds the complete integrability bound")
    fs={s:w**r,s+1:S.Integer(0)}
    for b in range(s-1,-1,-1):
        F=2*(b+1)*w*w*S.diff(R*fs[b+1]/w,w)
        F+=2*k*(b+1)*(b+2)*R*fs[b+1]
        F+=(b+2)*(b+1)*R*R*fs.get(b+2,0)
        fs[b]=red(sum(-co*red(1/C(a,b))*w**a for a,co in laurent_terms(F).items()))
    lam=red(sum(fs[b]*z**b for b in range(s+1)))
    rho=(k*z+R)/w
    x1z=red(-S.diff(lam,w)-S.diff(rho*lam,z))
    x1w=red(S.diff(rho*lam,w)+S.diff(rho*rho*lam,z))
    X1=red(S.integrate(x1z,z)+S.integrate(red(x1w.subs(z,0)),w))
    X2=red(S.integrate(w*S.diff(X1,z)+lam,z)+S.integrate(red((w*S.diff(X1,w)-rho*lam).subs(z,0)),w))
    X3=red(S.integrate(z*S.diff(X1,z)+rho*lam,z)+S.integrate(red((-rho*rho*lam).subs(z,0)),w))
    X=S.Matrix([X1,X2,X3]).applyfunc(S.cancel)
    V=red(X2-w*X1);W=red(X3-z*X1)
    rest={a:co for a,co in terms.items() if a not in [0,1]}
    j=max(rest) if rest else None
    if any(a in terms for a in [0,1]) and rest:
        raise ValueError("strip affine R before using the degree formula")
    N=r+s+1 if not rest else (s+2)*(r+(s+1)*j)
    from math import gcd
    e=(s+1)//gcd(s+1,r-1)
    count=gcd(s+1,r-1)
    group=("S" if (r+s+1)%2==0 else "A") if not rest else ("S" if r%2==0 and s%2==1 else "A")
    return dict(r=r,s=s,R=R,w=w,z=z,k=k,minimal=minimal,red=red,C=C,
      coefficients=fs,lambda_polynomial=lam,rho=rho,X=X,V=V,W=W,
      j=j,N=N,geometric_group=group+str(N),boundary_ramification_index=e,boundary_branch_count=count)
