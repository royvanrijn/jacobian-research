"""Finite exact controls for the written critical-normalization theorem.

Jets never certify rank, Gaussian realizability, arithmetic spreading or
algebraicity. Kummer certificates below check supplied rational data only.
"""
import sympy as s

t=s.Symbol('t')

def jet(a,n):
    return s.series(a,t,0,n+1).removeO().expand()

def wick(poly,X):
    ans=0
    for powers,c in s.Poly(s.expand(poly),*X).terms():
        if any(k%2 for k in powers):
            continue
        ans+=c*s.prod(s.factorial2(k-1) if k else 1 for k in powers)
    return s.expand(ans)

def period(P,Q,X,n):
    return s.expand(sum(wick(Q*P**k,X)*t**k/s.factorial(k) for k in range(n+1)))

def critical(P,X,n):
    g=s.zeros(len(X),1)
    for _ in range(n):
        sub=dict(zip(X,g))
        g=s.Matrix([jet(t*s.diff(P,x).subs(sub,simultaneous=True),n) for x in X])
    sub=dict(zip(X,g))
    psi=jet(t*P.subs(sub,simultaneous=True)-(g.T*g)[0]/2,n)
    B=s.eye(len(X))-t*s.hessian(P,X).subs(sub,simultaneous=True)
    A=B-s.eye(len(X));L=s.zeros(len(X));power=s.eye(len(X))
    for k in range(n+1):
        L+=s.binomial(-s.Rational(1,2),k)*power
        power=(power*A).applyfunc(lambda a:jet(a,n))
    return g,psi,B.applyfunc(lambda a:jet(a,n)),L.applyfunc(lambda a:jet(a,n))

def normalized_via_saddle(P,Q,X,n):
    g,psi,B,L=critical(P,X,n)
    Y=s.symbols('y0:'+str(len(X)))
    gv=dict(zip(X,g));Yv=s.Matrix(Y)
    shifted=P.subs(dict(zip(X,g+Yv)),simultaneous=True)
    grad=s.Matrix([s.diff(P,x).subs(gv,simultaneous=True) for x in X])
    Hess=s.hessian(P,X).subs(gv,simultaneous=True)
    residual=s.expand(shifted-P.subs(gv,simultaneous=True)-(grad.T*Yv)[0]-(Yv.T*Hess*Yv)[0]/2)
    residual=jet(residual.subs(dict(zip(Y,L*Yv)),simultaneous=True),n)
    amp=jet(Q.subs(dict(zip(X,g+L*Yv)),simultaneous=True),n)
    expjet=s.Integer(1);power=s.Integer(1)
    for k in range(1,n+1):
        power=jet(power*t*residual,n)
        expjet+=power/s.factorial(k)
    return jet(s.det(L)*wick(jet(amp*expjet,n),Y),n)

def scaled_integral(a,scale,marker_denominator=1):
    b=s.Poly(s.expand(marker_denominator*a.subs(t,scale*t)),t)
    return all(c.is_Integer for c in b.all_coeffs())

def occupation(alpha,counts,gamma):
    if any(sum(a)<2 or any(x<0 for x in a) for a in alpha):
        raise ValueError('phase vertices must have degree at least two')
    if len(alpha)!=len(counts) or any(k<0 for k in counts):
        raise ValueError('invalid occupation')
    beta=tuple(gamma[i]+sum(k*a[i] for k,a in zip(counts,alpha)) for i in range(len(gamma)))
    if any(k%2 for k in beta):
        return s.Integer(0)
    return s.Rational(s.prod(s.factorial2(k-1) if k else 1 for k in beta),s.prod(s.factorial(k) for k in counts))

def kummer_certificate(action,log_derivative):
    """QQ rational input; exact residue-resultant and Kummer witness.

    The action is supplied; this does not recover a phase's exact action
    from jets or certify a Gaussian image. Reject nonrational residues.
    """
    R=s.cancel(action);r=s.cancel(log_derivative)
    for val in [R,r]:
        num,den=s.fraction(val)
        s.Poly(num,t,domain=s.QQ);s.Poly(den,t,domain=s.QQ)
        if den.subs(t,0)==0:
            raise ValueError('data must be regular at origin')
    if R.subs(t,0)!=0:
        raise ValueError('critical action must vanish at origin')
    a=s.cancel(r-s.diff(R,t));num,den=s.fraction(a)
    num=s.Poly(num,t,domain=s.QQ);den=s.Poly(den,t,domain=s.QQ).monic()
    # Recompute numerator after making the denominator monic.
    num=s.Poly(s.cancel(a*den.as_expr()),t,domain=s.QQ)
    if a==0:
        return dict(action=str(R),log_derivative=str(r),normalized_log_derivative='0',N=1,q='1',residue_resultant='1',residue_factors=[],Gaussian_realization_inferred=False)
    if num.degree()>=den.degree() or s.gcd(den,den.diff()).degree()>0:
        raise ValueError('normalized connection has irregular pole')
    z=s.Symbol('z')
    res=s.Poly(s.resultant(den.as_expr(),num.as_expr()-z*den.diff().as_expr(),t),z,domain=s.QQ).monic()
    facts=s.factor_list(res.as_expr(),z)[1]
    if any(s.degree(f,z)!=1 for f,m in facts):
        raise ValueError('normalized residues are not all rational')
    residues=[s.solve(f,z)[0] for f,m in facts]
    N=int(s.ilcm(*[int(c.q) for c in residues],1))
    q=s.Integer(1);parts=[]
    for c in residues:
        divisor=s.gcd(den,s.Poly(num.as_expr()-c*den.diff().as_expr(),t,domain=s.QQ)).monic()
        if divisor.degree()==0:
            continue
        factor=divisor.as_expr()/divisor.eval(0)
        q*=factor**int(N*c)
        parts.append(dict(residue=str(c),factor=str(s.expand(factor)),exponent=int(N*c)))
    q=s.cancel(q)
    if s.cancel(s.diff(q,t)/q-N*a)!=0 or q.subs(t,0)!=1:
        raise ValueError('Kummer witness failed')
    return dict(action=str(R),log_derivative=str(r),normalized_log_derivative=str(a),N=N,q=str(q),residue_resultant=str(res.as_expr()),residue_factors=parts,Gaussian_realization_inferred=False)

def mod_rational(a,p):
    num,den=s.fraction(s.cancel(a))
    n=s.Poly(num,t,modulus=p);d=s.Poly(den,t,modulus=p)
    if d.is_zero:
        raise ValueError('bad reduction')
    g=s.gcd(n,d);n=n.exquo(g);d=d.exquo(g)
    inv=pow(int(d.LC())%p,-1,p)
    return s.cancel((n.mul_ground(inv)).as_expr()/(d.mul_ground(inv)).as_expr())

def curvature(a,p):
    """D^p h/h for h'=a h, in characteristic p."""
    return mod_rational(a**p+s.diff(a,t,p-1),p)
