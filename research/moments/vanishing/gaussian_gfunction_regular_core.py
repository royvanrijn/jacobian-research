"""Exact finite Fuchsian certificates, not a Gaussian realization oracle."""
import sympy as s

t=s.Symbol('t')

def _domain(parameters):
    return s.QQ.frac_field(*parameters) if parameters else s.QQ

def fuchsian_certificate(coefficients,parameters=()):
    """Coefficients are [p_0,...,p_r] for sum p_i(t) D^i.

    QQ(parameter)-polynomials in t only. Retains every divisibility witness
    and a sufficient open specialization guard. Does not verify annihilation,
    minimality, Gaussian realization, or analytic convergence of a supplied jet.
    """
    domain=_domain(parameters)
    polys=[s.Poly(s.cancel(v),t,domain=domain) for v in coefficients]
    r=len(polys)-1
    if r<1 or polys[-1].is_zero:
        raise ValueError('nonzero positive differential order required')
    pr=polys[-1];degree=pr.degree();factors=s.factor_list(pr.as_expr(),t,domain=domain)[1]
    factors=[(s.Poly(q,t,domain=domain).monic(),int(m)) for q,m in factors]
    witnesses=[];scalars=[pr.LC()]
    all_polys=list(polys)
    for q,m in factors:
        scalars.append(s.discriminant(q.as_expr(),t))
        all_polys.append(q)
        for i,pi in enumerate(polys[:-1]):
            e=max(m-r+i,0)
            quot,rem=s.div(pi,q**e,domain=domain)
            if not rem.is_zero:
                raise ValueError('finite singularity violates Fuchs bound')
            witnesses.append(dict(factor=str(q.as_expr()),multiplicity=m,derivative_order=i,required_power=e,quotient=str(quot.as_expr())))
            all_polys.append(quot)
    for j,(q,m) in enumerate(factors):
        for qq,mm in factors[j+1:]:scalars.append(s.resultant(q.as_expr(),qq.as_expr(),t))
    infinity=[]
    for i,pi in enumerate(polys[:-1]):
        bound=degree-r+i
        if not pi.is_zero and pi.degree()>bound:
            raise ValueError('infinity violates Fuchs bound')
        infinity.append(dict(derivative_order=i,degree=None if pi.is_zero else int(pi.degree()),maximum_degree=int(bound)))
    for poly in all_polys:
        for c in poly.all_coeffs():scalars.append(s.fraction(s.cancel(c))[1])
    guards=[]
    if parameters:
        for v in scalars:
            num,den=s.fraction(s.cancel(v))
            for a in [num,den]:
                if a==0:raise ValueError('zero specialization guard')
                for q,m in s.factor_list(a,*parameters)[1]:
                    q=s.Poly(q,*parameters,domain=s.QQ).monic().as_expr()
                    if q not in guards:guards.append(q)
    guard=s.prod(guards)
    return dict(order=r,coefficients=[str(z.as_expr()) for z in polys],parameters=list(map(str,parameters)),finite_witnesses=witnesses,infinity_bounds=infinity,specialization_guard=str(s.factor(guard)),guard_factors=list(map(str,guards)),Fuchsian=True,Gaussian_realization_inferred=False,minimality_inferred=False,annihilation_inferred=False,convergence_inferred_from_jets=False)

def indicial_at(coefficients,point):
    """Exact local indicial polynomial, after verifying the Fuchs bounds."""
    fuchsian_certificate(coefficients)
    x=s.Symbol('rho');r=len(coefficients)-1;pr=s.sympify(coefficients[-1]);ans=s.ff(x,r)
    if point==s.oo:
        # A solution t^-rho has derivatives (-1)^i (rho)^rising_i.
        ans=s.Integer(0)
        degree=s.degree(pr,t)
        for i,pi in enumerate(coefficients):
            c=s.limit(s.sympify(pi)*t**(r-i)/t**degree,t,s.oo)
            ans+=c*(-1)**i*s.rf(x,i)
    else:
        for i,pi in enumerate(coefficients[:-1]):
            c=s.limit((t-point)**(r-i)*s.sympify(pi)/pr,t,point)
            ans+=c*s.ff(x,i)
    return s.factor(ans)

def derivative_recurrence(coefficients,series,order):
    """Finite residual only; never used to certify an all-order identity."""
    expr=sum(s.sympify(pi)*s.diff(series,t,i) for i,pi in enumerate(coefficients))
    return s.series(expr,t,0,order+1).removeO().expand()
