"""Exact degree-free scalar multiplier test for a fixed homogeneous direction."""
import sympy as S
def feedback_form(V,ps,h):
    V=S.Matrix(V);r=len(ps);A=V.jacobian(ps);T=S.Symbol('T')
    C=S.Matrix.hstack(*[A**j*V for j in range(r)])
    detC=S.factor(C.det(method='domain-ge'))
    if detC==0:raise ValueError('primitive direction is not cyclic')
    chi=S.Poly(A.charpoly(T).as_expr(),T);cs=chi.all_coeffs()
    div=[S.Integer(1)]
    for k in range(1,r+1):
        nk=(-1)**k*S.binomial(r,k)*h**k
        div.append(S.expand(nk-sum(cs[j]*div[k-j] for j in range(1,k+1))))
    mu=S.Matrix([[-div[j] for j in range(1,r+1)]])
    eta=(mu*C.inv()).applyfunc(S.factor)
    N=(A+V*eta-h*S.eye(r)).applyfunc(S.cancel)
    rho=S.Matrix([[S.diff(h,p) for p in ps]])+h*eta
    return dict(V=V,A=A,C=C,detC=detC,T=T,chi=chi.as_expr(),moments=mu,
                eta=eta,N=N,rho=rho)
def logarithmic_residues(data,ps,parameter=None):
    """Exact coefficient matching; consistency is retained, never solved generically away."""
    eta=data['eta'];D=data['detC']
    factors=[f for f,e in S.factor_list(D,*ps)[1]]
    columns=[]
    for f in factors:columns.append([S.cancel(D*S.diff(f,p)/f) for p in ps])
    rhs=[S.cancel(D*e) for e in eta]
    equations=[];B=[]
    # Match source monomials, keeping the multiplier parameter in the coefficients.
    polys=[S.Poly(e,*ps) for e in rhs]
    cpolys=[[S.Poly(v,*ps) for v in col] for col in columns]
    for i in range(len(ps)):
        monoms=set(polys[i].monoms())
        for col in cpolys:monoms.update(col[i].monoms())
        for mon in sorted(monoms):
            equations.append([col[i].coeff_monomial(mon) for col in cpolys])
            B.append(polys[i].coeff_monomial(mon))
    AA=S.Matrix(equations);bb=S.Matrix(B)
    if not factors:
        residues=S.zeros(0,1);residual=list(bb)
    else:
        pivots=AA.T.rref()[1]
        if len(pivots)!=len(factors):raise AssertionError('distinct logarithmic forms not independent')
        selected=list(pivots);square=AA.extract(selected,list(range(len(factors))))
        residues=(square.inv()*bb.extract(selected,[0])).applyfunc(S.factor)
        residual=list((AA*residues-bb).applyfunc(S.factor))
    nonzero=[S.fraction(S.cancel(e))[0] for e in residual if e!=0]
    if parameter is None:
        compatibility=not nonzero;condition=None
    else:
        condition=S.Poly(0,parameter)
        for e in nonzero:
            condition=S.gcd(condition,S.Poly(e,parameter))
        if not condition.is_zero:condition=condition.monic()
        compatibility=None
    return dict(factors=factors,residues=list(residues),residual=residual,
       compatible=compatibility,parameter_condition=condition,
       matrix_rows=AA.rows,matrix_columns=AA.cols)
def quadratic_shear():
    u,v,w,t,a,b,c,q0,q1,q2=S.symbols('u v w t a b c q0 q1 q2')
    ps=[u,v,w,t];Q=q0*v*v+q1*v*w+q2*w*w
    V=S.Matrix([u*u,a*u*v,b*u*w,c*u*t+Q]);A=V.jacobian(ps)
    kk=c*(c-1)*(c-a)*(c-b)
    rr=(c-2)*(c-a)*(c-b)*Q+(c-1)*(a*(c-b)*v*S.diff(Q,v)+b*(c-a)*w*S.diff(Q,w))
    return dict(parameters=ps,a=a,b=b,c=c,q=[q0,q1,q2],Q=Q,V=V,A=A,K=kk,R=rr)

def direction_shape(V,ps):
    V=S.Matrix(V);r=len(ps)
    degrees={sum(mon) for entry in V for mon,co in S.Poly(entry,*ps).terms() if co!=0}
    if len(degrees)!=1 or next(iter(degrees))<1:raise ValueError('positive homogeneous direction required')
    degree=next(iter(degrees));det=S.factor(V.jacobian(ps).det(method='domain-ge'))
    if det==0:return dict(degree=degree,determinant=det,power_compatible=False)
    constant,factors=S.factor_list(det,*ps)
    possible=all(e%r==0 for f,e in factors)
    out=dict(degree=degree,determinant=det,power_compatible=possible)
    if possible:out.update(delta=constant,h0=S.prod(f**(e//r) for f,e in factors))
    return out
