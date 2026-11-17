"""Original-marker Stein action and exact QQ rank-one lattice classification.

The lattice routine consumes a faithful rank-one presentation. It does not
certify that arbitrary supplied rational data come from a Gaussian phase.
Universal and actual-image statements use STEIN_ACTION_AND_RANK_ONE_LATTICES.md.
"""
import sympy as s
t,xi=s.symbols('t xi')


def heat(poly,X,C,sign=1):
    X=tuple(X);C=s.Matrix(C);value=s.expand(poly);answer=value;k=0
    while value:
        k+=1
        value=s.expand(sign*sum(C[i,j]*s.diff(value,X[i],X[j])
                               for i in range(len(X)) for j in range(len(X)))/(2*k))
        answer+=value
    return s.expand(answer)


def inverse_number(A,X,C):
    X=tuple(X);C=s.Matrix(C)
    transformed=heat(A,X,C)
    without_constant=s.Add(*(c*s.prod(x**e for x,e in zip(X,alpha))/sum(alpha)
                             for alpha,c in s.Poly(transformed,*X).terms() if sum(alpha)))
    return heat(without_constant,X,C,-1)


def stein(P,A,X,C):
    X=tuple(X);C=s.Matrix(C)
    if not X or len(set(X))!=len(X) or C.shape!=(len(X),len(X)) or C!=C.T or C.det()==0:
        raise ValueError('Distinct coordinates and symmetric nonsingular covariance required')
    if any(c.is_Rational is not True for c in C):
        raise ValueError('Fixed rational covariance required')
    for value in [P,A]:
        try:s.Poly(value,*sorted(s.sympify(value).free_symbols|set(X),key=str))
        except s.PolynomialError as exc:raise ValueError('Polynomial phase and marker required') from exc
    inverse=inverse_number(A,X,C)
    return s.expand((s.Matrix([s.diff(P,x) for x in X]).T*C*
                     s.Matrix([s.diff(inverse,x) for x in X]))[0])


def shift(h,f=1):
    h=s.sympify(h);f=s.sympify(f)
    if s.cancel(f.subs(t,0)-1)!=0:
        raise ValueError('Normalized unit period required')
    value=h.subs(t,0)
    if value.has(s.zoo,s.nan,s.oo):raise ValueError('Series must be regular at zero')
    return s.cancel((h-value*f)/t)


def seed_matrix(f):
    f=s.Poly(f,t,domain=s.QQ)
    if f.nth(0)!=1 or f.degree()<1:raise ValueError('Nonconstant QQ polynomial with f(0)=1 required')
    m=f.degree()
    matrix=s.Matrix(m,m,lambda i,j:s.Poly(shift(t**j,f.as_expr()),t).nth(i))
    expected=s.Poly(sum(f.nth(j)*xi**(m-j) for j in range(m+1)),xi)
    assert matrix.charpoly(xi).as_expr()==expected.as_expr()
    return dict(basis=[str(t**j) for j in range(m)],matrix=[[str(x) for x in row] for row in matrix.tolist()],
                characteristic_polynomial=str(expected.as_expr()),dimension=m,
                nonzero_spectrum_is_reciprocal_seed_zeroes=True)


def rational(value):
    value=s.cancel(value)
    a,b=s.fraction(value)
    s.Poly(a,t,domain=s.QQ);s.Poly(b,t,domain=s.QQ)
    if b.subs(t,0)==0:raise ValueError('Rational connection data must be regular at normalization origin')
    return value


def factors(value):
    result={}
    for q,m in s.factor_list(s.denom(s.cancel(value)),t)[1]:
        result[s.Poly(q,t,domain=s.QQ).monic().as_expr()]=m
    return result


def classify_rank_one(log_derivative,normalized_generators):
    """Classify D_r closure of 1 and finitely many faithful rational ratios.

    f'=r*f,f(0)=1. Input normalized generators represent h_i/f, not h_i.
    Output M=(f/p)*QQ[t,d^-1]. No Gaussian realization or rank is inferred.
    """
    r=rational(log_derivative);generators=[s.Integer(1)]+[rational(v) for v in normalized_generators]
    rf=factors(r);allf=set(rf)
    for g in generators:allf.update(factors(g))
    p=d=s.Integer(1);profiles=[]
    num,den=s.fraction(r)
    for q in sorted(allf,key=str):
        order=rf.get(q,0);bound=0;residue=None;forced=False
        if order>=2:forced=True
        elif order==1:
            other=s.cancel(den/q)
            residue=s.rem(num*s.invert(s.diff(q,t)*other,q,t),q,t)
            if residue.is_Integer and residue>0:bound=int(residue)
            else:forced=True
        marker_order=max(factors(g).get(q,0) for g in generators)
        infinite=forced or marker_order>bound
        factor=q/q.subs(t,0)
        if infinite:d*=factor
        elif bound:p*=factor**bound
        profiles.append(dict(factor=str(q),connection_pole_order=order,
            residue=None if residue is None else str(residue),generator_pole_order=marker_order,
            permitted_finite_order=bound,infinite_pole_tower=infinite,
            final_pole_order='infinity' if infinite else bound))
    p=s.expand(p);d=s.expand(d);new_r=s.cancel(r-s.diff(p,t)/p)
    proper=(r==0 or s.degree(s.numer(r),t)<s.degree(s.denom(r),t))
    regular=all(m<=1 for m in rf.values()) and bool(proper)
    return dict(log_derivative=str(r),normalized_generators=[str(g) for g in generators],
        finite_pole_polynomial=str(p),infinite_pole_divisor=str(d),gauge_log_derivative=str(new_r),
        pole_profiles=profiles,normal_form='(f/p)*QQ[t,d^-1]',
        safe=(d==1 and new_r==0),pure=(r==0),globally_regular=regular,
        scope='Exact D_r closure of supplied faithful rank-one rational generators;not a Gaussian realization/rank certificate',
        Gaussian_realization_inferred=False,actual_rank_inferred_from_finite_jets=False)
