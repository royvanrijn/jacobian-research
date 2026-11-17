"""Exact Padé reductions and symbolic shifted separation certificates.

Universal evidence is PADE_ADJACENT_MONOMIALS.tex. A failed certificate is
unresolved by this sufficient test, not a counterexample to FC or SFC.
"""
from math import factorial, comb, prod
import sympy as s
N=s.Symbol('n',integer=True,positive=True)
Z=s.Symbol('z')
U=s.Symbol('u',nonnegative=True)


def integer(value,minimum=0):
    if isinstance(value,bool) or not isinstance(value,int) or value<minimum:
        raise ValueError('integer outside the specified range')
    return value


def pade_pair(M,Ndegree,variable=Z):
    M,Ndegree=integer(M),integer(Ndegree)
    S=M+Ndegree
    P=sum(s.Rational(factorial(S-j)*factorial(M),factorial(S)*factorial(M-j)*factorial(j))*variable**j for j in range(M+1))
    Q=sum(s.Rational((-1)**j*factorial(S-j)*factorial(Ndegree),factorial(S)*factorial(Ndegree-j)*factorial(j))*variable**j for j in range(Ndegree+1))
    return s.Poly(P,variable),s.Poly(Q,variable)


def auxiliary_symbolic(a):
    a=integer(a,1)
    coefficients=[]
    for k in range(a):
        value=sum(s.rf(N+1,i)*s.ff(a*N,i)/s.factorial(i)*(-1)**(k-i)*s.rf(a*(N-1)+1,k-i)*s.ff(N-1,k-i)/s.factorial(k-i) for i in range(k+1))
        if k==a-1:value+=s.ff(a*N,a)/N
        coefficients.append(s.factor(s.expand(value)))
    return s.Poly(sum(c*Z**(a-1-k) for k,c in enumerate(coefficients)),Z)


def schwarz_parameters(poly):
    degree=poly.degree()
    if degree<1 or poly.LC()!=1:raise ValueError('positive-degree monic polynomial required')
    b0=poly.nth(degree-1)
    if b0==0:return None
    E=s.Poly(sum(poly.nth(j)*Z**j for j in range(degree+1) if j%2==degree%2),Z)
    O=s.Poly(sum(poly.nth(j)*Z**j for j in range(degree+1) if j%2!=degree%2)/b0,Z)
    params=[s.factor(b0)]
    while O.degree()>0:
        remainder=E-s.Poly(Z,Z)*O
        if remainder.is_zero:return None
        b=s.factor(remainder.LC())
        if b==0:return None
        params.append(b)
        E,O=O,s.Poly(s.cancel(remainder.as_expr()/b),Z)
    if O!=s.Poly(1,Z):return None
    return params


def positive_ratio_data(expr):
    numerator,denominator=s.cancel(expr.subs(N,U+2)).as_numer_denom()
    p,q=s.Poly(numerator,U),s.Poly(denominator,U)
    if q.nth(0)<0:p,q=-p,-q
    accepted=all(c>=0 for c in p.all_coeffs()+q.all_coeffs()) and p.nth(0)>0 and q.nth(0)>0
    return {'accepted':bool(accepted),'numerator':[str(p.nth(j)) for j in range(p.degree()+1)],
            'denominator':[str(q.nth(j)) for j in range(q.degree()+1)],
            'variable':'u=n-2','coefficient_order':'increasing degree'}


def reconstruct(parameters):
    d0,d1=s.Poly(1,Z),s.Poly(Z+parameters[0],Z)
    for b in parameters[1:]:d0,d1=d1,s.Poly(s.cancel(Z*d1.as_expr()+b*d0.as_expr()),Z)
    return d1


def shifted_all_index_certificate(a):
    a=integer(a,2)
    H=auxiliary_symbolic(a);c=(a-1)*N+2
    shifted=s.Poly(s.expand(H.as_expr().subs(Z,Z+c)),Z)
    parameters=schwarz_parameters(shifted)
    if parameters is None:
        return {'a':a,'status':'unresolved_by_shift_certificate','counterexample_claimed':False}
    identity=s.simplify(reconstruct(parameters).as_expr()-shifted.as_expr())==0
    ratios=[positive_ratio_data(b) for b in parameters]
    accepted=bool(identity and all(row['accepted'] for row in ratios))
    return {'a':a,'status':'all_indices_certified' if accepted else 'unresolved_by_shift_certificate',
            'shift':str(c),'H_coefficients':[str(v) for v in H.all_coeffs()],
            'b_parameters':[str(s.factor(v)) for v in parameters],
            'positive_ratio_certificates':ratios,'reconstruction_identity':bool(identity),
            'all_n_ge2_certified':accepted,'all_a_claimed':False,'counterexample_claimed':False,
            'scope':'ordinary X^a(A+B*X); a supplied integer, arbitrary complex coefficients'}


def indexed_strip_certificate(a,n):
    a,n=integer(a,1),integer(n,2)
    accepted=a>=2*n-3
    return {'a':a,'terminal_index':n,'moment_indices':[n-1,n],
            'current_denominator_real_part_lower':(a-1)*n+2,
            'lower_bound_strict':True,
            'preceding_denominator_real_part_upper':(a+1)*(n-1),
            'accepted':accepted,'counterexample_claimed':False,
            'all_indices_for_this_a_claimed':False,'accepted_is_only_the_strip_test':True,'other_certificates_may_apply':True}


def auxiliary_at_integer(a,n):
    """Literal integer coefficient formula without expanding symbolic n."""
    a,n=integer(a,1),integer(n,2)
    coefficients=[]
    for k in range(a):
        value=0
        for j in range(min(k,n-1)+1):
            i=k-j
            value+=prod(range(n+1,n+1+i))*comb(a*n,i)*(-1)**j*prod(range(a*(n-1)+1,a*(n-1)+1+j))*comb(n-1,j)
        if k==a-1:
            top=prod(range(a*n-a+1,a*n+1));assert top%n==0
            value+=top//n
        coefficients.append(value)
    return s.Poly(sum(c*Z**(a-1-k) for k,c in enumerate(coefficients)),Z)


def finite_shift_certificate(a,n):
    a,n=integer(a,1),integer(n,2)
    H=auxiliary_at_integer(a,n)
    if H.degree()==0:
        return {'a':a,'terminal_index':n,'accepted':H!=s.Poly(0,Z),
                'status':'nonzero_auxiliary_constant','all_indices_claimed':False,'counterexample_claimed':False}
    shifted=s.Poly(H.as_expr().subs(Z,Z+(a-1)*n+2),Z)
    parameters=schwarz_parameters(shifted)
    accepted=bool(parameters is not None and all(b>0 for b in parameters)
                  and s.simplify(reconstruct(parameters).as_expr()-shifted.as_expr())==0)
    return {'a':a,'terminal_index':n,'moment_indices':[n-1,n],
            'accepted':accepted,'status':'terminal_window_certified' if accepted else 'unresolved_by_shift_criterion',
            'H_coefficients':[str(c) for c in H.all_coeffs()],
            'b_parameters_hex':[{'numerator':hex(int(s.numer(c))),'denominator':hex(int(s.denom(c)))} for c in parameters] if parameters is not None else None,
            'parameter_encoding':'exact hexadecimal integer numerator/denominator',
            'all_indices_claimed':False,'all_a_claimed':False,'counterexample_claimed':False}


def pencil_symbolic(a):
    a=integer(a,2);dimension=a-1;lower=a*(N-1)
    A=s.zeros(dimension);B=s.zeros(dimension)
    for j in range(dimension):
        A[j,j]=lower+j+1+N
        B[j,j]=1 if j<dimension-1 else s.Rational(1,a+1)
        if j<dimension-1:A[j,j+1]=-(lower+j+2)
        if j:B[j,j-1]=-1
    return A,B
