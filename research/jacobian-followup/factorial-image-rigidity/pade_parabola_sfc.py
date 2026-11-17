"""Exact instantiations of the uniform parabola proof for adjacent supports.

The universal proof is UNIFORM_PADE_PARABOLA_SFC.tex. No spectral or
finite-index extrapolation is used by the structural certificate.
"""
import sympy as s
Z=s.Symbol('z')

def R(t,degree,variable=Z):
    if not isinstance(degree,int) or degree<0:raise ValueError('nonnegative integer degree required')
    t=s.sympify(t)
    return s.Poly(sum((-1)**k*s.binomial(degree,k)*s.rf(t+1,k)*variable**(degree-k) for k in range(degree+1)),variable)

def pencil(l,a,N):
    if not isinstance(a,int) or a<2:raise ValueError('integer a>=2 required')
    l,N=s.sympify(l),s.sympify(N);M=l+a;D=a-1
    A=s.diag(*[l+N+j for j in range(1,D+1)]);B=s.eye(D)
    for j in range(D-1):A[j,j+1]=-(l+j+2);B[j+1,j]=-1
    B[-1,-1]=N/(M+N)
    return A,B

def normalized_continuants(l,a,N,variable=Z):
    l,N=s.sympify(l),s.sympify(N);M=l+a;D=a-1
    if D<1:raise ValueError('a>=2 required')
    result=[s.Poly(1,variable)]
    for j in range(1,D+1):
        L=l+N+j;beta=s.Integer(1) if j<D else N/(M+N)
        expr=(1+beta*variable/L)*result[-1].as_expr()
        if j>=2:expr-=variable*(l+j)/(L*(L-1))*result[-2].as_expr()
        result.append(s.Poly(s.cancel(expr),variable))
    return result

def adjacent_monomial_certificate(a,n,beta=1):
    if not isinstance(a,int) or a<0 or not isinstance(n,int) or n<1:
        raise ValueError('integer a>=0,n>=1 required')
    beta=s.Rational(beta)
    if beta<=0:raise ValueError('positive Gamma shape required')
    result={'accepted':True,'a':a,'terminal_index':n,'Gamma_shape':str(beta),
        'proof':'UNIFORM_PADE_PARABOLA_SFC.tex','criterion':'shared parabola; no uniform vertical bound assumed',
        'unrestricted_SFC_closed':False,'general_FC_closed':False,
        'independent_replay':False,'formal_verification':False,'external_review':False}
    if n==1:result['case']='order-zero moment is one'
    elif a==0:result['case']='backward nonzero three-term recurrence'
    elif a==1:result['case']='two consecutive parameter zeros force an identically zero polynomial'
    else:
        l=beta+a*(n-1)-1;M=l+a;D=a-1
        first=l+n+1 if D>=2 else (M+n)*(M+n-1)/n
        result.update({'case':'complementary parabola inclusions','l':str(l),'M':str(M),
            'first_b_minus_n':str(first-n),'interior_parabola_parameter':n,
            'last_parabola_parameter':str(M+n) if D>=2 else str(first),
            'auxiliary_degree':D,'root_tau_range':[0,n-1]})
        assert first>n and M+n>n and l>-1
    return result


def gap_R(t,gap,degree,variable=Z):
    if not isinstance(gap,int) or gap<1 or not isinstance(degree,int) or degree<0:
        raise ValueError('positive integer gap, nonnegative integer degree required')
    t=s.sympify(t)
    return s.Poly(sum((-1)**k*s.binomial(degree,k)*s.rf(t+1,gap*k)*variable**(degree-k) for k in range(degree+1)),variable)
