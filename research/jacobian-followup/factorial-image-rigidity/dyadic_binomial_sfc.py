"""Exact witnesses for X^(2^r)*(A+B*X^(2^r)) Strong Factorial.

The theorem is DYADIC_BINOMIAL_SFC.tex. Structural certificates avoid
constructing enormous exponents or factorials; literal polynomial helpers
are intended for explicitly supplied, manageable replay inputs.
"""
from math import comb,factorial,prod
import sympy as s
Z=s.Symbol('z')

def bessel_coefficients(n):
    return [comb(n+k,2*k)*prod(range(1,2*k,2)) for k in range(n+1)]

def bessel_polynomial(n,h=0):
    if not isinstance(n,int) or n<0 or h not in (0,1):raise ValueError('n>=0 and h=0,1 required')
    return s.Poly(sum(c*(n+k+1 if h else 1)*Z**k for k,c in enumerate(bessel_coefficients(n))),Z)

def factorial_unit(q,L,delta=0):
    if not isinstance(q,int) or q<1 or q&(q-1):raise ValueError('dyadic q required')
    if L<0 or not 0<=delta<q:raise ValueError('L>=0,0<=delta<q required')
    num=factorial(q*L+delta);den=2**((q-1)*L)*factorial(L)*factorial(delta)
    quotient,remainder=divmod(num,den)
    assert not remainder and quotient%2
    return quotient

def normalized_moment(q,beta,n):
    if not isinstance(beta,int) or not 1<=beta<=2*q or n<0:raise ValueError('1<=beta<=2q,n>=0 required')
    h,delta=divmod(beta-1,q)
    coeff=[]
    for k,c in enumerate(bessel_coefficients(n)):
        coeff.append(factorial_unit(q,n+k+h,delta)*c*(n+k+1 if h else 1))
    return s.Poly(sum(c*Z**k for k,c in enumerate(coeff)),Z)

def structural_certificate(r,beta=1):
    if not isinstance(r,int) or r<0 or not isinstance(beta,int) or beta<1:
        raise ValueError('integer r>=0,beta>=1 required')
    if (beta-1).bit_length()>r+1:raise ValueError('shape exceeds the proved range beta<=2q')
    h=int((beta-1).bit_length()>r)
    return {'accepted':True,'gap':'2^'+str(r),'lower_exponent':'2^'+str(r),
        'Gamma_shape_hex':hex(beta),'base_shift_h':h,'all_starting_indices':True,
        'proof':'DYADIC_BINOMIAL_SFC.tex','rational_unit_sequence_transfer':True,
        'degree_preserving_mod2_reduction':True,'finite_index_extrapolation':False,
        'all_lower_exponents_or_gaps':False,'all_positive_Gamma_shapes':False,
        'unrestricted_SFC_closed':False,'general_FC_closed':False,
        'independent_replay':False,'formal_verification':False,'external_review':False}


def exact_resultant_valuation(q,n):
    if not isinstance(q,int) or q<1 or q&(q-1) or not isinstance(n,int) or n<1:
        raise ValueError('dyadic q and n>=1 required')
    return 3*q*n*(n-1)-(n-1)*n.bit_count()-n*(n-1).bit_count()


def cubic_residual_kernel(d):
    if not isinstance(d,int) or d<1:raise ValueError('positive integer offset required')
    return s.Poly(sum((-1)**j*factorial(d-1+j)//(factorial(j)*factorial(d-1-2*j))*Z**j for j in range((d-1)//2+1)),Z)

def cubic_prime_window(n,q):
    if not isinstance(n,int) or n<6 or not isinstance(q,int) or not s.isprime(q) or not n<q or 2*q>3*(n-1):
        raise ValueError('n>=6,prime n<q<=3(n-1)/2 required')
    d=q-n;L=(d-1)//2
    A=cubic_residual_kernel(d);B=cubic_residual_kernel(d+1)
    residual=s.gcd(s.Poly(A,Z,modulus=q),s.Poly(B,Z,modulus=q))
    even_bound=2*((d-1)//4)
    return {'terminal_index':n,'prime':q,'offset':d,
        'horizontal_degree_bound':L,'real_even_degree_bound':even_bound,
        'residual_gcd_degree':residual.degree(),
        'rational_square_degree_bound':d//4,'real_even_square_degree_bound':2*(d//8),
        'pair_certified_coprime':min(even_bound,residual.degree())<=1 or 2*(d//8)==0,
        'all_indices_covered':False,'counterexample_claimed':False}


def lower_one_tilted(n,gap):
    if not isinstance(n,int) or n<1 or not isinstance(gap,int) or gap<1:
        raise ValueError('positive integers n,gap required')
    return s.Poly(sum(s.binomial(n,k)*s.rf(n,gap*k)*Z**k for k in range(n+1)),Z)

def cubic_square_kernel(d):
    if not isinstance(d,int) or d<1:raise ValueError('positive integer offset required')
    return s.Poly(sum((-1)**j*d*factorial(d+j-1)//(factorial(j)*factorial(d-2*j))*Z**j for j in range(d//2+1)),Z)
