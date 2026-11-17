"""Exact finite general-exponent local and resultant certificates.

The analytic rate is proved in the note using published theorems, never
inferred from a finite squarefreeness sample. The exponent is a fixed integer.
"""
from math import comb,factorial
import sympy as s
import prime_power_windows as p


def coefficients(m,n):
    m,n=p.integer(m,1),p.integer(n)
    return [comb(n,j)*factorial(m*j) for j in range(n+1)]


def _poly(coeffs,prime=None):
    t=s.Symbol('t');expression=sum(a*t**j for j,a in enumerate(coeffs))
    return s.Poly(expression,t,modulus=prime) if prime else s.Poly(expression,t,domain=s.ZZ)


def _mod_coeffs(f,prime):
    return [0] if f.is_zero else [int(f.nth(j)) % prime for j in range(f.degree()+1)]


def shifted_profile(m,a,prime,offset):
    m,a,prime,i=p.integer(m,1),p.integer(a,1),p.integer(prime,2),p.integer(offset)
    if not s.isprime(prime) or prime<=m*max(a,i):
        raise ValueError('prime>m*max(a,offset) required')
    n=a*prime+i
    vertices=[[0,0]]+([[i,0]] if i else [])+[[n,m*a]]
    return {'exponent':m,'multiplier':a,'prime':prime,'offset':i,'terminal_index':n,
            'factorial_range':m*n,'factorial_range_below_p_squared':m*n<prime*prime,
            'vertices':vertices,'sloped_reduced_denominator':prime,'sloped_numerator':m,
            'horizontal_square_degree_bound':i//2,'sloped_square_units_bound':a//2,
            'parity_used':False,'theorem':'FIRALLMLOCAL20261116'}


def preceding_prime_bound(m,n):
    m,n=p.integer(m,1),p.integer(n,2)
    q=int(s.prevprime(n+1));gap=n-q
    if q<=m*max(1,gap):raise ValueError('preceding prime must exceed m*max(1,gap)')
    profile=shifted_profile(m,1,q,gap)
    return {'exponent':m,'terminal_index':n,'preceding_prime':q,'gap':gap,
            'gap_below_sqrt_n':gap*gap<n,'square_divisor_degree_bound':gap//2,
            'profile':profile,'parity_used':False,'theorem':'FIRALLMLOCAL20261116'}


def local_squarefree_certificate(m,n,prime,offset):
    m,n,prime,i=p.integer(m,1),p.integer(n,2),p.integer(prime,2),p.integer(offset)
    if n<=i or (n-i) % prime:raise ValueError('n=a*p+offset with positive integer a required')
    profile=shifted_profile(m,(n-i)//prime,prime,i)
    gap=preceding_prime_bound(m,n)
    f=_poly(coefficients(m,i),prime)
    left,right,divisor=s.gcdex(f,f.diff())
    accepted=gap['square_divisor_degree_bound']<prime and divisor.degree()==0
    return {'status':'general_exponent_local_certified' if accepted else 'this_local_witness_unresolved',
            'exponent':m,'terminal_index':n,'prime':prime,'offset':i,'shifted_profile':profile,
            'preceding_prime_bound':gap,'residual_coefficients':_mod_coeffs(f,prime),
            'residual_gcd_coefficients':_mod_coeffs(divisor,prime),'residual_gcd_degree':divisor.degree(),
            'bezout_left':_mod_coeffs(left,prime),'bezout_right':_mod_coeffs(right,prime),
            'squarefree_certified':accepted,'counterexample_claimed':False,'parity_used':False,
            'numerical_density_evidence':False,'all_index_SFC_claimed':False,'gamma_shape':'1',
            'theorem':'FIRALLMLOCAL20261116'}


def resultant_size_certificate(m,offset):
    """One finite offset, with no irreducibility assumption or large-H sweep."""
    m,i=p.integer(m,1),p.integer(offset,1)
    coeffs=coefficients(m,i);f=_poly(coeffs);R=int(s.resultant(f,f.diff()))
    A=2**i*factorial(m*i)
    bound=i**i*A**(2*i-1)
    return {'exponent':m,'offset':i,'coefficients':[str(a) for a in coeffs],
            'resultant':str(R),'resultant_nonzero':R!=0,
            'A':str(A),'coefficient_l2_squared':str(sum(a*a for a in coeffs)),
            'derivative_l2_squared':str(sum(j*j*coeffs[j]**2 for j in range(1,i+1))),
            'sylvester_row_counts':{'polynomial':i-1,'derivative':i},
            'hadamard_absolute_bound':str(bound),'bound_formula':'i^i*(2^i*(m*i)!)^(2i-1)',
            'irreducibility_claimed':False,'all_offsets_squarefree_claimed':False,
            'numerical_density_evidence':False,'theorem':'FIRQUANTDENSITY20261116'}
