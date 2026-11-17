"""Finite witnesses for the density-one proof; the asymptotic uses analytic imports."""
from math import comb,factorial
import sympy as s
import prime_power_windows as p
import prime_gap_windows as g


def _coefficients(offset):
    return [comb(offset,j)*factorial(6*j) for j in range(offset+1)]


def _poly(coefficients,prime=None):
    t=s.Symbol('t')
    expression=sum(a*t**j for j,a in enumerate(coefficients))
    return s.Poly(expression,t,modulus=prime) if prime else s.Poly(expression,t,domain=s.ZZ)


def _mod_coefficients(poly,prime):
    if poly.is_zero:return [0]
    return [int(poly.nth(j)) % prime for j in range(poly.degree()+1)]


def fixed_offset_certificate(offset):
    """One supplied finite offset; no implicit growing-H computation."""
    offset=p.integer(offset)
    f=_poly(_coefficients(offset))
    squarefree=s.gcd(f,f.diff()).degree()==0
    disc=int(s.discriminant(f)) if offset else 1
    assert squarefree==(disc!=0)
    return {'offset':offset,'characteristic_zero_squarefree':squarefree,
            'integer_discriminant':str(disc),
            'safe_all_larger_primes_threshold':max(6*offset,abs(disc)) if squarefree else None,
            'irreducibility_claimed':False,'all_offsets_claimed':False,
            'numerical_density_evidence':False}


def shifted_good_offset_certificate(n,prime,offset):
    """One exact local witness, including offsets beyond the uniform prefix."""
    n,prime,offset=p.integer(n,30),p.integer(prime,2),p.integer(offset)
    if not s.isprime(prime) or n<=offset or (n-offset) % prime:
        raise ValueError('prime p and positive integer a with n=a*p+offset required')
    a=(n-offset)//prime
    if prime<=6*max(a,offset):raise ValueError('prime>6max(a,offset) required')
    coefficients=[a % prime for a in _coefficients(offset)]
    f=_poly(coefficients,prime)
    left,right,divisor=s.gcdex(f,f.diff())
    gap=g.prime_gap_certificate(n)
    degree_bound=gap['square_divisor_degree_bound']
    residual_squarefree=divisor.degree()==0
    accepted=degree_bound<prime and residual_squarefree
    return {'status':'shifted_good_offset_certified' if accepted else 'this_local_witness_unresolved',
            'terminal_index':n,'multiplier':a,'prime':prime,'offset':offset,
            'preceding_prime':gap['preceding_prime'],'gap':gap['gap'],
            'square_divisor_degree_bound':degree_bound,'degree_bound_below_shifted_prime':degree_bound<prime,
            'residual_coefficients':coefficients,'residual_degree':f.degree(),
            'residual_gcd_coefficients':_mod_coefficients(divisor,prime),
            'residual_gcd_degree':divisor.degree(),'residual_squarefree':residual_squarefree,
            'bezout_left':_mod_coefficients(left,prime),'bezout_right':_mod_coefficients(right,prime),
            'squarefree_certified':accepted,'counterexample_claimed':False,
            'degree_preservation_reason':'g(0)=+-1 and all factor slopes horizontal when degree<p',
            'all_index_SFC_claimed':False,'numerical_density_evidence':False,
            'gamma_shape':'1','theorem':'FIRDENSITYONE20261116'}
