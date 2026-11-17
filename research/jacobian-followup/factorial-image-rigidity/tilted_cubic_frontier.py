"""Tilted cubic Newton certificates and uniform prime-power families.

Universal statements are proved in TILTED_CUBIC_PRIME_POWER_AND_FOLD.tex.
Literal polynomial functions require explicitly supplied manageable n;
structural family certificates do not construct huge prime powers.
"""
from math import comb,factorial,gcd
import sympy as s
from prime_power_windows import factorial_valuation
from newton_slope_coverage import slope_square_divisor_degrees
Z=s.Symbol('t')
FAMILIES={(1,0),(1,1),(1,2),(1,3),(2,0),(2,1),(2,2),(3,0),(3,1),(4,0)}

def coefficients(n):
    if not isinstance(n,int) or n<1:raise ValueError('positive integer n required')
    return [comb(n,j)*factorial(n+2*j-1)//factorial(n-1) for j in range(n+1)]

def tilted(n):return s.Poly(sum(c*Z**j for j,c in enumerate(coefficients(n))),Z)

def valuation_edges(n,p):
    if not isinstance(n,int) or n<1 or not s.isprime(p):raise ValueError('positive n and prime p required')
    hull=[]
    for j in range(n+1):
        v=factorial_valuation(n,p)-factorial_valuation(j,p)-factorial_valuation(n-j,p)+factorial_valuation(n+2*j-1,p)-factorial_valuation(n-1,p)
        while len(hull)>1:
            x,y=hull[-2];xx,yy=hull[-1]
            if (xx-x)*(v-yy)-(yy-y)*(j-xx)>0:break
            hull.pop()
        hull.append((j,v))
    edges=[]
    for left,right in zip(hull,hull[1:]):
        width,height=right[0]-left[0],right[1]-left[1];g=gcd(width,height)
        edges.append({'left':list(left),'right':list(right),'e':width//g,'h':height//g})
    return edges

def expected_shifted_edges(a,p,r,i):
    if a<1 or r<1 or i<0 or not s.isprime(p) or p<=3*max(a,i):raise ValueError('shifted theorem hypotheses required')
    P=p**r;S=(P-1)//(p-1);n=a*P+i
    edges=[]
    if i:edges.append({'left':[0,0],'right':[i,0],'e':1,'h':0,'residual_coefficients':[int(s.binomial(i,j)*s.rf(i,2*j))%p for j in range(i+1)]})
    shape=a if i==0 else a+1
    edges.append({'left':[i,0],'right':[n,2*a*S],'e':P,'h':2*S,
        'residual_coefficients':[int(s.binomial(a,j)*s.rf(shape,2*j))%p for j in range(a+1)]})
    return edges

def prime_power_certificate(p,r):
    if not isinstance(r,int) or r<1 or not s.isprime(p):raise ValueError('prime p and r>=1 required')
    slope='2-1/P' if p==2 else ('1-1/P' if p==3 else '2*(P-1)/(P*(p-1))')
    return {'terminal_index':str(p)+'^'+str(r),'P':'p^r','prime':int(p),'power':r,
        'ordinary_Gamma_shape':1,'single_primitive_edge':True,'slope':slope,'irreducible_over_Qp':True,
        'cubic_pair_at_this_terminal_certified':True,'all_starting_indices_claimed':False,
        'proof':'TILTED_CUBIC_PRIME_POWER_AND_FOLD.tex',
        'unrestricted_SFC_closed':False,'general_FC_closed':False,
        'independent_replay':False,'formal_verification':False,'external_review':False}

def uniform_shift_family(a,p,i):
    if not s.isprime(p) or a<1 or i<0 or p<=3*max(a,i):raise ValueError('p>3max(a,i) and integer a>=1,i>=0 required')
    accepted=(a,i) in FAMILIES
    return {'multiplier':a,'prime':int(p),'offset':i,'terminal_indices':str(a)+'*'+str(p)+'^r+'+str(i)+',all integer r>=1',
        'ordinary_Gamma_shape':1,'all_terminals_in_family_certified':accepted,'generic_polygon_proved':True,
        'all_lower_exponents_or_all_starting_indices_claimed':False,
        'unresolved_does_not_mean_counterexample':True}

def slope_certificate(n,primes):
    candidates=set(range(2,n//2+1,2));profiles=[];witness={}
    for p in primes:
        p=int(p);edges=valuation_edges(n,p);allowed=slope_square_divisor_degrees(edges)
        for k in candidates-set(allowed):witness[str(k)]=p
        candidates.intersection_update(allowed)
        profiles.append({'prime':p,'edges':edges,'square_factor_degrees':allowed})
    return {'terminal_index':n,'finite_slope_prime_cap':3*n-1,'primes':list(map(int,primes)),
        'profiles':profiles,'surviving_even_degrees':sorted(candidates),'excluding_prime_by_degree':witness,
        'squarefree_at_this_index_certified':not candidates,'all_indices_claimed':False,
        'residuals_used':False,'counterexample_claimed':False}
