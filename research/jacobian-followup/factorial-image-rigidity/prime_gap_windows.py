"""Ordinary factorial prime-gap and retained small-kernel certificates.

All-index SFC remains open. The retained kernel prefix is0..12;
eventual fixed-multiplier claims do not provide an effective cutoff.
"""
from functools import lru_cache
from math import comb,factorial,gcd,prod
import sympy as s
import prime_power_windows as p

MAX_KERNEL_OFFSET=12


def _prime(q):
    q=p.integer(q,2)
    if not s.isprime(q):raise ValueError('prime required')
    return q


def exposed_edge(n,q):
    n,q=p.integer(n,1),_prime(q)
    if not (q>6 and 2*q>n and q<=n and q*q>6*n):
        raise ValueError('q>6, n/2<q<=n, and q^2>6n required')
    d=n-q
    h=max(range(d+1),key=lambda j:6*j % q)
    residue=6*h % q
    return {'terminal_index':n,'prime':q,'gap':d,'support_index':h,
            'edge_left':[h,6*h//q],'edge_right':[h+q,6+6*h//q],
            'width':q,'height':6,'reduced_denominator':q,
            'minimum_adjusted_height_numerator':-residue,
            'adjusted_height_denominator':q,'remaining_polygon_width':d,
            'square_divisor_degree_bound':d//2,'even_square_divisor_degree_bound':2*(d//4),
            'full_polygon_claimed':False,'gamma_shape':'1','theorem':'FIRPRIMEGAP20261116'}


def prime_gap_certificate(n):
    n=p.integer(n,30)
    q=int(s.prevprime(n+1))
    cert=exposed_edge(n,q)
    cert.update(preceding_prime=q,candidate_even_degrees=list(range(2,2*((n-q)//4)+1,2)),
                all_n_ge_30_theorem=True,all_index_SFC_claimed=False)
    return cert


def fourteenth_certificate(n):
    n=p.integer(n,30)
    q=int(s.nextprime(5*n//6))
    assert 6*q>5*n and q<n, 'Nagura interval witness failed'
    d=n-q
    if 7*q>6*n:
        vertices=[[0,0],[d,0],[n,6]]
        bound=d//2
        branch='q>6n/7: exposed width-q edge'
    else:
        L=q//6
        vertices=[[0,0],[L,0],[L+q,6],[n,7]]
        bound=L//2
        branch='q<=6n/7: two primitive nonhorizontal edges'
    edges=[]
    for left,right in zip(vertices,vertices[1:]):
        width,height=right[0]-left[0],right[1]-left[1]
        common=gcd(width,height)
        edges.append({'left':left,'right':right,'e':width//common,'h':height//common})
    assert 14*bound<n
    return {'terminal_index':n,'prime':q,'gap':d,'branch':branch,'expected_slope_edges':edges,
            'square_divisor_degree_bound':bound,'even_square_divisor_degree_bound':2*(bound//2),
            'uniform_bound':'degree<n/14 for every n>=30','published_prime_input':'Nagura1952 at5n/6>=25',
            'gamma_shape':'1','theorem':'FIRFOURTEENTH20261116'}


def _egcd(a,b):
    x,y,u,v=1,0,0,1
    aa,bb=abs(a),abs(b)
    while bb:
        q,rem=divmod(aa,bb)
        aa,bb=bb,rem
        x,u=u,x-q*u
        y,v=v,y-q*v
    return aa,x*(1 if a>=0 else -1),y*(1 if b>=0 else -1)


@lru_cache(maxsize=13)
def _kernel_values(d):
    t=s.Symbol('t')
    coefficients=tuple(comb(d,j)*factorial(6*j) for j in range(d+1))
    f=s.Poly(sum(a*t**j for j,a in enumerate(coefficients)),t,domain=s.ZZ)
    seq=s.subresultants(f,f.diff())
    # The retained prefix is regular, so its degree-one member is Sres_1.
    assert [int(a.degree()) for a in seq]==list(range(d,-1,-1)), 'nonregular PRS needs determinant normalization'
    S1=seq[-2]
    return coefficients,int(seq[-1].nth(0)),int(S1.nth(1)),int(S1.nth(0))


def small_kernel_certificate(d):
    """Exact retained integer certificate, universal over primes q>6d."""
    d=p.integer(d)
    if d>MAX_KERNEL_OFFSET:raise ValueError('only the retained kernel prefix0..12 is certified')
    if d<2:
        return {'offset':d,'status':'kernel_square_degree_at_most_one_certified',
                'reason':'constant or linear polynomial has no nonconstant square divisor',
                'kernel_square_degree_bound':0,'all_offsets_claimed':False,'theorem':'FIRSMALLKERNEL20261116'}
    coefficients,R,A,B=_kernel_values(d)
    first,x,y=_egcd(R,A)
    C,u,v=_egcd(first,B)
    bezout=(u*x,u*y,v)
    assert sum(a*b for a,b in zip((R,A,B),bezout))==C
    remainder=C
    factors={}
    if C:
        for q in s.primerange(2,6*d+1):
            power=0
            while remainder % q==0:
                remainder//=int(q);power+=1
            if power:factors[str(q)]=power
    accepted=C>0 and remainder==1
    return {'offset':d,'status':'kernel_square_degree_at_most_one_certified' if accepted else 'kernel_integer_test_unresolved',
            'coefficients':[str(a) for a in coefficients],'resultant':str(R),
            'first_subresultant_linear':str(A),'first_subresultant_constant':str(B),
            'integer_gcd':str(C),'integer_gcd_bezout':[str(x) for x in bezout],
            'small_prime_factorization':factors,'remainder_after_small_primes':str(remainder),
            'prime_condition':f'q>{6*d}','modular_gcd_degree_bound':1 if accepted else None,
            'kernel_square_degree_bound':1 if accepted else None,
            'all_offsets_claimed':False,'counterexample_claimed':False,'theorem':'FIRSMALLKERNEL20261116'}


def near_prime_window(n,q):
    n,q=p.integer(n,1),_prime(q)
    d=n-q
    if d<0 or d>12 or q<=6*max(1,d):
        raise ValueError('n=q+d, offset0..12, and q>6*max(1,d) required')
    kernel=small_kernel_certificate(d)
    assert kernel['status']=='kernel_square_degree_at_most_one_certified'
    return {'terminal_index':n,'prime':q,'offset':d,'horizontal_residual':'F_offset modq',
            'kernel_certificate':kernel,'squarefree_certified':True,
            'moment_indices':[n-1,n],'gamma_shape':'1','theorem':'FIRSMALLKERNEL20261116'}


def shifted_family_certificate(a,prime,offset):
    """Fixed kernels and n/14 rules; qualitative eventual scope is separate."""
    a,prime,i=p.integer(a,1),_prime(prime),p.integer(offset)
    if i>12 or prime<=6*max(a,i):
        raise ValueError('offset0..12 and prime>6*max(a,offset) required')
    kernel=small_kernel_certificate(i)
    assert kernel['status']=='kernel_square_degree_at_most_one_certified'
    raw=[(h,t) for h in range(i//2+1) for t in range(a//2+1)
         if (h or t) and (h+t) % 2==0]
    kernel_removed=[(h,t) for h,t in raw if h>kernel['kernel_square_degree_bound']]
    full_removed=[(h,t) for h,t in raw if 2*h==i and 2*t==a]
    band_removed=[]
    if a*prime+i>=30:
        for h,t in raw:
            coefficient=14*t-a
            constant=14*h-i
            if coefficient>=0 and coefficient*prime+constant>=0:
                band_removed.append((h,t))
    remaining=[pair for pair in raw if pair not in kernel_removed and pair not in full_removed and pair not in band_removed]
    accepted=not remaining
    encode=lambda pairs:[{'horizontal_degree':h,'sloped_residual_degree':t} for h,t in pairs]
    return {'status':'all_r_shifted_family_certified' if accepted else 'all_r_test_unresolved',
            'multiplier':a,'prime':prime,'offset':i,'gamma_shape':'1',
            'terminal_indices':f'{a}*{prime}^r+{i}, EVERY integer r>=1',
            'raw_even_patterns':encode(raw),'small_kernel_removed_patterns':encode(kernel_removed),
            'full_square_removed_patterns':encode(full_removed),'fourteenth_removed_patterns':encode(band_removed),
            'remaining_patterns':encode(remaining),'all_powers_r_ge1_certified':accepted,
            'eventual_fixed_multiplier_squarefree_theorem':True,
            'eventual_statement_uses':'preceding-prime gap o(n) from PNT, plus retained horizontal kernel',
            'effective_cutoff':None,'all_indices_SFC_claimed':False,'counterexample_claimed':False,
            'theorem':'FIRGAPKERNELFAMILY20261116'}
