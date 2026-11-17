"""Slope-only diagnostics and proved shifted-index families, factorial shape one.

Finite coverage never implies all-index coverage. Residual tests retain
information even above the slope-only prime cap d*n.
"""
from math import comb, factorial, gcd, prod
import sympy as s
import prime_power_windows as p
import mixed_prime_windows as m


def _prime(value):
    value=p.integer(value,2)
    if not s.isprime(value):
        raise ValueError('prime required')
    return value


def slope_edges(n, degree, prime):
    """Literal lower hull from exact Legendre valuations; no residual tests."""
    n,degree,prime=p.integer(n,1),p.integer(degree,1),_prime(prime)
    vn=p.factorial_valuation(n,prime)
    hull=[]
    for j in range(n+1):
        v=vn-p.factorial_valuation(j,prime)-p.factorial_valuation(n-j,prime)+p.factorial_valuation(degree*j,prime)
        point=(j,v)
        while len(hull)>=2:
            a,b=hull[-2:]
            if (b[0]-a[0])*(v-b[1])-(b[1]-a[1])*(j-b[0])>0:
                break
            hull.pop()
        hull.append(point)
    edges=[]
    for left,right in zip(hull,hull[1:]):
        width,height=right[0]-left[0],right[1]-left[1]
        common=gcd(width,height)
        edges.append({'left':list(left),'right':list(right),'e':width//common,'h':height//common})
    return edges


def slope_square_divisor_degrees(edges):
    """Doubled edge occupancy; this generic budget uses no real-place filter."""
    possible={0}
    for edge in edges:
        width=edge['right'][0]-edge['left'][0]
        e=edge['e']
        assert width>0 and e>0 and width % e==0
        possible={value+e*k for value in possible for k in range(width//(2*e)+1)}
    return sorted(possible)


def slope_coverage_certificate(n,degree=6):
    """Exhaust the entire useful slope prime range for ONE supplied index."""
    n,degree=p.integer(n,2),p.integer(degree,1)
    candidates=set(range(2,n//2+1,2))
    original=sorted(candidates)
    profiles=[]
    exclusions={}
    for prime in s.primerange(2,degree*n+1):
        prime=int(prime)
        edges=slope_edges(n,degree,prime)
        allowed=slope_square_divisor_degrees(edges)
        removed=candidates-set(allowed)
        for k in removed:
            exclusions[str(k)]=prime
        candidates.intersection_update(allowed)
        profiles.append({'prime':prime,'edges':edges,'possible_square_divisor_degrees':allowed})
    return {'status':'slope_parity_window_certified' if not candidates else 'full_slope_method_has_survivors',
            'terminal_index':n,'moment_indices':[n-1,n],'exponent':degree,'gamma_shape':'1',
            'candidate_even_degrees':original,'surviving_even_degrees':sorted(candidates),
            'prime_cap':degree*n,'all_primes_through_cap_checked':True,
            'excluding_prime_by_degree':exclusions,'profiles':profiles,
            'slopes_and_lengths_only':True,'residual_multiplicities_used':False,
            'no_larger_prime_improves_slope_method':True,
            'slope_method_limitation_at_this_index':bool(candidates),
            'counterexample_claimed':False,'all_starting_indices_claimed':False,
            'theorem':'FIRSLOPECAP20261116'}


def candidate_prime_set(n,k,degree=6):
    """A necessary finite shortlist of primes that could exclude degree k."""
    n,k,degree=p.integer(n,2),p.integer(k,1),p.integer(degree,1)
    if 2*k>n:
        raise ValueError('a square-divisor degree must satisfy 2*k<=n')
    small=list(map(int,s.primerange(2,2*degree*k+1)))
    large={}
    for offset in range(2*k):
        for prime in s.factorint(n-offset):
            prime=int(prime)
            if prime>2*degree*k:
                large.setdefault(str(prime),[]).append(offset)
    return {'terminal_index':n,'exponent':degree,'candidate_degree':k,
            'small_prime_bound':2*degree*k,'small_primes':small,
            'large_prime_divisor_offsets':large,
            'candidate_primes':sorted(set(small)|set(map(int,large))),
            'uniform_exclusion_claimed':False,'residual_tests_have_same_prime_restriction':False,
            'theorem':'FIRSLOPECAP20261116'}


def large_degree_certificate(n):
    """Uniform n>=9 exclusion: Nagura supplies 3n<q<18n/5."""
    n=p.integer(n,9)
    prime=int(s.nextprime(3*n))
    assert 5*prime<18*n, 'Nagura interval witness failed'
    length=(prime-1)//6
    edges=[{'left':[0,0],'right':[length,0],'e':1,'h':0},
           {'left':[length,0],'right':[n,1],'e':n-length,'h':1}]
    return {'terminal_index':n,'exponent':6,'gamma_shape':'1','prime':prime,
            'prime_interval':'3n<q<18n/5','expected_slope_edges':edges,
            'square_divisor_degree_bound':length//2,'uniform_strict_bound':'degree<3n/10',
            'published_prime_input':'Nagura1952, Theorem pp180-181, x>=25',
            'shape_extension_claimed':False,'theorem':'FIRLARGESQUARE20261116'}


def shifted_index_profile(degree,multiplier,prime,power,offset):
    """Exact two-kernel formula at N=multiplier*prime**power+offset."""
    d,u,prime,power,i=p.integer(degree,1),p.integer(multiplier,1),_prime(prime),p.integer(power,1),p.integer(offset)
    if prime<=d*max(u,i):
        raise ValueError('prime>d*max(multiplier,offset) required')
    e=prime**power
    c=(e-1)//(prime-1)
    n=u*e+i
    horizontal=[comb(i,j)*factorial(d*j) % prime for j in range(i+1)]
    kernel=[comb(u,j)*factorial(d*j) % prime for j in range(u+1)]
    epsilon=(-1)**(d*c) % prime
    edges=[]
    if i:
        edges.append({'left':[0,0],'right':[i,0],'e':1,'h':0,'residual_coefficients':horizontal})
    edges.append({'left':[i,0],'right':[n,d*u*c],'e':e,'h':d*c,
                  'residual_coefficients':[a*pow(epsilon,j,prime) % prime for j,a in enumerate(kernel)]})
    local=m.square_divisor_degrees(edges,prime)
    even=[k for k in local['possible_degrees'] if k % 2==0]
    full_square=[k for k in even if k>0 and 2*k==n]
    large_removed=[k for k in even if d==6 and n>=9 and 10*k>=3*n]
    surviving=[k for k in even if k not in full_square and k not in large_removed]
    large_witness=large_degree_certificate(n) if large_removed else None
    return {'terminal_index':n,'exponent':d,'multiplier':u,'prime':prime,'power':power,'offset':i,
            'gamma_shape':'1','edges':edges,'horizontal_kernel':horizontal,'nonhorizontal_kernel':kernel,
            'slope_only_degrees':slope_square_divisor_degrees(edges),
            'residual_square_divisor_degrees':local['possible_degrees'],
            'real_parity_filtered_degrees':even,
            'full_square_removed_degrees':full_square,
            'large_degree_removed_degrees':large_removed,'large_degree_witness':large_witness,
            'possible_square_divisor_degrees':surviving,
            'adjacent_window_certified':surviving==[0],
            'all_starting_indices_claimed':False,'theorem':'FIRSHIFTEDINDEX20261116'}


def shifted_family_certificate(degree,multiplier,prime,offset):
    """Certify every r>=1 using fixed kernels and uniform global arithmetic."""
    data=shifted_index_profile(degree,multiplier,prime,1,offset)
    d,u,prime,i=data['exponent'],data['multiplier'],data['prime'],data['offset']
    def kernel_degrees(coefficients):
        degree=len(coefficients)-1
        if not degree:
            return [0]
        edge={'left':[0,0],'right':[degree,0],'e':1,'h':0,'residual_coefficients':coefficients}
        return m.square_divisor_degrees([edge],prime)['possible_degrees']
    hs=kernel_degrees(data['horizontal_kernel'])
    ts=kernel_degrees(data['nonhorizontal_kernel'])
    even_patterns=[{'horizontal_degree':h,'nonhorizontal_residual_degree':t}
                   for h in hs for t in ts if (h or t) and (h+(prime % 2)*t) % 2==0]
    full_square=[pattern for pattern in even_patterns
                 if 2*pattern['horizontal_degree']==i and 2*pattern['nonhorizontal_residual_degree']==u]
    large_removed=[]
    if d==6 and u*prime+i>=9:
        for pattern in even_patterns:
            h,t=pattern['horizontal_degree'],pattern['nonhorizontal_residual_degree']
            coefficient=10*t-3*u
            constant=10*h-3*i
            if coefficient>=0 and coefficient*prime+constant>=0:
                large_removed.append(pattern)
    surviving=[pattern for pattern in even_patterns if pattern not in full_square and pattern not in large_removed]
    accepted=not surviving
    return {'status':'unbounded_shifted_family_certified' if accepted else 'shifted_kernel_test_unresolved',
            'exponent':d,'multiplier':u,'prime':prime,'offset':i,'gamma_shape':'1',
            'terminal_indices':f'{u}*{prime}^r+{i}, EVERY integer r>=1',
            'moment_indices':['N-1','N'],'horizontal_square_degrees':hs,
            'nonhorizontal_residual_square_degrees':ts,'raw_even_patterns':even_patterns,
            'full_square_removed_patterns':full_square,'large_degree_removed_patterns':large_removed,
            'surviving_even_patterns':surviving,
            'global_arithmetic_filters_used':bool(full_square or large_removed),
            'large_degree_filter_valid_for_every_r':True,
            'full_square_argument':'constant1 and ordinary factorial leading coefficient forbid a rational full square',
            'degree_pattern':'h+prime^r*t; kernel budgets and parity are independent of r>=1',
            'all_powers_r_ge_1_certified':accepted,'bounded_search_used_for_universal_claim':False,
            'all_starting_indices_claimed':False,'counterexample_claimed':False,
            'theorem':'FIRSHIFTEDINDEX20261116'}
