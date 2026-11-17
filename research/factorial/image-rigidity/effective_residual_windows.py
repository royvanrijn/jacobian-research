"""Certified Dusart cutoffs and full nearby-prime residual degree budgets.

Ordinary factorial shape one. Finite exceptions to a cutoff are not
counterexamples; residual-degree intersections remain an all-index target.
"""
from fractions import Fraction
from functools import lru_cache
from math import isqrt,factorial
import hashlib,json
import sympy as s
import prime_power_windows as p
import prime_gap_windows as g
import mixed_prime_windows as m


def _ceil_div(a,b):return (a+b-1)//b


def _log_series(t,terms):
    power=t;total=Fraction(0);square=t*t
    for j in range(terms):
        total+=power/(2*j+1)
        power*=square
    lower=2*total
    upper=lower+2*power/((2*terms+1)*(1-square))
    return lower,upper


@lru_cache(maxsize=None)
def _log_two(terms):return _log_series(Fraction(1,3),terms)


def log_integer_interval(value,terms=16):
    """Exact atanh-series bounds; no transcendental library rounding."""
    value,terms=p.integer(value,1),p.integer(terms,1)
    k=value.bit_length()-1;power=1<<k
    t=Fraction(value-power,value+power)
    lo,hi=_log_series(t,terms)
    l2,u2=_log_two(terms)
    return k*l2+lo,k*u2+hi


def _pack_fraction(value):
    return {'numerator_hex':hex(value.numerator),'denominator_hex':hex(value.denominator)}


def _log_condition(a,value):
    terms=8
    while True:
        lo,hi=log_integer_interval(value,terms)
        if 50*lo*lo>=a:decision=True;break
        if 50*hi*hi<a:decision=False;break
        terms*=2
    return decision,{'integer_argument':value,'series_terms':terms,
                     'lower':_pack_fraction(lo),'upper':_pack_fraction(hi),
                     'condition':'50*log(argument)^2>=a','condition_holds':decision}


def cutoff_certificate(a):
    """Compute the exact proposed P0 using rational comparisons and integers."""
    a=p.integer(a,15);b=a-2
    k=isqrt(a//50)
    if 50*k*k<a:k+=1
    constant_lower=_ceil_div(400000,b)
    integer_upper=_ceil_div(max(400000,3**k),b)
    if integer_upper==constant_lower:
        cutoff=constant_lower
        true_witness={'method':'exp(sqrt(a/50))<3^ceil_sqrt<=b*cutoff',
                      'ceil_sqrt':k,'integer_exp_upper':3**k}
    else:
        passes,witness=_log_condition(a,b*constant_lower)
        if passes:
            cutoff=constant_lower;true_witness=witness
        else:
            lo,hi=constant_lower,integer_upper
            while hi-lo>1:
                mid=(lo+hi)//2
                passes,_=_log_condition(a,b*mid)
                if passes:hi=mid
                else:lo=mid
            cutoff=hi
            passes,true_witness=_log_condition(a,b*cutoff)
            assert passes
    previous_argument=b*(cutoff-1)
    if previous_argument<400000:
        false_witness={'method':'previous integer fails400000 lower bound','integer_argument':previous_argument}
    else:
        passes,false_witness=_log_condition(a,previous_argument)
        assert not passes
    minimum_prime=int(s.nextprime(6*a))
    max_exponent=0;power=1
    while power*minimum_prime<cutoff:
        power*=minimum_prime;max_exponent+=1
    return {'multiplier':a,'exact_cutoff':cutoff,
            'formula':'ceil(max(400000,exp(sqrt(a/50)))/(a-2))',
            'exact_ceiling_certified':True,'constant_lower_cutoff':constant_lower,
            'true_witness':true_witness,'previous_integer_false_witness':false_witness,
            'minimum_allowed_prime':minimum_prime,
            'prime_requirement_implies_cutoff':minimum_prime>=cutoff,
            'finite_exceptions_empty':minimum_prime>=cutoff,
            'finite_exception_prime_base_upper':cutoff-1 if minimum_prime<cutoff else None,
            'finite_exception_exponent_upper':max_exponent,
            'exception_condition':'prime p>6a and p^r<exact_cutoff; not a counterexample claim',
            'all_offsets':list(range(13)),'gamma_shape':'1','theorem':'FIREFFECTIVECUTOFF20261116'}


@lru_cache(maxsize=1)
def multiplier_log_endpoint_certificate():
    """Exact rational witnesses; concavity extends endpoints through a=24000."""
    e_upper=sum((Fraction(1,factorial(j)) for j in range(9)),Fraction(0))+Fraction(10,9*factorial(9))
    e_comparison=Fraction(2719,1000)
    power_upper=e_comparison**22*Fraction(12,13)
    assert e_upper<e_comparison<3
    assert 3**6<6*15*13 and power_upper<6*24000*23998
    right_H=50*Fraction(263,12)**2-24000
    assert right_H==Fraction(1225,72)>0
    return {'multiplier_interval':[15,24000],
            'exponential_partial_sum_last_index':8,
            'exponential_tail_upper':_pack_fraction(Fraction(10,9*factorial(9))),
            'e_rational_upper':_pack_fraction(e_upper),
            'e_comparison_upper':_pack_fraction(e_comparison),
            'left_log_strict_lower':_pack_fraction(Fraction(6)),
            'left_log_integer_argument':1170,'left_H_strict_lower':_pack_fraction(Fraction(1785)),
            'right_log_strict_lower':_pack_fraction(Fraction(263,12)),
            'right_log_integer_argument':3455712000,
            'right_exponential_rational_upper':_pack_fraction(power_upper),
            'right_H_strict_lower':_pack_fraction(right_H),
            'analytic_extension':'H(a)=50*log(6a(a-2))^2-a is strictly concave for real a>=3',
            'all_multipliers_in_interval_log_condition':True,'outside_interval_claimed':False}


@lru_cache(maxsize=1)
def initial_prime_gap_certificate():
    """One explicit finite Eratosthenes calculation; no factorial-polynomial scan."""
    limit=500000
    flags=bytearray(b'\x01')*(limit+1);flags[:2]=b'\x00\x00'
    for q in range(2,isqrt(limit)+1):
        if flags[q]:flags[q*q:limit+1:q]=b'\x00'*((limit-q*q)//q+1)
    primes=[j for j in range(2,limit+1) if flags[j]]
    largest_gap=max(b-a for a,b in zip(primes,primes[1:]))
    endpoints=[[a,b] for a,b in zip(primes,primes[1:]) if b-a==largest_gap]
    largest_distance=max(largest_gap-1,limit-primes[-1])
    return {'limit':limit,'algorithm':'Eratosthenes, starting each crossing at q^2',
            'prime_count':len(primes),'primes':primes,
            'prime_list_sha256':hashlib.sha256(json.dumps(primes,separators=(',',':')).encode()).hexdigest(),
            'maximum_consecutive_prime_gap':largest_gap,'maximum_gap_endpoints':endpoints,
            'maximum_n_minus_preceding_prime':largest_distance,
            'last_prime':primes[-1],'terminal_distance':limit-primes[-1],
            'scope':'finite prime-gap evidence for2<=n<=500000; not squarefreeness of every such F_n',
            'all_index_SFC_claimed':False}


def bounded_multiplier_family_certificate(a,prime,offset):
    """All r>=1, all admissible p, multiplier1..24000 and offset0..12."""
    a,prime,i=p.integer(a,1),p.integer(prime,2),p.integer(offset)
    if a>24000 or i>12 or not s.isprime(prime) or prime<=6*max(a,i):
        raise ValueError('a1..24000, offset0..12, prime>6max(a,i) required')
    return {'status':'complete_bounded_multiplier_family_certified',
            'multiplier':a,'prime':prime,'offset':i,
            'terminal_indices':f'{a}*{prime}^r+{i}, EVERY integer r>=1',
            'multiplier_range':[1,24000],'all_offsets':list(range(13)),
            'all_powers_r_ge1_certified':True,'finite_exceptions_remaining':False,
            'method':'retained uniform family' if a<=14 else 'logarithmic concavity, Dusart and finite prime gaps',
            'initial_range':{'condition':'x=(a-2)P+i<400000, a>=15',
                             'n_strict_upper':461539,'prime_gap_certificate_limit':500000,
                             'preceding_prime_gap_bound':113,'twice_P_strict_lower':180},
            'large_range':{'condition':'x>=400000, a15..24000',
                           'prime_interval':'n-2P<q<n, Dusart Proposition6.8'},
            'square_divisor_degree':'<P, then horizontal degree<=1, then parity',
            'gamma_shape':'1','all_multipliers_or_offsets_claimed':False,
            'all_index_SFC_claimed':False,'theorem':'FIRBOUNDEDMULTIPLIER20261116'}


def effective_shifted_certificate(a,prime,power,offset):
    a,prime,power,i=p.integer(a,1),p.integer(prime,2),p.integer(power,1),p.integer(offset)
    if not s.isprime(prime) or i>12 or prime<=6*max(a,i):
        raise ValueError('prime>6max(a,i), offset0..12 and positive power required')
    P=prime**power
    if a<=24000:
        family=bounded_multiplier_family_certificate(a,prime,i)
        return {'status':'uniform_family_certified','multiplier':a,'prime':prime,'power':power,
                'offset':i,'P':P,'terminal_index':a*P+i,'cutoff':1,'squarefree_certified':True,
                'cutoff_needed':False,'finite_exceptions_remaining':False,
                'family_certificate':family,'all_index_SFC_claimed':False,
                'theorem':'FIRBOUNDEDMULTIPLIER20261116'}
    threshold=cutoff_certificate(a)
    accepted=P>=threshold['exact_cutoff']
    return {'status':'effective_shifted_window_certified' if accepted else 'below_effective_cutoff_no_decision',
            'multiplier':a,'prime':prime,'power':power,'offset':i,'P':P,'terminal_index':a*P+i,
            'cutoff':threshold['exact_cutoff'],'cutoff_certificate':threshold,
            'squarefree_certified':accepted,'counterexample_claimed':False,
            'finite_exceptions_explicitly_bounded':True,'gamma_shape':'1',
            'all_index_SFC_claimed':False,'theorem':'FIREFFECTIVECUTOFF20261116'}


def residual_degree_data(coefficients,prime):
    """Full generic square-part divisor degrees; no moment parity assumed."""
    prime=p.integer(prime,2)
    if not s.isprime(prime):raise ValueError('prime required')
    coefficients=list(coefficients)
    if any(isinstance(a,bool) or not isinstance(a,(int,s.Integer)) for a in coefficients):
        raise ValueError('exact integer coefficients required')
    coefficients=[int(a) % prime for a in coefficients]
    if not coefficients or not coefficients[0] or not coefficients[-1]:
        raise ValueError('nonzero constant and leading coefficients required')
    degree=len(coefficients)-1
    if not degree:
        return {'prime':prime,'coefficients':coefficients,'factor_branches':[],
                'square_part_coefficients':[1],'possible_divisor_degrees':[0],
                'parity_applied':False,'characteristic_zero_squarefree_claimed':False}
    edge={'left':[0,0],'right':[degree,0],'e':1,'h':0,'residual_coefficients':coefficients}
    data=m.square_divisor_degrees([edge],prime)
    x=s.Symbol('x');K=s.Poly(1,x,modulus=prime)
    for branch in data['branches']:
        phi=s.Poly(sum(a*x**j for j,a in enumerate(branch['residual_factor'])),x,modulus=prime)
        K*=phi**branch['max_copies_in_square_divisor']
    return {'prime':prime,'coefficients':coefficients,'factor_branches':data['branches'],
            'square_part_coefficients':[int(K.nth(j)) % prime for j in range(K.degree()+1)],
            'possible_divisor_degrees':data['possible_degrees'],'parity_applied':False,
            'characteristic_zero_squarefree_claimed':False}


def moment_kernel_degree_data(offset,prime):
    offset,prime=p.integer(offset),p.integer(prime,2)
    if not s.isprime(prime) or prime<=6*max(1,offset):
        raise ValueError('prime>6max(1,offset) required')
    coefficients=[1];binomial=1;factorial_unit=1
    for j in range(1,offset+1):
        binomial=binomial*(offset-j+1)*pow(j,-1,prime) % prime
        for k in range(6*(j-1)+1,6*j+1):factorial_unit=factorial_unit*k % prime
        coefficients.append(binomial*factorial_unit % prime)
    data=residual_degree_data(coefficients,prime)
    data['moment_kernel_offset']=offset
    return data


def nearby_residual_certificate(n,primes):
    """Explicit finite nearby primes; full degree sets, with global even bounds."""
    n=p.integer(n,1)
    primes=list(dict.fromkeys(p.integer(q,2) for q in primes))
    if not primes:raise ValueError('provide at least one finite nearby prime')
    profiles=[]
    if n>=30:
        gap=g.prime_gap_certificate(n);band=g.fourteenth_certificate(n)
        maximum=min(gap['even_square_divisor_degree_bound'],band['even_square_divisor_degree_bound'])
    else:
        maximum=2*(n//4)
    allowed={0}|set(range(2,maximum+1,2))
    for q in primes:
        d=n-q
        if d<0 or not s.isprime(q) or q<=6*max(1,d):
            raise ValueError('each prime must satisfy n=q+d and q>6max(1,d)')
        data=moment_kernel_degree_data(d,q)
        allowed.intersection_update(data['possible_divisor_degrees'])
        profiles.append({'nearby_prime':q,'offset':d,**data})
    accepted=allowed=={0}
    return {'status':'nearby_residual_window_certified' if accepted else 'residual_degree_intersection_unresolved',
            'terminal_index':n,'moment_indices':[n-1,n],'profiles':profiles,
            'global_even_degree_bound':maximum,'remaining_degrees':sorted(allowed),
            'squarefree_certified':accepted,'counterexample_claimed':False,
            'uniform_all_indices_coverage_claimed':False,'gamma_shape':'1',
            'theorem':'FIRRESIDUALDEGREES20261116'}


def fixed_prime_offset_family(prime,offset):
    """Full residual degrees certify every r>=1 at one fixed prime/offset."""
    data=moment_kernel_degree_data(offset,prime)
    even=[d for d in data['possible_divisor_degrees'] if d % 2==0]
    accepted=even==[0]
    return {'status':'fixed_prime_all_r_family_certified' if accepted else 'fixed_prime_degree_test_unresolved',
            'prime':prime,'offset':offset,'kernel_data':data,'remaining_even_degrees':even,
            'all_r_ge1_certified':accepted,'all_primes_or_offsets_claimed':False,
            'uniform_C_smoothness_claimed':False,'counterexample_claimed':False,
            'theorem':'FIRRESIDUALDEGREES20261116'}
