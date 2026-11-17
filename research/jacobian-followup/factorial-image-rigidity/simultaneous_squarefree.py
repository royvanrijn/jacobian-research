"""Exact parameters and finite common offsets; analytic inputs are source-qualified."""
from fractions import Fraction as Q
from math import factorial
import sympy as s
import prime_power_windows as p


def rational(value):return {'numerator':value.numerator,'denominator':value.denominator}


def parameter_certificate():
    nu1,nu2,eta,rho,delta=Q(1,100000),Q(3,20),Q(1,150),Q(1,10),Q(1,2000)
    h_power,m_power,p1_power=Q(1,5),Q(1,16),Q(1,10**8)
    rho1_lower=Q(1,9)  # pi>3 gives rho1>1/9.
    sqrt_ratio_upper=Q(49,6000)
    sieve_upper=Q(6,5)*(Q(1,1000)+2*sqrt_ratio_upper)
    complement_loss=sieve_upper+2*delta
    e_upper=sum((Q(1,factorial(j)) for j in range(9)),Q(0))+Q(10,9*factorial(9))
    gap=Q(1,20);witness_margin=gap-complement_loss
    saving1=p1_power*(Q(1,2)-nu2-3*eta);saving2=nu1**3/200
    incidence_power=Q(1,2)+Q(3,2)*m_power+2*h_power
    eps=Q(1,10**18)
    return {'nu1':rational(nu1),'nu2':rational(nu2),'eta':rational(eta),'rho':rational(rho),
            'delta':rational(delta),'rho1_strict_lower':rational(rho1_lower),'J_MR':1,
            'M_power':rational(m_power),'H_power':rational(h_power),'P1_power':rational(p1_power),
            'middle_interval_power_squared':rational(nu1*nu2),
            'source_P1_log_power':rational(2/eta),'dense_error_log_power':rational(rho/12),
            'sqrt_nu1_over_nu2_squared':rational(nu1/nu2),
            'sqrt_ratio_strict_upper':rational(sqrt_ratio_upper),
            'sieve_complement_upper':rational(sieve_upper),'total_loss_upper':rational(complement_loss),
            'e_strict_upper':rational(e_upper),'e_comparison':rational(Q(2719,1000)),
            'exp_one_twentieth_strict_upper':rational(Q(20,19)),
            'Delta_strict_lower':rational(gap),'witness_margin_strict_lower':rational(witness_margin),
            'usable_offset_fraction':rational(Q(1,50)),
            'dense_first_saving':rational(saving1),'dense_second_saving':rational(saving2),
            'incidence_error_power':rational(incidence_power),'final_epsilon':rational(eps),
            'J_density_log_argument':rational(Q(32,21)),'smooth_complement_log_argument':rational(Q(32,17)),
            'margin_log_argument':rational(Q(1024,357)),
            'published_2016_first_saving':rational(h_power*Q(1,100)/25),
            'analytic_exception_multiplied_by_M':False,'numeric_starting_threshold_extracted':False,
            'power_saving_input':'MR II arXiv2007.04290v1 Theorem9.2(ii), dense; preprint',
            'log_saving_input':'published MR2016 Theorem1; independent of MR II',
            'preprint_Theorem1_7_advertised_exponent_used':False,'Section10_sparse_sieve_used':False}


def dyadic_parameters(X):
    X=p.integer(X,2)
    M=int(s.integer_nthroot(2*X,16)[0]);H=int(s.integer_nthroot(X,5)[0])
    return {'X':X,'M':M,'H':H,'Y_squared':2*M*X,'Z_squared':M*H,
            'finite_parameters_only':True,'analytic_asymptotic_conditions_claimed_at_this_X':False}


def common_offset_certificate(H,M,offset):
    H,M,i=p.integer(H,1),p.integer(M,1),p.integer(offset,1)
    if i>H or i<=(M+1)**2:raise ValueError('offset<=H and offset>(M+1)^2 required')
    q=int(s.prevprime(i+1));d=i-q
    largest=int(max(s.factorint(i)))
    if d*d>=i or largest*largest<=M*H:
        raise ValueError('gap<sqrt(offset) and Pplus(offset)>sqrt(MH) required')
    assert q>M*max(1,d) and largest>d
    return {'H':H,'maximum_exponent':M,'offset':i,'largest_prime_factor':largest,
            'preceding_prime':q,'preceding_gap':d,'factor_degree_divisibility':largest,
            'large_factor_degree_at_least':q,'complement_degree_at_most':d,
            'irreducible_for_every_exponent_1_through_M':True,
            'uniform_local_condition':q>M*max(1,d),
            'parity_used':False,'numerical_density_evidence':False,
            'all_offsets_claimed':False,'theorem':'FIRCOMMONOFFSETS20261117'}
