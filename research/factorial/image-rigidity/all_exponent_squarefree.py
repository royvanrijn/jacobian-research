"""Deterministic large exponents and a finite necessary candidate list per index."""
from fractions import Fraction
from math import comb,factorial
import sympy as s
import prime_power_windows as p


def large_exponent_certificate(n,m):
    n,m=p.integer(n,1),p.integer(m,1)
    K=4*n*n;threshold=(K-1).bit_length()
    accepted=m>=threshold
    return {'status':'large_exponent_squarefree_certified' if accepted else 'criterion_not_met_no_decision',
            'terminal_index':n,'exponent':m,'least_large_exponent':threshold,
            'criterion':'2^m>=4*n^2, evaluated by integer bit length',
            'squarefree_certified':accepted,'equality_allowed':True,
            'weighted_mass_bound':'A,B<exp(1/2)-1<1','counterexample_claimed':False,
            'gamma_shape':'1','theorem':'FIRLARGEM20261117'}


def prime_power_class_certificate(m):
    m=p.integer(m,1)
    if m==1:
        return {'exponent':1,'prime':2,'prime_exponent':0,'scaling_power':0,
                'all_indices_n_ge1_certified':True,'theorem':'FIRPADICDEFORM20261116',
                'mechanism':'falling-factorial unit deformation with every unit1; r0 case'}
    factors=s.factorint(m)
    if len(factors)!=1:
        return {'exponent':m,'all_indices_n_ge1_certified':False,'counterexample_claimed':False,
                'reason':'outside exponent-one or prime-power class'}
    prime,power=next(iter(factors.items()));prime,power=int(prime),int(power)
    previous=p.class_certificate(m)
    return {'exponent':m,'prime':prime,'prime_exponent':power,'scaling_power':(m-1)//(prime-1),
            'all_indices_n_ge1_certified':True,'theorem':previous['theorem'],
            'mechanism':previous['all_index_mechanism'],'prime_power_result_already_registered':True}


def prime_power_residual_certificate(n,m):
    """An explicitly requested finite replay, including exponent one."""
    n=p.integer(n,1);data=prime_power_class_certificate(m)
    if not data['all_indices_n_ge1_certified']:raise ValueError('m=1 or a prime power required')
    m,prime,c=data['exponent'],data['prime'],data['scaling_power']
    coefficients=[Fraction(comb(n,k)*factorial(m*k),prime**(c*k)) for k in range(n+1)]
    rows=p.newton_residuals(coefficients,prime)
    expected=p.expected_edges(n,prime,(-1)**c % prime)
    assert rows==expected and all(row['residual_squarefree'] for row in rows)
    return {**data,'terminal_index':n,'newton_edges':rows,
            'squarefree_certified':True,'finite_replay_is_not_an_all_index_scan':True}


def remaining_exponent_candidates(n):
    """Necessary candidates only; no moment-polynomial or index sweep."""
    n=p.integer(n,1);threshold=(4*n*n-1).bit_length()
    candidates=[m for m in range(2,threshold) if len(s.factorint(m))>=2]
    return {'terminal_index':n,'least_large_exponent':threshold,'candidate_exponents':candidates,
            'necessary_only':True,'counterexamples_claimed':False,
            'all_positive_exponents_certified_by_these_two_classes':not candidates,
            'other_registered_index_certificates_applied':False,
            'criterion':'2^m<4*n^2 and at least two distinct prime divisors',
            'theorem':'FIRLOGEXPCANDIDATES20261117'}


def pointwise_certificate(n,m):
    n,m=p.integer(n,1),p.integer(m,1)
    large=large_exponent_certificate(n,m)
    if large['squarefree_certified']:return large
    prime_power=prime_power_class_certificate(m)
    accepted=prime_power['all_indices_n_ge1_certified']
    return {'status':'prime_power_or_one_squarefree_certified' if accepted else 'remaining_logarithmic_exponent_candidate',
            'terminal_index':n,'exponent':m,'squarefree_certified':accepted,
            'class_certificate':prime_power,'least_large_exponent':large['least_large_exponent'],
            'counterexample_claimed':False,'other_registered_index_certificates_applied':False,
            'all_index_SFC_closed':False,
            'theorem':prime_power.get('theorem','FIRLOGEXPCANDIDATES20261117')}


def overlap_certificate():
    return {'overlap_threshold':'2^144','overlap_threshold_hex':hex(1<<144),
            'overlap_binary_exponent':144,'u_threshold':512,'g_at_u_threshold':222,
            'g_formula':'u-32*log2(u)-2','derivative_strict_lower':{'numerator':7,'denominator':8},
            'log2_strict_lower':{'numerator':1,'denominator':2},
            'all_and_capped_bad_sets_agree_for_n_at_least_threshold':True,
            'count_transfer':'D(X)<=D_all(X)<=D(X)+2^144',
            'numeric_analytic_starting_threshold_extracted':False,
            'small_index_pointwise_audit_claimed':False}
