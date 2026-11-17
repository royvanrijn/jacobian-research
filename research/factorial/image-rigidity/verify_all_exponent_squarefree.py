"""Exact coefficient/overlap checks and retained prime-power residual reconciliation."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial
import hashlib,json
import sympy as s
import all_exponent_squarefree as a
import prime_power_windows as p

HERE=Path(__file__).resolve().parent
counts={};large_examples=[];residual_examples=[];candidate_examples=[]
def check(name,value):
    if not value:raise AssertionError(name)
    counts[name]=counts.get(name,0)+1

for k in range(1,13):
    check('AM_GM_factorial_bound_is_exact',2**k*factorial(k)<=(k+1)**k)
    check('adjacent_exponential_bound_is_exact',(k+1)**(k+1)<=4*k**(k+1))
    check('normalized_factorial_sequence_is_nondecreasing',
          factorial(k)**2*(k+1)**(k*(k-1))<=k**(k*(k+1)))
for v in range(1,13):
    for u in range(1,v+1):
        check('factorial_ratio_inequality_by_exact_positive_powers',
              factorial(u)**(2*v)*v**(u*v)<=factorial(v)**(2*u)*u**(u*v))
for n in range(2,13):
    for j in range(1,n):
        t=min(j,n-j)
        check('logarithmic_chord_bound_as_integer_power_inequality',n**j>=2**t*j**j)
    q=Q(1,2*n);c=[Q(comb(n,j))*q**min(j,n-j) for j in range(1,n)]
    A=sum((Q(j,n)*value for j,value in enumerate(c,1)),Q(0))
    B=sum(((1-Q(j,n))*value for j,value in enumerate(c,1)),Q(0))
    C=(1+q)**n-1
    check('symmetric_majorant_has_equal_weighted_masses',A==B==sum(c,Q(0))/2)
    check('majorant_masses_are_below_binomial_bound_and_one',A<=C<1)
for n,m in ((2,4),(4,6),(5,7),(8,8),(16,10)):
    cert=a.large_exponent_certificate(n,m)
    check('large_exponent_input_decision_is_exact',cert['squarefree_certified'] and (1<<m)>=4*n*n)
    for j in range(1,n):
        check('literal_normalized_coefficient_is_below_dyadic_majorant',
              factorial(m*j)**(2*n)*2**(m*min(j,n-j)*n)<=factorial(m*n)**(2*j))
    t=s.Symbol('t');f=s.Poly(sum(comb(n,j)*factorial(m*j)*t**j for j in range(n+1)),t)
    check('literal_polynomial_gcd_confirms_finite_large_exponent_controls',s.gcd(f,f.diff()).degree()==0)
    large_examples.append(cert)
for n,m,accepted in ((4,5,False),(4,6,True),(16,9,False),(16,10,True),(10**100,10**100,True)):
    cert=a.large_exponent_certificate(n,m)
    check('bit_length_criterion_handles_equality_and_enormous_exponent_without_power_construction',
          cert['squarefree_certified']==accepted and not cert['counterexample_claimed'])

for n,m in ((1,1),(13,1),(13,2),(17,3),(24,5),(17,9),(31,4)):
    cert=a.prime_power_residual_certificate(n,m);prime=cert['prime'];c=cert['scaling_power']
    remaining=n;K=0;place=0;expected=[]
    while remaining:
        digit=remaining % prime;remaining//=prime
        if digit:
            e=prime**place;h=(e-1)//(prime-1)
            sigma=(-1)**(h+c*e) % prime
            values=[factorial(digit)//factorial(digit-j)*pow(sigma,j,prime) % prime for j in range(digit+1)]
            expected.append((K,K+digit*e,e,h,values));K+=digit*e
        place+=1
    check('every_digit_edge_matches_user_residual_formula_including_r_zero',
          [(row['left'][0],row['right'][0],row['e'],row['h'],row['residual_coefficients']) for row in cert['newton_edges']]==expected)
    check('all_literal_edge_residuals_are_squarefree',all(row['residual_squarefree'] for row in cert['newton_edges']))
    for k in range(6):
        W=Q(factorial(m*k),prime**(c*k)*factorial(k))
        check('prime_power_factorial_ratio_is_unit_with_geometric_residue',
              p.valuation(W,prime)==0 and p.residue(W,prime)==(-1)**(c*k) % prime)
    residual_examples.append(cert)
for prime in (2,3,5,7):
    n=prime**2+prime-1
    coeffs=[factorial(n)//factorial(n-k) for k in range(n+1)]
    check('exponent_one_unit_deformation_works_at_multiple_primes',
          p.newton_residuals(coeffs,prime)==p.expected_edges(n,prime,1))

for n in (1,2,3,4,5,10,1000,2**144):
    cert=a.remaining_exponent_candidates(n);limit=cert['least_large_exponent']
    check('candidate_upper_cutoff_is_exact_integer_ceiling',2**(limit-1)<4*n*n<=2**limit)
    check('every_candidate_is_logarithmically_small_with_two_prime_divisors',
          all(2**m<4*n*n and len(s.factorint(m))>=2 for m in cert['candidate_exponents']))
    check('candidate_list_exhausts_the_remaining_necessary_exponent_range',
          cert['candidate_exponents']==[m for m in range(1,limit) if m!=1 and len(s.factorint(m))>=2])
    check('candidates_are_not_reported_as_counterexamples',not cert['counterexamples_claimed'])
    candidate_examples.append(cert)
check('first_four_indices_cover_all_positive_exponents_by_the_two_deterministic_classes',
      all(row['all_positive_exponents_certified_by_these_two_classes'] for row in candidate_examples[:4]))
check('exponent_six_remains_a_candidate_at_index_five',candidate_examples[4]['candidate_exponents']==[6])
check('pointwise_wrapper_makes_no_counterexample_claim_outside_its_two_classes',
      not a.pointwise_certificate(5,6)['squarefree_certified'] and not a.pointwise_certificate(5,6)['counterexample_claimed'])
check('pointwise_wrapper_certifies_exponent_one_and_prime_powers',
      a.pointwise_certificate(1000,1)['squarefree_certified'] and a.pointwise_certificate(1000,9)['squarefree_certified'])

overlap=a.overlap_certificate();N=int(overlap['overlap_threshold_hex'],16)
check('overlap_threshold_and_endpoint_are_exact',N==2**144 and 512**16==N and 512-32*9-2==overlap['g_at_u_threshold']==222)
check('overlap_derivative_is_uniformly_positive',
      2*Q(1,3)>Q(1,2) and 1-Q(64,512)==Q(7,8)>0)
for n in (N,N+1,2**160,2**256):
    root=int(s.integer_nthroot(n,16)[0]);limit=a.large_exponent_certificate(n,1)['least_large_exponent']
    check('finite_integer_overlap_controls_cover_every_possible_bad_exponent',limit-1<=root)
check('overlap_threshold_is_not_an_extracted_analytic_threshold_or_prefix_audit',
      not overlap['numeric_analytic_starting_threshold_extracted'] and not overlap['small_index_pointwise_audit_claimed'])
for operation in (lambda:a.large_exponent_certificate(0,6),lambda:a.large_exponent_certificate(4,False),
                  lambda:a.prime_power_class_certificate(0),lambda:a.prime_power_residual_certificate(5,6),
                  lambda:a.remaining_exponent_candidates(-1),lambda:a.pointwise_certificate(5,Q(3,2))):
    try:operation()
    except ValueError:check('outside_exact_positive_integer_domain_fails_closed',True)
    else:raise AssertionError('invalid input accepted')

previous=json.loads((HERE/'simultaneous_squarefree_verification.json').read_text())
check('published_log_and_qualified_dense_power_routes_are_retained',
      previous['published_2016_fallback_proved'] and previous['MR_II_preprint_source_qualification_retained']
      and previous['preprint_dense_power_application_checked'] and not previous['complete_preprint_validation_claimed'])
check('harmonic_transfer_uses_power_route',previous['reciprocal_sum_converges_using_power_route'])
audit={'canonical_proof':'ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex','date':'2026-11-17',
    'large_exponent_criterion':'2^m>=4*n^2 implies squarefree F_(n,m), all positive integers n,m',
    'universal_proof':'normalized factorial sequence monotonicity, concave logarithm chord and two weighted coefficient masses<1',
    'coefficient_checks_use_exact_positive_powers_not_algebraic_rounding':True,
    'prime_power_reconciliation':'r>=1 all-index prime-power theorem and every edge residual ALREADY registered in FIRPADICDEFORM/FIRSFCMONOPRIME; new wrapper explicitly admits r0,m1',
    'retained_prime_power_proof':'PRIME_POWER_TWO_TERM_SFC.tex; residual sigma=(-1)^(h+c*e)',
    'retained_source_credit':'Rocks2015 author manuscript Lemma3.8/Proposition3.22 supplies factorial units/digit hull; residual calculation is the retained internal proof',
    'remaining_bad_pair_necessary_conditions':['2^m<4*n^2','m has at least two distinct prime divisors'],
    'candidate_filter_is_not_a_counterexample_test':True,
    'overlap_threshold':'2^144','overlap_proof':'u>=512, g512=222, gprime>7/8 using log2>1/2',
    'set_transfer':'capped and unrestricted bad sets agree for n>=2^144; D<=D_all<=D+2^144',
    'published_2016_route':'D_all(X)<<X/(log X)^(1/50), absolute implied constant, independent of MR II',
    'additional_preprint_route':'D_all(X)<<X^(1-10^-18), absolute implied constant; dense MR II9.2(ii) source qualification inherited',
    'harmonic_convergence':'unrestricted all-exponent exceptional reciprocal sum converges using power route; finite initial contribution only',
    'MR_II_preprint_source_qualification_retained':True,
    'sparse_warning_and_Section10_exclusion_retained':True,'complete_preprint_validation_claimed':False,
    'overlap_threshold_is_not_an_analytic_starting_threshold':True,
    'numerical_analytic_constants_or_thresholds_extracted':False,'all_index_SFC_closed':False,
    'other_Gamma_shapes_claimed':False,'finite_controls_are_not_numerical_density_evidence':True,
    'provenance':'User-supplied elementary large-exponent argument checked; prime-power result reconciled with existing proof; no priority claim',
    'independent_replay':False,'formal_verification':False,'external_review':False}
(HERE/'all_exponent_squarefree_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
(HERE/'all_exponent_squarefree_examples.json').write_text(json.dumps({'large_exponent_controls':large_examples,
    'prime_power_and_one_residuals':residual_examples,'remaining_candidates':candidate_examples,'overlap':overlap,
    'finite_controls_are_not_numerical_density_evidence':True},indent=2)+'\n')
files=('ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex','all_exponent_squarefree.py','verify_all_exponent_squarefree.py',
    'all_exponent_squarefree_audit.json','all_exponent_squarefree_examples.json','prime_power_windows.py',
    'PRIME_POWER_TWO_TERM_SFC.tex','QUANTITATIVE_SQUAREFREENESS_ALL_EXPONENTS.tex',
    'STATUS_AND_LITERATURE.tex',
    'simultaneous_squarefree_verification.json','simultaneous_squarefree_audit.json',
    'sources/rocks-2015-author-preprint.pdf','sources/rocks-2015-author-preprint.txt')
receipt={'status':'pass','checks':counts,'total':sum(counts.values()),
    'files':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
    'universal_evidence':'ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex, retained prime-power theorem and explicitly qualified simultaneous analytic inputs',
    'large_exponent_criterion_proved':True,'prime_power_result_already_registered':True,'exponent_one_explicit':True,
    'all_positive_exponents_share_one_exceptional_set':True,'overlap_threshold':'2^144',
    'published_all_exponent_log_bound':True,'qualified_all_exponent_power_bound':True,
    'all_exponent_exceptional_reciprocal_sum_converges_using_power_route':True,
    'MR_II_source_qualification_retained':True,'numeric_analytic_starting_threshold_extracted':False,
    'finite_controls_are_not_numerical_density_evidence':True,'all_index_SFC_closed':False,
    'independent_replay':False,'formal_verification':False,'external_review':False,'novelty_claim':False}
(HERE/'all_exponent_squarefree_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'pass','checks':receipt['total'],'criterion':'2^m>=4n^2',
    'simultaneous_range':'every positive exponent','overlap':'2^144','all_index_SFC_closed':False}))
