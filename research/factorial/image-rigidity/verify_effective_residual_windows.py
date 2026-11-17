"""Exact cutoff witnesses and finite full-residual controls; no all-index scan."""
from pathlib import Path
from fractions import Fraction
from math import comb,factorial
import hashlib,json
import sympy as s
import effective_residual_windows as e

HERE=Path(__file__).resolve().parent
counts={};cutoffs=[];examples=[]
def check(name,value):
    if not value:raise AssertionError(name)
    counts[name]=counts.get(name,0)+1

def unpack(row):
    return Fraction(int(row['numerator_hex'],16),int(row['denominator_hex'],16))

def independent_log_bracket(value,terms):
    k=value.bit_length()-1
    def series(t):
        lo=2*sum((t**(2*j+1)/Fraction(2*j+1) for j in range(terms)),Fraction(0))
        tail=2*t**(2*terms+1)/((2*terms+1)*(1-t*t))
        return lo,lo+tail
    l2,u2=series(Fraction(1,3))
    lr,ur=series(Fraction(value-2**k,value+2**k))
    return k*l2+lr,k*u2+ur

# Exact endpoint certificate plus a second, structurally different finite sieve.
endpoints=e.multiplier_log_endpoint_certificate()
check('e_upper_bound_is_finite_sum_plus_positive_tail',
      unpack(endpoints['e_rational_upper'])==sum((Fraction(1,factorial(j)) for j in range(9)),Fraction(0))+Fraction(10,9*factorial(9)))
check('e_upper_strictly_below_rational_comparison',
      unpack(endpoints['e_rational_upper'])<Fraction(2719,1000)<3)
check('left_endpoint_log_strict_bound',3**6<6*15*13==endpoints['left_log_integer_argument'])
check('right_endpoint_log_strict_bound',
      Fraction(2719,1000)**22*Fraction(12,13)==unpack(endpoints['right_exponential_rational_upper'])<6*24000*23998)
check('endpoint_H_bounds_are_positive_exact_rationals',
      50*6**2-15==unpack(endpoints['left_H_strict_lower'])==1785
      and 50*Fraction(263,12)**2-24000==unpack(endpoints['right_H_strict_lower'])==Fraction(1225,72)>0)
A=s.Symbol('a',positive=True);L=s.log(6*A*(A-2));S=A**-2+(A-2)**-2
check('concavity_second_derivative_formula',s.simplify(s.diff(50*L**2-A,A,2)-100*(s.diff(L,A)**2-L*S))==0)
check('concavity_Cauchy_defect_is_a_negative_square',
      s.simplify(s.diff(L,A)**2-2*S+(1/A-1/(A-2))**2)==0)
check('L_exceeds_two_on_a_ge_three',3**2<6*3*(3-2))

finite=e.initial_prime_gap_certificate();limit=finite['limit']
least=[0]*(limit+1);replay=[]
for j in range(2,limit+1):
    if not least[j]:least[j]=j;replay.append(j)
    for q in replay:
        if q>least[j] or j*q>limit:break
        least[j*q]=q
check('Eratosthenes_and_Euler_sieves_agree_on_every_prime',finite['primes']==replay)
check('prime_count_through_500000',len(replay)==finite['prime_count']==41538)
check('prime_list_hash_is_reproducible',
      hashlib.sha256(json.dumps(replay,separators=(',',':')).encode()).hexdigest()==finite['prime_list_sha256'])
current=2;maximum=0;witnesses=[]
for j in range(2,limit+1):
    if least[j]==j:current=j
    distance=j-current
    if distance>maximum:maximum=distance;witnesses=[j]
    elif distance==maximum:witnesses.append(j)
check('preceding_prime_distances_checked_at_every_integer',
      maximum==finite['maximum_n_minus_preceding_prime']==113 and witnesses==[492226])
check('largest_consecutive_gap_and_unique_endpoints',
      max(b-a for a,b in zip(replay,replay[1:]))==finite['maximum_consecutive_prime_gap']==114
      and [[a,b] for a,b in zip(replay,replay[1:]) if b-a==114]==finite['maximum_gap_endpoints']==[[492113,492227]])
check('finite_gap_tail_included',finite['last_prime']==499979 and finite['terminal_distance']==21)
check('small_branch_index_strict_upper_and_domain',
      Fraction(400000*15,13)<461539<=limit and 15*97>30 and 113<180)
check('finite_sieve_is_not_an_all_index_squarefree_claim',not finite['all_index_SFC_claimed'])
for a,q,i in ((1,7,0),(1,73,12),(14,89,12),(15,97,0),(15,97,12),
               (100,601,12),(1000,6007,12),(24000,144013,12)):
    family=e.bounded_multiplier_family_certificate(a,q,i)
    check('complete_bounded_multiplier_family_has_no_initial_exceptions',
          family['all_powers_r_ge1_certified'] and not family['finite_exceptions_remaining']
          and not family['all_index_SFC_claimed'])

for value in (1,2,3,400000,425905,10**20):
    previous=None
    for terms in (8,16):
        actual=e.log_integer_interval(value,terms)
        check('log_series_replayed_from_independent_finite_sum',actual==independent_log_bracket(value,terms))
        check('log_interval_order',0<=actual[0]<=actual[1])
        if previous:check('log_intervals_nest_under_refinement',previous[0]<=actual[0]<=actual[1]<=previous[1])
        previous=actual

for a in (15,100,1000,10000,100000):
    row=e.cutoff_certificate(a);m=row['exact_cutoff'];b=a-2
    check('constant_cutoff_ceiling',row['constant_lower_cutoff']==(400000+b-1)//b)
    check('cutoff_constant_true_bound',b*m>=400000)
    for truth,witness in ((True,row['true_witness']),(False,row['previous_integer_false_witness'])):
        if 'series_terms' in witness:
            arg=witness['integer_argument'];terms=witness['series_terms']
            lo,hi=independent_log_bracket(arg,terms)
            check('retained_log_witness_reconstructed',lo==unpack(witness['lower']) and hi==unpack(witness['upper']))
            check('retained_log_comparison_is_rational_and_strict_when_false',
                  50*lo*lo>=a if truth else 50*hi*hi<a)
            check('retained_log_argument_at_ceiling_boundary',arg==b*(m if truth else m-1))
        elif truth:
            k=witness['ceil_sqrt']
            check('integer_upper_exponent_is_ceiling_sqrt',50*(k-1)**2<a<=50*k*k)
            check('integer_exponential_upper_dominates_exact_threshold',witness['integer_exp_upper']==3**k<=b*m)
        else:check('previous_integer_fails_constant_threshold',witness['integer_argument']==b*(m-1)<400000)
    q=row['minimum_allowed_prime'];r=row['finite_exception_exponent_upper']
    check('minimum_allowed_prime_is_exact',s.isprime(q) and q==s.nextprime(6*a))
    check('finite_exception_base_and_exponent_bounds',q**(r+1)>=m and (r==0 or q**r<m))
    check('automatic_prime_condition_case_is_explicit',row['finite_exceptions_empty']==(q>=m))
    cutoffs.append(row)
check('three_proposed_cutoffs_reproduced_exactly',[r['exact_cutoff'] for r in cutoffs[:3]]==[30770,4082,401])
check('nonconstant_exponential_branch_is_exercised',cutoffs[-1]['exact_cutoff']>cutoffs[-1]['constant_lower_cutoff'])
check('multiplier_thousand_has_no_cutoff_exceptions',cutoffs[2]['finite_exceptions_empty'])
for a,q,r,i,accepted in ((15,97,1,12,True),(15,181,2,12,True),(100,4091,1,12,True),
                         (1000,6007,1,12,True),(1,7,1,0,True),(14,89,1,12,True),
                         (24000,144013,1,12,True),(100000,600011,1,12,False),
                         (100000,600011,3,12,True)):
    row=e.effective_shifted_certificate(a,q,r,i)
    check('exported_cutoff_decision_and_domain',row['squarefree_certified']==accepted and not row['all_index_SFC_claimed'])
    examples.append(row)
P=181**2;a=15;i=12;n=a*P+i;x=(a-2)*P+i;q=int(s.nextprime(x))
lo,hi=independent_log_bracket(x,16)
check('accepted_cutoff_sample_meets_strict_prime_interval',x<q<n and n-q<2*P)
check('Dusart_sample_width_verified_with_rational_log_upper',25*(q-x)*hi*hi<=x)
check('cutoff_proof_strict_width_chain',50*lo*lo>=a and x<a*P and 2*Fraction(x,a)<2*P)

X=s.Symbol('x')
def poly(coeffs,q):return s.Poly(sum(c*X**j for j,c in enumerate(coeffs)),X,modulus=q)
def verify_residual(data):
    q=data['prime'];f=poly(data['coefficients'],q);K=poly(data['square_part_coefficients'],q)
    check('square_part_reconstructs_its_degree_and_divides_twice',f.rem(K*K).is_zero and max(data['possible_divisor_degrees'])==K.degree())
    degrees={0};reconstruction=s.Poly(1,X,modulus=q)
    for branch in data['factor_branches']:
        phi=poly(branch['residual_factor'],q)
        check('residual_factor_is_verified_irreducible',branch['factor_irreducibility_verified'] and phi.is_irreducible)
        reconstruction*=phi**branch['multiplicity']
        degrees={d+k*phi.degree() for d in degrees for k in range(branch['max_copies_in_square_divisor']+1)}
    check('all_residual_factors_reconstruct_up_to_scalar',reconstruction.monic()==f.monic())
    check('residual_degree_set_is_full_bounded_sum',sorted(degrees)==data['possible_divisor_degrees'])
    check('generic_residual_does_not_invent_real_parity_or_charzero_claim',
          not data['parity_applied'] and not data['characteristic_zero_squarefree_claimed'])

phi=X**3+X+1;psi=X**2+X+1
synthetic=[]
for expression,expected in ((phi**2,[0,3]),(phi**2*psi**2,[0,2,3,5])):
    f=s.Poly(expression,X,modulus=2)
    row=e.residual_degree_data([int(f.nth(j)) for j in range(f.degree()+1)],2)
    verify_residual(row)
    check('synthetic_residual_preserves_cubic_degree_information',row['possible_divisor_degrees']==expected)
    synthetic.append({'fixture_only_not_a_moment_kernel':True,**row})
actual=[]
for d in (0,1,2,3,12,13,17):
    q=int(s.nextprime(6*max(1,d)));row=e.moment_kernel_degree_data(d,q)
    check('moment_kernel_recurrence_matches_literal_coefficients',
          row['coefficients']==[comb(d,j)*factorial(6*j)%q for j in range(d+1)])
    verify_residual(row);actual.append(row)
nearby=[]
for n,primes in ((15,[13]),(118,[113]),(126,[109,113])):
    row=e.nearby_residual_certificate(n,primes)
    check('nearby_degree_intersection_certifies_sample',row['squarefree_certified'] and row['remaining_degrees']==[0])
    check('nearby_certificate_scope_is_finite',not row['uniform_all_indices_coverage_claimed'])
    for data in row['profiles']:
        q=data['prime'];d=n-q
        check('actual_horizontal_moment_reduction_equals_kernel',
              [comb(n,j)*factorial(6*j)%q for j in range(n+1)]==data['coefficients']+[0]*q)
        verify_residual(data)
    nearby.append(row)
family=e.fixed_prime_offset_family(113,13)
verify_residual(family['kernel_data'])
check('fixed_prime_offset13_family_is_all_r_only',family['all_r_ge1_certified']
      and not family['all_primes_or_offsets_claimed'] and not family['uniform_C_smoothness_claimed'])
for operation in (lambda:e.bounded_multiplier_family_certificate(24001,144013,12),
                  lambda:e.bounded_multiplier_family_certificate(24000,144001,13),
                  lambda:e.bounded_multiplier_family_certificate(15,89,12),
                  lambda:e.cutoff_certificate(14),lambda:e.cutoff_certificate(True),
                  lambda:e.effective_shifted_certificate(1,7,1,12),lambda:e.effective_shifted_certificate(15,91,1,0),
                  lambda:e.effective_shifted_certificate(15,97,0,0),lambda:e.effective_shifted_certificate(15,97,1,13),
                  lambda:e.residual_degree_data([1,True],7),lambda:e.residual_degree_data([1,Fraction(1,2)],7),
                  lambda:e.residual_degree_data([0,1],7),lambda:e.residual_degree_data([1,0],7),
                  lambda:e.residual_degree_data([1,1],9),lambda:e.moment_kernel_degree_data(13,73),
                  lambda:e.nearby_residual_certificate(126,[]),lambda:e.nearby_residual_certificate(126,[107])):
    try:operation()
    except ValueError:check('outside_proved_input_domain_fails_closed',True)
    else:raise AssertionError('invalid input accepted')

(HERE/'initial_prime_gaps_500000.json').write_text(json.dumps(finite,indent=2)+'\n')
(HERE/'multiplier_endpoint_certificate.json').write_text(json.dumps(endpoints,indent=2)+'\n')
(HERE/'effective_cutoff_examples.json').write_text(json.dumps({'cutoff_certificates':cutoffs,'window_examples':examples,
    'formula_exactly_evaluated':True,'finite_exception_domains_are_not_counterexamples':True},indent=2)+'\n')
(HERE/'residual_degree_examples.json').write_text(json.dumps({'actual_moment_kernels':actual,'nearby_index_certificates':nearby,
    'synthetic_residual_controls':synthetic,'fixed_prime_offset_family':family,
    'all_offset_C_smoothness_claimed':False,'uniform_nearby_prime_coverage_claimed':False},indent=2)+'\n')
files=('EFFECTIVE_CUTOFFS_AND_RESIDUAL_DEGREES.tex','effective_residual_windows.py','verify_effective_residual_windows.py',
       'prime_power_windows.py','prime_gap_windows.py','mixed_prime_windows.py','newton_slope_coverage.py',
       'sources/dusart-1002.0442.pdf','sources/dusart-1002.0442.txt',
       'effective_cutoff_examples.json','residual_degree_examples.json',
       'initial_prime_gaps_500000.json','multiplier_endpoint_certificate.json')
receipt={'status':'pass','checks':counts,'total':sum(counts.values()),
         'files':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
         'universal_evidence':'EFFECTIVE_CUTOFFS_AND_RESIDUAL_DEGREES.tex; imported pinned Dusart Proposition6.8 and retained prime-gap/kernel proofs',
         'complete_bounded_multiplier_family':{'multipliers':[1,24000],'offsets':[0,12],
             'prime_condition':'p>6max(a,i)','all_r_ge1':True,'finite_exceptions_remaining':False},
         'finite_prime_gap_replay':{'limit':500000,'algorithms':['Eratosthenes','Euler least-prime-factor'],
             'prime_count':41538,'maximum_consecutive_gap':114,'maximum_preceding_distance':113},
         'exact_endpoint_comparisons_verified':True,
         'effective_fixed_multiplier_cutoff_computed':True,'cutoff_examples':{str(r['multiplier']):r['exact_cutoff'] for r in cutoffs},
         'fixed_prime_offset13_example':{'prime':113,'offset':13,'all_r_ge1_certified':True,'uniform_C13_smoothness_claimed':False},
         'C_positive_iff_characteristic_zero_squarefree_for_this_sequence':True,
         'full_residual_divisor_degree_sets_retained':True,'uniform_nearby_prime_coverage_proved':False,
         'finite_initial_range_still_required_for_an_eventual_coverage_proof':True,
         'all_index_SFC_closed':False,'independent_replay':False,'formal_verification':False,
         'external_review':False,'novelty_claim':False}
(HERE/'effective_residual_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'checks':receipt['total'],'cutoffs':receipt['cutoff_examples'],
                  'complete_multiplier_range':[1,24000],'all_index_SFC_closed':False}))
