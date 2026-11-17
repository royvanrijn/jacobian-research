"""Exact finite controls for the general-m quantitative proof; no density sampling."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial,prod
import hashlib,json
import sympy as s
import quantitative_squarefree as a
import prime_power_windows as p

HERE=Path(__file__).resolve().parent
counts={};profiles=[];resultants=[];local=[]
def check(name,value):
    if not value:raise AssertionError(name)
    counts[name]=counts.get(name,0)+1

cases=((1,1,2,0),(1,1,2,1),(1,2,3,2),(2,1,3,1),(2,3,11,4),
       (3,4,13,2),(4,1,5,1),(6,1,7,1),(6,3,19,2),(10,1,11,1),(10,2,31,3))
for m,b,q,i in cases:
    cert=a.shifted_profile(m,b,q,i);n=cert['terminal_index']
    coeffs=[comb(n,j)*factorial(m*j) for j in range(n+1)]
    check('factorial_range_is_strictly_below_prime_square',m*n<q*q and cert['factorial_range_below_p_squared'])
    for j,value in enumerate(coeffs):
        block,residue=divmod(j,q)
        expected=m*block+m*residue//q+int(residue>i)
        check('general_exponent_valuation_matches_literal_integer_coefficients',p.valuation(value,q)==expected)
    rows=p.newton_residuals(coeffs,q)
    check('literal_hull_matches_horizontal_then_general_slope',
          [rows[0]['left']]+[row['right'] for row in rows]==cert['vertices']
          and rows[-1]['e']==q and rows[-1]['h']==m)
    check('full_reduction_is_exactly_offset_kernel',
          [value % q for value in coeffs]==[value % q for value in a.coefficients(m,i)]+[0]*(n-i))
    check('occupancy_certificate_uses_no_parity',not cert['parity_used'])
    profiles.append(cert)

T=s.Symbol('t')
def poly(coeffs,q):return s.Poly(sum(value*T**j for j,value in enumerate(coeffs)),T,modulus=q)
for m,n,q,i,accepted in ((1,3,2,1,True),(2,4,3,1,True),(6,59,19,2,True),(3,40,19,2,False)):
    cert=a.local_squarefree_certificate(m,n,q,i)
    f=poly(cert['residual_coefficients'],q)
    witness=poly(cert['bezout_left'],q)*f+poly(cert['bezout_right'],q)*f.diff()
    check('local_residual_Bezout_reconstructs_normalized_gcd',
          witness==poly(cert['residual_gcd_coefficients'],q)==s.gcd(f,f.diff()))
    gap=cert['preceding_prime_bound'];qn=gap['preceding_prime'];d=n-qn
    check('general_preceding_prime_bound_has_valid_local_hypotheses',
          s.isprime(qn) and qn==s.prevprime(n+1) and qn>m*max(1,d)
          and gap['square_divisor_degree_bound']==d//2)
    check('local_bad_reduction_fails_closed',cert['squarefree_certified']==accepted and not cert['counterexample_claimed'])
    check('local_controls_do_not_infer_density_or_apply_parity',not cert['numerical_density_evidence'] and not cert['parity_used'])
    local.append(cert)

for m,i in ((1,1),(1,3),(2,2),(3,2),(6,3),(10,3)):
    cert=a.resultant_size_certificate(m,i);values=list(map(int,cert['coefficients']))
    f=s.Poly(sum(value*T**j for j,value in enumerate(values)),T,domain=s.ZZ);df=f.diff()
    columns=[s.Poly(T**j*f.as_expr(),T) for j in range(i-2,-1,-1)]
    columns += [s.Poly(T**j*df.as_expr(),T) for j in range(i-1,-1,-1)]
    sylvester=s.Matrix([[column.nth(k) for column in columns] for k in range(2*i-2,-1,-1)])
    R=int(cert['resultant']);A=int(cert['A']);bound=int(cert['hadamard_absolute_bound'])
    check('independent_Sylvester_determinant_matches_resultant',abs(int(sylvester.det()))==abs(R))
    check('Hadamard_row_count_and_bound_are_exact',
          cert['sylvester_row_counts']=={'polynomial':i-1,'derivative':i}
          and A==2**i*factorial(m*i) and bound==i**i*A**(2*i-1))
    check('coefficient_and_derivative_norm_bounds_hold_as_integer_squares',
          sum(value*value for value in values)==int(cert['coefficient_l2_squared'])<=A*A
          and sum(j*j*values[j]**2 for j in range(1,i+1))==int(cert['derivative_l2_squared'])<=i*i*A*A)
    check('literal_resultant_obeys_Hadamard',abs(R)<=bound)
    check('resultant_test_does_not_assume_all_offsets_or_irreducibility',
          not cert['irreducibility_claimed'] and not cert['all_offsets_squarefree_claimed'] and not cert['numerical_density_evidence'])
    resultants.append(cert)

# One small, nonempty bad-prime progression union, not a density experiment.
m=3;X=40;offset=2;R=int(a.resultant_size_certificate(m,offset)['resultant']);Y2=2*m*X
large=[int(q) for q in s.factorint(abs(R)) if int(q)**2>Y2]
check('finite_bad_prime_fixture_is_nonempty',large==[19])
check('large_prime_log_count_bound_follows_from_integer_products',
      prod(large)<=abs(R) and Y2**len(large)<=R*R)
union=set();pairs=[]
for q in large:
    values=[n for n in range(X,2*X+1) if (n-offset) % q==0]
    by_floors=(2*X-offset)//q-(X-1-offset)//q
    check('bad_prime_progression_count_matches_exact_floor_formula',len(values)==by_floors)
    check('bad_prime_progression_count_obeys_X_over_p_plus_one',len(values)<=Q(X,q)+1)
    union.update(values);pairs.append({'offset':offset,'prime':q,'indices':values,'count':len(values)})
check('bad_prime_union_obeys_pair_count_sum',len(union)<=sum(row['count'] for row in pairs) and sorted(union)==[40,59,78])

# Verify the prime-gap exceptional-set containment used before importing its bound.
N=200;bad=[];q=2
for n in range(2,N+1):
    if s.isprime(n):q=n
    if (n-q)**2>=n:
        following=int(s.nextprime(q))
        check('bad_preceding_distance_lies_in_a_long_consecutive_prime_gap',
              following>n and (following-q)**2>q and q<=N)
        bad.append(n)
long_gap_sum=sum(int(s.nextprime(q))-int(q) for q in s.primerange(2,N+1)
                 if (int(s.nextprime(q))-int(q))**2>q)
check('finite_gap_containment_respects_total_gap_length',len(bad)<=long_gap_sum)

h_power=Q(1,8);delta=Q(1,4);rate=Q(1,50)
check('fixed_mean_margin_allows_delta_one_quarter',1-Q(56,81)<Q(1,3) and Q(1,3)+delta+Q(1,24)==Q(5,8)<Q(2,3))
check('MR_polynomial_exception_exponent_is_799_over_800',1-h_power*delta/25==Q(799,800))
check('bad_reduction_union_exponent_is_seven_eighths',Q(1,2)+3*h_power==Q(7,8))
check('remaining_error_terms_have_positive_power_savings',all(v<1 for v in (Q(799,800),Q(7,8),Q(2,3),h_power)))
check('MR_logarithmic_ratio_exponent_for_absorption',Q(1,3)+rate==Q(53,150))
for operation in (lambda:a.shifted_profile(0,1,7,0),lambda:a.shifted_profile(2,1,2,0),
                  lambda:a.shifted_profile(3,4,12,2),lambda:a.shifted_profile(1,2,3,3),
                  lambda:a.preceding_prime_bound(10,10),lambda:a.local_squarefree_certificate(3,40,19,3),
                  lambda:a.resultant_size_certificate(3,0),lambda:a.resultant_size_certificate(True,2)):
    try:operation()
    except ValueError:check('outside_proved_finite_input_domain_fails_closed',True)
    else:raise AssertionError('invalid input accepted')

sources={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
    'sources/filaseta-rocks-2016-published.pdf','sources/filaseta-rocks-2016-published.txt',
    'sources/matomaki-radziwill-2016-annals.pdf','sources/matomaki-radziwill-2016-annals.txt')}
audit={'claim':'FIRQUANTDENSITY20261116','canonical_proof':'QUANTITATIVE_SQUAREFREENESS_ALL_EXPONENTS.tex',
    'scope':'E_m(X)<<_m X/(log X)^(1/50), every FIXED positive integer m, ordinary Gamma shape1',
    'analytic_inputs':[
        {'source':'Filaseta-Rocks2016 Theorem1.1','scope':'irreducible offset lower density>=log2 for each fixed m','imported_not_reproved':True},
        {'source':'Filaseta-Rocks2016 equation3.1 page6','scope':'sum of consecutive gaps exceeding sqrt(left endpoint), with left endpoint<=Z, is O(Z^(2/3)); yields O(X^(2/3)) exceptional indices','imported_not_reproved':True},
        {'source':'Matomaki-Radziwill2016 Annals Theorem1 pp1015-1016','scope':'uniform h,delta,f; delta1/4 and H~X^(1/8) give X^(799/800)*(log X)^(1/3)+X/(log X)^(1/50)','imported_not_reproved':True},
        {'source':'classical PNT, retained DLMF27.2.3 reference','scope':'eventual smooth-indicator mean1-log2, fixed m; no quantitative PNT error needed','imported_not_reproved':True}],
    'local_valuation':'m*b+floor(m*c/p)+indicator(c>i), j=b*p+c; m*n<p^2',
    'local_occupancy':'degree=u+v*p, 2u<=i,2v<=a; no parity assumption',
    'preceding_prime_local_condition':'q_n>m*max(1,n-q_n); eventually follows from n-q_n<sqrt(n)',
    'H':'floor(X^(1/8))','Y':'sqrt(2mX)','MR_delta':'1/4',
    'resultant_bound':'abs R_i<=i^i*(2^i*(m*i)!)^(2i-1)',
    'bad_prime_union':'sum over i in I_H of primes p>Y dividing R_i, n=i modp; O_m(X^(7/8))',
    'error_terms':['X/(log X)^(1/50)','X^(799/800)*(log X)^(1/3)','X^(7/8)','X^(2/3)','X^(1/8)'],
    'dyadic_summation':'split at sqrt(X), freeze log on the larger blocks, bound the lower indices trivially',
    'old_density_one_theorem_used':False,'bootstrap_used':False,'parity_used':False,
    'all_offsets_squarefree_assumed':False,'all_C_smoothness_required':False,
    'uniform_in_growing_m_claimed':False,'numerical_constant_or_starting_threshold_extracted':False,
    'all_index_SFC_closed':False,'finite_controls_are_not_numerical_density_evidence':True,
    'source_pins':sources,'provenance':'User-supplied quantitative refinement checked and integrated; no novelty claim',
    'independent_replay':False,'formal_verification':False,'external_review':False}
(HERE/'quantitative_squarefree_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
(HERE/'quantitative_squarefree_examples.json').write_text(json.dumps({'shifted_profiles':profiles,'local_certificates':local,
    'resultant_certificates':resultants,'finite_bad_prime_fixture':{'exponent':m,'X':X,'Y_squared':Y2,'pairs':pairs,'union':sorted(union)},
    'finite_controls_are_not_numerical_density_evidence':True},indent=2)+'\n')
files=('QUANTITATIVE_SQUAREFREENESS_ALL_EXPONENTS.tex','quantitative_squarefree.py','verify_quantitative_squarefree.py',
       'quantitative_squarefree_audit.json','quantitative_squarefree_examples.json','prime_power_windows.py',*sources)
receipt={'status':'pass','checks':counts,'total':sum(counts.values()),
    'files':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
    'universal_evidence':'QUANTITATIVE_SQUAREFREENESS_ALL_EXPONENTS.tex and credited primary analytic inputs, not finite tests',
    'all_fixed_positive_integer_exponents':True,'quantitative_rate':'E_m(X)<<_m X/(log X)^(1/50)',
    'explicit_logarithmic_exponent':'1/50','numerical_constant_or_starting_threshold_extracted':False,
    'old_density_one_theorem_used':False,'bootstrap_used':False,'parity_used':False,
    'growing_offsets_bad_prime_union_controlled':True,'uniform_in_growing_m_claimed':False,
    'finite_controls_are_not_numerical_density_evidence':True,'all_index_SFC_closed':False,
    'independent_replay':False,'formal_verification':False,'external_review':False,'novelty_claim':False}
(HERE/'quantitative_squarefree_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'pass','checks':receipt['total'],'fixed_m':'every positive integer',
                  'rate':receipt['quantitative_rate'],'all_index_SFC_closed':False}))
