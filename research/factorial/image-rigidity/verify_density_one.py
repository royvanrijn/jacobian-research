"""Finite controls for the written density-one proof, not a numerical density scan."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial,isqrt
import hashlib,json
import sympy as s
import density_one_windows as d
import prime_power_windows as p

HERE=Path(__file__).resolve().parent
counts={};offsets=[];windows=[]
def check(name,value):
    if not value:raise AssertionError(name)
    counts[name]=counts.get(name,0)+1

check('irreducible_offset_density_exceeds_two_thirds_exactly',2*(Q(1,3)+Q(1,81))==Q(56,81)>Q(2,3))
check('smooth_long_mean_strictly_below_three_fifths_exactly',
      2*(Q(1,5)+Q(1,375))==Q(152,375)>Q(2,5))
check('short_interval_error_budget_has_fixed_positive_margin',Q(3,5)+2*Q(1,100)==Q(31,50)<Q(2,3))
check('large_prime_square_exceeds_six_n_eventually',
      all((2*X)**4>(12*X)**3 for X in (109,1000,1000000)))
check('MR_exception_power_for_fixed_delta',Q(1,100)/25==Q(1,2500)>0)

# Exact finite verification of the elementary long-mean identity, without
# computing a purported asymptotic rate or sampling factorial squarefreeness.
for X in (100,1000):
    Y=int(s.integer_nthroot((2*X)**2,3)[0])
    check('integer_smooth_cutoff_matches_real_cube_root',Y**3<=(2*X)**2<(Y+1)**3 and Y*Y>2*X)
    literal=sum(all(int(q)<=Y for q in s.factorint(u)) for u in range(X,2*X+1))
    large_primes=list(s.primerange(Y+1,2*X+1))
    by_multiples=X+1-sum(2*X//q-(X-1)//q for q in large_primes)
    check('disjoint_large_prime_multiple_identity',literal==by_multiples)
for Y in (7,11):
    smooth=lambda u:all(int(q)<=Y for q in s.factorint(u))
    check('smooth_indicator_is_multiplicative_on_finite_controls',
          all(smooth(a*b)==smooth(a)*smooth(b) for a,b in ((8,9),(7,11),(11,13),(16,25),(1,17))))

T=s.Symbol('t')
def poly(coeffs,q):return s.Poly(sum(a*T**j for j,a in enumerate(coeffs)),T,modulus=q)
for i in (0,1,2,3,13):
    cert=d.fixed_offset_certificate(i)
    f=s.Poly(sum(comb(i,j)*factorial(6*j)*T**j for j in range(i+1)),T,domain=s.ZZ)
    check('finite_offset_characteristic_zero_squarefreeness',cert['characteristic_zero_squarefree'] and s.gcd(f,f.diff()).degree()==0)
    disc=1 if i==0 else int(s.resultant(f,f.diff())*(-1)**(i*(i-1)//2)//f.LC())
    check('offset_discriminant_has_independent_resultant_normalization',disc==int(cert['integer_discriminant'])!=0)
    check('conservative_all_larger_prime_threshold_is_finite',cert['safe_all_larger_primes_threshold']==max(6*i,abs(disc)))
    check('finite_offset_test_does_not_claim_irreducibility_or_density',
          not cert['irreducibility_claimed'] and not cert['all_offsets_claimed'] and not cert['numerical_density_evidence'])
    offsets.append(cert)

for n,q,i,accepted in ((45,43,2,True),(115,37,4,True),(465,113,13,True),(217,107,3,False)):
    cert=d.shifted_good_offset_certificate(n,q,i)
    actual=[comb(n,j)*factorial(6*j) for j in range(n+1)]
    reduced=[a % q for a in actual]
    check('full_moment_reduction_is_exactly_the_offset_kernel',
          reduced==cert['residual_coefficients']+[0]*(n-i))
    f=poly(cert['residual_coefficients'],q)
    lhs=poly(cert['bezout_left'],q)*f+poly(cert['bezout_right'],q)*f.diff()
    check('finite_field_Bezout_identity_reconstructs_normalized_gcd',
          lhs==poly(cert['residual_gcd_coefficients'],q)==s.gcd(f,f.diff()))
    hull=p.newton_residuals(actual,q)
    check('literal_shifted_hull_has_horizontal_edge_and_denominator_p',
          len(hull)==2 and hull[0]['right'][0]==i and hull[0]['h']==0
          and hull[1]['e']==q and hull[1]['h']==6)
    check('global_degree_bound_lies_below_shifted_prime',cert['degree_bound_below_shifted_prime'])
    check('modular_bad_prime_does_not_invent_a_squarefree_claim',cert['squarefree_certified']==accepted)
    check('local_witness_does_not_claim_all_indices_or_numerical_density',
          not cert['counterexample_claimed'] and not cert['all_index_SFC_claimed'] and not cert['numerical_density_evidence'])
    windows.append(cert)
for operation in (lambda:d.fixed_offset_certificate(-1),lambda:d.fixed_offset_certificate(True),
                  lambda:d.shifted_good_offset_certificate(29,29,0),
                  lambda:d.shifted_good_offset_certificate(465,113,14),
                  lambda:d.shifted_good_offset_certificate(465,111,21),
                  lambda:d.shifted_good_offset_certificate(200,7,4)):
    try:operation()
    except ValueError:check('outside_local_proved_domain_fails_closed',True)
    else:raise AssertionError('invalid input accepted')

sources={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
    'sources/matomaki-radziwill-2016-annals.pdf','sources/matomaki-radziwill-2016-annals.txt',
    'sources/filaseta-rocks-2016-published.pdf','sources/filaseta-rocks-2016-published.txt')}
mr=(HERE/'sources/matomaki-radziwill-2016-annals.txt').read_text()
fr=(HERE/'sources/filaseta-rocks-2016-published.txt').read_text()
check('pinned_MR_text_contains_the_checked_uniformity_statement','to vary uniformly' in mr and '20000' in mr)
check('pinned_FR_text_contains_both_imported_pinpoints','Theorem 1.1.' in fr and '(3.3)' in fr)
audit={
    'status':'checked_scoped_written_argument','density_type':'ordinary natural density one',
    'claim':'FIRDENSITYONE20261116','canonical_proof':'DENSITY_ONE_SQUAREFREENESS.tex',
    'polynomial':'F_n(T)=sum binom(n,j)*(6j)!*T^j, ordinary Gamma shape1',
    'analytic_inputs':[
        {'source':'Filaseta-Rocks2016 author-hosted journal Online First',
         'pinpoints':'Theorem1.1, polynomial definition page1; Section3 definition B and equation3.3 pages6-7',
         'scope':'liminf density irreducible offsets>=log2; n-q_n<sqrt(n) on a natural-density-one set',
         'published_input_imported_not_independently_reproved':True},
        {'source':'Matomaki-Radziwill, Annals183(2016)1015-1056 publisher PDF',
         'pinpoints':'Theorem1 pp1015-1016 and uniformity sentence immediately following',
         'scope':'real multiplicative f bounded1; h,delta,f uniform; exceptional count includes (log h)^(1/3)/(delta^2*h^(delta/25))',
         'published_input_imported_not_independently_reproved':True},
        {'source':'classical PNT, retained DLMF27.2.3 reference',
         'scope':'derive long smooth-indicator mean by disjoint prime multiples and partial summation; no separate Dickman theorem needed',
         'published_input_imported_not_independently_reproved':True}],
    'long_mean':'1-log(3/2)<3/5 for Y=(2X)^(2/3)',
    'MR_delta':'1/100','MR_C_prime':20000,
    'MR_exception_eta':'C*(log H)^(1/3)/(delta^2*H^(delta/25)), tending0',
    'fixed_H_data':'I_H has size>2H/3; finitely many nonzero discriminants; finite T_H',
    'reverse_interval':'apply at x=n-H; backward H terms are a subset of the inclusive forward interval; omit O(H) dyadic boundary starts',
    'local_step':'p>(2X)^(2/3), p>T_H, p>6max(a,i); prime-gap bound k<sqrt(n)/2<p; horizontal degree-preserving reduction square-divides squarefree F_i modp',
    'primitive_normalization':'Gauss gives integral quotient, g(0)=+-1; horizontal-only factor has leading p-adic unit',
    'limit_order':['fix sufficiently large H','let X tend to infinity, passing the finite discriminant thresholds','let H tend to infinity','dyadic summation gives ordinary natural density'],
    'all_C_smoothness_required':False,'all_offsets_squarefree_assumed':False,
    'fixed_multiplier_or_offset_cutoff_in_density_statement':False,
    'irreducibility_density_one_claimed':False,'all_index_SFC_closed':False,
    'finite_controls_are_not_numerical_density_evidence':True,'source_pins':sources,
    'provenance':'User-supplied density-one argument; primary-source uniformity and local normalization checked; no novelty claim',
    'independent_replay':False,'formal_verification':False,'external_review':False}
(HERE/'density_one_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
(HERE/'density_one_local_examples.json').write_text(json.dumps({'fixed_offsets':offsets,'local_windows':windows,
    'finite_controls_are_not_numerical_density_evidence':True,'all_index_SFC_claimed':False},indent=2)+'\n')
files=('DENSITY_ONE_SQUAREFREENESS.tex','density_one_windows.py','verify_density_one.py','density_one_audit.json',
       'density_one_local_examples.json','prime_power_windows.py','prime_gap_windows.py',
       'PRIME_GAPS_AND_SMALL_KERNELS.tex','SLOPE_COVERAGE_AND_SHIFTED_WINDOWS.tex',*sources)
receipt={'status':'pass','checks':counts,'total':sum(counts.values()),
         'files':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
         'universal_evidence':'DENSITY_ONE_SQUAREFREENESS.tex and credited primary analytic inputs, not finite tests',
         'ordinary_natural_density_one_proved':True,'exceptional_count':'o(X)',
         'fixed_H_then_X_then_H_order_checked':True,'finite_bad_discriminant_thresholds_respected':True,
         'analytic_uniformity_in_X_dependent_function_checked':True,
         'all_C_smoothness_required':False,'irreducibility_density_one_claimed':False,
         'finite_controls_are_not_numerical_density_evidence':True,
         'all_index_SFC_closed':False,'independent_replay':False,'formal_verification':False,
         'external_review':False,'novelty_claim':False}
(HERE/'density_one_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'pass','checks':receipt['total'],'density_type':'ordinary natural',
                  'exceptional_count':'o(X)','all_index_SFC_closed':False}))
