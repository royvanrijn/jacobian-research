"""Exact exposed-edge and integer-subresultant controls; no all-index scan."""
from pathlib import Path
from math import comb,factorial,gcd,prod
import hashlib,json
import sympy as s
import prime_gap_windows as g
import prime_power_windows as p
import newton_slope_coverage as c

HERE=Path(__file__).resolve().parent
counts={};examples=[]
def check(name,value):
    if not value:raise AssertionError(name)
    counts[name]=counts.get(name,0)+1

for n in (30,31,37,42,64,100):
    coefficients=[comb(n,j)*factorial(6*j) for j in range(n+1)]
    for q in s.primerange(n//2+1,n+1):
        if q<=6 or q*q<=6*n:continue
        q=int(q);cert=g.exposed_edge(n,q)
        rows=p.newton_residuals(coefficients,q)
        check('exposed_edge_matches_independent_literal_hull',any(
              row['left']==cert['edge_left'] and row['right']==cert['edge_right']
              and row['e']==q and row['h']==6 for row in rows))
        check('complete_literal_slope_budget_respects_gap_bound',
              max(c.slope_square_divisor_degrees(rows))<=cert['square_divisor_degree_bound'])
        minimum=cert['minimum_adjusted_height_numerator']
        touched=[]
        for j,a in enumerate(coefficients):
            adjusted=q*p.valuation(a,q)-6*j
            check('all_literal_adjusted_heights_respect_support',adjusted>=minimum)
            if adjusted==minimum:touched.append(j)
        check('exactly_two_support_points_produce_primitive_edge',
              touched==[cert['support_index'],cert['support_index']+q])
        examples.append(cert)
for n in (30,31,35,37,42,64,100):
    cert=g.fourteenth_certificate(n)
    rows=p.newton_residuals([comb(n,j)*factorial(6*j) for j in range(n+1)],cert['prime'])
    check('both_fourteenth_polygon_branches_match_literal_coefficients',
          cert['expected_slope_edges']==[{k:row[k] for k in ('left','right','e','h')} for row in rows])
    check('strict_one_fourteenth_bound',14*cert['square_divisor_degree_bound']<n)
    check('fourteenth_budget_removes_both_nonhorizontal_edges',
          max(c.slope_square_divisor_degrees(rows))==cert['square_divisor_degree_bound'])
    examples.append(cert)
million=g.prime_gap_certificate(10**6)
check('million_index_prime_gap_example_without_large_polynomial',
      million['preceding_prime']==999983 and million['gap']==17
      and million['candidate_even_degrees']==[2,4,6,8])
examples.append(million)

x=s.Symbol('x');kernels=[]
for d in range(2,13):
    cert=g.small_kernel_certificate(d)
    coefficients=list(map(int,cert['coefficients']))
    f=s.Poly(sum(a*x**j for j,a in enumerate(coefficients)),x,domain=s.ZZ)
    derivative=f.diff()
    R,A,B,C=map(int,(cert['resultant'],cert['first_subresultant_linear'],cert['first_subresultant_constant'],cert['integer_gcd']))
    check('integer_gcd_Bezout_identity',sum(a*b for a,b in zip((R,A,B),map(int,cert['integer_gcd_bezout'])))==C)
    check('integer_gcd_matches_actual_three_integers',gcd(gcd(R,A),B)==C)
    check('smooth_factorization_reconstructs_C_exactly',
          prod(int(q)**power for q,power in cert['small_prime_factorization'].items())==C
          and all(s.isprime(int(q)) and int(q)<=6*d for q in cert['small_prime_factorization']))
    check('all_retained_kernel_smoothness_remainders_are_one',
          cert['status']=='kernel_square_degree_at_most_one_certified'
          and cert['remainder_after_small_primes']=='1')
    # Independent division-free cofactor construction of the first subresultant.
    columns=[s.Poly(x**j*f.as_expr(),x) for j in range(d-3,-1,-1)]
    columns += [s.Poly(x**j*derivative.as_expr(),x) for j in range(d-2,-1,-1)]
    high=s.Matrix([[column.nth(k) for column in columns] for k in range(2*d-3,1,-1)])
    vector=[(-1)**j*high[:,[k for k in range(high.cols) if k!=j]].det(method='domain-ge')
            for j in range(high.cols)] if high.rows else [s.Integer(1)]
    polynomial=s.Poly(sum(v*column.as_expr() for v,column in zip(vector,columns)),x)
    check('independent_integer_determinant_Bezout_reconstructs_Sres1',
          polynomial.as_expr() in (A*x+B,-A*x-B))
    # A separate full Sylvester determinant checks the resultant integer.
    sylvester_columns=[s.Poly(x**j*f.as_expr(),x) for j in range(d-2,-1,-1)]
    sylvester_columns += [s.Poly(x**j*derivative.as_expr(),x) for j in range(d-1,-1,-1)]
    sylvester=s.Matrix([[column.nth(k) for column in sylvester_columns] for k in range(2*d-2,-1,-1)])
    check('independent_full_Sylvester_resultant',abs(sylvester.det(method='domain-ge'))==abs(R))
    for q in (int(s.nextprime(6*d)),int(s.nextprime(12*d))):
        reduced=s.Poly(f.as_expr(),x,modulus=q)
        check('sample_modular_gcd_respects_universal_integer_certificate',s.gcd(reduced,reduced.diff()).degree()<=1)
    kernels.append(cert)
check('user_C2_exact_normalization',int(kernels[0]['integer_gcd'])==2**5*3**2*5)
check('user_C3_exact_normalization',int(kernels[1]['integer_gcd'])==2**30*3**17*5**6*7**3*11**2*13*17)
# Non-squarefree horizontal residuals are still covered by the stronger criterion.
for d,q in ((2,13),(3,107)):
    f=s.Poly(sum(comb(d,j)*factorial(6*j)*x**j for j in range(d+1)),x,modulus=q)
    check('residual_discriminant_zero_but_gcd_degree_one_is_allowed',s.gcd(f,f.diff()).degree()==1)
    n=q+d
    certificate=g.near_prime_window(n,q)
    literal=p.newton_residuals([comb(n,j)*factorial(6*j) for j in range(n+1)],q)
    check('near_prime_actual_horizontal_residual_matches_kernel',
          literal[0]['residual_coefficients']==[int(f.nth(j)) % q for j in range(d+1)]
          and certificate['squarefree_certified'])
for i in range(13):
    q=int(s.nextprime(6*max(1,i)))
    cert=g.shifted_family_certificate(1,q,i)
    check('all_prime_power_offsets_through_twelve_have_all_r_certificate',
          cert['all_powers_r_ge1_certified'] and not cert['all_indices_SFC_claimed'])
for a in range(1,15):
    for i in range(13):
        q=int(s.nextprime(6*max(a,i)))
        cert=g.shifted_family_certificate(a,q,i)
        check('combined_fourteen_by_thirteen_uniform_rectangle',cert['all_powers_r_ge1_certified'])
        if a in (1,14) and i in (0,12):examples.append(cert)
for a in range(1,28):
    for i in (0,1):
        q=int(s.nextprime(6*a))
        check('extended_offsets_zero_one_multiplier_blocks',g.shifted_family_certificate(a,q,i)['all_powers_r_ge1_certified'])
check('boundary_multiplier28_offset0_certified',g.shifted_family_certificate(28,173,0)['all_powers_r_ge1_certified'])
for a,q,i in ((28,173,1),(15,97,12)):
    cert=g.shifted_family_certificate(a,q,i)
    check('qualitative_eventual_statement_does_not_invent_all_r_or_cutoff',
          not cert['all_powers_r_ge1_certified'] and cert['eventual_fixed_multiplier_squarefree_theorem']
          and cert['effective_cutoff'] is None and not cert['counterexample_claimed'])
    examples.append(cert)
for operation in (lambda:g.prime_gap_certificate(29),lambda:g.fourteenth_certificate(29),
                  lambda:g.exposed_edge(100,47),lambda:g.exposed_edge(100,51),
                  lambda:g.small_kernel_certificate(13),lambda:g.small_kernel_certificate(-1),
                  lambda:g.shifted_family_certificate(1,73,13),lambda:g.shifted_family_certificate(2,11,0),
                  lambda:g.near_prime_window(100,83)):
    try:operation()
    except ValueError:check('outside_proved_input_domain_fails_closed',True)
    else:raise AssertionError('invalid input accepted')

(HERE/'small_kernel_certificates.json').write_text(json.dumps({'scope':'exact retained offsets2..12; universal q>6d, not all offsets',
    'certificates':kernels,'all_remainders_one':True,'all_offsets_claimed':False},indent=2)+'\n')
(HERE/'prime_gap_examples.json').write_text(json.dumps({'certificates':examples,'million_index_has_no_large_polynomial':True},indent=2)+'\n')
files=('PRIME_GAPS_AND_SMALL_KERNELS.tex','prime_gap_windows.py','verify_prime_gap_windows.py',
       'prime_power_windows.py','mixed_prime_windows.py','newton_slope_coverage.py',
       'sources/nagura-1952-prime-interval.pdf','sources/nagura-1952-prime-interval.txt')
receipt={'status':'pass','checks':counts,'total':sum(counts.values()),
         'files':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
         'universal_evidence':'PRIME_GAPS_AND_SMALL_KERNELS.tex; exposed-edge proofs, Nagura/PNT imports, integer Bezout/determinant certificates',
         'retained_kernel_offsets':list(range(2,13)),'all_kernel_remainders_one':True,
         'uniform_bound':'degree<n/14 for n>=30; preceding-prime even bound2*floor((n-q_n)/4)',
         'uniform_shifted_rectangle':'a1..14,i0..12,p>6max(a,i),r>=1',
         'eventual_fixed_multiplier_cutoff_computed':False,'all_offsets_claimed':False,
         'all_index_SFC_closed':False,'independent_replay':False,'formal_verification':False,
         'external_review':False,'novelty_claim':False}
(HERE/'prime_gap_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
