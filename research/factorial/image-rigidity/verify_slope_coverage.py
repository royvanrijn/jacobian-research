"""Finite controls for the proved slope cap and shifted-index formulas."""
from pathlib import Path
from math import factorial,comb
import hashlib,json
import sympy as s
import newton_slope_coverage as c
import prime_power_windows as p
import mixed_prime_windows as m

HERE=Path(__file__).resolve().parent
counts={}
examples=[]
def check(group,condition):
    if not condition:
        raise AssertionError(group)
    counts[group]=counts.get(group,0)+1

for n in (2,7,21,23,31,64):
    for d in (1,2,6):
        coefficients=[comb(n,j)*factorial(d*j) for j in range(n+1)]
        for prime in (2,3,7,17,int(s.nextprime(d*n))):
            direct=p.newton_residuals(coefficients,prime)
            slopes=c.slope_edges(n,d,prime)
            check('Legendre_hull_matches_literal_factorial_coefficients',
                  slopes==[{k:row[k] for k in ('left','right','e','h')} for row in direct])
            residual=m.square_divisor_degrees(direct,prime)['possible_degrees']
            budget=c.slope_square_divisor_degrees(slopes)
            check('residual_multiplicities_refine_slope_budget',set(residual)<=set(budget))
            if prime>d*n:
                check('above_prime_cap_slope_budget_is_unrestricted',budget==list(range(n//2+1)))
# Above the cap a residual test remains useful, although slopes do not.
slopes=c.slope_square_divisor_degrees(c.slope_edges(7,6,43))
check('above_cap_residual_test_can_certify_while_slopes_admit_even_degree',
      2 in slopes and m.mixed_window_certificate(7,6,[43])['status']=='adjacent_window_certified')

records=[]
for n in range(2,65):
    cert=c.slope_coverage_certificate(n,6)
    check('bounded_full_prime_cap_slope_diagnostic_covers_supplied_index',
          cert['status']=='slope_parity_window_certified' and cert['all_primes_through_cap_checked']
          and not cert['all_starting_indices_claimed'] and not cert['residual_multiplicities_used'])
    check('diagnostic_really_exhausts_every_prime_to_cap',
          [row['prime'] for row in cert['profiles']]==list(s.primerange(2,6*n+1)))
    for k,prime in cert['excluding_prime_by_degree'].items():
        k=int(k)
        check('exported_candidate_exclusion_replays',
              k not in c.slope_square_divisor_degrees(c.slope_edges(n,6,prime)))
    records.append(cert)
for n in (7,21,23,28,64):
    for k in range(1,min(n//2,4)+1):
        shortlist=c.candidate_prime_set(n,k,6)
        for prime in s.primerange(2,6*n+1):
            if prime not in shortlist['candidate_primes']:
                coefficients=[comb(n,j)*factorial(6*j) for j in range(2*k+1)]
                check('omitted_prime_has_first_2k_coefficients_units',all(a % prime for a in coefficients))
                check('omitted_prime_cannot_exclude_candidate_by_slopes',
                      k in c.slope_square_divisor_degrees(c.slope_edges(n,6,int(prime))))
        for prime,offsets in shortlist['large_prime_divisor_offsets'].items():
            check('shortlisted_large_prime_has_literal_falling_factor_divisor',
                  int(prime)>12*k and all((n-i) % int(prime)==0 for i in offsets))
        examples.append(shortlist)

x=s.Symbol('x')
fixtures=[(6,1,19,1,i) for i in range(4)]
fixtures += [(6,1,19,2,i) for i in (0,2,3)]
fixtures += [(6,2,13,1,2),(6,3,19,1,3),(3,2,11,2,1),(1,3,5,2,2)]
for d,u,prime,r,i in fixtures:
    cert=c.shifted_index_profile(d,u,prime,r,i)
    n=cert['terminal_index']
    coefficients=[comb(n,j)*factorial(d*j) for j in range(n+1)]
    direct=p.newton_residuals(coefficients,prime)
    check('shifted_exact_polygon_and_both_normalized_residuals',
          cert['edges']==[{k:row[k] for k in ('left','right','e','h','residual_coefficients')} for row in direct])
    e=prime**r
    expected={h+e*t for h in range(i//2+1) for t in range(u//2+1)}
    check('shifted_doubled_occupancy_is_exact_degree_budget',cert['slope_only_degrees']==sorted(expected))
    family=c.shifted_family_certificate(d,u,prime,i)
    if family['all_powers_r_ge_1_certified']:
        check('accepted_family_matches_finite_exact_profile',cert['adjacent_window_certified'])
    if n<=24:
        f=s.Poly(sum(a*x**j for j,a in enumerate(coefficients)),x)
        check('separate_rational_gcd_in_shifted_family',s.gcd(f,f.diff()).degree()==0)
    if cert['adjacent_window_certified']:
        check('shifted_exclusion_contains_no_positive_even_degree',cert['possible_square_divisor_degrees']==[0])
    examples.append(cert)
for prime in (19,23,71):
    for i in range(4):
        family=c.shifted_family_certificate(6,1,prime,i)
        check('all_r_shifted_cubic_offset_family_exported_without_scan',
              family['status']=='unbounded_shifted_family_certified'
              and family['all_powers_r_ge_1_certified']
              and not family['bounded_search_used_for_universal_claim']
              and not family['all_starting_indices_claimed'])
        examples.append(family)
# This accepted family genuinely needs real parity at a repeated horizontal residual.
family=c.shifted_family_certificate(6,1,71,2)
check('shifted_family_handles_repeated_horizontal_residual_using_real_place',
      family['horizontal_square_degrees']==[0,1] and family['all_powers_r_ge_1_certified'])
# The former local survivor is closed by the full-square obstruction.
full=c.shifted_family_certificate(6,2,13,2)
check('full_square_obstruction_closes_former_local_survivor',
      full['raw_even_patterns']==[{'horizontal_degree':1,'nonhorizontal_residual_degree':1}]
      and full['full_square_removed_patterns']==full['raw_even_patterns']
      and full['all_powers_r_ge_1_certified'])
examples.append(full)
# A different small shifted kernel still survives below the uniform large band.
bad=c.shifted_family_certificate(6,4,137,3)
check('surviving_small_even_shifted_pattern_is_unknown',
      bad['status']=='shifted_kernel_test_unresolved' and not bad['counterexample_claimed']
      and bad['surviving_even_patterns']==[{'horizontal_degree':1,'nonhorizontal_residual_degree':1}]
      and not bad['large_degree_removed_patterns'] and not bad['full_square_removed_patterns'])
examples.append(bad)
# Independent small-kernel arithmetic uniformly removes the 2P+3 and 3P+2 patterns.
f2=p.moment_polynomial(2,6,1,x)
f3=p.moment_polynomial(3,6,1,x)
d2=s.discriminant(f2.as_expr(),x)
d3=s.discriminant(f3.as_expr(),x)
check('independent_F2_F3_discriminant_gcd',s.gcd(d2,d3)==2**10*3**4*5**2)
for f,disc in [(f2,d2),(f3,d3)]:
    resultant=s.polys.subresultants(f,f.diff())[-1].as_expr()
    check('discriminant_matches_separate_resultant_computation',
          s.simplify((-1)**(f.degree()*(f.degree()-1)//2)*resultant/f.LC()-disc)==0)
for n in (9,10,28,41,59,60,324):
    band=c.large_degree_certificate(n)
    prime=band['prime']
    literal=p.newton_residuals([comb(n,j)*factorial(6*j) for j in range(n+1)],prime)
    check('Nagura_prime_interval_is_literal',3*n<prime and 5*prime<18*n and s.isprime(prime))
    check('large_factor_two_edge_formula_matches_literal_coefficients',
          band['expected_slope_edges']==[{k:row[k] for k in ('left','right','e','h')} for row in literal])
    check('large_factor_degree_bound_is_strictly_below_three_tenths',10*band['square_divisor_degree_bound']<3*n)
    check('literal_slope_budget_has_the_uniform_large_factor_bound',
          max(c.slope_square_divisor_degrees(band['expected_slope_edges']))==band['square_divisor_degree_bound'])
    examples.append(band)
# Both user blocks, then the additional offset and multiplier extensions.
for u in range(1,7):
    for i in range(4 if u<=3 else 2):
        for prime in (13,19,71,107,137):
            if prime<=6*max(u,i):
                continue
            family=c.shifted_family_certificate(6,u,prime,i)
            check('enlarged_all_r_blocks_have_no_remaining_pattern',family['all_powers_r_ge_1_certified'])
            examples.append(family)
large=c.shifted_family_certificate(6,3,107,3)
check('uniform_large_factor_rule_closes_repeated_kernel_survivor',
      large['raw_even_patterns']==[{'horizontal_degree':1,'nonhorizontal_residual_degree':1}]
      and not large['full_square_removed_patterns']
      and large['large_degree_removed_patterns']==large['raw_even_patterns']
      and large['all_powers_r_ge_1_certified'])
for r in (1,2,4):
    prof=c.shifted_index_profile(6,3,107,r,3)
    check('accepted_symbolic_band_removal_matches_each_supplied_power',
          prof['adjacent_window_certified'] and prof['real_parity_filtered_degrees']!=[0]
          and prof['possible_square_divisor_degrees']==[0])
for n in range(2,41):
    check('ordinary_factorial_leading_coefficient_is_not_square',not s.integer_nthroot(factorial(6*n),2)[1])
for value in (8,9.0,True):
    try:
        c.large_degree_certificate(value)
    except ValueError:
        check('large_band_domain_fails_closed',True)
    else:
        raise AssertionError('invalid large-band input accepted')
for args in ((6,1,18,1,0),(6,2,11,1,0),(6,1,13,1,3),(6,1,19,0,1)):
    try:
        c.shifted_index_profile(*args)
    except ValueError:
        check('outside_shifted_theorem_domain_fails_closed',True)
    else:
        raise AssertionError('invalid shifted input accepted')
for operation in (lambda:c.slope_coverage_certificate(1),lambda:c.slope_coverage_certificate(True),
                  lambda:c.candidate_prime_set(7,4),lambda:c.slope_edges(7,6,4)):
    try:
        operation()
    except ValueError:
        check('invalid_finite_diagnostic_inputs_rejected',True)
    else:
        raise AssertionError('invalid finite diagnostic accepted')

probe={'scope':'FULL slope-only prime range p<=6n, plus real parity, at finite indices2..64',
       'all_indices_claimed':False,'residual_factorization_used':False,
       'survivor_indices':[row['terminal_index'] for row in records if row['surviving_even_degrees']],
       'certificates':records}
(HERE/'slope_coverage_probe.json').write_text(json.dumps(probe,indent=2)+'\n')
(HERE/'shifted_window_examples.json').write_text(json.dumps(examples,indent=2)+'\n')
files=('SLOPE_COVERAGE_AND_SHIFTED_WINDOWS.tex','newton_slope_coverage.py','verify_slope_coverage.py',
       'prime_power_windows.py','mixed_prime_windows.py','sources/nagura-1952-prime-interval.pdf',
       'sources/nagura-1952-prime-interval.txt')
receipt={'status':'pass','checks':counts,'total':sum(counts.values()),
         'files':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
         'universal_evidence':'SLOPE_COVERAGE_AND_SHIFTED_WINDOWS.tex; mathematical proofs, not finite scans',
         'finite_scope':'full slope cap at2..64, independent literal hulls/residuals, candidate-prime restrictions, shifted formulas and finite rational gcd controls',
         'infinite_shifted_family':'n=a*p^r+i, a1..3,i0..3,p>6*max(a,i),r>=1; also a4..6,i0..1,p>6a',
         'uniform_large_square_degree_bound':'degree<3*n/10 for EVERY n>=9, ordinary exponent6',
         'all_index_exponent_six_closed':False,'independent_replay':False,'external_review':False,
         'formal_verification':False,'novelty_claim':False}
(HERE/'slope_coverage_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
